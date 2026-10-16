#!/usr/bin/env python3
"""Electron-density time series from a set of X-ray datasets (model + structure factors per dataset).

For each dataset (e.g. the CypA 278 K series 5KUZ, 5KV0 ... 5KV7) this:
  1. reads the structure factors (intensities, F_calc, phase_calc) and the model,
  2. converts I -> Fo = sqrt(max(I, 0)) and scales Fo onto |Fc| (resolution-dependent: ln k + b/d^2),
  3. computes a 2Fo-Fc map (phases from the deposited model; no sigma-A weighting),
  4. computes an isomorphous difference map  Fo(n) - Fo(1)  with the phases of dataset 1, so changes
     in the *data* show up even though the refined models barely differ,
  5. writes CCP4 maps (cropped around the molecule), per-dataset statistics (R_iso etc.),
     and ranks residues by how much difference density sits on them (to choose a zoom region).

Output (default assets/density/): maps/<ID>_2fofc.ccp4, maps/<ID>_diff.ccp4, meta.json
Used by scripts/render_density.py.
"""
import argparse
import json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np

SERIES_278K = ["5KUZ", "5KV0", "5KV1", "5KV2", "5KV3", "5KV4", "5KV5", "5KV6", "5KV7"]
SERIES_100K = ["5KUL", "5KUN", "5KUO", "5KUQ", "5KUR", "5KUS", "5KUU", "5KUV", "5KUW"]


def fetch(data_dir: Path, pdb_id: str):
    data_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("", "-sf"):
        f = data_dir / f"{pdb_id}{suffix}.cif"
        if not f.exists():
            print(f"Downloading {f.name}")
            urlretrieve(f"https://files.rcsb.org/download/{f.name}", f)
    return data_dir / f"{pdb_id}.cif", data_dir / f"{pdb_id}-sf.cif"


def load_dataset(model_path: Path, sf_path: Path):
    import gemmi

    st = gemmi.read_structure(str(model_path))
    st.setup_entities()
    rb = gemmi.as_refln_blocks(gemmi.cif.read(str(sf_path)))[0]
    mtz = gemmi.CifToMtz().convert_block_to_mtz(rb)
    cols = mtz.column_labels()
    a = mtz.array
    need = ("IMEAN", "SIGIMEAN", "FC", "PHIC")
    if any(c not in cols for c in need):
        raise SystemExit(f"{sf_path.name}: needs intensities + F_calc/phase_calc columns, found {cols}")
    d = {c: a[:, cols.index(c)].astype(float) for c in cols if c not in ("H", "K", "L")}
    hkl = a[:, :3].round().astype(int)
    cell = st.cell if st.cell.is_crystal() else mtz.cell
    inv_d2 = np.array([cell.calculate_1_d2(h) for h in hkl.tolist()])
    return {"st": st, "cell": cell, "sg": mtz.spacegroup, "hkl": hkl, "I": d["IMEAN"], "sig": d["SIGIMEAN"],
            "fc": d["FC"], "phic": d["PHIC"], "inv_d2": inv_d2}


def scaled_fo(ds):
    """Fo scaled onto the |Fc| scale with a resolution-dependent factor exp(a + b/d^2)."""
    fo = np.sqrt(np.clip(ds["I"], 0, None))
    ok = (fo > 0) & np.isfinite(ds["fc"]) & (ds["fc"] > 0)   # some reflections have no F_calc ('?')
    bins = np.quantile(ds["inv_d2"][ok], np.linspace(0, 1, 16))
    xs, ys = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = ok & (ds["inv_d2"] >= lo) & (ds["inv_d2"] <= hi)
        if m.sum() > 20:
            xs.append(ds["inv_d2"][m].mean())
            ys.append(np.log(ds["fc"][m].sum() / fo[m].sum()))
    b, a = np.polyfit(xs, ys, 1)
    return fo * np.exp(a + b * ds["inv_d2"])


def make_grid(ds, hkl, amp, phase_deg, sample_rate=3.0):
    import gemmi

    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = ds["sg"]
    mtz.set_cell_for_all(ds["cell"])
    mtz.add_dataset("d")
    mtz.add_column("F", "F")
    mtz.add_column("PHI", "P")
    mtz.set_data(np.column_stack([hkl, amp, phase_deg]).astype(np.float32))
    return mtz.transform_f_phi_to_map("F", "PHI", sample_rate=sample_rate)


def write_map(grid, st, path, margin=6.0):
    import gemmi

    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header()
    ccp4.set_extent(st.calculate_fractional_box(margin))
    ccp4.write_ccp4_map(str(path))


def residue_scores(st, diff_grids, noise):
    """Per residue: mean over difference maps of the strongest |difference density| (in noise sigmas) at its atoms."""
    scores = {}
    for g, nz in zip(diff_grids, noise):
        for ch in st[0]:
            for res in ch:
                if res.is_water() or res.het_flag != "A":
                    continue
                peak = max(abs(g.interpolate_value(a.pos)) for a in res) / nz
                scores.setdefault((ch.name, res.seqid.num, res.name), []).append(float(peak))
    return sorted(((float(np.mean(v)), k) for k, v in scores.items()), reverse=True)


# CypA active site / catalytic network (human numbering); override with --focus for other proteins
CYPA_ACTIVE_SITE = "55,60,61,63,72,101,102,103,111,113,121,122,126"


def prepare(data_dir, ids, out, dmin, diff_dmin, focus=CYPA_ACTIVE_SITE, nsigma=3.0, fetch_missing=True):
    import gemmi

    out = Path(out)
    (out / "maps").mkdir(parents=True, exist_ok=True)
    data = []
    for i in ids:
        mp, sp = fetch(Path(data_dir), i) if fetch_missing else (Path(data_dir) / f"{i}.cif", Path(data_dir) / f"{i}-sf.cif")
        print(f"Loading {i}")
        ds = load_dataset(mp, sp)
        ds["id"] = i
        ds["fo_s"] = scaled_fo(ds)
        data.append(ds)

    ref = data[0]
    key = lambda h: (int(h[0]), int(h[1]), int(h[2]))
    ref_index = {key(h): j for j, h in enumerate(ref["hkl"])}
    sigma_ref, grids2, diffs = None, [], []
    stats = []
    for n, ds in enumerate(data):
        sel = (ds["inv_d2"] <= 1.0 / dmin**2) & np.isfinite(ds["fc"]) & np.isfinite(ds["phic"]) & np.isfinite(ds["fo_s"])
        amp = 2 * ds["fo_s"] - ds["fc"]
        g2 = make_grid(ds, ds["hkl"][sel], amp[sel], ds["phic"][sel])
        arr = np.array(g2, copy=False)
        if n == 0:
            sigma_ref = float(arr.std())     # common density unit for ALL maps: std of set 1's 2Fo-Fc
        arr[:] = (arr - arr.mean()) / sigma_ref
        grids2.append(g2)

        with np.errstate(divide="ignore", invalid="ignore"):
            isig = np.where(ds["sig"] > 0, ds["I"] / ds["sig"], np.nan)
        st_ = {"id": ds["id"], "set": n + 1, "cell": [round(x, 3) for x in ds["cell"].parameters[:3]],
               "n_refl": int(len(ds["I"])), "frac_I_over_sig2": float(np.nanmean(isig > 2)),
               "mean_I_over_sig": float(np.nanmean(isig)),
               "n_atoms": int(sum(len(r) for ch in ds["st"][0] for r in ch))}
        if n > 0:
            idx = np.array([ref_index.get(key(h), -1) for h in ds["hkl"]])
            ok = (idx >= 0) & (ds["inv_d2"] <= 1.0 / diff_dmin**2) & np.isfinite(ds["fo_s"])
            ok[ok] &= np.isfinite(ref["phic"][idx[ok]]) & np.isfinite(ref["fo_s"][idx[ok]])
            ri = idx[ok]
            f_n = ds["fo_s"][ok]
            f_1 = ref["fo_s"][ri]
            # bring Fo(n) onto Fo(1)'s scale (single scalar over the common reflections)
            f_n = f_n * (f_1.sum() / f_n.sum())
            st_["r_iso"] = float(np.abs(f_n - f_1).sum() / f_1.sum())
            gd = make_grid(ref, ref["hkl"][ri], f_n - f_1, ref["phic"][ri])
            darr = np.array(gd, copy=False)
            darr[:] = (darr - darr.mean()) / sigma_ref
            st_["diff_rms"] = float(darr.std())
            diffs.append(gd)
        stats.append(st_)

    # One fixed contour for every frame: +/- nsigma x sigma of the LAST (largest-change) difference map.
    contour = nsigma * stats[-1]["diff_rms"]
    for n, ds in enumerate(data):
        write_map(grids2[n], ds["st"], out / "maps" / f"{ds['id']}_2fofc.ccp4")
        if n > 0:
            write_map(diffs[n - 1], ref["st"], out / "maps" / f"{ds['id']}_diff.ccp4")
    ranked = residue_scores(ref["st"], diffs, [s["diff_rms"] for s in stats[1:]])
    score_of = {k[1]: (k[0], k[2], sc) for sc, k in ranked}
    if focus == "auto":
        nums = [k[1] for _, k in ranked[:8]]
    else:
        nums = [int(x) for x in focus.split(",") if x.strip()]
    chain = ref["st"][0][0].name
    focus_res = [{"chain": chain, "resi": n, "resn": score_of[n][1], "score": round(score_of[n][2], 2)}
                 for n in nums if n in score_of]
    if not focus_res:
        raise SystemExit("None of the --focus residues exist in the model.")

    meta = {"ids": ids, "dmin": dmin, "diff_dmin": diff_dmin, "sigma_ref": sigma_ref,
            "diff_contour": round(contour, 4), "diff_nsigma": nsigma, "twofofc_contour": 1.0, "datasets": stats,
            "focus_residues": focus_res, "top_ranked": [{"resi": k[1], "resn": k[2], "score": round(sc, 2)} for sc, k in ranked[:10]]}
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    print(f"\nR_iso vs {ids[0]}:  " + "  ".join(f"{s['id']}={s['r_iso']:.3f}" for s in stats[1:]))
    print(f"Difference-map rms (units of set-1 2Fo-Fc sigma): " + "  ".join(f"{s['diff_rms']:.3f}" for s in stats[1:]))
    print("Zoom focus:", ", ".join(f"{f['resn']}{f['resi']}" for f in focus_res))
    return meta


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="ProteinTimeSeries", help="folder with <ID>.cif and <ID>-sf.cif (downloaded if missing)")
    ap.add_argument("--ids", default=",".join(SERIES_278K), help="comma-separated PDB IDs in series order")
    ap.add_argument("--dmin", type=float, default=1.8, help="resolution limit for 2Fo-Fc maps (Å)")
    ap.add_argument("--diff-dmin", type=float, default=2.0, help="resolution limit for difference maps (Å)")
    ap.add_argument("--focus", default=CYPA_ACTIVE_SITE,
                    help="residue numbers for the zoom panel (default: CypA active site) or 'auto' (strongest difference peaks)")
    ap.add_argument("--diff-nsigma", type=float, default=3.0,
                    help="difference-map contour = this x sigma of the last difference map (same for all frames)")
    ap.add_argument("--out", default="assets/density")
    a = ap.parse_args()
    prepare(a.dir, [s.strip().upper() for s in a.ids.split(",")], a.out, a.dmin, a.diff_dmin, a.focus, a.diff_nsigma)
