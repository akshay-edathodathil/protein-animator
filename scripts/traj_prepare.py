#!/usr/bin/env python3
"""Turn several PDB/mmCIF files (or multi-model files) into one aligned, interpolated trajectory.

Steps
  1. Every file (and every MODEL inside a file) is a *keyframe*, in natural filename order.
  2. Only atoms present in ALL keyframes are kept (matched by chain, residue number, residue name,
     atom name). Waters, hydrogens, ligands and alternate locations other than A are dropped.
  3. Each keyframe is superimposed (Kabsch, on CA atoms) onto the first one, so the camera can stay fixed.
  4. `--between N` interpolated frames are inserted between consecutive keyframes (straight-line
     Cartesian interpolation with smooth in/out easing).
  5. Writes <out>/morph.pdb (multi-model), <out>/disp.f32 (per-atom displacement from keyframe 0, Å,
     float32, frames x atoms) and <out>/meta.json (labels, keyframe flags, CA RMSD, ...).

Used by scripts/render_traj.py; can also be run on its own.
"""
import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

STRUCT_EXT = (".pdb", ".ent", ".cif", ".mmcif", ".pdb.gz", ".ent.gz", ".cif.gz")


def natural_key(p):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(p))]


def collect_inputs(paths):
    files = []
    for p in map(Path, paths):
        if p.is_dir():
            files += sorted((f for f in p.iterdir() if f.name.lower().endswith(STRUCT_EXT)), key=natural_key)
        elif p.exists():
            files.append(p)
        else:
            raise SystemExit(f"Input not found: {p}")
    if not files:
        raise SystemExit("No structure files found.")
    return files


@dataclass
class Keyframe:
    label: str
    atoms: dict  # key -> np.array(3)


def read_keyframes(files, chains=None):
    import gemmi

    frames = []
    for f in files:
        st = gemmi.read_structure(str(f))
        st.setup_entities()
        for mi, model in enumerate(st):
            atoms = {}
            for ch in model:
                if chains and ch.name not in chains:
                    continue
                for res in ch:
                    if res.entity_type != gemmi.EntityType.Polymer and res.het_flag != "A":
                        continue
                    for a in res:
                        if a.element.name == "H" or a.altloc not in ("\0", " ", "A"):
                            continue
                        key = (ch.name, res.seqid.num, res.seqid.icode.strip(), res.name, a.name)
                        atoms.setdefault(key, np.array([a.pos.x, a.pos.y, a.pos.z]))
            label = f.name.split(".")[0] if len(st) == 1 else f"{f.name.split('.')[0]} · model {mi + 1}"
            if atoms:
                frames.append(Keyframe(label, atoms))
    return frames


def kabsch(mobile, ref):
    """Rotation R and translations so that (mobile - mc) @ R.T + rc best fits ref."""
    mc, rc = mobile.mean(0), ref.mean(0)
    h = (mobile - mc).T @ (ref - rc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1, 1, d]) @ u.T
    return r, mc, rc


def prepare(files, between=24, labels=None, chains=None, align=True, ease=True):
    frames = read_keyframes(files, chains)
    if len(frames) < 2:
        raise SystemExit("Need at least 2 structures (files or models) to build a trajectory.")
    if labels:
        if len(labels) != len(frames):
            raise SystemExit(f"--labels has {len(labels)} entries but there are {len(frames)} structures.")
        for fr, lab in zip(frames, labels):
            fr.label = lab

    keys = [k for k in frames[0].atoms if all(k in fr.atoms for fr in frames[1:])]
    frac = len(keys) / len(frames[0].atoms)
    print(f"{len(frames)} structures; {len(keys)} atoms common to all ({frac:.0%} of the first structure)")
    if len(keys) < 20 or frac < 0.5:
        print("WARNING: few atoms in common. Check residue numbering/chain IDs (try --chains A).")
    ca = np.array([k[4] == "CA" for k in keys])
    if ca.sum() < 3:
        raise SystemExit("Fewer than 3 CA atoms in common; cannot align.")

    xyz = np.stack([np.array([fr.atoms[k] for k in keys]) for fr in frames])  # (K, N, 3)
    if align:
        ref_ca = xyz[0][ca]
        for i in range(1, len(xyz)):
            r, mc, rc = kabsch(xyz[i][ca], ref_ca)
            xyz[i] = (xyz[i] - mc) @ r.T + rc

    out_xyz, meta_frames = [], []
    for i in range(len(xyz) - 1):
        n_steps = between + 1
        for s in range(n_steps):
            t = s / n_steps
            if ease:
                t = t * t * (3 - 2 * t)
            out_xyz.append((1 - t) * xyz[i] + t * xyz[i + 1])
            meta_frames.append({"label": frames[i].label if s == 0 else f"{frames[i].label} → {frames[i + 1].label}",
                                "key": s == 0, "key_index": i, "t": (i + s / n_steps)})
    out_xyz.append(xyz[-1])
    meta_frames.append({"label": frames[-1].label, "key": True, "key_index": len(frames) - 1, "t": float(len(frames) - 1)})

    out_xyz = np.stack(out_xyz)
    disp = np.linalg.norm(out_xyz - out_xyz[0][None], axis=2).astype(np.float32)   # (F, N)
    ref_ca = out_xyz[0][ca]
    for k, mf in enumerate(meta_frames):
        mf["state"] = k + 1
        mf["rmsd"] = float(np.sqrt(((out_xyz[k][ca] - ref_ca) ** 2).sum(1).mean()))
    return keys, out_xyz, disp, meta_frames, [fr.label for fr in frames]


def _atom_line(serial, key, xyz, element):
    chain, resnum, icode, resname, name = key
    nm = f" {name:<3}" if len(element) == 1 and len(name) < 4 else f"{name:<4}"
    return (f"ATOM  {serial % 100000:>5} {nm} {resname:>3} {chain[:1]:1}{resnum:>4}{icode[:1]:1}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2}")


def write_outputs(out, keys, out_xyz, disp, meta_frames, key_labels, title):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    elements = [re.sub(r"[^A-Za-z]", "", k[4])[:1].upper() or "C" for k in keys]
    with open(out / "morph.pdb", "w") as fh:
        for fi, xyz in enumerate(out_xyz):
            fh.write(f"MODEL     {fi + 1:>4}\n")
            for si, (k, p, el) in enumerate(zip(keys, xyz, elements), 1):
                fh.write(_atom_line(si, k, p, el) + "\n")
            fh.write("ENDMDL\n")
        fh.write("END\n")
    disp.astype("<f4").tofile(out / "disp.f32")
    vmax = float(np.percentile(disp, 99))
    meta = {"title": title, "n_atoms": len(keys), "n_frames": len(out_xyz), "keyframes": key_labels,
            "disp_max": max(vmax, 0.5), "frames": meta_frames}
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    print(f"Wrote {out}/morph.pdb ({len(out_xyz)} states, {len(keys)} atoms), disp.f32, meta.json "
          f"(max CA RMSD vs first: {max(f['rmsd'] for f in meta_frames):.2f} Å)")
    return meta


def add_prepare_args(ap):
    ap.add_argument("inputs", nargs="+", help="structure files, directories, or multi-model files (in time order)")
    ap.add_argument("--between", type=int, default=24, help="interpolated frames between keyframes (0 = none)")
    ap.add_argument("--labels", help="comma-separated labels, one per structure (default: file names)")
    ap.add_argument("--chains", help="comma-separated chain IDs to keep (default: all polymer chains)")
    ap.add_argument("--no-align", action="store_true", help="do not superimpose structures onto the first")
    ap.add_argument("--no-ease", action="store_true", help="linear instead of smooth in/out interpolation")
    ap.add_argument("--title", default="Protein conformational change")
    ap.add_argument("--out", default="assets/trajectory")


def run_prepare(a):
    files = collect_inputs(a.inputs)
    labels = [s.strip() for s in a.labels.split(",")] if a.labels else None
    chains = set(a.chains.split(",")) if a.chains else None
    keys, xyz, disp, mf, kl = prepare(files, a.between, labels, chains, not a.no_align, not a.no_ease)
    return write_outputs(a.out, keys, xyz, disp, mf, kl, a.title)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_prepare_args(ap)
    run_prepare(ap.parse_args())
