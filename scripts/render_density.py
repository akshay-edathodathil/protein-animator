#!/usr/bin/env python3
"""End-to-end: X-ray datasets (model + structure factors) -> density maps -> PyMOL frames -> Manim video.

    python scripts/render_density.py                    # CypA 278 K series 5KUZ,5KV0..5KV7 (downloads missing files)
    python scripts/render_density.py --ids 5KUL,5KUN,5KUO,5KUQ,5KUR,5KUS,5KUU,5KUV,5KUW --title "CypA, 100 K"

For each dataset ID the folder --dir (default ProteinTimeSeries/) must hold a model (<ID>.cif/.pdb) and structure
factors (<ID>-sf.cif or <ID>.mtz). Files for real PDB IDs are downloaded from RCSB if missing.
Your own data: see docs/CUSTOM_DATA.md.
"""
import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from frame_cache import ensure_fresh  # noqa: E402
from density_prepare import CYPA_ACTIVE_SITE, SERIES_100K, SERIES_278K, discover, prepare  # noqa: E402
from fetch_and_convert import ROOT, find_pymol  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", default=None,
                help="folder with your data. With --dir alone, model + reflection pairs are found automatically")
ap.add_argument("--ids", default=None,
                help="dataset names/PDB IDs in series order (default: discovered in --dir; the built-in CypA 278 K series if no --dir)")
ap.add_argument("--reverse", action="store_true", help="reverse the automatically discovered order")
ap.add_argument("--title", default=None, help="video title (default: generic, or CypA text for the built-in series)")
ap.add_argument("--labels", help="comma-separated label per dataset, e.g. \"0 min,5 min,10 min\" (default: Data set k)")
ap.add_argument("--subtitle", default=None, help="line under the title")
ap.add_argument("--zoom-label", default=None, help="caption of the zoom panel (default: 'Focus residues')")
ap.add_argument("--focus", default=None,
                help="residue numbers for the zoom panel, e.g. 55,60,113, or 'auto' (strongest difference peaks; default for your own data)")
ap.add_argument("--no-download", action="store_true", help="never download missing files from RCSB")
ap.add_argument("--dmin", type=float, default=1.8)
ap.add_argument("--diff-dmin", type=float, default=2.0)
ap.add_argument("--diff-nsigma", type=float, default=3.0)
ap.add_argument("--size", type=int, default=900, help="square panel size in pixels")
ap.add_argument("--antialias", type=int, default=2)
ap.add_argument("--ray-mode", type=int, default=0)
ap.add_argument("--frames-per-set", type=int, default=10, help="frames ray-traced per dataset (each frame = 2 panel images)")
ap.add_argument("--total-frames", type=int, default=None,
                help="instead of --frames-per-set: about this many frames in total, spread evenly over the datasets")
ap.add_argument("--panels", choices=["both", "zoom", "overview"], default="both",
                help="which panel(s) to render: both (default), zoom (only the zoom on the focus residues, shown larger), overview")
ap.add_argument("--fast", action="store_true",
                help="draft preset for slow machines: --size 500 --antialias 1 --frames-per-set 2 (other options still override)")
ap.add_argument("--rotate", type=float, default=20.0, help="total slow rotation in degrees (0 = fixed view)")
ap.add_argument("--set-seconds", type=float, default=2.5, help="video seconds per dataset")
ap.add_argument("--manim-args", default="-r 1280,720 --fps 30")
ap.add_argument("--out", default="assets/density")
ap.add_argument("--no-video", action="store_true")
ap.add_argument("--prepare-only", action="store_true")
pre, _ = ap.parse_known_args()
if pre.fast:                                  # preset values become the defaults; explicit options still win
    ap.set_defaults(size=500, antialias=1, frames_per_set=2)
a = ap.parse_args()

out = Path(a.out).resolve()
norm = lambda x: x.strip().upper() if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", x.strip()) else x.strip()
files = None
if a.ids:                                    # explicit list of dataset names
    data_dir = Path(a.dir or "ProteinTimeSeries").resolve()
    ids = [norm(x) for x in a.ids.split(",")]
elif a.dir:                                  # only --dir: find the model + reflection pairs ourselves
    data_dir = Path(a.dir).resolve()
    found = discover(data_dir, a.reverse)
    ids = [k for k, _, _ in found]
    files = {k: (m, r) for k, m, r in found}
    print(f"Found {len(found)} datasets in {data_dir} (series order; the first is the reference):")
    for n, (k, m, r) in enumerate(found, 1):
        print(f"  {n:2d}. {k:20s} model: {m.name:28s} reflections: {r.name}")
else:                                        # nothing given: built-in CypA 278 K demo series
    data_dir = Path("ProteinTimeSeries").resolve()
    ids = list(SERIES_278K)
labels = [x.strip() for x in a.labels.split(",")] if a.labels else None
builtin = {tuple(SERIES_278K): "278 K", tuple(SERIES_100K): "100 K"}.get(tuple(ids))   # built-in CypA series
focus = a.focus or (CYPA_ACTIVE_SITE if builtin else "auto")
title = a.title or (f"Human cyclophilin A, {builtin}: data sets 1 to 9" if builtin else f"Electron density change across {len(ids)} datasets")
subtitle = a.subtitle or ("Refined models barely differ; the structure factors change across the series" if builtin == "278 K"
                          else "Density from each dataset's structure factors; difference density vs the first dataset")
zoom_label = a.zoom_label or ("Active site" if builtin else "Focus residues")
prepare(data_dir, ids, out, a.dmin, a.diff_dmin, focus, a.diff_nsigma, not a.no_download, labels, files)
if a.prepare_only:
    raise SystemExit(0)

pymol = find_pymol()
if not pymol:
    raise SystemExit("PyMOL not found (see README: ./setup.sh --with-pymol, or export PYMOL=/path/to/pymol)")
if a.total_frames:
    a.frames_per_set = max(1, round(a.total_frames / len(ids)))
panels = ["overview", "zoom"] if a.panels == "both" else [a.panels]
print(f"{len(ids)} datasets x {a.frames_per_set} frame(s) = {len(ids) * a.frames_per_set} frames "
      f"({len(ids) * a.frames_per_set * len(panels)} panel images: {', '.join(panels)}) at {a.size} px, antialias {a.antialias}")
cfg = {"density_dir": str(out), "data_dir": str(data_dir), "frames_dir": str(out / "frames"), "size": a.size,
       "antialias": a.antialias, "ray_mode": a.ray_mode, "frames_per_set": a.frames_per_set, "rotate": a.rotate,
       "panels": panels}
(out / "pymol_config.json").write_text(json.dumps(cfg))
meta = json.loads((out / "meta.json").read_text())
ensure_fresh(out / "frames", {**cfg, "ids": ids, "focus": [r["resi"] for r in meta["focus_residues"]],
                              "diff_contour": meta["diff_contour"], "dmin": a.dmin, "diff_dmin": a.diff_dmin})
print("Ray-tracing with PyMOL (progress lines are per frame)...")
subprocess.run([pymol, "-cq", str(Path(__file__).with_name("pymol_density.py")), "--", str(out / "pymol_config.json")], check=True)

if not a.no_video:
    env = dict(os.environ, DENSITY_DIR=str(out), DENSITY_SET_SECONDS=str(a.set_seconds),
               DENSITY_FRAMES_PER_SET=str(a.frames_per_set), DENSITY_TITLE=title, DENSITY_PANELS=a.panels,
               DENSITY_SUBTITLE=subtitle, DENSITY_ZOOM_LABEL=zoom_label)
    cmd = [sys.executable, "-m", "manim", *a.manim_args.split(), "--media_dir", str(out / "media"),
           str(ROOT / "scenes" / "density_scene.py"), "DensityScene"]
    print("Composing video with Manim...")
    subprocess.run(cmd, check=True, env=env)
    vids = sorted((out / "media").rglob("DensityScene.mp4"))
    print(f"\nVideo: {vids[-1] if vids else out / 'media'}")
