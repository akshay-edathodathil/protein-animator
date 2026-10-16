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
from density_prepare import CYPA_ACTIVE_SITE, SERIES_100K, SERIES_278K, prepare  # noqa: E402
from fetch_and_convert import ROOT, find_pymol  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", default="ProteinTimeSeries")
ap.add_argument("--ids", default=",".join(SERIES_278K), help="dataset names/PDB IDs in series order (default: CypA 278 K series)")
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
ap.add_argument("--frames-per-set", type=int, default=10)
ap.add_argument("--rotate", type=float, default=20.0, help="total slow rotation in degrees (0 = fixed view)")
ap.add_argument("--set-seconds", type=float, default=2.5, help="video seconds per dataset")
ap.add_argument("--manim-args", default="-r 1280,720 --fps 30")
ap.add_argument("--out", default="assets/density")
ap.add_argument("--no-video", action="store_true")
ap.add_argument("--prepare-only", action="store_true")
a = ap.parse_args()

out = Path(a.out).resolve()
data_dir = Path(a.dir).resolve()
norm = lambda x: x.strip().upper() if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", x.strip()) else x.strip()
ids = [norm(x) for x in a.ids.split(",")]
labels = [x.strip() for x in a.labels.split(",")] if a.labels else None
builtin = {tuple(SERIES_278K): "278 K", tuple(SERIES_100K): "100 K"}.get(tuple(ids))   # built-in CypA series
focus = a.focus or (CYPA_ACTIVE_SITE if builtin else "auto")
title = a.title or (f"Human cyclophilin A, {builtin}: data sets 1 to 9" if builtin else f"Electron density change across {len(ids)} datasets")
subtitle = a.subtitle or ("Refined models barely differ; the structure factors change across the series" if builtin == "278 K"
                          else "Density from each dataset's structure factors; difference density vs the first dataset")
zoom_label = a.zoom_label or ("Active site" if builtin else "Focus residues")
prepare(data_dir, ids, out, a.dmin, a.diff_dmin, focus, a.diff_nsigma, not a.no_download, labels)
if a.prepare_only:
    raise SystemExit(0)

pymol = find_pymol()
if not pymol:
    raise SystemExit("PyMOL not found (see README: ./setup.sh --with-pymol, or export PYMOL=/path/to/pymol)")
cfg = {"density_dir": str(out), "data_dir": str(data_dir), "frames_dir": str(out / "frames"), "size": a.size,
       "antialias": a.antialias, "ray_mode": a.ray_mode, "frames_per_set": a.frames_per_set, "rotate": a.rotate}
(out / "pymol_config.json").write_text(json.dumps(cfg))
meta = json.loads((out / "meta.json").read_text())
ensure_fresh(out / "frames", {**cfg, "ids": ids, "focus": [r["resi"] for r in meta["focus_residues"]],
                              "diff_contour": meta["diff_contour"], "dmin": a.dmin, "diff_dmin": a.diff_dmin})
print(f"Ray-tracing {len(ids) * a.frames_per_set * 2} panel images with PyMOL...")
subprocess.run([pymol, "-cq", str(Path(__file__).with_name("pymol_density.py")), "--", str(out / "pymol_config.json")], check=True)

if not a.no_video:
    env = dict(os.environ, DENSITY_DIR=str(out), DENSITY_SET_SECONDS=str(a.set_seconds),
               DENSITY_FRAMES_PER_SET=str(a.frames_per_set), DENSITY_TITLE=title,
               DENSITY_SUBTITLE=subtitle, DENSITY_ZOOM_LABEL=zoom_label)
    cmd = [sys.executable, "-m", "manim", *a.manim_args.split(), "--media_dir", str(out / "media"),
           str(ROOT / "scenes" / "density_scene.py"), "DensityScene"]
    print("Composing video with Manim...")
    subprocess.run(cmd, check=True, env=env)
    vids = sorted((out / "media").rglob("DensityScene.mp4"))
    print(f"\nVideo: {vids[-1] if vids else out / 'media'}")
