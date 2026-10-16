#!/usr/bin/env python3
"""End-to-end: X-ray datasets (model + structure factors) -> density maps -> PyMOL frames -> Manim video.

    python scripts/render_density.py                    # CypA 278 K series 5KUZ,5KV0..5KV7 (downloads missing files)
    python scripts/render_density.py --ids 5KUL,5KUN,5KUO,5KUQ,5KUR,5KUS,5KUU,5KUV,5KUW --title "CypA, 100 K"

Files are looked up in --dir (default ProteinTimeSeries/) as <ID>.cif and <ID>-sf.cif and downloaded from RCSB if missing.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from density_prepare import CYPA_ACTIVE_SITE, SERIES_278K, prepare  # noqa: E402
from fetch_and_convert import ROOT, find_pymol  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", default="ProteinTimeSeries")
ap.add_argument("--ids", default=",".join(SERIES_278K), help="PDB IDs in series order")
ap.add_argument("--title", default="Human cyclophilin A, 278 K: data sets 1 to 9")
ap.add_argument("--focus", default=CYPA_ACTIVE_SITE, help="residue numbers for the zoom panel, or 'auto'")
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
ids = [s.strip().upper() for s in a.ids.split(",")]
prepare(data_dir, ids, out, a.dmin, a.diff_dmin, a.focus, a.diff_nsigma)
if a.prepare_only:
    raise SystemExit(0)

pymol = find_pymol()
if not pymol:
    raise SystemExit("PyMOL not found (see README: ./setup.sh --with-pymol, or export PYMOL=/path/to/pymol)")
cfg = {"density_dir": str(out), "data_dir": str(data_dir), "frames_dir": str(out / "frames"), "size": a.size,
       "antialias": a.antialias, "ray_mode": a.ray_mode, "frames_per_set": a.frames_per_set, "rotate": a.rotate}
(out / "pymol_config.json").write_text(json.dumps(cfg))
print(f"Ray-tracing {len(ids) * a.frames_per_set * 2} panel images with PyMOL...")
subprocess.run([pymol, "-cq", str(Path(__file__).with_name("pymol_density.py")), "--", str(out / "pymol_config.json")], check=True)

if not a.no_video:
    env = dict(os.environ, DENSITY_DIR=str(out), DENSITY_SET_SECONDS=str(a.set_seconds),
               DENSITY_FRAMES_PER_SET=str(a.frames_per_set), DENSITY_TITLE=a.title)
    cmd = [sys.executable, "-m", "manim", *a.manim_args.split(), "--media_dir", str(out / "media"),
           str(ROOT / "scenes" / "density_scene.py"), "DensityScene"]
    print("Composing video with Manim...")
    subprocess.run(cmd, check=True, env=env)
    vids = sorted((out / "media").rglob("DensityScene.mp4"))
    print(f"\nVideo: {vids[-1] if vids else out / 'media'}")
