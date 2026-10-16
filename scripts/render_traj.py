#!/usr/bin/env python3
"""End-to-end: several PDB files -> aligned/interpolated trajectory -> PyMOL frames -> Manim video.

    python scripts/render_traj.py examples/adenylate_kinase/4AKE.pdb examples/adenylate_kinase/1AKE.pdb \\
        --labels "Open (4AKE),Closed (1AKE)" --chains A --title "Adenylate kinase: lid closing"

Inputs may be files, a directory (natural filename order), or multi-model files (NMR/MD snapshots).
Add --no-video to stop after the frames; --prepare-only to stop after the alignment/interpolation.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_and_convert import ROOT, find_pymol  # noqa: E402
from frame_cache import ensure_fresh  # noqa: E402
from traj_prepare import add_prepare_args, run_prepare  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
add_prepare_args(ap)
ap.add_argument("--width", type=int, default=1280)
ap.add_argument("--height", type=int, default=720)
ap.add_argument("--antialias", type=int, default=2, help="1 = fast, 2 = smoother")
ap.add_argument("--ray-mode", type=int, default=1, help="PyMOL ray_trace_mode: 0 plain, 1 black outlines")
ap.add_argument("--color-by", choices=["displacement", "chain", "none"], default="displacement")
ap.add_argument("--disp-max", type=float, default=None, help="Å mapped to full red (default: 99th percentile)")
ap.add_argument("--rotate", type=float, default=20.0, help="total slow camera rotation in degrees (0 = locked view)")
ap.add_argument("--ss", choices=["first", "each"], default="first",
                help="secondary structure: assigned once (stable) or recomputed for every frame")
ap.add_argument("--play-fps", type=float, default=15, help="video speed through interpolated frames")
ap.add_argument("--hold", type=float, default=1.2, help="seconds to pause on each real structure")
ap.add_argument("--manim-args", default="-r 1280,720 --fps 30", help="flags passed to manim")
ap.add_argument("--no-video", action="store_true")
ap.add_argument("--prepare-only", action="store_true")
a = ap.parse_args()

out = Path(a.out)
if not out.is_absolute():
    out = Path.cwd() / out
a.out = str(out)

meta = run_prepare(a)
meta["color_by"] = a.color_by
(out / "meta.json").write_text(json.dumps(meta, indent=1))
if a.prepare_only:
    raise SystemExit(0)

pymol = find_pymol()
if not pymol:
    raise SystemExit("PyMOL not found (see README: ./setup.sh --with-pymol, or export PYMOL=/path/to/pymol)")
cfg = {"traj_dir": str(out), "frames_dir": str(out / "frames"), "width": a.width, "height": a.height,
       "antialias": a.antialias, "ray_mode": a.ray_mode, "color_by": a.color_by, "ss": a.ss,
       "rotate": a.rotate, "disp_max": a.disp_max}
(out / "pymol_config.json").write_text(json.dumps(cfg))
ensure_fresh(out / "frames", {**cfg, "n_frames": meta["n_frames"], "frames": [round(f["rmsd"], 3) for f in meta["frames"]],
                              "disp_max": meta["disp_max"]})
print(f"Ray-tracing {meta['n_frames']} frames with PyMOL...")
subprocess.run([pymol, "-cq", str(Path(__file__).with_name("pymol_traj.py")), "--", str(out / "pymol_config.json")], check=True)

if not a.no_video:
    env = dict(os.environ, TRAJ_DIR=str(out), TRAJ_PLAY_FPS=str(a.play_fps), TRAJ_HOLD=str(a.hold))
    cmd = [sys.executable, "-m", "manim", *a.manim_args.split(), "--media_dir", str(out / "media"),
           str(ROOT / "scenes" / "trajectory_scene.py"), "TrajectoryScene"]
    print("Composing video with Manim...")
    subprocess.run(cmd, check=True, env=env)
    vids = sorted((out / "media").rglob("TrajectoryScene.mp4"))
    print(f"\nVideo: {vids[-1] if vids else out / 'media'}")
