#!/usr/bin/env python3
"""Launcher: ray-trace a transparent PNG turntable of the protein with PyMOL.

    python scripts/render_pymol_frames.py --frames 120 --width 1280 --height 720

Output: assets/pymol_frames/protein_frame0001.png ... (resumable: existing frames are skipped).
Ray tracing is slow (tens of seconds per frame for 6QL9); lower --frames/--width to go faster.
"""
import argparse
import subprocess
from pathlib import Path

from fetch_and_convert import ROOT, find_pymol

ap = argparse.ArgumentParser()
ap.add_argument("--pdb", default="6QL9")
ap.add_argument("--frames", type=int, default=120)
ap.add_argument("--width", type=int, default=1280)
ap.add_argument("--height", type=int, default=720)
ap.add_argument("--antialias", type=int, default=2, help="1 = fast, 2 = smoother (4x rays)")
ap.add_argument("--surface-quality", type=int, default=0, help="-1 fewer triangles (faster), 0 default")
ap.add_argument("--out", default=str(ROOT / "assets" / "pymol_frames"))
a = ap.parse_args()

exe = find_pymol()
if not exe:
    raise SystemExit("PyMOL not found (see README: conda create -n pymol -c conda-forge pymol-open-source)")
cif = ROOT / "assets" / f"{a.pdb.upper()}.cif"
script = Path(__file__).with_name("pymol_turntable.py")
subprocess.run([exe, "-cq", str(script), "--", str(cif), a.out, str(a.frames), str(a.width), str(a.height),
                str(a.antialias), str(a.surface_quality)], check=True)
