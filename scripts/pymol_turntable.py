"""Runs INSIDE PyMOL:  pymol -cq scripts/pymol_turntable.py -- <cif> <out_dir> <frames> <width> <height> [antialias] [surface_quality]

Ray-traces a 360° turntable of a chain-colored molecular surface with a transparent background.
Use scripts/render_pymol_frames.py as the launcher.
"""
import os
import sys
import time

from pymol import cmd

cif, out_dir, n_frames, width, height = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
antialias = int(sys.argv[6]) if len(sys.argv) > 6 else 2
quality = int(sys.argv[7]) if len(sys.argv) > 7 else 0
os.makedirs(out_dir, exist_ok=True)

cmd.load(cif, "m")
cmd.remove("solvent or hydrogens")
cmd.hide("everything")
cmd.set("surface_quality", quality)
cmd.show("surface", "m")
cmd.util.cbc()  # one color per chain

cmd.orient("m")
cmd.zoom("m", buffer=12, complete=1)

cmd.bg_color("black")
cmd.set("ray_opaque_background", "off")   # transparent PNG background
cmd.set("antialias", antialias)
cmd.set("ray_shadow", "on")
cmd.set("ambient", 0.4)
cmd.set("spec_reflect", 0.3)
cmd.set("ray_trace_mode", 0)

step = 360.0 / n_frames
t0 = time.time()
for i in range(n_frames):
    path = os.path.join(out_dir, f"protein_frame{i + 1:04d}.png")
    if not os.path.exists(path):  # resumable
        cmd.png(path, width, height, dpi=150, ray=1)
        done = i + 1
        print(f"frame {done}/{n_frames}  elapsed {time.time() - t0:.0f}s", flush=True)
    cmd.turn("y", step)
