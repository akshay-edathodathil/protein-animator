"""Runs INSIDE PyMOL:  pymol -cq scripts/pymol_traj.py -- <config.json>

Renders one ray-traced transparent PNG per trajectory state from a FIXED camera (plus an optional
slow rotation), coloring each residue by how far it has moved from the first structure.
Launched by scripts/render_traj.py. Resumable: existing frames are skipped.
"""
import array
import json
import os
import sys
import time

from pymol import cmd

cfg = json.load(open(sys.argv[-1]))
out_dir = cfg["frames_dir"]
os.makedirs(out_dir, exist_ok=True)
meta = json.load(open(os.path.join(cfg["traj_dir"], "meta.json")))
n_atoms, n_frames = meta["n_atoms"], meta["n_frames"]
disp = array.array("f")
with open(os.path.join(cfg["traj_dir"], "disp.f32"), "rb") as fh:
    disp.fromfile(fh, n_atoms * n_frames)

cmd.load(os.path.join(cfg["traj_dir"], "morph.pdb"), "m")  # one state per frame
cmd.dss("m", state=1)   # morph.pdb has no HELIX/SHEET records: assign secondary structure
cmd.hide("everything")
cmd.show("cartoon", "m")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)
cmd.set("cartoon_loop_radius", 0.25)
if cfg["color_by"] == "chain":
    cmd.util.cbc("m")
elif cfg["color_by"] == "none":
    cmd.color("marine", "m")

# Fixed camera that frames ALL states, so nothing leaves the view while the protein moves.
cmd.orient("m", state=0)
cmd.zoom("m", buffer=0, state=0)

cmd.bg_color("black")
cmd.set("ray_opaque_background", "off")  # transparent PNG
cmd.set("antialias", cfg["antialias"])
cmd.set("ray_shadow", "on")
cmd.set("ambient", 0.45)
cmd.set("ray_trace_mode", cfg["ray_mode"])  # 1 = black outlines (illustrative cartoon look)

vmax = cfg["disp_max"] or meta["disp_max"]
step = cfg["rotate"] / max(n_frames - 1, 1)  # total slow rotation spread over the movie
t0 = time.time()
for i in range(n_frames):
    state = i + 1
    path = os.path.join(out_dir, f"frame{state:05d}.png")
    if not os.path.exists(path):
        cmd.set("state", state)
        if cfg["ss"] == "each":
            cmd.dss("m", state=state)
        if cfg["color_by"] == "displacement":
            vals = iter(disp[i * n_atoms:(i + 1) * n_atoms])
            cmd.alter("m", "b=next(it)", space={"it": vals})
            cmd.spectrum("b", "blue_white_red", "m", minimum=0, maximum=vmax)
        cmd.png(path, cfg["width"], cfg["height"], dpi=150, ray=1)
        print(f"frame {state}/{n_frames}  elapsed {time.time() - t0:.0f}s", flush=True)
    cmd.turn("y", step)
