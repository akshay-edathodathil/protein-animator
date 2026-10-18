"""Runs INSIDE PyMOL:  pymol -cq scripts/pymol_density.py -- <config.json>

Per dataset, ray-traces two square panels (or only one, see cfg['panels']) from fixed cameras (+ slow rotation):
  overview/  whole molecule (cartoon) with difference density  Fo(n)-Fo(1)  (green +, red -)
  zoom/      focus residues as sticks with 2Fo-Fc density (blue) and the difference density
Frames are ordered dataset by dataset, `frames_per_set` frames each. Resumable. Launched by render_density.py.
"""
import json
import os
import sys
import time

from pymol import cmd

cfg = json.load(open(sys.argv[-1]))
meta = json.load(open(os.path.join(cfg["density_dir"], "meta.json")))
maps = os.path.join(cfg["density_dir"], "maps")
ids = meta["ids"]
nm = {i: f"s{k}" for k, i in enumerate(ids)}   # PyMOL-safe object names (dataset names may contain dots, dashes, "and", ...)
F = cfg["frames_per_set"]
size = cfg["size"]
focus_sel = " or ".join(f"(chain {r['chain']} and resi {r['resi']})" for r in meta["focus_residues"])
tight = cfg.get("tight", False)               # single-residue mode: tight view, context sticks, tighter mesh
lvl_e, lvl_d = meta["twofofc_contour"], meta["diff_contour"]

# Keep the absolute scale written by density_prepare.py (PyMOL would otherwise rescale each map to its own sigma,
# making contours incomparable between datasets).
cmd.set("normalize_ccp4_maps", 0)
for i in ids:
    cmd.load(os.path.join(cfg["data_dir"], meta["models"][i]), f"m_{nm[i]}")
    cmd.load(os.path.join(maps, f"{i}_2fofc.ccp4"), f"e_{nm[i]}", format="ccp4")
    if i != ids[0]:
        cmd.load(os.path.join(maps, f"{i}_diff.ccp4"), f"d_{nm[i]}", format="ccp4")
cmd.remove("solvent or hydrogens")
cmd.hide("everything")
cmd.bg_color("black")
cmd.set("ray_opaque_background", "off")
cmd.set("antialias", cfg["antialias"])
cmd.set("ray_shadow", "off")
cmd.set("ray_trace_mode", cfg["ray_mode"])
cmd.set("mesh_width", 1.2)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("stick_radius", 0.16)

first = ids[0]
# Views are set up on the first model; later models/maps share (nearly) the same crystal frame.
cmd.disable("all")
cmd.enable(f"m_{nm[first]}")
cmd.show("cartoon", f"m_{nm[first]}")
cmd.orient(f"m_{nm[first]}")
cmd.zoom(f"m_{nm[first]}", buffer=2)
view_over = cmd.get_view()
cmd.show("sticks", f"m_{nm[first]} and ({focus_sel})")
if tight:                                      # look at the residue along its own principal axes
    cmd.orient(f"m_{nm[first]} and ({focus_sel})")
cmd.zoom(f"m_{nm[first]} and ({focus_sel})", buffer=cfg.get("focus_buffer", 9))
view_zoom = cmd.get_view()
if cfg.get("view"):                           # user-supplied camera (PyMOL get_view) replaces the automatic one
    view_zoom = tuple(cfg["view"])
cmd.hide("everything", f"m_{nm[first]}")
cmd.disable("all")

n_frames = len(ids) * F
step = cfg["rotate"] / max(n_frames - 1, 1)
t0 = time.time()


def meshes(i, sele, carve, with_e, carve_d=None):
    carve_d = carve if carve_d is None else carve_d
    if with_e:
        cmd.isomesh("mesh_e", f"e_{nm[i]}", lvl_e, sele, buffer=max(2, carve), carve=carve)
        cmd.color("slate", "mesh_e")
        cmd.enable("mesh_e")
    if i != first:
        cmd.isomesh("mesh_dp", f"d_{nm[i]}", lvl_d, sele, buffer=max(2, carve_d), carve=carve_d)
        cmd.isomesh("mesh_dn", f"d_{nm[i]}", -lvl_d, sele, buffer=max(2, carve_d), carve=carve_d)
        cmd.color("green", "mesh_dp")
        cmd.color("red", "mesh_dn")
        cmd.enable("mesh_dp")
        cmd.enable("mesh_dn")


for j in range(n_frames):
    si = j // F
    i = ids[si]
    angle = j * step
    for panel in cfg.get("panels", ["overview", "zoom"]):
        os.makedirs(os.path.join(cfg["frames_dir"], panel), exist_ok=True)
        path = os.path.join(cfg["frames_dir"], panel, f"frame{j + 1:05d}.png")
        if os.path.exists(path):
            continue
        cmd.disable("all")
        cmd.delete("mesh_*")
        cmd.hide("everything")
        cmd.enable(f"m_{nm[i]}")
        cmd.show("cartoon", f"m_{nm[i]}")
        cmd.color("grey70", f"m_{nm[i]}")
        cmd.set("cartoon_transparency", (0.8 if tight else 0.65) if panel == "zoom" else 0.0, f"m_{nm[i]}")
        if panel == "overview":
            cmd.set_view(view_over)
            meshes(i, f"m_{nm[i]}", 2.6, False)
        else:
            cmd.show("sticks", f"m_{nm[i]} and ({focus_sel})")
            cmd.util.cnc(f"m_{nm[i]}")
            cmd.color("wheat", f"m_{nm[i]} and elem C and ({focus_sel})")
            if tight:                              # neighbouring residues as thin sticks for context
                nb = f"byres ((m_{nm[i]} and not ({focus_sel})) within 4.5 of (m_{nm[i]} and ({focus_sel}))) and not ({focus_sel})"
                cmd.show("sticks", nb)
                cmd.set("stick_radius", 0.07, nb)
                cmd.color("grey60", f"({nb}) and elem C")
                cmd.set("stick_radius", 0.2, f"m_{nm[i]} and ({focus_sel})")
            cmd.set_view(view_zoom)
            meshes(i, f"m_{nm[i]} and ({focus_sel})", *((1.8, True, 3.0) if tight else (2.2, True)))
        cmd.turn("y", angle)
        cmd.png(path, size, size, dpi=150, ray=1)
    print(f"frame {j + 1}/{n_frames}  elapsed {time.time() - t0:.0f}s", flush=True)
