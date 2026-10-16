#!/usr/bin/env python3
"""Fetch a PDB entry and convert it to a per-chain-colored OBJ surface mesh.

Surface geometry backends, tried in order:
  1. PyMOL headless (solvent-excluded surface):
         pymol -cq -d "load 6ql9.cif; show surface; save out.obj"
  2. ChimeraX headless:
         chimerax --nogui --exit --cmd "open 6ql9.cif; surface; save out.obj"
  3. Built-in fallback: Gaussian density grid (gemmi) + marching cubes (scikit-image).

PyMOL's OBJ export carries no colors, so every backend's mesh is welded, decimated and
then colored per chain (nearest atom) here. Output: vertex-colored OBJ.

PyMOL lookup: $PYMOL, `pymol` on PATH, then a conda env named `pymol` under ~/miniforge3, ~/miniconda3, etc.
(install: conda create -n pymol -c conda-forge pymol-open-source).

Usage:
    python scripts/fetch_and_convert.py                       # 6QL9 -> assets/6QL9_surface.obj
    python scripts/fetch_and_convert.py --backend builtin --voxel 4.0
"""
import argparse
import colorsys
import os
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ASSETS = ROOT / "assets"


def fetch(pdb_id: str) -> Path:
    cif = ASSETS / f"{pdb_id}.cif"
    if not cif.exists():
        ASSETS.mkdir(exist_ok=True)
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, cif)
    return cif


def read_atoms(cif: Path):
    """Heavy-atom coordinates (N,3) and the chain index of each atom (N,)."""
    import gemmi

    st = gemmi.read_structure(str(cif))
    st.remove_hydrogens()
    pts, chain_ids = [], []
    for ci, ch in enumerate(st[0]):
        for r in ch:
            if r.is_water():
                continue
            for a in r:
                pts.append((a.pos.x, a.pos.y, a.pos.z))
                chain_ids.append(ci)
    print(f"{len(pts)} atoms, {max(chain_ids) + 1} chains")
    return np.array(pts, dtype=np.float32), np.array(chain_ids)


def find_pymol():
    cands = [os.environ.get("PYMOL"), shutil.which("pymol")]
    for base in ("miniforge3", "mambaforge", "miniconda3", "anaconda3", "micromamba", ".conda"):
        cands.append(str(Path.home() / base / "envs/pymol/bin/pymol"))
    cands.append(str(Path.home() / "micromamba/envs/pymol/bin/pymol"))
    for cand in cands:
        if cand and Path(cand).exists():
            return cand
    return None


def run_pymol(cif: Path, quality: int):
    """Returns a trimesh.Trimesh (raw PyMOL surface) or None."""
    import trimesh

    exe = find_pymol()
    if not exe:
        return None
    with tempfile.TemporaryDirectory() as tmp:
        raw = Path(tmp) / "raw.obj"
        cmd = (f"load {cif}, m; remove solvent or hydrogens; hide everything; "
               f"set surface_quality, {quality}; show surface, m; save {raw}")
        print("Running PyMOL (about a minute)...")
        if subprocess.run([exe, "-cq", "-d", cmd]).returncode != 0 or not raw.exists():
            return None
        return trimesh.load(raw, force="mesh", process=False)


def run_chimerax(cif: Path):
    import trimesh

    exe = shutil.which("chimerax") or shutil.which("ChimeraX")
    if not exe:
        return None
    with tempfile.TemporaryDirectory() as tmp:
        raw = Path(tmp) / "raw.obj"
        cmd = f"open {cif}; surface; save {raw}"
        if subprocess.run([exe, "--nogui", "--exit", "--cmd", cmd]).returncode != 0 or not raw.exists():
            return None
        return trimesh.load(raw, force="mesh", process=False)


def run_builtin(pts, voxel: float, sigma: float, level: float):
    import trimesh
    from scipy.ndimage import gaussian_filter
    from skimage import measure

    pad = 4 * sigma * voxel + 6.0
    lo = pts.min(0) - pad
    shape = np.ceil((pts.max(0) + pad - lo) / voxel).astype(int) + 1
    idx = np.round((pts - lo) / voxel).astype(int)
    grid = np.zeros(shape, dtype=np.float32)
    np.add.at(grid, tuple(idx.T), 1.0)
    grid = gaussian_filter(grid, sigma)
    grid /= grid.max()
    verts, faces, _, _ = measure.marching_cubes(grid, level=level)
    return trimesh.Trimesh(verts * voxel + lo, faces, process=False)


def keep_exterior(mesh, pts, voxel=1.0, probe_reach=3.1):
    """Drop interior faces (chain interfaces, buried cavities) via a solvent flood fill.

    A voxel is blocked if within `probe_reach` Å of an atom (vdW + probe radius); free voxels
    connected to the grid corner are 'outside'. A face is kept if the point just off its
    surface, along the normal, lies in (or next to) an outside voxel.
    """
    from scipy import ndimage

    lo = pts.min(0) - 8.0
    shape = np.ceil((pts.max(0) + 8.0 - lo) / voxel).astype(int) + 1
    occ = np.zeros(shape, dtype=bool)
    occ[tuple(np.round((pts - lo) / voxel).astype(int).T)] = True
    r = int(np.ceil(probe_reach / voxel))
    g = np.mgrid[-r:r + 1, -r:r + 1, -r:r + 1]
    ball = (g ** 2).sum(0) * voxel ** 2 <= probe_reach ** 2
    occ = ndimage.binary_dilation(occ, structure=ball)
    labels, _ = ndimage.label(~occ)
    outside = ndimage.binary_dilation(labels == labels[0, 0, 0], iterations=1)

    centers = mesh.triangles_center + mesh.face_normals * (probe_reach - 1.7)
    idx = np.clip(np.round((centers - lo) / voxel).astype(int), 0, np.array(shape) - 1)
    keep = outside[tuple(idx.T)]
    print(f"Exterior faces: {keep.sum()} / {len(keep)}")
    mesh.update_faces(keep)
    mesh.remove_unreferenced_vertices()
    return mesh


def finalize(mesh, pts, chain_ids, max_faces: int, out: Path, exterior_only: bool = False):
    """Weld, drop interior faces, decimate, color per chain by nearest atom, export."""
    from scipy.spatial import cKDTree

    mesh.merge_vertices()
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    print(f"Welded mesh: {len(mesh.vertices)} verts, {len(mesh.faces)} faces")
    if exterior_only:
        mesh = keep_exterior(mesh, pts)
    if len(mesh.faces) > max_faces:
        mesh = mesh.simplify_quadric_decimation(face_count=max_faces)
        print(f"Decimated to {len(mesh.faces)} faces")

    _, nearest = cKDTree(pts).query(mesh.vertices)
    n_chains = chain_ids.max() + 1
    palette = np.array([colorsys.hsv_to_rgb((i * 0.382) % 1.0, 0.65, 0.85) for i in range(n_chains)])
    rgb = (palette[chain_ids[nearest]] * 255).astype(np.uint8)
    mesh.visual.vertex_colors = np.hstack([rgb, np.full((len(rgb), 1), 255, np.uint8)])
    mesh.export(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", default="6QL9")
    ap.add_argument("--backend", choices=["auto", "pymol", "chimerax", "builtin"], default="auto")
    ap.add_argument("--max-faces", type=int, default=200000)
    ap.add_argument("--exterior-only", action="store_true",
                    help="drop interior faces (helps the Cairo scene; the OpenGL scene doesn't need it)")
    ap.add_argument("--pymol-quality", type=int, default=-1, help="PyMOL surface_quality (-1..4)")
    ap.add_argument("--voxel", type=float, default=3.0, help="voxel size in Å (builtin)")
    ap.add_argument("--sigma", type=float, default=1.5, help="Gaussian sigma in voxels (builtin)")
    ap.add_argument("--level", type=float, default=0.05, help="iso-level, fraction of max density (builtin)")
    args = ap.parse_args()

    pdb = args.pdb.upper()
    cif = fetch(pdb)
    out = ASSETS / f"{pdb}_surface.obj"
    pts, chain_ids = read_atoms(cif)

    runners = {
        "pymol": lambda: run_pymol(cif, args.pymol_quality),
        "chimerax": lambda: run_chimerax(cif),
        "builtin": lambda: run_builtin(pts, args.voxel, args.sigma, args.level),
    }
    order = ["pymol", "chimerax", "builtin"] if args.backend == "auto" else [args.backend]
    for name in order:
        mesh = runners[name]()
        if mesh is not None:
            finalize(mesh, pts, chain_ids, args.max_faces, out, args.exterior_only)
            print(f"[{name}] wrote {out}")
            return
        print(f"[{name}] unavailable or failed, trying next")
    raise SystemExit("No backend produced a mesh")


if __name__ == "__main__":
    main()
