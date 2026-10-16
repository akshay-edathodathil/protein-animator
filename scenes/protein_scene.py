"""3D surface animation of PDB 6QL9 (fatty acid synthase) with Manim (Cairo renderer)."""
from pathlib import Path

import numpy as np
import trimesh
from manim import *

OBJ_PATH = Path(__file__).resolve().parent.parent / "assets" / "6QL9_surface.obj"
TARGET_FACES = 15000  # Cairo depth-sorts every face each frame; raise for quality, lower for speed
TARGET_WIDTH = 4.0    # longest mesh extent in Manim units
LIGHT_DIR = np.array([-0.4, -0.5, 0.75])


def load_mesh(path=OBJ_PATH, target_faces=TARGET_FACES, width=TARGET_WIDTH):
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run: python scripts/fetch_and_convert.py")
    mesh = trimesh.load(path, force="mesh")
    src_colors = None
    if mesh.visual.kind == "vertex":
        src_colors = (np.asarray(mesh.visual.vertex_colors)[:, :3] / 255.0, np.asarray(mesh.vertices))
    if len(mesh.faces) > target_faces:
        mesh = mesh.simplify_quadric_decimation(face_count=target_faces)
    colors = None
    if src_colors is not None:  # decimation drops colors: copy from nearest original vertex
        from scipy.spatial import cKDTree
        _, nn = cKDTree(src_colors[1]).query(np.asarray(mesh.vertices))
        colors = src_colors[0][nn]
    verts = np.asarray(mesh.vertices, dtype=float)
    verts -= (verts.max(0) + verts.min(0)) / 2          # center at origin
    verts *= width / np.ptp(verts, axis=0).max()          # fit to ~width units
    mesh = trimesh.Trimesh(verts, mesh.faces, process=False)
    mesh.vertex_chain_colors = colors  # None -> fall back to Z gradient
    return mesh


def mesh_to_mobject(mesh, c_low=BLUE_E, c_high=ORANGE):
    """One flat-shaded Polygon per face, colored by chain (or Z gradient) with Lambert-style shade."""
    verts, faces = mesh.vertices, mesh.faces
    normals = mesh.face_normals
    light = LIGHT_DIR / np.linalg.norm(LIGHT_DIR)
    zc = verts[faces][:, :, 2].mean(axis=1)
    t = (zc - zc.min()) / np.ptp(zc)
    vcol = getattr(mesh, "vertex_chain_colors", None)
    lo, hi = ManimColor(c_low), ManimColor(c_high)
    polys = []
    for f, n, ti in zip(faces, normals, t):
        shade = 0.55 + 0.45 * abs(float(n @ light))
        if vcol is not None:  # per-chain color from the OBJ vertex colors
            col = ManimColor(vcol[f].mean(axis=0))
        else:
            col = interpolate_color(lo, hi, float(ti))
        col = interpolate_color(BLACK, col, shade)
        p = Polygon(*verts[f], stroke_width=0)
        p.set_fill(col, opacity=1)
        p.set_stroke(col, width=0.3, opacity=1)  # hides hairline seams between faces
        polys.append(p)
    return VGroup(*polys)


class ProteinVolumeScene(ThreeDScene):
    def construct(self):
        protein = mesh_to_mobject(load_mesh())
        self.set_camera_orientation(phi=70 * DEGREES, theta=-45 * DEGREES)

        title = Text("Fatty Acid Synthase  ·  PDB 6QL9", font_size=34).to_corner(UL)
        self.add_fixed_in_frame_mobjects(title)
        self.play(FadeIn(title), FadeIn(protein, scale=0.6), run_time=2)

        # Slow ambient spin
        self.begin_ambient_camera_rotation(rate=0.25)

        notes = [
            "Megasynthase that builds palmitate (C16) from acetyl- and malonyl-CoA",
            "Multi-domain complex: KS, KR, DH, ER, MAT/PAT, ACP",
            "The mobile acyl carrier protein (ACP) shuttles intermediates between active sites",
        ]
        note = None
        for i, text in enumerate(notes):
            new = Text(text, font_size=24).to_edge(DOWN)
            self.add_fixed_in_frame_mobjects(new)
            if note is None:
                self.play(FadeIn(new))
            else:
                self.play(FadeOut(note), FadeIn(new))
            note = new
            if i == 1:  # zoom in on the structure
                self.move_camera(zoom=1.6, run_time=3)
            elif i == 2:
                self.move_camera(zoom=1.0, run_time=3)
            else:
                self.wait(3)
            self.wait(1)

        self.stop_ambient_camera_rotation()
        self.move_camera(phi=20 * DEGREES, theta=-90 * DEGREES, run_time=3)
        self.play(FadeOut(note), FadeOut(title))
        self.wait(1)
