"""OpenGL-renderer version of the 6QL9 scene: full-resolution mesh, real depth buffer.

Run:  manim -pql --renderer=opengl scenes/protein_scene_gl.py ProteinVolumeSceneGL
"""
from pathlib import Path

import numpy as np
import trimesh
from manim import *
from manim.mobject.opengl.opengl_surface import OpenGLSurface

OBJ_PATH = Path(__file__).resolve().parent.parent / "assets" / "6QL9_surface.obj"
TARGET_WIDTH = 4.0


class MeshSurface(OpenGLSurface):
    """Triangle mesh drawn with Manim's OpenGL surface shader (per-vertex color + smooth normals)."""

    def __init__(self, vertices, faces, vertex_colors, vertex_normals, **kwargs):
        self._verts = np.asarray(vertices, dtype=np.float32)
        self._faces = np.asarray(faces, dtype=np.int64)
        self._vcolors = np.asarray(vertex_colors, dtype=np.float32)
        self._vnormals = np.asarray(vertex_normals, dtype=np.float32)
        super().__init__(**kwargs)

    def init_points(self):
        # The surface shader recovers normals from two points nudged along tangent directions.
        n = self._vnormals
        helper = np.where(np.abs(n[:, :1]) < 0.9, [[1, 0, 0]], [[0, 1, 0]])
        t1 = np.cross(n, helper)
        t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
        t2 = np.cross(n, t1)
        eps = self.epsilon
        self.set_points(np.vstack([self._verts, self._verts + eps * t1, self._verts + eps * t2]))
        return self

    def compute_triangle_indices(self):
        self.triangle_indices = self._faces.reshape(-1).astype(np.uint32)
        return self

    def get_shader_data(self):
        s, du, dv = self.get_surface_points_and_nudged_points()
        data = np.zeros(len(s), dtype=self.shader_dtype)
        data["point"], data["du_point"], data["dv_point"] = s, du, dv
        data["color"] = self._vcolors
        return data


def load_mesh_mobject(path=OBJ_PATH, width=TARGET_WIDTH):
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run: python scripts/fetch_and_convert.py")
    mesh = trimesh.load(path, force="mesh")
    v = np.asarray(mesh.vertices, dtype=float)
    v -= (v.max(0) + v.min(0)) / 2
    v *= width / np.ptp(v, axis=0).max()
    mesh = trimesh.Trimesh(v, mesh.faces, vertex_colors=mesh.visual.vertex_colors, process=False)
    colors = np.asarray(mesh.visual.vertex_colors, dtype=np.float32) / 255.0
    # epsilon must stay small relative to the scene scale for stable normals
    return MeshSurface(mesh.vertices, mesh.faces, colors, mesh.vertex_normals,
                       gloss=0.2, shadow=0.3, epsilon=1e-4)


class ProteinVolumeSceneGL(ThreeDScene):
    def construct(self):
        protein = load_mesh_mobject()
        self.set_camera_orientation(phi=70 * DEGREES, theta=-45 * DEGREES)

        title = Text("Fatty Acid Synthase  ·  PDB 6QL9", font_size=34).to_corner(UL)
        title.fix_in_frame()
        self.add(protein)
        self.play(FadeIn(title), run_time=1.5)

        notes = [
            "Megasynthase that builds palmitate (C16) from acetyl- and malonyl-CoA",
            "Hetero 12-mer (α6β6, D3 symmetry): each chain shown in its own color",
            "Domains KS, KR, DH, ER, MAT/PAT and the mobile acyl carrier protein (ACP)",
        ]
        note = None
        for i, text in enumerate(notes):
            new = Text(text, font_size=24).to_edge(DOWN)
            new.fix_in_frame()
            if note is None:
                self.play(FadeIn(new))
            else:
                self.play(FadeOut(note), FadeIn(new))
                self.remove(note)
            note = new
            if i == 1:
                self.play(self.camera.animate.scale(0.6).increment_theta(60 * DEGREES), run_time=4)
            elif i == 2:
                self.play(self.camera.animate.scale(1 / 0.6).increment_theta(60 * DEGREES), run_time=4)
            else:
                self.play(self.camera.animate.increment_theta(60 * DEGREES), run_time=4)

        self.play(self.camera.animate.set_phi(20 * DEGREES), FadeOut(note), FadeOut(title), run_time=3)
        self.remove(note, title)
        self.wait(1)
