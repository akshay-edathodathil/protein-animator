"""Hybrid scene: PyMOL ray-traced turntable frames composited with Manim titles/callouts.

Prereq:  python scripts/render_pymol_frames.py
Run:     manim -pqh scenes/protein_scene_hybrid.py ProteinHybridScene
"""
from pathlib import Path

import numpy as np
from manim import *
from PIL import Image

import os

FRAME_DIR = Path(os.environ.get("FRAME_DIR") or Path(__file__).resolve().parent.parent / "assets" / "pymol_frames")
REV_SECONDS = 30.0   # time for one full 360° turn
TOTAL_SECONDS = 20.0


def load_frames():
    files = sorted(FRAME_DIR.glob("protein_frame*.png"))
    if not files:
        raise FileNotFoundError(f"No frames in {FRAME_DIR}. Run: python scripts/render_pymol_frames.py")
    return np.stack([np.asarray(Image.open(f).convert("RGBA")) for f in files])


class ProteinHybridScene(Scene):
    def construct(self):
        frames = load_frames()
        n = len(frames)

        img = ImageMobject(frames[0]).set(height=config.frame_height * 0.78).shift(UP * 0.3)  # leave room for caption
        clock = {"t": 0.0}

        def spin(mob, dt):
            clock["t"] += dt
            pos = (clock["t"] / REV_SECONDS * n) % n          # fractional frame index
            i0 = int(pos)
            w = pos - i0
            a, b = frames[i0], frames[(i0 + 1) % n]
            # cross-fade neighboring frames so a low-fps turntable still looks smooth
            mob.pixel_array = ((1 - w) * a + w * b).astype(np.uint8)

        img.add_updater(spin)
        self.add(img)

        title = Text("Fatty Acid Synthase  ·  PDB 6QL9", font_size=34).to_corner(UL)
        self.play(FadeIn(title), run_time=1.5)

        notes = [
            "Megasynthase that builds palmitate (C16) from acetyl- and malonyl-CoA",
            "Hetero 12-mer (α6β6, D3 symmetry): each chain shown in its own color",
            "Domains KS, KR, DH, ER, MAT/PAT and the mobile acyl carrier protein (ACP)",
        ]
        note = None
        for i, text in enumerate(notes):
            new = Text(text, font_size=24).to_edge(DOWN)
            if note is None:
                self.play(FadeIn(new))
            else:
                self.play(FadeOut(note), FadeIn(new))
            note = new
            if i == 1:
                self.play(img.animate.scale(1.15), run_time=4)
            elif i == 2:
                self.play(img.animate.scale(1 / 1.15), run_time=4)
            else:
                self.wait(4)

        self.play(FadeOut(note), FadeOut(title), run_time=3)
        self.wait(1)
