"""Electron-density time series scene: PyMOL panels (whole molecule and/or zoom) + R_iso bar chart.

Normally run via scripts/render_density.py. Manually:
    DENSITY_DIR=assets/density manim -r 1280,720 --fps 30 scenes/density_scene.py DensityScene
Reads $DENSITY_DIR/meta.json and $DENSITY_DIR/frames/{overview,zoom}/frame00001.png ...
"""
import json
import os
from pathlib import Path

import numpy as np
from manim import *
from PIL import Image

DENSITY_DIR = Path(os.environ.get("DENSITY_DIR") or Path(__file__).resolve().parent.parent / "assets" / "density")
SET_SECONDS = float(os.environ.get("DENSITY_SET_SECONDS", 2.5))   # seconds per dataset
FRAMES_PER_SET = int(os.environ.get("DENSITY_FRAMES_PER_SET", 10))
TITLE = os.environ.get("DENSITY_TITLE", "Electron density change across the series")
SUBTITLE = os.environ.get("DENSITY_SUBTITLE", "Density from each dataset's structure factors; difference density vs the first dataset")
ZOOM_LABEL = os.environ.get("DENSITY_ZOOM_LABEL", "Focus residues")
PANELS = os.environ.get("DENSITY_PANELS", "both")   # both | zoom | overview
PANEL_LIST = ["overview", "zoom"] if PANELS == "both" else [PANELS]


def load(panel, n):
    return np.asarray(Image.open(DENSITY_DIR / "frames" / panel / f"frame{n:05d}.png").convert("RGBA"))


class DensityScene(Scene):
    def construct(self):
        meta = json.loads((DENSITY_DIR / "meta.json").read_text())
        sets = meta["datasets"]
        n_sets = len(sets)
        n_frames = n_sets * FRAMES_PER_SET
        if not (DENSITY_DIR / "frames" / PANEL_LIST[0] / f"frame{n_frames:05d}.png").exists():
            raise FileNotFoundError(f"Expected {n_frames} frames in {DENSITY_DIR / 'frames'}; run scripts/render_density.py")
        total = n_sets * SET_SECONDS

        def keep_live(m):  # Cairo caches mobjects without updaters as a static background
            m.add_updater(lambda _m: None)
            return m

        # ---- layout ----
        title = Text(TITLE, font_size=32).to_corner(UL)
        sub = Text(SUBTITLE, font_size=18, color=GREY_A).next_to(title, DOWN, buff=0.12).align_to(title, LEFT)
        captions = {"overview": "Whole molecule: difference density Fo(n) − Fo(1)",
                    "zoom": f"{ZOOM_LABEL}: 2Fo−Fc and difference density"}
        if len(PANEL_LIST) == 2:
            side, centers = 4.2, {"overview": [-3.6, -0.05, 0], "zoom": [3.6, -0.05, 0]}
        else:                                      # a single panel is shown larger, in the middle
            side, centers = 4.7, {PANEL_LIST[0]: [0, 0.0, 0]}
        imgs = {pn: keep_live(ImageMobject(load(pn, 1)).set(height=side).move_to(centers[pn])) for pn in PANEL_LIST}
        caps = [Text(captions[pn], font_size=17, color=GREY_A).next_to(imgs[pn], UP, buff=0.08) for pn in PANEL_LIST]

        # legend
        def swatch(color, text):
            return VGroup(Line(LEFT * 0.18, RIGHT * 0.18, color=color, stroke_width=5),
                          Text(text, font_size=15, color=GREY_A)).arrange(RIGHT, buff=0.12)
        legend = VGroup(
            swatch("#6A6AFF", f"2Fo−Fc at {meta['twofofc_contour']:.1f} σ (set 1 scale)"),
            swatch(GREEN, f"+ Fo(n)−Fo(1), contour +{meta['diff_nsigma']:g} σ of final set"),
            swatch(RED, "− Fo(n)−Fo(1), same contour for every frame")).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        legend.to_corner(DR, buff=0.3)

        # R_iso bar chart (set 1 is the reference -> 0)
        vals = [0.0] + [100 * s["r_iso"] for s in sets[1:]]
        vmax = max(vals) * 1.15
        bw, gap, ch, x0, y0 = 0.5, 0.2, 0.8, -5.5, -3.45
        bars, labels = [], []
        for k, v in enumerate(vals):
            h = max(ch * v / vmax, 0.02)
            b = Rectangle(width=bw, height=h, stroke_width=0, fill_opacity=1, fill_color=GREY_C)
            b.move_to([x0 + k * (bw + gap) + bw / 2, y0 + h / 2, 0])
            bars.append(keep_live(b))
            labels.append(Text(str(k + 1), font_size=15, color=GREY_A).move_to([b.get_x(), y0 - 0.22, 0]))
        axis_cap = Text("Change in diffraction data vs data set 1 (R_iso, %)", font_size=16, color=GREY_A).move_to([x0 + 2.0, y0 + ch + 0.3, 0])
        set_cap = Text("data set", font_size=15, color=GREY_A).move_to([x0 - 0.85, y0 - 0.22, 0])

        # per-set readouts (pre-created, toggled by opacity)
        reads = []
        for s in sets:
            name = s.get("label") or f"Data set {s['set']}"
            txt = name + ("" if s["id"] == name else f"  ·  {s['id']}") + ("  (reference)" if s["set"] == 1 else f"  ·  R_iso {100 * s['r_iso']:.1f} %")
            reads.append(keep_live(Text(txt, font_size=24, color=YELLOW).to_corner(UR).set_opacity(0)))

        state = {"t": 0.0, "set": -1, "key": None}
        DISSOLVE = 0.2     # last 20% of each dataset's time is a smooth dissolve into the next dataset

        def panel_image(panel, cur, nxt, w):
            a = load(panel, cur + 1)
            if w <= 0:
                return a
            b = load(panel, nxt + 1)
            return ((1 - w) * a + w * b).astype(np.uint8)

        def update(dt):
            state["t"] = min(state["t"] + dt, total - 1e-6)
            s_idx = min(int(state["t"] / SET_SECONDS), n_sets - 1)
            u = state["t"] / SET_SECONDS - s_idx                     # 0..1 within this dataset
            cur = s_idx * FRAMES_PER_SET + min(int(u * FRAMES_PER_SET), FRAMES_PER_SET - 1)
            w = 0.0
            if s_idx < n_sets - 1 and u > 1 - DISSOLVE:
                x = (u - (1 - DISSOLVE)) / DISSOLVE
                w = x * x * (3 - 2 * x)                                # smoothstep
            key = (cur, round(w, 2))
            if key != state["key"]:
                state["key"] = key
                for pn in PANEL_LIST:
                    imgs[pn].pixel_array = panel_image(pn, cur, (s_idx + 1) * FRAMES_PER_SET, w)
            if s_idx != state["set"]:
                if state["set"] >= 0:
                    reads[state["set"]].set_opacity(0)
                    bars[state["set"]].set_fill(GREY_C)
                reads[s_idx].set_opacity(1)
                bars[s_idx].set_fill(YELLOW)
                state["set"] = s_idx

        self.add(*imgs.values(), *caps, legend, *bars, *labels, axis_cap, set_cap, *reads)
        self.play(FadeIn(title), FadeIn(sub), run_time=1.2)
        self.add_updater(update)
        self.wait(total)
        self.remove_updater(update)
        self.play(FadeOut(title), FadeOut(sub), run_time=1.0)
