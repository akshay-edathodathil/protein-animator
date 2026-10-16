"""Trajectory scene: PyMOL ray-traced frames of a changing structure + Manim title, timeline, RMSD, legend.

Normally run via scripts/render_traj.py. Manually:
    TRAJ_DIR=assets/trajectory manim -r 1280,720 --fps 30 scenes/trajectory_scene.py TrajectoryScene
Reads $TRAJ_DIR/meta.json and $TRAJ_DIR/frames/frame00001.png ...
"""
import bisect
import json
import os
from pathlib import Path

import numpy as np
from manim import *
from PIL import Image

TRAJ_DIR = Path(os.environ.get("TRAJ_DIR") or Path(__file__).resolve().parent.parent / "assets" / "trajectory")
PLAY_FPS = float(os.environ.get("TRAJ_PLAY_FPS", 15))   # trajectory frames per second of video between keyframes
HOLD = float(os.environ.get("TRAJ_HOLD", 1.2))          # seconds to pause on each real (keyframe) structure


class TrajectoryScene(Scene):
    def construct(self):
        meta = json.loads((TRAJ_DIR / "meta.json").read_text())
        frames = meta["frames"]
        n_keys = len(meta["keyframes"])
        frame_files = [TRAJ_DIR / "frames" / f"frame{f['state']:05d}.png" for f in frames]
        if not frame_files[0].exists():
            raise FileNotFoundError(f"No frames in {TRAJ_DIR / 'frames'}. Run scripts/render_traj.py first.")

        # Playback schedule: dwell on keyframes, run through interpolated frames at PLAY_FPS.
        starts, t = [], 0.0
        for f in frames:
            starts.append(t)
            t += HOLD if f["key"] else 1.0 / PLAY_FPS
        total = t

        # ---- static layout ----
        title = Text(meta["title"], font_size=34).to_corner(UL)
        img = ImageMobject(np.asarray(Image.open(frame_files[0]).convert("RGBA")))
        img.set(height=config.frame_height * 0.72).shift(UP * 0.35)

        left, right, y = -4.5, 3.5, -3.15
        line = Line([left, y, 0], [right, y, 0], color=GREY_B, stroke_width=3)
        tx = lambda tt: left + (right - left) * tt / max(n_keys - 1, 1)
        ticks = VGroup(*[Line([tx(i), y - 0.1, 0], [tx(i), y + 0.1, 0], color=WHITE, stroke_width=3) for i in range(n_keys)])
        show_all = n_keys <= 6
        tick_labels = VGroup(*[
            Text(lab, font_size=18, color=GREY_A).next_to([tx(i), y, 0], DOWN, buff=0.18)
            for i, lab in enumerate(meta["keyframes"]) if show_all or i in (0, n_keys - 1)])
        marker = Dot([left, y, 0], color=YELLOW, radius=0.09)
        static = VGroup(line, ticks, tick_labels)

        # Legend (matches PyMOL's blue_white_red spectrum used for displacement coloring)
        legend = VGroup()
        if meta.get("color_by", "displacement") == "displacement":
            n = 40
            bar = VGroup(*[Rectangle(width=3.2 / n, height=0.16, stroke_width=0, fill_opacity=1,
                                     fill_color=interpolate_color(ManimColor('#0000FF'), WHITE, 2 * k / (n - 1)) if k < n / 2
                                     else interpolate_color(WHITE, ManimColor('#FF0000'), 2 * (k - n / 2) / (n - 1)))
                           for k in range(n)]).arrange(RIGHT, buff=0)
            cap = Text("Distance moved from first structure", font_size=16, color=GREY_A)
            lo = Text("0 Å", font_size=16, color=GREY_A)
            hi = Text(f"{meta['disp_max']:.1f} Å", font_size=16, color=GREY_A)
            legend = VGroup(cap, bar, VGroup(lo, hi))
            lo.next_to(bar, DOWN, buff=0.08).align_to(bar, LEFT)
            hi.next_to(bar, DOWN, buff=0.08).align_to(bar, RIGHT)
            cap.next_to(bar, UP, buff=0.08)
            legend.to_corner(UR).shift(DOWN * 1.15)

        # NOTE: Manim's Cairo renderer caches mobjects WITHOUT updaters as a static background, so
        # anything that changes during playback must carry an updater. Labels/readouts are therefore
        # all created up front, kept in the scene, and shown/hidden by opacity.
        def keep_live(m):
            m.add_updater(lambda _m: None)
            return m

        labels = {text: keep_live(Text(text, font_size=26, color=YELLOW).to_corner(UR).set_opacity(0))
                  for text in dict.fromkeys(f["label"] for f in frames)}
        rmsd_keys = sorted({round(f["rmsd"], 1) for f in frames})
        rmsds = {k: keep_live(Text(f"CA RMSD vs first: {k:.1f} Å", font_size=20, color=GREY_A)
                              .to_corner(UR).shift(DOWN * 0.55).set_opacity(0)) for k in rmsd_keys}
        keep_live(img)
        keep_live(marker)

        state = {"t": 0.0, "idx": -1, "label": None, "rmsd": None}

        def show(old, new):
            if old is not None and old is not new:
                old.set_opacity(0)
            new.set_opacity(1)
            return new

        def update(dt):
            state["t"] = min(state["t"] + dt, total - 1e-6)
            idx = bisect.bisect_right(starts, state["t"]) - 1
            if idx != state["idx"]:
                state["idx"] = idx
                img.pixel_array = np.asarray(Image.open(frame_files[idx]).convert("RGBA"))
                f = frames[idx]
                state["label"] = show(state["label"], labels[f["label"]])
                state["rmsd"] = show(state["rmsd"], rmsds[round(f["rmsd"], 1)])
                marker.move_to([tx(f["t"]), y, 0])

        self.add(img, static, marker, legend, *labels.values(), *rmsds.values())
        self.play(FadeIn(title), run_time=1.2)
        self.add_updater(update)
        self.wait(total)
        self.remove_updater(update)
        self.play(FadeOut(title), run_time=1.0)
