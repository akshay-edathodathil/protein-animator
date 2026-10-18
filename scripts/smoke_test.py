#!/usr/bin/env python3
"""Self-test: draw one frame containing an image (what every PyMOL-composited video does) with Manim.

Catches dependency drift at setup time, e.g. Pillow 13 breaking Manim 0.22's image drawing
("ValueError: Invalid mode for data type"). Run by ./setup.sh; can be run on its own.
"""
import sys
import tempfile

import numpy as np


def main():
    from manim import ImageMobject, Scene, config

    with tempfile.TemporaryDirectory() as tmp:
        config.media_dir = tmp
        config.save_last_frame = True
        config.quality = "low_quality"
        config.verbosity = "ERROR"
        config.progress_bar = "none"

        class ImageSmoke(Scene):
            def construct(self):
                self.add(ImageMobject(np.full((16, 16, 4), 200, dtype=np.uint8)))

        try:
            ImageSmoke().render()
        except Exception as e:  # noqa: BLE001
            import PIL
            print(f"SMOKE TEST FAILED: {type(e).__name__}: {e}\n"
                  f"Pillow {PIL.__version__}. Known fix: pip install 'Pillow>=11,<13' (Pillow 13 breaks Manim 0.22).")
            return 1
    print("Smoke test OK (Manim can draw images).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
