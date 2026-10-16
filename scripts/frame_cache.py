"""Guard for resumable PyMOL frame renders.

The PyMOL scripts skip frames that already exist, which is only safe if those frames were made with the
same settings. `ensure_fresh` stores the settings next to the frames and wipes them when they change.
"""
import json
import shutil
from pathlib import Path


def ensure_fresh(frames_dir, settings: dict):
    frames_dir = Path(frames_dir)
    stamp = frames_dir / ".settings.json"
    new = json.dumps(settings, sort_keys=True, default=str)
    if frames_dir.exists():
        old = stamp.read_text() if stamp.exists() else None
        has_frames = any(frames_dir.rglob("*.png"))
        if has_frames and old != new:
            reason = "settings changed" if old else "no settings record"
            print(f"Existing frames in {frames_dir} are stale ({reason}); deleting and re-rendering.")
            shutil.rmtree(frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)
    stamp.write_text(new)
