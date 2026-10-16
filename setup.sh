#!/usr/bin/env bash
# One-shot setup: venv + Python deps (+ optional PyMOL via conda).
#   ./setup.sh              # core install (enough for the OpenGL / Cairo scenes)
#   ./setup.sh --with-pymol # also create a conda env "pymol" (needed for the hybrid scene)
set -euo pipefail
cd "$(dirname "$0")"

PY="${PYTHON:-python3}"
if ! command -v "$PY" >/dev/null; then echo "python3 not found. Install Python 3.10+."; exit 1; fi
"$PY" - <<'PYEOF' || { echo "Python 3.10 or newer is required."; exit 1; }
import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)
PYEOF

missing=()
command -v ffmpeg >/dev/null || missing+=(ffmpeg)
if [[ "$(uname)" == "Darwin" ]]; then
  command -v brew >/dev/null && for p in cairo pango pkg-config; do brew list "$p" >/dev/null 2>&1 || missing+=("$p"); done
  [[ ${#missing[@]} -gt 0 ]] && echo "Missing system packages: ${missing[*]}  ->  brew install ${missing[*]}"
else
  [[ ${#missing[@]} -gt 0 ]] && echo "Missing: ${missing[*]}"
  echo "Linux: sudo apt install ffmpeg libcairo2-dev libpango1.0-dev pkg-config python3-dev  (or your distro's equivalent)"
fi
if ! command -v ffmpeg >/dev/null; then echo "ffmpeg is required. Install it and re-run."; exit 1; fi

"$PY" -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt

if [[ "${1:-}" == "--with-pymol" ]]; then
  if command -v mamba >/dev/null; then C=mamba; elif command -v conda >/dev/null; then C=conda; else
    echo "conda/mamba not found. Install Miniforge (https://github.com/conda-forge/miniforge) then re-run with --with-pymol."; exit 1; fi
  "$C" create -y -n pymol -c conda-forge pymol-open-source
fi

.venv/bin/python -c "import manim, trimesh, gemmi; print('Python deps OK, manim', manim.__version__)"
echo
echo "Setup done. Next:  ./render.sh gl      (OpenGL scene, ~seconds)"
echo "Hybrid (PyMOL):    ./setup.sh --with-pymol  &&  ./render.sh frames-fast  &&  ./render.sh hybrid"
