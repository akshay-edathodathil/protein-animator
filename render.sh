#!/usr/bin/env bash
# Render helper. Extra arguments are passed to manim (e.g. -qh, -qk, -p to preview).
#   ./render.sh mesh            # (re)generate assets/6QL9_surface.obj (PyMOL if available, else built-in)
#   ./render.sh gl [flags]      # OpenGL scene, full-res mesh          (default 720p30)
#   ./render.sh cairo [flags]   # Cairo scene, decimated mesh           (default -ql)
#   ./render.sh frames-fast     # PyMOL ray-traced turntable, 60 frames, fast settings
#   ./render.sh frames-hq       # PyMOL ray-traced turntable, 60 frames, high quality (slow)
#   ./render.sh hybrid [flags]  # PyMOL frames + Manim titles           (default 720p30)
#   ./render.sh traj <pdb files...> [options]   # structure CHANGE over time from several PDBs (see README)
#   ./render.sh traj-demo       # adenylate kinase open -> closed example
set -euo pipefail
cd "$(dirname "$0")"
[[ -x .venv/bin/python ]] || { echo "Run ./setup.sh first."; exit 1; }
source .venv/bin/activate
mode="${1:-}"; shift || true
default_flags=(-r 1280,720 --fps 30)

case "$mode" in
  mesh)        python scripts/fetch_and_convert.py "$@" ;;
  gl)          manim "${@:-${default_flags[@]}}" --renderer=opengl scenes/protein_scene_gl.py ProteinVolumeSceneGL ;;
  cairo)       manim "${@:--ql}" scenes/protein_scene.py ProteinVolumeScene ;;
  frames-fast) python scripts/render_pymol_frames.py --frames 60 --width 1280 --height 720 --antialias 1 --surface-quality -1 --out assets/pymol_frames ;;
  frames-hq)   python scripts/render_pymol_frames.py --frames 60 --width 1280 --height 720 --antialias 2 --surface-quality 0 --out assets/pymol_frames_hq ;;
  hybrid)
    if [[ -z "${FRAME_DIR:-}" ]]; then
      if [[ -d assets/pymol_frames_hq && -n "$(ls assets/pymol_frames_hq 2>/dev/null)" ]]; then FRAME_DIR=assets/pymol_frames_hq
      else FRAME_DIR=assets/pymol_frames; fi
    fi
    echo "Using frames from $FRAME_DIR"
    [[ -n "$(ls "$FRAME_DIR" 2>/dev/null)" ]] || { echo "No frames. Run ./render.sh frames-fast first."; exit 1; }
    FRAME_DIR="$FRAME_DIR" manim "${@:-${default_flags[@]}}" scenes/protein_scene_hybrid.py ProteinHybridScene ;;
  traj)        python scripts/render_traj.py "$@" ;;
  traj-demo)   python scripts/render_traj.py examples/adenylate_kinase/4AKE.pdb examples/adenylate_kinase/1AKE.pdb \
                 --labels "Open (4AKE),Closed (1AKE)" --chains A --title "Adenylate kinase: lid closing" "$@" ;;
  *) sed -n '2,12p' "$0"; exit 1 ;;
esac
