# 6QL9 Fatty Acid Synthase — 3D protein animation with Manim

Animated 3D views of PDB **6QL9** (yeast fatty acid synthase, hetero 12-mer α6β6), each of the 12 chains in its own color, with titles and captions. Two ways to render it:

| | **A. OpenGL scene** | **B. Hybrid scene (PyMOL + Manim)** |
|---|---|---|
| How | Manim draws a 200k-face surface mesh itself on the GPU | PyMOL ray-traces a turntable of PNG frames; Manim adds titles/captions on top |
| Look | Flat shading, slightly noisy surface | Smooth, shadowed, publication-style (close to the RCSB image) |
| Camera | Free 3D camera (rotate, zoom, anything) | Fixed turntable (frames are pre-rendered images) |
| Needs PyMOL | No (mesh is included in the repo) | **Yes** |
| Time | seconds | minutes to ~half an hour (ray tracing) |

There is also an older **Cairo scene** (`scenes/protein_scene.py`), kept because it matches the original spec (`ThreeDScene`, default renderer). It draws a decimated mesh and looks noisier than A.

> Tested on macOS (Apple Silicon), Python 3.12, Manim 0.22. Linux should work with the packages listed below. Windows is untested.

---

## 1. Setup

### 1a. System packages (once)

**macOS** (with [Homebrew](https://brew.sh)):
```bash
brew install ffmpeg cairo pango pkg-config
```
**Ubuntu/Debian:**
```bash
sudo apt install ffmpeg libcairo2-dev libpango1.0-dev pkg-config python3-dev python3-venv
```
Python 3.10 or newer is required. `sox` is *not* needed (audio only).

### 1b. Project setup

```bash
git clone <repo-url> protein_manim_project
cd protein_manim_project
./setup.sh
```
`setup.sh` checks your system packages, creates `.venv`, and installs `requirements.txt`. That is all you need for the OpenGL scene.

<details>
<summary>Manual equivalent (what setup.sh does)</summary>

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```
`requirements-lock.txt` has the exact versions the project was developed with (macOS-specific packages included); use it only if the plain `requirements.txt` gives you trouble.
</details>

### 1c. PyMOL (only for the hybrid scene / regenerating the mesh)

PyMOL is installed in a **separate conda environment**, not in `.venv`. (The `pymol-open-source` wheel on PyPI is broken on macOS, so use conda-forge.) Install [Miniforge](https://github.com/conda-forge/miniforge) if you don't have conda, then:
```bash
./setup.sh --with-pymol            # or: conda create -n pymol -c conda-forge pymol-open-source
```
The scripts find PyMOL automatically (`$PYMOL`, `pymol` on PATH, or a conda env named `pymol` under `~/miniforge3`, `~/miniconda3`, `~/mambaforge`, ...). If yours lives elsewhere: `export PYMOL=/path/to/pymol`.

---

## 2. Render

Use the helper `./render.sh <mode>` (extra arguments go to `manim`), or run the commands shown underneath.

### A. OpenGL scene (no PyMOL needed)
```bash
./render.sh gl                     # 720p30
./render.sh gl -qh                 # 1080p60
./render.sh gl -qk                 # 4K60
```
Equivalent: `manim -r 1280,720 --fps 30 --renderer=opengl scenes/protein_scene_gl.py ProteinVolumeSceneGL`

### B. Hybrid scene (PyMOL ray-traced frames + Manim overlay)
Step 1 — ray-trace the turntable (60 frames of a 360° spin, transparent PNGs):
```bash
./render.sh frames-fast            # fast settings  (~6 min)
./render.sh frames-hq              # high quality   (~35 min)
```
Step 2 — composite with Manim:
```bash
./render.sh hybrid                 # uses frames-hq if present, otherwise frames-fast
FRAME_DIR=assets/pymol_frames ./render.sh hybrid     # force a specific frame set
```
Frames are resumable: if you stop the PyMOL step, re-running it skips frames already written.

### C. Cairo scene (original spec)
```bash
./render.sh cairo                  # 480p15 preview
./render.sh cairo -qh              # 1080p60
```

Videos are written to `media/videos/...` (ignored by git). Add `-p` to open the result when finished.

### Regenerate the surface mesh (optional)
`assets/6QL9_surface.obj` is included in the repo (200k faces, per-chain vertex colors), so you can skip this. To rebuild it (downloads `6QL9.cif` from RCSB if missing):
```bash
./render.sh mesh                                   # PyMOL if installed, else built-in fallback
./render.sh mesh --backend builtin                 # force the no-PyMOL fallback (blobbier surface)
```

---

## 3. Computation times

Measured on an Apple Silicon Mac, 10 CPU cores. Your numbers will differ; ray tracing scales with core count.

| Step | Settings | Time |
|---|---|---|
| Generate mesh | PyMOL backend | ~35 s |
| Generate mesh | built-in fallback (no PyMOL) | ~2 s |
| **A. OpenGL render** | 480p15 | ~8 s |
| **A. OpenGL render** | 720p30 | ~23 s |
| **B. PyMOL frames, fast** | 60 frames, 1280×720, antialias 1, reduced-triangle surface | ~6 min (~6 s/frame after a ~25 s load) |
| **B. PyMOL frames, high quality** | 60 frames, 1280×720, antialias 2, full surface | ~35 min (~35 s/frame) |
| **B. Hybrid compose (Manim)** | 720p30, 20.5 s video | ~55 s |
| C. Cairo render | 480p15, 4,000-face mesh | ~1 min |

Not measured: OpenGL/Cairo at 1080p/4K, and Cairo with the current 15,000-face setting (expect it to be noticeably slower than the 4,000-face figure above). The hybrid scene always composites 1280×720 frames; rendering it at 1080p+ just upscales them (re-run PyMOL with larger `--width/--height` for true 1080p; time grows roughly with pixel count).

**Why PyMOL is slow:** it ray-traces every frame on the CPU (no GPU path), with a ~2.4M-triangle surface, shadows, and antialiasing (antialias 2 = ~4 rays/pixel). Speed-ups: fewer `--frames`, smaller `--width/--height`, `--antialias 1`, `--surface-quality -1`.

---

## 4. Troubleshooting

| Problem | Fix |
|---|---|
| `pycairo`/`pangocairo` build errors on pip install | Install the system packages in 1a, then re-run `./setup.sh` |
| `ffmpeg: command not found` | `brew install ffmpeg` / `apt install ffmpeg` |
| `PyMOL not found` | Do 1c, or `export PYMOL=/path/to/pymol` |
| `No frames in assets/pymol_frames` (hybrid) | Run `./render.sh frames-fast` first |
| OpenGL window/context errors on a headless Linux box | Use a machine with a display/GPU, or fall back to the Cairo scene |
| Hybrid spin looks too fast/slow | `REV_SECONDS` in `scenes/protein_scene_hybrid.py` |

---

## 5. More

- **How the code works:** [`docs/CODE_GUIDE.md`](docs/CODE_GUIDE.md)
- Project layout:
```
assets/6QL9_surface.obj        200k-face colored surface mesh (committed)
assets/6QL9.cif                downloaded on demand (ignored by git)
assets/pymol_frames*/          ray-traced turntable frames (ignored by git)
scenes/protein_scene_gl.py     A. OpenGL scene
scenes/protein_scene_hybrid.py B. PyMOL frames + Manim overlay
scenes/protein_scene.py        C. Cairo scene
scripts/fetch_and_convert.py   download CIF -> surface mesh (PyMOL / ChimeraX / built-in)
scripts/render_pymol_frames.py launcher for the PyMOL turntable
scripts/pymol_turntable.py     runs inside PyMOL
setup.sh, render.sh            setup and render helpers
```
