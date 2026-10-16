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

## 3. Show how a protein changes over time (multiple PDB files)

Give the pipeline several structures in time order (separate PDB/mmCIF files, a folder, or a multi-model file such as an NMR ensemble or MD snapshots). It produces a video with a (nearly) fixed camera where you watch the structure move.

```bash
./setup.sh --with-pymol                       # PyMOL is required
./render.sh traj-demo                         # example: adenylate kinase, open (4AKE) -> closed (1AKE)

# your own data (files are used in natural filename order: t2.pdb before t10.pdb)
./render.sh traj snapshots/*.pdb --labels "0 ns,10 ns,20 ns" --title "My protein"
./render.sh traj my_folder/ --chains A
./render.sh traj ensemble.pdb                 # every MODEL in the file is one frame
```
The video path is printed at the end (`assets/trajectory/media/videos/trajectory_scene/...`).

**What it does**
1. Keeps only atoms present in *all* structures (matched by chain, residue number, residue name, atom name; waters, hydrogens and ligands are dropped).
2. Superimposes every structure on the first (CA atoms), so only real conformational change is visible.
3. Inserts `--between N` interpolated frames between consecutive structures (default 24) so the motion is smooth. Use `--between 0` if your files are already closely spaced MD snapshots.
4. PyMOL ray-traces every frame as a cartoon from the same viewpoint (plus a slow `--rotate` 20° drift so it feels like video; `--rotate 0` locks the view). Residues are colored blue → white → red by **how far they have moved from the first structure**.
5. Manim adds the title, current structure label, a timeline with a tick for each real structure, a live CA-RMSD readout and a color legend. The video pauses (`--hold`, default 1.2 s) on each real structure.

**Common options** (`python scripts/render_traj.py --help` for all)

| Option | Meaning |
|---|---|
| `--labels "a,b,c"` | one label per structure (default: file names) |
| `--chains A,B` | keep only these chains (needed if chain IDs differ between files) |
| `--between N` | interpolated frames between structures |
| `--color-by displacement\|chain\|none` | coloring |
| `--disp-max Å` | displacement that maps to full red (default: 99th percentile) |
| `--ss first\|each` | secondary structure assigned once (stable) or recomputed every frame |
| `--width/--height/--antialias` | PyMOL render quality (default 1280×720, antialias 2) |
| `--play-fps`, `--hold` | playback speed through interpolated frames / pause on real structures |
| `--prepare-only`, `--no-video` | stop after alignment, or after the frames |

**Requirements and limitations**
- The structures must be the *same protein* with matching residue numbering and chain IDs (different constructs work, but only residues common to all are shown).
- Interpolation is straight-line in Cartesian coordinates, so for large domain motions the in-between frames are an approximation (bonds can briefly look compressed). It is a visual guide, not a simulation. With real MD/NMR snapshots, use `--between 0` or a small value.
- Chain IDs must be one character (PDB format).
- Resume: frames already rendered are skipped if you re-run with the same `--out`.

**Computation times** (this Mac, 10 cores): adenylate kinase (214 residues), 2 structures, 30 interpolated frames, 1280×720, antialias 2 → **~1 min total** (~0.4 s/frame ray tracing). NMR ensemble (76 residues), 10 models, 55 frames at 960×540, antialias 1 → ~37 s. Time scales with number of frames × protein size (a huge complex like 6QL9 as *cartoon* would be far slower per frame; not measured).

---

## 4. Show how the *data* change across a series (structure factors -> electron density)

Use this when the refined PDB models of a series barely differ but the diffraction data do (e.g. successive datasets on one crystal, temperature series). Each dataset needs a model (`<ID>.cif`) **and** its structure-factor file (`<ID>-sf.cif`), both from RCSB.

```bash
./setup.sh --with-pymol
./render.sh density                           # CypA at 278 K, data sets 1-9 (5KUZ, 5KV0 ... 5KV7)
./render.sh density --ids 5KUL,5KUN,5KUO,5KUQ,5KUR,5KUS,5KUU,5KUV,5KUW --title "CypA, 100 K"   # the 100 K series
./render.sh density --ids ID1,ID2,ID3 --focus 12,45,67 --title "My series"                     # your own entries
```
Files are read from `ProteinTimeSeries/` (`<ID>.cif`, `<ID>-sf.cif`) and downloaded from RCSB when missing. The video path is printed at the end (`assets/density/media/...`).

**What it does** (`scripts/density_prepare.py`)
1. Reads intensities, F_calc and phases from each `-sf.cif`; Fo = sqrt(max(I, 0)), scaled onto |Fc| with a resolution-dependent factor.
2. **2Fo−Fc map** per dataset (phases from the deposited model).
3. **Difference map Fo(n) − Fo(1)** with set 1's phases (isomorphous difference density), plus **R_iso** per dataset: the overall size of the change.
4. PyMOL ray-traces two panels per frame: the whole molecule with the difference density (green = gain, red = loss) and a zoom on the active-site residues with 2Fo−Fc (blue) plus difference density. The **same absolute contour** is used for every frame (±3 σ of the *last* difference map), so growth of the change is visible rather than normalised away.
5. Manim composites the panels with a title, the current dataset, and a bar chart of R_iso across the series.

**Notes and limitations**
- The deposited structure-factor files contain intensities + F_calc/phase, not ready-made map coefficients, so the maps here are *not* sigma-A weighted (no refinement program is run). They are good for visualising change; use Phenix/Refmac maps for publication-quality density.
- The default zoom region (`--focus`) is the CypA active site; use `--focus auto` to pick the residues with the strongest difference peaks, or give your own residue numbers.
- Difference maps are in each entry's own crystal frame; cell edges drift by <0.3 Å across the CypA series, which is negligible for this visualisation.
- For 5KUZ-5KV7 the difference density is real but diffuse (peaks up to ~5 σ spread over many residues; R_iso rises 3% -> 11%), not a single dramatic site. The video shows exactly that.

**Computation times** (this Mac, 10 cores, CypA 278 K, 9 datasets): density maps + R_iso ~1 s (plus ~5 s the first time, to download the missing entries); PyMOL ray tracing of 144 panel images (8 frames per dataset x 2 panels, 800 px, antialias 2) **~2.7 min**; Manim compositing at 1080p30 (24.7 s video) **~1 min**. Total ~4 min.

---

## 5. Computation times (6QL9 videos)

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

## 6. Troubleshooting

| Problem | Fix |
|---|---|
| `pycairo`/`pangocairo` build errors on pip install | Install the system packages in 1a, then re-run `./setup.sh` |
| `ffmpeg: command not found` | `brew install ffmpeg` / `apt install ffmpeg` |
| `PyMOL not found` | Do 1c, or `export PYMOL=/path/to/pymol` |
| `No frames in assets/pymol_frames` (hybrid) | Run `./render.sh frames-fast` first |
| OpenGL window/context errors on a headless Linux box | Use a machine with a display/GPU, or fall back to the Cairo scene |
| Hybrid spin looks too fast/slow | `REV_SECONDS` in `scenes/protein_scene_hybrid.py` |

---

## 7. More

- **How the code works:** [`docs/CODE_GUIDE.md`](docs/CODE_GUIDE.md)
- Project layout:
```
assets/6QL9_surface.obj        200k-face colored surface mesh (committed)
assets/6QL9.cif                downloaded on demand (ignored by git)
assets/pymol_frames*/          ray-traced turntable frames (ignored by git)
scenes/protein_scene_gl.py     A. OpenGL scene
scenes/protein_scene_hybrid.py B. PyMOL frames + Manim overlay
scenes/protein_scene.py        C. Cairo scene
scenes/trajectory_scene.py     structure-change video (overlay for the PyMOL frames)
scripts/fetch_and_convert.py   download CIF -> surface mesh (PyMOL / ChimeraX / built-in)
scripts/render_pymol_frames.py launcher for the PyMOL turntable
scripts/pymol_turntable.py     runs inside PyMOL
scripts/render_traj.py         several PDBs -> change-over-time video (orchestrates the two below)
scripts/traj_prepare.py        align + interpolate structures
scripts/pymol_traj.py          runs inside PyMOL: fixed-camera frames colored by displacement
scripts/render_density.py      X-ray datasets (model + structure factors) -> density-change video (orchestrator)
scripts/density_prepare.py     structure factors -> 2Fo-Fc / difference maps, R_iso
scripts/pymol_density.py       runs inside PyMOL: density-mesh panels
scenes/density_scene.py        Manim overlay for the density video
examples/adenylate_kinase/     4AKE (open) and 1AKE (closed) demo inputs
setup.sh, render.sh            setup and render helpers
```
