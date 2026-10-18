# Using your own protein time-series data

Two pipelines turn a series of your own files into a Manim video. Pick by **what changes** in your data:

| Your data | What you want to show | Use |
|---|---|---|
| Several **coordinate files** (PDB/mmCIF; MD snapshots, NMR models, different conformations) whose atoms move | The protein *moving* (cartoon, colored by how far each residue moved) | **A. Conformational change** (`render.sh traj`) |
| Several **crystallographic datasets** (model + structure factors) where the refined models barely differ but the data do (dose series, temperature series, time-resolved) | Change in **electron density** across the series | **B. Density change** (`render.sh density`) |

**New here? Run the matching example first** (`./render.sh example run adk --quick`, or `cypa278`); [`examples/README.md`](../examples/README.md) shows the equivalent command for each pipeline, ready to copy. Both pipelines need PyMOL (ray tracing) and the Python environment. If you haven't set up yet, do [README §1](../README.md#1-setup) (`./setup.sh --with-pymol`) first. All commands below are run from the project folder.

---

## A. Conformational change from coordinate files

### 1. Prepare the files
- Put the structures in one folder (or list them on the command line), **in time order**. When you pass a **folder**, files are sorted in *natural* order (`snap_2.pdb` before `snap_10.pdb`). When you list **files** (including a shell glob such as `snap_*.pdb`), your order is used as given, and the shell sorts globs alphabetically (`snap_10` before `snap_2`), so for globs use zero-padded names (`snap_002.pdb`) or just pass the folder.
- Formats: `.pdb`, `.ent`, `.cif` (+ `.gz`), or **multi-model files** (every `MODEL` is one frame, e.g. an NMR ensemble or an MD trajectory exported as a multi-model PDB).
- Requirements: the **same protein** in every file with the **same residue numbering and chain IDs** (constructs that differ at the ends are fine; only residues/atoms present in *all* structures are shown). Chain IDs must be a single character.
- Waters, hydrogens, ligands and alternate locations other than A are ignored. Frames are superimposed on the first one automatically.

### 2. Run
```bash
# a folder of snapshots
./render.sh traj my_snapshots/ --title "My protein" --labels "0 ns,10 ns,20 ns,30 ns"

# explicit files, one chain only
./render.sh traj run/t0.pdb run/t1.pdb run/t2.pdb --chains A --labels "t0,t1,t2"

# one multi-model file (every MODEL is a frame)
./render.sh traj ensemble.pdb --between 4
```
At the end it prints the path of the video (`assets/trajectory/media/videos/trajectory_scene/<quality>/TrajectoryScene.mp4`).

### 3. Options you will actually use
| Option | Effect |
|---|---|
| `--labels "a,b,c"` | One label per structure (default: file names; for multi-model files "name · model k") |
| `--between N` | Interpolated frames between consecutive structures (default 24). **Use `--between 0`** (or a small number) if your snapshots are already closely spaced, e.g. MD frames; keep 20–30 for just a few distant structures |
| `--chains A,B` | Keep only these chains (needed when chain IDs differ between files) |
| `--color-by displacement\|chain\|none` | Blue→red by distance moved from the first structure (default), by chain, or plain |
| `--disp-max 8` | Å that map to full red (default: 99th percentile of all movement) |
| `--rotate 20` / `0` | Slow camera drift in degrees over the whole video (`0` = locked view) |
| `--ss first\|each` | Secondary structure assigned once (stable) or recomputed every frame (use `each` if helices/sheets form or melt) |
| `--hold 1.2`, `--play-fps 15` | Seconds to pause on each real structure; playback speed through interpolated frames |
| `--width 1920 --height 1080` | PyMOL render size (default 1280×720) |
| `--antialias 1` | Faster, rougher edges (default 2) |
| `--manim-args "-r 1920,1080 --fps 30"` | Video resolution/frame rate (default `-r 1280,720 --fps 30`) |
| `--out my_run/` | Output folder (default `assets/trajectory`); use one folder per data set so runs don't overwrite each other |
| `--prepare-only` / `--no-video` | Stop after alignment / after the PNG frames |

### 4. Check before a long render
```bash
./render.sh traj my_snapshots/ --prepare-only
```
Look at the printout: **"N atoms common to all"** should be most of your protein (the pipeline warns below 50%) and the **maximum CA RMSD** should be plausible. If few atoms are common, residue numbering or chain IDs differ between files; use `--chains` or renumber.

### 5. Re-make the video without re-rendering frames
Frames are cached with their settings. Re-running the *same* command with changed `--hold`, `--play-fps`, `--labels`, `--title` or `--manim-args` reuses the PNGs and only recomposes the video (about a minute). Changing anything that affects the picture (size, `--between`, `--color-by`, `--rotate`, files, ...) clears the cache and re-renders automatically.

### 6. How long will it take?
PyMOL time ≈ number of frames × time per frame; per-frame cost grows with protein size and render size. Measured on a 10-core Mac: adenylate kinase (214 residues), 2 structures, 30 interpolated frames, 1280×720: ~1 min total (~0.4 s/frame); a 76-residue NMR ensemble, 55 frames at 960×540: ~37 s. Large complexes (thousands of residues) cost much more per frame, so lower `--between`, `--width/--height`, or use `--antialias 1`.

---

## B. Electron-density change from structure factors

### Quickest way: just give it the folder
```bash
./render.sh density --dir /path/to/my_folder --prepare-only     # check what it found (about a second)
./render.sh density --dir /path/to/my_folder                     # make the video
```
With only `--dir`, the pipeline **finds the datasets itself**:
- **Models:** `.pdb`, `.ent`, `.cif`. **Reflections:** MTZ files, recognised by their content, so any or no file extension works (Refmac's `HKLOUT_...` files have none), or `*-sf.cif` files.
- **Pairing:** a leading `XYZOUT_` / `HKLOUT_` is ignored, as are the extension and a trailing `-sf`. So `XYZOUT_0.1and0.9.pdb` + `HKLOUT_0.1and0.9` form the dataset `0.1and0.9`, and `t00.pdb` + `t00.mtz` form `t00`.
- **Order:** by the numbers in the names, compared numerically including decimals (`0.15and0.85` comes before `0.2and0.8`, `t2` before `t10`). The **first dataset is the reference** for all difference maps. Add `--reverse` to flip the order.
- **It prints what it paired**, in series order, before doing anything else. Check this list: wrong pairing or order is the most common mistake. Unpaired files are listed in a warning and ignored; other files (README, `.DS_Store`) are skipped.

Add `--labels "10/90,20/80,..."` (one per dataset, in the printed order) if you want nicer text in the video than the file names, plus `--title`, `--focus` etc. as below. If the automatic pairing is not what you want, name the datasets explicitly with `--ids` (next sections).

### 1. What you need per dataset
1. A **model**: `<name>.cif`, `<name>.pdb` or `<name>.ent` (the refined structure for that dataset).
2. **Structure factors** from the refinement of that model: `<name>-sf.cif`, `<name>.mtz` or `<name>-sf.mtz`.
   The file must contain **observed data**, **calculated structure factors** and **calculated phases**. Accepted column names (first match wins):

   | | Accepted columns |
   |---|---|
   | observed intensity + sigma | `IMEAN`/`SIGIMEAN`, `I`/`SIGI`, `I-obs`/`SIGI-obs` |
   | *or* observed amplitude + sigma | `F`/`SIGF`, `FP`/`SIGFP`, `F-obs`/`SIGF-obs`, `F-obs-filtered`/`SIGF-obs-filtered` |
   | calculated amplitude + phase (degrees), used as a **matched pair**, best first | `FC_ALL`/`PHIC_ALL`, `FC_ALL_LS`/`PHIC_ALL_LS`, `FC`/`PHIC`, `F-model`/`PHIF-model` |

   Where they come from: an RCSB `-sf.cif` already has `IMEAN, SIGIMEAN, FC, PHIC`. A **Phenix** refinement `*_data.mtz` has `F-obs-filtered`, `SIGF-obs-filtered`, `F-model`, `PHIF-model`. A **Refmac** output MTZ (the file you give as `HKLOUT`) has `FP, SIGFP, FREE, FC, PHIC, FC_ALL, PHIC_ALL, FWT, PHWT, ...`; it works as is. When both pairs are present the pipeline uses `FC_ALL`/`PHIC_ALL` (includes the bulk-solvent contribution) rather than the model-only `FC`/`PHIC`. A raw data-reduction MTZ (only `I, SIGI`) is *not* enough: it has no model phases; run a refinement first. If a column name is not recognised the error lists the columns found; rename them (e.g. with CCP4's `cad` or a few lines of gemmi/Python).
3. Files for **real PDB IDs** (e.g. `5KUZ`) are downloaded from RCSB automatically if missing. Use `--no-download` to forbid that.

#### Using Refmac output directly
Refmac's `XYZOUT` (coordinates) and `HKLOUT` (reflections) normally share a name prefix, e.g. `t00_refmac1.pdb` and `t00_refmac1.mtz`, or `XYZOUT_<name>.pdb` and `HKLOUT_<name>`. With `--dir` alone these are paired automatically (see above). To name the datasets explicitly, use the shared prefix as the dataset name:
```bash
./render.sh density --dir my_refmac_runs --ids t00_refmac1,t10_refmac1,t20_refmac1 --labels "0 min,10 min,20 min"
```
If your two files do not share a stem, or the reflection file has another extension, copy or rename it to `<stem>.mtz` next to the PDB. The pipeline does not care how the MTZ file is called otherwise; it reads columns by name, so unmeasured reflections (empty `FP`) are fine.

### 2. Rules for a valid series
- **Same crystal form**: same space group (checked; the run stops otherwise) and a similar cell (warns above 3% difference).
- **Same indexing/origin**: datasets should be indexed consistently with the first (true for datasets from the same crystal processed the same way). Models must be in the same orientation/origin and have the same residue numbering. If you refined each dataset independently, start each refinement from the previous model without moving the molecule.
- The **order of `--ids` is the time order**; the first is the reference for all difference maps.

### 3. Run
Put the files in `ProteinTimeSeries/` (or another folder via `--dir`), named `<id>.pdb/.cif` + `<id>.mtz/-sf.cif`:
```
my_series/
  t00.pdb   t00.mtz
  t10.pdb   t10.mtz
  t20.cif   t20-sf.cif
```
```bash
./render.sh density --dir my_series --ids t00,t10,t20 \
    --labels "0 min,10 min,20 min" --title "My enzyme: time-resolved series" \
    --focus 55,60,113        # residues for the zoom panel (see below)
```
The video path is printed at the end (`assets/density/media/videos/density_scene/<quality>/DensityScene.mp4`).

### 4. Options you will actually use
| Option | Effect |
|---|---|
| `--ids a,b,c` | Dataset names in time order (file stems; real PDB IDs also download) |
| `--labels "0 min,10 min,..."` | Text shown for each dataset instead of "Data set k" |
| `--title`, `--subtitle`, `--zoom-label` | Text on the video (defaults are generic for custom data) |
| `--view "<get_view text>"` or `--view view.txt` | **Your own camera for the zoom panel**: the output of PyMOL's `get_view` (the `set_view ( ... )` block, or just its 18 numbers), pasted in quotes or saved in a text file whose path you give. Implies zoom panel only and no rotation. See "Choosing the camera yourself" below |
| `--residue A:113` | **Focus on one amino acid**: chain `A`, residue number `113`. Shows a tight view of that residue (with its neighbours as thin sticks), its 2Fo−Fc density and the difference density, **no rotation**, only the zoom panel, one frame per dataset (so it is also fast). The caption reads e.g. "Phe113 (chain A)". Residue numbers are those of your model (PDB numbering); a wrong chain/number stops immediately and lists the chains and residue ranges that exist |
| `--focus 55,60,A:113` | Residues for the zoom panel (several allowed; `CHAIN:NUMBER` or just `NUMBER`, which means the first chain that has that residue number). **Default for your own data is `auto`**: residues with the strongest difference peaks. When the change is diffuse, `auto` picks scattered residues; choosing a region you care about (active site, a loop) usually gives a better zoom |
| `--diff-nsigma 3` | Difference-map contour = this × σ of the *last* difference map, **the same absolute level in every frame**, so growth of the change is visible. Lower (2–2.5) to show weaker changes (more noise), raise to show only the strongest |
| `--dmin 1.8`, `--diff-dmin 2.0` | Resolution limits (Å) for 2Fo−Fc and difference maps. Use a coarser `--diff-dmin` (e.g. 2.5) if your high-resolution shell is noisy |
| `--size 800`, `--antialias 2` | PyMOL panel size (pixels) and smoothing |
| `--frames-per-set 8`, `--rotate 20` | Frames rendered per dataset and total slow rotation in degrees (`0` = fixed view) |
| `--fast` | **Draft preset for slow machines**: `--size 500 --antialias 1 --frames-per-set 2` (other options still override it) |
| `--total-frames 20` | Instead of `--frames-per-set`: about 20 frames in total, spread evenly over the datasets (e.g. 9 datasets -> 2 per dataset) |
| `--panels zoom` | Render **only the zoom panel** (the right-hand one), shown larger. Halves the work; `overview` renders only the whole-molecule panel; default `both` |
| `--set-seconds 2.5` | Video seconds per dataset |
| `--manim-args "-r 1920,1080 --fps 30"` | Video resolution/frame rate |
| `--out my_density/` | Output folder (default `assets/density`) |
| `--prepare-only`, `--no-video` | Stop after the maps/statistics, or after the PNG frames |

### 5. Check the numbers first (about a second)
```bash
./render.sh density --dir my_series --ids t00,t10,t20 --prepare-only
```
It prints **R_iso** of every dataset versus the first (overall size of the change), the difference-map rms, and the zoom residues. Before rendering, check that R_iso behaves as you expect for your experiment. Very large values (tens of percent) can indicate that datasets are not consistently indexed or that models sit in different origins/orientations (a rule of thumb, not a hard limit). The maps are written to `assets/density/maps/*.ccp4` and can be opened in Coot or PyMOL.

### 6. Re-make the video without re-ray-tracing
Re-run the same command with a different `--set-seconds`, `--title`, `--labels`, `--subtitle`, `--zoom-label` or `--manim-args`; the PNG frames are reused (settings are cached) and only the Manim step repeats (~1 min at 1080p). Changing `--size`, `--frames-per-set`, `--rotate`, `--focus`, `--diff-nsigma`, `--dmin`, the datasets, etc. re-renders the frames.

### 7. Slow ray tracing? Speed it up
PyMOL ray tracing dominates the run time, and it is proportional to **(number of images) x (pixels) x (antialiasing rays)**. Levers, biggest first:
1. **`--panels zoom`**: one image per frame instead of two.
2. **Fewer frames**: `--frames-per-set 1` or `2`, or `--total-frames 20`. The datasets are discrete anyway; frames within a dataset only add a slow rotation. (`--rotate 0` for a fixed view.)
3. **`--antialias 1`**: about 4x fewer rays than the default 2.
4. **Smaller panels**: `--size 500` instead of 900 (about 3x fewer pixels). 500 px is enough for a 720p video.
5. **Fewer mesh lines**: a noisy difference map draws very many mesh segments. Raise `--diff-nsigma` (e.g. 4) and/or pick a small `--focus` region.

Easiest: `--fast` combines 2-4. A good first try on a slow machine: `./render.sh density --dir my_folder --panels zoom --fast`. If a single frame still takes minutes, run with `--total-frames 9` first to see the per-frame time, then scale up. On a 10-core Mac, 36 images at 500 px / antialias 1 take about 10 s; the default settings (180 images at 900 px / antialias 2) take about 3 minutes. Slower or older CPUs, and datasets that produce dense difference meshes, can be much slower per image.

### 8. How long will it take?
Measured on a 10-core Mac with the 9-dataset CypA series (165 residues): maps ~1 s; PyMOL 144 panel images at 800 px, antialias 2: ~2.7 min (≈2.2 s per frame pair); Manim at 1080p30 (25 s video): ~1 min. A synthetic 3-dataset test at 500 px finished in 20 s in total. Cost scales with the number of frames and the size of the molecule/mesh.

#### Choosing the camera yourself (`--view`)
1. Open **PyMOL** and load the model of the **first** dataset (the reference), unchanged (do not `align`, `translate` or `rotate` the object, since the view is stored in the model's own coordinates).
2. Turn, move and zoom the molecule with the mouse until the frame looks the way you want.
3. Type `get_view` in PyMOL's command line and copy the printed block (`set_view ( ... )`).
4. Save it in a text file (e.g. `view.txt`), or paste it in quotes, and run:
```bash
./render.sh density --dir my_folder --residue A:113 --view view.txt
```
The meshes are drawn around the focus residue(s) (`--residue` or `--focus`), so point those at what is in your view; `--view` only replaces the camera. The video panel is square: PyMOL keeps the *vertical* extent of your window, so keep the object of interest near the middle of your PyMOL window. The same camera is used for every dataset, and there is no rotation unless you add `--rotate`.

### 9. How to read the video (and what not to over-read)
- **Bar chart = how much the diffraction data changed.** The height of each bar is **R_iso** of that data set compared with data set 1: `R_iso = Σ|Fo(n) − Fo(1)| / ΣFo(1)`, summed over all reflections the two data sets share (to the `--diff-dmin` resolution), after scaling Fo(n) to Fo(1) with one overall factor. So 0% = identical amplitudes; 10% = on average the structure-factor amplitudes differ by 10%. Data set 1 is the reference, so its bar is 0 by definition. The highlighted bar is the data set currently shown. It is an overall, whole-crystal number; it does not say *where* in the molecule the change is (that is what the difference density shows).
- Even two identical crystals measured twice give a non-zero R_iso (measurement noise, a few %), so compare bars with each other and look for a steady rise, not for the absolute value.
- Blue mesh = 2Fo−Fc at 1 σ (set 1 scale). Green/red = positive/negative Fo(n)−Fo(1) at the fixed contour.
- The maps come from the deposited/refined phases and are **not σ-A weighted**; they are for visualising change, not for model building or publication-quality figures.
- Early datasets can look empty by design: the contour is fixed from the *last* dataset. If everything looks empty, lower `--diff-nsigma`.

---

## Customising the Manim look
- Text, colours, layout: edit `scenes/trajectory_scene.py` (A) or `scenes/density_scene.py` (B). Both are plain Manim `Scene`s; positions are in Manim units (frame is ~14.2 × 8).
- Gotcha: Manim's default (Cairo) renderer caches anything **without an updater** as a static background. If you add an object that must change during playback, give it an updater (the scenes use `keep_live(...)` for this).
- To render the scene yourself without the launcher:
  ```bash
  TRAJ_DIR=assets/trajectory manim -r 1920,1080 --fps 30 scenes/trajectory_scene.py TrajectoryScene
  DENSITY_DIR=assets/density manim -r 1920,1080 --fps 30 scenes/density_scene.py DensityScene
  ```
  (set `DENSITY_FRAMES_PER_SET` and `DENSITY_SET_SECONDS` to match how the frames were rendered).
- For true 4K, render the PNGs larger as well (`--width 3840 --height 2160` for A, `--size 1600` for B) and use `--manim-args "-r 3840,2160 --fps 30"`; expect several times the PyMOL time.

## Troubleshooting
| Message / symptom | Meaning and fix |
|---|---|
| `Dataset 'x': missing files ...` | The folder does not contain `x.cif/.pdb` and `x-sf.cif/.mtz`; the message lists the expected names. Check `--dir` and `--ids` (names are case-sensitive) |
| `no ... column found ... Columns present: [...]` | The structure-factor file lacks Fo/σ, F_calc or phases, or uses other names. Use refinement output, or rename columns |
| `N atoms common to all (x%)` with a low x / `WARNING: few atoms in common` | Numbering/chain IDs differ between files. Use `--chains`, or renumber |
| `... is in space group ... but ...` | Datasets must share a crystal form for difference maps |
| `Existing frames ... are stale; deleting and re-rendering` | Normal: settings changed since the last run, so the frame cache was reset |
| `PyMOL not found` | `./setup.sh --with-pymol`, or `export PYMOL=/path/to/pymol` |
| Difference map is all noise | Raise `--diff-nsigma`, use a coarser `--diff-dmin`, or check indexing consistency |
| Zoom panel looks empty or off-target | Pick residues explicitly with `--focus` |
