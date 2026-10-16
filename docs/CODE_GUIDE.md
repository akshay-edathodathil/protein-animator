# Code guide

How the pieces fit together:

```
6QL9.cif ──► fetch_and_convert.py ──► assets/6QL9_surface.obj ──► protein_scene_gl.py  (A: OpenGL)
   │            (PyMOL surface,         (200k faces, colored        protein_scene.py     (C: Cairo)
   │             chain coloring)         per chain)
   │
   └──► render_pymol_frames.py ──► pymol_turntable.py (inside PyMOL) ──► assets/pymol_frames*/*.png
                                                                              │
                                                                  protein_scene_hybrid.py  (B: hybrid)
```

## `scripts/fetch_and_convert.py` — structure → colored surface mesh

1. `fetch()` uses `assets/6QL9.cif`, or downloads it from `files.rcsb.org`.
2. `read_atoms()` (gemmi) returns heavy-atom coordinates and the **chain index of every atom**. The file already contains the whole biological assembly (12 chains: α = A–F, β = G–L), so no symmetry expansion is needed.
3. A backend produces raw triangles; the first one available is used (`--backend auto`):
   - `run_pymol()` — headless PyMOL: `show surface` then `save raw.obj`. This is a real solvent-excluded molecular surface. Output is a ~2.4M-triangle "soup" with no colors.
   - `run_chimerax()` — same idea with ChimeraX (untested here; kept as an alternative).
   - `run_builtin()` — no external tools: atoms → Gaussian-blurred density grid → marching cubes. Much blobbier, but works anywhere.
4. `finalize()` is shared by all backends:
   - welds duplicate vertices, drops degenerate faces,
   - optionally `keep_exterior()` (`--exterior-only`) removes interior faces with a solvent flood fill — off by default because with a GPU depth buffer it isn't needed and it left holes,
   - decimates to `--max-faces` (default 200k),
   - **colors each vertex by the chain of its nearest atom** (KD-tree; hue wheel, one color per chain) — PyMOL's OBJ export carries no colors, so this is done here,
   - writes an OBJ with per-vertex RGB (`v x y z r g b`).

## A. `scenes/protein_scene_gl.py` — OpenGL scene

Rendered with `--renderer=opengl`. Manim's OpenGL renderer has a real depth buffer, so the full 200k-face mesh renders correctly in seconds.

- `MeshSurface(OpenGLSurface)` adapts Manim's *parametric* surface class to an arbitrary triangle mesh:
  - `triangle_indices` come straight from the mesh faces,
  - Manim's surface shader derives normals from two extra points nudged along the surface tangent; `init_points()` builds those from the mesh's vertex normals (tangents orthogonal to the normal),
  - `get_shader_data()` supplies per-vertex colors (the chain colors).
- `load_mesh_mobject()` loads the OBJ with trimesh, centers it at the origin and scales the longest side to `TARGET_WIDTH` (4 units).
- `ProteinVolumeSceneGL` sets the camera (φ=70°, θ=−45°), fades in a fixed-in-frame title, then steps through three captions while the camera orbits (`camera.animate.increment_theta`) and zooms (`camera.animate.scale`), ending with a top-down view.

Tweak: caption text (`notes`), camera moves/timings in `construct`, `TARGET_WIDTH`, shader `gloss`/`shadow`.

## B. Hybrid scene

### `scripts/pymol_turntable.py` (runs *inside* PyMOL)
Loads the CIF, removes solvent/hydrogens, shows a surface, colors by chain (`util.cbc`), fits the view, then loops `png(..., ray=1)` + `turn("y", 360/N)` to write `protein_frame0001.png…` Transparent background comes from `ray_opaque_background off`. It skips existing frames, so it is resumable. Quality knobs (positional args): frames, width, height, antialias, surface_quality.

### `scripts/render_pymol_frames.py` (launcher)
Finds PyMOL (`find_pymol()` in `fetch_and_convert.py`), then runs `pymol -cq pymol_turntable.py -- <args>`. Flags: `--frames --width --height --antialias --surface-quality --out`.

### `scenes/protein_scene_hybrid.py`
A plain 2D `Scene`; the 3D look is baked into the PNGs.
- `load_frames()` reads all PNGs into one NumPy array (from `$FRAME_DIR`, default `assets/pymol_frames`).
- An **updater** on a single `ImageMobject` advances a clock each rendered frame, computes a fractional frame index (`t / REV_SECONDS * n`) and **cross-fades the two neighboring PNGs**, so a 60-frame (6° per step) turntable still looks smooth. One revolution takes `REV_SECONDS` (30 s); the 20-second video shows about 240°.
- Title, captions and the zoom (`img.animate.scale`) are ordinary Manim animations. The image is sized to 78% of the frame height and shifted up so captions don't overlap the protein.

Why not one `ImageMobject` per frame with `remove/add` + `wait(1/30)`? Every `wait` becomes its own animation; with 60 fps output that duplicates frames and runs the spin at the wrong speed.

## C. `scenes/protein_scene.py` — Cairo scene

The default Manim renderer has **no depth buffer**: it paints triangles back-to-front by their centers. That is why this scene (a) decimates the mesh to `TARGET_FACES` (15,000) and builds one `Polygon` per face, and (b) still looks noisier than A, because interior surfaces interleave with the outer shell. It colors per face from the OBJ vertex colors (falls back to a Z-height gradient if the OBJ has none), with simple Lambert-style shading. Kept mainly because it follows the original spec.

## Helper scripts

- `scripts/examples.py` — runs the built-in examples (`adk`, `nmr`, `cypa278`, `cypa100`): downloads missing inputs from RCSB (`urlretrieve`, never overwriting existing files), builds the command for `render_traj.py` / `render_density.py` (`--quick` swaps in small, fast settings and a separate output folder), and forwards unknown options to the render script (`parse_known_args`, so later options override the example's defaults).
- `setup.sh` — verifies Python ≥ 3.10 and system packages, creates `.venv`, installs `requirements.txt`; `--with-pymol` also creates the conda env.
- `render.sh` — one entry point for `mesh | gl | cairo | frames-fast | frames-hq | hybrid`; extra args pass through to `manim`.

## D. Structure change over time (`render_traj.py` pipeline)

```
PDB files / models ──► traj_prepare.py ──► morph.pdb + disp.f32 + meta.json ──► pymol_traj.py ──► frames/*.png ──► trajectory_scene.py ──► video
```

- **`scripts/traj_prepare.py`** (runs in `.venv`, uses gemmi + numpy)
  - `collect_inputs()` expands files/directories in natural order; `read_keyframes()` makes one keyframe per file *or per MODEL*, keyed by `(chain, resnum, icode, resname, atom name)`.
  - `prepare()` keeps atoms common to all keyframes, superimposes each on the first (`kabsch()` on CA atoms), then builds `--between` interpolated frames per segment with smoothstep easing. It also computes per-atom displacement from keyframe 0 and per-frame CA RMSD.
  - `write_outputs()` writes `morph.pdb` (one MODEL per frame, identical atom order), `disp.f32` (raw float32, frames × atoms — read in PyMOL without numpy), and `meta.json` (labels, keyframe flags, RMSD, `disp_max`).
- **`scripts/pymol_traj.py`** (runs inside PyMOL): loads `morph.pdb` as multi-state, assigns secondary structure, frames the camera on **all states at once** (`orient/zoom state=0`) so the protein never leaves the view, then per frame: `set state`, write that frame's displacements into the B-factor column (`alter ... b=next(it)`; B is per-atom, not per-state, so it must be reset each frame), `spectrum b, blue_white_red`, `png(..., ray=1)`, and a small `turn y` for the slow rotation.
- **`scenes/trajectory_scene.py`**: a playback schedule gives each real structure `HOLD` seconds and each interpolated frame `1/PLAY_FPS`; an updater picks the current PNG and updates the label, RMSD text and timeline marker. **Gotcha:** Manim's Cairo renderer caches mobjects without updaters as a static background image, so anything that changes during playback (image, marker, labels) must carry an updater; labels and RMSD texts are pre-created and toggled by opacity instead of added/removed.
- **`scripts/render_traj.py`** runs the three stages in order and prints the video path.

## E. Density change across X-ray datasets (`render_density.py` pipeline)

```
<ID>.cif + <ID>-sf.cif (per dataset) ──► density_prepare.py ──► maps/*.ccp4 + meta.json ──► pymol_density.py ──► frames/{overview,zoom}/*.png ──► density_scene.py ──► video
```

- **`scripts/density_prepare.py`** (`.venv`; gemmi + numpy)
  - `load_dataset()` converts the structure-factor mmCIF to an MTZ-like table with gemmi (`IMEAN, SIGIMEAN, FC, PHIC`; some reflections have no F_calc, so NaNs are masked).
  - `scaled_fo()` sets Fo = sqrt(max(I,0)) and fits `ln(Fc/Fo) = a + b/d²` in 15 resolution bins to put every dataset on the |Fc| scale.
  - `make_grid()` Fourier-transforms coefficients to a map (`Mtz.transform_f_phi_to_map`). 2Fo−Fc uses each dataset's own phases; the difference map uses `Fo(n) − Fo(1)` (after one scalar scale) with **set 1's phases**.
  - All maps are written in one absolute unit (σ of set 1's 2Fo−Fc map) and cropped around the molecule (`write_map`). `diff_contour` = n·σ of the last difference map and is reused for every frame.
  - `R_iso = Σ|Fo(n)−Fo(1)| / ΣFo(1)` per dataset is the overall change measure used for the bar chart.
- **`scripts/pymol_density.py`** (inside PyMOL): loads all models and maps with **`normalize_ccp4_maps 0`** (otherwise PyMOL silently rescales every map to its own σ and the contours stop being comparable — an early version of this pipeline showed noise everywhere because of that). Per frame it draws `isomesh` objects for the current dataset only: overview = difference ± contour around all atoms; zoom = 2Fo−Fc at 1σ + difference around the focus residues, with a semi-transparent cartoon behind.
- **`scenes/density_scene.py`**: shows `frames_per_set` frames per dataset for `SET_SECONDS`, dissolving only between datasets, plus a bar chart with the current dataset highlighted. As in the trajectory scene, anything that changes during playback (images, bars, readouts) carries an updater so Cairo does not cache it as background.
