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

- `setup.sh` — verifies Python ≥ 3.10 and system packages, creates `.venv`, installs `requirements.txt`; `--with-pymol` also creates the conda env.
- `render.sh` — one entry point for `mesh | gl | cairo | frames-fast | frames-hq | hybrid`; extra args pass through to `manim`.
