# Examples

Four ready-made examples exercise every pipeline with real data from the [RCSB PDB](https://www.rcsb.org). **You do not have to download anything by hand**: each example fetches the files it needs, the first time you run it.

```bash
./render.sh example list                    # show the examples
./render.sh example run adk --quick         # try one in seconds (small, fast render)
./render.sh example run adk                 # full-quality render
```
(Needs `./setup.sh --with-pymol` once; see the [README](../README.md#1-setup).)

| Example | What it shows | Input data (from RCSB) | Pipeline |
|---|---|---|---|
| `adk` | **Shape change.** Adenylate kinase closing its lid: apo/open [4AKE](https://www.rcsb.org/structure/4AKE) → inhibitor-bound/closed [1AKE](https://www.rcsb.org/structure/1AKE) (the Ap5A ligand itself is ignored) | 2 PDB files, 0.66 MB, **committed** in `examples/adenylate_kinase/` (re-downloaded if missing) | conformational change (`render_traj.py`) |
| `nmr` | **Shape change.** Ubiquitin NMR ensemble: each of the 10 models is one frame ([1D3Z](https://www.rcsb.org/structure/1D3Z)) | 1 PDB file, 1 MB, downloaded to `examples/data/nmr_ubiquitin/` | conformational change |
| `cypa278` | **Data change.** Cyclophilin A at 278 K, 9 successive datasets ([5KUZ … 5KV7](https://www.rcsb.org/structure/5KV6)): electron density and difference density | 9 models + 9 structure-factor files, 15 MB, downloaded to `examples/data/cypa_278K/` | electron density (`render_density.py`) |
| `cypa100` | **Data change.** The same experiment at 100 K ([5KUL … 5KUW](https://www.rcsb.org/structure/5KUL)) | 9 + 9 files, 14 MB, downloaded to `examples/data/cypa_100K/` | electron density |

Only the small adenylate kinase files are committed to git. Everything under `examples/data/` is downloaded on demand and git-ignored, so the repository stays small. To download without rendering: `python scripts/examples.py fetch all` (or one example name).

## What to expect

| Example | Quick (`--quick`) | Full quality | Result to look for |
|---|---|---|---|
| `adk` | ~7 s (640×360) | ~36 s | 1,656 atoms common to both files; CA RMSD rises to **7.1 Å**; the lid domain turns red as it closes |
| `nmr` | ~25 s (640×360) | ~54 s | 602 atoms in all 10 models; max CA RMSD **1.7 Å**; the core barely moves, the flexible tail end moves most (up to 5.6 Å) |
| `cypa278` | ~22 s (640×360) | ~4 min 38 s | R_iso vs data set 1 rises **3.2% → 11.4%**; difference density (green/red) appears in the late datasets while 2Fo−Fc stays stable |
| `cypa100` | ~22 s (640×360) | ~4 min 4 s | R_iso rises **3.3% → 16.5%** (stronger than at 278 K) |

Times are from a 10-core Apple Silicon Mac with the data already downloaded; the one-time download adds a few seconds (up to ~15 s for the CypA series). Ray tracing dominates and scales with your CPU cores. "Full quality" means the default settings of the pipeline scripts (1280×720 video, antialias 2).

## Where the output goes

```
assets/examples/<name>/media/videos/<scene>/<quality>/*.mp4          full quality
assets/examples/<name>_quick/media/videos/...                         --quick
```
The path is also printed at the end of every run; open that file with your video player. These folders are git-ignored; delete `assets/examples/` to reclaim space.

## Changing an example

Any extra option is passed to the underlying render script, and later options override the example's defaults:
```bash
./render.sh example run adk --title "My title" --rotate 0 --between 40
./render.sh example run cypa278 --set-seconds 4 --manim-args "-r 1920,1080 --fps 30"
```
Options of the two pipelines: `python scripts/render_traj.py --help` (adk, nmr) and `python scripts/render_density.py --help` (cypa278, cypa100). Re-running with only video-level changes (title, durations, labels, resolution) reuses the ray-traced frames, so it takes seconds, not minutes.

## The examples are templates for your own data

Each example is just a call to a render script. These are the equivalent commands, which you can copy and point at your own files (step-by-step in [`docs/CUSTOM_DATA.md`](../docs/CUSTOM_DATA.md)):

```bash
# adk  (shape change from two coordinate files)
./render.sh traj examples/adenylate_kinase/4AKE.pdb examples/adenylate_kinase/1AKE.pdb \
    --labels "Open (4AKE),Closed (1AKE)" --chains A --title "Adenylate kinase: lid closing"

# nmr  (one multi-model file; every model is a frame)
./render.sh traj examples/data/nmr_ubiquitin/1D3Z.pdb --between 5 --title "Ubiquitin: NMR ensemble (PDB 1D3Z)"

# cypa278  (data change from model + structure-factor files; <ID>.cif and <ID>-sf.cif per dataset)
./render.sh density --dir examples/data/cypa_278K --ids 5KUZ,5KV0,5KV1,5KV2,5KV3,5KV4,5KV5,5KV6,5KV7 \
    --focus 55,60,61,63,72,101,102,103,111,113,121,122,126 --zoom-label "Active site"
```
Your own data only needs the same layout: a folder with one model and one structure-factor file per dataset (or one coordinate file per snapshot). Details and accepted formats: [`docs/CUSTOM_DATA.md`](../docs/CUSTOM_DATA.md).

## Offline use / troubleshooting
- **No internet:** run an example once with internet to cache the files, or place the files yourself in the folders above (`<ID>.pdb`; `<ID>.cif` + `<ID>-sf.cif`). Existing files are never re-downloaded.
- **A download fails or a file is corrupt:** delete the file (or `examples/data/`) and run again.
- **`PyMOL not found`:** run `./setup.sh --with-pymol` or set `PYMOL=/path/to/pymol`.
- **Cyclophilin A zoom region:** the examples zoom on the catalytic residues (Arg55, Phe60, Met61, Gln63, Gly72, Ala101, Asn102, Ala103, Gln111, Phe113, Trp121, Leu122, His126). That choice is a modelling decision of this project, not something derived from the data; change it with `--focus`.
