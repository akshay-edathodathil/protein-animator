#!/usr/bin/env python3
"""Run the built-in examples. Input data are downloaded from the RCSB PDB automatically when missing.

    python scripts/examples.py list
    python scripts/examples.py run adk            # protein shape change, open -> closed
    python scripts/examples.py run nmr --quick    # --quick = small, fast render to try the pipeline in ~1 min
    python scripts/examples.py run adk --title "My title"   # other options are passed to the render script
    python scripts/examples.py run cypa278        # electron-density change, CypA 278 K series
    python scripts/examples.py fetch all          # only download the data

Examples (details and expected output: examples/README.md):
    adk       4AKE -> 1AKE adenylate kinase; coordinates differ (shape change)       data: examples/adenylate_kinase (committed)
    nmr       ubiquitin NMR ensemble (1D3Z, 10 models); every model is a frame        data: examples/data/ (downloaded)
    cypa278   cyclophilin A at 278 K, 9 datasets (5KUZ..5KV7); structure factors      data: examples/data/ (downloaded)
    cypa100   cyclophilin A at 100 K, 9 datasets (5KUL..5KUW); structure factors      data: examples/data/ (downloaded)
Outputs go to assets/examples/<name>[_quick]/media/videos/... (git-ignored).
"""
import argparse
import subprocess
import sys
from pathlib import Path
from urllib.request import urlretrieve

sys.path.insert(0, str(Path(__file__).resolve().parent))
from density_prepare import SERIES_100K, SERIES_278K  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
EX = ROOT / "examples"
OUT = ROOT / "assets" / "examples"
CYPA_ACTIVE_SITE = "55,60,61,63,72,101,102,103,111,113,121,122,126"

QUICK_TRAJ = ["--width", "640", "--height", "360", "--antialias", "1", "--between", "8",
              "--manim-args", "-r 640,360 --fps 15"]
QUICK_DENSITY = ["--size", "400", "--antialias", "1", "--frames-per-set", "2",
                 "--manim-args", "-r 640,360 --fps 15"]


def download(url: str, dest: Path):
    if dest.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url}")
    urlretrieve(url, dest)


def pdb_url(pdb_id):
    return f"https://files.rcsb.org/download/{pdb_id}.pdb"


# ---- data fetchers (density datasets are fetched by render_density.py itself) ----
def fetch_adk():
    for i in ("4AKE", "1AKE"):
        download(pdb_url(i), EX / "adenylate_kinase" / f"{i}.pdb")
    return EX / "adenylate_kinase"


def fetch_nmr():
    download(pdb_url("1D3Z"), EX / "data" / "nmr_ubiquitin" / "1D3Z.pdb")
    return EX / "data" / "nmr_ubiquitin"


def fetch_cypa(ids, folder):
    for i in ids:
        download(f"https://files.rcsb.org/download/{i}.cif", EX / "data" / folder / f"{i}.cif")
        download(f"https://files.rcsb.org/download/{i}-sf.cif", EX / "data" / folder / f"{i}-sf.cif")
    return EX / "data" / folder


def command(name, quick):
    py, scripts = sys.executable, ROOT / "scripts"
    out = OUT / (name + ("_quick" if quick else ""))
    if name == "adk":
        d = fetch_adk()
        return [py, scripts / "render_traj.py", d / "4AKE.pdb", d / "1AKE.pdb",
                "--labels", "Open (4AKE),Closed (1AKE)", "--chains", "A",
                "--title", "Adenylate kinase: lid closing", "--out", out] + (QUICK_TRAJ if quick else [])
    if name == "nmr":
        d = fetch_nmr()
        return [py, scripts / "render_traj.py", d / "1D3Z.pdb", "--between", "5",
                "--title", "Ubiquitin: NMR ensemble (PDB 1D3Z)", "--out", out] + (QUICK_TRAJ if quick else [])
    if name in ("cypa278", "cypa100"):
        ids, folder, temp = (SERIES_278K, "cypa_278K", "278 K") if name == "cypa278" else (SERIES_100K, "cypa_100K", "100 K")
        d = fetch_cypa(ids, folder)
        return [py, scripts / "render_density.py", "--dir", d, "--ids", ",".join(ids), "--out", out,
                "--focus", CYPA_ACTIVE_SITE, "--zoom-label", "Active site",
                "--title", f"Human cyclophilin A, {temp}: data sets 1 to 9"] + (QUICK_DENSITY if quick else [])
    raise SystemExit(f"Unknown example '{name}'. Run: python scripts/examples.py list")


NAMES = ["adk", "nmr", "cypa278", "cypa100"]
DESC = {"adk": "adenylate kinase open -> closed (coordinates change)",
        "nmr": "ubiquitin NMR ensemble, 10 models (coordinates change)",
        "cypa278": "cyclophilin A 278 K, 9 datasets (electron density change)",
        "cypa100": "cyclophilin A 100 K, 9 datasets (electron density change)"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    r = sub.add_parser("run")
    r.add_argument("name", choices=NAMES)
    r.add_argument("--quick", action="store_true", help="small/fast render (about a minute) to try the pipeline")
    f = sub.add_parser("fetch")
    f.add_argument("name", choices=NAMES + ["all"])
    a, extra = ap.parse_known_args()   # unknown options are passed through to the render script
    if extra and a.cmd != "run":
        ap.error(f"unrecognized arguments: {' '.join(extra)}")

    if a.cmd == "list":
        for n in NAMES:
            print(f"  {n:8s} {DESC[n]}")
    elif a.cmd == "fetch":
        for n in (NAMES if a.name == "all" else [a.name]):
            {"adk": fetch_adk, "nmr": fetch_nmr,
             "cypa278": lambda: fetch_cypa(SERIES_278K, "cypa_278K"),
             "cypa100": lambda: fetch_cypa(SERIES_100K, "cypa_100K")}[n]()
        print("Data ready.")
    else:
        cmd = [str(c) for c in command(a.name, a.quick)] + extra
        print("Running:", " ".join(cmd[1:3]), "...")
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
