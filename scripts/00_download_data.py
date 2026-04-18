"""
Download raw scRNA-seq / snRNA-seq count matrices from GEO/SRA based on the
dataset registry (configs/dataset_registry.yaml).

Strategy:
  - For each dataset with `status: confirmed`, fetch the processed expression
    matrix if a 10x-style or h5ad supplementary file is available on GEO.
  - GEOparse fetches per-sample (GSM) supplementary files; this script also
    parses the downloaded soft.gz to pull series-level (GSE) supplementary
    files that GEOparse misses (e.g. GSE216146_chemo_brain.h5ad.gz).
  - Skip datasets marked `needs_verification` or `candidate`; print a TODO.

Usage:
    python scripts/00_download_data.py --registry configs/dataset_registry.yaml
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import requests
import yaml


def load_registry(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def _geo_stub(accession: str) -> str:
    """GSE216146 -> GSE216nnn  (GEO FTP directory convention)."""
    return accession[:-3] + "nnn" if len(accession) > 6 else accession + "nnn"


def _download_file(url: str, dest: Path) -> bool:
    """Stream-download url to dest. Returns True on success."""
    if dest.exists():
        print(f"    [skip] already have {dest.name}")
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        r = requests.get(url, stream=True, timeout=120)
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = 100 * downloaded // total
                    print(f"\r    {dest.name}  {pct}%", end="", flush=True)
        print()
        return True
    except Exception as e:
        print(f"\n    [FAIL] {url}: {e}")
        dest.unlink(missing_ok=True)
        return False


def _series_suppl_urls_from_soft(soft_gz: Path) -> list[str]:
    """
    Parse a GEO soft.gz and return URLs on !Series_supplementary_file lines.
    These are the series-level files that GEOparse.download_supplementary_files
    doesn't download (it only fetches per-sample GSM files).
    """
    urls = []
    try:
        with gzip.open(soft_gz, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("!Series_supplementary_file"):
                    parts = line.split("=", 1)
                    if len(parts) == 2:
                        url = parts[1].strip()
                        if url.startswith("ftp://"):
                            # convert ftp:// -> https:// for requests
                            url = "https://" + url[len("ftp://"):]
                        if url:
                            urls.append(url)
    except Exception as e:
        print(f"    [warn] could not parse {soft_gz.name}: {e}")
    return urls


def download_geo(accession: str, outdir: Path) -> None:
    """
    Download supplementary files for a GEO Series accession.

    1. GEOparse: fetches the soft file + per-sample (GSM) supplementary files.
    2. Soft parser: finds series-level (GSE) supplementary file URLs embedded
       in the soft file and downloads any we don't have yet.
    """
    try:
        import GEOparse  # type: ignore
    except ImportError:
        print(f"  [WARN] GEOparse not installed; skipping {accession}")
        return

    outdir.mkdir(parents=True, exist_ok=True)
    gse = GEOparse.get_GEO(geo=accession, destdir=str(outdir), silent=True)
    gse.download_supplementary_files(directory=str(outdir))

    # Also grab series-level supplementary files (GEOparse skips these)
    soft_files = list(outdir.glob("*.soft.gz"))
    for soft in soft_files:
        urls = _series_suppl_urls_from_soft(soft)
        for url in urls:
            fname = url.split("/")[-1]
            dest = outdir / fname
            if not dest.exists():
                print(f"    [series suppl] {fname}")
                _download_file(url, dest)

    print(f"  [OK]  {accession} -> {outdir}")


def _print_allen_instructions(ds_id: str, accession: str, outdir: Path) -> None:
    print(f"[manual] {ds_id} ({accession}) -- not a GEO accession, needs manual download.")
    print(f"  Allen Mouse Brain Cell Atlas (Yao et al. 2021):")
    print(f"  https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq")
    print(f"  Place the h5ad at: {outdir}/allen_ctx_hpf.h5ad")
    print(f"  Or skip for now: set cell_type_key: leiden in configs/default.yaml")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", type=Path, default="configs/dataset_registry.yaml")
    ap.add_argument("--out", type=Path, default="data/raw")
    ap.add_argument("--only", type=str, default=None,
                    help="Comma-separated dataset IDs to download (default: all confirmed)")
    args = ap.parse_args()

    reg = load_registry(args.registry)
    to_get = set(args.only.split(",")) if args.only else None

    for ds in reg["datasets"]:
        ds_id = ds["id"]
        status = ds.get("status", "candidate")
        if to_get is not None:
            if ds_id not in to_get:
                continue
        else:
            if status != "confirmed":
                print(f"[skip] {ds_id}: status={status}  (accession={ds.get('accession')})")
                continue

        accession = ds["accession"]
        if not accession or accession.startswith("TBD"):
            print(f"[skip] {ds_id}: no accession set")
            continue

        outdir = args.out / ds_id
        source = ds.get("source", "geo")

        if source == "allen_brain_map" or not accession.startswith("GSE"):
            _print_allen_instructions(ds_id, accession, outdir)
            continue

        print(f"[get ] {ds_id} ({accession}) ...")
        download_geo(accession, outdir)

    print("\nDone. Next: python scripts/01_preprocess.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
