"""
Preprocess raw per-dataset AnnData files into a single merged AnnData ready
for contrastive modeling.

Steps per dataset:
  1. Load raw counts (matrix market or h5ad)
  2. QC filter (min genes/cell, max %mito, min cells/gene)
  3. Basic doublet removal (optional; Scrublet)
  4. Annotate obs with: drug, dataset_id, region, assay, sex (from registry)

After all datasets are per-processed:
  5. Concatenate into a single AnnData
  6. Normalize total + log1p
  7. Select HVGs (batch-aware: flavor="seurat_v3", batch_key="dataset_id")
  8. Cell-type label transfer from an Allen mouse brain reference (scANVI / symphony-like)
  9. Save merged AnnData to data/processed/merged_chemobrain.h5ad

Notes
-----
Cross-study integration is notoriously sensitive. We keep the integration
step SEPARATE from the contrastive model because the model itself learns a
batch-aware decoder. Here we only harmonize cell-type labels, not expression.

If cell-type labels are unavailable for some datasets, run scanpy.tl.ingest
from the Allen reference as a fallback (implemented below).

Usage:
    python scripts/01_preprocess.py --registry configs/dataset_registry.yaml \\
        --raw-dir data/raw --out data/processed/merged_chemobrain.h5ad
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import yaml


# Import scanpy lazily so this module is importable without it (CI, tests)
def _imports():
    import anndata as ad
    import numpy as np
    import scanpy as sc

    return ad, np, sc


def _read_prefixed_10x_dir(mtx_dir: Path):
    """
    Read a 10x directory where files have a sample-specific prefix, e.g.:
      GSM8367861_CNT_barcodes.tsv.gz
      GSM8367861_CNT_features.tsv.gz
      GSM8367861_CNT_matrix.mtx.gz
    scanpy.read_10x_mtx requires exact filenames; we rename into a temp dir.
    """
    import shutil
    import tempfile

    ad_mod, np, sc = _imports()

    matrix = next(mtx_dir.glob("*matrix.mtx*"), None)
    barcodes = next(mtx_dir.glob("*barcodes.tsv*"), None)
    features = next((p for p in mtx_dir.glob("*features.tsv*") or
                     mtx_dir.glob("*genes.tsv*")), None)
    if features is None:
        features = next(mtx_dir.glob("*genes.tsv*"), None)

    if matrix is None or barcodes is None or features is None:
        return None

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        shutil.copy(matrix,   tmp / ("matrix.mtx.gz"   if matrix.suffix == ".gz"   else "matrix.mtx"))
        shutil.copy(barcodes, tmp / ("barcodes.tsv.gz"  if barcodes.suffix == ".gz" else "barcodes.tsv"))
        shutil.copy(features, tmp / ("features.tsv.gz"  if features.suffix == ".gz" else "features.tsv"))
        return sc.read_10x_mtx(str(tmp), var_names="gene_symbols", make_unique=True)


def _condition_from_dirname(dirname: str) -> str:
    """Infer condition label from a GEO Supp_* subdirectory name."""
    d = dirname.lower()
    if "control" in d or "_cnt" in d or "vehicle" in d or "ctrl" in d:
        return "control"
    # Check rescue arms before plain drug names
    if "acy" in d or "hdac6i" in d or "rescue" in d:
        return "rescue"
    if "dox" in d:
        return "doxorubicin"
    if "cisplat" in d or "cispl" in d:
        return "cisplatin"
    if "pbs" in d:
        return "control"
    if "flucy" in d or "flu" in d:
        return "fludarabine_cyclophosphamide"
    return "unknown"


def _read_supp_dirs(ds_dir: Path, drug: str):
    """
    GEO often deposits per-sample data in Supp_GSM* subdirectories with
    prefixed filenames.  Read each, label with condition, and concatenate.
    Returns None if no Supp_* dirs found.
    """
    ad_mod, np, sc = _imports()

    supp_dirs = sorted(ds_dir.glob("Supp_*"))
    if not supp_dirs:
        return None

    adatas = []
    for d in supp_dirs:
        adata = _read_prefixed_10x_dir(d)
        if adata is None:
            continue
        condition = _condition_from_dirname(d.name)
        if condition == "unknown":
            raise ValueError(f"Unmapped sample condition {d.name}")
        adata.obs["rescue"] = condition == "rescue"
        adata.obs["sample_id"] = d.name
        adata.obs["condition"] = condition
        # map to drug label used by the model
        adata.obs["drug"] = "control" if condition == "control" else drug
        adatas.append(adata)

    if not adatas:
        return None

    import anndata as ad
    return ad.concat(adatas, join="outer", index_unique="-")


def load_dataset(raw_dir: Path, ds_entry: dict):
    """
    Load one dataset's raw counts into an AnnData.
    Supports (in priority order):
      1. h5ad / h5ad.gz
      2. 10x filtered .h5
      3. Standard 10x mtx directory (matrix.mtx[.gz])
      4. GEO Supp_* dirs with sample-prefixed mtx files
    """
    ad_mod, np, sc = _imports()
    ds_id = ds_entry["id"]
    drug = str(ds_entry.get("drug", "unknown"))
    ds_dir = raw_dir / ds_id
    if not ds_dir.exists():
        raise FileNotFoundError(f"No raw data directory for {ds_id}: {ds_dir}")

    adata = None

    # 1. h5ad (plain or gzipped — h5py can't open .gz directly so decompress first)
    h5ads = sorted(ds_dir.glob("**/*.h5ad")) + sorted(ds_dir.glob("**/*.h5ad.gz"))
    if h5ads:
        src = h5ads[0]
        if src.suffix == ".gz":
            import gzip, shutil
            dest = src.with_suffix("")  # strip .gz -> .h5ad
            if not dest.exists():
                print(f"    [decompress] {src.name} ...")
                with gzip.open(src, "rb") as f_in, open(dest, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            src = dest
        adata = sc.read_h5ad(src)

    # 2. 10x filtered h5
    if adata is None:
        h5s = list(ds_dir.glob("**/*filtered*feature_bc_matrix.h5"))
        if h5s:
            adata = sc.read_10x_h5(str(h5s[0]))

    # 3. Standard 10x mtx (file named exactly matrix.mtx[.gz])
    if adata is None:
        mtx_dirs = [p.parent for p in ds_dir.glob("**/matrix.mtx*")
                    if p.name.startswith("matrix")]
        if mtx_dirs:
            adata = sc.read_10x_mtx(str(mtx_dirs[0]), var_names="gene_symbols", make_unique=True)

    # 4. GEO per-sample Supp_* dirs with prefixed filenames
    if adata is None:
        adata = _read_supp_dirs(ds_dir, drug)

    if adata is None:
        raise FileNotFoundError(
            f"No recognized count matrix under {ds_dir}. "
            "Expect .h5ad, 10x .h5, or 10x mtx directory."
        )

    # Apply condition_map if the registry entry specifies one (e.g. cisplatin h5ad)
    condition_col = ds_entry.get("condition_col")
    condition_map = ds_entry.get("condition_map")  # dict: drug_label -> [obs_values]
    if condition_col and condition_map and condition_col in adata.obs.columns:
        adata.obs["rescue"] = adata.obs[condition_col].astype(str).isin(ds_entry.get("rescue_values", []))
        if "biosample" in adata.obs:
            adata.obs["sample_id"] = adata.obs["biosample"].astype(str)
        value_to_drug = {}
        for drug_label, values in condition_map.items():
            for v in values:
                value_to_drug[str(v)] = drug_label
        adata.obs["drug"] = (
            adata.obs[condition_col].astype(str).map(value_to_drug).fillna("unknown")
        )
        unknown_mask = adata.obs["drug"] == "unknown"
        if unknown_mask.any():
            print(f"    [warn] {unknown_mask.sum()} cells had unmapped {condition_col} values; dropping")
            adata = adata[~unknown_mask].copy()

    # Annotate obs with dataset metadata (don't overwrite per-cell condition set above)
    for key in ("id", "drug", "region", "assay", "species"):
        if key in ds_entry and ds_entry[key] is not None:
            col = "dataset_id" if key == "id" else key
            if col not in adata.obs.columns:
                adata.obs[col] = str(ds_entry[key])
    adata.obs["dataset_id"] = ds_id

    return adata


def qc_filter(adata, min_genes=500, min_cells=10, max_pct_mt=10):
    _, np, sc = _imports()
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-", "Mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()
    return adata


def concat_and_normalize(per_ds_adatas: list, n_hvgs: int = 3000):
    ad, np, sc = _imports()
    adata = ad.concat(per_ds_adatas, join="inner", label="batch_concat", index_unique="-")
    # preserve raw
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # batch-aware HVG selection
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_hvgs, flavor="seurat_v3",
        batch_key="dataset_id", subset=True, layer="counts",
    )
    return adata


def annotate_cell_types(adata, reference_h5ad: Path | None):
    """
    Optional: transfer cell-type labels from an Allen reference using scanpy.tl.ingest.
    If no reference provided, leave 'cell_type' column set to 'unannotated'.
    """
    ad_mod, np, sc = _imports()
    if reference_h5ad is None or not reference_h5ad.exists():
        adata.obs["cell_type"] = "unannotated"
        print("[warn] No reference provided; cell_type = 'unannotated'. "
              "Run Azimuth/scANVI separately and write back the labels.")
        return adata

    ref = sc.read_h5ad(reference_h5ad)
    # basic alignment
    common = adata.var_names.intersection(ref.var_names)
    query = adata[:, common].copy()
    sc.tl.ingest(query, ref[:, common].copy(), obs="cell_type")
    adata.obs["cell_type"] = query.obs["cell_type"].reindex(adata.obs_names)
    return adata


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", type=Path, default="configs/dataset_registry.yaml")
    ap.add_argument("--raw-dir", type=Path, default="data/raw")
    ap.add_argument("--out", type=Path, default="data/processed/merged_chemobrain.h5ad")
    ap.add_argument("--n-hvgs", type=int, default=3000)
    ap.add_argument("--reference", type=Path, default=None,
                    help="Allen mouse brain reference h5ad for label transfer.")
    ap.add_argument("--include-rescue", action="store_true", help="Explicitly retain rescue arms in preprocessing output")
    args = ap.parse_args()

    with open(args.registry) as f:
        reg = yaml.safe_load(f)

    per_ds = []
    for ds in reg["datasets"]:
        if ds.get("status") != "confirmed":
            continue
        if ds["id"] == "allen_mousebrain_ctx_hpf":
            continue  # reference only, not a condition
        print(f"[load] {ds['id']}")
        try:
            a = load_dataset(args.raw_dir, ds)
        except FileNotFoundError as e:
            print(f"  [skip] {e}")
            continue
        if "rescue" in a.obs and not args.include_rescue:
            a = a[~a.obs["rescue"]].copy()
        a = qc_filter(a)
        per_ds.append(a)
        print(f"  -> {a.n_obs} cells, {a.n_vars} genes after QC")

    if not per_ds:
        print("No datasets loaded. Confirm accessions in the registry.")
        return 1

    print(f"\n[merge] {len(per_ds)} datasets")
    adata = concat_and_normalize(per_ds, n_hvgs=args.n_hvgs)
    print(f"  merged: {adata.n_obs} cells, {adata.n_vars} HVGs")

    adata = annotate_cell_types(adata, args.reference)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    print(f"\n[done] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
