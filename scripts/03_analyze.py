"""
Interpret the learned salient latents by attributing gene contributions to
each axis, comparing drugs, and highlighting cell-type-specific vulnerability.

Produces:
  - out/top_genes_shared_axis.csv
  - out/top_genes_per_drug.csv
  - out/celltype_shared_loadings.csv
  - out/figures/*.pdf (UMAPs of z_bg vs z_shared vs z_drug, per-cell-type heatmaps)

Usage:
    python scripts/03_analyze.py --run runs/<timestamp>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.models.mc_contrastive_vi import MCContrastiveVI, MCContrastiveVIConfig  # noqa: E402


# ---------------------------------------------------------------- attribution

def integrated_gradients_to_axis(
    model: MCContrastiveVI,
    x: torch.Tensor,
    drug_idx: torch.Tensor,
    axis: str,          # "shared" or "drug_k"
    drug_k: int = 1,    # only used when axis starts with "drug"
    n_steps: int = 32,
    device: str = "cpu",
) -> np.ndarray:
    """
    Integrated gradients of the mean of the salient latent axis wrt input genes.
    Returns (G,) attribution vector averaged over cells in `x`.

    For "shared": attributes to the shared-toxicity encoder mean.
    For "drug_k": attributes to the k-th drug-specific encoder mean.
    """
    model.eval()
    x = x.to(device)
    baseline = torch.zeros_like(x)  # zero baseline; swap to median if preferred
    attrs = torch.zeros_like(x)
    for alpha in torch.linspace(0, 1, n_steps):
        xi = baseline + alpha * (x - baseline)
        xi.requires_grad_(True)
        enc = model.encode(xi, drug_idx.to(device))
        if axis == "shared":
            scalar = enc["mu_shared"].sum()
        elif axis.startswith("drug_"):
            k = int(axis.split("_")[1])
            scalar = enc["mu_drug_list"][k].sum()
        else:
            raise ValueError(f"Unknown axis: {axis}")
        grad = torch.autograd.grad(scalar, xi)[0]
        attrs = attrs + grad.detach()
    ig = (x - baseline) * attrs / n_steps
    return ig.mean(dim=0).cpu().numpy()


def top_genes(attr: np.ndarray, gene_names, n: int = 50) -> pd.DataFrame:
    order = np.argsort(-np.abs(attr))
    return pd.DataFrame({
        "gene": [gene_names[i] for i in order[:n]],
        "attribution": attr[order[:n]],
    })


# ---------------------------------------------------------------- decomposition

def shared_vs_drug_variance(latents_h5ad_path: Path) -> pd.DataFrame:
    """
    How much of each cell's treated-state latent is captured by the shared vs
    drug-specific axis? Useful for the "convergent vs idiosyncratic" claim.
    """
    import anndata as ad
    adata = ad.read_h5ad(latents_h5ad_path)
    z_sh = adata.obsm["z_shared"]
    z_dr = adata.obsm["z_drug"]
    treated = (adata.obs["drug"].astype(str) != "control").values
    if treated.sum() == 0:
        return pd.DataFrame()
    v_sh = np.var(z_sh[treated], axis=0).sum()
    v_dr = np.var(z_dr[treated], axis=0).sum()
    return pd.DataFrame({
        "component": ["shared_toxicity", "drug_specific"],
        "total_variance": [v_sh, v_dr],
        "fraction": [v_sh / (v_sh + v_dr), v_dr / (v_sh + v_dr)],
    })


def celltype_shared_loadings(latents_h5ad_path: Path) -> pd.DataFrame:
    """
    Per-cell-type mean magnitude along the shared-toxicity axis.
    High value = that cell type is strongly engaged by the convergent signature.
    """
    import anndata as ad
    adata = ad.read_h5ad(latents_h5ad_path)
    treated = (adata.obs["drug"].astype(str) != "control").values
    cts = adata.obs["cell_type"].astype(str).values
    z_sh = adata.obsm["z_shared"]
    rows = []
    for ct in np.unique(cts):
        mask = treated & (cts == ct)
        if mask.sum() < 10:
            continue
        mag = np.linalg.norm(z_sh[mask], axis=1).mean()
        rows.append({"cell_type": ct, "n_cells": int(mask.sum()),
                     "mean_shared_magnitude": float(mag)})
    return pd.DataFrame(rows).sort_values("mean_shared_magnitude", ascending=False)


# ---------------------------------------------------------------- driver

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--n-top", type=int, default=50)
    args = ap.parse_args()

    run_dir = args.run
    outdir = args.outdir or (run_dir / "analysis")
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. variance decomposition
    latents_path = run_dir / "latents.h5ad"
    var_df = shared_vs_drug_variance(latents_path)
    var_df.to_csv(outdir / "variance_decomposition.csv", index=False)
    print(var_df)

    # 2. per-cell-type engagement of shared axis
    ct_df = celltype_shared_loadings(latents_path)
    ct_df.to_csv(outdir / "celltype_shared_loadings.csv", index=False)
    print(ct_df.head(15))

    # 3. gene attribution (load model + compute IG)
    import anndata as ad
    import scipy.sparse as sp
    adata = ad.read_h5ad(latents_path)
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()

    ckpt = torch.load(run_dir / "model.pt", map_location="cpu", weights_only=False)
    cfg = MCContrastiveVIConfig(**ckpt["model_cfg"])
    model = MCContrastiveVI(cfg)
    model.load_state_dict(ckpt["model_state"])

    drug_to_idx = ckpt["drug_to_idx"]
    # use treated cells only for shared-axis attribution
    drugs = adata.obs["drug"].astype(str).values
    treated_mask = drugs != "control"
    if treated_mask.sum() > 2000:
        # random subsample for speed
        idx = np.random.default_rng(0).choice(np.where(treated_mask)[0], 2000, replace=False)
    else:
        idx = np.where(treated_mask)[0]

    x_t = torch.tensor(X[idx], dtype=torch.float32)
    d_t = torch.tensor([drug_to_idx[d] for d in drugs[idx]], dtype=torch.long)

    print("[attribute] shared-toxicity axis ...")
    attr_sh = integrated_gradients_to_axis(model, x_t, d_t, axis="shared")
    top_sh = top_genes(attr_sh, list(adata.var_names), n=args.n_top)
    top_sh.to_csv(outdir / "top_genes_shared_axis.csv", index=False)
    print(top_sh.head(15))

    # per-drug
    rows = []
    for drug, k in drug_to_idx.items():
        if k == 0:
            continue
        print(f"[attribute] drug-specific axis: {drug}")
        mask = drugs == drug
        if mask.sum() < 50:
            continue
        sub = np.where(mask)[0]
        if len(sub) > 1000:
            sub = np.random.default_rng(0).choice(sub, 1000, replace=False)
        x_d = torch.tensor(X[sub], dtype=torch.float32)
        d_d = torch.tensor([k] * len(sub), dtype=torch.long)
        attr = integrated_gradients_to_axis(model, x_d, d_d, axis=f"drug_{k-1}")
        top = top_genes(attr, list(adata.var_names), n=args.n_top)
        top["drug"] = drug
        rows.append(top)
    if rows:
        pd.concat(rows).to_csv(outdir / "top_genes_per_drug.csv", index=False)

    print(f"\n[done] wrote analysis artifacts to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
