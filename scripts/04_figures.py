"""
Generate poster figures from trained latents + analysis CSVs.

Outputs (all PDF, 300 DPI):
  figures/fig1_variance_decomposition.pdf
  figures/fig2_umap_zbg_by_dataset.pdf
  figures/fig3_umap_zshared_by_drug.pdf
  figures/fig4_top_genes_shared_axis.pdf
  figures/fig5_drug_specific_genes.pdf

Usage:
    python scripts/04_figures.py --run runs/<timestamp>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
        "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    return plt


def fig_variance_decomposition(csv_path: Path, out: Path):
    plt = _setup_matplotlib()
    df = pd.read_csv(csv_path)
    fig, ax = plt.subplots(figsize=(4, 3))
    colors = ["#E15759", "#4E79A7"]
    bars = ax.bar(df["component"], df["fraction"] * 100, color=colors, edgecolor="white")
    for bar, frac in zip(bars, df["fraction"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{frac*100:.1f}%", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("% of salient variance")
    ax.set_title("Shared vs. drug-specific toxicity variance")
    ax.set_ylim(0, 110)
    ax.set_xticklabels(["Shared\ntoxicity", "Drug-\nspecific"])
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out)
    plt.close(fig)
    print(f"  [fig] {out.name}")


def fig_top_genes(csv_path: Path, out: Path, title: str, n: int = 20):
    plt = _setup_matplotlib()
    df = pd.read_csv(csv_path)
    # Remove pseudogenes (Gm*) — they appear as noise at the top of attribution lists
    df = df[~df["gene"].str.match(r"^Gm\d+$")].head(n).sort_values("attribution")
    colors = ["#E15759" if v > 0 else "#4E79A7" for v in df["attribution"]]
    fig, ax = plt.subplots(figsize=(5, 0.35 * n + 1))
    ax.barh(df["gene"], df["attribution"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Integrated gradient attribution")
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out)
    plt.close(fig)
    print(f"  [fig] {out.name}")


def fig_drug_specific_genes(csv_path: Path, out: Path, n: int = 15):
    plt = _setup_matplotlib()
    df = pd.read_csv(csv_path)
    drugs = df["drug"].unique()
    fig, axes = plt.subplots(1, len(drugs), figsize=(5 * len(drugs), 0.35 * n + 1),
                             sharey=False)
    if len(drugs) == 1:
        axes = [axes]
    for ax, drug in zip(axes, drugs):
        sub = df[df["drug"] == drug].head(n).sort_values("attribution")
        colors = ["#E15759" if v > 0 else "#4E79A7" for v in sub["attribution"]]
        ax.barh(sub["gene"], sub["attribution"], color=colors)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title(f"{drug}\n(drug-specific axis)")
        ax.set_xlabel("Attribution")
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out)
    plt.close(fig)
    print(f"  [fig] {out.name}")


def fig_umap(latents_path: Path, embedding_key: str, color_key: str,
             title: str, out: Path, palette: dict | None = None):
    import anndata as ad
    import scanpy as sc
    plt = _setup_matplotlib()

    adata = ad.read_h5ad(latents_path)
    adata.obsm["X_latent"] = adata.obsm[embedding_key]

    print(f"  [umap] computing neighbors + UMAP for {embedding_key} ...")
    sc.pp.neighbors(adata, use_rep="X_latent", n_neighbors=15, n_pcs=None)
    sc.tl.umap(adata, min_dist=0.3)

    fig, ax = plt.subplots(figsize=(5, 4))
    categories = adata.obs[color_key].astype(str).unique()
    if palette is None:
        import matplotlib.cm as cm
        cmap = cm.get_cmap("tab10", len(categories))
        palette = {c: cmap(i) for i, c in enumerate(sorted(categories))}

    for cat in sorted(categories):
        mask = adata.obs[color_key].astype(str) == cat
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=1, alpha=0.4, label=cat, color=palette.get(cat, "gray"), rasterized=True)

    ax.legend(markerscale=5, frameon=False, loc="upper left",
              bbox_to_anchor=(1, 1), fontsize=9)
    ax.set_title(title)
    ax.axis("off")
    fig.savefig(out)
    plt.close(fig)
    print(f"  [fig] {out.name}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True)
    args = ap.parse_args()

    run_dir = args.run
    analysis_dir = run_dir / "analysis"
    fig_dir = run_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    latents = run_dir / "latents.h5ad"

    # Fig 1: variance decomposition bar
    fig_variance_decomposition(
        analysis_dir / "variance_decomposition.csv",
        fig_dir / "fig1_variance_decomposition.pdf",
    )

    # Fig 2a: UMAP of z_bg colored by Leiden cluster (cell-type proxy)
    import anndata as ad, scanpy as sc
    adata_tmp = ad.read_h5ad(latents)
    adata_tmp.obsm["X_latent"] = adata_tmp.obsm["z_bg"]
    sc.pp.neighbors(adata_tmp, use_rep="X_latent", n_neighbors=15, n_pcs=None)
    sc.tl.leiden(adata_tmp, resolution=0.5, key_added="leiden_bg")
    sc.tl.umap(adata_tmp, min_dist=0.3)
    # save UMAP coords back so fig_umap reuses them
    adata_tmp.write_h5ad(latents)

    fig_umap(latents, "z_bg", "leiden_bg",
             "Background latent (z_bg) — Leiden clusters (cell-type proxy)",
             fig_dir / "fig2a_umap_zbg_leiden.pdf")

    # Fig 2b: same UMAP colored by dataset (batch check)
    drug_palette = {
        "control": "#999999",
        "doxorubicin": "#E15759",
        "cisplatin": "#4E79A7",
        "fludarabine_cyclophosphamide": "#F28E2B",
    }
    fig_umap(latents, "z_bg", "dataset_id",
             "Background latent (z_bg) — colored by dataset",
             fig_dir / "fig2_umap_zbg_by_dataset.pdf")

    # Fig 3: UMAP of z_shared colored by drug
    fig_umap(latents, "z_shared", "drug",
             "Shared-toxicity latent (z_shared) — colored by drug",
             fig_dir / "fig3_umap_zshared_by_drug.pdf",
             palette=drug_palette)

    # Fig 4: top genes on shared axis
    fig_top_genes(
        analysis_dir / "top_genes_shared_axis.csv",
        fig_dir / "fig4_top_genes_shared_axis.pdf",
        title="Top genes: shared-toxicity axis (IG attribution)",
    )

    # Fig 5: drug-specific genes (side by side)
    drug_csv = analysis_dir / "top_genes_per_drug.csv"
    if drug_csv.exists() and drug_csv.stat().st_size > 0:
        fig_drug_specific_genes(drug_csv, fig_dir / "fig5_drug_specific_genes.pdf")

    print(f"\n[done] figures saved to {fig_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
