# Methodology

Detailed methods for the chemo-brain contrastive analysis project. This is the source document for the **Methods** section of the poster/paper.

## 1. Data

### 1.1 Source datasets

We aggregate publicly deposited mouse brain single-cell or single-nucleus RNA-seq datasets in which a chemotherapeutic agent is a primary experimental variable and where a matched vehicle control is available. Candidate datasets and their GEO accessions are tracked in `configs/dataset_registry.yaml`. Datasets marked `status: confirmed` have been independently verified against the source publication's Data Availability statement.

### 1.2 Inclusion / exclusion criteria

- **Include**: mouse brain tissue, scRNA-seq or snRNA-seq, a chemotherapeutic agent or chemotherapy-like immunotherapy (e.g. CAR-T) condition with a matched vehicle or mock arm, and a count matrix publicly available on GEO or a comparable repository.
- **Exclude**: bulk RNA-seq-only studies (used only for pathway-level consistency checks, not model input), cell-line studies, human-only studies at model-training stage (reserved for future projection-based validation).

### 1.3 QC

Per-cell QC is applied per dataset before merging, using scanpy:

- Minimum 500 detected genes per cell
- Minimum 10 cells per gene
- Maximum 10% mitochondrial read content
- Doublet removal via Scrublet (optional; enabled when dataset size permits)

### 1.4 Normalization and feature selection

Counts are normalized per cell to 10,000 reads and log1p-transformed. Highly variable genes are selected with `scanpy.pp.highly_variable_genes` using `flavor="seurat_v3"` and `batch_key="dataset_id"` to obtain a consistent HVG set across studies. The default target count is 3,000 HVGs.

### 1.5 Cross-study integration and cell-type annotation

Cell-type labels are transferred from an Allen Mouse Brain reference (cortex + hippocampus taxonomy, Yao et al. 2021) using `scanpy.tl.ingest` or, for more robust labels, scANVI / Azimuth. Batch effects at the expression level are handled jointly by the contrastive model's batch-aware decoder (a one-hot encoding of `dataset_id` is concatenated to the latent code before decoding), which avoids removing biological variation we care about (the drug signal) during a preprocessing integration step.

## 2. Model

### 2.1 Architecture

We introduce **Multi-Condition Contrastive Variational Inference (MC-ContrastiveVI)**, extending the contrastiveVI framework of Weinberger et al. to the multi-drug setting.

Let `x ∈ ℝᴳ` denote the log-normalized expression of `G` HVGs for one cell, `d ∈ {0, 1, …, K}` its drug label (0 = control), and `b` its dataset-of-origin index. The model learns three latent spaces:

- **Background** `z_bg` ∈ ℝ^{d_bg}`: captures biological variance shared across treated and control cells (cell type, state, sex, age, etc.).
- **Shared-toxicity** `z_shared` ∈ ℝ^{d_sh}`: captures variation enriched in any treated cell relative to controls.
- **Drug-specific** `z_drug[k]` ∈ ℝ^{d_dr}`, one per non-control drug: captures variation unique to drug `k`.

Each latent has its own encoder: a small MLP trunk (2 hidden layers × 256 units by default, LayerNorm + ReLU + Dropout) followed by mean and log-variance heads producing a diagonal Gaussian posterior, sampled via the reparameterization trick.

### 2.2 Hard gating

To ensure control cells occupy only `z_bg`, the salient latents are hard-gated by drug label at the decoder input:

- For control cells (`d = 0`): `z_shared_gated = 0`, `z_drug[k]_gated = 0` for all `k`.
- For drug-`k` cells (`d = k`, `k ≥ 1`): `z_shared_gated = z_shared`, `z_drug[k]_gated = z_drug[k]`, and `z_drug[j≠k]_gated = 0`.

This gating is the key structural difference from standard multi-treatment analyses and is what makes the learned axes *interpretable* as shared-vs-idiosyncratic toxicity rather than arbitrary treatment-covarying directions.

### 2.3 Decoder

The decoder concatenates `[z_bg, z_shared_gated, concat_k(z_drug[k]_gated), one_hot(b)]` and passes it through an MLP trunk with a linear head producing the reconstructed log-expression `x̂`. The likelihood is Gaussian with a learnable per-gene log-sigma; a ZINB head on raw counts is a drop-in replacement for production and can be enabled by changing `model.input_type` in the config.

### 2.4 Loss

The full objective combines reconstruction, KL, and independence terms:

```
L = recon_nll(x, x̂)
  + β_kl · [ KL(q(z_bg)‖𝒩(0,I))
           + 1{d>0} · KL(q(z_shared)‖𝒩(0,I))
           + Σ_k 1{d=k} · KL(q(z_drug[k])‖𝒩(0,I)) ]
  + λ_bs · HSIC(z_bg, z_shared)
  + λ_bd · HSIC(z_bg, z_drug_cat)
  + λ_sd · HSIC(z_shared, z_drug_cat)
```

The HSIC (Hilbert-Schmidt Independence Criterion) terms, computed on each mini-batch with an RBF kernel using the median heuristic for bandwidth, push the three latent groups toward statistical independence. The KL on salient latents is masked per-cell so only cells actually using that latent contribute — a design choice that follows the contrastiveVI formulation and prevents the salient latents from collapsing under KL pressure when their gate is zero.

Default weights (tune per dataset in `configs/default.yaml`): `β_kl=1`, `λ_bs=λ_bd=10`, `λ_sd=5`. The shared–drug independence is weaker because shared and drug-specific toxicity may share some mechanistic directions by construction.

### 2.5 Training

- Optimizer: AdamW, learning rate 1e-3, weight decay 1e-5.
- Batch size 512.
- Gradient clipping at L2 norm 5.0.
- 80/10/10 train/val/test split by cell (no cell crosses splits; biological replicates within a drug stay within a split to prevent info leak).
- Early stopping on validation loss with patience 20 epochs.
- Seed fixed for reproducibility.

## 3. Interpretation

### 3.1 Gene attribution via integrated gradients

For each salient axis (`z_shared` overall, and each `z_drug[k]`), we compute **integrated gradients** of the axis's mean encoder output with respect to input gene expression, averaged over a held-out set of treated cells. The top-|attribution| genes for each axis define its molecular signature. We intentionally prefer integrated gradients over simple differential expression here because it reads the *model's* internal coordinates, not marginal per-gene statistics.

### 3.2 Variance decomposition

For treated cells only, we compute the trace of the per-cell covariance of `z_shared_gated` versus `z_drug_cat`. The ratio `trace(cov z_shared) / [trace(cov z_shared) + trace(cov z_drug_cat)]` quantifies how much of the treatment-induced latent variance is *shared across drugs*. A large fraction supports the convergent toxicity hypothesis.

### 3.3 Cell-type engagement

For each cell type `c`, we compute the mean L2 norm of `z_shared_gated` across cells of that type in the treated pool. High-magnitude cell types are those most engaged by the convergent toxicity manifold and are therefore the primary candidates for cell-type-of-origin of the cognitive phenotype.

### 3.4 Rescue-arm validation

For datasets with a rescue intervention (HDAC6i with doxorubicin; 40 Hz gamma with cisplatin), we encode rescued cells through the trained model. If the rescue is biologically reversing the shared-toxicity program, `‖z_shared_gated‖` in rescued cells should shift down toward control magnitudes. This is reported as a single scalar per rescue dataset with a bootstrap confidence interval.

### 3.5 Pathway enrichment

Top-attribution genes per axis are tested against Reactome, KEGG, and MSigDB Hallmark collections using `gseapy` with an FDR threshold of 0.05.

## 4. Baselines and ablations

To support the claim that MC-ContrastiveVI adds value over existing tools, we report:

- **PCA on treated cells only** — tests whether a simple linear decomposition finds the shared axis.
- **scVI latent + post-hoc treatment classifier** — standard deep-generative baseline without disentanglement.
- **Vanilla contrastiveVI (pairwise)** — run separately for each drug-vs-control pair, then measure overlap between the learned salient axes. This is the closest existing method and is the one to beat.
- **MC-ContrastiveVI without HSIC penalties** — tests whether disentanglement regularization actually matters.
- **MC-ContrastiveVI without hard gating** — tests whether the gating (vs. soft weighting) is necessary for interpretability.

For each, we report: reconstruction nats, latent-space treatment-classification accuracy, and qualitative UMAP separation.

## 5. Reproducibility

- All random seeds set in `configs/default.yaml`.
- Exact package versions pinned in `requirements.txt`.
- Every run is saved to `runs/<timestamp>/` with a frozen copy of the config, the checkpoint, per-epoch logs, and the latent-annotated AnnData.
- Synthetic-data smoke test (`tests/test_smoke.py`) must pass before any real run.
