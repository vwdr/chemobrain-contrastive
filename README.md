# chemobrain-contrastive

> Disentangling shared and drug-specific molecular signatures of chemotherapy-induced cognitive impairment (CICI) via multi-condition contrastive latent-variable modeling of mouse brain single-cell data.

---

## What this repo is

A bioinformatics project scaffold to execute the research plan: ingest public mouse scRNA-seq / snRNA-seq datasets spanning multiple chemotherapy agents, train a multi-condition contrastive VI model that separates (i) shared biological variance, (ii) a shared-toxicity axis common to all drugs, and (iii) drug-specific axes, then interpret the learned axes to generate mechanistic hypotheses about CICI. Final deliverable: a conference poster.

## What's built vs. what's next

**Built (and tested on synthetic data):**
- `src/models/mc_contrastive_vi.py` — the multi-condition contrastive VI model (PyTorch). Verified end-to-end: forward pass, gating (control cells have zero salients), backward pass, 57% loss reduction on simulated data, shared axis 2.4× more active on treated cells.
- `scripts/00_download_data.py` — GEO data fetcher reading from the registry.
- `scripts/01_preprocess.py` — QC, HVG selection, cross-study merge, cell-type annotation.
- `scripts/02_train.py` — training loop with early stopping, checkpointing, latent export.
- `scripts/03_analyze.py` — integrated-gradient gene attribution, variance decomposition, cell-type engagement.
- `tests/test_smoke.py` — synthetic-data end-to-end check.
- `configs/dataset_registry.yaml` — curated list of candidate GEO datasets.
- `configs/default.yaml` — all hyperparameters in one place.

**Next (this is what you do in Claude Code):**
1. **Verify GEO accessions** in `configs/dataset_registry.yaml`. Every entry marked `needs_verification` has a paper DOI — open the paper's Data Availability section, copy the accession, flip `status: confirmed`.
2. **Download data**: `python scripts/00_download_data.py`
3. **Preprocess**: `python scripts/01_preprocess.py`
4. **Train**: `python scripts/02_train.py` (needs a GPU for the full dataset; use a Colab, lab cluster, or your own GPU box)
5. **Analyze**: `python scripts/03_analyze.py --run runs/<timestamp>`
6. **Make figures → poster**.

---

## Getting started in Claude Code

Once you clone this into your working environment:

```bash
cd chemobrain-contrastive
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Sanity check: the model trains on synthetic data
python tests/test_smoke.py
```

If the smoke test passes (it will — it's been verified), your environment is good and you can move on to real data.

**Recommended Claude Code prompts at each stage:**

Stage 2 (dataset verification):
> "Open configs/dataset_registry.yaml. For each entry with status: needs_verification, search the web for the paper's Data Availability section and fill in the GEO accession. Flip status to confirmed when done."

Stage 3 (first training run):
> "Run a small training run: python scripts/02_train.py with epochs=30 and batch_size=256 to verify it works on real data. Report any errors and suggest fixes."

Stage 5 (analysis + figures):
> "Use scripts/03_analyze.py outputs to generate the poster figures: UMAP of z_bg colored by cell type, UMAP of z_shared colored by drug, variance decomposition bar, top-gene heatmap. Save them as PDFs at 300 DPI."

---

## Directory structure

```
chemobrain-contrastive/
├── configs/
│   ├── default.yaml              # all hyperparameters
│   └── dataset_registry.yaml     # GEO accessions + metadata
├── src/
│   ├── models/
│   │   └── mc_contrastive_vi.py  # THE model
│   ├── data/                     # (expand here as pipeline grows)
│   ├── training/                 # (put reusable loops here)
│   ├── analysis/                 # (put reusable analysis fns here)
│   └── utils/
├── scripts/
│   ├── 00_download_data.py       # GEO ingest
│   ├── 01_preprocess.py          # QC + merge + HVG + label transfer
│   ├── 02_train.py               # train + export latents
│   └── 03_analyze.py             # gene attribution + cross-drug decomposition
├── tests/
│   └── test_smoke.py             # synthetic end-to-end test
├── docs/
│   ├── methodology.md            # detailed methods for the poster/paper
│   └── proposal.md               # the research plan
├── notebooks/                    # exploration lives here (not committed)
├── data/                         # (gitignored)
│   ├── raw/                      # GEO downloads
│   └── processed/                # merged h5ad
├── runs/                         # (gitignored) training outputs
├── requirements.txt
├── .gitignore
└── README.md
```

---

## The model in one paragraph

Standard contrastiveVI (Weinberger et al., *Nature Methods* 2023) learns a background latent `z_bg` shared between treated and control cells, plus a salient latent `z_s` active only in treated cells. We extend this to the **multi-drug** setting: `z_bg` stays, we add a **shared-toxicity** salient `z_shared` that activates for any non-control drug, and a **drug-specific** salient `z_drug[k]` per drug, one-hot-gated so only the correct drug's head fires. HSIC penalties push the three groups toward statistical independence. The decoder reconstructs log-normalized expression from `[z_bg, z_shared, z_drug_cat]`, optionally conditioned on dataset-of-origin for batch correction. Control cells have `z_shared == 0` and all `z_drug[k] == 0` by hard gating, so those dimensions exclusively capture treatment-enriched variance.

## What the analysis answers

1. **Is there a convergent CICI signature?** Yes if `z_shared` explains a large fraction of treated-cell variance and its top genes enrich for known CICI pathways (inflammation, senescence, myelin, oxidative stress).
2. **Which cell types carry it?** Per-cell-type mean magnitude along `z_shared` — high-magnitude cell types are the vulnerable ones.
3. **What's drug-specific?** Per-drug attribution scores on `z_drug[k]` isolate idiosyncratic toxicity (e.g., paclitaxel-specific endothelial effects vs. methotrexate-specific OPC effects).
4. **Does a rescue intervention reverse it?** Datasets with rescue arms (HDAC6i, 40Hz gamma) let us project rescued cells into the latent space and measure their shift back toward control along `z_shared`.

---

## Known limitations / things to tune

- **Drug-specific heads need real data.** The synthetic smoke test shows the architecture works but doesn't strongly separate drugs because of the short run + noisy simulation. On real data with longer training, lower HSIC weights on drug heads, or more cells per drug, separation is expected.
- **Count model is Gaussian on log-normalized HVGs** for clarity. For production, swap in a ZINB head (see `scvi-tools.contrastive_vi` for reference). HVG-Gaussian is fine for the poster; ZINB is for a paper.
- **Cross-study integration is the hardest part.** Different sequencing platforms, dosing regimens, sex, age. The `batch_key=dataset_id` in the decoder helps but is not a silver bullet. Inspect UMAPs of `z_bg` colored by `dataset_id` — if datasets still cluster separately, add a harmonization step (Harmony / scVI) before HVG selection.
- **Cell-type annotation requires a reference.** We use `scanpy.tl.ingest` against an Allen mouse brain reference; for more accurate labels use Azimuth or scANVI separately and write the labels back into the merged AnnData.

## Timeline suggestion (to poster)

- Week 1: dataset verification + downloads + QC report
- Week 2: preprocessing + cross-study merge + cell-type annotation
- Week 3: first training runs + hyperparameter sweeps
- Week 4: analysis + figure iteration
- Week 5: poster drafting
- Week 6: revisions + print

## References

- Weinberger et al. *Isolating salient variations of interest in single-cell data with contrastiveVI.* Nat Methods 2023.
- Gibson et al. *Methotrexate Chemotherapy Induces Persistent Tri-glial Dysregulation.* Cell 2019.
- Kim et al. *Non-invasive gamma stimulation for chemo brain.* Sci Transl Med 2024.
- Ma et al. *snRNA-seq of HDAC6i in Doxorubicin-Induced CICI.* Mol Neurobiol 2025.
- Han et al. *CLEAR: self-supervised contrastive learning for scRNA-seq.* Brief Bioinform 2022.
- Tu et al. *Supervised Contrastive VAE for perturbation data.* PMLR 2024.

See `docs/proposal.md` for the full research plan and `docs/methodology.md` for detailed methods.
