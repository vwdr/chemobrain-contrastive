# Poster Outline

Target: standard academic conference poster (36" × 48" portrait or 48" × 36" landscape).

## Layout (4 columns, landscape)

### Column 1 — The Problem

**Title:** Disentangling Shared and Drug-Specific Molecular Signatures of Chemotherapy-Induced Cognitive Impairment

**Authors / Affiliations**

**Background panel** (short text + 1 figure):
- "Chemo brain" (CICI) affects a majority of cancer survivors
- Mechanisms proposed across many studies: neuroinflammation, oligodendrocyte/myelin damage, endothelial senescence, reduced hippocampal neurogenesis, oxidative stress
- But studies are single-drug, single-region, and differential-expression-based — signal gets lost in inter-animal and cell-state variance
- **The gap:** no one has asked whether structurally different chemo agents converge on a shared molecular signature

**Figure 1:** schematic — 4 chemo agents (cisplatin, doxorubicin, methotrexate, fludarabine/cyclophosphamide) pointing to "???" molecular program → cognitive impairment

### Column 2 — The Method

**Question:** Is there a convergent CICI signature across drugs, and which cell types carry it?

**Approach:** Multi-Condition Contrastive Variational Inference (MC-ContrastiveVI)

**Figure 2:** model schematic
- 3 encoders: background, shared-toxicity, drug-specific (one per drug)
- Hard gating: control cells have zero salients
- HSIC penalties enforce independence
- Batch-aware decoder handles cross-study integration

**Datasets table** (Figure 2 inset):
| Drug | Region | Assay | Source | Has rescue arm? |
| --- | --- | --- | --- | --- |
| Doxorubicin | Hippocampus | snRNA-seq | Ma 2025 | ✓ (HDAC6i) |
| Cisplatin | Cortex | scRNA-seq | Kim/Tsai 2024 | ✓ (40 Hz γ) |
| Fludarabine/Cyclophosphamide | Whole brain | snRNA-seq | Monje 2025 | (CAR-T arm) |
| Methotrexate | — | bulk-only | Gibson 2019 | pathway check only |

**Key innovation:** the architecture separates *what's common across all chemo agents* from *what's drug-specific*, in a single training run.

### Column 3 — The Results

**Figure 3:** UMAPs side-by-side
- (A) UMAP of `z_bg` colored by cell type — confirms biological structure is preserved
- (B) UMAP of `z_bg` colored by dataset_id — confirms batch is well-mixed (integration worked)
- (C) UMAP of `z_shared` for treated cells, colored by drug — convergence: do drugs mix or separate? (Expect mixing if the hypothesis holds)

**Figure 4:** variance decomposition bar
- % of treated-cell latent variance explained by shared-toxicity axis vs drug-specific axes
- Quantifies the convergence claim

**Figure 5:** cell-type engagement heatmap
- Rows: cell types; Columns: drugs; Color: mean magnitude along `z_shared`
- Identifies which cell types carry the convergent signature

**Figure 6:** top genes + pathway enrichment
- Panel A: Top 20 attribution genes for `z_shared`
- Panel B: Reactome / MSigDB enrichments for the top-100 gene set
- Panel C: Drug-specific axes' top genes side-by-side — do they match known drug-specific biology?

**Figure 7 (if rescue-arm data delivers):**
- Rescue validation: `‖z_shared‖` distribution in control vs treated vs rescued (HDAC6i and 40 Hz)
- If rescue cells shift toward control, the learned axis is biologically meaningful

### Column 4 — The Interpretation

**Key findings** (4–6 concise bullets):
- (Hypothesis-dependent — fill in after training. Examples below.)
- "Shared-toxicity axis explains X% of treatment-induced variance, supporting convergent CICI mechanism."
- "Microglia and oligodendrocyte-lineage cells carry the strongest signal along the shared axis."
- "Top shared-axis genes enrich for senescence-associated secretory phenotype (SASP) and TNF signaling."
- "Rescue interventions shift `‖z_shared‖` toward control, validating the axis as a functional readout."
- "Drug-specific axes recover known biology: cisplatin → DNA damage response; paclitaxel → microtubule/endothelial."

**Candidate targets:**
- A ranked table of 5–10 genes/pathways with cell-type, drug-consistency, and druggability annotations.

**Limitations** (honest, brief):
- Cross-study batch effects are the largest confound; mitigated but not eliminated.
- Mouse-to-human translation not yet validated.
- Rescue-arm datasets are a strong but imperfect functional proxy.

**Future directions:**
- Project human post-chemo snRNA-seq (where available) onto the mouse-learned space
- Spatial transcriptomics integration for region-specific targeting
- In-vitro validation of top candidates

**References** (small font at bottom) — ~5–8 key citations.

**QR code** linking to the public GitHub repo and interactive latent-space viewer (CELLxGENE or a custom Plotly Dash app).

## Design notes

- Keep text per panel under ~120 words.
- Two high-impact figures per column.
- A consistent color palette across all figures (drugs: categorical; cell types: viridis subset; latent magnitude: sequential).
- Title in one very clear sentence + a short subtitle.
- Big result callout in Column 3 — one sentence in large font that captures the main finding, e.g. "A single latent axis, shared across 3 chemically unrelated chemotherapies, explains X% of treatment-induced single-cell variance and engages microglia and oligodendrocytes most strongly."
