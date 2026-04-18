# Research Proposal

## Disentangling Shared and Drug-Specific Molecular Signatures of Chemotherapy-Induced Cognitive Impairment via Contrastive Latent-Variable Modeling of Mouse Brain scRNA-seq

### Summary

Chemotherapy-induced cognitive impairment (CICI, "chemo brain") is a clinically significant side effect affecting a large fraction of cancer survivors, yet its molecular etiology is fragmented across drugs (cisplatin, methotrexate, doxorubicin, paclitaxel, 5-FU, fludarabine/cyclophosphamide) and across cell types (microglia, oligodendrocytes, endothelial cells, neurons). We will build a unified contrastive deep-learning framework that leverages existing mouse brain scRNA-seq / snRNA-seq datasets to disentangle chemo-induced variation from shared biological variation across animals, then perform a cross-drug analysis to identify (a) a convergent core toxicity signature shared across agents and (b) drug-specific toxicity axes. The model's salient latent factors will be interpreted mechanistically and, where behavioral data are available from the source studies, tested for association with cognitive phenotypes, producing a prioritized list of candidate mediators of cognitive decline.

### Background and gap

Existing CICI scRNA-seq studies are typically single-drug, single-region, and analyzed with standard differential expression pipelines that conflate subtle treatment effects with inter-animal variance and cell-state heterogeneity. Contrastive latent-variable methods such as contrastiveVI (Weinberger et al., *Nat Methods* 2023), scDisInFact (Zhang & Zhang, *Nat Commun* 2024), and SC-VAE (Tu et al., PMLR 2024) can disentangle these sources within a single study but have not been applied systematically across chemo agents. No current work asks: is there a shared molecular chemo-brain latent dimension across chemically unrelated drugs, and if so, which cell types and pathways define it?

### Central hypothesis

There exists a low-dimensional, cell-type-specific convergent toxicity manifold in brain scRNA-seq space that is shared across structurally diverse chemotherapeutic agents. Position along this manifold predicts cognitive impairment severity better than gene-level differential expression signatures computed in isolation.

### What makes this novel

1. First cross-drug contrastive latent-variable meta-analysis of CICI scRNA-seq data.
2. Methodological contribution: extending contrastive analysis models to handle multiple simultaneous target conditions (drugs) rather than just treatment-vs-control, so shared and drug-specific axes can be separated within one model.
3. Integration of learned latent factors with rescue-arm data (HDAC6 inhibition in doxorubicin-CICI; 40 Hz gamma stimulation in cisplatin-CICI) to validate that the shared-toxicity axis responds to independently validated interventions.
4. An openly released, annotated latent-space atlas of chemo brain that other groups can project new data onto.

### Specific Aims

**Aim 1 — Build a harmonized cross-study mouse chemo-brain scRNA-seq corpus.**
Collect public datasets from GEO/SRA covering at least 3 chemotherapeutic agents with matched vehicle controls (target: doxorubicin, cisplatin, fludarabine/cyclophosphamide; stretch: methotrexate, paclitaxel). Harmonize using scanpy-based QC and batch-aware HVG selection. Annotate cell types consistently against the Allen Brain Atlas reference taxonomy via label transfer.

**Aim 2 — Develop a multi-condition contrastive disentanglement model.**
Extend contrastiveVI to learn: (i) a shared background latent space capturing biological variance common to treated and untreated mice; (ii) a shared-toxicity salient space capturing variation common across drugs relative to controls; (iii) drug-specific salient subspaces capturing what is unique to each agent. Disentanglement is enforced with HSIC independence penalties. Validate that drug-specific dimensions recover known drug-specific biology as positive controls (e.g. paclitaxel → endothelial senescence; methotrexate → oligodendrocyte lineage defects).

**Aim 3 — Identify and functionally interpret convergent toxicity factors.**
For each cell type, extract the top-contributing genes along the shared-toxicity axis using integrated gradients on the encoder. Enrich for pathways (oxidative stress, neuroinflammation, senescence-associated secretory phenotype, mitochondrial dysfunction, myelin biology). For datasets with rescue arms or reported behavioral data, correlate per-animal pseudo-bulk latent position with the rescue/cognitive outcome. Prioritize candidate mediators by (a) cross-drug consistency, (b) cell-type specificity, (c) rescue-arm reversal, (d) druggability.

### Expected deliverables

- A publicly available annotated chemo-brain latent-space atlas.
- A multi-condition contrastive analysis model (open-source code) generalizable beyond CICI.
- A prioritized, interpretable list of 5–20 candidate molecular mediators of CICI with cell-type- and drug-specific annotations.
- Falsifiable mechanistic hypotheses suitable for follow-up validation.
- A conference poster summarizing methodology, findings, and candidate targets.

### Limitations and mitigations

- **Dataset heterogeneity** (platforms, dosing, timepoints) is the largest risk — mitigated by rigorous batch modeling and by treating timepoint/dose as covariates rather than ignoring them.
- **Confounding with sickness behavior/cachexia** — include datasets with pair-fed or sham-stressed controls where available.
- **Mouse-to-human translation** is limited; future work can project human post-chemo snRNA-seq onto the mouse-learned latent space.
- **Contrastive model identifiability** is theoretically non-trivial; include ablations (no-HSIC, contrastiveVI baseline) and synthetic-data benchmarks to verify disentanglement.

### Feasibility

All candidate datasets are publicly available on GEO; the methodological primitives (contrastiveVI, scVI, HSIC regularization) have open-source implementations. The project is tractable on a single GPU (24 GB VRAM recommended) for the datasets at hand. Bottleneck is curation and integration, not compute.
