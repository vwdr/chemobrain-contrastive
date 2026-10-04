# R3 — One consistent set of numbers from the fresh canonical fits

## Methods

- Full model: fresh canonical fits from the previous run (`runs/peerreview_20261003/canonical/full_{0,1,2}`; see
  `analysis/peerreview_20261003/A1_fresh_canonical_fits.md`).
- Comparison fits (`scripts/36_r3_comparison_fits.py`): the unmodified `train()` / `pca()` of `scripts/06_validation.py` (modes
  `no_hsic`, `no_gating`, `gaussian_vae`; seeds 0–2; PCA(32)) and `train()` of `scripts/12_nb_sensitivity.py` (seeds 0–2), run by
  importing the scripts and redirecting their module-level `RUN`/`O` globals to `runs/peerreview_20261004/comparison/`. Thread counts
  are those set by the scripts (06: 2 torch threads; 12: 1). Checkpoints, histories, metrics JSON and latents saved there; NB latents
  were computed from the NB checkpoints with the 06 gated-means logic (`scripts/40_r3_analyses.py --part comparison`). 13 jobs,
  13:20–13:35 (ablation fits 129–161 s; NB fits about 5 min).
- `scripts/39_r3_rerun_canonical_scripts.py` executes the unmodified source of `scripts/09_interpretation.py`, `scripts/11_diagnostics.py`
  and `scripts/16_pseudobulk_paired.py` with only the protected path expressions replaced (each replacement asserted to occur once):
  09 → `R3_09_interpretation/` (both IG targets, both baselines, adaptive quadrature, attribution stability incl. cross-seed top-50 Jaccard,
  complete attribution table); 11 → `R3_11_diagnostics/` (variance_robustness with 200-resample cell bootstrap, latent_dependence with
  100-permutation null, latent probes, ablation diagnostics, rescue projection and the 8-assignment sign-flip reference, cross-seed CKA);
  16 → `runs/peerreview_20261004/pseudobulk_task4/` (paired PyDESeq2). Inputs staged by copying the regenerated canonical inputs, the
  fresh full-model checkpoints/latents and the R3 ablation latents into `runs/peerreview_20261004/r3_canonical/`.
- `scripts/40_r3_analyses.py`: comparison table (`R3_comparison_models*.csv`), Figure-5 data (`R3_consensus_top12_figure5_data.csv`),
  enrichment with the functions of `scripts/19_gene_pattern_audit.py` (`R3_enrichment_*.csv`; frozen lists recomputed alongside for
  reference; direction-consistent cisplatin set derived with `derive_cisplatin_de_support` from the regenerated pseudobulk results; the
  scDisInFact top-100 is the frozen set hard-coded in script 19, derived from Task 6 — the committed Task 6 gene-score table contains only
  seed-pair stability statistics, not per-gene scores), pseudobulk comparison (`R3_pseudobulk_*.csv`).
- `scripts/45_r3_paper_numbers.py`: `R3_paper_numbers.csv` (116 rows; every quote checked to occur in the text extraction after
  whitespace normalization).

## Comparison models (`R3_comparison_models_summary.csv`)

| model | n | held-out MSE range | background study accuracy | background source cell-type accuracy | F_sh range |
|---|---|---|---|---|---|
| full | 3 | 0.168616–0.168857 | 0.760–0.782 | 0.985–0.994 | 0.924–0.970 |
| no_hsic | 3 | 0.168448–0.169154 | 0.717–0.856 | 0.986–0.995 | 0.801–0.949 |
| no_gating | 3 | 0.167669–0.167821 | 0.794–0.919 | 0.787–0.862 | 0.302–0.334 |
| gaussian_vae | 3 | 0.167623–0.167943 | 0.802–0.869 | 0.814–0.914 | 0.184–0.330 |
| negative_binomial | 3 | 0.210833–0.211570 | 0.684–0.819 | 0.974–0.991 | 0.864–0.881 |
| pca32 | 1 | 0.155247–0.155247 | 0.999–0.999 | — | — |

Notes: background study accuracy is the metrics-JSON probe of `06_validation.py` (max_iter 500); source cell type follows
`11_diagnostics.py` (`ablation_diagnostics.csv`; NB computed from NB latents with the same definition); PCA study accuracy is from
`06_validation.pca` (PCA scores). F_sh for no_gating / gaussian_vae is computed by the unchanged 06 code on their ungated means.

## Pseudobulk regeneration (`R3_pseudobulk_committed_vs_regenerated.csv`, `R3_pseudobulk_per_celltype_contrast.csv`)

| quantity | committed | regenerated | changed |
|---|---|---|---|
| n_BH_findings_cisplatin_GENUS_vs_cisplatin | 55 | 55 | no |
| n_BH_findings_cisplatin_vs_control | 297 | 293 | yes |
| n_BH_findings_control_GENUS_vs_control | 2 | 2 | no |
| n_cell_type_contrasts_analyzed | 16 | 16 | no |
| n_findings_all3_pairs_same_direction | 353 | 349 | yes |
| n_welch_fdr_005_total | 1 | 1 | no |
| pooled_BH_findings_all_contrasts | 354 | 350 | yes |
| spearman_deseq2_vs_paired_max_all_rows | 0.996953 | 0.996955 | no |
| spearman_deseq2_vs_paired_max_rows_with_findings | 0.996953 | 0.996955 | no |
| spearman_deseq2_vs_paired_min_all_rows | 0.932491 | 0.932491 | no |
| spearman_deseq2_vs_paired_min_rows_with_findings | 0.932491 | 0.932491 | no |

## Changed manuscript values (`R3_paper_numbers.csv`, rows with changed = yes)

| section | quantity | manuscript | fresh |
|---|---|---|---|
| abstract | MC-ContrastiveVI held-out MSE range | 0.1687–0.1688 | 0.1686–0.1689 |
| abstract | pooled F_sh range (three primary fits) | 92.1–95.7% | 92.4–97.0% |
| 2.2 | MC-ContrastiveVI held-out MSE range | 0.1687–0.1688 | 0.1686–0.1689 |
| 2.2 | no HSIC held-out MSE range (qualitative statement; values from Figure 3A data) | committed 0.1683–0.1689 | 0.1684–0.1692 |
| 2.2 | no gating held-out MSE range (qualitative statement; values from Figure 3A data) | committed 0.1676–0.1678 | 0.1677–0.1678 |
| 2.2 | noncontrastive Gaussian VAE held-out MSE range (qualitative statement; values from Figure 3A data) | committed 0.1674–0.1678 | 0.1676–0.1679 |
| 2.2 | NB sensitivity held-out MSE range | 0.2109–0.2123 | 0.2108–0.2116 |
| 2.2 | background study balanced accuracy range | 0.726–0.856 | 0.760–0.782 |
| 2.2 | background source cell-type balanced accuracy range | 0.985–0.991 | 0.985–0.994 |
| 2.3 | pooled F_sh range | 92.1–95.7% | 92.4–97.0% |
| 2.3 | active shared coordinates per fit | 0, 0, 1 | 0, 0, 0 |
| 2.3 | NB F_sh range | 82.9–89.2% | 86.4–88.1% |
| 2.4 | shared top-50 cross-seed Jaccard (norm squared, zero baseline) | 0.190–0.639 | 0.449–0.562 |
| 2.4 | consensus list sizes | 34 / 20 / 33 | 40 / 27 / 31 |
| Figure 5 | top-12 consensus genes, shared | Mrc1, Pf4, Ms4a7, Cbr2, Cd163, F13a1, Gpx3, Dab2, mt-Co3, Cd52, Stab1, C5ar1 | Mrc1, Cbr2, Pf4, Ms4a7, F13a1, Cd163, Apoe, Gpx3, Cd52, mt-Nd4, Ttr, mt-Co2 |
| Figure 5 | top-12 consensus genes, doxorubicin | Snhg11, Kcnq1ot1, Meg3, Ttr, Gria2, Rtn1, Grin2b, C1ql3, Opcml, Pcdh9, Nrcam, Dclk1 | Snhg11, Meg3, Ttr, Apoe, Atp1b1, Kcnq1ot1, Nrxn1, Opcml, Arpp21, Sptbn1, Grm5, Ptprn |
| Figure 5 | top-12 consensus genes, cisplatin | Ccn3, Col25a1, Igfbp6, Mgp, Foxc2, Rspo3, Col6a2, mt-Nd4, Ogn, Col6a1, Adamtsl3, Foxd1 | Mrc1, mt-Co2, Ttr, Apoe, mt-Nd4, Ccn3, Mgp, Glul, Ctsd, Fos, Bsg, Junb |
| 2.5 | shared consensus: M8 / GO BP terms at FDR 0.05 | 46 / 11 | 73 / 6 |
| 2.5 | top M8 term (ZHANG_UTERUS_C5_MACROPHAGE): overlap / q | 14 of 34; 4.18e-10 | 19 of 40; 7.35e-16 |
| 2.5 | Reactome immune system q (shared) | 0.092 | 0.000596 (overlap 13) |
| 2.5 | cisplatin consensus: M8 / Reactome terms at FDR 0.05 | 77 / 5 | 56 / 0 |
| 2.5 | Reactome ECM degradation (cisplatin): overlap / fold / q | 5 of 33; 7.58; 0.034 | 4 of 31; 6.45; 0.321 |
| 2.5 | doxorubicin consensus: terms at FDR 0.05 (size) | 0 (20 genes) | 2 (27 genes): Reactome 2, GO 0, M8 0 |
| 2.5 | cisplatin vs control BH findings | 297 | 293 |
| 2.5 | pooled BH findings | 354 | 350 |
| 2.5 | findings with same direction in all 3 pairs | 353 | 349 |
| 2.5 | cisplatin consensus ∩ direction-consistent pseudobulk set | 5: Apoe, Mbp, Mt1, Tmsb10, Ttr | fresh list 5: Apoe, Junb, Mt1, Ttr, mt-Co2 |
| 2.6 | cisplatin rescue-minus-treatment mean distance difference | 0.0003 | -0.0032 |

Totals: 116 rows; changed = yes 28; no 77; not-recomputable 11.
Rows marked "no" include values taken from committed files for analyses that do not depend on the archived fits (cohort, Task 1
reruns, Task 2, Task 5, Task 6, annotation) and methods settings; each such row says so in `notes`.

## Deviations

- The comparison fits use the regenerated canonical inputs (content-identical to the committed inputs; X max abs diff 9.5e-7).
- Figure 4 and Figure 6 values are not stated numerically in the text; they are listed as not-recomputable with the fresh data file.
- The Figure 5 gene order for the fresh lists is the order of `A1_fresh_consensus.json` (mean of the six ranking vectors).
- No refit of the full model: the full-model numbers come from the previous run's fresh fits.
