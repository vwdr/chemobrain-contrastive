# C3 — Matched baselines: scVI and pairwise contrastiveVI

Script: `scripts/29_baselines_extended.py` (imports `scripts/20_peerreview_matched_baselines.py` unmodified and reuses its
`canonical_split`, `make_adata`, `modeled_gene_log`, `probe`, `linear_cka`, `EPOCHS` and model settings; outputs redirected).
Fit outputs: `runs/peerreview_20261003/baselines/`. Logs: `runs/peerreview_20261003/logs/29_*.log`. Design record:
`C3_design.json`.

## Settings

- scvi-tools 1.5.1 (pinned), CPU, 3 threads per fit, 3 fits concurrently; seeds 0, 1, 2; 100 epochs, no early
  stopping, no hyperparameter search (as in script 20).
- scVI: canonical 6,000 training cells, 1,500 committed genes, batch = study, n_latent 32, n_hidden 128, 2 layers,
  dropout 0.1, NB likelihood, gene dispersion, observed library size; batch size 256. Fit time 144 s.
- contrastiveVI: pairwise within study (GSE216146 control vs cisplatin; GSE271055 control vs doxorubicin), 3,000
  canonical training cells per pair, background latent 16, salient latent 8, n_hidden 128, 2 layers, dropout 0.1,
  wasserstein_penalty 0, batch size 128. Test cells: canonical test cells of the pair (883 cisplatin pair; 1,998
  doxorubicin pair). Fit time 102–106 s.
- The canonical split reproduced by script 20's `canonical_split` was asserted identical to the regenerated canonical
  inputs.

## Reconstruction MSE on both scales (`C3_mse_both_scales.csv`, `C3_mse_summary.csv`)

- Modeled-gene scale (script 20's metric): observed and predicted expression normalized to 1e4 within the 1,500
  modeled genes, log1p.
- Canonical scale (manuscript): log1p(count / all-gene total × 1e4). scvi-tools predictions (ρ·1e4 within modeled
  genes) were converted as log1p(ρ·1e4 · L_model / L_full), L_model = observed modeled-gene total (the library size
  the models condition on), L_full = all-gene total.
- MC-ContrastiveVI (fresh Phase 2 fits; gated posterior means decoded) and PCA(32) predictions p on the canonical
  scale were converted to the modeled-gene scale as log1p(e / Σe × 1e4), e = max(expm1(p), 0). PCA(32) was also fitted
  directly on modeled-gene-scale training data.

| method | test cells | n fits | canonical MSE mean (SD) | modeled-gene MSE mean (SD) |
|---|---|---|---|---|
| MC-ContrastiveVI (fresh) | all canonical test (2,881) | 3 | 0.168723 (0.000123) | 0.904167 (0.006085) |
| MC-ContrastiveVI (fresh) | cisplatin-pair test (883) | 3 | 0.091716 (0.000290) | 0.449616 (0.000906) |
| MC-ContrastiveVI (fresh) | doxorubicin-pair test (1,998) | 3 | 0.202756 (0.000305) | 1.105051 (0.008469) |
| scVI | all canonical test | 3 | 0.206108 (0.002009) | 1.080648 (0.005539) |
| contrastiveVI | cisplatin-pair test | 3 | 0.126042 (0.001984) | 0.714153 (0.005729) |
| contrastiveVI | doxorubicin-pair test | 3 | 0.262185 (0.001555) | 1.504420 (0.010913) |
| PCA(32), fit on canonical scale | all canonical test | 1 | 0.155247 | 0.873536 |
| PCA(32), fit on canonical scale | cisplatin-pair test | 1 | 0.085719 | 0.454030 |
| PCA(32), fit on canonical scale | doxorubicin-pair test | 1 | 0.185974 | 1.058934 |
| PCA(32), fit on modeled-gene scale | all canonical test | 1 | — | 0.550180 |
| PCA(32), fit on modeled-gene scale | cisplatin-pair test | 1 | — | 0.331658 |
| PCA(32), fit on modeled-gene scale | doxorubicin-pair test | 1 | — | 0.646755 |

The canonical-scale PCA(32) value reproduces the committed `runs/corrected_20260920/pca_metrics.json`
(0.15524699 vs 0.15524694).

## scVI probes (`C3_scvi_metrics.csv`)

Pooled probes as in script 20 (`LogisticRegression(max_iter=1500, class_weight='balanced')`, train → test). Within-study
treatment probes use the `latent_probes.csv` definition (`max_iter=800`, balanced; canonical training cells of the study
→ canonical test cells of the study; y = treated).

| seed | 3-class condition | study | cisplatin vs PBS within GSE216146 | doxorubicin vs control within GSE271055 |
|---|---|---|---|---|
| 0 | 0.6564 | 0.9457 | 0.7924 | 0.6219 |
| 1 | 0.6993 | 0.9295 | 0.8427 | 0.6394 |
| 2 | 0.6814 | 0.9287 | 0.7931 | 0.6224 |
| mean | 0.6791 | 0.9346 | 0.8094 | 0.6279 |

## contrastiveVI probes and salient-space utilization (`C3_contrastivevi_metrics.csv`, `C3_contrastivevi_salient_usage_per_dimension.csv`)

Utilization computed from `qs_m`, `qs_v` (module `_generic_inference`) on held-out pair test cells, split into
target (treated) and background (control) cells; active unit = Var[qs_m] > 0.01; KL to N(0, I) per dimension.

| pair | seed | treatment probe: background | treatment probe: salient | target AU (of 8) | target dims KL>0.01 | target total KL (nats) | target Σ Var[qs_m] | background AU | background total KL |
|---|---|---|---|---|---|---|---|---|---|
| cisplatin | 0 | 0.7111 | 0.7245 | 8 | 8 | 7.471 | 4.784 | 8 | 5.700 |
| cisplatin | 1 | 0.6896 | 0.6955 | 8 | 8 | 8.073 | 4.986 | 8 | 6.587 |
| cisplatin | 2 | 0.6496 | 0.6797 | 8 | 8 | 7.297 | 4.650 | 8 | 5.953 |
| doxorubicin | 0 | 0.5684 | 0.5246 | 8 | 8 | 7.511 | 4.229 | 8 | 7.211 |
| doxorubicin | 1 | 0.5600 | 0.6128 | 8 | 8 | 7.310 | 4.052 | 8 | 7.089 |
| doxorubicin | 2 | 0.5769 | 0.5597 | 8 | 8 | 7.036 | 4.102 | 8 | 6.693 |

Per-dimension ranges over seeds (target test cells): cisplatin pair KL 0.539–1.478 nats, posterior-mean variance
0.341–1.203; doxorubicin pair KL 0.771–1.104, variance 0.445–0.621.

## Cross-seed linear CKA on test cells (`C3_cross_seed_cka.csv`; mean, min–max over the 3 seed pairs)

| method | test cells | representation | mean | min–max |
|---|---|---|---|---|
| scVI | all | latent (32) | 0.774 | 0.766–0.791 |
| contrastiveVI | cisplatin pair | background | 0.846 | 0.838–0.858 |
| contrastiveVI | cisplatin pair | salient | 0.444 | 0.418–0.461 |
| contrastiveVI | cisplatin pair | salient, target cells only | 0.535 | 0.500–0.553 |
| contrastiveVI | doxorubicin pair | background | 0.729 | 0.727–0.731 |
| contrastiveVI | doxorubicin pair | salient | 0.124 | 0.111–0.134 |
| contrastiveVI | doxorubicin pair | salient, target cells only | 0.130 | 0.120–0.141 |
| MC-ContrastiveVI (fresh) | all | bg | 0.644 | 0.600–0.700 |
| MC-ContrastiveVI (fresh) | all | shared (gated) | 0.579 | 0.497–0.691 |
| MC-ContrastiveVI (fresh) | all | drug (gated) | 0.377 | 0.177–0.697 |

## Deviations

- None from script 20's model settings. The modeled-gene-scale PCA variant and the per-pair MC-ContrastiveVI MSE rows
  are additions for like-for-like comparison.
