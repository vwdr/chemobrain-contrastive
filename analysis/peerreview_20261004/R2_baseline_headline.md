# R2 — Baseline headline table

Script: `scripts/37_r2_baseline_headline.py`. All values are taken from existing result files (no refits); sources are
listed per row in `R2_baseline_headline.csv` (`source_files`). Means over seeds; per-seed values in the CSV.

Definitions: canonical-scale MSE = log1p(count / all-gene total × 1e4); modeled-gene MSE = normalization to 1e4 within the
1,500 modeled genes, log1p (see `analysis/peerreview_20261003/C3_baselines.md`). Within-study treatment probes:
`LogisticRegression(max_iter=800, class_weight='balanced')`, canonical training cells of the study → canonical test cells
of the study, y = treated; **ungated representations only**. Active units: Var[posterior mean] > 0.01 on held-out treated
(or own-group) cells. Cross-seed CKA: linear CKA on test cells, mean [min–max] over the three seed pairs.

| model | test cells | canonical MSE | modeled-gene MSE | probe representation | cisplatin vs PBS (GSE216146) | doxorubicin vs control (GSE271055) | secondary representation: cis / dox | treatment-block active units (per seed) | treatment-block total KL, nats (per seed) | cross-seed CKA |
|---|---|---|---|---|---|---|---|---|---|---|
| MC-ContrastiveVI (fresh fits) | all canonical test cells | 0.1687 | 0.9042 | drug block, ungated posterior means (8 dims) | 0.5303 | 0.5537 | shared block, ungated posterior means (8 dims): 0.5649 / 0.5441 | sh 0/dox 0/cis 0; sh 0/dox 0/cis 0; sh 0/dox 0/cis 0 | sh 0.0125/dox 0.0003/cis 0.0034; sh 0.0109/dox 0.0004/cis 0.0012; sh 0.0182/dox 0.0005/cis 0.0028 | shared(gated) 0.579 [0.497-0.691]; drug(gated) 0.377 [0.177-0.697] |
| PCA(32) | all canonical test cells | 0.1552 | 0.8735 | not available (no within-study probe in existing files) | — | — | —: — / — | not applicable | not applicable | not applicable (deterministic) |
| scVI | all canonical test cells | 0.2061 | 1.0806 | latent (32 dims; single block) | 0.8094 | 0.6279 | —: — / — | not computed in existing files | not computed in existing files | latent 0.774 [0.766-0.791] |
| contrastiveVI (cisplatin pair) | cisplatin-pair canonical test cells | 0.1260 | 0.7142 | salient (8 dims, ungated posterior means) | 0.6999 | — | background (16 dims): 0.6834 / — | 8/8; 8/8; 8/8 | 7.471; 8.073; 7.297 | salient 0.444 [0.418-0.461]; salient target cells 0.535 [0.500-0.553] |
| contrastiveVI (doxorubicin pair) | doxorubicin-pair canonical test cells | 0.2622 | 1.5044 | salient (8 dims, ungated posterior means) | — | 0.5657 | background (16 dims): — / 0.5685 | 8/8; 8/8; 8/8 | 7.511; 7.310; 7.036 | salient 0.124 [0.111-0.134]; salient target cells 0.130 [0.120-0.141] |
| multiGroupVI | all canonical test cells | 0.2454 | 1.4625 | group-specific latents, ungated posterior means (3 x 10 dims) | 0.8297 | 0.5698 | shared latent (10 dims): 0.6442 / 0.5623 | cis 10/10, dox 7/10; cis 10/10, dox 6/10; cis 10/10, dox 4/10 | cis 7.963, dox 6.485; cis 8.812, dox 5.889; cis 8.331, dox 5.950 | group-specific ungated 0.349 [0.306-0.381]; shared 0.762 [0.755-0.767] |
| scDisInFact (committed Task 6) | all canonical test cells | 0.2313 | — | unshared-bio (condition) latent; within-study probes not available in committed files | — | — | —: — / — | not available in committed files | not available in committed files | unshared_bio 0.842 [0.781-0.881]; shared_bio 0.815 [0.794-0.854] |

## gated_label_masked_probe (multiGroupVI only)

Within-study treatment probes on the multiGroupVI group-specific latents **after** multiplication by the group-label mask: cis 0.9503/0.9627/0.9534 (mean 0.9555); dox 0.9006/0.8943/0.8808 (mean 0.8919).
Gating by group label makes these values structurally label-informative (the mask itself encodes the group); they are listed
only in this separate column and are not used in the probe columns above.

## Not available in existing files

- PCA(32): no within-study treatment probe; no latent-utilization quantity (not a probabilistic model).
- scVI: active units / KL were not computed in Task C3.
- scDisInFact (committed Task 6): no within-study probes, no modeled-gene MSE, no active units/KL; only pooled 3-class condition probes (in `notes`).
- MC-ContrastiveVI cross-seed CKA is available only for gated blocks in the existing files.
- contrastiveVI probe values are the treatment probes of Task C3, computed with script 20's `probe()` (`max_iter=1500`, balanced) on the pair's training and test cells; they are within-study by construction.
