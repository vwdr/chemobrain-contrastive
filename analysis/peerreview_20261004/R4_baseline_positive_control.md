# R4 — Positive control on baseline models (semi-synthetic data reused)

Scripts: `scripts/41_r4_prepare.py` (injected raw counts), `scripts/42_r4_multigroupvi.py` (.venv-mgvi), `scripts/43_r4_contrastivevi.py`
(.venv), `scripts/44_r4_mccvi_and_collect.py` (MC-ContrastiveVI attribution on existing grid checkpoints; collection),
`src/peerreview/ig_generic.py` (integrated gradients). Fit outputs: `runs/peerreview_20261004/r4/`. Logs: `runs/peerreview_20261004/logs/4[234]_*.log`,
queue log `r4_queue.log`.

## Design

- Data: Phase 4 base cells, split, responders and gene sets (`runs/peerreview_20261003/semisynthetic/base.npz`, `D1_gene_sets.csv`). Injected
  raw counts rebuilt with the Phase 4 programs and `default_rng(20261003)` thinning draws; log-normalizing them reproduces the Phase 4
  model inputs exactly (asserted). Configurations: null, shared/specific/mixed × δ ∈ {1, 2}, responder fraction 1.0; seeds 0–2.
- multiGroupVI: groups pseudo-control (both studies), pseudo-cisplatin (GSE216146 pseudo-treated), pseudo-doxorubicin (GSE271055
  pseudo-treated); Phase 4 training (6,000) and test (1,429) cells; Task C4 settings (package defaults; 400 epochs). 21 fits.
- contrastiveVI: per study, background = pseudo-control, target = pseudo-treated; the study's Phase 4 training (3,000) and test cells;
  Task C3 settings (100 epochs). 42 fits.
- MC-ContrastiveVI: the existing Phase 4 grid fits (no refit); probes from `R1_probes_per_fit.csv` (`raw_drug`), active units and KL
  from the Phase 4 metrics (shared block + matching drug block, 8 + 4 dimensions).
- Metrics: canonical-scale MSE (log1p(count / all-gene total × 1e4); scvi-tools predictions converted with the observed modeled-gene
  library as in C3); active units (Var[posterior mean] > 0.01) and total KL of the treatment-related latent on the treated group's
  test cells (multiGroupVI: private block of the treated group, 10 dims; contrastiveVI: salient, 8 dims; MC-ContrastiveVI: shared +
  matching drug block, 12 dims); within-study probe pseudo-treated vs pseudo-control on the ungated treatment-related latent
  (LogisticRegression balanced, max_iter 800; for MC-ContrastiveVI the ungated drug block).
- Attribution: integrated gradients of the squared norm of the posterior mean of the treatment-related latent with respect to the
  encoder input (multiGroupVI group encoder: raw counts, as fed by the package; contrastiveVI salient encoder: log1p counts;
  MC-ContrastiveVI: log-normalized expression, concatenated ungated shared + matching drug block), zero and control-median baselines
  (median over the study's training pseudo-controls, in the same input space), 128 test pseudo-treated cells per group drawn with
  `default_rng(314)` (09 scheme), adaptive Gauss–Legendre quadrature with the 09 completeness rule. Top 50 genes by mean |IG| vs the
  union of genes injected into that group (GSE216146 treated: S ∪ A; GSE271055 treated: S ∪ B; 80 genes; random expectation
  50 × 80 / 1,500 = 2.67 hits), and per-set breakdown.

## Comparison table (`R4_comparison_table.csv`; means over seeds 0–2)

| model | configuration | group | probe (within study) | active units / dims | total KL (nats) | top-50 hits vs union, zero baseline (min–max) | control-median baseline | S / A / B hits (zero) | random expectation |
|---|---|---|---|---|---|---|---|---|---|
| multiGroupVI | mixed_d1_rf1.0 | cisplatin | 0.642 | 10.0 / 10 | 9.212 | 10.00 (9–11) | 9.67 | 4.67 / 5.33 / — | 2.67 |
| multiGroupVI | mixed_d1_rf1.0 | doxorubicin | 0.661 | 9.3 / 10 | 6.965 | 7.33 (7–8) | 8.00 | 1.33 / — / 6.00 | 2.67 |
| multiGroupVI | mixed_d2_rf1.0 | cisplatin | 0.768 | 10.0 / 10 | 9.113 | 11.00 (10–12) | 10.67 | 4.67 / 6.33 / — | 2.67 |
| multiGroupVI | mixed_d2_rf1.0 | doxorubicin | 0.804 | 9.0 / 10 | 6.275 | 6.33 (6–7) | 7.67 | 2.00 / — / 4.33 | 2.67 |
| multiGroupVI | null | cisplatin | 0.485 | 10.0 / 10 | 9.894 | 4.67 (4–5) | 4.33 | 2.67 / 2.00 / — | 2.67 |
| multiGroupVI | null | doxorubicin | 0.489 | 9.0 / 10 | 7.667 | 2.33 (2–3) | 2.00 | 1.00 / — / 1.33 | 2.67 |
| multiGroupVI | shared_d1_rf1.0 | cisplatin | 0.587 | 10.0 / 10 | 9.873 | 6.67 (5–8) | 6.33 | 5.00 / 1.67 / — | 2.67 |
| multiGroupVI | shared_d1_rf1.0 | doxorubicin | 0.582 | 9.7 / 10 | 6.924 | 3.67 (3–4) | 3.67 | 2.33 / — / 1.33 | 2.67 |
| multiGroupVI | shared_d2_rf1.0 | cisplatin | 0.671 | 10.0 / 10 | 8.982 | 6.00 (5–7) | 6.33 | 4.33 / 1.67 / — | 2.67 |
| multiGroupVI | shared_d2_rf1.0 | doxorubicin | 0.691 | 8.3 / 10 | 6.320 | 5.00 (4–6) | 5.67 | 3.67 / — / 1.33 | 2.67 |
| multiGroupVI | specific_d1_rf1.0 | cisplatin | 0.615 | 10.0 / 10 | 8.929 | 9.00 (8–10) | 7.33 | 2.00 / 7.00 / — | 2.67 |
| multiGroupVI | specific_d1_rf1.0 | doxorubicin | 0.581 | 9.0 / 10 | 7.369 | 6.33 (6–7) | 6.00 | 1.00 / — / 5.33 | 2.67 |
| multiGroupVI | specific_d2_rf1.0 | cisplatin | 0.718 | 10.0 / 10 | 9.226 | 9.67 (9–10) | 9.67 | 2.33 / 7.33 / — | 2.67 |
| multiGroupVI | specific_d2_rf1.0 | doxorubicin | 0.697 | 8.7 / 10 | 7.142 | 8.33 (8–9) | 8.33 | 1.00 / — / 7.33 | 2.67 |
| contrastiveVI | mixed_d1_rf1.0 | cisplatin | 0.565 | 8.0 / 8 | 6.267 | 6.33 (6–7) | 7.00 | 2.33 / 4.00 / — | 2.67 |
| contrastiveVI | mixed_d1_rf1.0 | doxorubicin | 0.589 | 8.0 / 8 | 7.178 | 1.67 (1–2) | 3.00 | 0.67 / — / 1.00 | 2.67 |
| contrastiveVI | mixed_d2_rf1.0 | cisplatin | 0.664 | 8.0 / 8 | 7.701 | 7.00 (6–8) | 7.67 | 3.33 / 3.67 / — | 2.67 |
| contrastiveVI | mixed_d2_rf1.0 | doxorubicin | 0.732 | 8.0 / 8 | 7.914 | 2.33 (1–3) | 3.67 | 0.67 / — / 1.67 | 2.67 |
| contrastiveVI | null | cisplatin | 0.520 | 8.0 / 8 | 5.195 | 4.33 (2–6) | 5.00 | 2.33 / 2.00 / — | 2.67 |
| contrastiveVI | null | doxorubicin | 0.489 | 8.0 / 8 | 6.915 | 2.33 (2–3) | 3.00 | 1.00 / — / 1.33 | 2.67 |
| contrastiveVI | shared_d1_rf1.0 | cisplatin | 0.518 | 8.0 / 8 | 5.816 | 5.67 (3–7) | 4.67 | 2.67 / 3.00 / — | 2.67 |
| contrastiveVI | shared_d1_rf1.0 | doxorubicin | 0.531 | 8.0 / 8 | 7.019 | 1.67 (0–3) | 1.67 | 0.67 / — / 1.00 | 2.67 |
| contrastiveVI | shared_d2_rf1.0 | cisplatin | 0.596 | 8.0 / 8 | 6.824 | 6.00 (5–8) | 6.00 | 3.00 / 3.00 / — | 2.67 |
| contrastiveVI | shared_d2_rf1.0 | doxorubicin | 0.626 | 8.0 / 8 | 7.336 | 2.33 (1–3) | 3.67 | 1.00 / — / 1.33 | 2.67 |
| contrastiveVI | specific_d1_rf1.0 | cisplatin | 0.545 | 8.0 / 8 | 6.002 | 5.00 (4–6) | 4.67 | 1.67 / 3.33 / — | 2.67 |
| contrastiveVI | specific_d1_rf1.0 | doxorubicin | 0.537 | 8.0 / 8 | 6.928 | 2.00 (2–2) | 3.33 | 1.00 / — / 1.00 | 2.67 |
| contrastiveVI | specific_d2_rf1.0 | cisplatin | 0.591 | 8.0 / 8 | 6.959 | 7.00 (6–8) | 7.00 | 2.33 / 4.67 / — | 2.67 |
| contrastiveVI | specific_d2_rf1.0 | doxorubicin | 0.628 | 8.0 / 8 | 7.598 | 2.33 (2–3) | 3.00 | 1.00 / — / 1.33 | 2.67 |
| MC-ContrastiveVI | null | cisplatin | 0.513 | 0.3 / 12 | 0.010 | 4.67 (4–6) | 4.00 | 3.67 / 1.00 / — | 2.67 |
| MC-ContrastiveVI | null | doxorubicin | 0.499 | 0.3 / 12 | 0.008 | 3.33 (2–5) | 2.00 | 1.33 / — / 2.00 | 2.67 |
| MC-ContrastiveVI | shared_d1_rf1.0 | cisplatin | 0.511 | 0.3 / 12 | 0.015 | 2.67 (1–4) | 3.33 | 2.00 / 0.67 / — | 2.67 |
| MC-ContrastiveVI | shared_d1_rf1.0 | doxorubicin | 0.489 | 0.3 / 12 | 0.014 | 2.33 (1–4) | 2.33 | 1.33 / — / 1.00 | 2.67 |
| MC-ContrastiveVI | specific_d1_rf1.0 | cisplatin | 0.518 | 0.3 / 12 | 0.017 | 3.00 (2–4) | 3.67 | 2.00 / 1.00 / — | 2.67 |
| MC-ContrastiveVI | specific_d1_rf1.0 | doxorubicin | 0.501 | 0.3 / 12 | 0.016 | 2.33 (2–3) | 1.33 | 1.33 / — / 1.00 | 2.67 |
| MC-ContrastiveVI | mixed_d1_rf1.0 | cisplatin | 0.512 | 1.3 / 12 | 0.028 | 1.67 (1–3) | 3.00 | 1.33 / 0.33 / — | 2.67 |
| MC-ContrastiveVI | mixed_d1_rf1.0 | doxorubicin | 0.493 | 1.3 / 12 | 0.028 | 2.33 (2–3) | 2.67 | 0.67 / — / 1.67 | 2.67 |
| MC-ContrastiveVI | shared_d2_rf1.0 | cisplatin | 0.512 | 0.0 / 12 | 0.015 | 3.33 (2–6) | 4.67 | 1.00 / 2.33 / — | 2.67 |
| MC-ContrastiveVI | shared_d2_rf1.0 | doxorubicin | 0.523 | 0.0 / 12 | 0.014 | 2.00 (1–3) | 2.67 | 1.00 / — / 1.00 | 2.67 |
| MC-ContrastiveVI | specific_d2_rf1.0 | cisplatin | 0.509 | 0.3 / 12 | 0.018 | 4.00 (3–5) | 4.33 | 3.00 / 1.00 / — | 2.67 |
| MC-ContrastiveVI | specific_d2_rf1.0 | doxorubicin | 0.511 | 0.3 / 12 | 0.016 | 2.00 (1–3) | 2.33 | 1.00 / — / 1.00 | 2.67 |
| MC-ContrastiveVI | mixed_d2_rf1.0 | cisplatin | 0.537 | 0.0 / 12 | 0.013 | 2.67 (2–4) | 4.00 | 1.33 / 1.33 / — | 2.67 |
| MC-ContrastiveVI | mixed_d2_rf1.0 | doxorubicin | 0.525 | 0.0 / 12 | 0.012 | 1.00 (0–2) | 2.33 | 0.33 / — / 0.67 | 2.67 |

## Canonical-scale MSE on test cells (mean over seeds; contrastiveVI: mean over the two study-specific models)

| configuration | MC-ContrastiveVI | contrastiveVI | multiGroupVI |
|---|---|---|---|
| mixed_d1_rf1.0 | 0.1620 | 0.1166 | 0.2315 |
| mixed_d2_rf1.0 | 0.1617 | 0.1161 | 0.2314 |
| null | 0.1634 | 0.1183 | 0.2317 |
| shared_d1_rf1.0 | 0.1625 | 0.1175 | 0.2318 |
| shared_d2_rf1.0 | 0.1626 | 0.1168 | 0.2324 |
| specific_d1_rf1.0 | 0.1626 | 0.1170 | 0.2317 |
| specific_d2_rf1.0 | 0.1625 | 0.1175 | 0.2318 |

## Attribution numerics

- MC-ContrastiveVI: maximum quadrature points 1024; maximum median scaled completeness error 0.0099.
- contrastiveVI: maximum quadrature points 256; maximum median scaled completeness error 0.0025.
- multiGroupVI: maximum quadrature points 256; maximum median scaled completeness error 0.0012.

Fit times: multiGroupVI 1068–1136 s per fit; contrastiveVI 87–93 s per fit.

## Notes

- The MC-ContrastiveVI rows use 12-dimensional treatment-related latents (8 shared + 4 drug); multiGroupVI 10; contrastiveVI 8.
- Probes and utilization for MC-ContrastiveVI are taken from existing files (no refit); its MSE is the Phase 4 per-fit test MSE.
- All fits completed on the first attempt unless listed in PROGRESS.md.
