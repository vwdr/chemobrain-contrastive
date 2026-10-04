# D1 — Semi-synthetic positive control (binomial thinning)

Scripts: `scripts/31_semisynthetic.py` (prepare / fit / perm / score), `src/peerreview/semisynth.py` (thinning),
`tests/test_peerreview_thinning.py` (unit test; passes), `scripts/35_d1_signal_check.py` (diagnostic),
`scripts/30_calibration_pseudobulk.py` (calibration), `scripts/33_d1_report_tables.py` (tables below).
Fits: `runs/peerreview_20261003/semisynthetic/` (per configuration: checkpoints, latents, histories, IG summaries,
`done.json`). Logs: `runs/peerreview_20261003/logs/31_*.log`, queue log `31_queue.log`.

## Design (as pre-specified)

1. Base cells: nonrescue controls only — GSE216146 PBS (PN, 4,830 cells; 3 libraries) and GSE271055 control (CNT,
   9,455 cells; 1 library): 14,285 cells. Within each library (`sample_id`, sorted) cells were split 50/50
   (`default_rng(20261003)`, floor(n/2) pseudo-treated). Pseudo-treated cells were coded cisplatin (d = 2) in GSE216146
   and doxorubicin (d = 1) in GSE271055; pseudo-controls d = 0; study one-hot as in the canonical model. Counts by
   library/split: `D1_base_cell_counts.csv`.
2. Genes: the 1,500 genes of `analysis/corrected_20260920/benchmark_gene_universe.csv` (no reselection).
3. Gene sets S, A, B (40 genes each, disjoint): deciles (`pandas.qcut`, 10 bins) of mean log-normalized expression over
   the uninjected base cells; per decile one permutation (`default_rng(20261003)`), genes 1–4 → S, 5–8 → A, 9–12 → B;
   within each set and decile the first two drawn are up and the next two down (20 up / 20 down per set). Saved to
   `D1_gene_sets.csv` and committed (96821c9) before any Phase 4 fit.
4. Injection: binomial thinning of raw integer counts with keep probability p = 2^(−δ). Up genes: thinned in all cells of
   the relevant study(ies) except responders. Down genes: thinned in responders. S in both studies, A in GSE216146 only,
   B in GSE271055 only. Responders: per study, `rng.choice` of round(rf × n_pseudo-treated) pseudo-treated cells
   (`default_rng(20261003)`, re-created per responder fraction; the same responder set is used by S and by A/B within a
   study and across δ). Thinning draws: `default_rng(20261003)` re-created per configuration, programs applied in the
   order S, A, B. Library sizes were reduced by the removed counts; X = log1p(count / library × 1e4) over all 18,271
   genes (max abs difference vs scanpy `normalize_total` + `log1p` on the uninjected base cells: 9.5e-7). Responder
   counts: `D1_responder_counts.csv` (rf 1.0: 2,415 / 4,727 responders; rf 0.3: 724 / 1,418; test-set pseudo-treated
   242 / 473). A data-level check of realized log2 fold changes is in `D1_injection_check.csv`.
5. Unit test: synthetic Poisson counts (100,000 cells × 30 genes), δ ∈ {0.25, 0.5, 1, 2}, 30 % responders; realized
   log2FC within 0.05 of +δ (up) and −δ (down); non-target genes, the non-target study, up genes in responders and down
   genes in non-responders are byte-identical; δ = 0 is the identity. Passes (run with the venv Python; pytest is not part
   of the pinned environment).
6. Scenarios: null (3 fits); shared (S), specific (A + B), mixed (S + A + B) × δ ∈ {0.25, 0.5, 1, 2} × responder fraction
   ∈ {1.0, 0.3} × seeds 0–2 (72 fits). Total 75 fits, all completed.
7. Fitting: canonical split on study × pseudo-label (80/10/10 stratified, ≤ 1,500 training cells per stratum, seed
   1729): 6,000 train / 1,428 validation / 1,429 test cells, identical for all configurations; full MC-ContrastiveVI with
   the canonical configuration and training loop (copy of `06_validation.train`, mode `full`); CPU, 3 threads per fit,
   3 fits concurrently.
8. Scoring per fit: held-out MSE; active units (Var[E q(z|x)] > 0.01) and per-dimension KL per block on pseudo-treated
   test cells; F_sh pooled and per drug (gated posterior means, pseudo-treated test cells); true shared fraction —
   primary energy fraction Σ‖e_S‖² / Σ(‖e_S‖² + ‖e_A‖² + ‖e_B‖²) and secondary trace-of-covariance fraction
   tr Cov(e_S) / (tr Cov(e_S) + tr Cov([e_A, e_B])), both over pseudo-treated test cells with e = ±δ·ln 2 per injected gene
   for responders and 0 otherwise (with every pseudo-treated cell a responder, e_S is constant and has zero variance, so the
   trace fraction is 0 for the specific and mixed scenarios and undefined (0/0, reported NA) for shared-only; it is undefined
   for the null; summed variances < 1e-12 are set to 0); cross-seed linear CKA per block on test cells; integrated gradients
   with Phase 2 settings (squared posterior-mean norm; zero and control-median baselines; `default_rng(314)` cells:
   64 + 64 pseudo-treated test cells for shared, 128 per drug block) with top-50 precision/recall of shared → S,
   cisplatin block → A, doxorubicin block → B, cross-routing counts, and the 4-of-6 consensus per configuration.
9. Timing and reductions: the first batch (null, seeds 0–2, 3 concurrent) took 197 s wall-clock including attribution;
   projected remaining Phase 4 time ≈ 4.8 h (< 48 h), so no reductions were applied (100 permutations; full grid).
10. Label-exchange reference: for null and shared_d1_rf1.0, the injected data were kept fixed and the pseudo-control /
    pseudo-treated labels permuted at the cell level within each library (one `default_rng(20261003)` stream; permutation
    k is the (k+1)-th pass over the sorted libraries), 100 permutations, one seed-0 fit each, with the canonical split
    held fixed; pooled F_sh computed on the permuted labels of the test cells. Observed statistic = the seed-0 fit of the
    unpermuted configuration; P = (1 + #{permuted F_sh ≥ observed}) / (1 + 100).
11. Calibration: Task 4 paired PyDESeq2 analysis (logic of `16_pseudobulk_paired.py`, outputs redirected to
    `runs/peerreview_20261003/pseudobulk_task4/`); rows with wald_padj < 0.05 and all three paired log2CPM differences in
    the direction of log2FoldChange (support rule of `19_gene_pattern_audit.py`).

## Per-configuration summary (`D1_config_summary.csv`; means over seeds 0–2)

| config | true shared (energy) | true shared (trace) | F_sh mean | F_sh SD | F_sh min–max | F_sh dox | F_sh cis | test MSE | AU sh/dox/cis | KL sh/dox/cis (nats) | CKA bg/sh/drug |
|---|---|---|---|---|---|---|---|---|---|---|---|
| null | NA | NA | 0.792 | 0.082 | 0.705–0.869 | 0.698 | 0.913 | 0.16339 | 0.33/0.00/0.00 | 0.0072/0.0010/0.0027 | 0.626/0.139/0.373 |
| shared_d0.25_rf1.0 | 1.000 | NA | 0.932 | 0.022 | 0.911–0.955 | 0.885 | 0.976 | 0.16291 | 0.33/0.00/0.00 | 0.0124/0.0007/0.0018 | 0.702/0.049/0.227 |
| shared_d0.25_rf0.3 | 1.000 | 1.000 | 0.748 | 0.314 | 0.389–0.973 | 0.584 | 0.915 | 0.16288 | 0.33/0.00/0.00 | 0.0117/0.0014/0.0026 | 0.710/0.093/0.558 |
| shared_d0.5_rf1.0 | 1.000 | NA | 0.775 | 0.064 | 0.717–0.843 | 0.560 | 0.909 | 0.16271 | 0.00/0.00/0.00 | 0.0051/0.0007/0.0020 | 0.703/0.274/0.371 |
| shared_d0.5_rf0.3 | 1.000 | 1.000 | 0.756 | 0.228 | 0.492–0.898 | 0.800 | 0.805 | 0.16268 | 0.00/0.00/0.00 | 0.0061/0.0004/0.0046 | 0.688/0.116/0.280 |
| shared_d1_rf1.0 | 1.000 | NA | 0.934 | 0.036 | 0.895–0.965 | 0.796 | 0.959 | 0.16249 | 0.33/0.00/0.00 | 0.0132/0.0006/0.0022 | 0.714/0.150/0.324 |
| shared_d1_rf0.3 | 1.000 | 1.000 | 0.833 | 0.049 | 0.804–0.890 | 0.641 | 0.941 | 0.16226 | 0.00/0.00/0.00 | 0.0058/0.0011/0.0013 | 0.723/0.291/0.299 |
| shared_d2_rf1.0 | 1.000 | NA | 0.898 | 0.041 | 0.872–0.944 | 0.796 | 0.959 | 0.16264 | 0.00/0.00/0.00 | 0.0129/0.0011/0.0019 | 0.699/0.093/0.152 |
| shared_d2_rf0.3 | 1.000 | 1.000 | 0.819 | 0.043 | 0.780–0.866 | 0.812 | 0.878 | 0.16218 | 0.00/0.00/0.00 | 0.0049/0.0006/0.0018 | 0.760/0.125/0.377 |
| specific_d0.25_rf1.0 | 0.000 | 0.000 | 0.684 | 0.203 | 0.474–0.879 | 0.630 | 0.771 | 0.16315 | 0.00/0.00/0.00 | 0.0058/0.0012/0.0023 | 0.655/0.218/0.388 |
| specific_d0.25_rf0.3 | 0.000 | 0.000 | 0.894 | 0.076 | 0.840–0.980 | 0.770 | 0.952 | 0.16293 | 1.00/0.00/0.00 | 0.0244/0.0013/0.0032 | 0.738/0.247/0.220 |
| specific_d0.5_rf1.0 | 0.000 | 0.000 | 0.929 | 0.025 | 0.910–0.957 | 0.910 | 0.948 | 0.16285 | 0.33/0.00/0.00 | 0.0085/0.0003/0.0017 | 0.688/0.230/0.257 |
| specific_d0.5_rf0.3 | 0.000 | 0.000 | 0.850 | 0.201 | 0.618–0.970 | 0.683 | 0.913 | 0.16269 | 0.33/0.00/0.00 | 0.0141/0.0010/0.0014 | 0.705/0.132/0.339 |
| specific_d1_rf1.0 | 0.000 | 0.000 | 0.919 | 0.055 | 0.858–0.966 | 0.785 | 0.969 | 0.16260 | 0.33/0.00/0.00 | 0.0150/0.0008/0.0017 | 0.702/0.102/0.364 |
| specific_d1_rf0.3 | 0.000 | 0.000 | 0.755 | 0.238 | 0.484–0.929 | 0.738 | 0.792 | 0.16246 | 0.00/0.00/0.00 | 0.0116/0.0009/0.0046 | 0.679/0.126/0.174 |
| specific_d2_rf1.0 | 0.000 | 0.000 | 0.889 | 0.030 | 0.856–0.913 | 0.771 | 0.928 | 0.16248 | 0.33/0.00/0.00 | 0.0148/0.0012/0.0034 | 0.746/0.257/0.108 |
| specific_d2_rf0.3 | 0.000 | 0.000 | 0.861 | 0.041 | 0.825–0.906 | 0.674 | 0.948 | 0.16205 | 0.00/0.00/0.00 | 0.0100/0.0011/0.0020 | 0.621/0.099/0.367 |
| mixed_d0.25_rf1.0 | 0.500 | 0.000 | 0.807 | 0.227 | 0.546–0.958 | 0.662 | 0.888 | 0.16274 | 0.00/0.00/0.00 | 0.0095/0.0009/0.0022 | 0.688/0.133/0.384 |
| mixed_d0.25_rf0.3 | 0.500 | 0.456 | 0.569 | 0.251 | 0.294–0.785 | 0.445 | 0.874 | 0.16261 | 0.00/0.00/0.00 | 0.0061/0.0028/0.0039 | 0.672/0.177/0.280 |
| mixed_d0.5_rf1.0 | 0.500 | 0.000 | 0.871 | 0.077 | 0.809–0.958 | 0.713 | 0.957 | 0.16241 | 0.00/0.00/0.00 | 0.0091/0.0006/0.0020 | 0.736/0.141/0.344 |
| mixed_d0.5_rf0.3 | 0.500 | 0.456 | 0.818 | 0.105 | 0.722–0.930 | 0.586 | 0.941 | 0.16254 | 0.67/0.00/0.00 | 0.0157/0.0017/0.0030 | 0.677/0.156/0.343 |
| mixed_d1_rf1.0 | 0.500 | 0.000 | 0.883 | 0.140 | 0.726–0.993 | 0.682 | 0.975 | 0.16195 | 1.33/0.00/0.00 | 0.0268/0.0014/0.0017 | 0.743/0.107/0.305 |
| mixed_d1_rf0.3 | 0.500 | 0.456 | 0.779 | 0.126 | 0.634–0.861 | 0.522 | 0.929 | 0.16187 | 0.00/0.00/0.00 | 0.0099/0.0021/0.0030 | 0.746/0.210/0.219 |
| mixed_d2_rf1.0 | 0.500 | 0.000 | 0.867 | 0.081 | 0.814–0.960 | 0.703 | 0.951 | 0.16173 | 0.00/0.00/0.00 | 0.0108/0.0013/0.0024 | 0.720/0.151/0.244 |
| mixed_d2_rf0.3 | 0.500 | 0.456 | 0.878 | 0.054 | 0.826–0.935 | 0.727 | 0.919 | 0.16135 | 0.00/0.00/0.00 | 0.0110/0.0007/0.0036 | 0.722/0.156/0.126 |

## Per-fit F_sh (`D1_fit_metrics.csv`)

| config | seed 0 | seed 1 | seed 2 |
|---|---|---|---|
| null | 0.7054 | 0.8688 | 0.8016 |
| shared_d0.25_rf1.0 | 0.9282 | 0.9555 | 0.9112 |
| shared_d0.25_rf0.3 | 0.3886 | 0.8823 | 0.9727 |
| shared_d0.5_rf1.0 | 0.7172 | 0.8433 | 0.7633 |
| shared_d0.5_rf0.3 | 0.8976 | 0.4923 | 0.8777 |
| shared_d1_rf1.0 | 0.9651 | 0.9409 | 0.8950 |
| shared_d1_rf0.3 | 0.8056 | 0.8041 | 0.8896 |
| shared_d2_rf1.0 | 0.8716 | 0.8766 | 0.9444 |
| shared_d2_rf0.3 | 0.8095 | 0.7804 | 0.8659 |
| specific_d0.25_rf1.0 | 0.6983 | 0.8790 | 0.4743 |
| specific_d0.25_rf0.3 | 0.9804 | 0.8397 | 0.8622 |
| specific_d0.5_rf1.0 | 0.9100 | 0.9205 | 0.9574 |
| specific_d0.5_rf0.3 | 0.6175 | 0.9703 | 0.9608 |
| specific_d1_rf1.0 | 0.8581 | 0.9337 | 0.9656 |
| specific_d1_rf0.3 | 0.8524 | 0.9294 | 0.4837 |
| specific_d2_rf1.0 | 0.8987 | 0.8560 | 0.9126 |
| specific_d2_rf0.3 | 0.8254 | 0.8503 | 0.9058 |
| mixed_d0.25_rf1.0 | 0.5461 | 0.9181 | 0.9578 |
| mixed_d0.25_rf0.3 | 0.6266 | 0.7852 | 0.2938 |
| mixed_d0.5_rf1.0 | 0.8469 | 0.8092 | 0.9583 |
| mixed_d0.5_rf0.3 | 0.8013 | 0.9301 | 0.7220 |
| mixed_d1_rf1.0 | 0.7257 | 0.9934 | 0.9302 |
| mixed_d1_rf0.3 | 0.8415 | 0.8613 | 0.6340 |
| mixed_d2_rf1.0 | 0.9600 | 0.8283 | 0.8138 |
| mixed_d2_rf0.3 | 0.8261 | 0.9347 | 0.8738 |

## Integrated-gradient recovery, per fit, mean over seeds (`D1_ig_per_fit.csv`, `D1_ig_per_config_mean.csv`)

Precision = hits/50, recall = hits/40 for the top 50 genes of each block. Cross-routing counts: A∪B genes in the shared block top 50; S genes in the cisplatin / doxorubicin block top 50.

| config | baseline | shared→S P/R | cis→A P/R | dox→B P/R | A∪B in shared top50 | S in cis top50 | S in dox top50 |
|---|---|---|---|---|---|---|---|
| mixed_d0.25_rf0.3 | control_median | 0.060/0.075 | 0.033/0.042 | 0.033/0.042 | 1.00 | 2.00 | 0.33 |
| mixed_d0.25_rf0.3 | zero | 0.040/0.050 | 0.033/0.042 | 0.020/0.025 | 0.67 | 1.33 | 1.00 |
| mixed_d0.25_rf1.0 | control_median | 0.067/0.083 | 0.027/0.033 | 0.027/0.033 | 2.00 | 2.33 | 0.33 |
| mixed_d0.25_rf1.0 | zero | 0.067/0.083 | 0.027/0.033 | 0.020/0.025 | 1.00 | 2.67 | 0.67 |
| mixed_d0.5_rf0.3 | control_median | 0.040/0.050 | 0.027/0.033 | 0.020/0.025 | 1.67 | 1.33 | 0.33 |
| mixed_d0.5_rf0.3 | zero | 0.033/0.042 | 0.027/0.033 | 0.020/0.025 | 2.33 | 2.00 | 1.00 |
| mixed_d0.5_rf1.0 | control_median | 0.033/0.042 | 0.027/0.033 | 0.027/0.033 | 1.67 | 3.00 | 0.33 |
| mixed_d0.5_rf1.0 | zero | 0.040/0.050 | 0.040/0.050 | 0.013/0.017 | 1.67 | 1.67 | 0.67 |
| mixed_d1_rf0.3 | control_median | 0.053/0.067 | 0.020/0.025 | 0.020/0.025 | 1.33 | 2.33 | 0.33 |
| mixed_d1_rf0.3 | zero | 0.053/0.067 | 0.027/0.033 | 0.013/0.017 | 1.67 | 1.33 | 1.00 |
| mixed_d1_rf1.0 | control_median | 0.027/0.033 | 0.040/0.050 | 0.020/0.025 | 1.33 | 1.00 | 0.33 |
| mixed_d1_rf1.0 | zero | 0.033/0.042 | 0.020/0.025 | 0.027/0.033 | 1.33 | 2.33 | 1.00 |
| mixed_d2_rf0.3 | control_median | 0.040/0.050 | 0.007/0.008 | 0.027/0.033 | 1.00 | 0.67 | 0.33 |
| mixed_d2_rf0.3 | zero | 0.040/0.050 | 0.007/0.008 | 0.033/0.042 | 1.33 | 0.67 | 0.33 |
| mixed_d2_rf1.0 | control_median | 0.020/0.025 | 0.040/0.050 | 0.033/0.042 | 3.00 | 2.33 | 1.00 |
| mixed_d2_rf1.0 | zero | 0.020/0.025 | 0.027/0.033 | 0.027/0.033 | 3.67 | 2.33 | 0.67 |
| null | control_median | 0.060/0.075 | 0.067/0.083 | 0.033/0.042 | 1.67 | 3.33 | 0.00 |
| null | zero | 0.060/0.075 | 0.040/0.050 | 0.027/0.033 | 1.33 | 2.00 | 0.67 |
| shared_d0.25_rf0.3 | control_median | 0.033/0.042 | 0.047/0.058 | 0.047/0.058 | 2.67 | 2.00 | 0.00 |
| shared_d0.25_rf0.3 | zero | 0.033/0.042 | 0.053/0.067 | 0.033/0.042 | 2.00 | 2.67 | 0.00 |
| shared_d0.25_rf1.0 | control_median | 0.047/0.058 | 0.053/0.067 | 0.020/0.025 | 2.67 | 1.67 | 0.67 |
| shared_d0.25_rf1.0 | zero | 0.047/0.058 | 0.047/0.058 | 0.013/0.017 | 2.67 | 1.33 | 0.33 |
| shared_d0.5_rf0.3 | control_median | 0.060/0.075 | 0.027/0.033 | 0.033/0.042 | 2.00 | 1.33 | 0.33 |
| shared_d0.5_rf0.3 | zero | 0.047/0.058 | 0.007/0.008 | 0.020/0.025 | 2.00 | 1.33 | 0.67 |
| shared_d0.5_rf1.0 | control_median | 0.053/0.067 | 0.033/0.042 | 0.033/0.042 | 1.33 | 2.00 | 0.00 |
| shared_d0.5_rf1.0 | zero | 0.040/0.050 | 0.040/0.050 | 0.013/0.017 | 1.33 | 2.00 | 0.33 |
| shared_d1_rf0.3 | control_median | 0.027/0.033 | 0.033/0.042 | 0.053/0.067 | 3.33 | 2.33 | 0.33 |
| shared_d1_rf0.3 | zero | 0.027/0.033 | 0.040/0.050 | 0.020/0.025 | 3.00 | 1.67 | 0.67 |
| shared_d1_rf1.0 | control_median | 0.047/0.058 | 0.020/0.025 | 0.013/0.017 | 2.67 | 2.33 | 1.00 |
| shared_d1_rf1.0 | zero | 0.033/0.042 | 0.033/0.042 | 0.027/0.033 | 2.00 | 2.00 | 1.33 |
| shared_d2_rf0.3 | control_median | 0.053/0.067 | 0.033/0.042 | 0.027/0.033 | 1.33 | 1.67 | 1.00 |
| shared_d2_rf0.3 | zero | 0.040/0.050 | 0.020/0.025 | 0.020/0.025 | 1.00 | 0.67 | 1.33 |
| shared_d2_rf1.0 | control_median | 0.013/0.017 | 0.033/0.042 | 0.013/0.017 | 3.33 | 2.33 | 1.33 |
| shared_d2_rf1.0 | zero | 0.027/0.033 | 0.027/0.033 | 0.020/0.025 | 2.67 | 1.00 | 0.67 |
| specific_d0.25_rf0.3 | control_median | 0.053/0.067 | 0.020/0.025 | 0.033/0.042 | 3.33 | 1.67 | 1.33 |
| specific_d0.25_rf0.3 | zero | 0.040/0.050 | 0.027/0.033 | 0.020/0.025 | 3.00 | 1.67 | 1.33 |
| specific_d0.25_rf1.0 | control_median | 0.027/0.033 | 0.040/0.050 | 0.033/0.042 | 1.67 | 3.33 | 0.00 |
| specific_d0.25_rf1.0 | zero | 0.040/0.050 | 0.047/0.058 | 0.027/0.033 | 1.67 | 2.67 | 1.00 |
| specific_d0.5_rf0.3 | control_median | 0.040/0.050 | 0.033/0.042 | 0.040/0.050 | 1.67 | 2.33 | 0.33 |
| specific_d0.5_rf0.3 | zero | 0.040/0.050 | 0.027/0.033 | 0.007/0.008 | 1.00 | 2.00 | 1.00 |
| specific_d0.5_rf1.0 | control_median | 0.060/0.075 | 0.027/0.033 | 0.027/0.033 | 1.67 | 1.67 | 0.33 |
| specific_d0.5_rf1.0 | zero | 0.060/0.075 | 0.027/0.033 | 0.027/0.033 | 2.33 | 2.33 | 1.00 |
| specific_d1_rf0.3 | control_median | 0.053/0.067 | 0.020/0.025 | 0.020/0.025 | 1.00 | 2.00 | 1.00 |
| specific_d1_rf0.3 | zero | 0.047/0.058 | 0.020/0.025 | 0.013/0.017 | 0.67 | 2.67 | 1.00 |
| specific_d1_rf1.0 | control_median | 0.040/0.050 | 0.013/0.017 | 0.020/0.025 | 1.00 | 2.00 | 0.67 |
| specific_d1_rf1.0 | zero | 0.033/0.042 | 0.020/0.025 | 0.007/0.008 | 1.00 | 1.00 | 1.33 |
| specific_d2_rf0.3 | control_median | 0.053/0.067 | 0.013/0.017 | 0.040/0.050 | 0.67 | 3.00 | 0.00 |
| specific_d2_rf0.3 | zero | 0.040/0.050 | 0.013/0.017 | 0.020/0.025 | 0.33 | 2.00 | 0.67 |
| specific_d2_rf1.0 | control_median | 0.047/0.058 | 0.033/0.042 | 0.047/0.058 | 2.33 | 3.00 | 0.67 |
| specific_d2_rf1.0 | zero | 0.053/0.067 | 0.020/0.025 | 0.007/0.008 | 2.33 | 3.00 | 0.67 |

## Integrated-gradient recovery, 4-of-6 consensus per configuration (`D1_ig_consensus.csv`)

Precision = hits / consensus-list size; recall = hits / 40.

| config | n shared/cis/dox | shared→S hits (P/R) | cis→A hits (P/R) | dox→B hits (P/R) | A∪B in shared | S in cis | S in dox |
|---|---|---|---|---|---|---|---|
| null | 28/22/24 | 2 (0.07/0.05) | 2 (0.09/0.05) | 1 (0.04/0.03) | 0 | 0 | 0 |
| shared_d0.25_rf1.0 | 23/16/20 | 1 (0.04/0.03) | 3 (0.19/0.07) | 0 (0.00/0.00) | 2 | 0 | 0 |
| shared_d0.25_rf0.3 | 25/27/21 | 0 (0.00/0.00) | 1 (0.04/0.03) | 2 (0.10/0.05) | 0 | 0 | 0 |
| shared_d0.5_rf1.0 | 28/13/26 | 2 (0.07/0.05) | 1 (0.08/0.03) | 1 (0.04/0.03) | 0 | 0 | 0 |
| shared_d0.5_rf0.3 | 28/22/25 | 2 (0.07/0.05) | 0 (0.00/0.00) | 1 (0.04/0.03) | 0 | 0 | 0 |
| shared_d1_rf1.0 | 28/15/20 | 1 (0.04/0.03) | 1 (0.07/0.03) | 0 (0.00/0.00) | 3 | 0 | 1 |
| shared_d1_rf0.3 | 24/18/24 | 1 (0.04/0.03) | 1 (0.06/0.03) | 0 (0.00/0.00) | 2 | 1 | 0 |
| shared_d2_rf1.0 | 27/11/23 | 0 (0.00/0.00) | 1 (0.09/0.03) | 0 (0.00/0.00) | 2 | 0 | 1 |
| shared_d2_rf0.3 | 22/23/19 | 1 (0.05/0.03) | 1 (0.04/0.03) | 0 (0.00/0.00) | 1 | 1 | 1 |
| specific_d0.25_rf1.0 | 25/26/27 | 1 (0.04/0.03) | 1 (0.04/0.03) | 1 (0.04/0.03) | 0 | 1 | 0 |
| specific_d0.25_rf0.3 | 26/21/28 | 2 (0.08/0.05) | 1 (0.05/0.03) | 0 (0.00/0.00) | 1 | 1 | 1 |
| specific_d0.5_rf1.0 | 28/19/22 | 2 (0.07/0.05) | 1 (0.05/0.03) | 0 (0.00/0.00) | 1 | 2 | 1 |
| specific_d0.5_rf0.3 | 25/21/27 | 0 (0.00/0.00) | 1 (0.05/0.03) | 1 (0.04/0.03) | 1 | 0 | 1 |
| specific_d1_rf1.0 | 19/19/24 | 1 (0.05/0.03) | 0 (0.00/0.00) | 0 (0.00/0.00) | 0 | 1 | 1 |
| specific_d1_rf0.3 | 28/23/26 | 2 (0.07/0.05) | 0 (0.00/0.00) | 0 (0.00/0.00) | 0 | 1 | 1 |
| specific_d2_rf1.0 | 29/14/21 | 2 (0.07/0.05) | 1 (0.07/0.03) | 0 (0.00/0.00) | 1 | 1 | 0 |
| specific_d2_rf0.3 | 33/22/21 | 2 (0.06/0.05) | 0 (0.00/0.00) | 1 (0.05/0.03) | 0 | 2 | 0 |
| mixed_d0.25_rf1.0 | 26/25/21 | 2 (0.08/0.05) | 1 (0.04/0.03) | 0 (0.00/0.00) | 0 | 1 | 0 |
| mixed_d0.25_rf0.3 | 22/20/25 | 2 (0.09/0.05) | 1 (0.05/0.03) | 1 (0.04/0.03) | 0 | 1 | 1 |
| mixed_d0.5_rf1.0 | 20/15/26 | 1 (0.05/0.03) | 1 (0.07/0.03) | 0 (0.00/0.00) | 0 | 0 | 0 |
| mixed_d0.5_rf0.3 | 23/14/26 | 1 (0.04/0.03) | 1 (0.07/0.03) | 0 (0.00/0.00) | 2 | 0 | 0 |
| mixed_d1_rf1.0 | 23/26/24 | 0 (0.00/0.00) | 1 (0.04/0.03) | 0 (0.00/0.00) | 0 | 0 | 1 |
| mixed_d1_rf0.3 | 36/20/17 | 3 (0.08/0.07) | 1 (0.05/0.03) | 0 (0.00/0.00) | 1 | 0 | 1 |
| mixed_d2_rf1.0 | 23/19/25 | 0 (0.00/0.00) | 1 (0.05/0.03) | 1 (0.04/0.03) | 1 | 2 | 0 |
| mixed_d2_rf0.3 | 28/37/30 | 1 (0.04/0.03) | 0 (0.00/0.00) | 1 (0.03/0.03) | 1 | 1 | 0 |

## Label-exchange reference (`D1_permutation_summary.csv`, `D1_permutation_values_*.csv`)

| config | observed F_sh (seed 0) | n perm | n perm ≥ observed | P | perm mean (SD) | perm 2.5% / 50% / 97.5% | perm min–max |
|---|---|---|---|---|---|---|---|
| null | 0.7054 | 100 | 77 | 0.7723 | 0.8065 (0.1809) | 0.3612 / 0.8685 / 0.9877 | 0.1317–0.9926 |
| shared_d1_rf1.0 | 0.9651 | 100 | 13 | 0.1386 | 0.8355 (0.1479) | 0.4539 / 0.8770 / 0.9813 | 0.2666–0.9892 |

## Calibration (`D1_calibration.json`)

```
{
 "n_supported_rows_all_genes": 292,
 "n_supported_unique_genes_all": 266,
 "median_abs_log2FC_supported_rows_all_genes": 0.7423859477427446,
 "n_supported_rows_in_1500_universe": 80,
 "n_supported_unique_genes_in_1500_universe": 69,
 "median_abs_log2FC_supported_rows_in_1500_universe": 0.9808300113324678,
 "median_of_per_gene_max_abs_log2FC_all": 0.7401869795822678,
 "median_of_per_gene_max_abs_log2FC_in_1500_universe": 0.9952179284199364,
 "quartiles_abs_log2FC_supported_rows_all_genes": [
  0.5692872936028686,
  0.7423859477427446,
  1.0764791266672646
 ],
 "delta_grid": [
  0.25,
  0.5,
  1.0,
  2.0
 ],
 "support_rule": "wald_padj < 0.05 and all three paired log2CPM differences share the sign of log2FoldChange",
 "pydeseq2_outputs": "runs\\peerreview_20261003\\pseudobulk_task4",
 "median_abs_log2FC_supported_rows_all_genes_grid_position": {
  "nearest_grid_value": 0.5,
  "between": [
   0.5,
   1.0
  ]
 },
 "median_abs_log2FC_supported_rows_in_1500_universe_grid_position": {
  "nearest_grid_value": 1.0,
  "between": [
   0.5,
   1.0
  ]
 }
}
```

## Injected gene sets (`D1_gene_sets.csv`)

| set | direction | genes (decile) |
|---|---|---|
| A | down | C3 (0), Twist1 (0), A830019P07Rik (1), 2010001K21Rik (1), Higd1b (2), Synpo2 (2), Dnah6 (3), Parp3 (3), Akna (4), Enpp6 (4), Baalc (5), AW047730 (5), Rnf182 (6), Psmb8 (6), Gm26532 (7), Sgms1 (7), Ctsz (8), Cacna2d3 (8), Snhg14 (9), Arpp21 (9) |
| A | up | Tcap (0), Folr2 (0), Crygn (1), Pla2g5 (1), Notch3 (2), Lrrc23 (2), Igfbp6 (3), Serpinb9 (3), Mirt1 (4), Gfod2 (4), G0s2 (5), G530011O06Rik (5), Phlda1 (6), Ets1 (6), Hist1h1e (7), Chst2 (7), Mobp (8), Arap2 (8), Fos (9), Tmem108 (9) |
| B | down | Gm20743 (0), Gm26737 (0), Cldn22 (1), Dnah11 (1), Apobr (2), Nts (2), Cped1 (3), Rnf213 (3), Hmcn1 (4), Tmem204 (4), Cpm (5), Nlrp3 (5), Plekha7 (6), Gad1 (6), Cited2 (7), Hspa1a (7), Ndrg2 (8), Icam5 (8), Syne1 (9), Ppfia2 (9) |
| B | up | Slc26a7 (0), Cyb5r2 (0), Lhx1 (1), Frmpd2 (1), Cfap43 (2), Col18a1 (2), Il34 (3), Fzd6 (3), Lims2 (4), Dlx1 (4), Slit2 (5), Psd2 (5), Mlc1 (6), Sorcs3 (6), Id2 (7), S100b (7), Ermn (8), Nrp1 (8), Oxr1 (9), Rims2 (9) |
| S | down | Cdc25c (0), 4933439K11Rik (0), Gm28729 (1), Sema3g (1), Mdfic (2), 1700008O03Rik (2), Pdlim1 (3), Edn1 (3), Acvrl1 (4), Hrk (4), Galntl6 (5), Rab33b (5), Sirt1 (6), Gsta4 (6), Klf4 (7), Plk5 (7), Adcy1 (8), Epb41l2 (8), Kcnq1ot1 (9), Sat1 (9) |
| S | up | Ccr2 (0), Pbk (0), Gm9946 (1), Defb11 (1), 9330185C12Rik (2), Gm26740 (2), Rsph1 (3), Ifi203 (3), Cd93 (4), Cfap54 (4), Grin2c (5), Atp13a4 (5), Mgat4c (6), Slc9a3r2 (6), Arl4a (7), Epha6 (7), Kif5a (8), Dcaf17 (8), Junb (9), Mal (9) |

## Diagnostic: injected signal in the inputs vs IG ranks, seed-0 fits (`D1_signal_check_seed0.csv`)

Added after inspecting the grid results, to check the injection pipeline. Input ranking: |mean(pseudo-treated) − mean(pseudo-control)| of the
log-normalized inputs over training cells (S: both studies; A: GSE216146; B: GSE271055), ranks 0–1,499. IG ranking: mean |IG|, zero baseline, in the
matching block (S → shared, A → cisplatin, B → doxorubicin). Probes: metrics-JSON condition probes (3-class) of the seed-0 fit.

| config | set | n in top 50 (input difference) | median rank (input) | n in top 50 (IG) | median rank (IG) | bg / shared-ungated / drug-ungated condition probe |
|---|---|---|---|---|---|---|
| null | S | 2 | 692 | 4 | 869 | 0.409 / 0.506 / 0.615 |
| null | A | 2 | 880 | 2 | 552 | 0.409 / 0.506 / 0.615 |
| null | B | 1 | 692 | 2 | 649 | 0.409 / 0.506 / 0.615 |
| shared_d0.25_rf1.0 | S | 6 | 735 | 3 | 752 | 0.511 / 0.487 / 0.600 |
| shared_d0.25_rf1.0 | A | 2 | 877 | 2 | 599 | 0.511 / 0.487 / 0.600 |
| shared_d0.25_rf1.0 | B | 1 | 694 | 1 | 675 | 0.511 / 0.487 / 0.600 |
| shared_d0.25_rf0.3 | S | 3 | 784 | 3 | 800 | 0.449 / 0.513 / 0.602 |
| shared_d0.25_rf0.3 | A | 2 | 878 | 1 | 678 | 0.449 / 0.513 / 0.602 |
| shared_d0.25_rf0.3 | B | 1 | 696 | 0 | 650 | 0.449 / 0.513 / 0.602 |
| shared_d0.5_rf1.0 | S | 10 | 340 | 2 | 812 | 0.436 / 0.527 / 0.600 |
| shared_d0.5_rf1.0 | A | 2 | 882 | 3 | 827 | 0.436 / 0.527 / 0.600 |
| shared_d0.5_rf1.0 | B | 1 | 698 | 0 | 657 | 0.436 / 0.527 / 0.600 |
| shared_d0.5_rf0.3 | S | 4 | 758 | 4 | 871 | 0.512 / 0.427 / 0.602 |
| shared_d0.5_rf0.3 | A | 2 | 876 | 0 | 756 | 0.512 / 0.427 / 0.602 |
| shared_d0.5_rf0.3 | B | 1 | 695 | 1 | 712 | 0.512 / 0.427 / 0.602 |
| shared_d1_rf1.0 | S | 18 | 109 | 2 | 914 | 0.462 / 0.425 / 0.598 |
| shared_d1_rf1.0 | A | 1 | 882 | 2 | 688 | 0.462 / 0.425 / 0.598 |
| shared_d1_rf1.0 | B | 0 | 700 | 2 | 658 | 0.462 / 0.425 / 0.598 |
| shared_d1_rf0.3 | S | 5 | 653 | 2 | 940 | 0.467 / 0.423 / 0.619 |
| shared_d1_rf0.3 | A | 2 | 880 | 0 | 666 | 0.467 / 0.423 / 0.619 |
| shared_d1_rf0.3 | B | 1 | 694 | 1 | 734 | 0.467 / 0.423 / 0.619 |
| shared_d2_rf1.0 | S | 20 | 53 | 0 | 945 | 0.432 / 0.457 / 0.595 |
| shared_d2_rf1.0 | A | 1 | 878 | 1 | 676 | 0.432 / 0.457 / 0.595 |
| shared_d2_rf1.0 | B | 0 | 704 | 0 | 788 | 0.432 / 0.457 / 0.595 |
| shared_d2_rf0.3 | S | 10 | 336 | 1 | 809 | 0.467 / 0.555 / 0.621 |
| shared_d2_rf0.3 | A | 2 | 876 | 0 | 639 | 0.467 / 0.555 / 0.621 |
| shared_d2_rf0.3 | B | 1 | 696 | 2 | 784 | 0.467 / 0.555 / 0.621 |
| specific_d0.25_rf1.0 | S | 1 | 693 | 3 | 884 | 0.492 / 0.424 / 0.617 |
| specific_d0.25_rf1.0 | A | 5 | 628 | 3 | 580 | 0.492 / 0.424 / 0.617 |
| specific_d0.25_rf1.0 | B | 7 | 744 | 0 | 686 | 0.492 / 0.424 / 0.617 |
| specific_d0.25_rf0.3 | S | 2 | 688 | 3 | 827 | 0.503 / 0.471 / 0.627 |
| specific_d0.25_rf0.3 | A | 2 | 968 | 1 | 550 | 0.503 / 0.471 / 0.627 |
| specific_d0.25_rf0.3 | B | 3 | 698 | 0 | 812 | 0.503 / 0.471 / 0.627 |
| specific_d0.5_rf1.0 | S | 1 | 703 | 2 | 768 | 0.489 / 0.440 / 0.609 |
| specific_d0.5_rf1.0 | A | 8 | 417 | 2 | 787 | 0.489 / 0.440 / 0.609 |
| specific_d0.5_rf1.0 | B | 9 | 550 | 2 | 712 | 0.489 / 0.440 / 0.609 |
| specific_d0.5_rf0.3 | S | 2 | 686 | 2 | 815 | 0.444 / 0.503 / 0.617 |
| specific_d0.5_rf0.3 | A | 4 | 683 | 1 | 605 | 0.444 / 0.503 / 0.617 |
| specific_d0.5_rf0.3 | B | 4 | 582 | 0 | 694 | 0.444 / 0.503 / 0.617 |
| specific_d1_rf1.0 | S | 1 | 712 | 2 | 830 | 0.527 / 0.497 / 0.597 |
| specific_d1_rf1.0 | A | 13 | 206 | 1 | 692 | 0.527 / 0.497 / 0.597 |
| specific_d1_rf1.0 | B | 12 | 404 | 1 | 712 | 0.527 / 0.497 / 0.597 |
| specific_d1_rf0.3 | S | 2 | 692 | 2 | 762 | 0.483 / 0.497 / 0.582 |
| specific_d1_rf0.3 | A | 6 | 662 | 1 | 804 | 0.483 / 0.497 / 0.582 |
| specific_d1_rf0.3 | B | 7 | 726 | 1 | 832 | 0.483 / 0.497 / 0.582 |
| specific_d2_rf1.0 | S | 1 | 713 | 4 | 885 | 0.470 / 0.541 / 0.607 |
| specific_d2_rf1.0 | A | 17 | 86 | 1 | 964 | 0.470 / 0.541 / 0.607 |
| specific_d2_rf1.0 | B | 14 | 211 | 1 | 862 | 0.470 / 0.541 / 0.607 |
| specific_d2_rf0.3 | S | 1 | 695 | 3 | 740 | 0.428 / 0.479 / 0.608 |
| specific_d2_rf0.3 | A | 8 | 484 | 0 | 866 | 0.428 / 0.479 / 0.608 |
| specific_d2_rf0.3 | B | 7 | 557 | 2 | 903 | 0.428 / 0.479 / 0.608 |
| mixed_d0.25_rf1.0 | S | 6 | 738 | 2 | 712 | 0.458 / 0.451 / 0.623 |
| mixed_d0.25_rf1.0 | A | 4 | 669 | 1 | 608 | 0.458 / 0.451 / 0.623 |
| mixed_d0.25_rf1.0 | B | 7 | 612 | 1 | 783 | 0.458 / 0.451 / 0.623 |
| mixed_d0.25_rf0.3 | S | 3 | 781 | 2 | 794 | 0.506 / 0.607 / 0.579 |
| mixed_d0.25_rf0.3 | A | 3 | 889 | 3 | 740 | 0.506 / 0.607 / 0.579 |
| mixed_d0.25_rf0.3 | B | 2 | 732 | 0 | 834 | 0.506 / 0.607 / 0.579 |
| mixed_d0.5_rf1.0 | S | 10 | 354 | 3 | 848 | 0.475 / 0.478 / 0.602 |
| mixed_d0.5_rf1.0 | A | 8 | 406 | 1 | 726 | 0.475 / 0.478 / 0.602 |
| mixed_d0.5_rf1.0 | B | 8 | 612 | 1 | 742 | 0.475 / 0.478 / 0.602 |
| mixed_d0.5_rf0.3 | S | 4 | 756 | 1 | 739 | 0.443 / 0.503 / 0.563 |
| mixed_d0.5_rf0.3 | A | 4 | 798 | 2 | 572 | 0.443 / 0.503 / 0.563 |
| mixed_d0.5_rf0.3 | B | 5 | 636 | 1 | 740 | 0.443 / 0.503 / 0.563 |
| mixed_d1_rf1.0 | S | 16 | 136 | 2 | 840 | 0.452 / 0.470 / 0.626 |
| mixed_d1_rf1.0 | A | 12 | 215 | 2 | 720 | 0.452 / 0.470 / 0.626 |
| mixed_d1_rf1.0 | B | 11 | 359 | 1 | 802 | 0.452 / 0.470 / 0.626 |
| mixed_d1_rf0.3 | S | 5 | 652 | 4 | 699 | 0.508 / 0.518 / 0.603 |
| mixed_d1_rf0.3 | A | 5 | 674 | 2 | 720 | 0.508 / 0.518 / 0.603 |
| mixed_d1_rf0.3 | B | 5 | 709 | 1 | 826 | 0.508 / 0.518 / 0.603 |
| mixed_d2_rf1.0 | S | 19 | 82 | 1 | 891 | 0.497 / 0.461 / 0.564 |
| mixed_d2_rf1.0 | A | 17 | 106 | 1 | 902 | 0.497 / 0.461 / 0.564 |
| mixed_d2_rf1.0 | B | 13 | 248 | 1 | 856 | 0.497 / 0.461 / 0.564 |
| mixed_d2_rf0.3 | S | 10 | 348 | 3 | 866 | 0.465 / 0.511 / 0.596 |
| mixed_d2_rf0.3 | A | 7 | 476 | 0 | 916 | 0.465 / 0.511 / 0.596 |
| mixed_d2_rf0.3 | B | 7 | 562 | 1 | 802 | 0.465 / 0.511 / 0.596 |

## Calibration summary

- Direction-consistent BH-significant cisplatin-vs-control rows (gene × cell type): 292 rows, 266 unique genes; median |log2FC| = 0.742 (quartiles 0.569 / 1.076); between δ = 0.5 and δ = 1 on the grid (nearest grid value 0.5).
- Restricted to the 1,500-gene universe: 80 rows, 69 unique genes; median |log2FC| = 0.981; between 0.5 and 1 (nearest grid value 1).
- Median over genes of the per-gene maximum |log2FC|: 0.740 (all), 0.995 (universe).
- Regenerated PyDESeq2 hit counts (`runs/peerreview_20261003/pseudobulk_task4/pseudobulk_task4_summary.csv`) equal the committed `pseudobulk_task4_summary.csv` for every cell type × contrast except Microglia cisplatin_vs_control (101 vs 109) and Oligodendrocyte cisplatin_vs_control (141 vs 137); the supported universe gene count is 69 (script 19 asserts 68 for the committed results). Rows used: `D1_calibration_supported_rows.csv`.
- The δ grid was not changed.

## Reference values

- Expected number of set genes in a random top-50 list from 1,500 genes: 50 × 40 / 1,500 = 1.33.

## Deviations and notes

- Pre-decided reductions: none applied (projection ≈ 4.8 h).
- Thinning random seed (not specified in the instructions): `default_rng(20261003)` re-created per configuration.
- Within-library 50/50 split uses floor(n/2) pseudo-treated cells per library.
- Permutation fits keep the canonical split fixed (not re-stratified on permuted labels) and do not compute IG.
- `scripts/35_d1_signal_check.py` is an added diagnostic, not part of the pre-specified scoring.
- The trace-fraction computation was corrected after the first scoring pass: floating-point variance of a constant e_S block (~1e-34) had produced 1.000 instead of 0/0 for shared-only configurations at responder fraction 1.0; summed variances < 1e-12 are now set to 0. No fit was affected.
- Runtime: grid 22:26:49–23:45 (75 fits; mean fit time 166 s); permutations 23:45–02:54:57 (200 fits).
