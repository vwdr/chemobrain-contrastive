# A1 — Fresh canonical full-model fits (seeds 0, 1, 2)

Script: `scripts/24_canonical_fits.py` (`--fit --seed S --threads=2`, then `--analyze`); helpers in
`src/peerreview/core.py` and `src/peerreview/fitjob.py`. Logs: `runs/peerreview_20261003/logs/24_*.log`.
Checkpoints, latents (all 46,857 cells incl. rescue; keys `bg, shared, drug, raw_shared, raw_drug, mse, nll`
plus `lv_shared, lv_drug`), histories, IG arrays: `runs/peerreview_20261003/canonical/`.

## Settings

- Inputs: regenerated canonical inputs (`runs/peerreview_20261003/canonical_inputs/inputs.npz`; content-identical
  to the committed inputs, see A0). 6,000 train / 2,881 validation / 2,881 test cells, 1,500 committed genes.
- Model and training: copy of `06_validation.train` mode `full` (hidden 128, dropout 0.1, latent 16/8/4+4, AdamW
  lr 1e-3 wd 1e-5, batch 256, ≤100 epochs, patience 15 on validation MSE, grad-clip 5, KL weight 1, HSIC
  weights 10/10/5, fast HSIC). CPU, `torch.set_num_threads(2)` (as in `06_validation.py`); three fits ran
  concurrently. Fit wall time 178–179 s each; IG 4–6 s.
- Diagnostics: F_sh (pooled/per drug; gated posterior means on treated test cells), active units
  (Var[E q(z|x)] > 0.01) and per-dimension KL on treated test cells (14_latent_usage logic), metrics-JSON probes
  (06), latent probes (11), cross-seed linear CKA on test cells (11).
- Integrated gradients: logic of `09_interpretation.py`; target = squared posterior-mean norm; baselines zero and
  control median (training controls; per-study controls for drug blocks); attribution cells drawn with
  `default_rng(314)` (64+64 treated test cells for shared; 128 doxorubicin; 128 cisplatin test cells); adaptive
  Gauss–Legendre (256 points sufficed in all 18 seed×axis×baseline cases; max median scaled completeness error
  0.0071). Ranking statistic: mean |IG| per gene. Consensus: top 50 in ≥ 4 of 6 seed×baseline rankings.

## S4 check

Held-out log-expression MSE: seed 0 = 0.168857, seed 1 = 0.168697, seed 2 = 0.168616 (all within
[0.165, 0.172]) → pass.

## Fresh vs committed metrics

Committed values: `runs/corrected_20260920/full_{seed}_metrics.json`, `full_{seed}_history.csv`,
`analysis/corrected_20260920/variance_robustness.csv`, `latent_usage_committed_summary.csv`, `latent_probes.csv`.
Full table: `A1_metric_comparison.csv`.

| metric | seed | committed | fresh | difference |
|---|---|---|---|---|
| best_epoch | 0 | 99 | 98 | -1 |
| best_epoch | 1 | 97 | 95 | -2 |
| best_epoch | 2 | 98 | 98 | +0 |
| test_mse | 0 | 0.168651 | 0.168857 | +0.0002059 |
| test_mse | 1 | 0.168801 | 0.168697 | -0.000104 |
| test_mse | 2 | 0.16868 | 0.168616 | -6.457e-05 |
| validation_mse | 0 | 0.166412 | 0.166632 | +0.0002203 |
| validation_mse | 1 | 0.166814 | 0.166664 | -0.00015 |
| validation_mse | 2 | 0.166695 | 0.166567 | -0.0001281 |
| test_plugin_nll | 0 | 181.08 | 184.389 | +3.309 |
| test_plugin_nll | 1 | 185.808 | 191.848 | +6.04 |
| test_plugin_nll | 2 | 185.557 | 184.565 | -0.9926 |
| shared_fraction | 0 | 0.920568 | 0.923721 | +0.003153 |
| shared_fraction | 1 | 0.957295 | 0.967862 | +0.01057 |
| shared_fraction | 2 | 0.922368 | 0.970478 | +0.04811 |
| shared_fraction_doxorubicin | 0 | 0.816703 | 0.907443 | +0.09074 |
| shared_fraction_doxorubicin | 1 | 0.922526 | 0.893718 | -0.02881 |
| shared_fraction_doxorubicin | 2 | 0.880545 | 0.973406 | +0.09286 |
| shared_fraction_cisplatin | 0 | 0.958817 | 0.941863 | -0.01695 |
| shared_fraction_cisplatin | 1 | 0.97556 | 0.983179 | +0.00762 |
| shared_fraction_cisplatin | 2 | 0.941598 | 0.986146 | +0.04455 |
| shared_active_units_var_gt_0p01 | 0 | 0 | 0 | +0 |
| shared_active_units_var_gt_0p01 | 1 | 0 | 0 | +0 |
| shared_active_units_var_gt_0p01 | 2 | 1 | 0 | -1 |
| dox_active_units_var_gt_0p01 | 0 | 0 | 0 | +0 |
| dox_active_units_var_gt_0p01 | 1 | 0 | 0 | +0 |
| dox_active_units_var_gt_0p01 | 2 | 0 | 0 | +0 |
| cis_active_units_var_gt_0p01 | 0 | 0 | 0 | +0 |
| cis_active_units_var_gt_0p01 | 1 | 0 | 0 | +0 |
| cis_active_units_var_gt_0p01 | 2 | 0 | 0 | +0 |
| shared_dims_kl_gt_0p01 | 0 | 0 | 0 | +0 |
| shared_dims_kl_gt_0p01 | 1 | 0 | 0 | +0 |
| shared_dims_kl_gt_0p01 | 2 | 0 | 0 | +0 |
| dox_dims_kl_gt_0p01 | 0 | 0 | 0 | +0 |
| dox_dims_kl_gt_0p01 | 1 | 0 | 0 | +0 |
| dox_dims_kl_gt_0p01 | 2 | 0 | 0 | +0 |
| cis_dims_kl_gt_0p01 | 0 | 0 | 0 | +0 |
| cis_dims_kl_gt_0p01 | 1 | 0 | 0 | +0 |
| cis_dims_kl_gt_0p01 | 2 | 0 | 0 | +0 |
| shared_mean_total_kl_nats | 0 | 0.0106912 | 0.0124993 | +0.001808 |
| shared_mean_total_kl_nats | 1 | 0.00790796 | 0.0109293 | +0.003021 |
| shared_mean_total_kl_nats | 2 | 0.0218533 | 0.0182141 | -0.003639 |
| dox_mean_total_kl_nats | 0 | 0.000408167 | 0.000300347 | -0.0001078 |
| dox_mean_total_kl_nats | 1 | 0.000590962 | 0.000353184 | -0.0002378 |
| dox_mean_total_kl_nats | 2 | 0.00117188 | 0.000488569 | -0.0006833 |
| cis_mean_total_kl_nats | 0 | 0.00309098 | 0.00338505 | +0.0002941 |
| cis_mean_total_kl_nats | 1 | 0.00136415 | 0.00122126 | -0.0001429 |
| cis_mean_total_kl_nats | 2 | 0.00677856 | 0.00279439 | -0.003984 |
| bg_study_balanced_accuracy | 0 | 0.777189 | 0.781957 | +0.004768 |
| bg_study_balanced_accuracy | 1 | 0.856462 | 0.760458 | -0.096 |
| bg_study_balanced_accuracy | 2 | 0.726082 | 0.770911 | +0.04483 |
| bg_condition_balanced_accuracy | 0 | 0.490264 | 0.483899 | -0.006365 |
| bg_condition_balanced_accuracy | 1 | 0.576802 | 0.521838 | -0.05496 |
| bg_condition_balanced_accuracy | 2 | 0.496888 | 0.51612 | +0.01923 |
| shared_ungated_study_balanced_accuracy | 0 | 0.774299 | 0.850741 | +0.07644 |
| shared_ungated_study_balanced_accuracy | 1 | 0.73084 | 0.811155 | +0.08031 |
| shared_ungated_study_balanced_accuracy | 2 | 0.656704 | 0.597438 | -0.05927 |
| shared_ungated_condition_balanced_accuracy | 0 | 0.548508 | 0.554493 | +0.005986 |
| shared_ungated_condition_balanced_accuracy | 1 | 0.496302 | 0.572563 | +0.07626 |
| shared_ungated_condition_balanced_accuracy | 2 | 0.466458 | 0.378951 | -0.08751 |
| drug_gated_study_balanced_accuracy | 0 | 0.721799 | 0.701717 | -0.02008 |
| drug_gated_study_balanced_accuracy | 1 | 0.690786 | 0.758364 | +0.06758 |
| drug_gated_study_balanced_accuracy | 2 | 0.720272 | 0.725368 | +0.005096 |
| drug_gated_condition_balanced_accuracy | 0 | 0.958527 | 0.955013 | -0.003514 |
| drug_gated_condition_balanced_accuracy | 1 | 0.934556 | 0.980482 | +0.04593 |
| drug_gated_condition_balanced_accuracy | 2 | 0.988333 | 0.990013 | +0.001679 |
| drug_ungated_study_balanced_accuracy | 0 | 0.842052 | 0.895057 | +0.053 |
| drug_ungated_study_balanced_accuracy | 1 | 0.903235 | 0.903854 | +0.0006192 |
| drug_ungated_study_balanced_accuracy | 2 | 0.782282 | 0.890795 | +0.1085 |
| drug_ungated_condition_balanced_accuracy | 0 | 0.590618 | 0.62399 | +0.03337 |
| drug_ungated_condition_balanced_accuracy | 1 | 0.603655 | 0.615329 | +0.01167 |
| drug_ungated_condition_balanced_accuracy | 2 | 0.528004 | 0.551256 | +0.02325 |
| probe:background:cis_source_cell_type | 0 | 0.990842 | 0.993572 | +0.00273 |
| probe:background:cis_source_cell_type | 1 | 0.985378 | 0.985378 | +0 |
| probe:background:cis_source_cell_type | 2 | 0.990486 | 0.987637 | -0.002849 |
| probe:background:doxorubicin_treatment | 0 | 0.531156 | 0.527263 | -0.003893 |
| probe:background:doxorubicin_treatment | 1 | 0.529259 | 0.542788 | +0.01353 |
| probe:background:doxorubicin_treatment | 2 | 0.541841 | 0.537412 | -0.004429 |
| probe:background:cisplatin_treatment | 0 | 0.645093 | 0.670797 | +0.0257 |
| probe:background:cisplatin_treatment | 1 | 0.621498 | 0.6356 | +0.0141 |
| probe:background:cisplatin_treatment | 2 | 0.637164 | 0.65779 | +0.02063 |
| probe:shared_ungated:cis_source_cell_type | 0 | 0.842083 | 0.788336 | -0.05375 |
| probe:shared_ungated:cis_source_cell_type | 1 | 0.755592 | 0.829901 | +0.07431 |
| probe:shared_ungated:cis_source_cell_type | 2 | 0.751295 | 0.74651 | -0.004785 |
| probe:shared_ungated:doxorubicin_treatment | 0 | 0.519155 | 0.523862 | +0.004707 |
| probe:shared_ungated:doxorubicin_treatment | 1 | 0.483603 | 0.558681 | +0.07508 |
| probe:shared_ungated:doxorubicin_treatment | 2 | 0.494952 | 0.549836 | +0.05488 |
| probe:shared_ungated:cisplatin_treatment | 0 | 0.593706 | 0.573841 | -0.01987 |
| probe:shared_ungated:cisplatin_treatment | 1 | 0.641227 | 0.626615 | -0.01461 |
| probe:shared_ungated:cisplatin_treatment | 2 | 0.636869 | 0.494125 | -0.1427 |
| probe:drug_ungated:cis_source_cell_type | 0 | 0.820508 | 0.904442 | +0.08393 |
| probe:drug_ungated:cis_source_cell_type | 1 | 0.797712 | 0.789126 | -0.008586 |
| probe:drug_ungated:cis_source_cell_type | 2 | 0.754309 | 0.768092 | +0.01378 |
| probe:drug_ungated:doxorubicin_treatment | 0 | 0.531839 | 0.544214 | +0.01237 |
| probe:drug_ungated:doxorubicin_treatment | 1 | 0.535813 | 0.55228 | +0.01647 |
| probe:drug_ungated:doxorubicin_treatment | 2 | 0.498235 | 0.564744 | +0.06651 |
| probe:drug_ungated:cisplatin_treatment | 0 | 0.622459 | 0.514754 | -0.1077 |
| probe:drug_ungated:cisplatin_treatment | 1 | 0.60654 | 0.50411 | -0.1024 |
| probe:drug_ungated:cisplatin_treatment | 2 | 0.636793 | 0.572045 | -0.06475 |
## Cross-seed linear CKA (test cells) — `A1_cross_seed_cka.csv`

| seeds | block | committed | fresh |
|---|---|---|---|
| 0–1 | bg | 0.5578 | 0.5997 |
| 0–1 | shared | 0.0818 | 0.6910 |
| 0–1 | drug | 0.3892 | 0.1770 |
| 0–2 | bg | 0.6306 | 0.6316 |
| 0–2 | shared | 0.3586 | 0.5474 |
| 0–2 | drug | 0.1146 | 0.2579 |
| 1–2 | bg | 0.7084 | 0.7004 |
| 1–2 | shared | 0.5247 | 0.4972 |
| 1–2 | drug | 0.2175 | 0.6968 |

## Fresh consensus lists vs the frozen lists in `scripts/19_gene_pattern_audit.py`

Files: `A1_fresh_consensus_genes.csv`, `A1_fresh_consensus.json`, `A1_consensus_jaccard.csv`.

| block | n fresh | n frozen | overlap | Jaccard |
|---|---|---|---|---|
| shared | 40 | 34 | 29 | 0.644 |
| doxorubicin | 27 | 20 | 8 | 0.205 |
| cisplatin | 31 | 33 | 17 | 0.362 |

- shared — fresh only: Cflar, Ctsc, Ctss, Ly6e, Lyz2, Maf, Mylip, Pik3ap1, Sat1, Sparcl1, Tyrobp; frozen only: C5ar1, Ccl4, Cpne2, Ptprc, Qk.
- doxorubicin — overlap: Kcnq1ot1, Meg3, Mir100hg, Opcml, Peg3, Snhg11, Ttr, Zeb2; fresh only: Apoe, Arpp21, Atp1b1, Bc1, Camk2a, Dlgap1, Grm5, Nrgn, Nrxn1, Ntng1, Pclo, Ptprn, R3hdm1, Slc17a7, Slc1a2, Snap25, Sparcl1, Sptbn1, Ttyh1; frozen only: C1ql3, Dclk1, Gria2, Grin2b, Miat, Nfib, Nrcam, Pcdh9, Pcsk2, Ptprd, Rtn1, Zbtb20.
- cisplatin — overlap: Adamtsl3, Apoe, Bsg, Ccn3, Col25a1, Col6a2, Fos, Foxc2, Foxd1, Igfbp6, Mgp, Mrc1, Mt1, Ogn, Rspo3, Ttr, mt-Nd4; fresh only: Chchd10, Ctsd, Glul, Igfbp5, Junb, Klf4, Mt3, Myl9, Nbl1, Qk, Slc26a7, Slc38a2, Zbtb20, mt-Co2; frozen only: Anxa1, Cbr2, Cd163, Col1a2, Col6a1, Efemp1, F13a1, Folr2, Hes1, Id3, Ildr2, Mbp, Pf4, S100a10, Tmsb10, Tsc22d1.

Pairwise top-50 Jaccard among the six fresh seed×baseline rankings (`A1_fresh_ranking_top50_jaccard.csv`):
shared mean 0.502 (range 0.370–0.818); doxorubicin mean 0.273 (0.124–0.471); cisplatin mean 0.312 (0.075–0.754).
The committed per-seed IG table (`complete_corrected_attributions.csv`) is git-ignored and not present in this
checkout, so per-ranking fresh-vs-committed comparisons were not possible.

## Fresh seed 0 vs `ved_archive/full_0.pt` — `A1_seed0_vs_ved.csv`

Both evaluated on the regenerated canonical inputs.

| quantity | ved full_0 | fresh seed 0 |
|---|---|---|
| test MSE | 0.168553 | 0.168857 |
| validation MSE | 0.166477 | 0.166632 |
| test plug-in NLL | 187.054 | 184.389 |
| pooled F_sh | 0.7542 | 0.9237 |
| doxorubicin F_sh | 0.9440 | 0.9074 |
| cisplatin F_sh | 0.7320 | 0.9419 |
| active units shared/dox/cis | 0/0/0 | 0/0/0 |
| total KL shared/dox/cis (nats) | 0.00379/0.00019/0.00302 | 0.01250/0.00030/0.00339 |

Linear CKA between ved full_0 and fresh seed 0 on test cells: bg 0.690; shared (gated) 0.188; drug (gated) 0.156;
ungated shared 0.116; ungated drug 0.192.

## Deviations

- Thread count 2 per fit (canonical setting) rather than the run's 3-thread limit, to match `06_validation.py`.
- Committed fits ran on different hardware/software; fresh fits are not bitwise reproductions (best epochs
  98/95/98 vs 99/97/98).
