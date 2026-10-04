# A0 — Content check of the archived collaborator files (`ved_archive/`)

Script: `scripts/23_prepare_and_ved_check.py` (log: `runs/peerreview_20261003/logs/23_ved_check.log`).
The archive was opened read-only; nothing in it was modified, moved, or copied into the repository.

## Method

1. Canonical inputs were regenerated from `data/processed/recovered_counts.h5ad` (SHA-256
   `62a5047e…8ac1`, identical to the manifest) using the logic of `scripts/06_validation.py --prepare`
   (split seed 1729; 80/10/10 stratified by study × drug among nonrescue cells; 1,500 training cells per
   stratum; normalize_total 1e4 over all 18,271 genes; log1p). The modeled genes were taken from
   `analysis/corrected_20260920/benchmark_gene_universe.csv` (in its order). Outputs:
   `runs/peerreview_20261003/canonical_inputs/{inputs.npz,count_inputs.npz,cells.csv}`.
2. `ved_archive/full_0.pt` was loaded and evaluated on (a) `ved_archive/inputs.npz` and (b) the
   regenerated inputs with the metric definitions of `06_validation.py` (gated posterior means; test MSE,
   plug-in NLL, pooled F_sh on treated test cells, metrics-JSON probes). Validation MSE was compared with
   `runs/corrected_20260920/full_0_history.csv` at the committed best epoch.
3. Latents were recomputed from `full_0.pt` and compared with `ved_archive/full_0_latents.npz`.
   Seed-0 rows of `latent_probes.csv` were recomputed with the logic of `11_diagnostics.py`.
4. `ved_archive/cells.csv` was compared column by column with the regenerated cell metadata.

Hashes of the archived files (none match `task1_committed_checkpoint_hashes.json`):

| file | SHA-256 |
|---|---|
| cells.csv | 9a1f064d15a78aa5cfd70b3f6650bff0eb09292ccae4145c4679c83f6cae50ea |
| full_0.pt | 5b983b2d0fbfba70d563f836fcaef0127dff4c9b613bd7fec8a1f9764a3927b4 |
| full_0_latents.npz | 38e3eca71d8aca52acd9fd06c305bc752a3e318d60e66286c6aec5fd261989be |
| inputs.npz | 06918885b1db72bd506b15335fa4954dc607e285d3a4c2b4088fcedddf3f8242 |

## Classification

| file | classification | evidence |
|---|---|---|
| `inputs.npz` | content-identical to the committed run (within tolerance) | `train`, `val`, `test`, `d`, `b`, `genes` exactly equal; `X` max abs diff 9.54e-7 (382,362 of 70.3M entries differ at float32 rounding level) |
| `cells.csv` | content-identical to the committed run (within tolerance) | 46,857 rows, identical index and columns; 0 mismatches in every string column; max abs diff 0 in every numeric column |
| `full_0.pt` | differs | checkpoint `best_epoch` = 97 (committed metrics JSON: 99); test plug-in NLL 187.054 vs 181.080; pooled F_sh 0.7542 vs 0.9206; probe accuracies differ by up to 0.22 (table below) |
| `full_0_latents.npz` | differs from the committed run; consistent with `ved_archive/full_0.pt` | latents recomputed from `full_0.pt` match the archived file (max abs diff ≤ 1.4e-6 for latent blocks, 1.5e-7 for per-cell MSE, 1.2e-3 for per-cell NLL of magnitude up to 4,226); seed-0 probes from these latents differ from committed `latent_probes.csv` by up to 0.092 |

### `full_0.pt` metrics vs `runs/corrected_20260920/full_0_metrics.json`

Evaluated on `ved_archive/inputs.npz` (results on the regenerated inputs agree to ≤ 2.2e-3 in every probe and ≤ 6e-8 in F_sh; see `A0_full0_metrics_comparison.csv`).

| metric | committed | recomputed from ved `full_0.pt` | difference |
|---|---|---|---|
| best_epoch | 99 | 97 | −2 |
| test_mse | 0.16865083575248718 | 0.16855283081531525 | −9.800e-05 |
| validation_mse (history at best epoch) | 0.16641180362721270 | 0.16647653281688690 | +6.473e-05 |
| test_plugin_nll | 181.07989501953125 | 187.05416870117188 | +5.974 |
| shared_fraction (pooled F_sh) | 0.92056787014007568 | 0.75418812036514282 | −0.16638 |
| bg_study_balanced_accuracy | 0.77718913704191173 | 0.73026877387013278 | −0.04692 |
| bg_condition_balanced_accuracy | 0.49026382927960649 | 0.48197639028950379 | −0.00829 |
| shared_ungated_study_balanced_accuracy | 0.77429921427656423 | 0.55551984600682225 | −0.21878 |
| shared_ungated_condition_balanced_accuracy | 0.54850774868933494 | 0.45900316723740797 | −0.08950 |
| drug_gated_study_balanced_accuracy | 0.72179880900152704 | 0.72065355276000798 | −0.00115 |
| drug_gated_condition_balanced_accuracy | 0.95852661596958189 | 0.96810202788339661 | +0.00958 |
| drug_ungated_study_balanced_accuracy | 0.84205213140660473 | 0.91059887747316970 | +0.06855 |
| drug_ungated_condition_balanced_accuracy | 0.59061774876915885 | 0.59876563108380543 | +0.00815 |

### Seed-0 `latent_probes.csv` rows (11_diagnostics logic) from ved `full_0.pt`

| representation | task | committed | recomputed | difference |
|---|---|---|---|---|
| background | cis_source_cell_type | 0.990842 | 0.988794 | −0.00205 |
| background | doxorubicin_treatment | 0.531156 | 0.539784 | +0.00863 |
| background | cisplatin_treatment | 0.645093 | 0.666832 | +0.02174 |
| shared_ungated | cis_source_cell_type | 0.842083 | 0.749643 | −0.09244 |
| shared_ungated | doxorubicin_treatment | 0.519155 | 0.539788 | +0.02063 |
| shared_ungated | cisplatin_treatment | 0.593706 | 0.602531 | +0.00883 |
| drug_ungated | cis_source_cell_type | 0.820508 | 0.744458 | −0.07605 |
| drug_ungated | doxorubicin_treatment | 0.531839 | 0.553125 | +0.02129 |
| drug_ungated | cisplatin_treatment | 0.622459 | 0.567629 | −0.05483 |

## Additional observation (regeneration)

Re-running the seurat_v3 HVG selection of `06_validation.prepare` on this machine (scanpy 1.11.5,
scikit-misc 0.5.2) returns 1,499 of the 1,500 committed genes: `Cebpd` is selected instead of
`Aspm`. As specified, all new analyses use the committed `benchmark_gene_universe.csv`.
(`runs/peerreview_20261003/canonical_inputs/hvg_reselection_check.json`)

## Output tables

`A0_inputs_comparison.csv`, `A0_full0_metrics_comparison.csv`, `A0_full0_latents_comparison.csv`,
`A0_latent_probes_seed0_comparison.csv`, `A0_cells_comparison.csv`.
