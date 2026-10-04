# R1 — Location of the injected semi-synthetic signal

Script: `scripts/38_r1_injected_signal.py` (log `runs/peerreview_20261004/logs/38_r1.log`). Inputs (read-only): the 75 Phase 4
grid fits in `runs/peerreview_20261003/semisynthetic/<config>/` (latents, checkpoints), `base.npz`, `D1_gene_sets.csv`;
inputs X rebuilt with `build_inputs` of `scripts/31_semisynthetic.py` (deterministic). No refits.

## Methods

- Probes: `LogisticRegression(class_weight='balanced', max_iter=800)`, Phase 4 training cells → test cells, within
  GSE216146 and within GSE271055 separately, on `bg` (16 dims), `raw_shared` (8, ungated) and `raw_drug` (8, ungated). Labels:
  (a) pseudo-treated vs pseudo-control; (b) responders vs all other cells of the study (responder set of the configuration;
  identical to (a) at responder fraction 1.0; for the null configuration (b) was computed with both responder sets).
- Reconstruction: checkpoint decoded on test cells with canonical gated decoding (gated posterior means, study one-hot).
  For each set × study × direction: observed = mean over the set's genes of [mean X over test responders − mean X over the
  other test cells of that study]; reconstructed = same on the decoder output; ratio = reconstructed / observed (per fit,
  then averaged). Sets not injected in a configuration are reported with `injected = False` (reference). `expected_effect_ln`
  = ±δ·ln 2 is the per-gene shift in log counts before log1p normalization, listed for reference only.

Files: `R1_probes_per_fit.csv` (918 rows), `R1_reconstruction_per_fit.csv`, `R1_summary.csv` (both summaries stacked),
`R1_summary_probes.csv`, `R1_summary_reconstruction.csv` (mean, SD, min, max over seeds 0–2).

## Probes, label (a) pseudo-treated vs pseudo-control (mean balanced accuracy over seeds)

| config | GSE216146 bg | GSE216146 raw_shared | GSE216146 raw_drug | GSE271055 bg | GSE271055 raw_shared | GSE271055 raw_drug |
|---|---|---|---|---|---|---|
| null | 0.492 | 0.498 | 0.513 | 0.483 | 0.496 | 0.499 |
| shared_d0.25_rf1.0 | 0.498 | 0.489 | 0.502 | 0.487 | 0.497 | 0.505 |
| shared_d0.25_rf0.3 | 0.503 | 0.513 | 0.511 | 0.488 | 0.505 | 0.484 |
| shared_d0.5_rf1.0 | 0.487 | 0.511 | 0.503 | 0.482 | 0.509 | 0.508 |
| shared_d0.5_rf0.3 | 0.506 | 0.497 | 0.511 | 0.485 | 0.495 | 0.495 |
| shared_d1_rf1.0 | 0.499 | 0.503 | 0.511 | 0.487 | 0.496 | 0.489 |
| shared_d1_rf0.3 | 0.495 | 0.513 | 0.520 | 0.489 | 0.501 | 0.509 |
| shared_d2_rf1.0 | 0.504 | 0.507 | 0.512 | 0.489 | 0.506 | 0.523 |
| shared_d2_rf0.3 | 0.493 | 0.511 | 0.504 | 0.481 | 0.504 | 0.513 |
| specific_d0.25_rf1.0 | 0.506 | 0.513 | 0.519 | 0.485 | 0.495 | 0.506 |
| specific_d0.25_rf0.3 | 0.499 | 0.499 | 0.507 | 0.482 | 0.513 | 0.506 |
| specific_d0.5_rf1.0 | 0.499 | 0.493 | 0.514 | 0.478 | 0.514 | 0.502 |
| specific_d0.5_rf0.3 | 0.506 | 0.529 | 0.515 | 0.473 | 0.505 | 0.509 |
| specific_d1_rf1.0 | 0.498 | 0.514 | 0.518 | 0.476 | 0.504 | 0.501 |
| specific_d1_rf0.3 | 0.505 | 0.511 | 0.519 | 0.477 | 0.499 | 0.501 |
| specific_d2_rf1.0 | 0.516 | 0.518 | 0.509 | 0.496 | 0.525 | 0.511 |
| specific_d2_rf0.3 | 0.501 | 0.506 | 0.509 | 0.486 | 0.499 | 0.501 |
| mixed_d0.25_rf1.0 | 0.502 | 0.500 | 0.507 | 0.488 | 0.486 | 0.507 |
| mixed_d0.25_rf0.3 | 0.496 | 0.502 | 0.517 | 0.482 | 0.492 | 0.496 |
| mixed_d0.5_rf1.0 | 0.497 | 0.518 | 0.508 | 0.481 | 0.506 | 0.511 |
| mixed_d0.5_rf0.3 | 0.498 | 0.494 | 0.502 | 0.485 | 0.522 | 0.502 |
| mixed_d1_rf1.0 | 0.511 | 0.499 | 0.512 | 0.481 | 0.509 | 0.493 |
| mixed_d1_rf0.3 | 0.497 | 0.498 | 0.513 | 0.479 | 0.499 | 0.499 |
| mixed_d2_rf1.0 | 0.527 | 0.541 | 0.537 | 0.495 | 0.513 | 0.525 |
| mixed_d2_rf0.3 | 0.505 | 0.524 | 0.515 | 0.482 | 0.487 | 0.497 |

## Probes, label (b) responders vs other cells of the study, responder fraction 0.3 configurations and null

| config | responder set | GSE216146 bg | GSE216146 raw_shared | GSE216146 raw_drug | GSE271055 bg | GSE271055 raw_shared | GSE271055 raw_drug |
|---|---|---|---|---|---|---|---|
| null | 1 | 0.492 | 0.498 | 0.513 | 0.483 | 0.496 | 0.499 |
| null | 0.3 | 0.456 | 0.425 | 0.476 | 0.512 | 0.509 | 0.497 |
| shared_d0.25_rf0.3 | 0.3 | 0.455 | 0.511 | 0.484 | 0.505 | 0.507 | 0.508 |
| shared_d0.5_rf0.3 | 0.3 | 0.468 | 0.473 | 0.437 | 0.523 | 0.528 | 0.494 |
| shared_d1_rf0.3 | 0.3 | 0.457 | 0.500 | 0.446 | 0.513 | 0.510 | 0.486 |
| shared_d2_rf0.3 | 0.3 | 0.488 | 0.466 | 0.470 | 0.500 | 0.500 | 0.507 |
| specific_d0.25_rf0.3 | 0.3 | 0.455 | 0.488 | 0.467 | 0.516 | 0.510 | 0.533 |
| specific_d0.5_rf0.3 | 0.3 | 0.478 | 0.481 | 0.482 | 0.524 | 0.522 | 0.503 |
| specific_d1_rf0.3 | 0.3 | 0.465 | 0.477 | 0.521 | 0.513 | 0.530 | 0.493 |
| specific_d2_rf0.3 | 0.3 | 0.482 | 0.548 | 0.484 | 0.519 | 0.488 | 0.521 |
| mixed_d0.25_rf0.3 | 0.3 | 0.456 | 0.469 | 0.464 | 0.523 | 0.509 | 0.488 |
| mixed_d0.5_rf0.3 | 0.3 | 0.472 | 0.521 | 0.443 | 0.502 | 0.498 | 0.501 |
| mixed_d1_rf0.3 | 0.3 | 0.485 | 0.533 | 0.450 | 0.510 | 0.517 | 0.501 |
| mixed_d2_rf0.3 | 0.3 | 0.472 | 0.523 | 0.456 | 0.509 | 0.505 | 0.542 |

Range of all label-(a) mean balanced accuracies over configurations, studies and representations: 0.473–0.541; label (b): 0.425–0.548.

## Reconstructed effect, injected sets (mean over seeds; log1p-normalized units)

| config | set | study | direction | observed | reconstructed | ratio (mean of per-fit ratios) |
|---|---|---|---|---|---|---|
| shared_d0.25_rf1.0 | S | GSE216146 | up | +0.0020 | -0.0159 | -7.983 |
| shared_d0.25_rf1.0 | S | GSE216146 | down | -0.0290 | -0.0020 | +0.070 |
| shared_d0.25_rf1.0 | S | GSE271055 | up | +0.0074 | +0.0020 | +0.274 |
| shared_d0.25_rf1.0 | S | GSE271055 | down | -0.0279 | -0.0021 | +0.074 |
| shared_d0.25_rf0.3 | S | GSE216146 | up | -0.0080 | -0.0331 | +4.152 |
| shared_d0.25_rf0.3 | S | GSE216146 | down | +0.0107 | +0.0300 | +2.811 |
| shared_d0.25_rf0.3 | S | GSE271055 | up | +0.0063 | +0.0001 | +0.012 |
| shared_d0.25_rf0.3 | S | GSE271055 | down | -0.0169 | +0.0031 | -0.183 |
| shared_d0.5_rf1.0 | S | GSE216146 | up | +0.0201 | -0.0156 | -0.775 |
| shared_d0.5_rf1.0 | S | GSE216146 | down | -0.0489 | -0.0022 | +0.045 |
| shared_d0.5_rf1.0 | S | GSE271055 | up | +0.0199 | +0.0019 | +0.095 |
| shared_d0.5_rf1.0 | S | GSE271055 | down | -0.0479 | -0.0016 | +0.034 |
| shared_d0.5_rf0.3 | S | GSE216146 | up | +0.0082 | -0.0320 | -3.927 |
| shared_d0.5_rf0.3 | S | GSE216146 | down | -0.0170 | +0.0276 | -1.623 |
| shared_d0.5_rf0.3 | S | GSE271055 | up | +0.0193 | +0.0008 | +0.042 |
| shared_d0.5_rf0.3 | S | GSE271055 | down | -0.0430 | +0.0019 | -0.045 |
| shared_d1_rf1.0 | S | GSE216146 | up | +0.0532 | -0.0146 | -0.275 |
| shared_d1_rf1.0 | S | GSE216146 | down | -0.0807 | -0.0021 | +0.026 |
| shared_d1_rf1.0 | S | GSE271055 | up | +0.0402 | +0.0019 | +0.048 |
| shared_d1_rf1.0 | S | GSE271055 | down | -0.0772 | -0.0018 | +0.023 |
| shared_d1_rf0.3 | S | GSE216146 | up | +0.0405 | -0.0278 | -0.685 |
| shared_d1_rf0.3 | S | GSE216146 | down | -0.0600 | +0.0278 | -0.463 |
| shared_d1_rf0.3 | S | GSE271055 | up | +0.0342 | +0.0005 | +0.015 |
| shared_d1_rf0.3 | S | GSE271055 | down | -0.0584 | +0.0014 | -0.024 |
| shared_d2_rf1.0 | S | GSE216146 | up | +0.1032 | -0.0122 | -0.118 |
| shared_d2_rf1.0 | S | GSE216146 | down | -0.1433 | -0.0032 | +0.022 |
| shared_d2_rf1.0 | S | GSE271055 | up | +0.0697 | +0.0015 | +0.022 |
| shared_d2_rf1.0 | S | GSE271055 | down | -0.1288 | -0.0010 | +0.008 |
| shared_d2_rf0.3 | S | GSE216146 | up | +0.0873 | -0.0240 | -0.275 |
| shared_d2_rf0.3 | S | GSE216146 | down | -0.1193 | +0.0253 | -0.212 |
| shared_d2_rf0.3 | S | GSE271055 | up | +0.0619 | +0.0005 | +0.008 |
| shared_d2_rf0.3 | S | GSE271055 | down | -0.1251 | +0.0018 | -0.015 |
| specific_d0.25_rf1.0 | A | GSE216146 | up | +0.0237 | -0.0102 | -0.431 |
| specific_d0.25_rf1.0 | A | GSE216146 | down | -0.0171 | +0.0011 | -0.067 |
| specific_d0.25_rf1.0 | B | GSE271055 | up | +0.0186 | +0.0028 | +0.151 |
| specific_d0.25_rf1.0 | B | GSE271055 | down | -0.0295 | -0.0003 | +0.010 |
| specific_d0.25_rf0.3 | A | GSE216146 | up | +0.0316 | -0.0094 | -0.299 |
| specific_d0.25_rf0.3 | A | GSE216146 | down | +0.0066 | +0.0250 | +3.815 |
| specific_d0.25_rf0.3 | B | GSE271055 | up | +0.0096 | +0.0002 | +0.017 |
| specific_d0.25_rf0.3 | B | GSE271055 | down | -0.0231 | -0.0015 | +0.067 |
| specific_d0.5_rf1.0 | A | GSE216146 | up | +0.0490 | -0.0093 | -0.190 |
| specific_d0.5_rf1.0 | A | GSE216146 | down | -0.0383 | +0.0011 | -0.027 |
| specific_d0.5_rf1.0 | B | GSE271055 | up | +0.0382 | +0.0029 | +0.075 |
| specific_d0.5_rf1.0 | B | GSE271055 | down | -0.0495 | +0.0005 | -0.010 |
| specific_d0.5_rf0.3 | A | GSE216146 | up | +0.0538 | -0.0084 | -0.155 |
| specific_d0.5_rf0.3 | A | GSE216146 | down | -0.0083 | +0.0255 | -3.071 |
| specific_d0.5_rf0.3 | B | GSE271055 | up | +0.0300 | +0.0013 | +0.042 |
| specific_d0.5_rf0.3 | B | GSE271055 | down | -0.0391 | -0.0013 | +0.033 |
| specific_d1_rf1.0 | A | GSE216146 | up | +0.0871 | -0.0092 | -0.106 |
| specific_d1_rf1.0 | A | GSE216146 | down | -0.0650 | +0.0003 | -0.005 |
| specific_d1_rf1.0 | B | GSE271055 | up | +0.0684 | +0.0029 | +0.042 |
| specific_d1_rf1.0 | B | GSE271055 | down | -0.0813 | +0.0002 | -0.003 |
| specific_d1_rf0.3 | A | GSE216146 | up | +0.0926 | -0.0090 | -0.097 |
| specific_d1_rf0.3 | A | GSE216146 | down | -0.0490 | +0.0256 | -0.522 |
| specific_d1_rf0.3 | B | GSE271055 | up | +0.0605 | +0.0013 | +0.022 |
| specific_d1_rf0.3 | B | GSE271055 | down | -0.0756 | -0.0009 | +0.012 |
| specific_d2_rf1.0 | A | GSE216146 | up | +0.1458 | -0.0077 | -0.053 |
| specific_d2_rf1.0 | A | GSE216146 | down | -0.1054 | +0.0002 | -0.002 |
| specific_d2_rf1.0 | B | GSE271055 | up | +0.1079 | +0.0029 | +0.027 |
| specific_d2_rf1.0 | B | GSE271055 | down | -0.1230 | +0.0006 | -0.005 |
| specific_d2_rf0.3 | A | GSE216146 | up | +0.1562 | -0.0057 | -0.036 |
| specific_d2_rf0.3 | A | GSE216146 | down | -0.0891 | +0.0232 | -0.260 |
| specific_d2_rf0.3 | B | GSE271055 | up | +0.1004 | +0.0020 | +0.020 |
| specific_d2_rf0.3 | B | GSE271055 | down | -0.1135 | -0.0003 | +0.003 |
| mixed_d0.25_rf1.0 | S | GSE216146 | up | +0.0020 | -0.0162 | -8.278 |
| mixed_d0.25_rf1.0 | S | GSE216146 | down | -0.0290 | -0.0022 | +0.076 |
| mixed_d0.25_rf1.0 | S | GSE271055 | up | +0.0074 | +0.0022 | +0.300 |
| mixed_d0.25_rf1.0 | S | GSE271055 | down | -0.0279 | -0.0019 | +0.069 |
| mixed_d0.25_rf1.0 | A | GSE216146 | up | +0.0225 | -0.0103 | -0.455 |
| mixed_d0.25_rf1.0 | A | GSE216146 | down | -0.0195 | +0.0019 | -0.100 |
| mixed_d0.25_rf1.0 | B | GSE271055 | up | +0.0160 | +0.0027 | +0.168 |
| mixed_d0.25_rf1.0 | B | GSE271055 | down | -0.0286 | +0.0004 | -0.013 |
| mixed_d0.25_rf0.3 | S | GSE216146 | up | -0.0080 | -0.0323 | +4.032 |
| mixed_d0.25_rf0.3 | S | GSE216146 | down | +0.0106 | +0.0284 | +2.679 |
| mixed_d0.25_rf0.3 | S | GSE271055 | up | +0.0063 | +0.0003 | +0.042 |
| mixed_d0.25_rf0.3 | S | GSE271055 | down | -0.0169 | +0.0022 | -0.128 |
| mixed_d0.25_rf0.3 | A | GSE216146 | up | +0.0317 | -0.0093 | -0.292 |
| mixed_d0.25_rf0.3 | A | GSE216146 | down | +0.0080 | +0.0257 | +3.209 |
| mixed_d0.25_rf0.3 | B | GSE271055 | up | +0.0102 | +0.0005 | +0.052 |
| mixed_d0.25_rf0.3 | B | GSE271055 | down | -0.0271 | -0.0018 | +0.068 |
| mixed_d0.5_rf1.0 | S | GSE216146 | up | +0.0201 | -0.0152 | -0.759 |
| mixed_d0.5_rf1.0 | S | GSE216146 | down | -0.0490 | -0.0030 | +0.061 |
| mixed_d0.5_rf1.0 | S | GSE271055 | up | +0.0199 | +0.0021 | +0.105 |
| mixed_d0.5_rf1.0 | S | GSE271055 | down | -0.0479 | -0.0024 | +0.049 |
| mixed_d0.5_rf1.0 | A | GSE216146 | up | +0.0462 | -0.0096 | -0.208 |
| mixed_d0.5_rf1.0 | A | GSE216146 | down | -0.0349 | +0.0012 | -0.035 |
| mixed_d0.5_rf1.0 | B | GSE271055 | up | +0.0376 | +0.0033 | +0.087 |
| mixed_d0.5_rf1.0 | B | GSE271055 | down | -0.0508 | -0.0001 | +0.003 |
| mixed_d0.5_rf0.3 | S | GSE216146 | up | +0.0081 | -0.0326 | -4.025 |
| mixed_d0.5_rf0.3 | S | GSE216146 | down | -0.0171 | +0.0280 | -1.633 |
| mixed_d0.5_rf0.3 | S | GSE271055 | up | +0.0193 | +0.0006 | +0.033 |
| mixed_d0.5_rf0.3 | S | GSE271055 | down | -0.0430 | +0.0030 | -0.070 |
| mixed_d0.5_rf0.3 | A | GSE216146 | up | +0.0532 | -0.0101 | -0.189 |
| mixed_d0.5_rf0.3 | A | GSE216146 | down | -0.0136 | +0.0247 | -1.815 |
| mixed_d0.5_rf0.3 | B | GSE271055 | up | +0.0286 | +0.0011 | +0.038 |
| mixed_d0.5_rf0.3 | B | GSE271055 | down | -0.0411 | -0.0007 | +0.016 |
| mixed_d1_rf1.0 | S | GSE216146 | up | +0.0531 | -0.0139 | -0.261 |
| mixed_d1_rf1.0 | S | GSE216146 | down | -0.0809 | -0.0026 | +0.032 |
| mixed_d1_rf1.0 | S | GSE271055 | up | +0.0402 | +0.0020 | +0.050 |
| mixed_d1_rf1.0 | S | GSE271055 | down | -0.0772 | -0.0016 | +0.021 |
| mixed_d1_rf1.0 | A | GSE216146 | up | +0.0816 | -0.0093 | -0.114 |
| mixed_d1_rf1.0 | A | GSE216146 | down | -0.0662 | +0.0011 | -0.016 |
| mixed_d1_rf1.0 | B | GSE271055 | up | +0.0683 | +0.0027 | +0.040 |
| mixed_d1_rf1.0 | B | GSE271055 | down | -0.0845 | -0.0000 | +0.000 |
| mixed_d1_rf0.3 | S | GSE216146 | up | +0.0404 | -0.0282 | -0.697 |
| mixed_d1_rf0.3 | S | GSE216146 | down | -0.0601 | +0.0279 | -0.464 |
| mixed_d1_rf0.3 | S | GSE271055 | up | +0.0343 | +0.0008 | +0.024 |
| mixed_d1_rf0.3 | S | GSE271055 | down | -0.0584 | +0.0025 | -0.042 |
| mixed_d1_rf0.3 | A | GSE216146 | up | +0.0942 | -0.0089 | -0.095 |
| mixed_d1_rf0.3 | A | GSE216146 | down | -0.0342 | +0.0237 | -0.694 |
| mixed_d1_rf0.3 | B | GSE271055 | up | +0.0558 | +0.0015 | +0.027 |
| mixed_d1_rf0.3 | B | GSE271055 | down | -0.0798 | -0.0011 | +0.013 |
| mixed_d2_rf1.0 | S | GSE216146 | up | +0.1031 | -0.0109 | -0.106 |
| mixed_d2_rf1.0 | S | GSE216146 | down | -0.1436 | -0.0042 | +0.029 |
| mixed_d2_rf1.0 | S | GSE271055 | up | +0.0697 | +0.0022 | +0.031 |
| mixed_d2_rf1.0 | S | GSE271055 | down | -0.1288 | -0.0011 | +0.009 |
| mixed_d2_rf1.0 | A | GSE216146 | up | +0.1467 | -0.0063 | -0.043 |
| mixed_d2_rf1.0 | A | GSE216146 | down | -0.1080 | -0.0008 | +0.007 |
| mixed_d2_rf1.0 | B | GSE271055 | up | +0.1097 | +0.0031 | +0.029 |
| mixed_d2_rf1.0 | B | GSE271055 | down | -0.1304 | +0.0003 | -0.002 |
| mixed_d2_rf0.3 | S | GSE216146 | up | +0.0872 | -0.0234 | -0.269 |
| mixed_d2_rf0.3 | S | GSE216146 | down | -0.1196 | +0.0254 | -0.212 |
| mixed_d2_rf0.3 | S | GSE271055 | up | +0.0619 | +0.0009 | +0.014 |
| mixed_d2_rf0.3 | S | GSE271055 | down | -0.1252 | +0.0015 | -0.012 |
| mixed_d2_rf0.3 | A | GSE216146 | up | +0.1577 | -0.0063 | -0.040 |
| mixed_d2_rf0.3 | A | GSE216146 | down | -0.0829 | +0.0226 | -0.273 |
| mixed_d2_rf0.3 | B | GSE271055 | up | +0.1013 | +0.0013 | +0.013 |
| mixed_d2_rf0.3 | B | GSE271055 | down | -0.1244 | -0.0017 | +0.014 |

## Reconstructed effect, null configuration (nothing injected; nominal responder sets)

| responder set | set | study | direction | observed | reconstructed |
|---|---|---|---|---|---|
| 1 | S | GSE216146 | up | -0.0188 | -0.0168 |
| 1 | S | GSE216146 | down | -0.0069 | -0.0024 |
| 1 | S | GSE271055 | up | -0.0070 | +0.0022 |
| 1 | S | GSE271055 | down | -0.0065 | -0.0019 |
| 1 | A | GSE216146 | up | -0.0014 | -0.0113 |
| 1 | A | GSE216146 | down | -0.0015 | +0.0020 |
| 1 | B | GSE271055 | up | -0.0040 | +0.0027 |
| 1 | B | GSE271055 | down | -0.0095 | -0.0002 |
| 0.3 | S | GSE216146 | up | -0.0279 | -0.0343 |
| 0.3 | S | GSE216146 | down | +0.0293 | +0.0303 |
| 0.3 | S | GSE271055 | up | -0.0086 | +0.0002 |
| 0.3 | S | GSE271055 | down | +0.0026 | +0.0029 |
| 0.3 | A | GSE216146 | up | +0.0078 | -0.0109 |
| 0.3 | A | GSE216146 | down | +0.0271 | +0.0273 |
| 0.3 | B | GSE271055 | up | -0.0104 | +0.0005 |
| 0.3 | B | GSE271055 | down | -0.0016 | -0.0010 |

Over all non-injected set × study × direction rows in all configurations (n = 80): observed effect range -0.0281 to +0.0293; reconstructed range -0.0346 to +0.0303.
Over injected rows (n = 128): observed range -0.1436 to +0.1577; reconstructed range -0.0331 to +0.0300.

Ratios are unstable when the observed effect is close to 0 (e.g. δ = 0.25 and non-injected rows).
