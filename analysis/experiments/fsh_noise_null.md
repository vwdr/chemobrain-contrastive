# F_sh on pure noise

Simulation: 10000 draws, 358 doxorubicin + 357 cisplatin treated test cells, collapsed blocks = iid Gaussian posterior means (variance 0.0001 in every block), drug blocks gated as in the model (a cell's non-matching drug coordinates are zero), F_sh from `src.peerreview.core.shared_fractions`. The cell counts are an assumption; the expected value does not depend on them, the spread does.

## Observed null fits versus matched-variance noise

**All three observed null values (0.705, 0.802, 0.869) exceed the matched-variance simulation.** With 8 shared + two 4-dim drug blocks at equal variance the simulated pooled F_sh has mean 0.6668 (analytic 8/(8+4) = 0.6667), 97.5th percentile 0.6809 and maximum 0.6932 over 10,000 draws, so none of the draws reaches even the smallest observed value.

Shared:drug per-dimension variance ratio r that would give each observed value (analytic, E[F] = 8r / (8r + 4) so r = 4F / (8(1 - F))); r = 1 is the matched-variance case simulated here:

| observed null F_sh | P(noise >= observed) | required r |
|---|---|---|
| 0.705 | 0.0000 | 1.19 |
| 0.802 | 0.0000 | 2.03 |
| 0.869 | 0.0000 | 3.32 |

The ratio says how much more variance per dimension the shared block would need than a drug dimension for noise alone to produce the observed value; it does not say the shared block has it.

## Dimension split alone: total latent budget 16, shared dimensions 2 to 14

Shared dimensions s and two drug blocks of (16 - s) / 2 dimensions each, equal per-dimension variance, pooled F_sh. Simulated mean with its analytic expectation s / (s + (16 - s) / 2) alongside:

| shared dims | drug dims per block | shared:drug dims | simulated mean | sd | 2.5% | 97.5% | analytic |
|---|---|---|---|---|---|---|---|
| 2 | 7 | 0.29 | 0.2223 | 0.0073 | 0.2080 | 0.2367 | 0.2222 |
| 4 | 6 | 0.67 | 0.4000 | 0.0082 | 0.3834 | 0.4161 | 0.4000 |
| 6 | 5 | 1.20 | 0.5455 | 0.0079 | 0.5299 | 0.5610 | 0.5455 |
| 8 | 4 | 2.00 | 0.6667 | 0.0072 | 0.6525 | 0.6807 | 0.6667 |
| 10 | 3 | 3.33 | 0.7692 | 0.0062 | 0.7568 | 0.7809 | 0.7692 |
| 12 | 2 | 6.00 | 0.8571 | 0.0050 | 0.8473 | 0.8667 | 0.8571 |
| 14 | 1 | 14.00 | 0.9334 | 0.0034 | 0.9265 | 0.9399 | 0.9333 |

The model's own split (8 shared, 4 + 4 drug) sits at 0.667. At equal variances the observed null values would need 8.7, 10.7, 12.3 shared dimensions out of 16 (F = 2s / (s + 16)), against the 8 the model has.

## Main case, all three F_sh definitions

| F_sh | mean | sd | 2.5% | median | 97.5% | analytic E |
|---|---|---|---|---|---|---|
| pooled | 0.6668 | 0.0072 | 0.6523 | 0.6669 | 0.6809 | 0.6667 |
| doxorubicin | 0.6667 | 0.0102 | 0.6468 | 0.6667 | 0.6866 | 0.6667 |
| cisplatin | 0.6668 | 0.0102 | 0.6464 | 0.6669 | 0.6866 | 0.6667 |

## Drug block size with shared fixed at 8 (pooled F_sh, equal variances)

| shared dims | drug dims per block | shared:drug | mean | 2.5% | 97.5% | analytic |
|---|---|---|---|---|---|---|
| 8 | 1 | 8.00 | 0.8890 | 0.8776 | 0.8998 | 0.8889 |
| 8 | 2 | 4.00 | 0.7999 | 0.7866 | 0.8129 | 0.8000 |
| 8 | 4 | 2.00 | 0.6665 | 0.6521 | 0.6806 | 0.6667 |
| 8 | 8 | 1.00 | 0.5001 | 0.4873 | 0.5129 | 0.5000 |
| 8 | 16 | 0.50 | 0.3334 | 0.3233 | 0.3434 | 0.3333 |
