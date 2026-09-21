# Task 3: rescue projection permutation-resolution audit

## Endpoint

The rescue endpoint is the Euclidean distance of the ungated shared posterior mean from the study-matched training-control centroid. Distances are summarized within each deposited library and then averaged across model seeds. For cisplatin, the contrast is rescue minus untreated cisplatin; a negative value would indicate movement toward the control centroid.

## Replicate structure

The GEO series record identifies three experimental replicate blocks. The cisplatin nonstim / cisplatin + GENUS pairs are:

- replicate 1: D20-6409 / D20-6410;
- replicate 2: D21-2750 / D21-2752;
- replicate 3: D21-2751 / D21-2753.

The previous D20/D21 permutation scheme incorrectly treated the four D21 libraries as one exchangeability stratum. The corrected sensitivity reference exchanges treatment/rescue labels only within each GEO replicate pair. This is equivalent to an exact paired sign-flip reference with 2^3 = 8 assignments.

## Result

Using the committed per-library rescue projections, the seed-averaged paired differences were:

- replicate 1: -0.019507;
- replicate 2: +0.016234;
- replicate 3: +0.004245.

The mean rescue-minus-treatment difference was +0.000324, indicating essentially no net movement toward the control centroid. All 8 sign assignments had an absolute mean difference at least as large as the observed value, so the exact two-sided reference probability is P = 1.000.

Because the 8-assignment sign-flip distribution is symmetric, the smallest attainable two-sided probability is 2/8 = 0.25. Even a one-sided exact reference has a minimum attainable probability of 1/8 = 0.125. This design therefore cannot support a conventional alpha = 0.05 significance claim regardless of the observed effect.

Seed-specific mean differences were +0.001021, -0.005827 and +0.005778 for seeds 0, 1 and 2, with exact two-sided reference probabilities 1.00, 0.50 and 0.75. These are computational sensitivity diagnostics, not additional biological replicates.

The doxorubicin study contributes one pooled doxorubicin library and one pooled doxorubicin + ACY-1083 library, so no biological inferential test is performed for that rescue comparison.

## Interpretation

The corrected replicate-pair analysis leaves the substantive conclusion unchanged: this latent-space endpoint does not establish rescue or reversal of a common treatment-associated state. The exact probability is conditional on pairwise exchangeability and should be treated as a low-resolution sensitivity reference, not a randomized-treatment-effect test. Failure of this endpoint does not refute rescue effects measured using behavioral, histologic or molecular endpoints in the source experiments.
