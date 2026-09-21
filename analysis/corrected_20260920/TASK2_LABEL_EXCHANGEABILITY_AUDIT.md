# Task 2: library-label exchangeability reference for the shared fraction

## Design

A separate label-independent sensitivity benchmark was constructed to ask whether the observed treatment labels produce an unusually large shared posterior-mean variance fraction. This analysis is deliberately separate from the canonical cell-level benchmark.

Eight non-rescue deposited libraries were used. Cell selection was fixed before treatment labels were permuted. Each GSE216146 library contributed 500 training, 50 validation and 50 test cells; each GSE271055 pooled library contributed 1,500 training, 150 validation and 150 test cells. This yields 6,000 training cells, 600 validation cells and 600 test cells, with 3,000 training cells per study and 1,500 training cells in each study-by-treatment stratum under every admissible label assignment. The 1,500 highly variable genes were selected once from the fixed training cells using study as the batch key and then held fixed.

The exact exchangeability set preserved the observed treatment counts within the GSE216146 D20 and D21 preparation batches and within the two GSE271055 pooled libraries. This gives 2 x 6 x 2 = 24 assignments. Seeds 0, 1 and 2 were refit for every assignment, for 72 fits total. The primary statistic was the mean pooled shared fraction across the three fixed seeds.

## Results

For the observed label assignment, mean F_sh across seeds was 0.776222. Sixteen of the 24 admissible assignments had a statistic at least this large, giving an exact upper-tail reference probability of 16/24 = 0.667. The 24-assignment reference distribution had median 0.866731 and range 0.232501-0.988365. The observed assignment ranked 16th from the top and was not in the upper tail.

With the doxorubicin labels fixed and only the cisplatin-study labels permuted within D20/D21, 8 of 12 assignments were at least as large as observed (P = 0.667; reference median 0.817119). With the cisplatin labels fixed and the two pooled doxorubicin-study labels exchanged, both assignments were at least as large as observed (P = 1.0). The latter comparison has only two possible assignments and a minimum attainable one-sided probability of 0.5.

Seed-specific observed shared fractions were 0.621400, 0.969542 and 0.737724. Their combined-24 upper-tail probabilities were 0.833, 0.167 and 0.625. Seed 1 reached the minimum attainable cisplatin-only probability of 1/12 = 0.0833, but this pattern was not reproduced by seeds 0 or 2.

The observed assignment's held-out MSE was 0.140131-0.140470 across seeds. This number is specific to the balanced Task 2 benchmark and should not be compared numerically with the canonical benchmark's MSE because the evaluation cells and gene universe differ.

## Interpretation

Under this model-specific library-label exchangeability reference, a large shared fraction is common for admissible alternative labels. The observed labels do not yield an unusually large shared fraction. This reinforces the Task 1 finding that F_sh is a sensitive internal allocation statistic rather than direct evidence of a stable cross-drug biological program.

The exact probabilities are conditional on the stated exchangeability assumptions. The source experiments were not randomized according to this relabeling scheme, and the doxorubicin study contributes only one pooled library per arm. These values are therefore sensitivity-reference probabilities, not randomized-trial treatment-effect P values and not evidence that the treatments have no biological effect.
