# Task 2: replicate-block library-label exchangeability reference for the shared fraction

## Audit correction

A post-analysis audit of the GEO series metadata showed that the initial 24-assignment reference grouped all four D21 non-rescue libraries together. GEO instead identifies three experimental replicates: replicate 1 contains D20-6407 (PBS nonstim) and D20-6409 (cisplatin nonstim), replicate 2 contains D21-2746 and D21-2750, and replicate 3 contains D21-2747 and D21-2751. The corrected primary reference therefore exchanges control/cisplatin labels only within each replicate pair. This gives 2^3 = 8 cisplatin-study assignments and, after the two-way GSE271055 pooled-library exchange, 16 combined assignments.

No model refitting was needed for this correction: all 16 replicate-preserving assignments were already contained in the original 24-assignment superset, and every assignment/seed fit reset the PyTorch and NumPy seeds independently. The corrected analysis therefore filters the existing 72-fit audit to the 48 fits belonging to the replicate-preserving reference set.

## Results

For the observed label assignment, mean F_sh across seeds was 0.776222. Twelve of the 16 replicate-preserving assignments had a statistic at least this large, giving an exact upper-tail reference probability of 12/16 = 0.750. The corrected reference distribution had median 0.887638 and range 0.292013-0.988365; the observed assignment ranked 12th from the top.

With doxorubicin labels fixed and only the three cisplatin-study replicate pairs exchanged, 5 of 8 assignments were at least as large as observed (P = 0.625; reference median 0.832065). With the cisplatin labels fixed and only the two pooled doxorubicin-study labels exchanged, both assignments were at least as large as observed (P = 1.0); this comparison has a minimum attainable one-sided probability of 0.5.

Seed-specific observed shared fractions were 0.621400, 0.969542 and 0.737724. Their corrected combined-reference probabilities were 0.8125, 0.125 and 0.750. The result therefore remains seed-dependent and does not show an unusually large observed-label shared fraction.

## Interpretation

The correction changes the exact reference set and numerical probabilities but not the qualitative conclusion. Under the replicate-preserving, model-specific sensitivity reference, large shared fractions remain common under admissible alternative library labels. These probabilities are conditional on the stated exchangeability assumptions: the source experiments were not documented as randomized according to these relabelings, and the doxorubicin study contributes one pooled library per arm. They are not randomized-treatment-effect P values and are not evidence that chemotherapy has no biological effect.
