# Task 4: paired pseudobulk differential-expression audit

## Primary design

The GSE216146 pseudobulk analysis aggregates integer counts within source cell type and deposited library. The three GEO-defined experimental replicate blocks are retained and each eligible 3-vs-3 comparison is fit with PyDESeq2 0.5.4 using

`~ replicate + condition`.

A cell type/contrast is analyzed only when all six required library-by-cell-type aggregates contain at least 20 cells. Genes require at least 10 raw pseudobulk counts in at least 3 of the 6 libraries. Benjamini-Hochberg adjusted Wald probabilities from the paired negative-binomial model are the primary model-based inferential summary.

Two sensitivity analyses are retained:

1. paired log2(CPM+1) effects with all 2^3 within-replicate sign flips;
2. the historical unpaired Welch test on the same six log2(CPM+1) library values.

GSE271055 has one pooled deposited library per arm and is descriptive only.

## Reproducibility check

Two independently implemented Task 4 scripts were run against the same recovered-data checksum. Across all 131,657 tested gene-by-cell-type/contrast rows, the PyDESeq2 base means, log2 fold changes, standard errors, Wald P values and adjusted P values were identical. The paired implementation is retained as canonical because it also exports per-replicate effects and exact sign-flip sensitivity results.

## Results

Sixteen cell-type/contrast combinations met the six-library eligibility rule. PyDESeq2 reported 354 gene-by-cell-type/contrast discoveries at adjusted P < 0.05:

- cisplatin vs control: 297;
- cisplatin + GENUS vs cisplatin: 55;
- control + GENUS vs control: 2.

The largest discovery counts occurred in oligodendrocytes (137 cisplatin-vs-control; 36 rescue-vs-cisplatin) and microglia (109 cisplatin-vs-control; 11 rescue-vs-cisplatin).

The PyDESeq2 log2 fold changes agreed closely with paired mean log2(CPM+1) effects across all eligible comparisons (Spearman rho 0.932-0.997). Of the 354 adjusted-P<0.05 PyDESeq2 discoveries, 353 had the same effect direction in all three GEO replicate pairs; the remaining discovery had the model direction in two pairs and a zero paired effect in the third.

Inference is nonetheless sensitive to the modeling assumptions. With only three replicate pairs, the exact paired two-sided sign-flip reference has only 2^3=8 assignments and a minimum attainable P value of 0.25. Consequently, no gene can reach P < 0.05 under that exact sensitivity reference. The historical unpaired Welch sensitivity analysis yielded only one FDR<0.05 gene across the 16 eligible comparisons. Thus the large difference in discovery counts reflects the paired negative-binomial model and cross-gene dispersion moderation, not disagreement in effect direction.

## Interpretation

The paired pseudobulk analysis supports reproducible within-GSE216146 expression shifts in several source cell types, particularly for cisplatin vs control. However, with three deposited replicate blocks, the adjusted Wald discoveries should be described as model-based exploratory differential-expression evidence rather than definitive population-level inference.

These results do not validate the latent gene-attribution rankings and do not address cross-drug convergence. The doxorubicin study remains unsuitable for gene-wise biological inference because each arm is represented by one pooled deposited library.
