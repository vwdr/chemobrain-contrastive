# Task 7: robust gene-pattern and marker-set audit

## Motivation and primary list definition

The original biological interpretation used top-50 integrated-gradient lists from one baseline. Because the corrected attribution rankings are sensitive to initialization and baseline choice, Task 7 replaces that narrative with a stricter seed-and-baseline consensus.

The attribution target is the squared norm of the posterior-mean latent block, which is invariant to sign flips and orthogonal rotations within a block. For each block, a gene enters the primary consensus only if it appears in the top 50 mean-absolute integrated-gradient ranking in at least four of six seed-by-baseline analyses (seeds 0, 1 and 2; zero and control-median baselines). The resulting frozen consensus sizes are:

- shared-treatment block: 34 genes;
- doxorubicin-associated block: 20 genes;
- cisplatin-associated block: 33 genes.

All enrichment tests use the exact 1,500-gene canonical model universe.

## Robust gene-set interpretation

The consensus lists were tested against mouse MSigDB 2025.1 Reactome, GO Biological Process and M8 Cell Type Signature collections, restricting tested sets to 10-500 model-universe genes and applying BH correction separately within each latent block and resource family.

The shared consensus retains immune/myeloid expression-pattern associations: 46 M8 signatures and 11 GO Biological Process terms pass FDR 0.05. The strongest M8 association is a macrophage signature (14 of 34 consensus genes; adjusted P approximately 4.18e-10). Reactome immune system does not remain below FDR under the stricter consensus (adjusted P approximately 0.092), demonstrating resource/list-definition dependence.

The cisplatin consensus retains stromal/pericyte and extracellular-matrix structure. Seventy-seven M8 signatures and five Reactome terms pass FDR 0.05. Reactome extracellular-matrix degradation contains five of 33 consensus genes, has approximately 7.58-fold enrichment and adjusted P approximately 0.034.

The doxorubicin consensus has no Reactome, GO or M8 term below FDR 0.05. Therefore the postsynaptic enrichment visible in the historical zero-baseline top-50 analysis is not robust to requiring agreement across seeds and baselines and should not be presented as a stable doxorubicin-associated pathway result.

## Cisplatin pseudobulk cross-check

Task 4 paired pseudobulk results are regenerated from the recovered-count checksum. Within the 1,500-gene model universe, 68 genes have PyDESeq2 adjusted P<0.05 for cisplatin versus control in at least one eligible source cell type and all three paired log2(CPM+1) effects agreeing with the fitted effect direction.

Five of the 33 cisplatin consensus genes overlap this robust support set: Apoe, Mbp, Mt1, Tmsb10 and Ttr. The hypergeometric probability is approximately 0.0146; after BH adjustment across the three latent-axis overlap tests, q is approximately 0.0219. This is limited within-study support for part of the cisplatin ranking, not validation of a complete latent program.

The shared consensus overlaps seven cisplatin-supported genes, but that one-drug overlap cannot establish a cross-drug shared mechanism. The doxorubicin study has one pooled deposited library per arm, so there is no matched replicated doxorubicin gene-wise validation.

## Cross-method recurrence

The Task 6 scDisInFact top-100 condition-associated gene ranking overlaps 6/34 shared, 6/20 doxorubicin and 9/33 cisplatin consensus genes. The corresponding BH-adjusted hypergeometric probabilities are approximately 0.0220, 0.00204 and 0.000549. This indicates ranking recurrence across architectures, but not independent biological validation because scDisInFact also retained strong study information under the same drug/study aliasing.

## Interpretation

The robust result is narrower than the historical biological narrative. The shared block is associated with immune/myeloid expression patterns and the cisplatin block with stromal/pericyte/ECM patterns; neither establishes cell of origin, activation, causality or cross-drug convergence. M8 labels aggregate signatures from heterogeneous tissues and experiments and are used only as expression-pattern annotations. The doxorubicin synaptic theme is not robust under the primary consensus.

The zero-baseline top-50 enrichment remains archived as a sensitivity analysis rather than the primary biological interpretation.
