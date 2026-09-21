# Task 8: related-work update

## Scope

The manuscript's related-work section was updated against primary literature available through September 2026. The goal was not to create a broad review of single-cell perturbation modeling, but to position the present analysis among methods that solve distinct tasks: integration, shared/group-specific decomposition, batch/condition disentanglement and counterfactual perturbation prediction.

## Added or corrected references

- scGen: Lotfollahi et al., Nature Methods 16, 715-721 (2019), DOI 10.1038/s41592-019-0494-8.
- trVAE: Lotfollahi et al., Bioinformatics 36, i610-i617 (2020), DOI 10.1093/bioinformatics/btaa800.
- multiGroupVI: Weinberger et al., Proceedings of Machine Learning Research 200, 16-32 (2022). The manuscript's previous bioRxiv-style entry was corrected to the peer-reviewed MLCB proceedings record.
- CPA: Lotfollahi et al., Molecular Systems Biology 19, e11517 (2023), DOI 10.15252/msb.202211517.
- biolord: Piran et al., Nature Biotechnology 42, 1678-1683 (2024), DOI 10.1038/s41587-023-02079-x.
- simple-baseline benchmark: Ahlmann-Eltze et al., Nature Methods 22, 1657-1661 (2025), DOI 10.1038/s41592-025-02772-6.
- broad perturbation benchmark: Wei et al., Nature Methods 23, 451-464 (2026), DOI 10.1038/s41592-025-02980-0.

Existing citations to scVI, contrastiveVI and scDisInFact were retained.

## Positioning changes

The revised introduction explicitly separates three related but non-equivalent goals:

1. latent-variable integration and representation learning;
2. decomposition or disentanglement of shared, group-specific, condition and batch-associated factors;
3. counterfactual perturbation-response prediction.

The manuscript now states that successful reconstruction or perturbation prediction does not by itself identify a latent component as a biological mechanism shared across drugs.

Recent benchmark literature is used only for the methodological lesson that evaluation must match the claimed generalization target and that simple baselines remain necessary. Numerical results from those external benchmarks are not transferred to the present chemotherapy datasets.

## Comparator rationale

scDisInFact remains the only external method fitted in this manuscript because it is the closest structural comparator for simultaneous batch/condition disentanglement. scGen, trVAE and CPA are discussed as neighboring perturbation-prediction methods but were not fit. The present two-study dataset does not cross both drugs with both studies, so a formal unseen-drug prediction benchmark would not resolve the core study/drug aliasing.

## Interpretation

The literature update does not change the empirical results from Tasks 1-7. It sharpens the claim boundary: this work is an audit of representation reproducibility and identifiability in a confounded retrospective design, not a benchmark of general perturbation-response prediction.
