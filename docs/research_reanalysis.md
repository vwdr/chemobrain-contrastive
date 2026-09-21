# Corrected research analysis

## Scope and provenance

The upstream source commit is `e6f2816220bb47d2b1c792603e46b65b64de9e40`. The completed analysis uses three seeds for each of the full model, no-HSIC ablation, no-gating ablation, noncontrastive Gaussian baseline and negative-binomial sensitivity model, plus PCA. Source and output paths have been reorganized for a research-only repository. No numerical findings were changed during publication of this update.

GSE216146 is hippocampal whole-cell expression with three deposited libraries per arm. Its specific series design and generic extraction fields disagree about pooling in the D20 preparations. Deposited libraries are the available independent units. GSE271055 is hippocampal nucleus expression with one library pooling four mice in each arm. GSE286221 is bulk microglial expression and is excluded. No independent external atlas was used for annotation.

## Preprocessing and benchmark

Within each study, QC requires at least 500 detected genes per observation, detection in at least ten observations per gene and mitochondrial fraction below 10%. No new doublet detection was performed. Intersection retains 18,271 genes. The cohort contains 28,808 nonrescue and 18,049 rescue observations after filtering.

Split seed 1729 produces stratified 80/10/10 candidate training, validation and test sets. Training is subsampled to 1,500 cells per study-by-condition stratum. Seurat-v3 HVGs use training counts only. Normalization uses 10,000 counts across common genes followed by log1p. Gene detection QC is transductive and libraries overlap between splits. Evaluation measures within-library interpolation, not animal-held-out generalization.

Models use background dimension 16, shared dimension eight and drug dimensions four each, two hidden layers of width 128, dropout 0.1, AdamW learning rate 0.001, weight decay 0.00001 and gradient clipping at five. Seeds are 0, 1 and 2. Training is capped at 100 epochs and checkpoint selection uses validation posterior-mean reconstruction MSE with patience 15 and improvement tolerance 0.00001.

The inherited loss sums conditional KL averages without prevalence weighting. HSIC weights are 10, 10 and five. The shared HSIC input is ungated. Removing gating also changes the KL restrictions and is not a pure one-factor ablation. The noncontrastive baseline removes gating and HSIC. The count model uses softmax decoder means scaled by observed modeled-gene counts and positive gene-specific dispersion. Cross-likelihood log-expression MSE does not assess count-distribution fidelity.

## Interpretation and uncertainty

Shared variance fractions are calculated on treated test-cell posterior means. Cell bootstrap intervals describe sampling sensitivity conditional on a fit, not biological uncertainty. Drug, study and sequencing modality are confounded among treated observations.

Integrated gradients use fixed held-out cells across seeds, coordinate-sum and squared-norm targets, and zero and control-median baselines. Adaptive Gauss-Legendre integration satisfied the median scaled completeness threshold below 0.01 for all 36 combinations. Earlier inadequate quadrature diagnostics are retained separately. Their rankings are not used. Attribution signs do not estimate expression direction.

Exploratory enrichment uses mouse MSigDB 2025.1 Reactome and GO Biological Process. Hypergeometric tests assess top-50 rankings against the modeled 1,500-gene universe, with Benjamini-Hochberg correction within each latent-block and resource family. Significant annotations are not evidence of causal pathway activation.

Cisplatin source labels are preserved. Doxorubicin transfer requires classifier probability at least 0.8 and marker agreement. Many nuclei remain uncertain. Reference feature selection is transductive in the annotation leave-library-out diagnostic, and probabilities are not independently calibrated.

Exploratory pseudobulk sums counts within source cell type and library, requires three libraries in both arms and uses Welch tests on log-CPM with within-comparison FDR control. It lacks moderated count-based dispersion estimates. No inferential doxorubicin differential-expression test is performed.

Rescue projections use ungated shared means and study-specific training-control centroids, summarized per library and averaged across seeds. The cisplatin test enumerates 12 preparation-batch-stratified allocations and yields a two-sided probability of 1.0. Doxorubicin is descriptive because its arms are unreplicated. These results do not establish or refute behavioral rescue.

## Remaining validation

Independent animal-held-out evaluation, tuned scVI and pairwise contrastiveVI comparisons, broader hyperparameter sensitivity, independent nucleus-reference annotation and replicated within-study cross-drug experiments remain outstanding. The completed results support exploratory representations, not an established common toxicity mechanism.
