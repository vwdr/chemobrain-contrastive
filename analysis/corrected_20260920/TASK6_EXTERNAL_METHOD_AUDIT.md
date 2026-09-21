# Task 6: external-method comparison with scDisInFact

## Comparator and design

scDisInFact was selected because it is an independently developed conditional variational model for multi-batch, multi-condition single-cell data with explicit shared-bio (condition-irrelevant) and unshared-bio (condition-associated) factors. The official upstream repository was pinned at commit `9466ecf257d923f879cb39712c530d24826d88bc` (Zhang et al., Nature Communications 15, 912, 2024).

The comparison used the exact corrected canonical benchmark: 6,000 training cells, 2,881 held-out test cells and the committed 1,500-gene universe. The treatment label was supplied as the single condition factor and study as the batch factor. The upstream one-condition robustness configuration was used without dataset-specific hyperparameter tuning: `Ks=[8,2]`, batch size 64, learning rate 5e-4, 50 epochs, default regularization weights, and seeds 0-2.

The original condition/study aliasing remains present. Cisplatin occurs only in GSE216146 and doxorubicin only in GSE271055, while control occurs in both. No method comparison can identify separate drug and study effects from that design alone.

## Reconstruction

Mean held-out log-expression MSE was 0.23134 for scDisInFact (SD 0.00217 across seeds). On the same canonical cells and 1,500-gene target, the archived MC-ContrastiveVI checkpoints had mean MSE 0.16871 and PCA(32) had MSE 0.15525. Thus scDisInFact's MSE was approximately 37% higher than MC-ContrastiveVI and 49% higher than PCA on this particular common output-scale metric.

This is a descriptive reconstruction comparison, not a likelihood ranking. scDisInFact was trained with its native negative-binomial objective and upstream defaults rather than tuned for log-expression MSE.

## Representation probes

The scDisInFact shared-bio factor, which is intended to be condition-irrelevant, retained balanced condition-prediction accuracy 0.560 on held-out cells and strong study predictability (0.895; chance for study is 0.5). The unshared-bio factor had condition balanced accuracy 0.756 and study balanced accuracy 0.981. The combined representation had condition balanced accuracy 0.761 and study balanced accuracy 0.986.

These results show that the external comparator also retained substantial study information under the present confounded design. High study predictability in the unshared factor is especially unsurprising because drug and study are partially aliased.

The latent semantics are not interchangeable with MC-ContrastiveVI. scDisInFact's shared-bio factor is intended to remove condition effects, whereas MC-ContrastiveVI's shared-treatment block is intended to contain treatment-associated signal common to drugs. The comparison therefore does not validate or refute the MC shared-fraction statistic directly.

## Cross-seed stability

scDisInFact representations were substantially more stable across the three computational seeds than the archived MC salient spaces. Pairwise linear CKA for scDisInFact was 0.794-0.854 for shared-bio (mean 0.815) and 0.781-0.881 for unshared-bio (mean 0.842). For reference, the archived MC-ContrastiveVI pairwise CKA values were 0.082-0.525 for the shared block (mean 0.322) and 0.115-0.389 for the drug block (mean 0.240).

Condition-associated scDisInFact gene scores were also stable across seeds: pairwise Spearman correlations were 0.843-0.882, with 67-79 genes shared between each pair of top-100 lists (mean 72.3; Jaccard 0.504-0.653).

This stability is a property of the fitted external method, not evidence that its condition-associated genes are chemotherapy-specific. Because the condition variable is confounded with study for the treated groups, those scores can contain technical and study-associated information.

## Interpretation

The external comparison provides two useful controls. First, the failure to remove study information is not unique to MC-ContrastiveVI; scDisInFact also retains strong study predictability when applied to this partially aliased two-study design. Second, stable latent factors and gene rankings are achievable with an independently developed disentanglement framework on the same cells, so the low cross-seed stability of MC-ContrastiveVI's salient allocation should not be treated as an inevitable property of the dataset alone.

At the same time, scDisInFact does not dominate all endpoints: its reconstruction MSE is worse on the manuscript's common log-expression metric. The appropriate conclusion is therefore not that one method is globally superior, but that architectural choice materially changes reconstruction, latent stability and where study/condition information is allocated. The experimental confounding remains the dominant barrier to causal cross-drug interpretation.
