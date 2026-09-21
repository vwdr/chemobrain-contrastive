# Task 5: unseen-library reconstruction validation

## Design

The canonical benchmark holds out cells from libraries that also contribute training cells. Task 5 instead removes one complete GSE216146 control/cisplatin GEO replicate pair at a time. The three folds hold out D20-6407/D20-6409, D21-2746/D21-2750 and D21-2747/D21-2751.

The held-out libraries contribute no cells to training, validation or feature selection. Each fold uses 6,000 balanced training cells (1,500 per study-by-treatment stratum), 400 validation cells (100 per stratum) and 1,500 HVGs selected from training cells only. All retained non-rescue cells from the held-out pair are used for testing. The full MC-ContrastiveVI and noncontrastive Gaussian VAE are fit for seeds 0-2. PCA(32) and a training gene-wise mean predictor are deterministic baselines.

The primary endpoint is balanced held-out-library reconstruction MSE: the mean of control-library and cisplatin-library MSE within each fold, with the three folds then weighted equally.

## Results

Across folds, full MC-ContrastiveVI had mean seen-library validation MSE 0.09099 and mean unseen-library MSE 0.09401, a mean absolute gap of 0.00302 and mean unseen/seen ratio 1.0348. The noncontrastive Gaussian VAE was similar (0.08957 seen, 0.09255 unseen; ratio 1.0348). PCA(32) had lower MSE on both seen and unseen libraries (0.08422 seen, 0.08596 unseen; ratio 1.0223). The training-mean baseline was much worse (0.31939 seen, 0.32678 unseen).

For the full model, fold-specific unseen/seen ratios were 1.1113, 1.0433 and 0.9498 for replicate pairs 1, 2 and 3. Thus the first held-out pair showed an approximately 11% reconstruction penalty, the second a roughly 4% penalty, and the third was reconstructed slightly better than the sampled seen-library validation cells. The direction of the gap therefore depends on which deposited libraries are held out.

PCA(32) had lower unseen-library MSE than the full model in every fold: 0.09092 versus 0.09920, 0.08443 versus 0.09309 and 0.08252 versus 0.08973. Averaged across folds, the full model's unseen MSE was approximately 9.4% higher than PCA; the noncontrastive Gaussian VAE's was approximately 7.7% higher. Full MC-ContrastiveVI was approximately 1.6% worse than the noncontrastive Gaussian VAE on this endpoint.

## Interpretation

The full contrastive model does not show catastrophic reconstruction failure when complete GSE216146 libraries are withheld: its mean unseen-library reconstruction error is only about 3.5% higher than its seen-library validation error. This is evidence that basic expression reconstruction transfers beyond cells from the exact deposited libraries used for fitting.

However, this sensitivity analysis does not establish a reconstruction advantage for MC-ContrastiveVI. PCA is better in every fold, and the noncontrastive Gaussian VAE is slightly better than the full model. Moreover, the test supplies the known study and treatment labels to the conditional model; it evaluates reconstruction of expression conditional on those labels, not prediction of treatment identity, recovery of a stable shared-treatment mechanism, or transport to an unseen study.

The analysis is leave-one-deposited-library-pair-out within GSE216146, not leave-one-animal-out. Some deposited libraries pool material from more than one animal. GSE271055 cannot support an analogous treatment-balanced holdout because there is only one pooled control library and one pooled doxorubicin library. These constraints should remain explicit in the manuscript.
