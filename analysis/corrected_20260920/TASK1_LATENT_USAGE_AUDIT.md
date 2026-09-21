# Task 1: corrected latent-utilization and posterior-collapse audit

## Scope

This audit evaluates the three corrected full MC-ContrastiveVI checkpoints used for the manuscript and fresh reruns of the same fixed benchmark. The benchmark uses the recovered corrected cohort, the committed 1,500-gene universe, the fixed train/validation/test split, and seeds 0, 1, and 2.

For each held-out treated cell and latent coordinate, the audit records:

- posterior-mean variance, `Var_x[E_q(z_j|x)]`;
- an active-unit flag defined for this audit as posterior-mean variance > 0.01;
- mean coordinate-wise KL divergence to the standard normal prior;
- a separate KL sensitivity flag at mean KL > 0.01 nats.

The KL flag is not used as the active-unit definition.

## Archived corrected checkpoints

The shared posterior-mean variance fraction recomputes exactly to the manuscript values: 0.920568, 0.957295, and 0.922368 for seeds 0, 1, and 2.

However, latent utilization is weak:

| Seed | Shared F_sh | Shared active units / 8 | Dox active units / 4 | Cis active units / 4 | Shared total mean KL (nats) | Dox total mean KL | Cis total mean KL |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.920568 | 0 | 0 | 0 | 0.010691 | 0.000408 | 0.003091 |
| 1 | 0.957295 | 0 | 0 | 0 | 0.007908 | 0.000591 | 0.001364 |
| 2 | 0.922368 | 1 | 0 | 0 | 0.021853 | 0.001172 | 0.006779 |

No coordinate in any archived checkpoint exceeded the separate mean-KL > 0.01-nat sensitivity threshold. Thus the large shared-fraction statistic occurs while most salient coordinates show very small posterior-mean variation and KL divergence. The ratio should therefore not be interpreted as evidence that a large, strongly utilized common biological program dominates the representation.

The exact archived checkpoint SHA-256 hashes used for this audit are recorded in `task1_committed_checkpoint_hashes.json`.

## Fresh CPU reruns

Three independent GitHub-hosted Ubuntu x86-64 jobs retrained seeds 0, 1, and 2 with the same data hash, committed gene universe, Python 3.12.14, PyTorch 2.5.1+cpu, NumPy 2.5.3, Scanpy 1.11.5, AnnData 0.12.4, and scikit-learn 1.9.1.

Across the nine reruns:

- held-out MSE remained in a narrow range, 0.168516-0.168816;
- pooled shared fraction ranged from 0.749563 to 0.981290;
- eight of nine fits had 0/8 active shared dimensions and one had 1/8;
- seven of nine fits had 0/8 active drug-specific dimensions in total and two had 1/8;
- only one fit had any coordinate above the mean-KL > 0.01-nat sensitivity threshold.

Two runners used AMD EPYC 7763 CPUs and produced identical results for every seed. A third runner used an AMD EPYC 9V45 CPU and produced different latent allocations despite nearly unchanged reconstruction error. This pattern is consistent with sensitivity of the optimization trajectory and latent partition to low-level numerical/execution context. It does not by itself prove that CPU model is the sole cause.

## Interpretation

The corrected models reconstruct expression reproducibly at the level of held-out MSE, but the internal shared-versus-condition-specific allocation is substantially less stable. The active-unit and KL diagnostics indicate that most salient coordinates are close to inactive under the evaluated thresholds. Consequently, a high `F_sh` can be a ratio formed from weakly utilized latent blocks and can change substantially without a commensurate change in reconstruction.

For the manuscript, `F_sh` should remain a descriptive coordinate-partition statistic. It should not be described as a stable estimate of common chemotherapy biology, percentage of toxicity, or evidence that most treatment-associated signal is shared across drugs.
