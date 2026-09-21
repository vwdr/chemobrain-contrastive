# Audit of Tasks 1 and 2

Date: 2026-09-21

## Task 1: latent utilization and execution-context sensitivity

**Verdict: PASS.**

The exact archived corrected checkpoints were reloaded and the Task 1 latent-utilization calculations were rerun. Checkpoint and input SHA-256 hashes matched the recorded provenance. The pooled shared fractions reproduced at 0.920568, 0.957295 and 0.922368 for seeds 0, 1 and 2. The active-unit audit reproduced 0/8, 0/8 and 1/8 active shared coordinates, with no active doxorubicin- or cisplatin-specific coordinate in the archived checkpoints and no salient coordinate with mean KL > 0.01 nats.

The Linux reproduction artifact was also checked against the archived canonical input matrix: the feature matrix, genes, labels and train/validation/test indices were identical. The fresh-run result therefore reflects optimization/execution sensitivity rather than an unnoticed split or gene-universe change.

Scope clarification: the posterior-collapse audit concerns the salient shared and condition-specific blocks. It does not establish collapse of the always-active background block. Manuscript wording should therefore use "strong under-utilization of the salient blocks" or "partial posterior collapse," not whole-model collapse.

## Task 2: library-label exchangeability reference

**Verdict: REQUIRED CORRECTION; PASS AFTER CORRECTION.**

The initial Task 2 reference used 24 assignments by grouping the four D21 non-rescue cisplatin-study libraries into one exchangeability stratum. Audit of the archived GEO series record showed that this was broader than the source experimental structure. GSE216146 documents three experimental replicate blocks:

- replicate 1: D20-6407 control / D20-6409 cisplatin;
- replicate 2: D21-2746 control / D21-2750 cisplatin;
- replicate 3: D21-2747 control / D21-2751 cisplatin.

The corrected reference exchanges labels only within those three replicate pairs and exchanges the two pooled GSE271055 control/doxorubicin libraries. The exact reference therefore contains 2^3 x 2 = 16 assignments, not 24.

All 16 corrected assignments were already present in the original 24-assignment superset. Because every assignment/seed fit resets the PyTorch and NumPy random seeds independently, filtering the existing computation is equivalent to fitting those 16 assignments directly; no model refitting is required.

Corrected primary result:

- observed mean pooled F_sh across seeds: 0.776222;
- combined replicate-preserving reference: 12/16 assignments at least as large, exact upper-tail reference P = 0.750;
- reference median 0.887638, range 0.292013-0.988365;
- cisplatin-only with doxorubicin fixed: 5/8, P = 0.625;
- doxorubicin-only with cisplatin fixed: 2/2, P = 1.0.

The qualitative conclusion is unchanged: the observed labels do not produce an unusually large shared fraction under this model-specific sensitivity reference. These values are not randomized-treatment-effect P values.

## Carryover to Task 3

The same source-metadata finding applies to the rescue permutation analysis. Any Task 3 calculation based on a D20/D21 12-allocation scheme must be re-audited before use. The rescue analysis should respect the three GEO replicate pairs rather than treating all four D21 libraries as one exchangeability stratum. No Task 3 inference should be finalized from the previously proposed 12-allocation resolution.
