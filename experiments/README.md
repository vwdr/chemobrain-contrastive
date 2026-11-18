# experiments/ — KL-utilization (free bits) sweep and calibration

This package tests one idea: that MC-ContrastiveVI's treatment latents are under-used (posterior collapse), and
that forcing them to be used with **free bits** changes how well integrated gradients recover genes we injected
into control cells. Everything here is a thin layer over `scripts/31` (the semi-synthetic benchmark) and
`scripts/41-44` (the baseline harness); no result in the original analysis is overwritten.

## Read this first: the data must exist

`data/processed/recovered_counts.h5ad` is **not in the repository and not in the Zenodo archive.** Every real
run (not `--smoke`) needs it, because the semi-synthetic base cells are built from it. Regenerate it:

```
python scripts/05_download_research_data.py        # downloads GEO data + MSigDB files (several GB)
sha256sum data/raw/GSE216146_chemo_brain.h5ad.gz   # 541822fb6983a4a07a2a4166479228078373c62153b4b9b237914e490bfaa04d
sha256sum data/raw/GSE271055_RAW.tar               # 4d31e82d33e04b842b09f2744fe9a13b2e49afb2b9bc8ae5091b34da4e1c1ba1
python scripts/05_recover_cohort.py
sha256sum data/processed/recovered_counts.h5ad     # 62a5047e74d4ea704f2a10b5b91b35a4836ccefc99186d4c72101db78cbd8ac1
```

(`shasum -a 256` on macOS.) These hashes are from `REPRODUCE.md` in the Zenodo archive (section 2). If the final
hash differs, stop: results would not be comparable with the archived ones. Then build the semi-synthetic base
(writes `runs/peerreview_20261003/semisynthetic/base.npz` and the `D1_*` files; takes minutes):

```
python scripts/31_semisynthetic.py --prepare
```

## Environments

* **Main** (everything except multiGroupVI): Python 3.12, `requirements-research.txt` + `scvi-tools==1.5.1` + `pytest`.
  That file pins `torch==2.5.1+cpu`, a Linux/Windows build; on macOS install `torch==2.5.1`.
* **multiGroupVI** (exp3 only): `.venv-mgvi`, built as in `REPRODUCE.md` section 1 (scvi-tools 0.18.0, a different
  stack). It cannot share the main environment.

## What each experiment answers

| script | question | size |
|---|---|---|
| `exp6_fsh_noise_null` | What shared fraction F_sh does pure noise give for a given shared:drug dimension split? Is the observed null F_sh (0.705, 0.802, 0.869) just a dimension-count artefact? Simulation, no training. | 13 cases x 10,000 draws |
| `exp5_auroc_rescore` | Re-scores every existing IG ranking with AUROC / average precision against the injected genes, instead of top-50 counts. No training. | one job per `*_ig_summary.csv.gz` under `runs/` |
| `exp1_freebits_sweep` | **Headline.** Within one architecture, does forcing KL utilization (free bits 0 .. 1) change active units and recovery, while F_sh stays put? Up to 300 epochs, early-stopping patience 30; every record has `epochs_run`, `best_epoch`, `converged`. | 6 free_bits x 3 configs x 3 seeds = 54 fits |
| `exp7_epoch_ablation` | **Why exp1 now trains up to 300 epochs.** All 75 archived fits stopped at the 100-epoch cap with the best validation MSE at epochs 96-99, i.e. still improving. Does the collapse persist once the model has converged? free_bits 0, `shared_d1_rf1.0`, epoch caps 100/200/300/500. | 4 caps x 3 seeds = 12 fits |
| `exp2_null_calibration` | How variable is F_sh across fits when nothing is injected? Gives the null to test the 72 injected fits against. | 30 null fits |
| `exp3_effect_grid` | Do contrastiveVI and multiGroupVI recover injected genes at small effect sizes (delta 0.25, 0.5)? | contrastiveVI 36 jobs + multiGroupVI 18 jobs |
| `exp4_more_methods` | Do scVI and PCA(32)+Wilcoxon, scored on the same benchmark, do any better? | 18 + 18 jobs |
| `collect.py` | One table (`master.csv`) of every fit, new and archived, plus `STATS.md` (dose-response, null calibration, noise null, cross-method AUROC tests, caveats). | seconds |
| `make_figures.py` | `fig7_freebits.pdf`, `fig8_nullcalib.pdf` (7.16 in wide, 300 dpi; needs `paperstyle.py`). | seconds |

Free bits (`src/models/mc_contrastive_vi.py`): the per-dimension KL of the shared and drug blocks (never the
background) is floored at `free_bits` nats, averaged over the cells of the block, before summing. Defaults
(`free_bits=0`, `kl_warmup_epochs=0`, `treatment_kl_beta=1`) reproduce the original loss exactly
(`tests/test_free_bits.py`). The sweep only varies `free_bits`.

## Commands, in order

Run from the repository root with the main environment active.

```
# 0. sanity (seconds, CPU, no data needed)
python -m pytest tests/test_free_bits.py -v
python experiments/run_all.py --smoke --tiers 0,1,2,3,4      # ~10 s; multiGroupVI jobs FAIL here by design (main env)

# 1. real run, most valuable first (see timings below); add --max-hours N to cap a session
python experiments/run_all.py --tiers 0
python experiments/run_all.py --tiers 1 --resume             # exp1 then exp7 (the convergence check)
python experiments/run_all.py --tiers 2 --resume
python experiments/run_all.py --tiers 3 --resume             # contrastiveVI part; multiGroupVI jobs FAIL until step 2
python experiments/run_all.py --tiers 4 --resume

# 2. multiGroupVI half of exp3, in the other environment
.venv-mgvi/bin/python experiments/exp3_effect_grid.py --method multigroupvi --resume

# 3. tables, statistics, figures
python experiments/collect.py
python experiments/make_figures.py
```

`run_all.py --tiers 0,1,2,3,4 --resume --max-hours 12` does steps 1 in one go. Outputs:
`runs/experiments/<exp>/` (per-fit files and `jobs/<exp>/*.json`; exp7 also writes `epoch_ablation*.csv`), `analysis/experiments/` (exp5, exp6, master.csv,
STATS.md, figures). Each experiment's `progress.log` has one line per job (timestamp, experiment, job id, status,
seconds); `errors.log` has the tracebacks.

For `collect.py` to compare against the archived fits, `runs/peerreview_20261003/semisynthetic/` (the 75 archived
fits' `done.json` and `ig_summary.csv.gz`) and `analysis/peerreview_20261004/R4_*` must be present. They are in
the repository; if a source is missing `collect.py` says so and STATS.md lists it under Caveats.

## Expected wall-clock

Estimated from the archived runs (Windows, Intel i5-10400, CPU only, 3 threads, fits run one to three at a time),
**not measured on this machine**. Time the first job of each tier and scale.

| tier | experiment | jobs | per job (archived) | sequential total |
|---|---|---|---|---|
| 0 | exp6 | 13 | ~2 s (measured here: 1 min total) | ~1 min |
| 0 | exp5 | ~75 | seconds | minutes |
| 1 | exp1 | 54 | archived 100-epoch fit ~165 s + ~7 s IG; up to 3x that if a fit runs the 300-epoch cap | 2.6 h (all stop at 100) to ~8 h (all hit 300) |
| 1 | exp7 | 12 | 100/200/300/500 epochs, ~1.7 s per epoch | ~1.7 h if none stops early |
| 2 | exp2 | 30 | ~172 s | ~1.4 h |
| 3 | exp3 contrastiveVI | 36 | ~90-105 s | ~1 h |
| 3 | exp3 multiGroupVI | 18 | ~1,080 s | ~5.5 h |
| 4 | exp4 scVI | 18 | ~145 s | ~45 min |
| 4 | exp4 PCA+Wilcoxon | 18 | tens of seconds | ~10 min |

Roughly 11 hours in total if exp1 plateaus near 100 epochs, up to ~17 hours if every exp1 fit needs 300; multiGroupVI is the largest single item. Fits are CPU-bound; three concurrent `--threads=N` runs
of different tiers (different `--out`) fit the pattern the archive used.

## Disk

Outputs are small because exp1/exp2 do not save checkpoints or latents: ~0.2 MB per fit (about 20 MB for 84 fits).
exp3/exp4 need injected input files (~15-20 MB each, 6 configurations each, shared between experiments when they
are rebuilt in the same `--out`) and exp3 saves model folders: budget **about 1 GB** for `runs/experiments/`.
Separately, the data download behind `recovered_counts.h5ad` is several GB, and the two environments take roughly
1.5 GB (main) plus the multiGroupVI one. Free disk before starting; a full disk makes job files fail to write.

## Resuming after a crash or kill

Re-run the same command with `--resume`. A job is done if `runs/experiments/<exp>/jobs/<exp>/<job>.json` exists;
done jobs are skipped, failed jobs (no JSON) are retried, and a fit that finished training but not scoring reuses
its `*_done.json` marker instead of retraining. Without `--resume`, an experiment with existing results refuses to
start rather than mixing runs. `kill -TERM <pid>` (or `--max-hours`) finishes the current job, writes it, and exits 0.
`kill -9` loses at most the job in progress (results are written atomically). A failed job never stops the run;
read `errors.log`.

## Known limits

* The multiGroupVI code path (`exp3 --method multigroupvi`) has not been executed by anyone yet: its smoke run
  needs `.venv-mgvi`, which was not built.
* Real-data (non-smoke) runs of exp1-exp4 have not been executed; the smoke runs use 200 synthetic cells and
  50 genes and prove the code paths, not the numbers. exp6 and `collect.py` on the archived fits have run on real inputs.
* **Regression gate.** `collect.py` scores recovery per study group with one convention for every method (two IG baselines
  averaged gene by gene, then scored once; truth = injected S u that group's own set; top-50 enrichment and AUROC always
  together; MC-ContrastiveVI group ranking = mean of its shared axis and the group's own drug axis) and refuses to write
  `STATS.md` unless the archived fits give these group-level medians (delta 1 and 2, responder fraction 1.0, n = 36 groups
  each), AUROC / enrichment: MC-ContrastiveVI 0.461 / 0.75, contrastiveVI 0.479 / 1.50, multiGroupVI 0.528 / 3.75.
  They reproduce. STATS.md also holds the paired one-sided Wilcoxon tests against MC-ContrastiveVI and a six-readout
  robustness table for MC-ContrastiveVI.
* `collect.py` reports archived fits that ran the full 100 epochs with the best validation MSE in the last
  epochs (all 75 archived semi-synthetic fits do). That is a convergence caveat for the whole comparison.
