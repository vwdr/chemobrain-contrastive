# Peer-review run summary (branch `peerreview-20261003`)

Branched from `main` at `7288dc88a3ebe95c5ab914d5228063eda44c786a`. Instructions:
`D:/chemobrain-contrastive-main/CLAUDE_CODE_PROMPT.md`. Run: 2026-10-03 21:49 to 2026-10-04 (EDT), Windows 11, Intel
i5-10400 (12 logical CPUs), CPU only. Step log with commands and times: `PROGRESS.md`. Environment: `environment.json`
(main venv) and `C4_mgvi_pip_freeze.txt` (multiGroupVI venv). No push, merge, rebase, reset, stash or remote change.
No stop condition was triggered.

## Task status

| task | status | output note | main tables |
|---|---|---|---|
| Phase 0 setup and data | done | `PROGRESS.md` | `environment.json` |
| Phase 1 ved_archive check | done | `A0_ved_archive_check.md` | `A0_*.csv` |
| Phase 2 fresh canonical fits | done | `A1_fresh_canonical_fits.md` | `A1_*.csv`, `A1_fresh_consensus.json` |
| C1 marker lists | done (one marker category limited, see deviations) | `C1_marker_comparison.md` | `C1_marker_list.csv`, `C1_marker_sources.json`, `C1_marker_overlap.csv`, `C1_consensus_gene_annotation.csv` |
| C2 shuffled-label attribution | done | `C2_shuffled_label_attribution.md` | `C2_*.csv` |
| C3 scVI / contrastiveVI | done | `C3_baselines.md` | `C3_*.csv`, `C3_design.json` |
| C4 multiGroupVI | done (third attempt; time box ≈ 59 min of effort) | `C4_multigroupvi.md` | `C4_*.csv` |
| C5 control-cell identity | done | `C5_control_cell_identity.md` | `C5_*.csv` |
| C6 Ttr distribution | done | `C6_ttr_distribution.md` | `C6_*.csv` |
| Phase 4 semi-synthetic positive control | done (full grid, 100 + 100 permutations, calibration) | `D1_semisynthetic.md` | `D1_*.csv`, `D1_design.json`, `D1_calibration.json` |

## Stop-condition checks

| condition | result |
|---|---|
| S1 environment | all pins in requirements-research / task4 / peerreview-baselines installed as pinned; no relaxed pin |
| S2 raw input SHA-256 | both match `data/evidence/source_file_manifest.csv` |
| S3 recovered cohort | SHA-256 of `recovered_counts.h5ad` identical to `62a5047e…8ac1`; regenerated CSVs identical to committed |
| S4 fresh MSE in [0.165, 0.172] | 0.168857 / 0.168697 / 0.168616 |
| S5 protected paths | no protected file modified by code written in this run (final `git status` check below) |
| S6 disk ≥ 20 GB | ≥ 377 GB free throughout |
| S7 Phase 4 projection ≤ 48 h | projected ≈ 4.8 h; no reductions |

## Headline numbers

| task | quantity | value |
|---|---|---|
| A0 | ved `inputs.npz`, `cells.csv` vs regenerated | content-identical (X max abs diff 9.5e-7) |
| A0 | ved `full_0.pt` vs committed `full_0_metrics.json` | differs: best epoch 97 vs 99; F_sh 0.7542 vs 0.9206; NLL 187.05 vs 181.08 |
| A1 | fresh test MSE (seeds 0/1/2) | 0.168857 / 0.168697 / 0.168616 (committed 0.168651 / 0.168801 / 0.168680) |
| A1 | fresh pooled F_sh | 0.9237 / 0.9679 / 0.9705 (committed 0.9206 / 0.9573 / 0.9224) |
| A1 | active units (Var > 0.01) shared/dox/cis | 0/0/0 in all three fresh seeds (committed: seed 2 shared 1) |
| A1 | fresh vs frozen consensus Jaccard (shared/dox/cis) | 0.644 / 0.205 / 0.362 (sizes 40/27/31 vs 34/20/33) |
| C1 | frozen shared list vs BAM markers | 9/20 overlap, BH q = 1.1e-9; vs PanglaoDB macrophages 7/45, q = 3.0e-4 |
| C1 | frozen doxorubicin list vs nucleus-enriched (Bakken 2018) | 6/32, q = 2.7e-5 |
| C1 | frozen cisplatin list vs BAM / fibroblasts | 4/20, q = 0.011 / 5/52, q = 0.035 |
| C2 | shuffled-label F_sh (seeds 0/1/2) | 0.1245 / 0.8674 / 0.9786 |
| C2 | shuffled vs fresh real consensus Jaccard (shared/dox/cis) | 0.015 / 0.188 / 0.133 |
| C2 | ZHANG_UTERUS_C5_MACROPHAGE (top M8 term of paper shared list) in shuffled shared list | 0 overlap, q = 1 |
| C3 | canonical-scale test MSE | MC-CVI 0.1687; scVI 0.2061; contrastiveVI cis pair 0.1260 (MC-CVI same cells 0.0917), dox pair 0.2622 (MC-CVI 0.2028); PCA(32) 0.1552 |
| C3 | scVI within-study treatment probes | cisplatin 0.809; doxorubicin 0.628 (means over seeds) |
| C3 | contrastiveVI salient active units | 8/8 in every fit; total KL 7.0–8.1 nats (target cells) |
| C5 | ungated shared → source cell type (control / cisplatin) | 0.819/0.780/0.717 and 0.760/0.718/0.700 (chance 0.083) |
| C6 | Ttr fraction nonzero, Choroid Plexus vs other (control arm) | 1.000 vs 0.428; mean log-norm 6.655 vs 0.499 |
| D1 | null F_sh (seeds 0/1/2) | 0.7054 / 0.8688 / 0.8016 |
| D1 | F_sh range across all 24 non-null configurations (config means) | 0.569 – 0.934 |
| D1 | IG top-50 hits of injected set, per fit (150 rankings), mean (max) | shared→S 2.14 (4); cis→A 1.48 (4); dox→B 1.23 (3); random expectation 1.33 |
| D1 | label exchange P | null 0.772; shared_d1_rf1.0 0.139 |
| D1 | calibration median abs log2FC | 0.742 (all genes) / 0.981 (1,500-gene universe) |
| C4 | multiGroupVI canonical / modeled-gene MSE (seeds 0/1/2) | 0.2453 / 0.2429 / 0.2480 and 1.476 / 1.440 / 1.471 |
| C4 | multiGroupVI within-study treatment probes, gated group-specific latents | cisplatin 0.950 / 0.963 / 0.953; doxorubicin 0.901 / 0.894 / 0.881 |
| C4 | multiGroupVI cross-seed CKA shared / gated group-specific | 0.755–0.767 / 0.364–0.440 |

## Deviations from the instructions and other decisions

1. Order of execution: C2 fits were run in the background before C1/C5 were finished (CPU otherwise idle); C6 was run
   during the Phase 4 queue (CPU-light); C4 installation (no training) was done during the Phase 4 queue. All C4 fits ran
   after Phase 4 finished. At most 3 model fits ran at any time.
2. Phase 2 fits used 2 threads (the canonical `torch.set_num_threads(2)`); all other fits used 3.
3. C2: the within-study shuffle was applied to training, validation and test cells (each split separately), so early
   stopping, gating, F_sh and attribution-cell selection use one consistent label assignment.
4. C1: no source was found that gave a meningeal/perivascular-fibroblast versus all-cell-type marker table other than
   the Zeisel et al. 2018 VLMC/ABC cluster markers (5 genes per cluster; 15 listed, 6 in the universe). Several sources
   were inspected and not used (listed in `C1_marker_sources.json`). The PMC copy of Dani et al. Table S1 sat behind a
   proof-of-work bot challenge, which was not bypassed; the publisher copy was used instead. BAM rule: intersection of the
   female-control and male-control BAM top-30 lists of Ochocka et al. 2021, plus Zeisel PVM1/PVM2 markers (rules fixed
   before any overlap was computed; list committed in 3690ba7 before the overlap).
5. C3: additions beyond script 20 (requested): canonical-scale MSE conversion, MC-ContrastiveVI and PCA(32) on both
   scales, plus a PCA(32) fitted directly on modeled-gene-scale data and MC-ContrastiveVI MSE restricted to each pair's
   test cells.
6. C4: see `C4_multigroupvi.md` for environment workarounds (pip 24.0 inside `.venv-mgvi` to install pytorch-lightning
   1.7.7; pinned older dependency versions; Python 3.12 is the only interpreter on the machine) and for the explicit
   `max_epochs=400` argument (the package's default rule value passed as a Python int).
7. Phase 4: thinning random generator `default_rng(20261003)` per configuration (seed not specified in the
   instructions); floor(n/2) pseudo-treated cells per library; permutation fits keep the canonical split fixed and skip
   IG. Added diagnostic `scripts/35_d1_signal_check.py`. Trace-of-covariance true shared fraction: summed variances below
   1e-12 set to 0 (corrected after the first scoring pass; no refit).
8. Calibration: regenerated PyDESeq2 counts differ from the committed Task 4 summary in two rows (Microglia 101 vs 109,
   Oligodendrocyte 141 vs 137 cisplatin_vs_control hits); 69 supported universe genes vs 68 asserted in script 19.
9. Phase 1: HVG reselection on this machine returned 1,499/1,500 committed genes (Cebpd instead of Aspm); the committed
   universe was used throughout.
10. Thinning unit test run with the venv Python directly (pytest is not part of the pinned environment).
11. Additional packages installed in `.venv` for reading supplementary spreadsheets: openpyxl 3.1.5, xlrd 2.0.2 (not in
    any requirements file; requirements files unchanged).
12. Times in early `PROGRESS.md` entries were corrected from file modification times.
13. 200 Phase 4 permutation log files were renamed: a trailing carriage return in the job list had become part of their
    file names (the arguments themselves were parsed correctly; output files are named by the script).

## Relaxed package pins

None in `.venv` (requirements-research.txt, requirements-task4.txt, requirements-peerreview-baselines.txt installed as
pinned). The separate `.venv-mgvi` follows multiGroupVI's own requirement (scvi-tools==0.18.0, protobuf<=3.20.1) with
compatible versions of transitive dependencies chosen by hand (`C4_mgvi_pip_freeze.txt`).

## Run times (wall clock)

| step | time |
|---|---|
| environment install | 21:49–21:56 |
| download + cohort recovery | 21:53–21:57 |
| Phase 1 | 21:56–21:58 |
| Phase 2 fits (3 concurrent) + analysis | 22:00–22:04 (≈178 s per fit, IG 4–6 s) |
| C2 fits | 22:06–22:09 (189 s per fit) |
| C1 retrieval + overlap | 22:05–22:14 |
| C3 fits (9) | 22:17–22:26 (scVI 144 s, contrastiveVI 102–106 s) |
| Phase 4 grid (75 fits) | 22:26:49–23:45 (mean 166 s per fit) |
| Phase 4 permutations (200 fits) | 23:45–02:54:57 |
| C4 | install 22:30:53–22:49:56; fits 02:55:08–03:34:38 (attempt 1 failed at start; attempt 2 failed after training in post-processing; attempt 3: 1,077–1,079 s per fit) |

## Final git status check

Before the final commit, `git status` showed no modified tracked file outside the new locations, and
`git diff main --name-only` lists only `analysis/peerreview_20261003/`, `scripts/23_*`–`scripts/35_*`, `src/peerreview/`
and `tests/test_peerreview_thinning.py`. Restore-procedure events: `data/evidence/GSE216146.soft`, `GSE271055.soft`,
`GSE286221.soft` were overwritten by `05_download_research_data.py` (identical except line endings) and restored with
`git checkout --`; the three CSVs rewritten by `05_recover_cohort.py` were identical to the committed versions.
Large or local-only outputs are excluded via `.git/info/exclude` (`.venv-mgvi/`, `runs/peerreview_20261003/scratch/`,
`runs/peerreview_20261003/c1_sources/`, two pseudobulk result files > 5 MB, `canonical_inputs/cells.csv`); `.pt` and
`.npz` files are git-ignored by the existing `.gitignore`.

