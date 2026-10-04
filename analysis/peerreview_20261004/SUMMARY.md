# Run A summary (outputs `analysis/peerreview_20261004/`, `runs/peerreview_20261004/`)

Instructions: `D:/chemobrain-contrastive-main/CLAUDE_CODE_PROMPT_RUN_A.md`. Branch `peerreview-20261003` (started at commit
e82fc38). Run 2026-10-04 13:18–16:20 EDT, Windows 11, Intel i5-10400, CPU only, at most 3 concurrent fits with ≤ 3 threads.
No push. No stop condition (T1–T4) was triggered. No figures were produced. Step log: `PROGRESS.md`; file index: `INDEX.md`;
environments: `environment_notes.md`.

## Task status

| task | status | note | main tables |
|---|---|---|---|
| R1 injected-signal location | done | `R1_injected_signal_location.md` | `R1_probes_per_fit.csv`, `R1_reconstruction_per_fit.csv`, `R1_summary*.csv` |
| R2 baseline headline table | done (from existing files; no refits) | `R2_baseline_headline.md` | `R2_baseline_headline.csv` |
| R3 consistent numbers | done | `R3_consistent_numbers.md` | `R3_paper_numbers.csv`, `R3_comparison_models*.csv`, `R3_09_interpretation/`, `R3_11_diagnostics/`, `R3_enrichment_*.csv`, `R3_pseudobulk_*.csv`, `R3_consensus_top12_figure5_data.csv` |
| R4 baseline positive control | done (63/63 fits; attribution implemented for all three models) | `R4_baseline_positive_control.md` | `R4_*.csv` |
| R5 archive package | done (prepared, not uploaded) | `REPRODUCE.md` | `MANIFEST.csv`, `ARCHIVE_ZIPS.csv`, `EXCLUDED_INPUT_HASHES.csv`, `pip_freeze_*.txt` |
| R6 wrap-up | done | this file | `INDEX.md` |

## Headline numbers

| task | quantity | value |
|---|---|---|
| R1 | within-study probes (pseudo-treated vs pseudo-control) on bg / ungated shared / ungated drug, all 25 configurations, mean over seeds | 0.473–0.541 |
| R1 | responder-vs-other probes (label b) | 0.425–0.548 |
| R1 | injected sets: observed responder-minus-other effect vs decoder-reconstructed effect (log1p units, means over seeds) | observed −0.144 to +0.158; reconstructed −0.033 to +0.030 |
| R2 | canonical MSE: MC-ContrastiveVI / PCA(32) / scVI / contrastiveVI cis pair / dox pair / multiGroupVI / scDisInFact | 0.1687 / 0.1552 / 0.2061 / 0.1260 / 0.2622 / 0.2454 / 0.2313 |
| R2 | ungated within-study treatment probes, cisplatin / doxorubicin: MC-ContrastiveVI drug block; scVI; multiGroupVI group-specific | 0.530 / 0.554; 0.809 / 0.628; 0.830 / 0.570 |
| R2 | multiGroupVI gated (label-masked) probes, separate column | cisplatin 0.955; doxorubicin 0.892 (means) |
| R3 | ablation / NB held-out MSE ranges (fresh refits) | no HSIC 0.1684–0.1692; no gating 0.1677–0.1678; Gaussian VAE 0.1676–0.1679; NB 0.2108–0.2116; PCA 0.1552 |
| R3 | NB shared fraction | 86.4–88.1% (manuscript 82.9–89.2%) |
| R3 | pooled F_sh, fresh fits (bootstrap intervals in `R3_11_diagnostics/variance_robustness.csv`) | 92.4–97.0% (manuscript 92.1–95.7%) |
| R3 | shared top-50 cross-seed Jaccard, squared norm, zero baseline | 0.449–0.562 (manuscript 0.190–0.639) |
| R3 | rescue-minus-treatment mean distance difference (cisplatin); sign-flip P | −0.0032; 1.000 (manuscript 0.0003; 1.000) |
| R3 | pseudobulk regenerated: cisplatin vs control / pooled / same-direction | 293 / 350 / 349 (manuscript 297 / 354 / 353); Welch 1; Spearman 0.932–0.997 (unchanged) |
| R3 | manuscript numbers listed in `R3_paper_numbers.csv` | 116 rows: 28 changed, 77 unchanged, 11 not recomputable |
| R4 | top-50 hits vs injected union (80 genes; random 2.67), non-null configurations, means over seeds, zero baseline | multiGroupVI 3.67–11.0; contrastiveVI 1.67–7.0; MC-ContrastiveVI 1.0–4.0 |
| R4 | same, null configuration | multiGroupVI 2.33–4.67; contrastiveVI 2.33–4.33; MC-ContrastiveVI 3.33–4.67 |
| R4 | within-study probe on treatment-related latent, non-null configurations | multiGroupVI 0.58–0.80; contrastiveVI 0.52–0.73; MC-ContrastiveVI 0.49–0.54 |
| R4 | active units of treatment-related latent (mean) | multiGroupVI 8.3–10 of 10; contrastiveVI 8 of 8; MC-ContrastiveVI 0–1.3 of 12 |
| R5 | archive | 7 zips, 1.2 GB (largest 540 MB), 1,931 files in `MANIFEST.csv` |

## Deviations and decisions

1. `.venv-baselines` (named in the instructions) does not exist; scvi-tools 1.5.1 is installed in `.venv` (previous run) and was
   used for contrastiveVI. Recorded in `environment_notes.md`. No environment was changed.
2. Order: the R3 comparison fits were launched first (13:20), as R3 permits early launch; R2 and R1 ran while they trained. The R3
   reruns of scripts 09 and 16, the R4 input preparation and the MC-ContrastiveVI attribution (no fitting) also ran while fits
   trained. R4 fits started when the comparison fits freed the slots (13:36).
3. R3 comparison fits use the thread counts set by the original scripts (06: 2 torch threads; 12: 1).
4. R3 reruns of scripts 09, 11 and 16 execute their unmodified source with only the protected path expressions replaced
   (`scripts/39_r3_rerun_canonical_scripts.py`). Inputs were staged by copying from the previous run (read-only sources).
5. R3 scDisInFact top-100 overlap uses the frozen top-100 set in `scripts/19_gene_pattern_audit.py`; the committed Task 6 gene-score
   file contains only seed-pair stability statistics, not per-gene scores.
6. R3 Figure 4 and Figure 6 values and Table/figure-only quantities are listed as not-recomputable in `R3_paper_numbers.csv` with the
   fresh data file; values from analyses independent of the archived fits are the committed values (marked in `notes`).
7. R1: for the null configuration the responder label (b) was computed with both responder sets; non-injected sets are reported with
   `injected = False`.
8. R2: contrastiveVI probes are the C3 values (`max_iter=1500`, script 20 helper); MC-ContrastiveVI cross-seed CKA is available only
   for gated blocks; PCA, scVI and scDisInFact lack some quantities (listed as not available).
9. R4: integrated gradients are taken with respect to each encoder's actual input (multiGroupVI group encoders: raw counts, as the
   package feeds them; contrastiveVI salient encoder: log1p counts; MC-ContrastiveVI: log-normalized expression). 2-epoch smoke-test
   fits were run before the queue and are not part of the results or the archive.
10. Two data-handling fixes: pandas parses the configuration name "null" as missing; R4 collection was re-run with
    `keep_default_na=False` (no refit).
11. Archive: raw supplementary marker-source files are not included (their SHA-256 values are in `C1_marker_sources.json`); pip
    freezes exist for `.venv` and `.venv-mgvi` only.

## Run times (wall clock, 2026-10-04)

| step | time |
|---|---|
| setup | 13:18–13:22 |
| R3 comparison fits (13) | 13:20–13:35 (ablations 129–161 s; NB 257–294 s) |
| R2 | 13:22–13:24 |
| R1 | 13:24–13:27 |
| R3 reruns (09, 16, 11) and analyses | 13:24–13:50 |
| R4 queue (63 fits) | 13:36–16:10 (multiGroupVI 1,068–1,136 s; contrastiveVI 87–93 s) |
| R5 archive | 16:12–16:14 |
| R6 | 16:14–16:20 |

## Final git status

All tracked changes are inside `analysis/peerreview_20261004/`, `runs/peerreview_20261004/` (small text outputs), `scripts/36_*`–`47_*` and `src/peerreview/ig_generic.py`; no file in
the protected folders, in the previous run's outputs, or in `src/peerreview/` / `scripts/23_*`–`35_*` was modified. Large outputs
are excluded via `.git/info/exclude` (`runs/peerreview_20261004/scratch/`, `runs/peerreview_20261004/archive/`); `.pt`/`.npz` are
git-ignored by `.gitignore`.
