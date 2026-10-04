# Run A progress (peerreview-20261003 branch, outputs peerreview_20261004)

Instructions: `D:/chemobrain-contrastive-main/CLAUDE_CODE_PROMPT_RUN_A.md`. Start 2026-10-04 13:18 EDT at commit e82fc38.

## Setup — done (13:18–13:22)
- Folders `analysis/peerreview_20261004/`, `runs/peerreview_20261004/{logs,scratch}`; scratch and archive added to `.git/info/exclude`.
- `.venv` imports torch/scanpy/anndata/sklearn/scvi-tools 1.5.1/pydeseq2; `.venv-mgvi` imports scvi-tools 0.18.0 + multigroup_vi.
- `.venv-baselines` does not exist (see environment_notes.md); `.venv` is used for scvi-tools baselines, as in Task C3.

## R3.1 comparison fits — launched 13:20 in background (scripts/36_r3_comparison_fits.py; 13 jobs, xargs -P 3)

## R2 — done (13:21)
- `scripts/37_r2_baseline_headline.py` → `R2_baseline_headline.csv`; note `R2_baseline_headline.md`. No refits.

## R1 — done (13:27–13:24)
- `scripts/38_r1_injected_signal.py` → `R1_probes_per_fit.csv`, `R1_reconstruction_per_fit.csv`, `R1_summary*.csv`, `R1_injected_signal_location.md`.
- Label-(a) probe means 0.473–0.541 across configs/studies/representations; injected-set reconstructed effects −0.033 to +0.030 vs observed −0.144 to +0.158.

## R3 — done (13:20–13:42)
- Comparison fits (13 jobs, 13:20–13:35; all exit 0) → runs/peerreview_20261004/comparison/.
- scripts/39 reruns of 09 (13:24–13:25), 16 (13:24–13:34), 11 (13:38) with redirected paths → R3_09_interpretation/, R3_11_diagnostics/, runs/.../pseudobulk_task4/.
- scripts/40 parts comparison/consensus/enrichment/pseudobulk; scripts/45 → R3_paper_numbers.csv (116 rows: 28 changed, 77 unchanged, 11 not recomputable); note R3_consistent_numbers.md.
- Pseudobulk regenerated: cisplatin vs control 293 (committed 297), pooled 350 (354), same-direction 349 (353); other quantities unchanged.

## R4 — queue launched 13:37 (21 multiGroupVI + 42 contrastiveVI fits, xargs -P 3); MC-ContrastiveVI attribution on grid checkpoints done 13:36 (scripts/44 --mccvi).
- Smoke tests (2 epochs) of scripts/42 and 43 passed (outputs in runs/peerreview_20261004/r4/*/smoke/).

## R4 — done (queue 13:36:15–16:10; 63/63 fits exit 0, first attempt)
- multiGroupVI 21 fits (1,068–1,136 s), contrastiveVI 42 fits (87–93 s), MC-ContrastiveVI attribution on existing grid checkpoints.
- scripts/44 --collect → R4_*.csv; note R4_baseline_positive_control.md. Fix in collect: configuration name "null" read as NaN by pandas; reads use keep_default_na=False.

## R5 — done (16:14)
- scripts/46_r5_archive.py → runs/peerreview_20261004/archive/ (7 zips, 1.2 GB; largest semisynthetic.zip 540 MB), MANIFEST.csv (1,931 rows), ARCHIVE_ZIPS.csv,
  EXCLUDED_INPUT_HASHES.csv, pip freezes (.venv, .venv-mgvi; .venv-baselines does not exist), REPRODUCE.md. Archive folder in .git/info/exclude; copies committed in this folder.
- Not archived: raw GEO data, recovered_counts.h5ad, MSigDB files (hashes recorded), smoke-test outputs, staging copies of inputs, raw supplementary marker-source files (SHA-256 in C1_marker_sources.json).
