# Peer-review run progress (peerreview-20261003)

Branch `peerreview-20261003` created from `main` at `7288dc88a3ebe95c5ab914d5228063eda44c786a`.
Instructions: `D:/chemobrain-contrastive-main/CLAUDE_CODE_PROMPT.md`.

## Phase 0 — Setup and data

### P0.1 Branch and folders — done
- Start/end: 2026-10-03
- Commands: `git checkout -b peerreview-20261003`; `mkdir -p analysis/peerreview_20261003 runs/peerreview_20261003/{logs,scratch} src/peerreview`

### P0.2 Environment — done (21:49–21:56)
- `py -3.12 -m venv .venv` (Python 3.12.10); `pip install -r requirements-research.txt` then
  `-r requirements-task4.txt -r requirements-peerreview-baselines.txt`. All pins installed as pinned
  (including `torch==2.5.1+cpu`); **no relaxed pins**. requirements-task6.txt not installed.
- Logs: `runs/peerreview_20261003/logs/pip_*.log`; environment: `analysis/peerreview_20261003/environment.json`.

### P0.3 Download — done (21:53)
- `.venv/Scripts/python scripts/05_download_research_data.py` (log `logs/05_download.log`).
- m8.all.v2025.1.Mm.symbols.gmt (not in manifest) SHA-256 `0b11177c438565361ce7dc917d216cf40dd71b12982d7f117b51b182be6a7788`; m2/m5 match manifest.

### P0.4 SHA-256 (S2) — done, pass
- GSE216146_chemo_brain.h5ad.gz `541822fb…a04d` = manifest; GSE271055_RAW.tar `4d31e82d…a1` = manifest.
- Restore procedure: `data/evidence/GSE216146.soft`, `GSE271055.soft`, `GSE286221.soft` showed as modified;
  `git diff --ignore-cr-at-eol` empty → identical except line endings (downloaded files' SHA-256 also equal
  the manifest). Restored with `git checkout --`.

### P0.5–6 Cohort recovery and content check (S3) — done, pass
- `.venv/Scripts/python scripts/05_recover_cohort.py` (log `logs/05_recover.log`): 46,857 × 18,271.
- `recovered_counts.h5ad` SHA-256 `62a5047e74d4ea704f2a10b5b91b35a4836ccefc99186d4c72101db78cbd8ac1` → exact match.
- `cohort_qc.csv`, `source_annotation_counts.csv`, `label_audit.csv` rewritten; `git status` reported no
  modification (identical to committed versions). `git checkout --` not needed (no-op).

## Phase 1 — ved_archive content check — done (21:56–21:58)
- `.venv/Scripts/python scripts/23_prepare_and_ved_check.py` (log `logs/23_ved_check.log`).
- Regenerated canonical inputs in `runs/peerreview_20261003/canonical_inputs/` (6,000/2,881/2,881).
- inputs.npz and cells.csv: content-identical (X max abs diff 9.54e-7). full_0.pt: differs (best_epoch 97 vs 99;
  F_sh 0.7542 vs 0.9206; NLL 187.05 vs 181.08). full_0_latents.npz: consistent with ved full_0.pt, differs from committed.
- HVG reselection on this machine: 1,499/1,500 overlap (Cebpd in, Aspm out); committed universe used.
- Output: `A0_ved_archive_check.md` + 5 CSVs.

## Phase 2 — Fresh canonical fits — done (22:00–22:04)
- `scripts/24_canonical_fits.py --fit --seed {0,1,2} --threads=2` (3 concurrent; 178–179 s/fit, IG 4–6 s), then `--analyze`.
- Outputs: `runs/peerreview_20261003/canonical/`; `A1_*.csv`, `A1_fresh_canonical_fits.md`.
- S4: test MSE 0.168857 / 0.168697 / 0.168616 → pass.
- Fresh F_sh pooled 0.9237/0.9679/0.9705; consensus sizes shared 40, dox 27, cis 31; Jaccard vs frozen 0.644/0.205/0.362.
- Deviation: openpyxl 3.1.5 installed into .venv (needed to read supplementary .xlsx for C1; not in requirements).

## Task C2 fits (started early, before C1 overlap) — fits done (22:06–22:09)
- Deviation from the order in section 6: the three C2 fits were run in the background while C1 sources were
  being retrieved (CPU otherwise idle). C2 analysis is run after C1/C5.
- `scripts/25_shuffled_label_attribution.py --fit --seed {0,1,2} --threads=3`; outputs `runs/peerreview_20261003/shuffled_labels/`.
- test MSE 0.168045/0.168701/0.168329; F_sh 0.1245/0.8674/0.9786.

## Task C1 — marker list built and committed before overlap (22:05–22:13)
- Sources retrieved into `runs/peerreview_20261003/c1_sources/` (see `C1_marker_sources.json` for URLs/SHA-256):
  PanglaoDB 27 Mar 2020 (11 types, mouse entries); Ochocka 2021 Suppl. Data 1 (BAM, female∩male control);
  Zeisel 2018 Table S4 (PVM1/2 → BAM; VLMC1/VLMC2/ABC → meningeal/perivascular fibroblasts; CHOR → CP epithelium);
  Dani 2021 Table S1 Epithelial (CP epithelium); Bakken 2018 S2 Table (159 nucleus-enriched genes).
- Inspected but not usable (no marker table): Van Hove 2019 Suppl. Tables 2-3; Pietilä 2023 mmc2-4; DeSisto 2020
  mmc2-6 (sub-cluster contrasts only); Vanlandewijck 2018 Suppl. Table 3/Fig.1 source data. PMC download of Dani
  Table S1 hit a proof-of-work bot challenge (not bypassed); the publisher CDN copy was used instead.
- Extra packages installed in .venv for reading supplementary files: openpyxl 3.1.5, xlrd 2.0.2.
- `scripts/26_c1_marker_list.py` → `C1_marker_list.csv` (1,421 rows; 15 cell-type sets), `C1_marker_sources.json`.

Note: times above are wall-clock (EDT) taken from file modification times; earlier draft times were corrected.

### C1 overlap — done (22:14)
- `scripts/27_c1_marker_comparison.py` → `C1_marker_overlap.csv`, `C1_consensus_gene_annotation.csv`,
  `C1_marker_set_sizes.csv`, `C1_marker_comparison.md`.
- Frozen shared list: BAM overlap 9/20 (q = 1.1e-9), Macrophages 7/45 (q = 3.0e-4), Microglia 4/24 (q = 8.5e-3);
  frozen dox: nucleus-enriched 6/32 (q = 2.7e-5), Neurons 4/47 (q = 0.021); frozen cis: BAM 4/20 (q = 0.011),
  Fibroblasts 5/52 (q = 0.035).

## Task C5 — done (22:15)
- `scripts/28_c5_control_identity.py` → `C5_*.csv`, `C5_control_cell_identity.md`.
- Ungated shared → source cell type balanced accuracy: control 0.819/0.780/0.717; cisplatin 0.760/0.718/0.700.

## Task C2 analysis — done (22:15)
- `scripts/25_shuffled_label_attribution.py --analyze` → `C2_*.csv`, `C2_shuffled_label_attribution.md`.
- Shuffled consensus sizes shared 27 / dox 30 / cis 20; Jaccard vs fresh real 0.015 / 0.188 / 0.133.
- Paper shared top M8 term ZHANG_UTERUS_C5_MACROPHAGE: shuffled shared overlap 0, q = 1 (not significant).
- Deviation: shuffle applied within study in train, validation and test cells (see C2 note).

## Task C3 — fits running (22:23)
- `scripts/29_baselines_extended.py` (imports scripts/20 unmodified; outputs runs/peerreview_20261003/baselines/), 9 jobs via xargs -P 3, 3 threads each. scVI seeds 0-2 done.

## Phase 4 — preparation (22:23)
- Unit test `tests/test_peerreview_thinning.py` passes (run with `.venv/Scripts/python tests/test_peerreview_thinning.py`; pytest is not installed in the pinned environment).
- `scripts/31_semisynthetic.py --prepare`: 14,285 base cells (PN 4,830; CNT 9,455); split 6,000/1,428/1,429;
  gene sets saved to `D1_gene_sets.csv` and committed before any fitting. Manual vs scanpy normalization max diff 9.5e-7.
- Calibration (`scripts/30_calibration_pseudobulk.py`, script 16 logic redirected) running in background.

## Task C3 — done (fits 22:17–22:26; collect 22:28)
- 9 fits (scVI ×3: 144 s each; contrastiveVI 2 pairs ×3: 102–106 s), `--collect` → `C3_*.csv`, `C3_design.json`, `C3_baselines.md`.
- Canonical MSE: MC-CVI 0.1687, scVI 0.2061, cVI cis-pair 0.1260 / dox-pair 0.2622, PCA(32) 0.1552 (= committed).
- contrastiveVI salient: 8/8 active units, all dims KL > 0.01 in every fit.

## Task C6 — done (22:28) (run early while CPU was otherwise used by fits; CPU-light)
- `scripts/32_c6_ttr_distribution.py` → `C6_ttr_by_celltype_arm.csv`, `C6_ttr_choroid_plexus_vs_other.csv`, `C6_ttr_distribution.md`.

## Phase 4 — timing batch started 22:26:49 (null seeds 0-2, 3 concurrent × 3 threads)
