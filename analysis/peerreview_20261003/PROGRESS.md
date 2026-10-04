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

## Phase 1 — ved_archive content check — done (22:00–22:10)
- `.venv/Scripts/python scripts/23_prepare_and_ved_check.py` (log `logs/23_ved_check.log`).
- Regenerated canonical inputs in `runs/peerreview_20261003/canonical_inputs/` (6,000/2,881/2,881).
- inputs.npz and cells.csv: content-identical (X max abs diff 9.54e-7). full_0.pt: differs (best_epoch 97 vs 99;
  F_sh 0.7542 vs 0.9206; NLL 187.05 vs 181.08). full_0_latents.npz: consistent with ved full_0.pt, differs from committed.
- HVG reselection on this machine: 1,499/1,500 overlap (Cebpd in, Aspm out); committed universe used.
- Output: `A0_ved_archive_check.md` + 5 CSVs.
