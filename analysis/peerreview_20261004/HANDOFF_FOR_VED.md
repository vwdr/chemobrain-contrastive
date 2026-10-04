# Handoff for Ved — peer-review analyses for the MC-ContrastiveVI paper

Prepared by Manan (with Claude), 2026-10-04. Everything below was produced by code on branch
`peerreview-20261003` of `vwdr/chemobrain-contrastive`. No manuscript text and no figures were written —
those are yours. This note explains what was done, what the results show, which decisions are yours, and
what we recommend.

---

## 1. Summary

- Every analysis the reviewer asked for has now been run.
- **Positive control.** On semi-synthetic data with a known injected treatment effect, MC-ContrastiveVI never
  represents the injected signal:
  - its treatment blocks stay collapsed;
  - F_sh is as high when *no* shared program exists as when everything is shared;
  - integrated gradients recover the injected genes no better than the null configuration does.
  
  multiGroupVI does recover the signal, and contrastiveVI partly does. So the test can tell a working model
  from a failing one.
- **The paper's gene lists look like cell identity, not treatment.** The shared list matches border-associated
  macrophage markers. The doxorubicin list matches transcripts enriched in single-nucleus data.
- **On baselines,** MC-ContrastiveVI loses to PCA on reconstruction. It also carries less treatment information
  than scVI or multiGroupVI.
- **The archived model files behind many manuscript numbers are lost.** We refit them. Most numbers barely move.
  The cisplatin gene list and its ECM enrichment do **not** reproduce.
- **Our recommendation:** present the paper as a methods study, with the diagnostics as the contribution and
  MC-ContrastiveVI as the worked example. Drop the remaining biological claims. Details are in section 6.

---

## 2. What was done

| Run | When | Instructions | Outputs |
|---|---|---|---|
| Run 1 — reviewer analyses | 2026-10-03/04 | `CLAUDE_CODE_PROMPT.md` | `analysis/peerreview_20261003/`, `runs/peerreview_20261003/`, scripts 23–35 |
| Run A — follow-ups, consistent numbers, baseline positive control, archive | 2026-10-04 | `CLAUDE_CODE_PROMPT_RUN_A.md` | `analysis/peerreview_20261004/`, `runs/peerreview_20261004/`, scripts 36–47 |

Where to start reading:
- **`analysis/peerreview_20261004/INDEX.md`** describes every result file from both runs.
- Each run has a `SUMMARY.md` (deviations, run times, headline numbers) and a `PROGRESS.md` (step log).

The safeguards held in both runs:
- No existing file or committed result was modified.
- The raw data and the recovered cohort match their recorded SHA-256 hashes, so these are exactly the paper's
  data.
- Nothing was pushed.

---

## 3. Status of each reviewer comment

| Reviewer comment | Status | Where |
|---|---|---|
| Semi-synthetic positive control (null / shared / specific / mixed; effect sizes; responder fraction; F_sh vs truth; active units/KL; seed stability; attribution precision/recall; label exchange) | Done | `analysis/peerreview_20261003/D1_semisynthetic.md`; follow-up `analysis/peerreview_20261004/R1_injected_signal_location.md` |
| Optional parts of the positive control: A/B gene overlap, injection into one cell type | Not done, deliberately. GSE271055 labels are too unreliable (72.6% unassigned) for cell-type-restricted injection. | — |
| Compare attributions with cell-type marker lists | Done | `C1_marker_comparison.md` |
| Attribution with shuffled treatment labels | Done (see caveat in section 7) | `C2_shuffled_label_attribution.md` |
| State the cell-identity interpretation; soften "most stable gene-level signal" | **Writing — yours** | — |
| Within-study treatment decodability (Methods promised it, Results omit it) | Numbers available | `analysis/corrected_20260920/latent_probes.csv` (archived fits); `A1_latent_probes_fresh.csv` (fresh fits) |
| Cross-seed CKA for MC-ContrastiveVI as numbers | Numbers available | `cross_seed_latent_stability.csv` (archived); `A1_cross_seed_cka.csv` (fresh) |
| Baselines: contrastiveVI, multiGroupVI, scVI | Done | `C3_baselines.md`, `C4_multigroupvi.md`, combined table `R2_baseline_headline.md` |
| (Ours) Run the positive control on the baselines, to show the test can pass a working model | Done | `R4_baseline_positive_control.md` |

Paths without a folder are in `analysis/peerreview_20261003/` (C*, D1, A*) or `analysis/peerreview_20261004/` (R*).

---

## 4. Key results (facts)

**Positive control on MC-ContrastiveVI** (75 fits; `D1_semisynthetic.md`, `R1_*`)
- **F_sh does not track the true shared fraction:**

  | Scenario | True shared fraction | F_sh |
  |---|---|---|
  | Null (nothing injected) | — | 0.71–0.87 |
  | Specific only | 0 | 0.68–0.93 (configuration means) |
  | Shared only | 1 | 0.75–0.93 |
  | Mixed | 0.5 | 0.57–0.88 |

- **Active units in the treatment blocks** are about 0 in every configuration, including δ = 2 with every
  pseudo-treated cell responding.
- **The injected signal is not represented anywhere in the model.** Within-study probes on the background,
  ungated shared and ungated drug latents give 0.47–0.54 balanced accuracy in all configurations. The decoder
  output shows essentially none of the injected effect (≈ 0 ± 0.03 log units, against observed effects up to
  ±0.16), so the signal is ignored rather than absorbed by the background block.
- **Label exchange cannot detect a real injected shared program:** P = 0.14 (null: P = 0.77), from 100
  permutations each.
- **Integrated gradients don't recover the injected genes:** top-50 hits are at the level of random and null
  configurations.
- **The injected effects are realistic in size.** The median |log2FC| of real cisplatin pseudobulk hits is 0.74
  (0.98 within the 1,500-gene universe), and the grid spans δ = 0.25–2.

**Positive control on the baselines** (`R4_baseline_positive_control.md`). Configurations: δ ∈ {1, 2}, all
pseudo-treated cells responding. "Hits" are injected genes in the top 50, out of 80 injected genes per group;
random expectation 2.67.

| Model | Top-50 hits (null configuration) | Within-study probe | Active units |
|---|---|---|---|
| multiGroupVI | 3.7–11.0 (2.3–4.7) | 0.58–0.80 | 8–10 of 10 |
| contrastiveVI | cisplatin side 5–7 (4.3); doxorubicin side ≈ null | 0.52–0.73 | 8 of 8 |
| MC-ContrastiveVI | 1.0–4.0 (3.3–4.7) | 0.49–0.54 | 0–1.3 of 12 |

Even the best model recovers only a small share of injected genes in its top 50. Attribution-derived gene lists
are weak evidence in general.

**Gene lists vs. cell identity** (`C1`, `C5`, `C6`)
- **Paper's shared list:** 9 of 20 border-associated macrophage markers (BH q = 1.1e-9).
- **Paper's doxorubicin list:** 6 of 32 nucleus-enriched transcripts (q = 2.7e-5).
- **Paper's cisplatin list:** BAM 4/20 (q = 0.011); fibroblasts 5/52 (q = 0.035).
- Marker lists came only from published tables (PanglaoDB, Ochocka 2021, Zeisel 2018, Dani 2021, Bakken 2018).
  They were frozen and committed before any overlap was computed.
- **The ungated shared block predicts source cell type** in control cells at 0.72–0.82 balanced accuracy
  (chance 0.08).
- **Ttr is a choroid plexus marker:** mean log-expression 6.66 in choroid plexus vs 0.50 in all other cells.

**Shuffled labels** (`C2`)
- F_sh for shuffled-label fits: 0.12 / 0.87 / 0.98.
- The shuffled consensus lists barely overlap the real-label ones (Jaccard 0.015 / 0.19 / 0.13).
- The macrophage M8 signature is absent from the shuffled shared list.

**Baselines on the real data** (`R2_baseline_headline.md`)

| Model | Canonical-scale MSE | Within-study probe, ungated (cisplatin / doxorubicin) |
|---|---|---|
| PCA(32) | 0.155 | — |
| MC-ContrastiveVI | 0.169 | 0.53 / 0.55 |
| scVI | 0.206 | 0.81 / 0.63 |
| scDisInFact | 0.231 | — |
| multiGroupVI | 0.245 | 0.83 / 0.57 |

- contrastiveVI (fit per drug pair): cisplatin pair 0.126 (MC-ContrastiveVI on the same cells 0.092);
  doxorubicin pair 0.262 (MC-ContrastiveVI 0.203).
- multiGroupVI's *gated* probes (0.96 / 0.89) are inflated by the label masking. Don't use them. They are kept
  in a separate column only for transparency.
- **The MSE metric favors MC-ContrastiveVI and PCA,** which optimize it directly. The scvi-tools models are
  trained on counts.
- **Ablations (fresh refits):**
  - Removing gating or the contrastive structure barely changes MSE (0.1676–0.1678).
  - The same change drops F_sh from 92–97% to 18–33%. The high F_sh depends on the gating architecture.

---

## 5. Decisions you need to make

### Decision 1 — Archived or fresh numbers (most important)

The archived full-model checkpoints are lost. The one file you had (`full_0.pt`) turned out to be a different
fit (F_sh 0.754 vs 0.921). We refit seeds 0–2 with identical data and settings.

`analysis/peerreview_20261004/R3_paper_numbers.csv` lists **all 116 numbers in the manuscript**. For each, it
gives the manuscript value, the fresh value, whether it changed, and the source file:

| Category | Count |
|---|---|
| Changed | 28 |
| Unchanged | 77 |
| Not recomputable (figure-only or source metadata) | 11 |

`R3_consistent_numbers.md` has the compact table. The main changes:

| Quantity | Manuscript | Fresh |
|---|---|---|
| Pooled F_sh | 92.1–95.7% | 92.4–97.0% |
| Active shared coordinates | 0, 0, 1 | 0, 0, 0 |
| Shared top-50 cross-seed Jaccard | 0.190–0.639 | 0.449–0.562 |
| Consensus list sizes | 34/20/33 | 40/27/31 |
| Shared list (macrophage) | — | Reproduces (Jaccard 0.64). Macrophage signature stronger (q 7e-16). Reactome immune system now significant. |
| Cisplatin list | — | Does **not** reproduce (Jaccard 0.36). Top genes now include Mrc1, mt-Co2, Ttr, Apoe, Fos, Junb. ECM degradation no longer significant (q 0.32 vs 0.034). 0 Reactome terms. |
| Doxorubicin list | — | Changed (Jaccard 0.21). Now 2 Reactome terms at FDR 0.05 (manuscript: none). |
| Pseudobulk findings (cisplatin / pooled / same direction) | 297 / 354 / 353 | 293 / 350 / 349. Platform-level PyDESeq2 drift, unrelated to the models. |
| Rescue difference | 0.0003 | −0.0032 (P still 1.000) |

**Recommendation: use the fresh numbers throughout.**
- They are the only ones anyone, including a reviewer, can reproduce now.
- All the new analyses use them.
- They can be archived with a DOI.

State in Methods that the original checkpoints were not retained and the models were refit. The instability of
the cisplatin list across refits supports the paper's own message, so report it rather than hiding it.

### Decision 2 — What the paper claims

**Recommendation: frame it as a methods / reproducibility study.**
- The contribution is a set of diagnostics that catch a failure the headline statistic hides: active units and
  KL, a semi-synthetic positive control, label-exchange references, ungated probes, and comparison against
  models that do use their latents.
- MC-ContrastiveVI is the worked example.
- Presenting it as a new method that improves on existing ones isn't supportable: it loses to PCA on
  reconstruction and to scVI and multiGroupVI on treatment information.

### Decision 3 — Biological claims

**Recommendation: drop treatment-biology claims.**
- Report what the attributions track (cell identity, assay modality), citing C1, C5 and C2.
- The pseudobulk results can stay as exploratory, within-study (GSE216146) results, with the three-pair
  limitation.
- The cross-drug question isn't answerable with these two studies. Drug, study and assay are fully confounded.
  That was already in the limitations.

### Decision 4 — Further work (optional)

Not done:
- **Attempt to fix MC-ContrastiveVI** (free-bits KL, a contrastiveVI-style background constraint, a count
  likelihood) and show the fix passes the positive control. This would turn the paper into a positive methods
  result. It needs a pre-registered plan and roughly a day of compute. Manan can run it if you want it.
- **A new dataset with several drugs under one protocol.** Only needed for biological conclusions. Not needed for
  publication.

### Decision 5 — Release

- Merge the branch, or keep it as a separate PR.
- Who uploads the archive to Zenodo (needs an account; gives a DOI), and when. Upload once the numbers in the
  manuscript are final.
- The administrative items your own `FINAL_SUBMISSION_READINESS_AUDIT.md` lists: corresponding author,
  contributions, funding, competing interests, ethics wording for public animal data.

---

## 6. What the results support adding to the manuscript

Suggestions only; the wording is yours.

**Results**
- **The positive control as a central result:** F_sh vs. truth, collapse even at a strong effect, label exchange
  P = 0.14 with a real injected program, and attribution recovery.
- **The baseline positive control:** multiGroupVI and contrastiveVI recover the signal, so the test can tell a
  working model from a failing one. This is the answer to "this is just one bad model."
- **The missing numbers the reviewer flagged:** within-study treatment probes and cross-seed CKA.
- **The marker comparison, the shuffled-label test and the control-cell identity probe,** supporting the
  cell-identity reading.
- **The ablation F_sh values** (18–33% without gating), showing F_sh depends on the architecture.

**Discussion**
- Replace "most stable gene-level signal" (immune/myeloid) with the cell-identity interpretation.
- Say why Ttr appears in two lists (choroid plexus marker, likely ambient RNA).

**Methods**
- Refit and reproducibility notes: Windows CPU; pins in `pip_freeze_venv.txt`; the multiGroupVI environment.
- Semi-synthetic design: fixed in advance and listed in `D1_design.json` and `D1_semisynthetic.md`.

**Data and code availability**
- The branch, plus the Zenodo DOI for the archive.

Data for figures (no figures were made):
- Figure 5 data: `R3_consensus_top12_figure5_data.csv`.
- Positive-control tables: `D1_config_summary.csv`, `R1_summary.csv`, `R4_comparison_table.csv`.
- Baseline table: `R2_baseline_headline.csv`.
- Figure 4 data: `R3_09_interpretation/` and `R3_11_diagnostics/`.

---

## 7. Caveats to disclose or keep in mind

- **The shuffled-label test (C2)** shuffles at the cell level, which breaks the link between label and library.
  It shows the attributions depend on the labels. It can't separate a treatment effect from library-level
  composition differences.
- **contrastiveVI handles two groups at a time,** so it was fit per drug pair. It is not a multi-drug
  shared/specific model.
- **multiGroupVI was run without a study covariate** (authors' setup). Its *gated* probes are label-masked; use
  the ungated ones.
- **Integrated gradients are taken in each model's own input space** (counts for multiGroupVI, log1p for
  contrastiveVI, log-normalized for MC-ContrastiveVI).
- **The original fits ran on Linux; the refits ran on Windows (CPU).** Small numerical drift is expected.
  Highly variable gene reselection gave 1,499 of the 1,500 genes, so the committed gene universe was used
  throughout.
- **The semi-synthetic responder-fraction-0.3 runs** use the same responder cells for all injected sets.

---

## 8. How to reproduce

`analysis/peerreview_20261004/REPRODUCE.md` lists every command in order:
1. Environments (pinned).
2. Data download and recovery, with hash checks.
3. Both runs' scripts.

The archive is 7 zips, 1.2 GB, listed in `MANIFEST.csv` with SHA-256 checksums. It holds all fitted
checkpoints and latents, the canonical inputs, the semi-synthetic data and the baseline outputs. Raw GEO data is
excluded (public; hashes in `EXCLUDED_INPUT_HASHES.csv`).
