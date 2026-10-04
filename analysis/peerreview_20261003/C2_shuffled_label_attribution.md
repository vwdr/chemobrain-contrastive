# C2 — Shuffled-label attribution

Script: `scripts/25_shuffled_label_attribution.py` (`--fit --seed S --threads=3`; `--analyze`). Fits:
`runs/peerreview_20261003/shuffled_labels/`. Logs: `runs/peerreview_20261003/logs/25_*.log`.

## Design

- Cells, genes, split and model: identical to Phase 2 (canonical inputs; 6,000 training cells; full model,
  canonical configuration, seeds 0, 1, 2; 3 CPU threads per fit).
- Labels: treatment labels permuted at the cell level within each study (GSE216146 control/cisplatin; GSE271055
  control/doxorubicin) with `numpy.random.default_rng(20261003)`, applied in turn to the training, validation and test
  cells of each study (so early stopping, gating, F_sh and attribution-cell selection all use the shuffled design).
  Class counts per study × split are unchanged (`C2_label_shuffle_design.csv`):

| split | study | n | treated (real = shuffled) | fraction of labels changed | shuffled-treated that are real-treated |
|---|---|---|---|---|---|
| train | GSE216146 | 3000 | 1500 | 0.501 | 749 |
| train | GSE271055 | 3000 | 1500 | 0.499 | 751 |
| validation | GSE216146 | 884 | 401 | 0.471 | 193 |
| validation | GSE271055 | 1997 | 1052 | 0.516 | 537 |
| test | GSE216146 | 883 | 400 | 0.503 | 178 |
| test | GSE271055 | 1998 | 1052 | 0.488 | 564 |

- Attribution and consensus: identical to Phase 2 step 4 (squared posterior-mean norm; zero and control-median
  baselines; `default_rng(314)` attribution cells drawn from shuffled-label test cells; top 50 in ≥ 4 of 6).
- Enrichment: `enrichment()` and `LIBRARIES` imported from `scripts/19_gene_pattern_audit.py` (1,500-gene universe;
  set size 10–500 after intersection; one-sided hypergeometric; BH within list × block × library).

## Fit metrics (`C2_shuffled_fit_metrics.csv`; F_sh and diagnostics on shuffled-label test cells)

| seed | best epoch | val MSE | test MSE | F_sh pooled | F_sh dox | F_sh cis | active units sh/dox/cis | dims KL>0.01 sh/dox/cis | total KL sh/dox/cis (nats) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 98 | 0.165987 | 0.168045 | 0.1245 | 0.3177 | 0.0618 | 0/0/3 | 0/0/3 | 0.00314/0.00372/0.03775 |
| 1 | 99 | 0.166636 | 0.168701 | 0.8674 | 0.5300 | 0.9868 | 0/0/0 | 0/0/0 | 0.01525/0.00317/0.00215 |
| 2 | 98 | 0.166282 | 0.168329 | 0.9786 | 0.8208 | 0.9946 | 1/0/0 | 0/0/0 | 0.02525/0.00092/0.00100 |

Fit time 189 s per fit; IG maximum median scaled completeness error 0.0088 (256 quadrature points).

## Shuffled consensus lists (`C2_shuffled_consensus_genes.csv`)

- shared (27): Ccdc153, Sncg, Dmkn, Enkur, S100a10, Anxa1, Gm19935, Odf3b, Crip2, Ttr, Tcp11, Vim, Cdhr4, Dnah10,
  Crocc2, Igfbp2, Lrrc36, Hydin, Frmpd2, Dnah3, Dnah11, Cfap46, Ak9, Igfbp6, Camk2a, Mbp, Cfap97d2
- doxorubicin (30): Gria2, Ttr, Snhg11, Meg3, Plp1, Atp1a2, Grm5, Nrxn1, Miat, Neat1, Syt1, Apoe, Map1b, Mbp, Brd9,
  Grin2a, Zeb2, Peg3, Rnf112, Dclk1, Gjc3, Gatm, Chn1, Camk2a, Ndrg4, C1ql3, Exoc6b, Gria1, Enpp2, Ryr2
- cisplatin (20): Ccdc153, mt-Co2, Tmsb10, Nnat, Ptgds, Mt3, Bc1, Ly6h, Plp1, mt-Co3, Sat1, Qk, Bsg, Apoe, Zbtb20,
  C1qc, Pltp, Meg3, Ly6e, Snhg11

## Overlap with real-label lists (`C2_overlap.csv`; hypergeometric, N = 1,500)

| comparison | block | n shuffled | n reference | overlap | Jaccard | hypergeom P | overlap genes |
|---|---|---|---|---|---|---|---|
| (a) vs fresh real-label consensus | shared | 27 | 40 | 1 | 0.015 | 0.521 | Ttr |
| (a) | doxorubicin | 30 | 27 | 9 | 0.188 | 5.1e-10 | Apoe, Camk2a, Grm5, Meg3, Nrxn1, Peg3, Snhg11, Ttr, Zeb2 |
| (a) | cisplatin | 20 | 31 | 6 | 0.133 | 1.5e-06 | Apoe, Bsg, Mt3, Qk, Zbtb20, mt-Co2 |
| (b) vs frozen paper lists | shared | 27 | 34 | 1 | 0.017 | 0.464 | Ttr |
| (b) | doxorubicin | 30 | 20 | 9 | 0.220 | 2.0e-11 | C1ql3, Dclk1, Gria2, Meg3, Miat, Peg3, Snhg11, Ttr, Zeb2 |
| (b) | cisplatin | 20 | 33 | 3 | 0.060 | 8.6e-03 | Apoe, Bsg, Tmsb10 |

Cross-block overlaps (shuffled block vs a different fresh real-label block) are also in `C2_overlap.csv`
(e.g. shuffled cisplatin vs fresh shared: 5 genes, P = 1.2e-4; shuffled cisplatin vs fresh doxorubicin: 4, P = 3.3e-4).

## Enrichment (`C2_enrichment_summary.csv`, `C2_top_enrichment_terms.csv`, `C2_enrichment_all.csv.gz`)

Terms with BH q < 0.05 (terms tested: GO_BP 1,625; M8 164; Reactome 114):

| list | block | GO_BP | M8 | Reactome |
|---|---|---|---|---|
| shuffled | shared | 0 | 1 | 0 |
| shuffled | doxorubicin | 31 | 3 | 4 |
| shuffled | cisplatin | 0 | 56 | 0 |
| frozen paper | shared | 11 | 46 | 0 |
| frozen paper | doxorubicin | 0 | 0 | 0 |
| frozen paper | cisplatin | 0 | 77 | 5 |
| fresh real | shared | 6 | 73 | 6 |
| fresh real | doxorubicin | 0 | 0 | 2 |
| fresh real | cisplatin | 0 | 56 | 0 |

Top terms of the shuffled lists: shared — M8 DESCARTES_ORGANOGENESIS_EPENDYMAL_CELL (10/72, q = 2.5e-5);
doxorubicin — GO_BP GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY (7/17, q = 1.7e-5), Reactome NEURONAL_SYSTEM (7/58,
q = 0.010), M8 TABULA_MURIS_SENIS_BRAIN_NON_MYELOID_NEURON_AGEING (8/93, q = 0.030); cisplatin — M8
TABULA_MURIS_SENIS_BRAIN_NON_MYELOID_OLIGODENDROCYTE_PRECURSOR_CELL_AGEING (7/47, q = 1.8e-4).

### Macrophage M8 signature

The top-ranked M8 term for the paper's frozen shared list (recomputed here with the 19_gene_pattern_audit procedure)
is `ZHANG_UTERUS_C5_MACROPHAGE` (58 genes in universe). Results for that term (`C2_paper_top_m8_term_lookup.csv`):

| list | block | overlap | fold enrichment | P | BH q |
|---|---|---|---|---|---|
| frozen paper | shared | 14 | 10.65 | 2.6e-12 | 4.2e-10 |
| fresh real | shared | 19 | 12.28 | 4.5e-18 | 7.3e-16 |
| **shuffled** | **shared** | **0** | 0 | 1.0 | 1.0 — not significant |
| shuffled | cisplatin | 5 | 6.47 | 7.2e-4 | 7.9e-3 |
| shuffled | doxorubicin | 1 | 0.86 | 0.70 | 1.0 |

## Deviation

The shuffle was applied to validation and test cells as well as training cells (the specification names the training
cells; the shuffled design was carried through so that early stopping, gating and attribution-cell selection use one
consistent label assignment). C2 fits were run before C1/C5 (background CPU use); analysis followed C1/C5.
