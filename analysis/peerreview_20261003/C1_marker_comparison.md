# C1 â€” Attribution consensus lists vs published cell-type marker lists

Scripts: `scripts/26_c1_marker_list.py` (list construction; list committed in 3690ba7 before any overlap was
computed), `scripts/27_c1_marker_comparison.py` (overlap). Log: `runs/peerreview_20261003/logs/27_c1.log`.

## Marker list (`C1_marker_list.csv`, sources and SHA-256 in `C1_marker_sources.json`)

Retrieved 2026-10-03 from the internet; no gene was added from memory. Symbols are stored as given by each source
and matched case-insensitively to the 1,500-gene universe.

| cell type | source | genes listed | in universe |
|---|---|---|---|
| Microglia, Macrophages, Astrocytes, Oligodendrocytes, OPCs (PanglaoDB "Oligodendrocyte progenitor cells"), Neurons, Endothelial cells, Pericytes, Fibroblasts, Ependymal cells, Choroid plexus cells | PanglaoDB markers file 27 Mar 2020 (FranzÃ©n et al. 2019, doi:10.1093/database/baz046), rows whose species field contains "Mm" | 79, 147, 63, 87, 28, 211, 188, 64, 177, 59, 22 | 24, 45, 35, 24, 21, 47, 78, 32, 52, 35, 8 |
| Border-associated macrophages | Ochocka et al. 2021 Nat Commun (doi:10.1038/s41467-021-21407-w) Suppl. Data 1: BAM-cluster top-30 DEGs present in both the "female control" and "male control" sheets; plus Zeisel et al. 2018 Cell (doi:10.1016/j.cell.2018.06.021) Table S4 marker genes of PVM1 and PVM2 | 25 (union) | 20 (Ochocka 20; Zeisel 5, all contained in the Ochocka set) |
| Meningeal/perivascular fibroblasts | Zeisel et al. 2018 Table S4 marker genes of VLMC1, VLMC2, ABC ("Vascular leptomeningeal cells") | 15 | 6 |
| Choroid plexus epithelium | Dani et al. 2021 Cell (doi:10.1016/j.cell.2021.04.003) Table S1 sheet "Cell Type" column "Epithelial"; plus Zeisel et al. 2018 Table S4 CHOR | 88 (union) | 37 (Dani 36; Zeisel 4) |
| Nucleus-enriched transcripts | Bakken et al. 2018 PLoS One (doi:10.1371/journal.pone.0209648) S2 Table: logFC < âˆ’log2(1.5) and adj.P.Val < 0.05 (159 genes, matching the count stated in the paper) | 159 | 32 |

Sources inspected and not used (no suitable marker table): Van Hove et al. 2019 Nat Neurosci Suppl. Tables 2â€“3;
PietilÃ¤ et al. 2023 Neuron mmc2â€“mmc4; DeSisto et al. 2020 Dev Cell mmc2â€“mmc6 (meningeal sub-cluster contrasts only);
Vanlandewijck et al. 2018 Nature Suppl. Table 3 and Fig. 1 source data. The PMC copy of Dani Table S1 was behind a
proof-of-work bot challenge (not bypassed); the publisher-hosted copy (mmc1.xlsx) was used.

## Test

Universe = the 1,500 canonical genes. For each consensus list and each marker set (union over sources):
overlap k, one-sided hypergeometric P(X â‰¥ k) with N = 1,500, K = markers in universe, n = list size; BH q across the
15 marker sets within each list. Per-source results for the two multi-source sets are in `C1_marker_overlap.csv`
(rows with `marker_source` â‰  `all_sources`). Primary lists: frozen lists in `scripts/19_gene_pattern_audit.py`
(34 shared / 20 doxorubicin / 33 cisplatin); secondary: fresh Phase 2 consensus (40 / 27 / 31).

## Results (`C1_marker_overlap.csv`)

### frozen_paper consensus lists

| block | cell type | n list | n markers (in universe) | overlap | expected | hypergeom P | BH q | overlap genes |
|---|---|---|---|---|---|---|---|---|
| shared | Astrocytes | 34 | 35 | 1 | 0.79 | 0.556 | 1 | Apoe |
| shared | Border-associated macrophages | 34 | 20 | 9 | 0.45 | 7.2e-11 | 1.08e-09 | Apoe|Dab2|F13a1|Ifi27l2a|Ifitm3|Lgals1|Mrc1|Ms4a7|Pf4 |
| shared | Choroid plexus cells | 34 | 8 | 1 | 0.18 | 0.168 | 0.611 | Ttr |
| shared | Choroid plexus epithelium | 34 | 37 | 2 | 0.84 | 0.204 | 0.611 | Atp5g1|Ttr |
| shared | Endothelial cells | 34 | 78 | 1 | 1.77 | 0.841 | 1 | Stab1 |
| shared | Ependymal cells | 34 | 35 | 0 | 0.79 | 1 | 1 |  |
| shared | Fibroblasts | 34 | 52 | 0 | 1.18 | 1 | 1 |  |
| shared | Macrophages | 34 | 45 | 7 | 1.02 | 3.97e-05 | 0.000298 | C5ar1|Cd163|Dab2|F13a1|Mrc1|Ms4a7|Stab1 |
| shared | Meningeal/perivascular fibroblasts | 34 | 6 | 0 | 0.14 | 1 | 1 |  |
| shared | Microglia | 34 | 24 | 4 | 0.54 | 0.0017 | 0.00849 | C5ar1|Cbr2|Ccl4|Pf4 |
| shared | Neurons | 34 | 47 | 0 | 1.07 | 1 | 1 |  |
| shared | Nucleus-enriched transcripts | 34 | 32 | 0 | 0.73 | 1 | 1 |  |
| shared | OPCs | 34 | 21 | 0 | 0.48 | 1 | 1 |  |
| shared | Oligodendrocytes | 34 | 24 | 0 | 0.54 | 1 | 1 |  |
| shared | Pericytes | 34 | 32 | 0 | 0.73 | 1 | 1 |  |
| doxorubicin | Astrocytes | 20 | 35 | 0 | 0.47 | 1 | 1 |  |
| doxorubicin | Border-associated macrophages | 20 | 20 | 0 | 0.27 | 1 | 1 |  |
| doxorubicin | Choroid plexus cells | 20 | 8 | 1 | 0.11 | 0.102 | 0.51 | Ttr |
| doxorubicin | Choroid plexus epithelium | 20 | 37 | 1 | 0.49 | 0.395 | 1 | Ttr |
| doxorubicin | Endothelial cells | 20 | 78 | 0 | 1.04 | 1 | 1 |  |
| doxorubicin | Ependymal cells | 20 | 35 | 0 | 0.47 | 1 | 1 |  |
| doxorubicin | Fibroblasts | 20 | 52 | 1 | 0.69 | 0.508 | 1 | Zeb2 |
| doxorubicin | Macrophages | 20 | 45 | 0 | 0.60 | 1 | 1 |  |
| doxorubicin | Meningeal/perivascular fibroblasts | 20 | 6 | 0 | 0.08 | 1 | 1 |  |
| doxorubicin | Microglia | 20 | 24 | 0 | 0.32 | 1 | 1 |  |
| doxorubicin | Neurons | 20 | 47 | 4 | 0.63 | 0.00284 | 0.0213 | Meg3|Miat|Pcsk2|Snhg11 |
| doxorubicin | Nucleus-enriched transcripts | 20 | 32 | 6 | 0.43 | 1.82e-06 | 2.72e-05 | Kcnq1ot1|Meg3|Miat|Pcdh9|Ptprd|Snhg11 |
| doxorubicin | OPCs | 20 | 21 | 0 | 0.28 | 1 | 1 |  |
| doxorubicin | Oligodendrocytes | 20 | 24 | 0 | 0.32 | 1 | 1 |  |
| doxorubicin | Pericytes | 20 | 32 | 0 | 0.43 | 1 | 1 |  |
| cisplatin | Astrocytes | 33 | 35 | 1 | 0.77 | 0.545 | 0.818 | Apoe |
| cisplatin | Border-associated macrophages | 33 | 20 | 4 | 0.44 | 0.000735 | 0.011 | Apoe|F13a1|Mrc1|Pf4 |
| cisplatin | Choroid plexus cells | 33 | 8 | 1 | 0.18 | 0.163 | 0.408 | Ttr |
| cisplatin | Choroid plexus epithelium | 33 | 37 | 2 | 0.81 | 0.195 | 0.417 | Bsg|Ttr |
| cisplatin | Endothelial cells | 33 | 78 | 3 | 1.72 | 0.244 | 0.457 | Id3|Mgp|S100a10 |
| cisplatin | Ependymal cells | 33 | 35 | 0 | 0.77 | 1 | 1 |  |
| cisplatin | Fibroblasts | 33 | 52 | 5 | 1.14 | 0.00468 | 0.0351 | Col1a2|Col6a2|Efemp1|Igfbp6|Mgp |
| cisplatin | Macrophages | 33 | 45 | 3 | 0.99 | 0.0735 | 0.276 | Cd163|F13a1|Mrc1 |
| cisplatin | Meningeal/perivascular fibroblasts | 33 | 6 | 1 | 0.13 | 0.125 | 0.375 | Mgp |
| cisplatin | Microglia | 33 | 24 | 3 | 0.53 | 0.0143 | 0.0717 | Cbr2|Fos|Pf4 |
| cisplatin | Neurons | 33 | 47 | 0 | 1.03 | 1 | 1 |  |
| cisplatin | Nucleus-enriched transcripts | 33 | 32 | 0 | 0.70 | 1 | 1 |  |
| cisplatin | OPCs | 33 | 21 | 0 | 0.46 | 1 | 1 |  |
| cisplatin | Oligodendrocytes | 33 | 24 | 1 | 0.53 | 0.416 | 0.694 | Mbp |
| cisplatin | Pericytes | 33 | 32 | 0 | 0.70 | 1 | 1 |  |

### fresh_phase2 consensus lists

| block | cell type | n list | n markers (in universe) | overlap | expected | hypergeom P | BH q | overlap genes |
|---|---|---|---|---|---|---|---|---|
| shared | Astrocytes | 40 | 35 | 1 | 0.93 | 0.616 | 1 | Apoe |
| shared | Border-associated macrophages | 40 | 20 | 10 | 0.53 | 8.44e-12 | 1.27e-10 | Apoe|Dab2|F13a1|Ifi27l2a|Ifitm3|Lgals1|Lyz2|Mrc1|Ms4a7|Pf4 |
| shared | Choroid plexus cells | 40 | 8 | 1 | 0.21 | 0.195 | 0.731 | Ttr |
| shared | Choroid plexus epithelium | 40 | 37 | 2 | 0.99 | 0.259 | 0.777 | Atp5g1|Ttr |
| shared | Endothelial cells | 40 | 78 | 2 | 2.08 | 0.626 | 1 | Sparcl1|Stab1 |
| shared | Ependymal cells | 40 | 35 | 0 | 0.93 | 1 | 1 |  |
| shared | Fibroblasts | 40 | 52 | 0 | 1.39 | 1 | 1 |  |
| shared | Macrophages | 40 | 45 | 9 | 1.20 | 1.18e-06 | 8.86e-06 | Cd163|Dab2|F13a1|Lyz2|Maf|Mrc1|Ms4a7|Stab1|Tyrobp |
| shared | Meningeal/perivascular fibroblasts | 40 | 6 | 0 | 0.16 | 1 | 1 |  |
| shared | Microglia | 40 | 24 | 3 | 0.64 | 0.0241 | 0.121 | Cbr2|Ctss|Pf4 |
| shared | Neurons | 40 | 47 | 0 | 1.25 | 1 | 1 |  |
| shared | Nucleus-enriched transcripts | 40 | 32 | 0 | 0.85 | 1 | 1 |  |
| shared | OPCs | 40 | 21 | 0 | 0.56 | 1 | 1 |  |
| shared | Oligodendrocytes | 40 | 24 | 0 | 0.64 | 1 | 1 |  |
| shared | Pericytes | 40 | 32 | 0 | 0.85 | 1 | 1 |  |
| doxorubicin | Astrocytes | 27 | 35 | 2 | 0.63 | 0.129 | 0.424 | Apoe|Slc1a2 |
| doxorubicin | Border-associated macrophages | 27 | 20 | 1 | 0.36 | 0.306 | 0.766 | Apoe |
| doxorubicin | Choroid plexus cells | 27 | 8 | 1 | 0.14 | 0.136 | 0.424 | Ttr |
| doxorubicin | Choroid plexus epithelium | 27 | 37 | 2 | 0.67 | 0.141 | 0.424 | Atp1b1|Ttr |
| doxorubicin | Endothelial cells | 27 | 78 | 1 | 1.40 | 0.767 | 1 | Sparcl1 |
| doxorubicin | Ependymal cells | 27 | 35 | 0 | 0.63 | 1 | 1 |  |
| doxorubicin | Fibroblasts | 27 | 52 | 1 | 0.94 | 0.618 | 1 | Zeb2 |
| doxorubicin | Macrophages | 27 | 45 | 0 | 0.81 | 1 | 1 |  |
| doxorubicin | Meningeal/perivascular fibroblasts | 27 | 6 | 0 | 0.11 | 1 | 1 |  |
| doxorubicin | Microglia | 27 | 24 | 0 | 0.43 | 1 | 1 |  |
| doxorubicin | Neurons | 27 | 47 | 8 | 0.85 | 7.15e-07 | 1.07e-05 | Dlgap1|Grm5|Meg3|Nrgn|Pclo|Slc17a7|Snap25|Snhg11 |
| doxorubicin | Nucleus-enriched transcripts | 27 | 32 | 3 | 0.58 | 0.0182 | 0.137 | Kcnq1ot1|Meg3|Snhg11 |
| doxorubicin | OPCs | 27 | 21 | 0 | 0.38 | 1 | 1 |  |
| doxorubicin | Oligodendrocytes | 27 | 24 | 0 | 0.43 | 1 | 1 |  |
| doxorubicin | Pericytes | 27 | 32 | 0 | 0.58 | 1 | 1 |  |
| cisplatin | Astrocytes | 31 | 35 | 1 | 0.72 | 0.523 | 1 | Apoe |
| cisplatin | Border-associated macrophages | 31 | 20 | 2 | 0.41 | 0.0623 | 0.312 | Apoe|Mrc1 |
| cisplatin | Choroid plexus cells | 31 | 8 | 1 | 0.17 | 0.154 | 0.463 | Ttr |
| cisplatin | Choroid plexus epithelium | 31 | 37 | 4 | 0.76 | 0.00613 | 0.0919 | Bsg|Chchd10|Mt3|Ttr |
| cisplatin | Endothelial cells | 31 | 78 | 1 | 1.61 | 0.812 | 1 | Mgp |
| cisplatin | Ependymal cells | 31 | 35 | 0 | 0.72 | 1 | 1 |  |
| cisplatin | Fibroblasts | 31 | 52 | 4 | 1.07 | 0.0202 | 0.152 | Col6a2|Igfbp6|Klf4|Mgp |
| cisplatin | Macrophages | 31 | 45 | 1 | 0.93 | 0.615 | 1 | Mrc1 |
| cisplatin | Meningeal/perivascular fibroblasts | 31 | 6 | 1 | 0.12 | 0.118 | 0.442 | Mgp |
| cisplatin | Microglia | 31 | 24 | 1 | 0.50 | 0.397 | 0.991 | Fos |
| cisplatin | Neurons | 31 | 47 | 0 | 0.97 | 1 | 1 |  |
| cisplatin | Nucleus-enriched transcripts | 31 | 32 | 0 | 0.66 | 1 | 1 |  |
| cisplatin | OPCs | 31 | 21 | 0 | 0.43 | 1 | 1 |  |
| cisplatin | Oligodendrocytes | 31 | 24 | 0 | 0.50 | 1 | 1 |  |
| cisplatin | Pericytes | 31 | 32 | 0 | 0.66 | 1 | 1 |  |

## Per-gene annotation, frozen lists (`C1_consensus_gene_annotation.csv` also covers the fresh lists)

| block | rank | gene | marker sets |
|---|---|---|---|
| shared | 1 | Mrc1 | Border-associated macrophages|Macrophages |
| shared | 2 | Pf4 | Border-associated macrophages|Microglia |
| shared | 3 | Ms4a7 | Border-associated macrophages|Macrophages |
| shared | 4 | Cbr2 | Microglia |
| shared | 5 | Cd163 | Macrophages |
| shared | 6 | F13a1 | Border-associated macrophages|Macrophages |
| shared | 7 | Gpx3 | — |
| shared | 8 | Dab2 | Border-associated macrophages|Macrophages |
| shared | 9 | mt-Co3 | — |
| shared | 10 | Cd52 | — |
| shared | 11 | Stab1 | Endothelial cells|Macrophages |
| shared | 12 | C5ar1 | Macrophages|Microglia |
| shared | 13 | mt-Nd4 | — |
| shared | 14 | mt-Co2 | — |
| shared | 15 | Apoe | Astrocytes|Border-associated macrophages |
| shared | 16 | Ttr | Choroid plexus cells|Choroid plexus epithelium |
| shared | 17 | Laptm5 | — |
| shared | 18 | Lilrb4a | — |
| shared | 19 | Ifi27l2a | Border-associated macrophages |
| shared | 20 | Bst2 | — |
| shared | 21 | Atp5g1 | Choroid plexus epithelium |
| shared | 22 | Psmb8 | — |
| shared | 23 | Cpne2 | — |
| shared | 24 | Ifitm3 | Border-associated macrophages |
| shared | 25 | Lgals1 | Border-associated macrophages |
| shared | 26 | Ctsd | — |
| shared | 27 | Tec | — |
| shared | 28 | Junb | — |
| shared | 29 | Cxcl12 | — |
| shared | 30 | Slc6a6 | — |
| shared | 31 | Ptprc | — |
| shared | 32 | Ccl4 | Microglia |
| shared | 33 | C1qa | — |
| shared | 34 | Qk | — |
| doxorubicin | 1 | Snhg11 | Neurons|Nucleus-enriched transcripts |
| doxorubicin | 2 | Kcnq1ot1 | Nucleus-enriched transcripts |
| doxorubicin | 3 | Meg3 | Neurons|Nucleus-enriched transcripts |
| doxorubicin | 4 | Ttr | Choroid plexus cells|Choroid plexus epithelium |
| doxorubicin | 5 | Gria2 | — |
| doxorubicin | 6 | Rtn1 | — |
| doxorubicin | 7 | Grin2b | — |
| doxorubicin | 8 | C1ql3 | — |
| doxorubicin | 9 | Opcml | — |
| doxorubicin | 10 | Pcdh9 | Nucleus-enriched transcripts |
| doxorubicin | 11 | Nrcam | — |
| doxorubicin | 12 | Dclk1 | — |
| doxorubicin | 13 | Peg3 | — |
| doxorubicin | 14 | Mir100hg | — |
| doxorubicin | 15 | Zeb2 | Fibroblasts |
| doxorubicin | 16 | Pcsk2 | Neurons |
| doxorubicin | 17 | Ptprd | Nucleus-enriched transcripts |
| doxorubicin | 18 | Zbtb20 | — |
| doxorubicin | 19 | Nfib | — |
| doxorubicin | 20 | Miat | Neurons|Nucleus-enriched transcripts |
| cisplatin | 1 | Ccn3 | — |
| cisplatin | 2 | Col25a1 | — |
| cisplatin | 3 | Igfbp6 | Fibroblasts |
| cisplatin | 4 | Mgp | Endothelial cells|Fibroblasts|Meningeal/perivascular fibroblasts |
| cisplatin | 5 | Foxc2 | — |
| cisplatin | 6 | Rspo3 | — |
| cisplatin | 7 | Col6a2 | Fibroblasts |
| cisplatin | 8 | mt-Nd4 | — |
| cisplatin | 9 | Ogn | — |
| cisplatin | 10 | Col6a1 | — |
| cisplatin | 11 | Adamtsl3 | — |
| cisplatin | 12 | Foxd1 | — |
| cisplatin | 13 | Efemp1 | Fibroblasts |
| cisplatin | 14 | S100a10 | Endothelial cells |
| cisplatin | 15 | Bsg | Choroid plexus epithelium |
| cisplatin | 16 | Hes1 | — |
| cisplatin | 17 | Fos | Microglia |
| cisplatin | 18 | Ttr | Choroid plexus cells|Choroid plexus epithelium |
| cisplatin | 19 | Ildr2 | — |
| cisplatin | 20 | Mt1 | — |
| cisplatin | 21 | Col1a2 | Fibroblasts |
| cisplatin | 22 | Tmsb10 | — |
| cisplatin | 23 | Apoe | Astrocytes|Border-associated macrophages |
| cisplatin | 24 | Tsc22d1 | — |
| cisplatin | 25 | Anxa1 | — |
| cisplatin | 26 | Id3 | Endothelial cells |
| cisplatin | 27 | Mbp | Oligodendrocytes |
| cisplatin | 28 | Cbr2 | Microglia |
| cisplatin | 29 | Folr2 | — |
| cisplatin | 30 | Mrc1 | Border-associated macrophages|Macrophages |
| cisplatin | 31 | Pf4 | Border-associated macrophages|Microglia |
| cisplatin | 32 | F13a1 | Border-associated macrophages|Macrophages |
| cisplatin | 33 | Cd163 | Macrophages |
