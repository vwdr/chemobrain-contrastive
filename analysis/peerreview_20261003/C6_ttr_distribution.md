# C6 — Ttr distribution in GSE216146

Script: `scripts/32_c6_ttr_distribution.py` (log `runs/peerreview_20261003/logs/32_c6.log`).
Data: `data/processed/recovered_counts.h5ad`, all GSE216146 cells after the canonical QC (14,965 cells).
Arms: `control` = PBS (source code PN), `control_GENUS` = PS, `cisplatin` = CN, `cisplatin_GENUS` = CS
(the GENUS arms are rescue arms, excluded from model training). Log-normalized = log1p(count / total counts over all
18,271 genes × 1e4), as in the canonical pipeline.

## Choroid Plexus vs all other cell types (`C6_ttr_choroid_plexus_vs_other.csv`)

| group | arm | n cells | fraction nonzero | mean raw counts | mean log-normalized |
|---|---|---|---|---|---|
| All other cell types | cisplatin | 2841 | 0.907 | 8.768 | 1.905 |
| All other cell types | cisplatin_GENUS | 1557 | 0.844 | 2.617 | 1.248 |
| All other cell types | control | 4577 | 0.428 | 0.988 | 0.499 |
| All other cell types | control_GENUS | 3362 | 0.617 | 1.526 | 0.949 |
| Choroid Plexus | cisplatin | 1164 | 1.000 | 1278.875 | 7.290 |
| Choroid Plexus | cisplatin_GENUS | 475 | 1.000 | 809.025 | 6.891 |
| Choroid Plexus | control | 253 | 1.000 | 632.581 | 6.655 |
| Choroid Plexus | control_GENUS | 736 | 1.000 | 543.758 | 6.559 |
| All other cell types | all_arms | 12337 | 0.642 | 3.132 | 1.040 |
| Choroid Plexus | all_arms | 2628 | 1.000 | 925.855 | 6.952 |

## By source cell type and arm (`C6_ttr_by_celltype_arm.csv`)

| cell type | arm | n cells | fraction nonzero | mean raw counts | mean log-normalized |
|---|---|---|---|---|---|
| Astrocyte | cisplatin | 296 | 0.916 | 8.287 | 1.893 |
| Astrocyte | cisplatin_GENUS | 204 | 0.848 | 2.608 | 1.288 |
| Astrocyte | control | 564 | 0.429 | 0.612 | 0.519 |
| Astrocyte | control_GENUS | 373 | 0.601 | 1.220 | 0.916 |
| Choroid Plexus | cisplatin | 1164 | 1.000 | 1278.875 | 7.290 |
| Choroid Plexus | cisplatin_GENUS | 475 | 1.000 | 809.025 | 6.891 |
| Choroid Plexus | control | 253 | 1.000 | 632.581 | 6.655 |
| Choroid Plexus | control_GENUS | 736 | 1.000 | 543.758 | 6.559 |
| Endothelial | cisplatin | 398 | 0.859 | 7.063 | 1.725 |
| Endothelial | cisplatin_GENUS | 207 | 0.778 | 2.227 | 1.156 |
| Endothelial | control | 399 | 0.388 | 0.594 | 0.540 |
| Endothelial | control_GENUS | 286 | 0.566 | 1.164 | 1.002 |
| Ependymal | cisplatin | 36 | 0.861 | 10.833 | 1.676 |
| Ependymal | cisplatin_GENUS | 11 | 0.909 | 2.545 | 1.299 |
| Ependymal | control | 49 | 0.408 | 0.510 | 0.382 |
| Ependymal | control_GENUS | 38 | 0.605 | 1.105 | 0.826 |
| Fibroblast | cisplatin | 9 | 1.000 | 5.444 | 1.640 |
| Fibroblast | cisplatin_GENUS | 9 | 1.000 | 2.444 | 1.548 |
| Fibroblast | control | 4 | 0.250 | 0.250 | 0.174 |
| Fibroblast | control_GENUS | 4 | 1.000 | 1.000 | 1.931 |
| Immune | cisplatin | 97 | 0.732 | 2.567 | 1.673 |
| Immune | cisplatin_GENUS | 58 | 0.707 | 1.397 | 1.230 |
| Immune | control | 140 | 0.379 | 0.521 | 0.670 |
| Immune | control_GENUS | 129 | 0.457 | 0.729 | 0.882 |
| Microglia | cisplatin | 907 | 0.936 | 8.451 | 2.074 |
| Microglia | cisplatin_GENUS | 411 | 0.886 | 2.927 | 1.322 |
| Microglia | control | 1494 | 0.432 | 0.664 | 0.504 |
| Microglia | control_GENUS | 1253 | 0.630 | 1.698 | 0.987 |
| Neuron | cisplatin | 22 | 1.000 | 6.636 | 2.080 |
| Neuron | cisplatin_GENUS | 21 | 0.857 | 2.524 | 1.662 |
| Neuron | control | 135 | 0.378 | 0.526 | 0.512 |
| Neuron | control_GENUS | 109 | 0.651 | 1.468 | 1.167 |
| OPC | cisplatin | 122 | 0.943 | 18.943 | 1.730 |
| OPC | cisplatin_GENUS | 72 | 0.819 | 3.056 | 1.225 |
| OPC | control | 280 | 0.371 | 0.621 | 0.433 |
| OPC | control_GENUS | 161 | 0.571 | 3.087 | 0.806 |
| Oligodendrocyte | cisplatin | 855 | 0.918 | 7.868 | 1.885 |
| Oligodendrocyte | cisplatin_GENUS | 486 | 0.848 | 2.566 | 1.186 |
| Oligodendrocyte | control | 1328 | 0.466 | 1.891 | 0.482 |
| Oligodendrocyte | control_GENUS | 870 | 0.644 | 1.299 | 0.897 |
| Pericyte | cisplatin | 88 | 0.818 | 22.818 | 1.778 |
| Pericyte | cisplatin_GENUS | 73 | 0.849 | 2.932 | 1.242 |
| Pericyte | control | 110 | 0.355 | 0.500 | 0.402 |
| Pericyte | control_GENUS | 75 | 0.613 | 2.400 | 0.781 |
| Progenitor | cisplatin | 11 | 0.909 | 9.273 | 1.842 |
| Progenitor | cisplatin_GENUS | 5 | 1.000 | 2.800 | 1.648 |
| Progenitor | control | 74 | 0.405 | 0.541 | 0.439 |
| Progenitor | control_GENUS | 64 | 0.672 | 1.688 | 1.184 |
