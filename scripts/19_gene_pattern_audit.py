"""Task 7: robust gene-pattern and marker-set enrichment audit.

Primary gene lists are frozen consensus lists derived from the complete corrected
integrated-gradient table in the manuscript source package
(SHA-256 c1ab9d911b6aad2d4f6058fb4c649c6ec338aef6aa2c58016d5c32af26b5d7dc).

Consensus derivation:
  * attribution target: squared posterior-mean norm (orthogonally invariant
    within a latent block);
  * baselines: zero and control median;
  * seeds: 0, 1, 2;
  * a gene must appear in the top 50 in at least 4 of the 6 seed-by-baseline
    rankings.

The resulting consensus list sizes are 34 shared, 20 doxorubicin-associated,
and 33 cisplatin-associated genes. This script tests those frozen lists against
the exact 1,500-gene benchmark universe using mouse MSigDB 2025.1 Reactome,
GO Biological Process, and M8 Cell Type Signature collections.

For cross-checking, the script also tests overlap with:
  * 69 unique benchmark genes reaching PyDESeq2 BH-adjusted P<0.05 in at least
    one GSE216146 cell type for cisplatin vs control (Task 4);
  * the top 100 scDisInFact condition-associated gene scores (Task 6).

These cross-checks do not convert the lists into causal or drug-specific
mechanisms. Doxorubicin has one pooled deposited library per arm, so there is
no independent gene-wise doxorubicin differential-expression test.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

R = Path(__file__).resolve().parents[1]
OUT = R / "analysis" / "corrected_20260920"
EVID = R / "data" / "evidence"
OUT.mkdir(parents=True, exist_ok=True)

ATTRIBUTION_SOURCE_SHA256 = (
    "c1ab9d911b6aad2d4f6058fb4c649c6ec338aef6aa2c58016d5c32af26b5d7dc"
)

CONSENSUS = {
    "shared": [
        "Mrc1","Pf4","Ms4a7","Cbr2","Cd163","F13a1","Gpx3","Dab2",
        "mt-Co3","Cd52","Stab1","C5ar1","mt-Nd4","mt-Co2","Apoe","Ttr",
        "Laptm5","Lilrb4a","Ifi27l2a","Bst2","Atp5g1","Psmb8","Cpne2",
        "Ifitm3","Lgals1","Ctsd","Tec","Junb","Cxcl12","Slc6a6","Ptprc",
        "Ccl4","C1qa","Qk",
    ],
    "doxorubicin": [
        "Snhg11","Kcnq1ot1","Meg3","Ttr","Gria2","Rtn1","Grin2b","C1ql3",
        "Opcml","Pcdh9","Nrcam","Dclk1","Peg3","Mir100hg","Zeb2","Pcsk2",
        "Ptprd","Zbtb20","Nfib","Miat",
    ],
    "cisplatin": [
        "Ccn3","Col25a1","Igfbp6","Mgp","Foxc2","Rspo3","Col6a2","mt-Nd4",
        "Ogn","Col6a1","Adamtsl3","Foxd1","Efemp1","S100a10","Bsg","Hes1",
        "Fos","Ttr","Ildr2","Mt1","Col1a2","Tmsb10","Apoe","Tsc22d1",
        "Anxa1","Id3","Mbp","Cbr2","Folr2","Mrc1","Pf4","F13a1","Cd163",
    ],
}

def derive_cisplatin_de_support():
    """Derive the Task 4 cisplatin support set from regenerated full results."""
    path = OUT / "pseudobulk_task4_paired_results.csv.gz"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing; run scripts/16_pseudobulk_paired.py first"
        )
    x = pd.read_csv(path)
    x = x[x.contrast == "cisplatin_vs_control"].copy()
    pair_cols = [
        "replicate_1_log2cpm_difference",
        "replicate_2_log2cpm_difference",
        "replicate_3_log2cpm_difference",
    ]
    pair_sign = np.sign(x[pair_cols].to_numpy())
    model_sign = np.sign(x.log2FoldChange.to_numpy())[:, None]
    x["all_three_pairs_same_direction"] = (pair_sign == model_sign).all(axis=1)
    x["supported_row"] = (
        (x.wald_padj < 0.05) & x.all_three_pairs_same_direction
    )
    support = set(x.loc[x.supported_row, "gene"].astype(str))

    gene = x.groupby("gene", as_index=False).agg(
        n_eligible_cell_types=("cell_type", "nunique"),
        n_fdr005_cell_types=("wald_padj", lambda v: int((v < 0.05).sum())),
        n_direction_consistent_fdr005_cell_types=(
            "supported_row", lambda v: int(np.sum(v))
        ),
        minimum_padj=("wald_padj", "min"),
        max_abs_log2FoldChange=(
            "log2FoldChange", lambda v: float(np.max(np.abs(v)))
        ),
    )
    gene["supported_de_gene"] = (
        gene.n_direction_consistent_fdr005_cell_types > 0
    )
    gene.to_csv(OUT / "task7_cis_pseudobulk_gene_support.csv", index=False)
    return support


SCD_TOP100 = set("""
Ttr mt-Co3 Plp1 Cldn11 Mobp Ptprd Bc1 Gstm1 Ttyh1 Ptgds Apoe Aldoc Mbp
Atp11b Neat1 Mt2 Serpinb1a Nrxn1 Gria2 Kcna1 Arhgap45 Elmo1 Etv1 mt-Co2
Fabp7 Tubb3 Ecrg4 Slc6a6 Mt1 Tmsb10 Ahcyl2 Adamts1 Sez6l2 Slc6a1 Pcdh9
Enpp2 mt-Nd4 Brd9 Slc38a1 Arhgap25 Prlr Hist1h1e Cdkn1a Peg3 Ptgs1 Npas3
Ntng1 Arl4a Slco1c1 Dclk1 Garnl3 Ntm Sema6d Sox2ot 4933429O19Rik Etnppl
Itgad Spint2 Rmst Bcan Slc8a1 A230006K03Rik Sat1 Kl Ndst4 Camk1d Plek Spef2
BC028528 Tsc22d1 Ifi203 S100a10 Kif5a Crip1 Eprs Gpcpd1 Cntnap2 Peg10 Mdh1
Col4a1 Gm43689 Chst2 Pcdh15 Spock2 6330403K07Rik 9330185C12Rik Slc1a2
4833427G06Rik Dbndd2 Col6a2 Prox1 Hspb1 Pcp4 Hist1h2ae Dscam Ppp1r15a
Snap25 Epha4 Zfp536 Homer2
""".split())

LIBRARIES = {
    "Reactome": "m2.cp.reactome.v2025.1.Mm.symbols.gmt",
    "GO_BP": "m5.go.bp.v2025.1.Mm.symbols.gmt",
    "M8_cell_type": "m8.all.v2025.1.Mm.symbols.gmt",
}


def read_gmt(path: Path):
    out = []
    for line in path.read_text().splitlines():
        fields = line.rstrip().split("\t")
        if len(fields) < 3:
            continue
        out.append((fields[0], set(fields[2:])))
    return out


def enrichment(selected: set[str], universe: set[str], library: str, path: Path):
    rows = []
    N = len(universe)
    n = len(selected)
    for term, genes in read_gmt(path):
        g = genes & universe
        if not 10 <= len(g) <= 500:
            continue
        overlap = selected & g
        k = len(overlap)
        p = float(hypergeom.sf(k - 1, N, len(g), n)) if k else 1.0
        rows.append(
            {
                "library": library,
                "term": term,
                "overlap": k,
                "set_size_in_universe": len(g),
                "list_size": n,
                "universe_size": N,
                "fold_enrichment": (
                    k / (n * len(g) / N) if k and len(g) else 0.0
                ),
                "p_value": p,
                "genes": "|".join(sorted(overlap)),
            }
        )
    if rows:
        q = multipletests([x["p_value"] for x in rows], method="fdr_bh")[1]
        for x, qv in zip(rows, q):
            x["q_value"] = float(qv)
    return rows


def overlap_test(selected, reference, universe_size):
    k = len(selected & reference)
    p = (
        float(
            hypergeom.sf(
                k - 1,
                universe_size,
                len(reference),
                len(selected),
            )
        )
        if k
        else 1.0
    )
    return k, p, "|".join(sorted(selected & reference))


def main():
    universe = set(
        pd.read_csv(OUT / "benchmark_gene_universe.csv")["gene"].astype(str)
    )
    assert len(universe) == 1500
    cisplatin_de = derive_cisplatin_de_support() & universe
    assert len(cisplatin_de) == 69, len(cisplatin_de)
    assert len(SCD_TOP100) == 100

    consensus_rows = []
    for axis, genes in CONSENSUS.items():
        assert len(genes) == len(set(genes))
        missing = set(genes) - universe
        assert not missing, (axis, missing)
        for rank, gene in enumerate(genes, start=1):
            consensus_rows.append(
                {
                    "axis": axis,
                    "consensus_rank_by_mean_abs_ig": rank,
                    "gene": gene,
                    "consensus_rule": "top50_in_at_least_4_of_6_seed_by_baseline_rankings",
                    "primary_target": "posterior_mean_block_norm_squared",
                    "baselines": "zero|control_median",
                    "seeds": "0|1|2",
                }
            )
    pd.DataFrame(consensus_rows).to_csv(
        OUT / "task7_consensus_genes.csv", index=False
    )

    cross = []
    for axis, genes in CONSENSUS.items():
        selected = set(genes)
        k_de, p_de, g_de = overlap_test(selected, cisplatin_de, len(universe))
        k_scd, p_scd, g_scd = overlap_test(selected, SCD_TOP100, len(universe))
        cross.append(
            {
                "axis": axis,
                "n_consensus_genes": len(selected),
                "n_overlap_cisplatin_DE_unique_genes": k_de,
                "cisplatin_DE_overlap_genes": g_de,
                "hypergeom_p_cisplatin_DE_overlap": p_de,
                "n_overlap_scdisinfact_top100": k_scd,
                "scdisinfact_top100_overlap_genes": g_scd,
                "hypergeom_p_scdisinfact_top100_overlap": p_scd,
            }
        )
    pd.DataFrame(cross).to_csv(
        OUT / "task7_cross_validation_summary.csv", index=False
    )

    all_enrich = []
    for axis, genes in CONSENSUS.items():
        selected = set(genes)
        for library, filename in LIBRARIES.items():
            path = EVID / filename
            if not path.exists():
                raise FileNotFoundError(path)
            rows = enrichment(selected, universe, library, path)
            for row in rows:
                row["axis"] = axis
            all_enrich.extend(rows)

    enrich = pd.DataFrame(all_enrich)
    enrich.to_csv(OUT / "task7_consensus_enrichment.csv", index=False)

    top = (
        enrich.sort_values(
            ["axis", "library", "q_value", "p_value", "term"],
            ascending=[True, True, True, True, True],
        )
        .groupby(["axis", "library"], as_index=False, group_keys=False)
        .head(5)
    )
    top.to_csv(OUT / "task7_top_enrichment_terms.csv", index=False)

    sig = (
        enrich.assign(significant=enrich.q_value < 0.05)
        .groupby(["axis", "library"], as_index=False)
        .agg(
            n_tested=("term", "size"),
            n_fdr_005=("significant", "sum"),
            minimum_q=("q_value", "min"),
        )
    )
    sig.to_csv(OUT / "task7_enrichment_summary.csv", index=False)

    design = {
        "attribution_source_sha256": ATTRIBUTION_SOURCE_SHA256,
        "primary_target": "squared norm of posterior-mean latent block",
        "target_rationale": (
            "invariant to sign flips and orthogonal rotations within a latent block"
        ),
        "ranking_statistic": "mean absolute integrated gradient within each seed/baseline",
        "consensus_rule": "top 50 in at least 4 of 6 seed-by-baseline rankings",
        "seeds": [0, 1, 2],
        "baselines": ["zero", "control_median"],
        "consensus_sizes": {k: len(v) for k, v in CONSENSUS.items()},
        "universe": "exact 1,500-gene canonical benchmark universe",
        "gene_set_resources": {
            "Reactome": "MSigDB mouse 2025.1 m2.cp.reactome",
            "GO_BP": "MSigDB mouse 2025.1 m5.go.bp",
            "M8_cell_type": "MSigDB mouse 2025.1 M8 cell type signatures",
        },
        "set_size_filter": "10-500 genes after intersection with model universe",
        "test": "one-sided hypergeometric overrepresentation",
        "multiple_testing": "BH separately within axis x resource family",
        "frozen_inputs": {
            "attribution_table_sha256": ATTRIBUTION_SOURCE_SHA256,
            "scdisinfact_top100": "derived from Task 6 across-seed mean condition-associated scores",
        },
        "cross_checks": {
            "cisplatin_DE": (
                "69 unique benchmark genes with paired PyDESeq2 adjusted P<0.05 "
                "in at least one GSE216146 cell type and all three paired "
                "logCPM effects agreeing with the fitted effect direction"
            ),
            "scDisInFact": (
                "top 100 mean condition-associated gene scores across Task 6 seeds"
            ),
        },
        "limitations": [
            "Consensus selection improves robustness to seed and baseline but does not create biological replication.",
            "Cisplatin differential expression is based on three paired deposited replicate blocks and is model-based exploratory evidence.",
            "No gene-wise inferential doxorubicin validation is possible with one pooled library per arm.",
            "Drug identity remains partially aliased with study and sequencing modality.",
            "Gene-set overlap does not imply pathway activation or causal mechanism.",
        ],
    }
    (OUT / "task7_gene_pattern_design.json").write_text(
        json.dumps(design, indent=2) + "\n"
    )

    print("===== TASK 7 CROSS VALIDATION =====")
    print(pd.DataFrame(cross).to_string(index=False))
    print("\n===== TASK 7 ENRICHMENT SUMMARY =====")
    print(sig.to_string(index=False))
    print("\n===== TASK 7 TOP TERMS =====")
    print(
        top[
            [
                "axis",
                "library",
                "term",
                "overlap",
                "fold_enrichment",
                "q_value",
                "genes",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()