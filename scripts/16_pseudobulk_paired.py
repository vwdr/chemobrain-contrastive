"""Task 4: paired library-level pseudobulk analysis.

Primary cisplatin-study analysis:
  * aggregate raw counts by source cell type and deposited library;
  * preserve the three GEO replicate blocks;
  * for each contrast, fit a paired negative-binomial model with PyDESeq2:
        ~ replicate + condition
  * report DESeq2-style Wald statistics and BH-adjusted p values;
  * additionally compute paired log2(CPM+1) effect sizes and the exact 2^3
    sign-flip sensitivity reference. With three pairs, the two-sided exact
    reference has a minimum attainable p value of 0.25.

Doxorubicin-study analysis:
  * one pooled library per arm, so report descriptive log2(CPM+1) differences
    only. No biological inferential p values are computed.

This script is intentionally separate from scripts/10_pseudobulk.py, which is
retained as the historical unpaired Welch sensitivity analysis.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

R = Path(__file__).resolve().parents[1]
T = R / "analysis" / "corrected_20260920"
T.mkdir(parents=True, exist_ok=True)

MIN_CELLS_PER_LIBRARY = 20
MIN_COUNT = 10
MIN_LIBRARIES_WITH_COUNT = 3

REPLICATE_ARM_TO_LIBRARY = {
    "replicate_1": {
        "control": "D20-6407",
        "control_GENUS": "D20-6408",
        "cisplatin": "D20-6409",
        "cisplatin_GENUS": "D20-6410",
    },
    "replicate_2": {
        "control": "D21-2746",
        "control_GENUS": "D21-2748",
        "cisplatin": "D21-2750",
        "cisplatin_GENUS": "D21-2752",
    },
    "replicate_3": {
        "control": "D21-2747",
        "control_GENUS": "D21-2749",
        "cisplatin": "D21-2751",
        "cisplatin_GENUS": "D21-2753",
    },
}

CONTRASTS = [
    ("cisplatin", "control"),
    ("control_GENUS", "control"),
    ("cisplatin_GENUS", "cisplatin"),
]

DOX_CONTRASTS = [
    ("doxorubicin", "control"),
    ("doxorubicin_ACY1083", "doxorubicin"),
]

DOX_LIBRARIES = {
    "control": "GSM8367861",
    "doxorubicin": "GSM8367862",
    "doxorubicin_ACY1083": "GSM8367863",
}


def aggregate_library_counts(a: ad.AnnData):
    groups = {}
    qc_rows = []
    for (cell_type, sample_id, arm), ix in a.obs.groupby(
        ["source_cell_type", "sample_id", "arm"], observed=True
    ).indices.items():
        n_cells = len(ix)
        qc_rows.append(
            {
                "study": str(a.obs.study.iloc[ix[0]]),
                "cell_type": str(cell_type),
                "sample_id": str(sample_id),
                "arm": str(arm),
                "n_cells": int(n_cells),
                "included_ge_20_cells": bool(n_cells >= MIN_CELLS_PER_LIBRARY),
            }
        )
        if n_cells >= MIN_CELLS_PER_LIBRARY:
            counts = np.asarray(a.X[ix].sum(axis=0)).ravel()
            if not np.allclose(counts, np.round(counts)):
                raise ValueError("Pseudobulk input contains non-integer counts.")
            groups[(str(cell_type), str(sample_id), str(arm))] = np.round(
                counts
            ).astype(np.int64)
    return groups, pd.DataFrame(qc_rows)


def library_rows_for_contrast(cell_type, tested, reference, groups):
    rows = []
    meta = []
    for rep in ["replicate_1", "replicate_2", "replicate_3"]:
        for condition, arm in [("reference", reference), ("tested", tested)]:
            sid = REPLICATE_ARM_TO_LIBRARY[rep][arm]
            key = (cell_type, sid, arm)
            if key not in groups:
                return None, None
            rows.append(groups[key])
            meta.append(
                {
                    "sample": sid,
                    "replicate": rep,
                    "condition": condition,
                    "arm": arm,
                }
            )
    return np.asarray(rows, dtype=np.int64), pd.DataFrame(meta).set_index("sample")


def exact_paired_signflip(logcpm, metadata):
    # metadata rows are generated as ref/test within replicate 1, 2, 3.
    diffs = []
    for rep in ["replicate_1", "replicate_2", "replicate_3"]:
        rr = np.flatnonzero(
            (metadata.replicate.to_numpy() == rep)
            & (metadata.condition.to_numpy() == "reference")
        )
        tt = np.flatnonzero(
            (metadata.replicate.to_numpy() == rep)
            & (metadata.condition.to_numpy() == "tested")
        )
        assert len(rr) == 1 and len(tt) == 1
        diffs.append(logcpm[tt[0]] - logcpm[rr[0]])
    diffs = np.asarray(diffs)

    observed = diffs.mean(axis=0)
    reference = np.asarray(
        [
            (diffs * np.asarray(signs)[:, None]).mean(axis=0)
            for signs in itertools.product([-1, 1], repeat=3)
        ]
    )
    p = np.mean(
        np.abs(reference) >= np.abs(observed)[None, :] - 1e-12,
        axis=0,
    )
    return observed, p, diffs


def fit_pydeseq2(counts, metadata, genes):
    counts_df = pd.DataFrame(
        counts,
        index=metadata.index,
        columns=genes,
    )
    meta = metadata[["replicate", "condition"]].copy()
    meta["replicate"] = pd.Categorical(
        meta["replicate"],
        categories=["replicate_1", "replicate_2", "replicate_3"],
    )
    meta["condition"] = pd.Categorical(
        meta["condition"],
        categories=["reference", "tested"],
    )

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=meta,
        design="~ replicate + condition",
        refit_cooks=False,
        n_cpus=1,
        quiet=True,
    )
    dds.deseq2()

    ds = DeseqStats(
        dds,
        contrast=["condition", "tested", "reference"],
        alpha=0.05,
        cooks_filter=False,
        independent_filter=True,
        n_cpus=1,
        quiet=True,
    )
    ds.summary()
    out = ds.results_df.copy()
    out.index = out.index.astype(str)
    out.index.name = "gene"
    return out.reset_index()


def analyze_cisplatin(a, groups):
    gene_names = np.asarray(a.var_names.astype(str))
    result_rows = []
    summary_rows = []
    skipped_rows = []

    cell_types = sorted(a.obs.source_cell_type.astype(str).unique())
    for cell_type in cell_types:
        for tested, reference in CONTRASTS:
            contrast = f"{tested}_vs_{reference}"
            raw, metadata = library_rows_for_contrast(
                cell_type, tested, reference, groups
            )
            if raw is None:
                skipped_rows.append(
                    {
                        "cell_type": cell_type,
                        "contrast": contrast,
                        "reason": "one_or_more_libraries_below_20_cells",
                    }
                )
                continue

            keep = (raw >= MIN_COUNT).sum(axis=0) >= MIN_LIBRARIES_WITH_COUNT
            if keep.sum() == 0:
                skipped_rows.append(
                    {
                        "cell_type": cell_type,
                        "contrast": contrast,
                        "reason": "no_genes_pass_count_filter",
                    }
                )
                continue

            raw_k = raw[:, keep]
            genes = gene_names[keep]

            libsize = raw_k.sum(axis=1).astype(float)
            logcpm = np.log2(
                raw_k / libsize[:, None] * 1e6 + 1.0
            )
            paired_effect, exact_p, pair_diffs = exact_paired_signflip(
                logcpm, metadata
            )
            exact_q = multipletests(exact_p, method="fdr_bh")[1]

            ref_ix = metadata.condition.to_numpy() == "reference"
            test_ix = metadata.condition.to_numpy() == "tested"
            welch_stat, welch_p = ttest_ind(
                logcpm[test_ix],
                logcpm[ref_ix],
                axis=0,
                equal_var=False,
            )
            welch_p = np.nan_to_num(welch_p, nan=1.0)
            welch_q = multipletests(welch_p, method="fdr_bh")[1]

            de = fit_pydeseq2(raw_k, metadata, genes)
            de = de.set_index("gene").reindex(genes).reset_index()

            for j, gene in enumerate(genes):
                row = {
                    "cell_type": cell_type,
                    "contrast": contrast,
                    "tested_arm": tested,
                    "reference_arm": reference,
                    "gene": str(gene),
                    "n_pairs": 3,
                    "paired_mean_log2cpm_difference": float(
                        paired_effect[j]
                    ),
                    "replicate_1_log2cpm_difference": float(
                        pair_diffs[0, j]
                    ),
                    "replicate_2_log2cpm_difference": float(
                        pair_diffs[1, j]
                    ),
                    "replicate_3_log2cpm_difference": float(
                        pair_diffs[2, j]
                    ),
                    "exact_two_sided_signflip_p": float(exact_p[j]),
                    "exact_signflip_q_bh": float(exact_q[j]),
                    "welch_t_sensitivity": float(welch_stat[j]),
                    "welch_p_sensitivity": float(welch_p[j]),
                    "welch_q_sensitivity": float(welch_q[j]),
                    "baseMean": float(de.loc[j, "baseMean"]),
                    "log2FoldChange": float(de.loc[j, "log2FoldChange"]),
                    "lfcSE": float(de.loc[j, "lfcSE"]),
                    "wald_stat": float(de.loc[j, "stat"]),
                    "wald_pvalue": float(de.loc[j, "pvalue"]),
                    "wald_padj": (
                        float(de.loc[j, "padj"])
                        if pd.notna(de.loc[j, "padj"])
                        else np.nan
                    ),
                }
                result_rows.append(row)

            padj = de["padj"].to_numpy(dtype=float)
            pval = de["pvalue"].to_numpy(dtype=float)
            lfc = de["log2FoldChange"].to_numpy(dtype=float)
            finite_padj = padj[np.isfinite(padj)]
            corr = np.corrcoef(
                lfc,
                paired_effect,
            )[0, 1] if len(lfc) > 1 else np.nan

            summary_rows.append(
                {
                    "cell_type": cell_type,
                    "contrast": contrast,
                    "n_pairs": 3,
                    "n_genes_tested": int(len(genes)),
                    "n_wald_p_lt_0p05": int(
                        np.sum(np.isfinite(pval) & (pval < 0.05))
                    ),
                    "n_wald_padj_lt_0p05": int(
                        np.sum(np.isfinite(padj) & (padj < 0.05))
                    ),
                    "minimum_wald_padj": (
                        float(np.min(finite_padj))
                        if len(finite_padj)
                        else np.nan
                    ),
                    "n_exact_signflip_p_lt_0p05": int(
                        np.sum(exact_p < 0.05)
                    ),
                    "minimum_exact_two_sided_p": float(
                        np.min(exact_p)
                    ),
                    "theoretical_exact_two_sided_p_floor": 0.25,
                    "n_welch_q_lt_0p05": int(
                        np.sum(welch_q < 0.05)
                    ),
                    "deseq2_vs_paired_logcpm_effect_correlation": float(
                        corr
                    ),
                }
            )

            print(
                f"{cell_type:24s} {contrast:36s} "
                f"genes={len(genes):5d} "
                f"padj<.05={np.sum(np.isfinite(padj) & (padj < 0.05)):4d}",
                flush=True,
            )

    return (
        pd.DataFrame(result_rows),
        pd.DataFrame(summary_rows),
        pd.DataFrame(skipped_rows),
    )


def analyze_dox(a, groups):
    genes = np.asarray(a.var_names.astype(str))
    rows = []
    summary = []
    cell_types = sorted(a.obs.source_cell_type.astype(str).unique())

    for cell_type in cell_types:
        for tested, reference in DOX_CONTRASTS:
            sid_t = DOX_LIBRARIES[tested]
            sid_r = DOX_LIBRARIES[reference]
            kt = (cell_type, sid_t, tested)
            kr = (cell_type, sid_r, reference)
            if kt not in groups or kr not in groups:
                continue

            raw = np.asarray([groups[kr], groups[kt]], dtype=np.int64)
            keep = np.max(raw, axis=0) >= MIN_COUNT
            if not keep.any():
                continue
            raw_k = raw[:, keep]
            gg = genes[keep]
            libsize = raw_k.sum(axis=1).astype(float)
            logcpm = np.log2(
                raw_k / libsize[:, None] * 1e6 + 1.0
            )
            effect = logcpm[1] - logcpm[0]

            for gene, e, rc, tc in zip(
                gg, effect, raw_k[0], raw_k[1]
            ):
                rows.append(
                    {
                        "cell_type": cell_type,
                        "contrast": f"{tested}_vs_{reference}",
                        "gene": str(gene),
                        "reference_library": sid_r,
                        "tested_library": sid_t,
                        "reference_count": int(rc),
                        "tested_count": int(tc),
                        "descriptive_log2cpm_difference": float(e),
                        "inferential_pvalue": np.nan,
                        "inference_note": "descriptive_only_one_pooled_library_per_arm",
                    }
                )

            summary.append(
                {
                    "cell_type": cell_type,
                    "contrast": f"{tested}_vs_{reference}",
                    "n_genes_described": int(len(gg)),
                    "reference_library": sid_r,
                    "tested_library": sid_t,
                    "inferential_test_performed": False,
                    "reason": "one_pooled_library_per_arm",
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(summary)


def main():
    a = ad.read_h5ad(R / "data" / "processed" / "recovered_counts.h5ad")
    if not np.allclose(a.X.data if hasattr(a.X, "data") else a.X,
                       np.round(a.X.data if hasattr(a.X, "data") else a.X)):
        raise ValueError("Recovered cohort matrix is not integer counts.")

    groups, qc = aggregate_library_counts(a)
    qc.to_csv(T / "pseudobulk_task4_library_qc.csv", index=False)

    cis = a[a.obs.study == "GSE216146"].copy()
    cis_groups = {
        k: v for k, v in groups.items()
        if k[1].startswith("D20-") or k[1].startswith("D21-")
    }
    results, summary, skipped = analyze_cisplatin(cis, cis_groups)
    results.to_csv(
        T / "pseudobulk_task4_paired_results.csv.gz",
        index=False,
        compression="gzip",
    )
    summary.to_csv(T / "pseudobulk_task4_summary.csv", index=False)
    skipped.to_csv(T / "pseudobulk_task4_skipped.csv", index=False)

    # Compact robustness summary for manuscript interpretation.
    pair_cols = [
        "replicate_1_log2cpm_difference",
        "replicate_2_log2cpm_difference",
        "replicate_3_log2cpm_difference",
    ]
    pair_sign = np.sign(results[pair_cols].to_numpy())
    model_sign = np.sign(results.log2FoldChange.to_numpy())[:, None]
    results["n_pairs_same_direction_as_model"] = (
        pair_sign == model_sign
    ).sum(axis=1)

    robustness_rows = []
    top_rows = []
    for (cell_type, contrast), g in results.groupby(
        ["cell_type", "contrast"], sort=True
    ):
        hits = g[g.wald_padj < 0.05].copy()
        effect_spearman = g[
            ["log2FoldChange", "paired_mean_log2cpm_difference"]
        ].corr(method="spearman").iloc[0, 1]
        effect_pearson = g[
            ["log2FoldChange", "paired_mean_log2cpm_difference"]
        ].corr(method="pearson").iloc[0, 1]
        robustness_rows.append(
            {
                "cell_type": cell_type,
                "contrast": contrast,
                "n_genes_tested": int(len(g)),
                "n_deseq2_fdr_005": int(len(hits)),
                "n_deseq2_fdr_010": int((g.wald_padj < 0.10).sum()),
                "n_exact_signflip_p_lt_005": int(
                    (g.exact_two_sided_signflip_p < 0.05).sum()
                ),
                "minimum_exact_two_sided_p": float(
                    g.exact_two_sided_signflip_p.min()
                ),
                "n_welch_fdr_005": int(
                    (g.welch_q_sensitivity < 0.05).sum()
                ),
                "effect_spearman_deseq2_vs_paired_logcpm": float(
                    effect_spearman
                ),
                "effect_pearson_deseq2_vs_paired_logcpm": float(
                    effect_pearson
                ),
                "n_deseq2_hits_all3_pairs_same_direction": int(
                    (hits.n_pairs_same_direction_as_model == 3).sum()
                ),
                "n_deseq2_hits_2of3_pairs_same_direction": int(
                    (hits.n_pairs_same_direction_as_model == 2).sum()
                ),
            }
        )
        if len(hits):
            top_rows.append(
                hits.nsmallest(5, "wald_padj")[
                    [
                        "cell_type",
                        "contrast",
                        "gene",
                        "log2FoldChange",
                        "wald_padj",
                        "paired_mean_log2cpm_difference",
                        *pair_cols,
                        "n_pairs_same_direction_as_model",
                        "welch_q_sensitivity",
                        "exact_two_sided_signflip_p",
                    ]
                ]
            )

    robustness = pd.DataFrame(robustness_rows)
    robustness.to_csv(
        T / "pseudobulk_task4_robustness.csv", index=False
    )
    if top_rows:
        pd.concat(top_rows, ignore_index=True).to_csv(
            T / "pseudobulk_task4_top_hits.csv", index=False
        )

    dox = a[a.obs.study == "GSE271055"].copy()
    dox_groups = {
        k: v for k, v in groups.items()
        if k[1].startswith("GSM836")
    }
    dox_results, dox_summary = analyze_dox(dox, dox_groups)
    dox_results.to_csv(
        T / "pseudobulk_task4_dox_descriptive.csv.gz",
        index=False,
        compression="gzip",
    )
    dox_summary.to_csv(
        T / "pseudobulk_task4_dox_summary.csv",
        index=False,
    )

    method = {
        "package": "PyDESeq2",
        "package_version": "0.5.4",
        "cisplatin_design": "~ replicate + condition",
        "cisplatin_replicates": 3,
        "contrasts": [
            f"{a}_vs_{b}" for a, b in CONTRASTS
        ],
        "minimum_cells_per_library_cell_type": MIN_CELLS_PER_LIBRARY,
        "gene_filter": (
            f"raw count >= {MIN_COUNT} in at least "
            f"{MIN_LIBRARIES_WITH_COUNT} of 6 libraries"
        ),
        "deseq2_cooks_filter": False,
        "deseq2_independent_filter": True,
        "paired_exact_sensitivity": (
            "all 2^3 within-replicate sign flips of paired "
            "log2(CPM+1) differences"
        ),
        "paired_exact_two_sided_p_floor": 0.25,
        "welch_analysis": (
            "reported only as a sensitivity column for continuity "
            "with the historical analysis; it is not the primary model"
        ),
        "doxorubicin_inference": (
            "descriptive only because one pooled deposited library "
            "is available per arm"
        ),
    }
    (T / "pseudobulk_task4_method.json").write_text(
        json.dumps(method, indent=2) + "\n"
    )

    print("\n===== TASK 4 CISPLATIN SUMMARY =====")
    print(summary.to_string(index=False))
    print("\n===== TASK 4 DOX SUMMARY =====")
    print(dox_summary.to_string(index=False))


if __name__ == "__main__":
    main()