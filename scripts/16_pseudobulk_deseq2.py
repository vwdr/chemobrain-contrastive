"""Library-level pseudobulk differential-expression analysis.

Primary inferential analysis for GSE216146 uses integer pseudobulk counts at the
deposited-library level and a paired replicate-block DESeq2-style model through
PyDESeq2.  The three GEO replicate blocks are included as a blocking factor:

    ~ replicate + condition

Only cell types with >=20 cells in every library required for a given 3-vs-3
paired contrast are analyzed. Genes must have >=10 pseudobulk counts in at least
3 of the 6 libraries before model fitting.

The prior Welch test on log2(CPM+1) is retained as a sensitivity analysis.  It is
not the primary inferential model.

GSE271055 contains one pooled library per arm.  It is summarized descriptively
only; no gene-wise P values are reported.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

R = Path(__file__).resolve().parents[1]
T = R / "analysis" / "corrected_20260920"
DATA = R / "data" / "processed" / "recovered_counts.h5ad"
T.mkdir(parents=True, exist_ok=True)

MIN_CELLS = 20
MIN_COUNT = 10
MIN_LIBRARIES_WITH_COUNT = 3
ALPHA = 0.05

REPLICATE = {
    "D20-6407": "replicate_1",
    "D20-6408": "replicate_1",
    "D20-6409": "replicate_1",
    "D20-6410": "replicate_1",
    "D21-2746": "replicate_2",
    "D21-2748": "replicate_2",
    "D21-2750": "replicate_2",
    "D21-2752": "replicate_2",
    "D21-2747": "replicate_3",
    "D21-2749": "replicate_3",
    "D21-2751": "replicate_3",
    "D21-2753": "replicate_3",
}

CONTRASTS = [
    ("cisplatin_vs_control", "cisplatin", "control"),
    ("cisplatin_GENUS_vs_cisplatin", "cisplatin_GENUS", "cisplatin"),
    ("control_GENUS_vs_control", "control_GENUS", "control"),
]


def aggregate_counts(a: ad.AnnData) -> tuple[dict, pd.DataFrame]:
    groups = {}
    qc = []
    for (cell_type, sample_id, arm), ix in a.obs.groupby(
        ["source_cell_type", "sample_id", "arm"], observed=True
    ).indices.items():
        vec = np.asarray(a.X[ix].sum(0)).ravel()
        assert np.allclose(vec, np.round(vec))
        groups[(str(cell_type), str(sample_id), str(arm))] = np.round(vec).astype(
            np.int64
        )
        qc.append(
            {
                "cell_type": str(cell_type),
                "sample_id": str(sample_id),
                "arm": str(arm),
                "n_cells": int(len(ix)),
                "library_counts": int(vec.sum()),
            }
        )
    return groups, pd.DataFrame(qc)


def samples_for_contrast(
    qc: pd.DataFrame, cell_type: str, tested: str, reference: str
) -> pd.DataFrame | None:
    q = qc[
        (qc.cell_type == cell_type) & qc.arm.isin([tested, reference])
    ].copy()
    q = q[q.sample_id.isin(REPLICATE)].copy()
    q["replicate"] = q.sample_id.map(REPLICATE)

    wanted = []
    for rep in ("replicate_1", "replicate_2", "replicate_3"):
        for condition in (reference, tested):
            hit = q[(q.replicate == rep) & (q.arm == condition)]
            if len(hit) != 1:
                return None
            if int(hit.n_cells.iloc[0]) < MIN_CELLS:
                return None
            wanted.append(hit.iloc[0])
    return pd.DataFrame(wanted).reset_index(drop=True)


def fit_deseq(
    count_df: pd.DataFrame,
    metadata: pd.DataFrame,
    tested: str,
    reference: str,
) -> pd.DataFrame:
    # Preserve explicit paired blocking; contrast controls tested/reference direction.
    meta = metadata[["replicate", "condition"]].copy()
    meta["replicate"] = pd.Categorical(
        meta["replicate"],
        categories=["replicate_1", "replicate_2", "replicate_3"],
        ordered=True,
    )
    meta["condition"] = pd.Categorical(
        meta["condition"], categories=[reference, tested], ordered=True
    )

    dds = DeseqDataSet(
        counts=count_df,
        metadata=meta,
        design="~ replicate + condition",
        refit_cooks=False,
        quiet=True,
        n_cpus=1,
    )
    dds.deseq2()

    stats = DeseqStats(
        dds,
        contrast=["condition", tested, reference],
        alpha=ALPHA,
        cooks_filter=True,
        independent_filter=True,
        quiet=True,
        n_cpus=1,
    )
    stats.summary()
    out = stats.results_df.copy()
    out.index.name = "gene"
    return out.reset_index()


def welch_sensitivity(
    raw: np.ndarray, genes: np.ndarray
) -> pd.DataFrame:
    cpm = raw / raw.sum(axis=1, keepdims=True) * 1e6
    v = np.log2(cpm + 1)
    # metadata order is reference rep1-3 interleaved with tested in caller;
    # select using explicit positions supplied there instead of relying on halves.
    return pd.DataFrame({"gene": genes})


def main() -> None:
    a = ad.read_h5ad(DATA)
    assert np.allclose(a.X.data if hasattr(a.X, "data") else a.X, np.round(a.X.data if hasattr(a.X, "data") else a.X))

    cis = a[a.obs.study == "GSE216146"].copy()
    groups, qc = aggregate_counts(cis)
    qc.to_csv(T / "pseudobulk_library_qc_v2.csv", index=False)

    result_rows = []
    summary_rows = []
    sensitivity_rows = []
    eligible_rows = []

    genes_all = np.asarray(cis.var_names.astype(str))

    for cell_type in sorted(cis.obs.source_cell_type.astype(str).unique()):
        for contrast_name, tested, reference in CONTRASTS:
            sm = samples_for_contrast(qc, cell_type, tested, reference)
            if sm is None:
                eligible_rows.append(
                    {
                        "cell_type": cell_type,
                        "contrast": contrast_name,
                        "eligible": False,
                        "reason": "missing library or <20 cells in at least one required library",
                    }
                )
                continue

            eligible_rows.append(
                {
                    "cell_type": cell_type,
                    "contrast": contrast_name,
                    "eligible": True,
                    "reason": "",
                }
            )

            raw = np.vstack(
                [
                    groups[(cell_type, str(r.sample_id), str(r.arm))]
                    for _, r in sm.iterrows()
                ]
            )
            keep = (raw >= MIN_COUNT).sum(axis=0) >= MIN_LIBRARIES_WITH_COUNT
            raw = raw[:, keep]
            genes = genes_all[keep]
            if len(genes) == 0:
                continue

            sample_names = sm.sample_id.astype(str).tolist()
            count_df = pd.DataFrame(raw, index=sample_names, columns=genes)
            meta = pd.DataFrame(
                {
                    "replicate": sm.replicate.astype(str).values,
                    "condition": sm.arm.astype(str).values,
                },
                index=sample_names,
            )

            res = fit_deseq(count_df, meta, tested, reference)
            res.insert(0, "contrast", contrast_name)
            res.insert(0, "cell_type", cell_type)
            result_rows.append(res)

            # Welch sensitivity using the same six libraries and same filtered genes.
            cond = meta.condition.to_numpy()
            tmask = cond == tested
            rmask = cond == reference
            cpm = raw / raw.sum(axis=1, keepdims=True) * 1e6
            lv = np.log2(cpm + 1)
            stat, p = ttest_ind(
                lv[tmask], lv[rmask], axis=0, equal_var=False, nan_policy="omit"
            )
            p = np.nan_to_num(p, nan=1.0)
            welch_effect = lv[tmask].mean(axis=0) - lv[rmask].mean(axis=0)
            sensitivity_rows.append(
                pd.DataFrame(
                    {
                        "cell_type": cell_type,
                        "contrast": contrast_name,
                        "gene": genes,
                        "mean_log2CPM_difference": welch_effect,
                        "welch_t": stat,
                        "welch_p_value": p,
                    }
                )
            )

            finite_padj = res["padj"].dropna()
            summary_rows.append(
                {
                    "cell_type": cell_type,
                    "contrast": contrast_name,
                    "n_libraries": 6,
                    "n_pairs": 3,
                    "n_genes_prefiltered": int(len(genes)),
                    "n_genes_with_padj": int(len(finite_padj)),
                    "n_fdr_005": int((finite_padj < 0.05).sum()),
                    "n_fdr_010": int((finite_padj < 0.10).sum()),
                    "min_padj": (
                        float(finite_padj.min()) if len(finite_padj) else np.nan
                    ),
                    "max_abs_log2FoldChange": float(
                        res.log2FoldChange.abs().max()
                    ),
                }
            )
            print(
                cell_type,
                contrast_name,
                "genes=", len(genes),
                "FDR<.05=", summary_rows[-1]["n_fdr_005"],
                flush=True,
            )

    results = pd.concat(result_rows, ignore_index=True)
    sensitivity = pd.concat(sensitivity_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    eligibility = pd.DataFrame(eligible_rows)

    merged = results.merge(
        sensitivity,
        on=["cell_type", "contrast", "gene"],
        how="left",
        validate="one_to_one",
    )
    merged.to_csv(T / "pseudobulk_deseq2_results.csv.gz", index=False)
    sensitivity.to_csv(T / "pseudobulk_welch_sensitivity_v2.csv.gz", index=False)
    summary.to_csv(T / "pseudobulk_deseq2_summary.csv", index=False)
    eligibility.to_csv(T / "pseudobulk_deseq2_eligibility.csv", index=False)

    # Contrast-level agreement between primary DESeq2 effect sizes and Welch effects.
    agree = []
    for (cell_type, contrast), g in merged.groupby(["cell_type", "contrast"]):
        ok = g.log2FoldChange.notna() & g.mean_log2CPM_difference.notna()
        rho = (
            g.loc[ok, ["log2FoldChange", "mean_log2CPM_difference"]]
            .corr(method="spearman")
            .iloc[0, 1]
            if ok.sum() > 2
            else np.nan
        )
        agree.append(
            {
                "cell_type": cell_type,
                "contrast": contrast,
                "n_genes": int(ok.sum()),
                "spearman_deseq2_vs_welch_effect": float(rho),
            }
        )
    pd.DataFrame(agree).to_csv(
        T / "pseudobulk_method_agreement.csv", index=False
    )

    # GSE271055: one pooled deposited library per arm -> descriptive only.
    dox = a[a.obs.study == "GSE271055"].copy()
    dox_groups, dox_qc = aggregate_counts(dox)
    dox_qc.to_csv(T / "pseudobulk_dox_library_qc.csv", index=False)

    desc_rows = []
    for cell_type in sorted(dox.obs.source_cell_type.astype(str).unique()):
        by_arm = {}
        for arm in ("control", "doxorubicin", "doxorubicin_ACY1083"):
            hits = dox_qc[
                (dox_qc.cell_type == cell_type) & (dox_qc.arm == arm)
            ]
            if len(hits) != 1 or int(hits.n_cells.iloc[0]) < MIN_CELLS:
                continue
            sid = str(hits.sample_id.iloc[0])
            by_arm[arm] = dox_groups[(cell_type, sid, arm)]
        if set(by_arm) != {
            "control",
            "doxorubicin",
            "doxorubicin_ACY1083",
        }:
            continue

        raw = np.vstack(
            [
                by_arm["control"],
                by_arm["doxorubicin"],
                by_arm["doxorubicin_ACY1083"],
            ]
        )
        keep = (raw >= MIN_COUNT).sum(axis=0) >= 1
        cpm = raw[:, keep] / raw[:, keep].sum(axis=1, keepdims=True) * 1e6
        lv = np.log2(cpm + 1)
        genes = np.asarray(dox.var_names.astype(str))[keep]
        for j, gene in enumerate(genes):
            desc_rows.append(
                {
                    "cell_type": cell_type,
                    "gene": gene,
                    "doxorubicin_vs_control_log2CPM_difference": float(
                        lv[1, j] - lv[0, j]
                    ),
                    "ACY1083_vs_doxorubicin_log2CPM_difference": float(
                        lv[2, j] - lv[1, j]
                    ),
                    "p_value": np.nan,
                    "padj": np.nan,
                    "inference": "descriptive_only_one_pooled_library_per_arm",
                }
            )
    pd.DataFrame(desc_rows).to_csv(
        T / "pseudobulk_dox_descriptive.csv.gz", index=False
    )

    design = {
        "software": "PyDESeq2 0.5.4",
        "primary_study": "GSE216146",
        "unit_of_replication": "deposited pseudobulk library",
        "design": "~ replicate + condition",
        "replicate_blocks": {
            "replicate_1": "D20",
            "replicate_2": "D21 replicate 2",
            "replicate_3": "D21 replicate 3",
        },
        "minimum_cells_per_library_cell_type": MIN_CELLS,
        "gene_prefilter": (
            f">={MIN_COUNT} counts in >={MIN_LIBRARIES_WITH_COUNT} of 6 libraries"
        ),
        "multiple_testing": "PyDESeq2 independent filtering + BH adjusted P values",
        "significance_threshold": ALPHA,
        "welch_sensitivity": "unpaired Welch t test on log2(CPM+1), same six libraries and genes",
        "GSE271055": "descriptive only; one pooled deposited library per arm",
    }
    (T / "pseudobulk_deseq2_design.json").write_text(
        json.dumps(design, indent=2) + "\n"
    )

    print("\n===== PSEUDOBULK DESEQ2 SUMMARY =====")
    print(summary.to_string(index=False))
    print("\n===== ELIGIBILITY =====")
    print(
        eligibility.groupby(["contrast", "eligible"]).size().to_string()
    )


if __name__ == "__main__":
    main()
