"""Task 7: gene-pattern and marker-set interpretation audit.

This analysis deliberately separates four questions that were previously mixed:

1. Which integrated-gradient genes are stable enough to recur among each
   seed's top 50?
2. Do attribution lists preferentially overlap source-cell-type marker sets?
3. Does the cisplatin-associated attribution list overlap independently
   summarized cisplatin-vs-control pseudobulk differential-expression genes?
4. Do broad pathway associations survive stricter seed-consensus definitions?

The fixed attribution lists in task7_attribution_lists.csv come from the
corrected zero-baseline, squared-posterior-norm integrated gradients already
audited for numerical completeness. No attribution values are recomputed here.

Marker analyses use two resources:
  * the small preregistered/canonical marker panel already used in the paper;
  * empirical source-label marker sets derived only from GSE216146 PBS control
    cells. Empirical markers are descriptive rather than independent
    validation. To reduce library-size effects, expression specificity is
    calculated within each of the three PBS libraries and then averaged
    equally across libraries.

For each source cell type represented by >=20 cells in all three PBS control
libraries, the empirical score for a model-universe gene is the mean across
libraries of:
    mean log-expression(cell type) - mean log-expression(all other cells)
Genes must have pooled fraction-positive >=0.10 in the focal cell type and a
positive specificity score. The top 50 become that cell type's empirical marker
set.

Pseudobulk support uses Task 4's paired PyDESeq2 cisplatin-vs-control results.
A supported DE row has adjusted P<0.05 and all three GEO replicate-pair logCPM
effects in the same direction as the DESeq2 log2 fold change. A gene is in the
DE-support set if at least one eligible source cell type has such a row.

All overrepresentation P values use the committed 1,500-gene model universe
and are BH adjusted within each axis/list definition across marker classes.
These tests describe overlap only; they do not establish cell-of-origin,
activation, causal regulation, or cross-drug convergence.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

R=Path(__file__).resolve().parents[1]
T=R/"analysis"/"corrected_20260920"
MODEL_GENES=pd.read_csv(T/"benchmark_gene_universe.csv")["gene"].astype(str).tolist()
UNIVERSE=set(MODEL_GENES)
ATTR=pd.read_csv(T/"task7_attribution_lists.csv")
CANON=pd.read_csv(T/"canonical_markers.csv")
PB_PATH=T/"pseudobulk_task4_paired_results.csv.gz"

CONTROL_LIBRARIES=["D20-6407","D21-2746","D21-2747"]
MIN_CELLS_PER_LIBRARY=20
EMPIRICAL_MARKER_COUNT=50
MIN_FRACTION_POSITIVE=0.10


def hypergeom_rows(axis,list_type,geneset,marker_sets,resource):
    rows=[]
    M=len(UNIVERSE); N=len(geneset)
    for marker_class, markers in sorted(marker_sets.items()):
        markers=set(markers)&UNIVERSE
        K=len(markers); overlap=geneset&markers; k=len(overlap)
        p=float(hypergeom.sf(k-1,M,K,N)) if k else 1.0
        expected=N*K/M if M else np.nan
        fold=float(k/expected) if expected>0 else np.nan
        rows.append({
            "axis":axis,
            "list_type":list_type,
            "list_size":N,
            "resource":resource,
            "marker_class":marker_class,
            "marker_set_size_in_universe":K,
            "overlap":k,
            "fold_enrichment":fold,
            "p_value":p,
            "genes":"|".join(sorted(overlap)),
        })
    q=multipletests([r["p_value"] for r in rows],method="fdr_bh")[1]
    for r,qq in zip(rows,q):
        r["q_value"]=float(qq)
    return rows


def empirical_marker_sets(a):
    obs=a.obs
    mask=(
        (obs.study.astype(str)=="GSE216146")
        & (obs.arm.astype(str)=="control")
        & obs.sample_id.astype(str).isin(CONTROL_LIBRARIES)
    ).to_numpy()
    x=a[mask,MODEL_GENES].copy()
    scounts=np.asarray(x.X.sum(axis=1)).ravel()
    # Match manuscript expression scale.
    mat=x.X.toarray().astype(np.float64)
    mat=np.log1p(mat/np.maximum(scounts[:,None],1.0)*1e4)
    meta=x.obs.copy().reset_index(drop=True)
    cell_types=sorted(meta.source_cell_type.astype(str).unique())

    eligible=[]
    qc=[]
    for ct in cell_types:
        counts=[]
        for sid in CONTROL_LIBRARIES:
            n=int(((meta.source_cell_type.astype(str)==ct)&(meta.sample_id.astype(str)==sid)).sum())
            counts.append(n)
        ok=min(counts)>=MIN_CELLS_PER_LIBRARY
        qc.append({
            "cell_type":ct,
            **{f"{sid}_cells":n for sid,n in zip(CONTROL_LIBRARIES,counts)},
            "eligible_all_three_libraries":ok,
        })
        if ok: eligible.append(ct)

    score_rows=[]
    marker_sets={}
    ct_arr=meta.source_cell_type.astype(str).to_numpy()
    sid_arr=meta.sample_id.astype(str).to_numpy()
    for ct in eligible:
        diffs=[]
        for sid in CONTROL_LIBRARIES:
            focal=(ct_arr==ct)&(sid_arr==sid)
            other=(ct_arr!=ct)&(sid_arr==sid)
            if focal.sum()<MIN_CELLS_PER_LIBRARY or other.sum()==0:
                raise RuntimeError((ct,sid,focal.sum(),other.sum()))
            diffs.append(mat[focal].mean(axis=0)-mat[other].mean(axis=0))
        specificity=np.vstack(diffs).mean(axis=0)
        focal_all=ct_arr==ct
        fraction_positive=(mat[focal_all]>0).mean(axis=0)
        df=pd.DataFrame({
            "cell_type":ct,
            "gene":MODEL_GENES,
            "specificity_score":specificity,
            "fraction_positive":fraction_positive,
        })
        df["eligible_marker"]=(
            (df.specificity_score>0)
            & (df.fraction_positive>=MIN_FRACTION_POSITIVE)
        )
        rank=df[df.eligible_marker].sort_values(
            ["specificity_score","fraction_positive"],
            ascending=[False,False],
        ).head(EMPIRICAL_MARKER_COUNT).copy()
        rank["empirical_marker_rank"]=np.arange(1,len(rank)+1)
        marker_sets[ct]=set(rank.gene)
        score_rows.append(rank)

    return marker_sets,pd.concat(score_rows,ignore_index=True),pd.DataFrame(qc)


def load_pseudobulk_support():
    if not PB_PATH.exists():
        raise FileNotFoundError(
            f"{PB_PATH} missing; run scripts/16_pseudobulk_paired.py first"
        )
    pb=pd.read_csv(PB_PATH)
    q=pb[pb.contrast=="cisplatin_vs_control"].copy()
    pair_cols=[
        "replicate_1_log2cpm_difference",
        "replicate_2_log2cpm_difference",
        "replicate_3_log2cpm_difference",
    ]
    signs=np.sign(q[pair_cols].to_numpy())
    model_sign=np.sign(q.log2FoldChange.to_numpy())[:,None]
    q["all_three_pairs_same_direction"]=(signs==model_sign).all(axis=1)
    q["supported_row"]=(
        (q.wald_padj<0.05)
        & q.all_three_pairs_same_direction
    )

    gene=q.groupby("gene",as_index=False).agg(
        n_eligible_cell_types=("cell_type","nunique"),
        n_fdr005_cell_types=("wald_padj",lambda x:int((x<0.05).sum())),
        n_direction_consistent_fdr005_cell_types=(
            "supported_row",lambda x:int(np.sum(x))
        ),
        minimum_padj=("wald_padj","min"),
        max_abs_log2FoldChange=("log2FoldChange",lambda x:float(np.max(np.abs(x)))),
    )
    gene["supported_de_gene"]=gene.n_direction_consistent_fdr005_cell_types>0
    support=set(gene.loc[gene.supported_de_gene,"gene"])&UNIVERSE
    return support,gene


def main():
    a=ad.read_h5ad(R/"data"/"processed"/"recovered_counts.h5ad")
    empirical_sets,empirical_table,empirical_qc=empirical_marker_sets(a)
    empirical_table.to_csv(T/"task7_empirical_source_markers.csv",index=False)
    empirical_qc.to_csv(T/"task7_empirical_source_marker_qc.csv",index=False)

    canonical_sets={
        ct:set(g.gene.astype(str))&UNIVERSE
        for ct,g in CANON.groupby("cell_type")
    }

    supported_de,de_table=load_pseudobulk_support()
    de_table.to_csv(T/"task7_cis_pseudobulk_gene_support.csv",index=False)

    enr=[]
    overlap_rows=[]
    gene_pattern_rows=[]

    for (axis,list_type),g in ATTR.groupby(["axis","list_type"],sort=True):
        geneset=set(g.gene.astype(str))&UNIVERSE
        enr.extend(hypergeom_rows(
            axis,list_type,geneset,canonical_sets,"canonical_marker_panel"
        ))
        enr.extend(hypergeom_rows(
            axis,list_type,geneset,empirical_sets,
            "GSE216146_PBS_empirical_markers"
        ))

        k=len(geneset&supported_de)
        p=float(hypergeom.sf(
            k-1,len(UNIVERSE),len(supported_de),len(geneset)
        )) if k else 1.0
        overlap_rows.append({
            "axis":axis,
            "list_type":list_type,
            "list_size":len(geneset),
            "cis_supported_de_universe_size":len(supported_de),
            "overlap":k,
            "fold_enrichment":(
                k/(len(geneset)*len(supported_de)/len(UNIVERSE))
                if supported_de else np.nan
            ),
            "p_value":p,
            "genes":"|".join(sorted(geneset&supported_de)),
            "interpretation":(
                "direct cisplatin attribution-vs-DE check"
                if axis=="drug1"
                else "context only; not a matched drug-specific validation"
            ),
        })

        for _,r in g.iterrows():
            gene=str(r.gene)
            canon_classes=[
                ct for ct,s in canonical_sets.items() if gene in s
            ]
            empirical_classes=[
                ct for ct,s in empirical_sets.items() if gene in s
            ]
            drow=de_table[de_table.gene==gene]
            gene_pattern_rows.append({
                "axis":axis,
                "list_type":list_type,
                "rank":int(r["rank"]),
                "gene":gene,
                "top50_seed_support":int(r.top50_seed_support),
                "across_seed_mean_abs_ig":float(r.across_seed_mean_abs_ig),
                "canonical_marker_classes":"|".join(canon_classes),
                "empirical_PBS_marker_classes":"|".join(empirical_classes),
                "cis_supported_de_gene":(
                    bool(drow.supported_de_gene.iloc[0])
                    if len(drow) else False
                ),
                "cis_n_supported_cell_types":(
                    int(drow.n_direction_consistent_fdr005_cell_types.iloc[0])
                    if len(drow) else 0
                ),
                "cis_minimum_padj":(
                    float(drow.minimum_padj.iloc[0])
                    if len(drow) else np.nan
                ),
            })

    enrichment=pd.DataFrame(enr)
    enrichment.to_csv(T/"task7_marker_set_enrichment.csv",index=False)

    overlap=pd.DataFrame(overlap_rows)
    # BH separately for the nine attribution-list definitions.
    overlap["q_value_bh"]=multipletests(
        overlap.p_value,method="fdr_bh"
    )[1]
    overlap.to_csv(T/"task7_cis_de_overlap.csv",index=False)
    pd.DataFrame(gene_pattern_rows).to_csv(
        T/"task7_gene_pattern_table.csv",index=False
    )

    stable_summary=[]
    for axis in ["shared","drug0","drug1"]:
        vals={}
        for lt in [
            "consensus_mean_top50",
            "top50_in_at_least_2_seeds",
            "top50_in_all_3_seeds",
        ]:
            vals[lt]=int(
                ((ATTR.axis==axis)&(ATTR.list_type==lt)).sum()
            )
        q=enrichment[
            (enrichment.axis==axis)
            & (enrichment.q_value<0.05)
        ]
        for lt in vals:
            qq=q[q.list_type==lt]
            stable_summary.append({
                "axis":axis,
                "list_type":lt,
                "list_size":vals[lt],
                "n_significant_canonical_marker_classes":int(
                    ((qq.resource=="canonical_marker_panel")).sum()
                ),
                "n_significant_empirical_marker_classes":int(
                    ((qq.resource=="GSE216146_PBS_empirical_markers")).sum()
                ),
                "significant_empirical_classes":"|".join(
                    qq.loc[
                        qq.resource=="GSE216146_PBS_empirical_markers",
                        "marker_class"
                    ].astype(str).tolist()
                ),
                "cis_DE_overlap":int(
                    overlap.loc[
                        (overlap.axis==axis)&(overlap.list_type==lt),
                        "overlap"
                    ].iloc[0]
                ),
                "cis_DE_overlap_q":float(
                    overlap.loc[
                        (overlap.axis==axis)&(overlap.list_type==lt),
                        "q_value_bh"
                    ].iloc[0]
                ),
            })
    summary=pd.DataFrame(stable_summary)
    summary.to_csv(T/"task7_gene_pattern_summary.csv",index=False)

    design={
        "attribution_target":"squared posterior-mean norm",
        "attribution_baseline":"zero",
        "list_definitions":{
            "consensus_mean_top50":"top 50 by across-seed mean absolute IG",
            "top50_in_at_least_2_seeds":"union genes present in seed top 50 for >=2/3 seeds",
            "top50_in_all_3_seeds":"intersection of all three seed top-50 lists",
        },
        "model_gene_universe":len(UNIVERSE),
        "empirical_marker_source":"GSE216146 PBS control only",
        "empirical_marker_library_weighting":"equal across D20-6407, D21-2746, D21-2747",
        "empirical_marker_min_cells_per_library":MIN_CELLS_PER_LIBRARY,
        "empirical_marker_min_fraction_positive":MIN_FRACTION_POSITIVE,
        "empirical_marker_top_n":EMPIRICAL_MARKER_COUNT,
        "cis_DE_support":"Task4 paired PyDESeq2 FDR<0.05 and all 3 paired logCPM effects agree with model direction",
        "limitations":[
            "Empirical source markers are descriptive and derived from the same GSE216146 study, so they are not independent validation.",
            "Canonical marker panel is intentionally small and therefore has limited enrichment power.",
            "Cisplatin pseudobulk support is available only for GSE216146; the pooled doxorubicin study cannot support replicated gene-level inference.",
            "Marker/pathway overlap does not identify cell of origin or establish activation, causality, or cross-drug convergence."
        ],
    }
    (T/"task7_gene_pattern_design.json").write_text(
        json.dumps(design,indent=2)+"\n"
    )

    print("===== TASK 7 GENE PATTERN SUMMARY =====")
    print(summary.to_string(index=False))
    print("\n===== TASK 7 CIS DE OVERLAP =====")
    print(overlap.to_string(index=False))
    print("\n===== TASK 7 SIGNIFICANT MARKER ENRICHMENTS =====")
    sig=enrichment[enrichment.q_value<0.05].sort_values(
        ["axis","list_type","resource","q_value"]
    )
    print(sig.to_string(index=False) if len(sig) else "None")


if __name__=="__main__":
    main()
