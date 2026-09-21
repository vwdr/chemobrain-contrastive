"""Task 6: external disentanglement comparison with scDisInFact.

scDisInFact is an independently developed conditional VAE for multi-batch,
multi-condition single-cell data. It learns a shared-bio factor intended to be
condition/batch invariant and an unshared-bio factor intended to encode the
specified condition effect.

This benchmark uses the exact corrected canonical cell split and committed
1,500-gene universe. scDisInFact is trained only on the 6,000 canonical
training cells with:
    condition_key = ["drug"]  (control, cisplatin, doxorubicin)
    batch_key = "study"       (GSE216146, GSE271055)

The condition/study structure is therefore the same confounded structure as in
the manuscript benchmark: cisplatin occurs only in GSE216146 and doxorubicin
only in GSE271055, while control occurs in both. This comparison must not be
interpreted as resolving that confounding.

We use the official one-condition robustness configuration from the upstream
repository: Ks=[8,2], batch_size=64, lr=5e-4, 50 epochs and the package-default
regularization weights. Upstream is pinned in requirements-task6.txt.

Evaluation:
  * reconstruction on the canonical held-out cells, transformed back to the
    manuscript's log1p(1e4-normalized) 1,500-gene target;
  * linear-probe balanced accuracy for drug condition and study from the
    shared-bio, unshared-bio and combined latent factors;
  * linear CKA across seeds for latent stability;
  * cross-seed stability of condition-associated gene scores.

The comparison is descriptive because latent semantics differ from
MC-ContrastiveVI: scDisInFact's shared-bio factor is condition-irrelevant,
whereas MC-ContrastiveVI's "shared" treatment block is intended to represent
treatment-associated signal common to drugs.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import train_test_split

import scDisInFact.model as scd

R = Path(__file__).resolve().parents[1]
OUT = R / "analysis" / "corrected_20260920"
RUN = R / "runs" / "corrected_20260920"
OUT.mkdir(parents=True, exist_ok=True)
RUN.mkdir(parents=True, exist_ok=True)

SPLIT_SEED = 1729
SEEDS = (0, 1, 2)
N_EPOCHS = 50
K_SHARED = 8
K_CONDITION = 2
BATCH_SIZE = 64
LR = 5e-4
UPSTREAM_COMMIT = "9466ecf257d923f879cb39712c530d24826d88bc"

torch.set_num_threads(2)


def canonical_split(obs: pd.DataFrame):
    pure = np.where(~obs.rescue.astype(bool).to_numpy())[0]
    strata = (
        obs.study.astype(str) + "_" + obs.drug.astype(str)
    ).to_numpy()

    tr, te = train_test_split(
        pure,
        test_size=0.2,
        stratify=strata[pure],
        random_state=SPLIT_SEED,
    )
    va, te = train_test_split(
        te,
        test_size=0.5,
        stratify=strata[te],
        random_state=SPLIT_SEED,
    )

    rng = np.random.default_rng(SPLIT_SEED)
    tr = np.concatenate(
        [
            rng.choice(
                tr[strata[tr] == s],
                min(1500, int(np.sum(strata[tr] == s))),
                replace=False,
            )
            for s in np.unique(strata[tr])
        ]
    )
    return tr.astype(np.int64), va.astype(np.int64), te.astype(np.int64)


def verify_split(obs: pd.DataFrame, tr, va, te):
    split = np.full(len(obs), "not_used", dtype=object)
    split[tr] = "train"
    split[va] = "validation"
    split[te] = "test"
    split[obs.rescue.astype(bool).to_numpy()] = "rescue_projection"

    current = (
        obs.assign(split=split)
        .groupby(["study", "arm", "sample_id", "split"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["study", "arm", "sample_id", "split"])
        .reset_index(drop=True)
    )
    expected = (
        pd.read_csv(OUT / "benchmark_split_counts.csv")
        .sort_values(["study", "arm", "sample_id", "split"])
        .reset_index(drop=True)
    )
    # Pandas preserves AnnData categorical dtypes in `current`, whereas
    # read_csv returns string-backed columns in `expected`. Compare the
    # actual split labels/counts after normalizing representation dtypes.
    key_cols = ["study", "arm", "sample_id", "split"]
    current_cmp = current[expected.columns].copy()
    expected_cmp = expected.copy()
    for col in key_cols:
        current_cmp[col] = current_cmp[col].astype(str)
        expected_cmp[col] = expected_cmp[col].astype(str)
    current_cmp["n_cells"] = current_cmp["n_cells"].astype(np.int64)
    expected_cmp["n_cells"] = expected_cmp["n_cells"].astype(np.int64)

    pd.testing.assert_frame_equal(
        current_cmp,
        expected_cmp,
        check_dtype=True,
        check_categorical=False,
    )


def normalize_for_scdisinfact(counts: np.ndarray, scaler):
    total = counts.sum(axis=1, keepdims=True).astype(np.float64)
    if np.any(total <= 0):
        raise RuntimeError("Zero HVG count total encountered.")
    size_factor = total / 100.0
    norm = np.log1p(counts / size_factor)
    standardized = scaler.transform(norm).astype(np.float32)
    return standardized, size_factor.astype(np.float32)


def batch_ids_from_meta(data_dict, meta: pd.DataFrame):
    batch_names = np.asarray(
        data_dict["matching_dict"]["batch_name"]
    ).astype(str)
    lookup = {name: i for i, name in enumerate(batch_names)}
    vals = meta["study"].astype(str).to_numpy()
    missing = sorted(set(vals) - set(lookup))
    if missing:
        raise RuntimeError(f"Unknown scDisInFact batch labels: {missing}")
    return np.asarray([lookup[x] for x in vals], dtype=np.float32)


def infer_external(model, data_dict, counts: np.ndarray, meta: pd.DataFrame):
    x_std, size_factor = normalize_for_scdisinfact(
        counts, data_dict["scaler"]
    )
    batch_ids = batch_ids_from_meta(data_dict, meta)

    x_t = torch.from_numpy(x_std).to(model.device)
    b_t = torch.from_numpy(batch_ids[:, None]).to(model.device)

    with torch.no_grad():
        inf = model.inference(
            counts=x_t,
            batch_ids=b_t,
            print_stat=False,
        )
        gen = model.generative(
            z_c=inf["mu_c"],
            z_d=inf["mu_d"],
            batch_ids=b_t,
        )

    common = inf["mu_c"].detach().cpu().numpy()
    unshared = inf["mu_d"][0].detach().cpu().numpy()
    predicted_counts = (
        gen["mu"].detach().cpu().numpy() * size_factor
    )
    return common, unshared, predicted_counts


def probe(train_rep, test_rep, y_train, y_test):
    clf = LogisticRegression(
        max_iter=1000,
        class_weight="balanced",
    )
    clf.fit(train_rep, y_train)
    return float(
        balanced_accuracy_score(y_test, clf.predict(test_rep))
    )


def linear_cka(x: np.ndarray, y: np.ndarray) -> float:
    x = x - x.mean(axis=0, keepdims=True)
    y = y - y.mean(axis=0, keepdims=True)
    xy = np.linalg.norm(x.T @ y, ord="fro") ** 2
    xx = np.linalg.norm(x.T @ x, ord="fro")
    yy = np.linalg.norm(y.T @ y, ord="fro")
    denom = xx * yy
    return float(xy / denom) if denom > 0 else np.nan


def target_log_expression(
    raw_hvg: np.ndarray,
    full_library: np.ndarray,
):
    return np.log1p(
        raw_hvg / full_library[:, None] * 1e4
    ).astype(np.float32)


def build_mc_reference():
    rows = []
    full_metrics = []
    for seed in SEEDS:
        p = RUN / f"full_{seed}_metrics.json"
        if not p.exists():
            raise FileNotFoundError(p)
        m = json.loads(p.read_text())
        full_metrics.append(m)

    def add(metric, values, note):
        rows.append(
            {
                "method": "MC-ContrastiveVI",
                "representation": "canonical",
                "metric": metric,
                "mean": float(np.mean(values)),
                "sd_across_seeds": float(
                    np.std(values, ddof=1)
                ) if len(values) > 1 else np.nan,
                "n_seeds": len(values),
                "note": note,
            }
        )

    add(
        "test_log_expression_mse",
        [x["test_mse"] for x in full_metrics],
        "Archived corrected full-model checkpoints; canonical held-out cells.",
    )
    add(
        "condition_balanced_accuracy",
        [x["bg_condition_balanced_accuracy"] for x in full_metrics],
        "Background factor; intended treatment-irrelevant component.",
    )
    add(
        "study_balanced_accuracy",
        [x["bg_study_balanced_accuracy"] for x in full_metrics],
        "Background factor; study leakage diagnostic.",
    )
    add(
        "raw_shared_condition_balanced_accuracy",
        [
            x["shared_ungated_condition_balanced_accuracy"]
            for x in full_metrics
        ],
        "Raw treatment-shared posterior means; semantics differ from scDisInFact shared-bio.",
    )
    add(
        "raw_drug_condition_balanced_accuracy",
        [
            x["drug_ungated_condition_balanced_accuracy"]
            for x in full_metrics
        ],
        "Raw drug-head posterior means before label gating.",
    )

    pca_path = RUN / "pca_metrics.json"
    if pca_path.exists():
        pca = json.loads(pca_path.read_text())
        rows.append(
            {
                "method": "PCA(32)",
                "representation": "canonical",
                "metric": "test_log_expression_mse",
                "mean": float(pca["test_mse"]),
                "sd_across_seeds": np.nan,
                "n_seeds": 1,
                "note": "Canonical PCA reconstruction baseline.",
            }
        )
    return rows


def main():
    raw = ad.read_h5ad(
        R / "data" / "processed" / "recovered_counts.h5ad"
    )
    obs = raw.obs.copy()
    train, val, test = canonical_split(obs)
    verify_split(obs, train, val, test)

    genes = (
        pd.read_csv(OUT / "benchmark_gene_universe.csv")["gene"]
        .astype(str)
        .tolist()
    )
    if len(genes) != 1500:
        raise RuntimeError(f"Expected 1500 genes, got {len(genes)}")

    counts = raw[:, genes].X.toarray().astype(np.float32)
    full_library = np.asarray(raw.X.sum(axis=1)).ravel().astype(
        np.float32
    )

    meta = pd.DataFrame(
        {
            "study": obs.study.astype(str).to_numpy(),
            "drug": obs.drug.astype(str).to_numpy(),
        },
        index=obs.index.astype(str),
    )
    meta_train = meta.iloc[train].reset_index(drop=True)
    meta_test = meta.iloc[test].reset_index(drop=True)

    y_condition_train = meta_train["drug"].to_numpy()
    y_condition_test = meta_test["drug"].to_numpy()
    y_study_train = meta_train["study"].to_numpy()
    y_study_test = meta_test["study"].to_numpy()

    target_test = target_log_expression(
        counts[test],
        full_library[test],
    )

    metric_rows = []
    common_test_by_seed = {}
    unshared_test_by_seed = {}
    gene_scores = {}

    for seed in SEEDS:
        print(f"\n===== scDisInFact seed {seed} =====", flush=True)

        data_dict = scd.create_scdisinfact_dataset(
            counts[train],
            meta_train.rename(
                columns={"study": "batch", "drug": "condition"}
            ),
            condition_key=["condition"],
            batch_key="batch",
        )

        model = scd.scdisinfact(
            data_dict=data_dict,
            Ks=[K_SHARED, K_CONDITION],
            batch_size=BATCH_SIZE,
            interval=10,
            lr=LR,
            reg_mmd_comm=1e-4,
            reg_mmd_diff=1e-4,
            reg_gl=1,
            reg_class=1,
            reg_kl_comm=1e-5,
            reg_kl_diff=1e-2,
            seed=seed,
            device=torch.device("cpu"),
        )
        model.train()
        losses = model.train_model(
            nepochs=N_EPOCHS,
            recon_loss="NB",
        )
        model.eval()

        train_meta_api = meta_train.rename(
            columns={"study": "batch", "drug": "condition"}
        )
        test_meta_api = meta_test.rename(
            columns={"study": "batch", "drug": "condition"}
        )

        common_train, unshared_train, _ = infer_external(
            model,
            data_dict,
            counts[train],
            train_meta_api.rename(columns={"batch": "study"}),
        )
        common_test, unshared_test, pred_counts_test = infer_external(
            model,
            data_dict,
            counts[test],
            test_meta_api.rename(columns={"batch": "study"}),
        )

        # infer_external expects a column named "study"; restoring the
        # original study label above keeps mapping explicit.
        combined_train = np.concatenate(
            [common_train, unshared_train], axis=1
        )
        combined_test = np.concatenate(
            [common_test, unshared_test], axis=1
        )

        pred_log = target_log_expression(
            np.clip(pred_counts_test, 0, None),
            full_library[test],
        )
        test_mse = float(np.mean((pred_log - target_test) ** 2))

        reps = {
            "shared_bio": (common_train, common_test),
            "unshared_bio": (unshared_train, unshared_test),
            "combined": (combined_train, combined_test),
        }
        row = {
            "seed": seed,
            "test_log_expression_mse": test_mse,
        }
        for name, (rtr, rte) in reps.items():
            row[
                f"{name}_condition_balanced_accuracy"
            ] = probe(
                rtr,
                rte,
                y_condition_train,
                y_condition_test,
            )
            row[
                f"{name}_study_balanced_accuracy"
            ] = probe(
                rtr,
                rte,
                y_study_train,
                y_study_test,
            )

        scores = model.extract_gene_scores()[0]
        if len(scores) != len(genes):
            raise RuntimeError("Gene-score length mismatch.")
        gene_scores[seed] = np.asarray(scores, dtype=float)

        common_test_by_seed[seed] = common_test
        unshared_test_by_seed[seed] = unshared_test
        metric_rows.append(row)

        pd.DataFrame(
            {
                "seed": seed,
                "gene": genes,
                "condition_gene_score": scores,
            }
        ).to_csv(
            RUN / f"task6_scdisinfact_gene_scores_seed{seed}.csv",
            index=False,
        )

        torch.save(
            model.state_dict(),
            RUN / f"task6_scdisinfact_seed{seed}.pt",
        )
        np.savez_compressed(
            RUN / f"task6_scdisinfact_seed{seed}_test_latents.npz",
            common=common_test,
            unshared=unshared_test,
        )

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(
        OUT / "task6_scdisinfact_metrics.csv", index=False
    )

    summary_rows = []
    for col in metrics.columns:
        if col == "seed":
            continue
        summary_rows.append(
            {
                "metric": col,
                "mean": float(metrics[col].mean()),
                "sd_across_seeds": float(metrics[col].std(ddof=1)),
                "min": float(metrics[col].min()),
                "max": float(metrics[col].max()),
                "n_seeds": len(metrics),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        OUT / "task6_scdisinfact_summary.csv", index=False
    )

    cka_rows = []
    score_stability_rows = []
    for a, b in combinations(SEEDS, 2):
        cka_rows.extend(
            [
                {
                    "seed_a": a,
                    "seed_b": b,
                    "representation": "shared_bio",
                    "linear_cka": linear_cka(
                        common_test_by_seed[a],
                        common_test_by_seed[b],
                    ),
                },
                {
                    "seed_a": a,
                    "seed_b": b,
                    "representation": "unshared_bio",
                    "linear_cka": linear_cka(
                        unshared_test_by_seed[a],
                        unshared_test_by_seed[b],
                    ),
                },
            ]
        )

        sa = gene_scores[a]
        sb = gene_scores[b]
        rho = float(spearmanr(sa, sb).statistic)
        top_a = set(np.argsort(sa)[-100:])
        top_b = set(np.argsort(sb)[-100:])
        jaccard = len(top_a & top_b) / len(top_a | top_b)
        score_stability_rows.append(
            {
                "seed_a": a,
                "seed_b": b,
                "spearman_gene_scores": rho,
                "top100_jaccard": float(jaccard),
                "top100_overlap": int(len(top_a & top_b)),
            }
        )

    pd.DataFrame(cka_rows).to_csv(
        OUT / "task6_scdisinfact_seed_cka.csv", index=False
    )
    pd.DataFrame(score_stability_rows).to_csv(
        OUT / "task6_scdisinfact_gene_score_stability.csv",
        index=False,
    )

    gene_score_table = pd.DataFrame({"gene": genes})
    for seed in SEEDS:
        gene_score_table[f"score_seed{seed}"] = gene_scores[seed]
    gene_score_table["mean_score"] = gene_score_table[
        [f"score_seed{x}" for x in SEEDS]
    ].mean(axis=1)
    gene_score_table["sd_score"] = gene_score_table[
        [f"score_seed{x}" for x in SEEDS]
    ].std(axis=1, ddof=1)
    gene_score_table.sort_values(
        "mean_score", ascending=False
    ).to_csv(
        OUT / "task6_scdisinfact_gene_scores.csv.gz",
        index=False,
        compression="gzip",
    )

    comparison = build_mc_reference()
    for _, x in summary.iterrows():
        comparison.append(
            {
                "method": "scDisInFact",
                "representation": (
                    x["metric"].split("_condition_balanced_accuracy")[0]
                    if "_condition_balanced_accuracy" in x["metric"]
                    else (
                        x["metric"].split("_study_balanced_accuracy")[0]
                        if "_study_balanced_accuracy" in x["metric"]
                        else "canonical"
                    )
                ),
                "metric": x["metric"],
                "mean": float(x["mean"]),
                "sd_across_seeds": float(x["sd_across_seeds"]),
                "n_seeds": int(x["n_seeds"]),
                "note": (
                    "Official upstream one-condition robustness configuration; "
                    "same canonical train/test cells and 1,500 genes."
                ),
            }
        )
    pd.DataFrame(comparison).to_csv(
        OUT / "task6_external_method_comparison.csv",
        index=False,
    )

    design = {
        "external_method": "scDisInFact",
        "upstream_repository": "ZhangLabGT/scDisInFact",
        "upstream_commit": UPSTREAM_COMMIT,
        "upstream_publication": (
            "Zhang et al., Nature Communications 15, 912 (2024)"
        ),
        "canonical_split_seed": SPLIT_SEED,
        "training_cells": int(len(train)),
        "validation_cells_not_used_by_scdisinfact": int(len(val)),
        "test_cells": int(len(test)),
        "genes": len(genes),
        "condition_key": "drug",
        "batch_key": "study",
        "condition_levels": sorted(meta_train.drug.unique().tolist()),
        "batch_levels": sorted(meta_train.study.unique().tolist()),
        "scdisinfact_Ks": [K_SHARED, K_CONDITION],
        "scdisinfact_epochs": N_EPOCHS,
        "scdisinfact_batch_size": BATCH_SIZE,
        "scdisinfact_lr": LR,
        "scdisinfact_seeds": list(SEEDS),
        "reconstruction_evaluation": (
            "Predicted raw HVG counts were transformed using each cell's "
            "original full-library size into the same log1p(1e4-normalized) "
            "gene-expression target used by the canonical benchmark."
        ),
        "interpretation_limitations": [
            (
                "Drug identity is partially confounded with study: cisplatin "
                "occurs only in GSE216146 and doxorubicin only in GSE271055."
            ),
            (
                "scDisInFact shared-bio is defined as condition-irrelevant; "
                "it is not semantically equivalent to MC-ContrastiveVI's "
                "treatment-shared latent block."
            ),
            (
                "The comparison is therefore descriptive and does not provide "
                "an external validation of the MC shared-fraction statistic."
            ),
        ],
    }
    (OUT / "task6_external_method_design.json").write_text(
        json.dumps(design, indent=2) + "\n"
    )

    print("\n===== TASK 6 SCDISINFACT METRICS =====")
    print(metrics.to_string(index=False))
    print("\n===== TASK 6 SCDISINFACT SUMMARY =====")
    print(summary.to_string(index=False))
    print("\n===== TASK 6 LATENT CKA =====")
    print(pd.DataFrame(cka_rows).to_string(index=False))
    print("\n===== TASK 6 GENE-SCORE STABILITY =====")
    print(pd.DataFrame(score_stability_rows).to_string(index=False))


if __name__ == "__main__":
    main()