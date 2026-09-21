"""Matched peer-review baselines using scvi-tools.

This analysis adds two widely used comparators that are closer to the manuscript
than PCA alone:

1. scVI on the exact canonical 6,000-cell training set and committed 1,500-gene
   universe, with study as the batch covariate.
2. Pairwise contrastiveVI within each source study, comparing control versus the
   study's chemotherapy condition. This avoids pretending that pairwise
   contrastiveVI is a multi-drug shared/specific model.

All models are run for seeds 0-2 with fixed architecture-scale choices and no
dataset-specific hyperparameter search. These comparisons are sensitivity
analyses rather than claims of exhaustive benchmarking.

For reconstruction, scvi-tools normalized expression is evaluated on a
modeled-gene scale: both predictions and observed counts are normalized to
10,000 within the 1,500 modeled genes and log1p transformed. This metric is
reported separately from the manuscript's canonical full-library normalization
MSE and must not be numerically conflated with it.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import train_test_split

R = Path(__file__).resolve().parents[1]
OUT = R / "analysis" / "corrected_20260920"
OUT.mkdir(parents=True, exist_ok=True)

SPLIT_SEED = 1729
SEEDS = (0, 1, 2)
EPOCHS = 100
N_GENES = 1500


def canonical_split(obs: pd.DataFrame):
    pure = np.where(~obs.rescue.astype(bool).to_numpy())[0]
    strata = (obs.study.astype(str) + "_" + obs.drug.astype(str)).to_numpy()
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


def probe(x_train, x_test, y_train, y_test):
    clf = LogisticRegression(max_iter=1500, class_weight="balanced")
    clf.fit(x_train, y_train)
    return float(balanced_accuracy_score(y_test, clf.predict(x_test)))


def linear_cka(x, y):
    x = x - x.mean(0, keepdims=True)
    y = y - y.mean(0, keepdims=True)
    num = np.linalg.norm(x.T @ y, "fro") ** 2
    den = np.linalg.norm(x.T @ x, "fro") * np.linalg.norm(y.T @ y, "fro")
    return float(num / den) if den > 0 else np.nan


def modeled_gene_log(counts):
    counts = np.asarray(counts, dtype=np.float64)
    total = counts.sum(1, keepdims=True)
    total[total <= 0] = 1.0
    return np.log1p(counts / total * 1e4).astype(np.float32)


def make_adata(raw, indices, genes):
    x = raw[indices, genes].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    x = np.asarray(np.rint(x), dtype=np.int64)
    obs = raw.obs.iloc[indices].copy()
    a = ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=pd.Index(genes)))
    a.obs_names_make_unique()
    return a


def main():
    raw = ad.read_h5ad(R / "data" / "processed" / "recovered_counts.h5ad")
    train, val, test = canonical_split(raw.obs)
    genes = pd.read_csv(OUT / "benchmark_gene_universe.csv")["gene"].astype(str).tolist()
    assert len(genes) == N_GENES

    # ---------- scVI ----------
    scvi_rows = []
    scvi_test_latents = {}
    for seed in SEEDS:
        print(f"===== scVI seed {seed} =====", flush=True)
        scvi.settings.seed = seed
        torch.manual_seed(seed)
        np.random.seed(seed)

        atr = make_adata(raw, train, genes)
        ate = make_adata(raw, test, genes)
        scvi.model.SCVI.setup_anndata(atr, batch_key="study")
        model = scvi.model.SCVI(
            atr,
            n_hidden=128,
            n_latent=32,
            n_layers=2,
            dropout_rate=0.1,
            dispersion="gene",
            gene_likelihood="nb",
            use_observed_lib_size=True,
        )
        model.train(
            max_epochs=EPOCHS,
            train_size=1.0,
            validation_size=0.0,
            early_stopping=False,
            accelerator="cpu",
            devices=1,
            batch_size=256,
            enable_progress_bar=False,
        )

        ztr = model.get_latent_representation(atr)
        zte = model.get_latent_representation(ate)
        pred = model.get_normalized_expression(
            ate,
            library_size=1e4,
            n_samples=1,
            return_mean=True,
            return_numpy=True,
        )
        obs_log = modeled_gene_log(ate.X)
        pred_log = np.log1p(np.asarray(pred, dtype=np.float64)).astype(np.float32)
        mse = float(np.mean((pred_log - obs_log) ** 2))

        row = {
            "method": "scVI",
            "study_or_pair": "combined",
            "seed": seed,
            "modeled_gene_log_mse": mse,
            "condition_balanced_accuracy": probe(
                ztr, zte,
                atr.obs["drug"].astype(str).to_numpy(),
                ate.obs["drug"].astype(str).to_numpy(),
            ),
            "study_balanced_accuracy": probe(
                ztr, zte,
                atr.obs["study"].astype(str).to_numpy(),
                ate.obs["study"].astype(str).to_numpy(),
            ),
        }
        scvi_rows.append(row)
        scvi_test_latents[seed] = zte

    # ---------- pairwise contrastiveVI ----------
    pair_specs = {
        "cisplatin": ("GSE216146", "cisplatin"),
        "doxorubicin": ("GSE271055", "doxorubicin"),
    }
    cvi_rows = []
    cvi_latents = {}
    for pair_name, (study, treated) in pair_specs.items():
        tr_mask = (
            (raw.obs.iloc[train].study.astype(str).to_numpy() == study)
            & np.isin(
                raw.obs.iloc[train].drug.astype(str).to_numpy(),
                ["control", treated],
            )
        )
        te_mask = (
            (raw.obs.iloc[test].study.astype(str).to_numpy() == study)
            & np.isin(
                raw.obs.iloc[test].drug.astype(str).to_numpy(),
                ["control", treated],
            )
        )
        pair_train = train[tr_mask]
        pair_test = test[te_mask]
        assert len(pair_train) == 3000

        for seed in SEEDS:
            print(f"===== contrastiveVI {pair_name} seed {seed} =====", flush=True)
            scvi.settings.seed = seed
            torch.manual_seed(seed)
            np.random.seed(seed)

            atr = make_adata(raw, pair_train, genes)
            ate = make_adata(raw, pair_test, genes)
            scvi.external.ContrastiveVI.setup_anndata(atr)

            background_idx = np.where(atr.obs["drug"].astype(str).to_numpy() == "control")[0]
            target_idx = np.where(atr.obs["drug"].astype(str).to_numpy() == treated)[0]

            model = scvi.external.ContrastiveVI(
                atr,
                n_hidden=128,
                n_background_latent=16,
                n_salient_latent=8,
                n_layers=2,
                dropout_rate=0.1,
                use_observed_lib_size=True,
                wasserstein_penalty=0,
            )
            model.train(
                background_indices=background_idx,
                target_indices=target_idx,
                max_epochs=EPOCHS,
                train_size=1.0,
                validation_size=0.0,
                early_stopping=False,
                accelerator="cpu",
                devices=1,
                batch_size=128,
                enable_progress_bar=False,
            )

            zbg_tr = model.get_latent_representation(atr, representation_kind="background")
            zsal_tr = model.get_latent_representation(atr, representation_kind="salient")
            zbg_te = model.get_latent_representation(ate, representation_kind="background")
            zsal_te = model.get_latent_representation(ate, representation_kind="salient")
            pred = model.get_normalized_expression(
                ate,
                library_size=1e4,
                n_samples=1,
                return_mean=True,
                return_numpy=True,
            )
            obs_log = modeled_gene_log(ate.X)
            pred_log = np.log1p(np.asarray(pred, dtype=np.float64)).astype(np.float32)
            mse = float(np.mean((pred_log - obs_log) ** 2))
            ytr = (atr.obs["drug"].astype(str).to_numpy() == treated).astype(int)
            yte = (ate.obs["drug"].astype(str).to_numpy() == treated).astype(int)

            cvi_rows.append(
                {
                    "method": "contrastiveVI",
                    "study_or_pair": pair_name,
                    "seed": seed,
                    "n_train": len(atr),
                    "n_test": len(ate),
                    "modeled_gene_log_mse": mse,
                    "background_condition_balanced_accuracy": probe(
                        zbg_tr, zbg_te, ytr, yte
                    ),
                    "salient_condition_balanced_accuracy": probe(
                        zsal_tr, zsal_te, ytr, yte
                    ),
                }
            )
            cvi_latents[(pair_name, seed, "background")] = zbg_te
            cvi_latents[(pair_name, seed, "salient")] = zsal_te

    scvi_df = pd.DataFrame(scvi_rows)
    cvi_df = pd.DataFrame(cvi_rows)
    scvi_df.to_csv(OUT / "peerreview_scvi_metrics.csv", index=False)
    cvi_df.to_csv(OUT / "peerreview_contrastivevi_metrics.csv", index=False)

    cka_rows = []
    for a, b in combinations(SEEDS, 2):
        cka_rows.append(
            {
                "method": "scVI",
                "study_or_pair": "combined",
                "representation": "latent",
                "seed_a": a,
                "seed_b": b,
                "linear_cka": linear_cka(scvi_test_latents[a], scvi_test_latents[b]),
            }
        )
        for pair_name in pair_specs:
            for rep in ("background", "salient"):
                cka_rows.append(
                    {
                        "method": "contrastiveVI",
                        "study_or_pair": pair_name,
                        "representation": rep,
                        "seed_a": a,
                        "seed_b": b,
                        "linear_cka": linear_cka(
                            cvi_latents[(pair_name, a, rep)],
                            cvi_latents[(pair_name, b, rep)],
                        ),
                    }
                )
    cka = pd.DataFrame(cka_rows)
    cka.to_csv(OUT / "peerreview_matched_baseline_cka.csv", index=False)

    summary = []
    for method, frame in [("scVI", scvi_df), ("contrastiveVI", cvi_df)]:
        group_cols = ["study_or_pair"]
        for key, g in frame.groupby(group_cols, dropna=False):
            if not isinstance(key, tuple):
                key = (key,)
            for col in frame.columns:
                if col in {"method", "study_or_pair", "seed", "n_train", "n_test"}:
                    continue
                if np.issubdtype(frame[col].dtype, np.number):
                    summary.append(
                        {
                            "method": method,
                            "study_or_pair": key[0],
                            "metric": col,
                            "mean": float(g[col].mean()),
                            "sd_across_seeds": float(g[col].std(ddof=1)),
                            "min": float(g[col].min()),
                            "max": float(g[col].max()),
                        }
                    )
    pd.DataFrame(summary).to_csv(
        OUT / "peerreview_matched_baseline_summary.csv", index=False
    )

    design = {
        "scvi_tools_version": scvi.__version__,
        "seeds": list(SEEDS),
        "epochs": EPOCHS,
        "canonical_training_cells": int(len(train)),
        "canonical_test_cells": int(len(test)),
        "modeled_genes": N_GENES,
        "scVI": {
            "n_latent": 32,
            "n_hidden": 128,
            "n_layers": 2,
            "dropout_rate": 0.1,
            "gene_likelihood": "nb",
            "batch_key": "study",
        },
        "contrastiveVI": {
            "pairwise_only": True,
            "pairs": pair_specs,
            "n_background_latent": 16,
            "n_salient_latent": 8,
            "n_hidden": 128,
            "n_layers": 2,
            "dropout_rate": 0.1,
        },
        "reconstruction_metric": (
            "log1p expression after normalization to 10,000 within the 1,500 "
            "modeled genes. This is intentionally reported separately from the "
            "manuscript canonical full-library normalization MSE."
        ),
        "interpretation": (
            "Sensitivity comparison only. No dataset-specific hyperparameter "
            "search was performed, and pairwise contrastiveVI does not provide "
            "a multi-drug shared-versus-specific decomposition."
        ),
    }
    (OUT / "peerreview_matched_baseline_design.json").write_text(
        json.dumps(design, indent=2) + "\n"
    )

    print("\n===== MATCHED BASELINE SUMMARY =====")
    print(pd.DataFrame(summary).to_string(index=False))
    print("\n===== MATCHED BASELINE CKA =====")
    print(cka.to_string(index=False))


if __name__ == "__main__":
    main()
