"""Task 5: cisplatin unseen-library validation.

The canonical benchmark holds out cells from libraries that also contribute
training cells. This sensitivity analysis instead holds out an entire
GSE216146 control/cisplatin GEO replicate pair at a time.

For each of the three replicate-pair folds:
  * rescue libraries are excluded;
  * the two libraries in the held-out pair contribute no training, validation,
    or feature-selection cells;
  * the remaining libraries contribute a balanced training set of 1,500 cells
    per study-by-treatment stratum (6,000 total) and 100 validation cells per
    stratum (400 total);
  * 1,500 HVGs are selected from the training cells only;
  * all retained cells from the held-out control/cisplatin pair form the test
    set.

The full MC-ContrastiveVI and the noncontrastive Gaussian-VAE sensitivity model
are fit for seeds 0, 1, and 2. PCA(32) and a training-mean predictor are included
as deterministic reconstruction baselines.

The primary metric is held-out-library balanced MSE: the mean of the control-
library MSE and cisplatin-library MSE within each fold. Results are also reported
for the seen-library GSE216146 validation cells so the within-study
generalization gap can be inspected.

This is unseen-deposited-library validation within GSE216146. It is not an
unseen-study test and, because some deposited libraries pool material from
multiple animals, it is not equivalent to leave-one-animal-out validation.
GSE271055 cannot support an analogous treatment-balanced library holdout because
there is only one pooled control library and one pooled doxorubicin library.
"""
from __future__ import annotations

import copy
import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn.decomposition import PCA

R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R))
import src.models.mc_contrastive_vi as mm

OUT = R / "analysis" / "corrected_20260920"
RUN = R / "runs" / "corrected_20260920"
OUT.mkdir(parents=True, exist_ok=True)
RUN.mkdir(parents=True, exist_ok=True)

SPLIT_SEED = 1729
TRAIN_PER_STRATUM = 1500
VAL_PER_STRATUM = 100
N_HVG = 1500
SEEDS = (0, 1, 2)
MODES = ("full", "gaussian_vae")
EPOCHS = 100

FOLDS = {
    "replicate_1": {
        "control": "D20-6407",
        "cisplatin": "D20-6409",
    },
    "replicate_2": {
        "control": "D21-2746",
        "cisplatin": "D21-2750",
    },
    "replicate_3": {
        "control": "D21-2747",
        "cisplatin": "D21-2751",
    },
}

torch.set_num_threads(2)


def fast_hsic(x, y):
    if len(x) < 4:
        return x.new_tensor(0.0)
    K = mm._rbf_kernel(x)
    L = mm._rbf_kernel(y)
    K = K - K.mean(0, keepdim=True) - K.mean(1, keepdim=True) + K.mean()
    L = L - L.mean(0, keepdim=True) - L.mean(1, keepdim=True) + L.mean()
    return (K * L).sum() / (len(x) - 1) ** 2


mm.hsic = fast_hsic


def deterministic_means(model, x, d, b, mode):
    e = model.encode(x, d)
    bg = e["mu_bg"]
    if mode == "gaussian_vae":
        sh = e["mu_shared"]
        dr = torch.cat(e["mu_drug_list"], dim=1)
    else:
        sh = e["mu_shared"] * (d > 0)[:, None]
        dr = torch.cat(
            [
                mu * (d == k + 1)[:, None]
                for k, mu in enumerate(e["mu_drug_list"])
            ],
            dim=1,
        )
    pred = model.decode(bg, sh, dr, b)
    return pred


def make_fold_design(a, fold_name):
    obs = a.obs
    held = FOLDS[fold_name]
    held_libraries = set(held.values())

    pure = ~obs.rescue.astype(bool).to_numpy()
    sample = obs.sample_id.astype(str).to_numpy()
    study = obs.study.astype(str).to_numpy()
    drug = obs.drug.astype(str).to_numpy()

    pool = np.flatnonzero(pure & ~np.isin(sample, list(held_libraries)))
    test = np.flatnonzero(pure & np.isin(sample, list(held_libraries)))

    strata = np.char.add(np.char.add(study.astype(str), "_"), drug.astype(str))
    expected = (
        "GSE216146_control",
        "GSE216146_cisplatin",
        "GSE271055_control",
        "GSE271055_doxorubicin",
    )

    rng = np.random.default_rng(SPLIT_SEED + list(FOLDS).index(fold_name))
    train_parts = []
    val_parts = []
    design_rows = []

    for stratum in expected:
        idx = pool[strata[pool] == stratum]
        need = TRAIN_PER_STRATUM + VAL_PER_STRATUM
        if len(idx) < need:
            raise RuntimeError(
                f"{fold_name} {stratum}: only {len(idx)} available; {need} required"
            )
        chosen = rng.permutation(idx)[:need]
        tr = chosen[:TRAIN_PER_STRATUM]
        va = chosen[TRAIN_PER_STRATUM:]
        train_parts.append(tr)
        val_parts.append(va)
        design_rows.append(
            {
                "fold": fold_name,
                "stratum": stratum,
                "n_available_training_libraries_cells": int(len(idx)),
                "n_train": int(len(tr)),
                "n_validation": int(len(va)),
            }
        )

    train = np.concatenate(train_parts)
    val = np.concatenate(val_parts)

    assert len(train) == 6000
    assert len(val) == 400
    assert not np.isin(sample[train], list(held_libraries)).any()
    assert not np.isin(sample[val], list(held_libraries)).any()
    assert set(sample[test]) == held_libraries
    assert len(set(train) & set(val)) == 0
    assert len(set(train) & set(test)) == 0
    assert len(set(val) & set(test)) == 0

    for arm, sid in held.items():
        n = int(np.sum(sample[test] == sid))
        design_rows.append(
            {
                "fold": fold_name,
                "stratum": f"heldout_{arm}",
                "n_available_training_libraries_cells": n,
                "n_train": 0,
                "n_validation": 0,
            }
        )

    return train, val, test, pd.DataFrame(design_rows)


def train_model(X, d, b, train, val, mode, seed):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)

    cfg = mm.MCContrastiveVIConfig(
        n_genes=X.shape[1],
        n_drugs=3,
        n_batches=2,
        hidden_dim=128,
        dropout=0.1,
    )
    model = mm.MCContrastiveVI(cfg)
    opt = torch.optim.AdamW(model.parameters(), lr=0.001, weight_decay=1e-5)

    best = float("inf")
    best_epoch = -1
    state = None
    bad = 0
    history = []
    started = time.time()

    for epoch in range(EPOCHS):
        model.train()
        vals = []
        for ix in np.array_split(
            rng.permutation(train), int(np.ceil(len(train) / 256))
        ):
            x = X[ix]
            dd = d[ix]
            bb = b[ix]
            e = model(x, dd, bb)

            if mode == "gaussian_vae":
                e = model.encode(x, dd)
                e["z_shared_gated"] = e["z_shared"]
                e["z_drug_cat"] = torch.cat(
                    [
                        mm.reparameterize(mu, lv)
                        for mu, lv in zip(
                            e["mu_drug_list"], e["lv_drug_list"]
                        )
                    ],
                    dim=1,
                )
                e["recon"] = model.decode(
                    e["z_bg"], e["z_shared_gated"], e["z_drug_cat"], bb
                )
                loss = mm.gaussian_nll(
                    x, e["recon"], model.log_sigma
                ).mean()
                loss = loss + sum(
                    mm.kl_standard_normal(mu, lv).mean()
                    for mu, lv in [
                        (e["mu_bg"], e["lv_bg"]),
                        (e["mu_shared"], e["lv_shared"]),
                        *list(zip(e["mu_drug_list"], e["lv_drug_list"])),
                    ]
                )
            else:
                loss = mm.mc_contrastive_loss(
                    model, e, x, dd, 1.0, 10.0, 10.0, 5.0
                )["loss"]

            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
            opt.step()
            vals.append(float(loss.detach()))

        model.eval()
        sq = 0.0
        n = 0
        with torch.no_grad():
            for ix in np.array_split(val, int(np.ceil(len(val) / 512))):
                pred = deterministic_means(
                    model, X[ix], d[ix], b[ix], mode
                )
                sq += float(((pred - X[ix]) ** 2).sum())
                n += X[ix].numel()
        val_mse = sq / n
        history.append(
            {
                "epoch": epoch,
                "train_objective": float(np.mean(vals)),
                "validation_mse": val_mse,
            }
        )

        if val_mse < best - 1e-5:
            best = val_mse
            best_epoch = epoch
            state = copy.deepcopy(model.state_dict())
            bad = 0
        else:
            bad += 1

        if epoch % 20 == 0:
            print(
                mode,
                seed,
                "epoch",
                epoch,
                "val_mse",
                round(val_mse, 6),
                "seconds",
                round(time.time() - started),
                flush=True,
            )
        if bad >= 15:
            break

    assert state is not None
    model.load_state_dict(state)
    model.eval()
    return model, best, best_epoch, pd.DataFrame(history)


def per_cell_mse(model, X, d, b, indices, mode):
    out = []
    model.eval()
    with torch.no_grad():
        for ix in np.array_split(indices, int(np.ceil(len(indices) / 512))):
            pred = deterministic_means(model, X[ix], d[ix], b[ix], mode)
            out.append(((pred - X[ix]) ** 2).mean(1).cpu().numpy())
    return np.concatenate(out)


def summarize_mse(
    per_cell,
    indices,
    obs,
    fold_name,
    mode,
    seed,
    split,
):
    rows = []
    frame = pd.DataFrame(
        {
            "mse": per_cell,
            "sample_id": obs.sample_id.astype(str).to_numpy()[indices],
            "arm": obs.arm.astype(str).to_numpy()[indices],
            "drug": obs.drug.astype(str).to_numpy()[indices],
            "cell_type": obs.source_cell_type.astype(str).to_numpy()[indices],
        }
    )

    for (sid, arm), g in frame.groupby(["sample_id", "arm"], sort=True):
        rows.append(
            {
                "fold": fold_name,
                "mode": mode,
                "seed": seed,
                "split": split,
                "sample_id": sid,
                "arm": arm,
                "cell_type": "__all__",
                "n_cells": len(g),
                "mse": float(g.mse.mean()),
            }
        )
    for (sid, arm, ct), g in frame.groupby(
        ["sample_id", "arm", "cell_type"], sort=True
    ):
        rows.append(
            {
                "fold": fold_name,
                "mode": mode,
                "seed": seed,
                "split": split,
                "sample_id": sid,
                "arm": arm,
                "cell_type": ct,
                "n_cells": len(g),
                "mse": float(g.mse.mean()),
            }
        )
    return frame, rows


def arm_balanced_mse(frame):
    by = frame.groupby("drug").mse.mean()
    required = {"control", "cisplatin"}
    if not required.issubset(by.index):
        raise RuntimeError(f"Missing required drugs: {set(by.index)}")
    return float((by["control"] + by["cisplatin"]) / 2)


def pca_baseline(X_np, train, val, test):
    pca = PCA(n_components=32, random_state=SPLIT_SEED)
    pca.fit(X_np[train])

    def mse(indices):
        z = pca.transform(X_np[indices])
        pred = pca.inverse_transform(z)
        return ((pred - X_np[indices]) ** 2).mean(axis=1)

    return mse(val), mse(test)


def mean_baseline(X_np, train, val, test):
    mu = X_np[train].mean(axis=0, keepdims=True)
    return (
        ((X_np[val] - mu) ** 2).mean(axis=1),
        ((X_np[test] - mu) ** 2).mean(axis=1),
    )


def main():
    raw = ad.read_h5ad(R / "data" / "processed" / "recovered_counts.h5ad")
    obs = raw.obs.copy()

    fold_specs = {}
    design_frames = []
    gene_rows = []

    for fold_name in FOLDS:
        train, val, test, design = make_fold_design(raw, fold_name)
        v = raw[train].copy()
        sc.pp.highly_variable_genes(
            v,
            flavor="seurat_v3",
            n_top_genes=N_HVG,
            batch_key="study",
        )
        genes = v.var_names[v.var.highly_variable].astype(str).tolist()
        if len(genes) != N_HVG:
            raise RuntimeError(
                f"{fold_name}: expected {N_HVG} HVGs, got {len(genes)}"
            )
        fold_specs[fold_name] = {
            "train": train,
            "val": val,
            "test": test,
            "genes": genes,
        }
        design_frames.append(design)
        gene_rows.extend(
            {
                "fold": fold_name,
                "rank_in_var_order": rank,
                "gene": gene,
            }
            for rank, gene in enumerate(genes)
        )

    pd.concat(design_frames, ignore_index=True).to_csv(
        OUT / "task5_fold_design.csv", index=False
    )
    pd.DataFrame(gene_rows).to_csv(
        OUT / "task5_gene_universes.csv", index=False
    )

    # Cell-wise normalization is independent of the fold. Gene subsetting occurs
    # after each fold-specific training-only HVG set has been fixed.
    norm = raw.copy()
    sc.pp.normalize_total(norm, target_sum=1e4)
    sc.pp.log1p(norm)

    d_all = obs.drug.map(
        {"control": 0, "doxorubicin": 1, "cisplatin": 2}
    ).to_numpy(dtype=np.int64)
    b_all = (obs.study == "GSE271055").to_numpy(dtype=np.int64)

    metric_rows = []
    celltype_rows = []

    for fold_name, spec in fold_specs.items():
        print("\n===== FOLD", fold_name, "=====", flush=True)
        genes = spec["genes"]
        train = spec["train"]
        val = spec["val"]
        test = spec["test"]

        X_np = norm[:, genes].X.toarray().astype(np.float32)
        X = torch.from_numpy(X_np)
        d = torch.from_numpy(d_all)
        b = torch.from_numpy(b_all)

        # Seen-library validation restricted to GSE216146 control/cisplatin,
        # matching the treatment classes present in the unseen test pair.
        val_cis_study = val[b_all[val] == 0]
        assert len(val_cis_study) == 200

        for mode in MODES:
            for seed in SEEDS:
                print(
                    f"Training {fold_name} {mode} seed={seed}",
                    flush=True,
                )
                model, best_val, best_epoch, history = train_model(
                    X, d, b, train, val, mode, seed
                )

                val_pc = per_cell_mse(
                    model, X, d, b, val_cis_study, mode
                )
                test_pc = per_cell_mse(model, X, d, b, test, mode)

                val_frame, val_rows = summarize_mse(
                    val_pc,
                    val_cis_study,
                    obs,
                    fold_name,
                    mode,
                    seed,
                    "seen_library_validation",
                )
                test_frame, test_rows = summarize_mse(
                    test_pc,
                    test,
                    obs,
                    fold_name,
                    mode,
                    seed,
                    "unseen_library_test",
                )
                celltype_rows.extend(val_rows)
                celltype_rows.extend(test_rows)

                val_bal = arm_balanced_mse(val_frame)
                test_bal = arm_balanced_mse(test_frame)
                test_by = test_frame.groupby("drug").mse.mean()

                metric_rows.append(
                    {
                        "fold": fold_name,
                        "heldout_control_library": FOLDS[fold_name][
                            "control"
                        ],
                        "heldout_cisplatin_library": FOLDS[fold_name][
                            "cisplatin"
                        ],
                        "mode": mode,
                        "seed": seed,
                        "best_epoch": best_epoch,
                        "training_selection_val_mse_all_strata": best_val,
                        "seen_library_validation_balanced_mse": val_bal,
                        "unseen_library_balanced_mse": test_bal,
                        "unseen_minus_seen_mse": test_bal - val_bal,
                        "unseen_over_seen_mse_ratio": test_bal / val_bal,
                        "unseen_control_mse": float(test_by["control"]),
                        "unseen_cisplatin_mse": float(
                            test_by["cisplatin"]
                        ),
                        "n_unseen_cells": int(len(test)),
                    }
                )

                history.to_csv(
                    RUN
                    / f"task5_{fold_name}_{mode}_{seed}_history.csv",
                    index=False,
                )

        for mode, fn in (
            ("pca32", pca_baseline),
            ("training_mean", mean_baseline),
        ):
            val_pc, test_pc = fn(X_np, train, val_cis_study, test)
            val_frame, val_rows = summarize_mse(
                val_pc,
                val_cis_study,
                obs,
                fold_name,
                mode,
                -1,
                "seen_library_validation",
            )
            test_frame, test_rows = summarize_mse(
                test_pc,
                test,
                obs,
                fold_name,
                mode,
                -1,
                "unseen_library_test",
            )
            celltype_rows.extend(val_rows)
            celltype_rows.extend(test_rows)
            val_bal = arm_balanced_mse(val_frame)
            test_bal = arm_balanced_mse(test_frame)
            test_by = test_frame.groupby("drug").mse.mean()
            metric_rows.append(
                {
                    "fold": fold_name,
                    "heldout_control_library": FOLDS[fold_name]["control"],
                    "heldout_cisplatin_library": FOLDS[fold_name][
                        "cisplatin"
                    ],
                    "mode": mode,
                    "seed": -1,
                    "best_epoch": np.nan,
                    "training_selection_val_mse_all_strata": np.nan,
                    "seen_library_validation_balanced_mse": val_bal,
                    "unseen_library_balanced_mse": test_bal,
                    "unseen_minus_seen_mse": test_bal - val_bal,
                    "unseen_over_seen_mse_ratio": test_bal / val_bal,
                    "unseen_control_mse": float(test_by["control"]),
                    "unseen_cisplatin_mse": float(
                        test_by["cisplatin"]
                    ),
                    "n_unseen_cells": int(len(test)),
                }
            )

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "task5_unseen_library_metrics.csv", index=False)
    pd.DataFrame(celltype_rows).to_csv(
        OUT / "task5_unseen_library_mse_by_celltype.csv", index=False
    )

    # Equal-weight folds. Neural modes are first averaged across seeds within a
    # fold, then across folds, so the large replicate-2 library does not
    # dominate the primary summary.
    fold_summary = (
        metrics.groupby(["fold", "mode"], as_index=False)
        .agg(
            n_seed_runs=("seed", "size"),
            mean_seen_library_validation_balanced_mse=(
                "seen_library_validation_balanced_mse",
                "mean",
            ),
            sd_seen_library_validation_balanced_mse=(
                "seen_library_validation_balanced_mse",
                "std",
            ),
            mean_unseen_library_balanced_mse=(
                "unseen_library_balanced_mse",
                "mean",
            ),
            sd_unseen_library_balanced_mse=(
                "unseen_library_balanced_mse",
                "std",
            ),
            mean_unseen_minus_seen_mse=(
                "unseen_minus_seen_mse",
                "mean",
            ),
            mean_unseen_over_seen_mse_ratio=(
                "unseen_over_seen_mse_ratio",
                "mean",
            ),
            mean_unseen_control_mse=("unseen_control_mse", "mean"),
            mean_unseen_cisplatin_mse=(
                "unseen_cisplatin_mse",
                "mean",
            ),
            n_unseen_cells=("n_unseen_cells", "first"),
        )
    )
    fold_summary.to_csv(
        OUT / "task5_unseen_library_fold_summary.csv", index=False
    )

    model_summary = (
        fold_summary.groupby("mode", as_index=False)
        .agg(
            n_folds=("fold", "size"),
            mean_seen_library_validation_balanced_mse=(
                "mean_seen_library_validation_balanced_mse",
                "mean",
            ),
            mean_unseen_library_balanced_mse=(
                "mean_unseen_library_balanced_mse",
                "mean",
            ),
            sd_unseen_library_balanced_mse_across_folds=(
                "mean_unseen_library_balanced_mse",
                "std",
            ),
            mean_generalization_gap=(
                "mean_unseen_minus_seen_mse",
                "mean",
            ),
            mean_unseen_over_seen_mse_ratio=(
                "mean_unseen_over_seen_mse_ratio",
                "mean",
            ),
        )
    )
    model_summary.to_csv(
        OUT / "task5_unseen_library_model_summary.csv", index=False
    )

    design = {
        "validation_type": (
            "leave-one-GEO-replicate-pair-out deposited-library validation"
        ),
        "heldout_pairs": FOLDS,
        "rescue_arms_used": False,
        "train_per_study_treatment_stratum": TRAIN_PER_STRATUM,
        "validation_per_study_treatment_stratum": VAL_PER_STRATUM,
        "training_cells_per_fold": 6000,
        "validation_cells_per_fold": 400,
        "hvg_count": N_HVG,
        "hvg_selection": "training cells only, seurat_v3, batch_key=study",
        "neural_modes": list(MODES),
        "neural_seeds": list(SEEDS),
        "baselines": ["PCA(32)", "training gene-wise mean"],
        "primary_metric": (
            "within-fold mean of heldout control-library MSE and "
            "heldout cisplatin-library MSE; folds weighted equally"
        ),
        "scope": (
            "unseen deposited libraries within GSE216146; not unseen study "
            "and not necessarily leave-one-animal-out"
        ),
        "GSE271055_limitation": (
            "no analogous treatment-balanced library holdout because only "
            "one pooled control and one pooled doxorubicin library exist"
        ),
    }
    (OUT / "task5_unseen_library_design.json").write_text(
        json.dumps(design, indent=2) + "\n"
    )

    print("\n===== TASK 5 FOLD SUMMARY =====")
    print(fold_summary.to_string(index=False))
    print("\n===== TASK 5 MODEL SUMMARY =====")
    print(model_summary.to_string(index=False))


if __name__ == "__main__":
    main()
