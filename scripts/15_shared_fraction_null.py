"""Exact library-label reference distribution for the shared latent variance fraction.

This analysis is intentionally separate from the manuscript's canonical cell-level
benchmark. It constructs a label-independent, balanced library benchmark so that
all treatment-label permutations use exactly the same cells, feature set, split,
and number of cells per library.

Eight non-rescue deposited libraries are used:
  * GSE216146: 6 libraries (cisplatin study), paired by the three GEO replicate blocks.
  * GSE271055: 2 pooled libraries (doxorubicin study).

Cell sampling is fixed before any label permutation and is balanced both within
study and across study-by-treatment strata for every assignment. Each
GSE216146 library contributes 500 train, 50 validation and 50 test cells; each
GSE271055 pooled library contributes 1,500 train, 150 validation and 150 test
cells. Consequently, every assignment contains 3,000 training cells per study
and 1,500 training cells in each study-by-treatment stratum, matching the
canonical benchmark's 6,000-cell training size without using treatment labels
to select cells. Highly variable genes are selected from these fixed training
cells using study as the batch key and are then held fixed across assignments.

The exact label-exchangeability reference set has:
  2^3 * C(2,1) = 16 assignments,
where each of the three GSE216146 GEO replicate blocks contributes one
control/cisplatin swap and the last factor exchanges the two pooled
GSE271055 control/doxorubicin libraries.

The primary statistic is the mean pooled shared posterior-mean variance fraction
across fixed model seeds 0, 1, and 2. Exact upper-tail P values include the
observed assignment in the reference distribution. They are conditional on the
stated library-label exchangeability assumption and should not be interpreted as
a randomized-trial biological significance test.
"""
from __future__ import annotations

import argparse
import copy
import itertools
import json
import math
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import torch

R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R))

import src.models.mc_contrastive_vi as mm

OUT = R / "analysis" / "corrected_20260920"
RUN = R / "runs" / "corrected_20260920"
OUT.mkdir(parents=True, exist_ok=True)
RUN.mkdir(parents=True, exist_ok=True)

SPLIT_SEED = 1729
SEEDS = (0, 1, 2)
SPLIT_COUNTS = {
    "GSE216146": {"train": 500, "validation": 50, "test": 50},
    "GSE271055": {"train": 1500, "validation": 150, "test": 150},
}
N_HVG = 1500
AU_THRESHOLD = 0.01
KL_THRESHOLD = 0.01

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


def _exchangeability_group(study: str, sample_id: str) -> str:
    """Return the experimental replicate block documented in the GEO series record."""
    if study == "GSE216146":
        replicate_map = {
            "D20-6407": "replicate_1",
            "D20-6409": "replicate_1",
            "D21-2746": "replicate_2",
            "D21-2750": "replicate_2",
            "D21-2747": "replicate_3",
            "D21-2751": "replicate_3",
        }
        if sample_id not in replicate_map:
            raise ValueError(f"Unexpected non-rescue GSE216146 sample: {sample_id}")
        return replicate_map[sample_id]
    if study == "GSE271055":
        return "pooled_pair"
    raise ValueError(f"Unexpected study: {study}")


def _enumerate_assignments(lib_meta: pd.DataFrame) -> pd.DataFrame:
    """Enumerate all assignments preserving observed treated counts per stratum."""
    cis = lib_meta[lib_meta.study == "GSE216146"].copy()
    dox = lib_meta[lib_meta.study == "GSE271055"].copy()

    cis_groups = []
    for replicate, g in cis.groupby("exchangeability_group", sort=True):
        libs = tuple(sorted(g.sample_id.astype(str)))
        n_treated = int((g.original_drug == "cisplatin").sum())
        assert len(libs) == 2 and n_treated == 1
        combos = list(itertools.combinations(libs, n_treated))
        cis_groups.append((replicate, combos))

    dox_libs = tuple(sorted(dox.sample_id.astype(str)))
    n_dox_treated = int((dox.original_drug == "doxorubicin").sum())
    dox_combos = list(itertools.combinations(dox_libs, n_dox_treated))

    observed_cis = set(
        cis.loc[cis.original_drug == "cisplatin", "sample_id"].astype(str)
    )
    observed_dox = set(
        dox.loc[dox.original_drug == "doxorubicin", "sample_id"].astype(str)
    )

    rows = []
    assignment_index = 0
    cis_products = itertools.product(*(combos for _, combos in cis_groups))
    for cis_choice_tuple in cis_products:
        cis_treated = set().union(*(set(x) for x in cis_choice_tuple))
        replicate_choice = {
            replicate: "|".join(choice)
            for (replicate, _), choice in zip(cis_groups, cis_choice_tuple)
        }
        for dox_choice in dox_combos:
            dox_treated = set(dox_choice)
            row = {
                "assignment_index": assignment_index,
                "cis_treated_libraries": "|".join(sorted(cis_treated)),
                "dox_treated_libraries": "|".join(sorted(dox_treated)),
                "is_observed": (
                    cis_treated == observed_cis and dox_treated == observed_dox
                ),
                "cis_matches_observed": cis_treated == observed_cis,
                "dox_matches_observed": dox_treated == observed_dox,
            }
            for replicate, choice in replicate_choice.items():
                row[f"{replicate}_cis_treated"] = choice
            rows.append(row)
            assignment_index += 1

    out = pd.DataFrame(rows)
    assert len(out) == 16, f"Expected 16 assignments, got {len(out)}"
    assert out.is_observed.sum() == 1
    assert out.cis_matches_observed.sum() == 2
    assert out.dox_matches_observed.sum() == 8
    return out


def prepare() -> None:
    a = ad.read_h5ad(R / "data" / "processed" / "recovered_counts.h5ad")
    obs = a.obs.copy()
    pure = obs.loc[~obs.rescue.astype(bool)].copy()

    lib_meta = (
        pure[["study", "sample_id", "drug"]]
        .drop_duplicates()
        .rename(columns={"drug": "original_drug"})
        .sort_values(["study", "sample_id"])
        .reset_index(drop=True)
    )
    lib_meta["exchangeability_group"] = [
        _exchangeability_group(study, sid)
        for study, sid in zip(lib_meta.study.astype(str), lib_meta.sample_id.astype(str))
    ]

    assert len(lib_meta) == 8, f"Expected 8 non-rescue libraries, got {len(lib_meta)}"
    assert set(lib_meta.study) == {"GSE216146", "GSE271055"}

    split_path = OUT / "task2_balanced_library_split.csv"
    if split_path.exists():
        split_df = pd.read_csv(split_path)
        required = {"global_index", "sample_id", "split"}
        assert required.issubset(split_df.columns)
        assert len(split_df) == 7200
        gi = split_df.global_index.to_numpy(dtype=np.int64)
        assert np.all(gi >= 0) and np.all(gi < len(obs))
        assert not obs.rescue.astype(bool).to_numpy()[gi].any()
        assert np.array_equal(
            obs.sample_id.astype(str).to_numpy()[gi],
            split_df.sample_id.astype(str).to_numpy(),
        )
        train_global = split_df.loc[
            split_df.split == "train", "global_index"
        ].to_numpy(dtype=np.int64)
        val_global = split_df.loc[
            split_df.split == "validation", "global_index"
        ].to_numpy(dtype=np.int64)
        test_global = split_df.loc[
            split_df.split == "test", "global_index"
        ].to_numpy(dtype=np.int64)
        split_source = "existing_task2_balanced_library_split"
    else:
        rng = np.random.default_rng(SPLIT_SEED)
        split_rows = []
        train_global, val_global, test_global = [], [], []

        for sid in sorted(lib_meta.sample_id.astype(str)):
            study_for_library = str(
                lib_meta.loc[lib_meta.sample_id.astype(str) == sid, "study"].iloc[0]
            )
            counts = SPLIT_COUNTS[study_for_library]
            n_required = sum(counts.values())
            idx = np.flatnonzero(
                (~obs.rescue.astype(bool).to_numpy())
                & (obs.sample_id.astype(str).to_numpy() == sid)
            )
            assert len(idx) >= n_required, (
                f"{sid} has only {len(idx)} retained pure cells; "
                f"{n_required} required"
            )
            chosen = rng.permutation(idx)[:n_required]
            ntr = counts["train"]
            nva = counts["validation"]
            parts = {
                "train": chosen[:ntr],
                "validation": chosen[ntr : ntr + nva],
                "test": chosen[ntr + nva : n_required],
            }
            assert len(parts["test"]) == counts["test"]

            train_global.extend(parts["train"].tolist())
            val_global.extend(parts["validation"].tolist())
            test_global.extend(parts["test"].tolist())

            for split_name, vals in parts.items():
                for gi in vals:
                    split_rows.append(
                        {
                            "global_index": int(gi),
                            "cell_id": str(a.obs_names[gi]),
                            "study": str(obs.study.iloc[gi]),
                            "sample_id": str(obs.sample_id.iloc[gi]),
                            "exchangeability_group": _exchangeability_group(
                                str(obs.study.iloc[gi]),
                                str(obs.sample_id.iloc[gi]),
                            ),
                            "original_drug": str(obs.drug.iloc[gi]),
                            "source_cell_type": str(
                                obs.source_cell_type.iloc[gi]
                            ),
                            "split": split_name,
                        }
                    )

        train_global = np.asarray(train_global, dtype=np.int64)
        val_global = np.asarray(val_global, dtype=np.int64)
        test_global = np.asarray(test_global, dtype=np.int64)
        split_df = pd.DataFrame(split_rows)
        split_df.to_csv(split_path, index=False)
        split_source = "computed_label_independent_split"

    assert len(train_global) == 6000
    assert len(val_global) == 600
    assert len(test_global) == 600
    assert len(set(train_global) & set(val_global)) == 0
    assert len(set(train_global) & set(test_global)) == 0
    assert len(set(val_global) & set(test_global)) == 0

    gene_path = OUT / "task2_gene_universe.csv"
    if gene_path.exists():
        genes = pd.read_csv(gene_path)["gene"].astype(str).tolist()
        missing = [g for g in genes if g not in a.var_names]
        assert not missing, f"Committed Task 2 genes missing from cohort: {missing}"
        assert len(genes) == N_HVG
        gene_source = "existing_task2_gene_universe"
    else:
        v = a[train_global].copy()
        sc.pp.highly_variable_genes(
            v,
            flavor="seurat_v3",
            n_top_genes=N_HVG,
            batch_key="study",
        )
        genes = v.var_names[v.var.highly_variable].astype(str).tolist()
        assert len(genes) == N_HVG
        pd.DataFrame({"gene": genes}).to_csv(gene_path, index=False)
        gene_source = "computed_from_label_independent_training_cells"

    selected_global = np.concatenate([train_global, val_global, test_global])
    sub = a[selected_global].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    X = sub[:, genes].X.toarray().astype("float32")

    selected_obs = obs.iloc[selected_global].copy()
    sample_id = selected_obs.sample_id.astype(str).to_numpy()
    study = selected_obs.study.astype(str).to_numpy()
    b = (study == "GSE271055").astype("int64")

    n_train = len(train_global)
    n_val = len(val_global)
    train = np.arange(0, n_train, dtype=np.int64)
    val = np.arange(n_train, n_train + n_val, dtype=np.int64)
    test = np.arange(n_train + n_val, len(selected_global), dtype=np.int64)

    np.savez_compressed(
        RUN / "task2_balanced_inputs.npz",
        X=X,
        b=b,
        sample_id=np.asarray(sample_id, dtype="U"),
        study=np.asarray(study, dtype="U"),
        train=train,
        val=val,
        test=test,
        genes=np.asarray(genes, dtype="U"),
        global_index=selected_global,
    )

    lib_meta.to_csv(OUT / "task2_library_design.csv", index=False)

    assignments = _enumerate_assignments(lib_meta)
    assignments.to_csv(OUT / "task2_label_assignments.csv", index=False)

    metadata = {
        "split_seed": SPLIT_SEED,
        "split_source": split_source,
        "per_library_split_counts": SPLIT_COUNTS,
        "n_libraries": int(len(lib_meta)),
        "n_train": int(len(train)),
        "n_validation": int(len(val)),
        "n_test": int(len(test)),
        "n_hvg": int(len(genes)),
        "gene_source": gene_source,
        "n_exact_label_assignments": int(len(assignments)),
        "primary_statistic": "mean pooled shared_fraction across seeds 0,1,2",
        "exchangeability_constraints": {
            "GSE216146_replicate_1": "exchange control/cisplatin labels only within replicate 1",
            "GSE216146_replicate_2": "exchange control/cisplatin labels only within replicate 2",
            "GSE216146_replicate_3": "exchange control/cisplatin labels only within replicate 3",
            "GSE271055": "exchange the control and doxorubicin pooled-library labels",
        },
        "inference_scope": (
            "Exact upper-tail label-exchangeability reference distribution; "
            "conditional on exchangeability assumptions and not a randomized-trial "
            "biological significance test."
        ),
    }
    (OUT / "task2_design.json").write_text(json.dumps(metadata, indent=2) + "\n")

    print("===== TASK 2 PREPARATION =====")
    print(json.dumps(metadata, indent=2))
    print("\nLibraries:")
    print(lib_meta.to_string(index=False))
    print("\nAssignments:", len(assignments))
    print("Observed assignment index:", int(assignments.loc[assignments.is_observed, "assignment_index"].iloc[0]))


def _labels_for_assignment(sample_ids: np.ndarray, row: pd.Series) -> np.ndarray:
    cis = set(str(row.cis_treated_libraries).split("|"))
    dox = set(str(row.dox_treated_libraries).split("|"))
    d = np.zeros(len(sample_ids), dtype=np.int64)
    d[np.isin(sample_ids, list(dox))] = 1
    d[np.isin(sample_ids, list(cis))] = 2
    return d


def _means(model, X, d, b):
    e = model.encode(X, d)
    sh = e["mu_shared"] * (d > 0)[:, None]
    dr = torch.cat(
        [
            mu * (d == k + 1)[:, None]
            for k, mu in enumerate(e["mu_drug_list"])
        ],
        dim=1,
    )
    pred = model.decode(e["mu_bg"], sh, dr, b)
    return e, sh, dr, pred


def _per_dim_kl(mu: np.ndarray, logvar: np.ndarray) -> np.ndarray:
    return 0.5 * (mu**2 + np.exp(logvar) - 1.0 - logvar)


def fit_assignment(
    X: torch.Tensor,
    b: torch.Tensor,
    d: torch.Tensor,
    train: np.ndarray,
    val: np.ndarray,
    test: np.ndarray,
    seed: int,
    epochs: int,
) -> dict:
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
    bad = 0
    best_epoch = -1
    state = None

    for epoch in range(epochs):
        model.train()
        for ix in np.array_split(
            rng.permutation(train), int(np.ceil(len(train) / 256))
        ):
            e = model(X[ix], d[ix], b[ix])
            loss = mm.mc_contrastive_loss(
                model,
                e,
                X[ix],
                d[ix],
                1.0,
                10.0,
                10.0,
                5.0,
            )["loss"]
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
            opt.step()

        model.eval()
        sq = 0.0
        n = 0
        with torch.no_grad():
            for ix in np.array_split(val, int(np.ceil(len(val) / 512))):
                _, _, _, pred = _means(model, X[ix], d[ix], b[ix])
                sq += float(((pred - X[ix]) ** 2).sum())
                n += X[ix].numel()
        val_mse = sq / n

        if val_mse < best - 1e-5:
            best = val_mse
            state = copy.deepcopy(model.state_dict())
            best_epoch = epoch
            bad = 0
        else:
            bad += 1
        if bad >= 15:
            break

    assert state is not None
    model.load_state_dict(state)
    model.eval()

    with torch.no_grad():
        e, sh, dr, pred = _means(model, X[test], d[test], b[test])

    test_np = test
    d_test = d[test].cpu().numpy()
    treated_pos = np.flatnonzero(d_test > 0)
    dox_pos = np.flatnonzero(d_test == 1)
    cis_pos = np.flatnonzero(d_test == 2)

    sh_np = sh.cpu().numpy()
    dr_np = dr.cpu().numpy()
    pred_np = pred.cpu().numpy()
    X_test_np = X[test].cpu().numpy()

    def frac(a, c):
        vs = np.var(a, axis=0).sum()
        vd = np.var(c, axis=0).sum()
        return float(vs / (vs + vd)), float(vs), float(vd)

    pooled_fsh, shared_var_sum, drug_var_sum = frac(
        sh_np[treated_pos], dr_np[treated_pos]
    )

    # Within-condition fractions use the relevant 4-dimensional drug block.
    shared_mu_all = e["mu_shared"].cpu().numpy()
    dox_mu_all = e["mu_drug_list"][0].cpu().numpy()
    cis_mu_all = e["mu_drug_list"][1].cpu().numpy()
    dox_fsh, _, _ = frac(shared_mu_all[dox_pos], dox_mu_all[dox_pos])
    cis_fsh, _, _ = frac(shared_mu_all[cis_pos], cis_mu_all[cis_pos])

    shared_mu = shared_mu_all[treated_pos]
    shared_lv = e["lv_shared"].cpu().numpy()[treated_pos]
    dox_mu = dox_mu_all[dox_pos]
    dox_lv = e["lv_drug_list"][0].cpu().numpy()[dox_pos]
    cis_mu = cis_mu_all[cis_pos]
    cis_lv = e["lv_drug_list"][1].cpu().numpy()[cis_pos]

    sh_var_dim = np.var(shared_mu, axis=0)
    dox_var_dim = np.var(dox_mu, axis=0)
    cis_var_dim = np.var(cis_mu, axis=0)

    sh_kl_dim = _per_dim_kl(shared_mu, shared_lv).mean(axis=0)
    dox_kl_dim = _per_dim_kl(dox_mu, dox_lv).mean(axis=0)
    cis_kl_dim = _per_dim_kl(cis_mu, cis_lv).mean(axis=0)

    return {
        "seed": seed,
        "best_epoch": best_epoch,
        "epochs_run": epoch + 1,
        "validation_mse": float(best),
        "test_mse": float(np.mean((pred_np - X_test_np) ** 2)),
        "shared_fraction": pooled_fsh,
        "doxorubicin_shared_fraction": dox_fsh,
        "cisplatin_shared_fraction": cis_fsh,
        "shared_var_sum": shared_var_sum,
        "drug_var_sum": drug_var_sum,
        "shared_active_units": int((sh_var_dim > AU_THRESHOLD).sum()),
        "dox_active_units": int((dox_var_dim > AU_THRESHOLD).sum()),
        "cis_active_units": int((cis_var_dim > AU_THRESHOLD).sum()),
        "shared_dims_kl_gt_0p01": int((sh_kl_dim > KL_THRESHOLD).sum()),
        "dox_dims_kl_gt_0p01": int((dox_kl_dim > KL_THRESHOLD).sum()),
        "cis_dims_kl_gt_0p01": int((cis_kl_dim > KL_THRESHOLD).sum()),
        "shared_mean_total_kl_nats": float(sh_kl_dim.sum()),
        "dox_mean_total_kl_nats": float(dox_kl_dim.sum()),
        "cis_mean_total_kl_nats": float(cis_kl_dim.sum()),
        "n_treated_test": int(len(treated_pos)),
        "n_doxorubicin_test": int(len(dox_pos)),
        "n_cisplatin_test": int(len(cis_pos)),
    }


def _exact_upper_tail(values: np.ndarray, observed: float) -> tuple[int, float]:
    count = int(np.sum(values >= observed - 1e-12))
    return count, count / len(values)


def summarize(fits: pd.DataFrame, assignments: pd.DataFrame) -> None:
    grouped = (
        fits.groupby(
            [
                "assignment_index",
                "is_observed",
                "cis_matches_observed",
                "dox_matches_observed",
            ],
            as_index=False,
        )
        .agg(
            mean_shared_fraction=("shared_fraction", "mean"),
            sd_shared_fraction=("shared_fraction", "std"),
            min_shared_fraction=("shared_fraction", "min"),
            max_shared_fraction=("shared_fraction", "max"),
            mean_test_mse=("test_mse", "mean"),
            mean_shared_var_sum=("shared_var_sum", "mean"),
            mean_drug_var_sum=("drug_var_sum", "mean"),
            mean_shared_active_units=("shared_active_units", "mean"),
            mean_drug_active_units=(
                "dox_active_units",
                lambda x: float(np.mean(x)),
            ),
        )
    )

    # Replace the temporary drug-active summary with the actual combined block count.
    drug_active = (
        fits.assign(
            drug_active_units_total=fits.dox_active_units + fits.cis_active_units
        )
        .groupby("assignment_index")
        .drug_active_units_total.mean()
    )
    grouped["mean_drug_active_units"] = grouped.assignment_index.map(drug_active)

    grouped.to_csv(OUT / "task2_assignment_summary.csv", index=False)

    observed = grouped[grouped.is_observed]
    assert len(observed) == 1
    obs_value = float(observed.mean_shared_fraction.iloc[0])

    scopes = {
        "combined_16": np.ones(len(grouped), dtype=bool),
        "cis_only_dox_fixed_8": grouped.dox_matches_observed.to_numpy(dtype=bool),
        "dox_only_cis_fixed_2": grouped.cis_matches_observed.to_numpy(dtype=bool),
    }

    null_rows = []
    for name, mask in scopes.items():
        sub = grouped.loc[mask].copy()
        values = sub.mean_shared_fraction.to_numpy()
        count, p = _exact_upper_tail(values, obs_value)
        null_rows.append(
            {
                "scope": name,
                "n_assignments": len(values),
                "statistic": "mean_shared_fraction_across_seeds_0_1_2",
                "observed_value": obs_value,
                "upper_tail_count": count,
                "exact_upper_tail_p": p,
                "minimum_attainable_p": 1 / len(values),
                "reference_min": float(np.min(values)),
                "reference_median": float(np.median(values)),
                "reference_max": float(np.max(values)),
                "observed_rank_descending": int(
                    1 + np.sum(values > obs_value + 1e-12)
                ),
            }
        )
    pd.DataFrame(null_rows).to_csv(OUT / "task2_null_test_summary.csv", index=False)

    seedwise_rows = []
    for seed in SEEDS:
        f = fits[fits.seed == seed].copy()
        obs = float(f.loc[f.is_observed, "shared_fraction"].iloc[0])
        seed_scopes = {
            "combined_16": np.ones(len(f), dtype=bool),
            "cis_only_dox_fixed_8": f.dox_matches_observed.to_numpy(dtype=bool),
            "dox_only_cis_fixed_2": f.cis_matches_observed.to_numpy(dtype=bool),
        }
        for name, mask in seed_scopes.items():
            vals = f.loc[mask, "shared_fraction"].to_numpy()
            count, p = _exact_upper_tail(vals, obs)
            seedwise_rows.append(
                {
                    "seed": seed,
                    "scope": name,
                    "n_assignments": len(vals),
                    "observed_shared_fraction": obs,
                    "upper_tail_count": count,
                    "exact_upper_tail_p": p,
                    "reference_min": float(vals.min()),
                    "reference_median": float(np.median(vals)),
                    "reference_max": float(vals.max()),
                }
            )
    pd.DataFrame(seedwise_rows).to_csv(
        OUT / "task2_seedwise_null_summary.csv", index=False
    )


def run(epochs: int) -> None:
    z = np.load(RUN / "task2_balanced_inputs.npz")
    X = torch.from_numpy(z["X"])
    b = torch.from_numpy(z["b"])
    train = z["train"]
    val = z["val"]
    test = z["test"]
    sample_id = z["sample_id"].astype(str)

    assignments = pd.read_csv(OUT / "task2_label_assignments.csv")
    rows = []
    total = len(assignments) * len(SEEDS)
    started = time.time()
    counter = 0

    for _, assignment in assignments.iterrows():
        d_np = _labels_for_assignment(sample_id, assignment)
        d = torch.from_numpy(d_np)

        # The label-independent split keeps every study-treatment stratum balanced.
        assert int((d[train] == 0).sum()) == 3000
        assert int((d[train] == 1).sum()) == 1500
        assert int((d[train] == 2).sum()) == 1500
        assert int((d[test] > 0).sum()) == 300
        assert int((d[test] == 1).sum()) == 150
        assert int((d[test] == 2).sum()) == 150

        for seed in SEEDS:
            counter += 1
            print(
                f"[{counter}/{total}] assignment={int(assignment.assignment_index)} "
                f"seed={seed} observed={bool(assignment.is_observed)}",
                flush=True,
            )
            result = fit_assignment(X, b, d, train, val, test, seed, epochs)
            result.update(
                {
                    "assignment_index": int(assignment.assignment_index),
                    "is_observed": bool(assignment.is_observed),
                    "cis_matches_observed": bool(assignment.cis_matches_observed),
                    "dox_matches_observed": bool(assignment.dox_matches_observed),
                    "cis_treated_libraries": str(
                        assignment.cis_treated_libraries
                    ),
                    "dox_treated_libraries": str(
                        assignment.dox_treated_libraries
                    ),
                }
            )
            rows.append(result)
            elapsed = time.time() - started
            print(
                f"  F_sh={result['shared_fraction']:.6f} "
                f"MSE={result['test_mse']:.6f} "
                f"elapsed={elapsed/60:.1f} min",
                flush=True,
            )

            # Checkpoint progress after every fit so partial results survive local runs.
            pd.DataFrame(rows).to_csv(
                OUT / "task2_permutation_fits_partial.csv", index=False
            )

    fits = pd.DataFrame(rows).sort_values(
        ["assignment_index", "seed"]
    ).reset_index(drop=True)
    fits.to_csv(OUT / "task2_permutation_fits.csv", index=False)
    summarize(fits, assignments)

    print("\n===== TASK 2 PRIMARY NULL TEST =====")
    print(pd.read_csv(OUT / "task2_null_test_summary.csv").to_string(index=False))
    print("\n===== TASK 2 OBSERVED ASSIGNMENT =====")
    print(
        fits[fits.is_observed][
            [
                "seed",
                "shared_fraction",
                "test_mse",
                "shared_var_sum",
                "drug_var_sum",
                "shared_active_units",
                "dox_active_units",
                "cis_active_units",
            ]
        ].to_string(index=False)
    )
    print("\nWrote Task 2 outputs to", OUT)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--epochs", type=int, default=100)
    args = ap.parse_args()

    if not args.prepare and not args.run:
        ap.error("Specify --prepare and/or --run")
    if args.prepare:
        prepare()
    if args.run:
        run(args.epochs)


if __name__ == "__main__":
    main()