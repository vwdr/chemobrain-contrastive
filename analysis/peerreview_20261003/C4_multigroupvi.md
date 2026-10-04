# C4 — multiGroupVI baseline

Script: `scripts/34_c4_multigroupvi.py` (run with `.venv-mgvi`). Fit outputs (saved models, latents, histories):
`runs/peerreview_20261003/multigroupvi/`. Logs: `runs/peerreview_20261003/logs/34_*.log` (failed attempts in
`logs/c4_attempt1/`, `logs/c4_attempt2/`), time-box log `logs/c4_timebox.txt`, pip logs `logs/c4_pip_attempt*.log`.

## Implementation and environment

- Official implementation: https://github.com/Genentech/multiGroupVI, commit `3d001dbe6905803195165f936f879a1a166a6fc8`
  (Weinberger, Lopez, Hütter & Regev, MLCB 2022 / PMLR 200). Its `setup.cfg` requires `scvi-tools==0.18.0`,
  `protobuf<=3.20.1`, `scanpy>=1.8.1`. Installed with `pip install --no-deps` into a separate venv `.venv-mgvi`
  (Python 3.12.10, the only interpreter available). Package source not modified.
- Environment steps (full freeze in `C4_mgvi_pip_freeze.txt`):
  1. `scvi-tools==0.18.0` with unpinned transitive dependencies resolved to torch 2.14.1, numpy 2.5.3, anndata 0.13.4,
     pytorch-lightning 1.6.5 → `import scvi` failed (`pkg_resources` missing; scvi 0.18 imports
     `pytorch_lightning.loggers.logger`, which does not exist in 1.6.5).
  2. Pinned: setuptools 69.5.1, torch 2.5.1+cpu, numpy 1.26.4, scipy 1.13.1, pandas 2.1.4, jax/jaxlib 0.4.30,
     flax 0.8.5, optax 0.2.3, chex 0.1.86, numpyro 0.15.3, pyro-ppl 1.9.1, torchmetrics 0.11.4, zarr 2.18.7,
     numcodecs 0.12.1, anndata 0.9.2, mudata 0.2.3, scanpy 1.9.8, scikit-learn 1.5.2, matplotlib 3.9.4, h5py 3.12.1.
  3. `pytorch-lightning==1.7.7` (within scvi 0.18's `>=1.6,<1.8`) could not be installed by pip 26 ("No matching
     distribution"); pip was downgraded to 24.0 inside `.venv-mgvi`, after which it installed. `import scvi` (0.18.0)
     and `multigroup_vi` then worked.
- Run-time compatibility fix (call argument only): the package computes its default `max_epochs` as
  `np.min([round(20000 / n_cells * 400), 400])`, a numpy integer, which pytorch-lightning 1.7.7 rejects
  (`assert isinstance(self.max_epochs, int)`). The same value (400 for 6,000 cells) was passed explicitly as a Python int.
- Failed attempts: (1) the `max_epochs` assertion above (all 3 seeds, before training); (2) training completed but
  this run's post-processing code called the data loader on unregistered test AnnData (fixed by
  `model._validate_anndata`; models are now saved right after training; a 2-epoch smoke test confirmed the full path).
  Attempt 3 succeeded for all seeds.
- Time box: first install attempt 22:30:53; import working 22:47:13; model initialisation 22:49:56; fits deferred while
  Phase 4 occupied the 3 fit slots; fit attempts 02:55:08–03:34:38. Effort ≈ 59 min (< 4 h).

## Settings

- Groups: control (both studies; 3,000 training cells), cisplatin (1,500), doxorubicin (1,500); `group_key` only (no
  batch key, as in the authors' notebook).
- Data: canonical 6,000 training cells, 1,500 committed genes, raw counts (`count_inputs.npz`); canonical test cells
  (2,881). Seeds 0, 1, 2 (`scvi.settings.seed`, `torch.manual_seed`, `np.random.seed`).
- Model: package defaults — n_hidden 128, n_layers 1, n_shared_latent 10, n_private_latent 10, wasserstein_penalty 1,
  dropout 0.1, ZINB likelihood; training defaults — batch size 128, train_size 0.9 (of the 6,000 training cells),
  400 epochs, no early stopping; CPU, 3 threads per fit, 3 fits concurrently (1,077–1,079 s per fit). No
  hyperparameter search.
- Outputs: `get_normalized_expression(library_size=1e4)` for MSE; posterior means/variances from
  `module.inference` (shared `qz_m/qz_v`; private `qt_m_all/qt_v_all`, gated by the package's group mask).
- MSE scales as in C3: modeled-gene scale (normalize to 1e4 within the 1,500 genes, log1p) and canonical scale
  log1p(ρ·1e4 · L_model / L_full). Note that multiGroupVI uses a latent (encoded) library size; the conversion uses the
  observed modeled-gene total L_model.
- Probes: `LogisticRegression(class_weight='balanced')`, canonical train → test; study and 3-class condition with
  max_iter 1500 (as in script 20); within-study treatment (cisplatin vs PBS in GSE216146; doxorubicin vs control in
  GSE271055) with max_iter 800 (`latent_probes.csv` definition).
- Utilization: active unit = Var[posterior mean] > 0.01, per-dimension KL to N(0, I), on held-out test cells (shared:
  all test and treated test cells; each private block: test cells of its own group).

## Results (`C4_multigroupvi_metrics.csv`)

| seed | canonical MSE | modeled-gene MSE | fit time (s) |
|---|---|---|---|
| 0 | 0.245305 | 1.476459 | 1077 |
| 1 | 0.242906 | 1.440070 | 1079 |
| 2 | 0.248022 | 1.471026 | 1077 |

Reference (C3, same test cells): MC-ContrastiveVI 0.168723 / 0.904167; scVI 0.206108 / 1.080648.

### Probes (balanced accuracy; seeds 0 / 1 / 2)

| representation | study | 3-class condition | cisplatin vs PBS (GSE216146) | doxorubicin vs control (GSE271055) |
|---|---|---|---|---|
| shared (10) | 0.983 / 0.974 / 0.973 | 0.631 / 0.637 / 0.625 | 0.637 / 0.654 / 0.641 | 0.562 / 0.563 / 0.561 |
| group-specific, gated (30) | 0.844 / 0.853 / 0.816 | 0.770 / 0.773 / 0.772 | 0.950 / 0.963 / 0.953 | 0.901 / 0.894 / 0.881 |
| group-specific, ungated (30) | 0.868 / 0.889 / 0.863 | 0.616 / 0.599 / 0.619 | 0.821 / 0.835 / 0.834 | 0.574 / 0.577 / 0.559 |

### Cross-seed linear CKA on test cells (`C4_multigroupvi_cross_seed_cka.csv`)

| representation | 0–1 | 0–2 | 1–2 |
|---|---|---|---|
| shared | 0.755 | 0.765 | 0.767 |
| group-specific, gated | 0.364 | 0.415 | 0.440 |
| group-specific, ungated | 0.306 | 0.381 | 0.359 |

### Latent utilization (`C4_multigroupvi_usage_per_dimension.csv`; seeds 0 / 1 / 2)

| block | population | active units | dims KL > 0.01 | total KL (nats) |
|---|---|---|---|---|
| shared (10 dims) | all test | 10 / 10 / 10 | 10 / 10 / 10 | 12.76 / 13.27 / 12.98 |
| shared | treated test | 10 / 10 / 10 | 10 / 10 / 10 | 12.50 / 13.04 / 12.76 |
| private control (10) | control test | 7 / 9 / 8 | 10 / 10 / 10 | 6.93 / 6.67 / 7.52 |
| private cisplatin (10) | cisplatin test | 10 / 10 / 10 | 10 / 10 / 10 | 7.96 / 8.81 / 8.33 |
| private doxorubicin (10) | doxorubicin test | 7 / 6 / 4 | 10 / 10 / 10 | 6.49 / 5.89 / 5.95 |

`epochs_run` in the metrics file is 401 because it records `trainer.current_epoch + 1` after the final epoch;
training ran for `max_epochs = 400`.
