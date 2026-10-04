# Reproduction commands (peer-review runs of 2026-10-03 and 2026-10-04)

Repository: `github.com/vwdr/chemobrain-contrastive`, branch `peerreview-20261003`. Platform used: Windows 11, Git Bash,
Python 3.12.10, CPU only (Intel i5-10400). Run all commands from the repository root. Thread counts are set inside the
scripts (`--threads=N` where supported); at most three fits were run concurrently.

Archived outputs, their SHA-256 values, producing scripts, commands and commits are listed in `MANIFEST.csv`; the zip files
in `ARCHIVE_ZIPS.csv`. Neural-network fits are not bitwise reproducible across hardware and library builds; the archived
checkpoints/latents are the outputs of the runs described here.

## 1. Environments

```
# main environment
py -3.12 -m venv .venv
.venv/Scripts/python -m pip install --upgrade pip
.venv/Scripts/python -m pip install -r requirements-research.txt
.venv/Scripts/python -m pip install -r requirements-task4.txt -r requirements-peerreview-baselines.txt
.venv/Scripts/python -m pip install openpyxl==3.1.5 xlrd==2.0.2      # reading supplementary spreadsheets (Task C1)
# exact versions: pip_freeze_venv.txt

# multiGroupVI environment (scvi-tools 0.18.0)
git clone https://github.com/Genentech/multiGroupVI <dir>; git -C <dir> checkout 3d001dbe6905803195165f936f879a1a166a6fc8
py -3.12 -m venv .venv-mgvi
.venv-mgvi/Scripts/python -m pip install "scvi-tools==0.18.0" "protobuf<=3.20.1" "scanpy>=1.8.1"
.venv-mgvi/Scripts/python -m pip install --extra-index-url https://download.pytorch.org/whl/cpu "setuptools<70" "torch==2.5.1+cpu" \
    "numpy==1.26.4" "scipy<1.14" "pandas<2.2" "jax==0.4.30" "jaxlib==0.4.30" "flax==0.8.5" "optax==0.2.3" "chex==0.1.86" \
    "numpyro==0.15.3" "pyro-ppl==1.9.1" "torchmetrics==0.11.4" "zarr<3" "numcodecs<0.13" "fast-array-utils<1.2" \
    "anndata==0.9.2" "mudata==0.2.3" "scanpy==1.9.8" "scikit-learn<1.6" "matplotlib<3.10" "h5py<3.13" "numba<0.61"
.venv-mgvi/Scripts/python -m pip install "pip==24.0"
.venv-mgvi/Scripts/python -m pip install --no-deps "pytorch-lightning==1.7.7"
.venv-mgvi/Scripts/python -m pip install --no-deps <dir>
# exact versions: pip_freeze_venv-mgvi.txt
```

The instructions for Run A mention `.venv-baselines`; it was never created. scvi-tools 1.5.1 (requirements-peerreview-baselines.txt)
is installed in `.venv` and was used for scVI and contrastiveVI.

## 2. Data

```
.venv/Scripts/python scripts/05_download_research_data.py
sha256sum data/raw/GSE216146_chemo_brain.h5ad.gz   # 541822fb6983a4a07a2a4166479228078373c62153b4b9b237914e490bfaa04d
sha256sum data/raw/GSE271055_RAW.tar               # 4d31e82d33e04b842b09f2744fe9a13b2e49afb2b9bc8ae5091b34da4e1c1ba1
git checkout -- data/evidence/GSE216146.soft data/evidence/GSE271055.soft data/evidence/GSE286221.soft   # only line endings change
.venv/Scripts/python scripts/05_recover_cohort.py
sha256sum data/processed/recovered_counts.h5ad      # 62a5047e74d4ea704f2a10b5b91b35a4836ccefc99186d4c72101db78cbd8ac1
```

MSigDB files (downloaded by the same script): m2 reactome `03fbe5cb…acb5`, m5 GO BP `343192bd…e865` (both as in
`data/evidence/source_file_manifest.csv`), m8 `0b11177c438565361ce7dc917d216cf40dd71b12982d7f117b51b182be6a7788`. Full hashes:
`EXCLUDED_INPUT_HASHES.csv`. `05_recover_cohort.py` rewrites three CSVs in `analysis/corrected_20260920/`; they were byte-identical
to the committed versions.

## 3. Previous run (outputs in analysis/ and runs/peerreview_20261003)

```
.venv/Scripts/python scripts/23_prepare_and_ved_check.py                       # canonical inputs; ved_archive comparison (needs ../ved_archive)
for s in 0 1 2; do .venv/Scripts/python scripts/24_canonical_fits.py --fit --seed $s --threads=2; done
.venv/Scripts/python scripts/24_canonical_fits.py --analyze
for s in 0 1 2; do .venv/Scripts/python scripts/25_shuffled_label_attribution.py --fit --seed $s --threads=3; done
.venv/Scripts/python scripts/26_c1_marker_list.py                              # needs the downloaded source files in runs/peerreview_20261003/c1_sources (see C1_marker_sources.json)
.venv/Scripts/python scripts/27_c1_marker_comparison.py
.venv/Scripts/python scripts/28_c5_control_identity.py
.venv/Scripts/python scripts/25_shuffled_label_attribution.py --analyze
for s in 0 1 2; do .venv/Scripts/python scripts/29_baselines_extended.py --job scvi --seed $s --threads=3; done
for p in cisplatin doxorubicin; do for s in 0 1 2; do .venv/Scripts/python scripts/29_baselines_extended.py --job cvi --pair $p --seed $s --threads=3; done; done
.venv/Scripts/python scripts/29_baselines_extended.py --collect --threads=2
.venv/Scripts/python scripts/30_calibration_pseudobulk.py
.venv/Scripts/python scripts/31_semisynthetic.py --prepare
.venv/Scripts/python tests/test_peerreview_thinning.py                         # prints PASS
.venv/Scripts/python scripts/31_semisynthetic.py --list-jobs grid             # 75 lines "fit <config> <seed>"
.venv/Scripts/python scripts/31_semisynthetic.py --fit --config <config> --seed <seed> --threads=3      # for every grid line
.venv/Scripts/python scripts/31_semisynthetic.py --list-jobs perm             # 200 lines "perm <config> <k>"
.venv/Scripts/python scripts/31_semisynthetic.py --perm --config <config> --perm-index <k> --threads=3  # for every perm line
.venv/Scripts/python scripts/31_semisynthetic.py --score --threads=1
.venv/Scripts/python scripts/35_d1_signal_check.py
.venv/Scripts/python scripts/33_d1_report_tables.py
for s in 0 1 2; do .venv-mgvi/Scripts/python scripts/34_c4_multigroupvi.py --fit --seed $s --threads=3; done
.venv-mgvi/Scripts/python scripts/34_c4_multigroupvi.py --collect
.venv/Scripts/python scripts/32_c6_ttr_distribution.py
```

## 4. Run A (outputs in analysis/ and runs/peerreview_20261004)

```
for s in 0 1 2; do for m in no_hsic no_gating gaussian_vae; do .venv/Scripts/python scripts/36_r3_comparison_fits.py --job mode --mode $m --seed $s; done; \
  .venv/Scripts/python scripts/36_r3_comparison_fits.py --job nb --seed $s; done
.venv/Scripts/python scripts/36_r3_comparison_fits.py --job pca
.venv/Scripts/python scripts/37_r2_baseline_headline.py
.venv/Scripts/python scripts/38_r1_injected_signal.py
.venv/Scripts/python scripts/39_r3_rerun_canonical_scripts.py --script 09
.venv/Scripts/python scripts/39_r3_rerun_canonical_scripts.py --script 16
.venv/Scripts/python scripts/39_r3_rerun_canonical_scripts.py --script 11      # after the comparison fits
for part in comparison consensus enrichment pseudobulk; do .venv/Scripts/python scripts/40_r3_analyses.py --part $part; done
.venv/Scripts/python scripts/45_r3_paper_numbers.py                            # needs ../manuscript_final_text.txt
.venv/Scripts/python scripts/41_r4_prepare.py
for c in null shared_d1_rf1.0 specific_d1_rf1.0 mixed_d1_rf1.0 shared_d2_rf1.0 specific_d2_rf1.0 mixed_d2_rf1.0; do for s in 0 1 2; do
  .venv-mgvi/Scripts/python scripts/42_r4_multigroupvi.py --config $c --seed $s --threads=3
  for st in GSE216146 GSE271055; do .venv/Scripts/python scripts/43_r4_contrastivevi.py --config $c --study $st --seed $s --threads=3; done
done; done
.venv/Scripts/python scripts/44_r4_mccvi_and_collect.py --mccvi
.venv/Scripts/python scripts/44_r4_mccvi_and_collect.py --collect
.venv/Scripts/python scripts/46_r5_archive.py
.venv/Scripts/python scripts/47_r6_index.py
```

## 5. Expected outputs

- Deterministic and checked: the three input SHA-256 values above; the regenerated canonical split and inputs are
  content-identical to the committed inputs (`analysis/peerreview_20261003/A0_inputs_comparison.csv`: index/label arrays identical,
  X max abs difference 9.5e-7); the semi-synthetic design (`D1_gene_sets.csv`, responders, split) and the injected counts
  (`scripts/41_r4_prepare.py` asserts equality with the Phase 4 inputs).
- Model fits, attributions and statistics derived from them: archived values and SHA-256 in `MANIFEST.csv`; refits on other
  hardware are expected to differ numerically.
- `.npz` files are zip containers whose bytes depend on write time; compare array contents rather than file hashes.
