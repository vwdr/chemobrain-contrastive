#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON=${PYTHON:-python}
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=2
mkdir -p runs/corrected_20260920/logs analysis/corrected_20260920
"$PYTHON" scripts/05_download_research_data.py
"$PYTHON" scripts/05_recover_cohort.py
"$PYTHON" scripts/06_validation.py --prepare
for seed in 0 1 2; do
  for mode in full no_hsic no_gating gaussian_vae; do
    "$PYTHON" scripts/06_validation.py --mode "$mode" --seed "$seed" > "runs/corrected_20260920/logs/${mode}_${seed}.log" 2>&1
  done
  "$PYTHON" scripts/12_nb_sensitivity.py --seed "$seed" > "runs/corrected_20260920/logs/nb_${seed}.log" 2>&1
done
"$PYTHON" scripts/06_validation.py --pca
"$PYTHON" scripts/08_annotation.py
"$PYTHON" scripts/09_interpretation.py
"$PYTHON" scripts/10_pseudobulk.py
"$PYTHON" scripts/11_diagnostics.py
"$PYTHON" tests/test_research_invariants.py
"$PYTHON" scripts/13_collect_results.py
