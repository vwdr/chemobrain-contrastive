# Environment notes (Run A)

- `.venv` (main): unchanged from the previous run (`analysis/peerreview_20261003/environment.json`). It contains
  scvi-tools 1.5.1 from `requirements-peerreview-baselines.txt`, which was installed into `.venv` in the previous run.
- `.venv-baselines`: named in the Run A instructions but never created; the previous run installed the baseline
  requirements into `.venv`. All scvi-tools 1.5.1 work in this run (contrastiveVI, R4) uses `.venv`. No new venv created.
- `.venv-mgvi`: unchanged (`analysis/peerreview_20261003/C4_mgvi_pip_freeze.txt`).
