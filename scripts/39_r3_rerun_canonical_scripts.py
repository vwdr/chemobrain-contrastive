"""R3: re-run canonical analysis scripts unmodified on the fresh fits, with redirected paths.

The source of scripts/09_interpretation.py, scripts/11_diagnostics.py and scripts/16_pseudobulk_paired.py is read,
only the path expressions that point to the protected folders are replaced (each replacement is asserted to occur
exactly once), and the code is executed with ``__file__`` set to the original script path. The scripts' own logic,
random seeds and thresholds are unchanged.

Staging folder runs/peerreview_20261004/r3_canonical/ holds copies of:
  inputs.npz, cells.csv (regenerated canonical inputs, previous run), full_{0,1,2}.pt and full_{0,1,2}_latents.npz
  (fresh canonical fits, previous run), {no_hsic,no_gating,gaussian_vae}_{0,1,2}_latents.npz (R3 comparison fits).

Usage:
  python scripts/39_r3_rerun_canonical_scripts.py --script 09
  python scripts/39_r3_rerun_canonical_scripts.py --script 11
  python scripts/39_r3_rerun_canonical_scripts.py --script 16
"""
import argparse
import shutil
import sys
from pathlib import Path

R = Path(__file__).resolve().parents[1]
PREV = R / 'runs' / 'peerreview_20261003'
STAGE = R / 'runs' / 'peerreview_20261004' / 'r3_canonical'
COMP = R / 'runs' / 'peerreview_20261004' / 'comparison'
OUTS = {'09': R / 'analysis' / 'peerreview_20261004' / 'R3_09_interpretation',
        '11': R / 'analysis' / 'peerreview_20261004' / 'R3_11_diagnostics',
        '16': R / 'runs' / 'peerreview_20261004' / 'pseudobulk_task4'}
STAGE_REL = "R/'runs'/'peerreview_20261004'/'r3_canonical'"


def stage(need_ablations):
    STAGE.mkdir(parents=True, exist_ok=True)
    pairs = [(PREV / 'canonical_inputs' / 'inputs.npz', 'inputs.npz'), (PREV / 'canonical_inputs' / 'cells.csv', 'cells.csv')]
    for s in (0, 1, 2):
        pairs += [(PREV / 'canonical' / f'full_{s}.pt', f'full_{s}.pt'), (PREV / 'canonical' / f'full_{s}_latents.npz', f'full_{s}_latents.npz')]
        if need_ablations:
            for m in ('no_hsic', 'no_gating', 'gaussian_vae'):
                pairs.append((COMP / f'{m}_{s}_latents.npz', f'{m}_{s}_latents.npz'))
    for src, name in pairs:
        if not src.exists():
            raise FileNotFoundError(src)
        if not (STAGE / name).exists():
            shutil.copy2(src, STAGE / name)


def run(script):
    fname = {'09': '09_interpretation.py', '11': '11_diagnostics.py', '16': '16_pseudobulk_paired.py'}[script]
    src = (R / 'scripts' / fname).read_text(encoding='utf-8')
    out = OUTS[script]; out.mkdir(parents=True, exist_ok=True)
    out_rel = "R/" + "/".join(f"'{p}'" for p in out.relative_to(R).parts)
    if script in ('09', '11'):
        reps = [("RUN=R/'runs/corrected_20260920'", f"RUN={STAGE_REL}"), ("T=R/'analysis/corrected_20260920'", f"T={out_rel}")]
    else:
        reps = [('T = R / "analysis" / "corrected_20260920"', f'T = {out_rel}')]
    for old, new in reps:
        assert src.count(old) == 1, (fname, old, src.count(old))
        src = src.replace(old, new)
    for forbidden in ('corrected_20260920',):
        assert forbidden not in src.replace('# ', ''), f'{fname}: protected path still referenced after redirection'
    sys.argv = [str(R / 'scripts' / fname)]
    g = {'__name__': '__main__', '__file__': str(R / 'scripts' / fname)}
    exec(compile(src, str(R / 'scripts' / fname), 'exec'), g)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--script', required=True); a = ap.parse_args()
    if a.script in ('09', '11'):
        stage(need_ablations=(a.script == '11'))
    run(a.script)
    print('finished', a.script)
