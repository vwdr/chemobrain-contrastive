"""Figures 7 and 8 from analysis/experiments/master.csv (run experiments/collect.py first) and fsh_noise_null_draws.csv.

  fig7_freebits.pdf   active units, recovery (AUROC and top-50 enrichment, always together) and shared_fraction against
                      free_bits; 2 x 2 panels, one point per seed, one colour per configuration, line = mean over seeds.
                      (The null configuration has no injected genes, so it is absent from the recovery panels.)
  fig8_nullcalib.pdf  (a) shared_fraction of the 30 null fits (exp2) with the injected archived fits overlaid;
                      (b) the simulated pure-noise F_sh null (exp6) on the same x axis, the 30 null fits as an outline and the
                      three observed null fits (0.705, 0.802, 0.869) as ticks.

Needs experiments/paperstyle.py exposing apply(), panel(ax, letter, title) and the colour constants NAVY, BRICK, FOREST,
PLUM, OCHRE, GREY (and METHOD / BLOCK dicts, unused here). 300 dpi PDFs, 7.16 inches wide, nothing is trained.

Usage:
  python experiments/make_figures.py [--out analysis/experiments]
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import argparse  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from experiments import paperstyle as ps  # noqa: E402

R = Path(__file__).resolve().parents[1]
WIDTH = 7.16
DPI = 300
OBSERVED_NULL = [0.705, 0.802, 0.869]


def config_style():
    return {'null': (ps.GREY, 'null (no injection)'),
            'shared_d1_rf1.0': (ps.NAVY, 'shared, delta 1'),
            'shared_d2_rf1.0': (ps.BRICK, 'shared, delta 2')}


def level_positions(levels):
    return {lv: i for i, lv in enumerate(levels)}


def dose_panel(ax, df, ycol, ylabel, levels, style, jitter_seed=0):
    pos = level_positions(levels); rng = np.random.default_rng(jitter_seed)
    for cfg, (colour, label) in style.items():
        d = df[(df.config == cfg)].dropna(subset=[ycol])
        if d.empty:
            continue
        x = d.free_bits.map(pos).to_numpy(dtype=float)
        ax.scatter(x + rng.uniform(-0.12, 0.12, len(x)), d[ycol], s=9, color=colour, alpha=0.75, linewidths=0, label=label, zorder=3)
        m = d.groupby('free_bits')[ycol].mean().reindex(levels)
        ax.plot([pos[lv] for lv in m.index if not np.isnan(m[lv])], m.dropna().to_numpy(), color=colour, lw=1.2, zorder=2)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f'{lv:g}' for lv in levels])
    ax.set_xlabel('free_bits (nats per dimension)')
    ax.set_ylabel(ylabel)


def fig7(master, out):
    e1 = master[master.experiment == 'exp1']
    if e1.empty:
        raise SystemExit('no exp1 rows in master.csv; run exp1 and experiments/collect.py first')
    levels = sorted(e1.free_bits.dropna().unique())
    style = config_style()
    fig, axs = plt.subplots(2, 2, figsize=(WIDTH, 4.6), constrained_layout=True)
    axes = axs.ravel()
    dose_panel(axes[0], e1, 'active_units', 'active units (of 16)', levels, style)
    axes[0].set_ylim(bottom=0)
    dose_panel(axes[1], e1, 'recovery_auroc', 'recovery AUROC (injected genes)', levels, style)
    axes[1].axhline(0.5, color=ps.GREY, lw=0.6, ls=':', zorder=1)
    dose_panel(axes[2], e1, 'recovery_top50_enrichment', 'recovery top-50 enrichment (x random)', levels, style)
    axes[2].axhline(1.0, color=ps.GREY, lw=0.6, ls=':', zorder=1)
    dose_panel(axes[3], e1, 'shared_fraction', r'shared fraction $F_{sh}$', levels, style)
    axes[3].set_ylim(0, 1)
    for ax, (letter, title) in zip(axes, [('a', 'Active units'), ('b', 'Recovery (AUROC)'), ('c', 'Recovery (top-50 enrichment)'),
                                          ('d', 'Shared fraction')]):
        ps.panel(ax, letter, title)
    h, l = axes[3].get_legend_handles_labels()
    fig.legend(h, l, loc='outside lower center', ncol=3, frameon=False)
    path = out / 'fig7_freebits.pdf'
    fig.savefig(path, dpi=DPI); plt.close(fig)
    return path


def hist(ax, v, bins, colour, label, **kw):
    ax.hist(v, bins=bins, density=True, color=colour, label=label, **kw)


def fig8(master, draws, out):
    null = master[master.experiment == 'exp2'].shared_fraction.dropna().to_numpy()
    inj = master[(master.experiment == 'archive_semisynthetic') & (master.scenario != 'null')].shared_fraction.dropna().to_numpy()
    if len(null) == 0:
        raise SystemExit('no exp2 rows in master.csv; run exp2 and experiments/collect.py first')
    allv = np.concatenate([null, inj, draws if draws is not None else [], OBSERVED_NULL])
    lo, hi = max(0.0, np.floor(allv.min() * 20) / 20 - 0.05), min(1.0, np.ceil(allv.max() * 20) / 20 + 0.05)
    bins = np.linspace(lo, hi, 31)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 2.7), sharex=True, constrained_layout=True)
    a, b = axes
    hist(a, null, bins, ps.GREY, f'null fits (n = {len(null)})', alpha=0.85)
    if len(inj):
        hist(a, inj, bins, ps.BRICK, f'injected fits (n = {len(inj)})', alpha=0.45)
    a.set_ylabel('density'); a.set_xlabel(r'shared fraction $F_{sh}$')
    a.legend(frameon=False, loc='best')
    if draws is not None:
        hist(b, draws, bins, ps.PLUM, f'simulated noise null (n = {len(draws):,})', alpha=0.6)
    b.hist(null, bins=bins, density=True, histtype='step', color=ps.GREY, lw=1.2, label='null fits')
    for i, o in enumerate(OBSERVED_NULL):
        b.plot([o, o], [0, 0.06], transform=b.get_xaxis_transform(), color=ps.OCHRE, lw=1.4, solid_capstyle='butt',
               label='observed null fits' if i == 0 else None, clip_on=False)
    b.set_xlabel(r'shared fraction $F_{sh}$'); b.set_xlim(lo, hi)
    b.legend(frameon=False, loc='best')
    ps.panel(a, 'a', 'Null vs injected fits'); ps.panel(b, 'b', 'Pure-noise F_sh')
    path = out / 'fig8_nullcalib.pdf'
    fig.savefig(path, dpi=DPI); plt.close(fig)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, default=R / 'analysis' / 'experiments')
    ap.add_argument('--smoke', action='store_true', help='use <out>/smoke')
    args = ap.parse_args(argv)
    if args.smoke:
        args.out = args.out / 'smoke'
    mp = args.out / 'master.csv'
    if not mp.exists():
        raise SystemExit(f'{mp} missing: run experiments/collect.py first')
    master = pd.read_csv(mp, keep_default_na=False, na_values=[''])
    ps.apply()
    plt.rcParams['savefig.bbox'] = 'standard'          # paperstyle sets 'tight', which would crop the 7.16 in width
    print('wrote', fig7(master, args.out))
    dp = args.out / 'fsh_noise_null_draws.csv'
    draws = pd.read_csv(dp)['pooled'].to_numpy() if dp.exists() else None
    if draws is None:
        print(f'note: {dp} not found, panel b of fig8 shows the null fits only')
    print('wrote', fig8(master, draws, args.out))
    return 0


if __name__ == '__main__':
    sys.exit(main())
