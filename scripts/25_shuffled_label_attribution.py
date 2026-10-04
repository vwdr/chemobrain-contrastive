"""Task C2: shuffled-label attribution.

Treatment labels are permuted at the cell level within each study (GSE216146: control/cisplatin;
GSE271055: control/doxorubicin), separately within the training, validation and test cells, with
shuffle seed 20261003. Class proportions within each study x split are preserved. The full model
is then fitted for seeds 0-2 exactly as in Phase 2 and the same integrated-gradient consensus is
derived.

Usage:
  python scripts/25_shuffled_label_attribution.py --fit --seed 0 [--threads=3]
  python scripts/25_shuffled_label_attribution.py --analyze
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure, strip_threads_arg  # noqa: E402
N_THREADS = configure(3)

import argparse  # noqa: E402
import json  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

from src.peerreview import core, fitjob  # noqa: E402

torch.set_num_threads(N_THREADS)
IN = core.PR_RUN / 'canonical_inputs'
RUN = core.PR_RUN / 'shuffled_labels'
OUT = core.PR_OUT
SHUFFLE_SEED = 20261003


def shuffled_inputs():
    z = np.load(IN / 'inputs.npz'); inp = {k: z[k] for k in z.files}
    d = inp['d'].copy(); b = inp['b']; rng = np.random.default_rng(SHUFFLE_SEED)
    for split in ('train', 'val', 'test'):
        ix = inp[split]
        for study in (0, 1):
            jj = ix[b[ix] == study]
            d[jj] = rng.permutation(inp['d'][jj])
    inp['d_real'] = inp['d']; inp['d'] = d
    return inp


def fit(seed):
    inp = shuffled_inputs()
    RUN.mkdir(parents=True, exist_ok=True)
    if not (RUN / 'shuffled_labels.npz').exists():
        np.savez_compressed(RUN / 'shuffled_labels.npz', d_shuffled=inp['d'], d_real=inp['d_real'])
    fitjob.run_fit(inp, seed, RUN, f'shuffled_{seed}', save_ig_arrays=True, extra={'threads': N_THREADS})


def overlap_row(name_a, A, name_b, B, N=1500):
    A, B = set(A), set(B); k = len(A & B)
    p = float(hypergeom.sf(k - 1, N, len(B), len(A))) if k else 1.0
    return dict(list_a=name_a, list_b=name_b, n_a=len(A), n_b=len(B), n_overlap=k, jaccard=core.jaccard(A, B),
                hypergeom_p=p, overlap_genes='|'.join(sorted(A & B)))


def analyze():
    inp = shuffled_inputs(); d, dr, b, tr, te = inp['d'], inp['d_real'], inp['b'], inp['train'], inp['test']
    # label design audit
    aud = []
    for split in ('train', 'val', 'test'):
        ix = inp[split]
        for study in (0, 1):
            jj = ix[b[ix] == study]
            aud.append(dict(split=split, study=['GSE216146', 'GSE271055'][study], n=len(jj),
                            n_treated_real=int((dr[jj] > 0).sum()), n_treated_shuffled=int((d[jj] > 0).sum()),
                            fraction_label_changed=float((dr[jj] != d[jj]).mean()),
                            n_shuffled_treated_that_are_real_treated=int(((d[jj] > 0) & (dr[jj] > 0)).sum())))
    pd.DataFrame(aud).to_csv(OUT / 'C2_label_shuffle_design.csv', index=False)
    mets = pd.DataFrame([json.loads((RUN / f'shuffled_{s}_metrics.json').read_text()) for s in (0, 1, 2)])
    mets.to_csv(OUT / 'C2_shuffled_fit_metrics.csv', index=False)
    rankings, genes = fitjob.load_ig_rankings(RUN, [f'shuffled_{s}' for s in (0, 1, 2)])
    shuf, counts = core.consensus_lists(rankings, genes)
    core.write_json(OUT / 'C2_shuffled_consensus.json', shuf)
    rows = []
    for axis, gl in shuf.items():
        for rank, g in enumerate(gl, 1):
            rows.append(dict(axis=axis, rank_by_mean_abs_ig=rank, gene=g, n_top50_of_6=int(counts[axis][list(genes).index(g)])))
    pd.DataFrame(rows).to_csv(OUT / 'C2_shuffled_consensus_genes.csv', index=False)
    fresh = json.loads((OUT / 'A1_fresh_consensus.json').read_text())
    g19 = core.load_script('gene_pattern_audit', '19_gene_pattern_audit.py'); frozen = g19.CONSENSUS
    ov = []
    for axis in ('shared', 'doxorubicin', 'cisplatin'):
        ov.append(dict(comparison='primary_vs_fresh_real', axis=axis, **overlap_row('shuffled', shuf[axis], 'fresh_real', fresh[axis])))
        ov.append(dict(comparison='secondary_vs_frozen_paper', axis=axis, **overlap_row('shuffled', shuf[axis], 'frozen', frozen[axis])))
    # cross-block overlaps of the shuffled lists with each real list (all 3x3) for completeness
    for a in ('shared', 'doxorubicin', 'cisplatin'):
        for c in ('shared', 'doxorubicin', 'cisplatin'):
            if a != c:
                ov.append(dict(comparison='cross_block_vs_fresh_real', axis=f'{a}|{c}', **overlap_row('shuffled_' + a, shuf[a], 'fresh_' + c, fresh[c])))
    ov = pd.DataFrame(ov); ov.to_csv(OUT / 'C2_overlap.csv', index=False)

    universe = set(pd.read_csv(core.CORR_OUT / 'benchmark_gene_universe.csv')['gene'].astype(str))
    allr = []
    for list_name, lists in [('shuffled', shuf), ('frozen_paper', frozen), ('fresh_real', fresh)]:
        for axis, gl in lists.items():
            for lib, fn in g19.LIBRARIES.items():
                for r in g19.enrichment(set(gl), universe, lib, core.R / 'data' / 'evidence' / fn):
                    r.update(list=list_name, axis=axis); allr.append(r)
    en = pd.DataFrame(allr); en.to_csv(OUT / 'C2_enrichment_all.csv.gz', index=False, compression='gzip')
    summ = (en.assign(sig=en.q_value < 0.05).groupby(['list', 'axis', 'library'], as_index=False)
            .agg(n_tested=('term', 'size'), n_fdr_005=('sig', 'sum'), minimum_q=('q_value', 'min')))
    summ.to_csv(OUT / 'C2_enrichment_summary.csv', index=False)
    top = (en.sort_values(['list', 'axis', 'library', 'q_value', 'p_value', 'term'])
           .groupby(['list', 'axis', 'library'], group_keys=False).head(5))
    top.to_csv(OUT / 'C2_top_enrichment_terms.csv', index=False)
    # the top-ranked M8 term for the paper's shared list, looked up in every list
    m8 = en[(en.library == 'M8_cell_type')]
    top_term = m8[(m8.list == 'frozen_paper') & (m8.axis == 'shared')].sort_values(['q_value', 'p_value', 'term']).iloc[0].term
    look = m8[m8.term == top_term][['list', 'axis', 'term', 'overlap', 'set_size_in_universe', 'list_size',
                                     'fold_enrichment', 'p_value', 'q_value', 'genes']]
    look.to_csv(OUT / 'C2_paper_top_m8_term_lookup.csv', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 80)
    print(mets[['seed', 'test_mse', 'shared_fraction', 'shared_fraction_doxorubicin', 'shared_fraction_cisplatin',
                'shared_active_units_var_gt_0p01', 'dox_active_units_var_gt_0p01', 'cis_active_units_var_gt_0p01',
                'shared_mean_total_kl_nats', 'dox_mean_total_kl_nats', 'cis_mean_total_kl_nats']].to_string(index=False))
    print(pd.DataFrame(aud).to_string(index=False))
    print(ov.drop(columns='overlap_genes').to_string(index=False))
    print(summ.to_string(index=False)); print('TOP M8 TERM (paper shared):', top_term); print(look.to_string(index=False))
    print({k: len(v) for k, v in shuf.items()}); print(shuf)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--fit', action='store_true'); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--analyze', action='store_true'); a = ap.parse_args(strip_threads_arg(sys.argv[1:]))
    if a.fit:
        fit(a.seed)
    if a.analyze:
        analyze()
