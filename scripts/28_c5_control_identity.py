"""Task C5: control-cell identity check on fresh full-model latents (seeds 0-2), GSE216146 only.

1. Balanced logistic classifiers (class_weight='balanced', max_iter=800, as in 11_diagnostics) predict the
   source cell type from ungated posterior means, separately in control (PBS, nonrescue) and cisplatin
   (nonrescue) cells, trained on canonical training cells and scored on canonical test cells.
2. Per-cell-type mean and SD of the ungated block norm in control vs cisplatin cells (all nonrescue
   GSE216146 cells, and test cells only).
Blocks: shared (8 dims); cisplatin-specific drug block (dims 4-7 of the ungated drug means); full ungated
drug block (8 dims, doxorubicin + cisplatin heads).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure  # noqa: E402
configure(3)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import balanced_accuracy_score  # noqa: E402

from src.peerreview import core  # noqa: E402

IN = core.PR_RUN / 'canonical_inputs'; RUN = core.PR_RUN / 'canonical'; OUT = core.PR_OUT


def main():
    z = np.load(IN / 'inputs.npz'); d, b, tr, te = z['d'], z['b'], z['train'], z['test']
    cells = pd.read_csv(IN / 'cells.csv', index_col=0)
    ct = cells.source_cell_type.astype(str).to_numpy(); rescue = cells.rescue.astype(bool).to_numpy()
    acc = []; norms = []; counts = []
    arms = {'control_PBS': (b == 0) & (d == 0) & ~rescue, 'cisplatin': (b == 0) & (d == 2) & ~rescue}
    split = np.full(len(d), 'not_used', dtype=object); split[tr] = 'train'; split[z['val']] = 'validation'; split[te] = 'test'
    for arm, mask in arms.items():
        for s in ('train', 'test'):
            ix = np.where(mask & (split == s))[0]
            for c, n in pd.Series(ct[ix]).value_counts().items():
                counts.append(dict(arm=arm, split=s, cell_type=c, n_cells=int(n)))
    for seed in (0, 1, 2):
        lat = np.load(RUN / f'full_{seed}_latents.npz')
        blocks = {'shared_ungated': lat['raw_shared'], 'cisplatin_drug_block_ungated': lat['raw_drug'][:, 4:],
                  'drug_block_ungated_all8': lat['raw_drug']}
        for bname, rep in blocks.items():
            for arm, mask in arms.items():
                trn = tr[mask[tr]]; tst = te[mask[te]]
                clf = LogisticRegression(max_iter=800, class_weight='balanced').fit(rep[trn], ct[trn])
                pred = clf.predict(rep[tst])
                acc.append(dict(seed=seed, block=bname, arm=arm, n_train=len(trn), n_test=len(tst),
                                n_classes_train=len(np.unique(ct[trn])), n_classes_test=len(np.unique(ct[tst])),
                                balanced_accuracy=float(balanced_accuracy_score(ct[tst], pred)),
                                chance_level=1.0 / len(np.unique(ct[tst]))))
                nr = np.linalg.norm(rep, axis=1)
                for scope, sm in [('all_nonrescue_cells', mask), ('test_cells', mask & (split == 'test'))]:
                    for c in np.unique(ct[sm]):
                        v = nr[sm & (ct == c)]
                        norms.append(dict(seed=seed, block=bname, arm=arm, scope=scope, cell_type=c, n_cells=len(v),
                                          mean_norm=float(v.mean()), sd_norm=float(v.std(ddof=1)) if len(v) > 1 else np.nan))
    acc = pd.DataFrame(acc); norms = pd.DataFrame(norms)
    acc.to_csv(OUT / 'C5_celltype_probe_accuracy.csv', index=False)
    norms.to_csv(OUT / 'C5_block_norms_by_celltype.csv', index=False)
    pd.DataFrame(counts).to_csv(OUT / 'C5_cell_counts.csv', index=False)
    wide = norms[norms.scope == 'all_nonrescue_cells'].pivot_table(index=['block', 'cell_type', 'seed'], columns='arm',
                                                                      values=['mean_norm', 'sd_norm', 'n_cells'])
    wide.columns = [f'{a}_{b}' for a, b in wide.columns]; wide = wide.reset_index()
    wide['cisplatin_minus_control_mean_norm'] = wide.mean_norm_cisplatin - wide.mean_norm_control_PBS
    wide.to_csv(OUT / 'C5_block_norms_control_vs_cisplatin.csv', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 400)
    print(acc.to_string(index=False)); print(wide.to_string(index=False)); print(pd.DataFrame(counts).to_string(index=False))


if __name__ == '__main__':
    main()
