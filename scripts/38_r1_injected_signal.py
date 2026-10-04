"""R1: where the injected semi-synthetic signal went (follow-up on the 75 Phase 4 grid fits; no refits).

1. Within-study probes: LogisticRegression(class_weight='balanced', max_iter=800) trained on training cells and scored
   (balanced accuracy) on test cells, separately within GSE216146 and GSE271055, on bg, raw_shared and raw_drug
   (from the saved *_latents.npz). Labels: (a) pseudo-treated vs pseudo-control; (b) responders vs all other cells of
   the study (responder set of the configuration's responder fraction; for the null configuration (b) is computed with
   both responder sets, labelled rf 1.0 and rf 0.3, although nothing was injected).
2. Reconstructed effect: checkpoint decoded on test cells with canonical gated decoding (06_validation.means).
   For each set (S in each study; A in GSE216146; B in GSE271055) and direction (up/down):
   observed = mean over set genes of [mean(X, test responders of the study) - mean(X, other test cells of the study)];
   reconstructed = the same on the decoder output; ratio = reconstructed / observed. Rows are produced for every set in
   every configuration with a flag `injected` (whether the set was injected in that configuration).
Inputs are read from runs/peerreview_20261003/semisynthetic (read-only).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure  # noqa: E402
N_THREADS = configure(2)

import importlib.util  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import balanced_accuracy_score  # noqa: E402

from src.peerreview import core  # noqa: E402

torch.set_num_threads(N_THREADS)
spec = importlib.util.spec_from_file_location('ss', core.R / 'scripts' / '31_semisynthetic.py')
ss = importlib.util.module_from_spec(spec); _argv = sys.argv; sys.argv = sys.argv[:1]; spec.loader.exec_module(ss); sys.argv = _argv
OUT = core.R / 'analysis' / 'peerreview_20261004'
STUDIES = {0: 'GSE216146', 1: 'GSE271055'}
SET_STUDIES = {'S': (0, 1), 'A': (0,), 'B': (1,)}


def probe(rep, y, tr, te):
    if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
        return np.nan
    c = LogisticRegression(class_weight='balanced', max_iter=800).fit(rep[tr], y[tr])
    return float(balanced_accuracy_score(y[te], c.predict(rep[te])))


def main():
    base = ss.load_base(); gs = pd.read_csv(ss.OUT / 'D1_gene_sets.csv')
    d, b, tr, te, treated = base['d'], base['b'], base['train'], base['test'], base['treated']
    resp_sets = {1.0: base['resp_1p0'], 0.3: base['resp_0p3']}
    prow, rrow = [], []
    for cfg in ss.all_configs():
        scen, delta, rf = ss.parse_config(cfg)
        inp, _ = ss.build_inputs(cfg, base); X = inp['X']
        rf_list = [rf] if scen != 'null' else [1.0, 0.3]
        for seed in (0, 1, 2):
            name = f'{cfg}_s{seed}'
            lat = np.load(ss.RUN / cfg / f'{name}_latents.npz')
            for sb, sname in STUDIES.items():
                trs = tr[b[tr] == sb]; tes = te[b[te] == sb]
                for rep in ('bg', 'raw_shared', 'raw_drug'):
                    y = treated.astype(int)
                    prow.append(dict(config=cfg, scenario=scen, delta=delta, responder_fraction=rf, seed=seed, study=sname, representation=rep,
                                     label='pseudo_treated_vs_pseudo_control', responder_set=None, balanced_accuracy=probe(lat[rep], y, trs, tes)))
                    for r in rf_list:
                        y = resp_sets[r].astype(int)
                        prow.append(dict(config=cfg, scenario=scen, delta=delta, responder_fraction=rf, seed=seed, study=sname, representation=rep,
                                         label='responders_vs_other_cells_of_study', responder_set=r, balanced_accuracy=probe(lat[rep], y, trs, tes)))
            # reconstruction
            m, _ = core.load_model(ss.RUN / cfg / f'{name}.pt')
            with torch.no_grad():
                xt = torch.from_numpy(X[te]); bg, sh, dr, _ = core.means(m, xt, torch.from_numpy(d[te]))
                rec = m.decode(bg, sh, dr, torch.from_numpy(b[te])).numpy()
            Xte = X[te]; bte = b[te]
            for r in rf_list:
                rte = resp_sets[r][te]
                for st, studies in SET_STUDIES.items():
                    injected = st in ss.SETS_BY_SCENARIO[scen]
                    for sb in studies:
                        msk_s = bte == sb; R_ = msk_s & rte; O_ = msk_s & ~rte
                        for direction in ('up', 'down'):
                            cols = gs[(gs.set == st) & (gs.direction == direction)].gene_index.to_numpy()
                            obs = float((Xte[np.ix_(R_, cols)].mean(0) - Xte[np.ix_(O_, cols)].mean(0)).mean())
                            rc = float((rec[np.ix_(R_, cols)].mean(0) - rec[np.ix_(O_, cols)].mean(0)).mean())
                            rrow.append(dict(config=cfg, scenario=scen, delta=delta, responder_fraction=rf, responder_set=r, seed=seed, set=st,
                                             study=STUDIES[sb], direction=direction, injected=injected, n_test_responders=int(R_.sum()),
                                             n_test_other=int(O_.sum()), observed_effect=obs, reconstructed_effect=rc,
                                             ratio_reconstructed_over_observed=rc / obs if obs != 0 else np.nan,
                                             expected_effect_ln=(1 if direction == 'up' else -1) * delta * np.log(2) if injected else 0.0))
        print('done', cfg, flush=True)
    P = pd.DataFrame(prow); Rr = pd.DataFrame(rrow)
    P.to_csv(OUT / 'R1_probes_per_fit.csv', index=False); Rr.to_csv(OUT / 'R1_reconstruction_per_fit.csv', index=False)
    keys = ['config', 'scenario', 'delta', 'responder_fraction']
    ps = (P.groupby(keys + ['study', 'representation', 'label', 'responder_set'], dropna=False, sort=False).balanced_accuracy
          .agg(['mean', 'std', 'min', 'max']).reset_index().assign(quantity='probe_balanced_accuracy'))
    rs = (Rr.groupby(keys + ['responder_set', 'set', 'study', 'direction', 'injected'], dropna=False, sort=False)
          .agg(observed_mean=('observed_effect', 'mean'), observed_sd=('observed_effect', 'std'),
               reconstructed_mean=('reconstructed_effect', 'mean'), reconstructed_sd=('reconstructed_effect', 'std'),
               ratio_mean=('ratio_reconstructed_over_observed', 'mean'), ratio_sd=('ratio_reconstructed_over_observed', 'std'),
               expected_effect_ln=('expected_effect_ln', 'first')).reset_index())
    ps.to_csv(OUT / 'R1_summary_probes.csv', index=False); rs.to_csv(OUT / 'R1_summary_reconstruction.csv', index=False)
    summ = pd.concat([ps.rename(columns={'mean': 'value_mean', 'std': 'value_sd'}).assign(table='probes'),
                      rs.assign(table='reconstruction')], ignore_index=True)
    summ.to_csv(OUT / 'R1_summary.csv', index=False)
    print(len(P), len(Rr))


if __name__ == '__main__':
    main()
