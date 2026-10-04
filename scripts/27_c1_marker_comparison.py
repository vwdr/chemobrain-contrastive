"""Task C1 step 3: overlap of consensus attribution lists with published cell-type marker sets.

Universe: the 1,500 canonical genes. Marker symbols are matched case-insensitively.
For each consensus list (primary: frozen lists in scripts/19_gene_pattern_audit.py; secondary:
fresh Phase 2 consensus) and each marker set: overlap count, one-sided hypergeometric P
(P[X >= k]), and BH q across the marker sets tested for that list.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

from src.peerreview import core  # noqa: E402

OUT = core.PR_OUT


def main():
    universe = pd.read_csv(core.CORR_OUT / 'benchmark_gene_universe.csv')['gene'].astype(str).tolist()
    N = len(universe); assert N == 1500
    lower = {g.lower(): g for g in universe}
    mk = pd.read_csv(OUT / 'C1_marker_list.csv')
    mk['gene_in_universe'] = mk.gene.astype(str).str.lower().map(lower)
    sets = {}
    for ct, g in mk.groupby('cell_type'):
        sets[(ct, 'all_sources')] = set(g.gene_in_universe.dropna())
        if g.source.nunique() > 1:
            for src, gg in g.groupby('source'):
                short = src.split(' (')[0].split(';')[0]
                sets[(ct, short)] = set(gg.gene_in_universe.dropna())
    size_rows = [dict(cell_type=ct, source=s, n_genes_listed=int(mk[(mk.cell_type == ct)].gene.nunique()) if s == 'all_sources' else None,
                      n_in_universe=len(v)) for (ct, s), v in sets.items()]
    pd.DataFrame(size_rows).to_csv(OUT / 'C1_marker_set_sizes.csv', index=False)

    g19 = core.load_script('gene_pattern_audit', '19_gene_pattern_audit.py')
    fresh = json.loads((OUT / 'A1_fresh_consensus.json').read_text())
    lists = {('frozen_paper', a): v for a, v in g19.CONSENSUS.items()}
    lists.update({('fresh_phase2', a): v for a, v in fresh.items()})
    rows = []
    for (ver, axis), genes in lists.items():
        sel = set(genes); assert sel <= set(universe)
        fam = []
        for (ct, src), ref in sets.items():
            k = len(sel & ref)
            p = float(hypergeom.sf(k - 1, N, len(ref), len(sel))) if k else 1.0
            fam.append(dict(consensus_version=ver, axis=axis, cell_type=ct, marker_source=src, n_consensus=len(sel),
                            n_markers_in_universe=len(ref), overlap=k,
                            expected_overlap=len(sel) * len(ref) / N,
                            fold_enrichment=(k / (len(sel) * len(ref) / N)) if len(ref) else np.nan,
                            hypergeom_p=p, overlap_genes='|'.join(sorted(sel & ref))))
        prim = [r for r in fam if r['marker_source'] == 'all_sources']
        q = multipletests([r['hypergeom_p'] for r in prim], method='fdr_bh')[1]
        for r, qq in zip(prim, q):
            r['bh_q_across_cell_types'] = float(qq)
        rows.extend(fam)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / 'C1_marker_overlap.csv', index=False)

    ann = []
    for (ver, axis), genes in lists.items():
        for rank, g in enumerate(genes, 1):
            member = sorted(ct for (ct, s), ref in sets.items() if s == 'all_sources' and g in ref)
            ann.append(dict(consensus_version=ver, axis=axis, rank=rank, gene=g, n_marker_sets=len(member),
                            marker_sets='|'.join(member)))
    pd.DataFrame(ann).to_csv(OUT / 'C1_consensus_gene_annotation.csv', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70); pd.set_option('display.max_rows', 300)
    p = res[res.marker_source == 'all_sources'][['consensus_version', 'axis', 'cell_type', 'n_consensus', 'n_markers_in_universe',
                                                  'overlap', 'expected_overlap', 'hypergeom_p', 'bh_q_across_cell_types', 'overlap_genes']]
    print(p.to_string(index=False))
    print(res[res.marker_source != 'all_sources'][['consensus_version', 'axis', 'cell_type', 'marker_source', 'n_markers_in_universe', 'overlap', 'hypergeom_p']].to_string(index=False))
    print(pd.DataFrame(size_rows).to_string(index=False))


if __name__ == '__main__':
    main()
