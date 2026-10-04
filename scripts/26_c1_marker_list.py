"""Task C1 step 1-2: build the cell-type marker list from published sources only.

All inputs were downloaded on 2026-10-03 into runs/peerreview_20261003/c1_sources/ (URLs and
SHA-256 recorded in analysis/peerreview_20261003/C1_marker_sources.json). No gene was added
from memory. Gene symbols are stored exactly as given by each source; matching to the mouse
1,500-gene universe is done case-insensitively in scripts/27_c1_marker_comparison.py.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
SRC = R / 'runs' / 'peerreview_20261003' / 'c1_sources'
OUT = R / 'analysis' / 'peerreview_20261003'

PANGLAO_URL = 'https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz'
PANGLAO_TYPES = {
    'Microglia': 'Microglia', 'Macrophages': 'Macrophages', 'Astrocytes': 'Astrocytes',
    'Oligodendrocytes': 'Oligodendrocytes', 'Oligodendrocyte progenitor cells': 'OPCs', 'Neurons': 'Neurons',
    'Endothelial cells': 'Endothelial cells', 'Pericytes': 'Pericytes', 'Fibroblasts': 'Fibroblasts',
    'Ependymal cells': 'Ependymal cells', 'Choroid plexus cells': 'Choroid plexus cells'}
OCHOCKA = ('Ochocka N et al. (2021) Nat Commun, doi:10.1038/s41467-021-21407-w; Supplementary Data 1 '
           '(41467_2021_21407_MOESM3_ESM.xlsx), BAM cluster top-30 DEGs, intersection of sheets "female control" '
           'and "male control"; retrieved via https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7895824/supplementaryFiles')
ZEISEL = ('Zeisel A et al. (2018) Cell, doi:10.1016/j.cell.2018.06.021; Table S4 (mmc4.xlsx, '
          'markers_spec_selec_rob), "Marker" rows for clusters {clusters}; '
          'https://ars.els-cdn.com/content/image/1-s2.0-S009286741830789X-mmc4.xlsx')
DANI = ('Dani N et al. (2021) Cell, doi:10.1016/j.cell.2021.04.003; Table S1 (mmc1.xlsx), sheet '
        '"Cell Type", column "Epithelial"; https://ars.els-cdn.com/content/image/1-s2.0-S0092867421004384-mmc1.xlsx')
BAKKEN = ('Bakken TE et al. (2018) PLoS One, doi:10.1371/journal.pone.0209648; S2 Table '
          '(pone.0209648.s007.xlsx), genes higher in nuclei: logFC < -log2(1.5) and adj.P.Val < 0.05 (159 genes, as '
          'stated in the paper); retrieved via https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6306246/supplementaryFiles')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    rows = []
    # 1. PanglaoDB (mouse entries)
    p = pd.read_csv(SRC / 'PanglaoDB_markers_27_Mar_2020.tsv.gz', sep='\t')
    p = p[p.species.astype(str).str.contains('Mm')]
    for src_type, label in PANGLAO_TYPES.items():
        genes = p[p['cell type'] == src_type]['official gene symbol'].astype(str).unique()
        assert len(genes) > 0, src_type
        for g in genes:
            rows.append(dict(cell_type=label, gene=g, source=f'PanglaoDB markers 27 Mar 2020 (Franzen O et al. 2019, Database, '
                             f'doi:10.1093/database/baz046), mouse entries, cell type "{src_type}"; {PANGLAO_URL}'))
    # 2. Border-associated macrophages
    f = pd.read_excel(SRC / '41467_2021_21407_MOESM3_ESM.xlsx', sheet_name='female control', header=1)
    m = pd.read_excel(SRC / '41467_2021_21407_MOESM3_ESM.xlsx', sheet_name='male control', header=1)
    fb = list(f[f.cluster == 'BAM']['Gene.name'].astype(str)); mb = set(m[m.cluster == 'BAM']['Gene.name'].astype(str))
    for g in [g for g in fb if g in mb]:
        rows.append(dict(cell_type='Border-associated macrophages', gene=g, source=OCHOCKA))
    z = pd.read_excel(SRC / 'zeisel2018_mmc4.xlsx', header=None)

    def zeisel(clusters):
        out = []
        for c in clusters:
            r = z[(z[0] == c) & (z[1] == 'Marker')]
            assert len(r) == 1, c
            out += [str(x) for x in r.iloc[0, 2:7].tolist() if pd.notna(x)]
        return list(dict.fromkeys(out))
    for g in zeisel(['PVM1', 'PVM2']):
        rows.append(dict(cell_type='Border-associated macrophages', gene=g, source=ZEISEL.format(clusters='PVM1, PVM2 (perivascular macrophages)')))
    # 3. Meningeal / perivascular fibroblasts (vascular leptomeningeal cells)
    for g in zeisel(['VLMC1', 'VLMC2', 'ABC']):
        rows.append(dict(cell_type='Meningeal/perivascular fibroblasts', gene=g,
                         source=ZEISEL.format(clusters='VLMC1, VLMC2, ABC (all described as "Vascular leptomeningeal cells")')))
    # 4. Choroid plexus epithelium
    d = pd.read_excel(SRC / 'dani2021_mmc1.xlsx', sheet_name='Cell Type')
    for g in d['Epithelial'].dropna().astype(str).unique():
        rows.append(dict(cell_type='Choroid plexus epithelium', gene=g, source=DANI))
    for g in zeisel(['CHOR']):
        rows.append(dict(cell_type='Choroid plexus epithelium', gene=g, source=ZEISEL.format(clusters='CHOR (choroid plexus epithelial cells)')))
    # 5. Nucleus-enriched transcripts
    b = pd.read_excel(SRC / 'pone.0209648.s007.xlsx')
    nb = b[(b.logFC < -np.log2(1.5)) & (b['adj.P.Val'] < 0.05)]
    assert len(nb) == 159, len(nb)
    for g in nb.gene.astype(str):
        rows.append(dict(cell_type='Nucleus-enriched transcripts', gene=g, source=BAKKEN))

    df = pd.DataFrame(rows).drop_duplicates()
    df.to_csv(OUT / 'C1_marker_list.csv', index=False)
    files = ['PanglaoDB_markers_27_Mar_2020.tsv.gz', '41467_2021_21407_MOESM3_ESM.xlsx', 'zeisel2018_mmc4.xlsx',
             'dani2021_mmc1.xlsx', 'pone.0209648.s007.xlsx']
    sources = {
        'retrieved': '2026-10-03',
        'files': {fn: {'sha256': sha(SRC / fn), 'bytes': (SRC / fn).stat().st_size} for fn in files},
        'symbol_matching': 'case-insensitive match of source symbols to the mouse 1,500-gene universe (PanglaoDB and '
                           'Bakken S2 symbols are upper-case)',
        'not_used_after_inspection': {
            'Van Hove H et al. 2019 Nat Neurosci, Supplementary Tables 2-3': 'statistics and dissociation-gene lists, no marker table',
            'Pietila R et al. 2023 Neuron, mmc2-mmc4': 'GO tables and top-variable genes, no marker table',
            'DeSisto J et al. 2020 Dev Cell, mmc2-mmc6': 'markers contrast meningeal sub-clusters against each other, not fibroblasts against other cell types',
            'Vanlandewijck M et al. 2018 Nature, Supplementary Table 3 / Source data Fig. 1': 'mural-cell and endothelial zonation genes, no fibroblast-vs-all table',
            'Dani N et al. 2021 Table S1 via PMC (NIHMS1690857-supplement-8.xlsx)': 'PMC returned a proof-of-work bot challenge; the same table was obtained from the publisher CDN (mmc1.xlsx)',
        },
        'counts': df.groupby(['cell_type', 'source']).size().rename('n').reset_index().to_dict('records'),
    }
    (OUT / 'C1_marker_sources.json').write_text(json.dumps(sources, indent=2) + '\n')
    print(df.groupby('cell_type').gene.nunique().to_string())


if __name__ == '__main__':
    sys.exit(main())
