"""R5: build the archive package (not uploaded) in runs/peerreview_20261004/archive/.

Creates one zip per group (each well below 10 GB), MANIFEST.csv (one row per archived file: path inside the
archive, bytes, SHA-256, producing script, exact command, commit of the producing script), pip freezes of the
environments, and EXCLUDED_INPUT_HASHES.csv for the public/regenerable inputs that are not archived. Smoke-test
outputs and staging copies of inputs are not archived. Raw supplementary marker-source files are not archived;
their SHA-256 values are in C1_marker_sources.json.
"""
import csv
import hashlib
import subprocess
import sys
import zipfile
from pathlib import Path

R = Path(__file__).resolve().parents[1]
A = R / 'runs' / 'peerreview_20261004' / 'archive'
P3 = R / 'runs' / 'peerreview_20261003'; P4 = R / 'runs' / 'peerreview_20261004'
AN3 = R / 'analysis' / 'peerreview_20261003'

GROUPS = {
    'canonical_fits.zip': [(P3 / 'canonical', 'canonical')],
    'canonical_inputs.zip': [(P3 / 'canonical_inputs', 'canonical_inputs')],
    'r3_comparison_fits.zip': [(P4 / 'comparison', 'comparison')],
    'c2_shuffled_label_fits.zip': [(P3 / 'shuffled_labels', 'shuffled_labels')],
    'semisynthetic.zip': [(P3 / 'semisynthetic', 'semisynthetic'), (AN3 / 'D1_gene_sets.csv', 'semisynthetic/D1_gene_sets.csv')],
    'baselines.zip': [(P3 / 'baselines', 'c3_baselines'), (P3 / 'multigroupvi', 'c4_multigroupvi'), (P4 / 'r4', 'r4')],
    'marker_list.zip': [(AN3 / 'C1_marker_list.csv', 'C1_marker_list.csv'), (AN3 / 'C1_marker_sources.json', 'C1_marker_sources.json')],
}
EXCLUDE_PARTS = {'smoke'}
EXCLUDE_NAMES = {('comparison', 'inputs.npz'), ('comparison', 'count_inputs.npz')}   # staging copies of canonical_inputs


def producer(inner: str):
    """(script, exact command) for a path inside the archive."""
    parts = inner.split('/'); top = parts[1] if len(parts) > 1 else parts[0]; name = parts[-1]
    if top == 'canonical':
        return 'scripts/24_canonical_fits.py', '.venv/Scripts/python scripts/24_canonical_fits.py --fit --seed {0,1,2} --threads=2'
    if top == 'canonical_inputs':
        return 'scripts/23_prepare_and_ved_check.py', '.venv/Scripts/python scripts/23_prepare_and_ved_check.py'
    if top == 'comparison':
        if name.startswith('negative_binomial') and name.endswith('_latents.npz'):
            return 'scripts/40_r3_analyses.py', '.venv/Scripts/python scripts/40_r3_analyses.py --part comparison'
        if name.startswith('negative_binomial'):
            return 'scripts/36_r3_comparison_fits.py', '.venv/Scripts/python scripts/36_r3_comparison_fits.py --job nb --seed {0,1,2}'
        if name.startswith('pca'):
            return 'scripts/36_r3_comparison_fits.py', '.venv/Scripts/python scripts/36_r3_comparison_fits.py --job pca'
        return 'scripts/36_r3_comparison_fits.py', '.venv/Scripts/python scripts/36_r3_comparison_fits.py --job mode --mode {no_hsic,no_gating,gaussian_vae} --seed {0,1,2}'
    if top == 'shuffled_labels':
        return 'scripts/25_shuffled_label_attribution.py', '.venv/Scripts/python scripts/25_shuffled_label_attribution.py --fit --seed {0,1,2} --threads=3'
    if top == 'semisynthetic':
        if name == 'D1_gene_sets.csv' or name.startswith('base'):
            return 'scripts/31_semisynthetic.py', '.venv/Scripts/python scripts/31_semisynthetic.py --prepare'
        if len(parts) > 2 and parts[2].startswith('perm_'):
            return 'scripts/31_semisynthetic.py', f'.venv/Scripts/python scripts/31_semisynthetic.py --perm --config {parts[2][5:]} --perm-index K --threads=3'
        return 'scripts/31_semisynthetic.py', f'.venv/Scripts/python scripts/31_semisynthetic.py --fit --config {parts[2] if len(parts) > 2 else "?"} --seed S --threads=3'
    if top == 'c3_baselines':
        return 'scripts/29_baselines_extended.py', '.venv/Scripts/python scripts/29_baselines_extended.py --job scvi --seed S | --job cvi --pair {cisplatin,doxorubicin} --seed S (--threads=3)'
    if top == 'c4_multigroupvi':
        return 'scripts/34_c4_multigroupvi.py', '.venv-mgvi/Scripts/python scripts/34_c4_multigroupvi.py --fit --seed {0,1,2} --threads=3'
    if top == 'r4':
        sub = parts[2] if len(parts) > 2 else ''
        if sub.startswith('inputs_'):
            return 'scripts/41_r4_prepare.py', '.venv/Scripts/python scripts/41_r4_prepare.py'
        if sub == 'multigroupvi':
            return 'scripts/42_r4_multigroupvi.py', '.venv-mgvi/Scripts/python scripts/42_r4_multigroupvi.py --config CFG --seed S --threads=3'
        if sub == 'contrastivevi':
            return 'scripts/43_r4_contrastivevi.py', '.venv/Scripts/python scripts/43_r4_contrastivevi.py --config CFG --study STUDY --seed S --threads=3'
        if sub == 'mccvi':
            return 'scripts/44_r4_mccvi_and_collect.py', '.venv/Scripts/python scripts/44_r4_mccvi_and_collect.py --mccvi'
    if top in ('C1_marker_list.csv', 'C1_marker_sources.json'):
        return 'scripts/26_c1_marker_list.py', '.venv/Scripts/python scripts/26_c1_marker_list.py'
    return '', ''


def sha(path, block=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(block), b''):
            h.update(chunk)
    return h.hexdigest()


def commit_of(script):
    if not script:
        return ''
    return subprocess.run(['git', '-C', str(R), 'log', '-1', '--format=%H', '--', script], capture_output=True, text=True).stdout.strip()


def files_of(src: Path, arc: str):
    if src.is_file():
        yield src, arc; return
    for p in sorted(src.rglob('*')):
        if not p.is_file():
            continue
        rel = p.relative_to(src)
        if EXCLUDE_PARTS & set(rel.parts) or (arc, rel.as_posix()) in EXCLUDE_NAMES:
            continue
        yield p, f'{arc}/{rel.as_posix()}'


def main():
    A.mkdir(parents=True, exist_ok=True)
    manifest = []; commits = {}
    for zname, sources in GROUPS.items():
        zpath = A / zname
        entries = [x for s, a in sources for x in files_of(s, a)]
        with zipfile.ZipFile(zpath, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as zf:
            for p, arc in entries:
                zf.write(p, arc)
                script, cmd = producer(f'{zname}/{arc}')
                if script not in commits:
                    commits[script] = commit_of(script)
                manifest.append(dict(archive_file=zname, path_in_archive=arc, bytes=p.stat().st_size, sha256=sha(p), producing_script=script,
                                     command=cmd, producing_script_commit=commits[script], source_path=p.relative_to(R).as_posix()))
        print(zname, len(entries), round(zpath.stat().st_size / 1e6, 1), 'MB', flush=True)
    with open(A / 'MANIFEST.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(manifest[0])); w.writeheader(); w.writerows(manifest)
    zips = [dict(archive_file=z, bytes=(A / z).stat().st_size, sha256=sha(A / z)) for z in GROUPS]
    with open(A / 'ARCHIVE_ZIPS.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(zips[0])); w.writeheader(); w.writerows(zips)
    # excluded inputs: hashes only
    ex = []
    for rel in ('data/raw/GSE216146_chemo_brain.h5ad.gz', 'data/raw/GSE271055_RAW.tar', 'data/processed/recovered_counts.h5ad',
                'data/evidence/m2.cp.reactome.v2025.1.Mm.symbols.gmt', 'data/evidence/m5.go.bp.v2025.1.Mm.symbols.gmt', 'data/evidence/m8.all.v2025.1.Mm.symbols.gmt'):
        p = R / rel
        ex.append(dict(path=rel, bytes=p.stat().st_size if p.exists() else '', sha256=sha(p) if p.exists() else 'missing',
                       reason='public GEO / MSigDB download or regenerable with scripts/05_recover_cohort.py; not archived'))
    with open(A / 'EXCLUDED_INPUT_HASHES.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(ex[0])); w.writeheader(); w.writerows(ex)
    for env, py in (('venv', R / '.venv' / 'Scripts' / 'python.exe'), ('venv-mgvi', R / '.venv-mgvi' / 'Scripts' / 'python.exe')):
        out = subprocess.run([str(py), '-m', 'pip', 'freeze'], capture_output=True, text=True).stdout
        (A / f'pip_freeze_{env}.txt').write_text(out, encoding='utf-8')
    print('manifest rows', len(manifest))


if __name__ == '__main__':
    sys.exit(main())
