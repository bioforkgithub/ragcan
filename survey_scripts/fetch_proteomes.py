#!/usr/bin/env python3
"""
Assemble proteomes/<Genus>/ for one genus: reuse what is already on disk, download
only what is genuinely missing.

CRASH SAFETY (this whole pipeline assumes it will be killed mid-run):
  * A proteome is only ever "present" if a non-empty .faa exists at its final path.
  * Downloads go to <name>.faa.part, are gzip-integrity-checked in full while being
    decompressed, and are then os.replace()d into place. os.replace is atomic on the
    same filesystem, so a killed download leaves a .part file that is ignored and
    overwritten next time. A truncated .faa can never appear.
  * Reused files are symlinks; a broken symlink counts as absent.
  * The per-genus manifest reports/fetch_<Genus>.tsv is written atomically at the end.

Species whose assembly carries no *_protein.faa.gz (no protein annotation) are NOT
silently dropped: they are recorded with status=no_protein_annotation in the manifest
and in reports/skipped_species.tsv, and excluded from the run for that genus.
"""

import argparse
import csv
import gzip

import os
import random
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROJ = os.path.dirname(ROOT)

SELECTION = os.path.join(ROOT, 'reports', 'genome_selection.tsv')
PROTEOMES = os.path.join(ROOT, 'proteomes')
REPORTS = os.path.join(ROOT, 'reports')

# Directories searched for proteomes already on disk. Matched by ACCESSION parsed
# out of the filename, not by name, so a differently-named copy is still reused.
REUSE_DIRS = [
    os.path.join(PROJ, 'Stutzerimonas_Validation', 'proteomes_pseudomonas_full'),
    os.path.join(PROJ, 'Stutzerimonas_Validation', 'proteomes_stutzerimonas'),
    os.path.join(PROJ, 'Stutzerimonas_Validation', 'proteomes_halopseudomonas'),
    os.path.join(PROJ, 'Stutzerimonas_Validation', 'proteomes_pseudomonadaceae_all'),
    os.path.join(PROJ, 'pygemini_work', 'PY-GEMINI', 'chryseo_all', 'proteomes'),
    os.path.join(PROJ, 'pygemini_work', 'PY-GEMINI', 'chryseo_all', 'proteomes_hq'),
    os.path.join(PROJ, 'QualityTiers', 'tiers', 'A_all138'),
    os.path.join(PROJ, 'QualityTiers', 'tiers', 'B_qc116'),
    os.path.join(PROJ, 'QualityTiers', 'tiers', 'C_complete30'),
    os.path.join(PROJ, 'QualityTiers', 'tiers', 'PsB_qc'),
    os.path.join(PROJ, 'QualityTiers', 'tiers', 'PsC_complete'),
]

ACC_RE = __import__('re').compile(r'(GC[AF]_\d+\.\d+)')

MAX_WORKERS = 4          # polite to NCBI
RETRIES = 5
TIMEOUT = 300


def build_reuse_index():
    """accession -> absolute path of an existing, non-empty .faa on disk."""
    index = {}
    for d in REUSE_DIRS:
        if not os.path.isdir(d):
            continue
        for fn in os.listdir(d):
            if not fn.endswith('.faa'):
                continue
            m = ACC_RE.search(fn)
            if not m:
                continue
            p = os.path.join(d, fn)
            try:
                if os.path.getsize(p) == 0:
                    continue
            except OSError:
                continue
            index.setdefault(m.group(1), p)
    return index


def load_selection(genus):
    with open(SELECTION, newline='') as fh:
        return [r for r in csv.DictReader(fh, delimiter='\t') if r['genus'] == genus]


def present(path):
    """True only for a real, non-empty file (or a symlink resolving to one)."""
    try:
        return os.path.getsize(path) > 0
    except OSError:
        return False


def download(row, dest, session):
    """Fetch <ftp_path>/<asm>_protein.faa.gz, decompress, atomically install.

    Returns 'downloaded', or 'no_protein_annotation' if the file does not exist.
    Raises on any other persistent failure.
    """
    ftp = row['ftp_path'].replace('ftp://', 'https://').rstrip('/')
    asm = ftp.rsplit('/', 1)[-1]
    url = '%s/%s_protein.faa.gz' % (ftp, asm)
    part = dest + '.part'
    gzpart = dest + '.gz.part'
    last = None
    for attempt in range(RETRIES):
        try:
            with session.get(url, timeout=TIMEOUT, stream=True) as r:
                if r.status_code == 404:
                    return 'no_protein_annotation'
                r.raise_for_status()
                with open(gzpart, 'wb') as out:
                    for chunk in r.iter_content(1 << 20):
                        out.write(chunk)
            # Decompress the whole file. A truncated transfer raises here, so an
            # incomplete download can never be installed as a .faa.
            with gzip.open(gzpart, 'rb') as gz, open(part, 'wb') as out:
                shutil.copyfileobj(gz, out, 1 << 20)
            if os.path.getsize(part) == 0:
                raise IOError('empty proteome after decompression')
            os.replace(part, dest)          # atomic
            os.unlink(gzpart)
            return 'downloaded'
        except Exception as exc:            # noqa: BLE001 - retried below
            last = exc
            for p in (part, gzpart):
                try:
                    os.unlink(p)
                except OSError:
                    pass
            if attempt < RETRIES - 1:
                time.sleep((2 ** attempt) + random.random())
    raise RuntimeError('%s: %s' % (url, last))


def fetch_genus(genus, log=print):
    rows = load_selection(genus)
    if not rows:
        raise SystemExit('no rows for genus %s in %s' % (genus, SELECTION))
    outdir = os.path.join(PROTEOMES, genus)
    os.makedirs(outdir, exist_ok=True)

    reuse = build_reuse_index()
    manifest = []
    todo = []

    for row in rows:
        dest = os.path.join(outdir, row['faa_name'])
        if present(dest):
            manifest.append((row, 'already_present', os.path.realpath(dest)))
            continue
        src = reuse.get(row['accession'])
        if src and present(src):
            # remove a stale/broken symlink before relinking
            if os.path.islink(dest) or os.path.exists(dest):
                try:
                    os.unlink(dest)
                except OSError:
                    pass
            os.symlink(src, dest)
            manifest.append((row, 'reused', src))
            continue
        todo.append((row, dest))

    log('[fetch] %s: %d selected, %d already present, %d reused from disk, '
        '%d to download' % (genus, len(rows),
                            sum(1 for m in manifest if m[1] == 'already_present'),
                            sum(1 for m in manifest if m[1] == 'reused'),
                            len(todo)))

    skipped = []
    if todo:
        session = requests.Session()
        session.headers['User-Agent'] = 'RaGCAn-LPSN-survey/1.0 (contact: 114657046+bioforkgithub@users.noreply.github.com)'

        def one(job):
            row, dest = job
            try:
                st = download(row, dest, session)
            except Exception as exc:                     # noqa: BLE001
                return (row, 'download_failed', str(exc))
            if st == 'no_protein_annotation':
                return (row, 'no_protein_annotation', '')
            return (row, 'downloaded', dest)

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
            for res in pool.map(one, todo):
                manifest.append(res)
                if res[1] != 'downloaded':
                    skipped.append(res)
                    log('[fetch] %s: %s %s -> %s %s'
                        % (genus, res[0]['species'], res[0]['accession'], res[1], res[2]))

    # per-genus manifest, written atomically
    mpath = os.path.join(REPORTS, 'fetch_%s.tsv' % genus)
    tmp = mpath + '.tmp'
    with open(tmp, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['genus', 'species', 'accession', 'assembly_level',
                    'protein_coding_gene_count', 'faa_name', 'status', 'detail'])
        for row, status, detail in sorted(manifest, key=lambda m: m[0]['species']):
            w.writerow([row['genus'], row['species'], row['accession'],
                        row['assembly_level'], row['protein_coding_gene_count'],
                        row['faa_name'], status, detail])
    os.replace(tmp, mpath)

    # append skipped species to the global record (append-only; never truncated)
    if skipped:
        spath = os.path.join(REPORTS, 'skipped_species.tsv')
        new = not os.path.exists(spath)
        with open(spath, 'a', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            if new:
                w.writerow(['genus', 'species', 'accession', 'reason', 'detail'])
            for row, status, detail in skipped:
                w.writerow([row['genus'], row['species'], row['accession'], status, detail])

    usable = [m for m in manifest if m[1] in ('reused', 'already_present', 'downloaded')]
    hard_fail = [m for m in manifest if m[1] == 'download_failed']
    if hard_fail:
        raise RuntimeError('%s: %d proteomes failed to download after %d retries: %s'
                           % (genus, len(hard_fail), RETRIES,
                              ', '.join(m[0]['accession'] for m in hard_fail[:5])))

    counts = {
        'selected': len(rows),
        'usable': len(usable),
        'reused': sum(1 for m in manifest if m[1] in ('reused', 'already_present')),
        'downloaded': sum(1 for m in manifest if m[1] == 'downloaded'),
        'no_annotation': sum(1 for m in manifest if m[1] == 'no_protein_annotation'),
    }
    log('[fetch] %s: done - %s' % (genus, counts))
    return outdir, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('genus')
    args = ap.parse_args()
    fetch_genus(args.genus)


if __name__ == '__main__':
    sys.exit(main())
