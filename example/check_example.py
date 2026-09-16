#!/usr/bin/env python3
"""Compare your run of the worked example with the expected result.

    python3 example/check_example.py [RESULT_DIR]     (default: example/result)

It passes if the same genomes are grouped together. It does not demand identical files:
a newer DIAMOND can move an AAI value by a fraction of a point without changing any group.
Tested on DIAMOND 2.1.8 and 2.2.6, where the largest AAI difference was 0.13 points.
"""
import csv
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
got_dir = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / 'result'


def groups(path):
    with open(path) as fh:
        return {frozenset(r['members'].split()) for r in csv.DictReader(fh, delimiter='\t')}


def aai(path):
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    names = rows[0][1:]
    return {(r[0], names[j]): float(v) for r in rows[1:]
            for j, v in enumerate(r[1:]) if v not in ('', 'NA', 'nan')}


got_bins = got_dir / 'reports' / 'proposed_bins.tsv'
if not got_bins.is_file():
    sys.exit(f'No result found at {got_bins}. Run the example first (see example/README.md).')

want, got = groups(HERE / 'expected' / 'proposed_bins.tsv'), groups(got_bins)
a, b = aai(HERE / 'expected' / 'aai_matrix.tsv'), aai(got_dir / 'reports' / 'aai_matrix.tsv')
diff = max(abs(a[k] - b[k]) for k in a if k in b)

print(f'groups expected {len(want)}, found {len(got)}')
for g in sorted(got, key=len, reverse=True):
    print(f'  {len(g)} genome(s): ' + ', '.join(sorted(x.split("__")[0] for x in g)))
print(f'largest AAI difference from the expected matrix: {diff:.2f} points')

if got != want:
    print('FAIL: the genomes are grouped differently from the expected result.')
    sys.exit(1)
if diff > 1.0:
    print('NOTE: groups match, but some AAI values differ by more than 1 point.')
print('PASS: same groups as the expected result.')
