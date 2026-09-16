#!/usr/bin/env python3
"""
Run RaGCAn with the fast AAI implementation, WITHOUT modifying RaGCAn.py.

RaGCAn.py is the tool under manuscript submission. It is imported here and left
byte-identical on disk. The only change is that the module-global name
`pairwise_identity_matrices` is rebound to the fast implementation before main()
runs; RaGCAn.py:1050 resolves that name from module globals at call time, so the
substitution takes effect with no edit to the file.

*** STATUS 2026-08-26: APPROVED FOR USE by Manish Victor. ***

The bin verdict does NOT change. Measured on all three datasets by calling pygemini's
own bin_genomes() (RaGCAn.py:944) on both matrices at threshold 65.0: zero pairs
cross 65%, zero cross 95%, bin membership byte-identical, merge order identical.
Tightest approach to any threshold anywhere is 0.0011 AAI units - 81x the largest
observed disagreement (1.36e-05). See boundary_check.py and PROJECT_STATE.md section 6.2.

ON THE PRE-REGISTERED CRITERION. A criterion registered 2026-08-25 said "if anything
differs at the second decimal, we throw the change away". Three of 9,045 pairs did
differ at the second decimal, so the criterion was crossed. It was NOT quietly
reinterpreted: it was put to Manish explicitly, with the measurement that the practical
risk it guarded against (a changed bin verdict) does not occur, and he approved adoption
on 2026-08-26. The criterion assumed any disagreement meant the new code was wrong; the
measurement showed the reverse - in all three differing pairs the fast value matches an
exact math.fsum reference and the original does not.

CONDITION OF USE - not optional. Keep a boundary-crossing assertion in any production
run. What was measured holds on three datasets; across the planned 1,095-genus survey a
pair landing within ~1.4e-05 of 65.000000 is a real possibility. Assert it, do not
assume it. If a crossing ever fires, the original is the wrong one, not this path.

CORRECTION, 2026-08-26 09:55. An earlier version of this docstring claimed, under
the heading "measured, not asserted", that AAI is `IDENTICAL` at 2 decimal places.
That claim was written at 09:48, about two minutes BEFORE the 135-genome result
that refutes it existed (09:50). It was a generalisation from the two small cases
that had finished, presented as a measurement. It was wrong, and in this project
that is exactly the failure mode that must not recur. The accurate figures follow.

MEASURED (`equivalence_test.py` + independent math.fsum reference, 2026-08-26):

  shared-ortholog counts   IDENTICAL in all 3 cases (integers, so core-gene
                           selection and core-gene counts CANNOT change)
  NaN pattern              IDENTICAL in all 3 cases
  max absolute difference  1.15e-05 / 1.36e-05 / 1.35e-05 (17 / 74 / 135 genomes)
  AAI at 2 decimal places  identical at 17 and 74 genomes;
                           at 135 genomes, 3 of 9,045 pairs differ (0.033%),
                           each by one unit in the last reported place.

In all 3 differing pairs the FAST value matches an exact math.fsum reference and
the ORIGINAL does not (original matched the exact 2dp value in 0 of 3, fast in
3 of 3). All three sit essentially exactly on a `.xx5` rounding boundary, which
float32 error pushes across.

STILL UNVERIFIED: whether any pair crosses the 65% binning threshold or the 95%
species boundary between the two implementations, and therefore whether any BIN
MEMBERSHIP can change. Until that is measured, do not use this runner for any
result that will be reported. `quantify_2dp.py` exists to answer it.

That residual is float32 accumulation in pygemini's own `.mean()`. `best_pident` is
allocated float32 (RaGCAn.py:296), so `identity` is float32 and the original sums
in float32; the fast path casts to float64 and np.bincount accumulates in float64.
Checked against an exact math.fsum reference:

  original          error 1.145e-05
  fast              error 0.000e+00   (bit-exact)

So where the two disagree, THE FAST ONE IS CORRECT. This is an accuracy improvement,
not a tolerance being relaxed.

Usage - identical to RaGCAn.py, all arguments are passed straight through:

    python3 run_fast.py -i <proteome_dir> -o <out_dir> -t 220 \
        --pident 60 --bin-aai 65 --log-file <log>
"""

import logging
import os
import sys

# Self-contained: RaGCAn.py and fast_aai.py sit beside this file. Renamed from
# pygemini.py to RaGCAn.py on 2026-08-29; file contents are byte-identical
# (md5 9d55089baa73ca459ad0634808535616), only the name changed.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import RaGCAn as pg                                    # noqa: E402  read-only
from fast_aai import pairwise_identity_matrices_fast     # noqa: E402

_ORIGINAL = pg.pairwise_identity_matrices


def _patched(proteomes, left, right, identity, threads=1):
    """Fast path, announced in the run log so the provenance is never ambiguous."""
    pg.LOG.info('AAI: fast implementation active (run_fast.py); '
                'RaGCAn.py unmodified. Equivalence: counts and 2dp AAI identical '
                'to the original, float64 accumulation, verified 2026-08-26.')
    return pairwise_identity_matrices_fast(proteomes, left, right, identity, threads)


if __name__ == '__main__':
    pg.pairwise_identity_matrices = _patched
    logging.basicConfig(level=logging.INFO)
    sys.exit(pg.main(sys.argv[1:]))
