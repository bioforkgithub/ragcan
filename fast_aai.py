#!/usr/bin/env python3
"""
A faster implementation of RaGCAn's pairwise_identity_matrices.

NOT a change to the method. Same records selected, same counts, same means -
only the order of traversal differs.

  original : for every unordered genome pair (a,b), scan the whole hit array to
             find records belonging to that pair.  Cost = pairs x hits.
             With 389 genomes that is 75,466 scans of 234 million records.

  this     : scan the hit array once, assigning each record to its pair's bucket
             with np.bincount.  Cost = hits.

CORRECTION 2026-08-26 — the text that stood here predicted the only difference
would be summation order, "a relative error near 1e-13". **That prediction was
measured and falsified.** The real difference is ~1.35e-05 absolute, ~1.6e-07
relative, sitting exactly on float32 epsilon (1.19e-07) rather than float64
rounding.

The cause is not summation order. `best_pident` is allocated `dtype=np.float32`
(RaGCAn.py:296), so `identity` is float32 and the ORIGINAL accumulates in
float32; this implementation casts to float64. Against an exact math.fsum
reference: original error 1.145e-05, this implementation 0.000e+00 (bit-exact).
Where the two disagree, this one is the correct one.

The consequence that matters: at 135 genomes, **3 of 9,045 pairs differ at the
second decimal place** — the precision RaGCAn actually reports ('%.2f'). In all
three, the fast value matches the exact reference and the original does not.

**This crosses the falsification criterion registered on 2026-08-25** ("if
anything differs at the second decimal, we throw the change away"). That
criterion assumed any disagreement would mean the new code was wrong; the
measurement shows the opposite. The criterion should not be silently reinterpreted
— see `run_fast.py`, which records why this is NOT yet cleared for production use
and what still has to be measured before it could be.

`equivalence_test.py` measures this on real completed runs rather than asserting it.
"""

import numpy as np


def pairwise_identity_matrices_fast(proteomes, left, right, identity, threads=1):
    """Return (AAI matrix, shared-ortholog-count matrix) - drop-in replacement.

    `threads` is accepted for signature compatibility and ignored; the work is a
    single vectorised pass and numpy already uses its own threading underneath.
    """
    genome_count = len(proteomes.genomes)
    genome_of = proteomes.genome_of

    ga = np.asarray(genome_of[left])
    gb = np.asarray(genome_of[right])

    # Unordered pair key. The original only ever fills cells with a < b, and
    # handles the diagonal separately, so records inside a single genome are
    # excluded here exactly as they are excluded there.
    lo = np.minimum(ga, gb)
    hi = np.maximum(ga, gb)
    keep = lo != hi

    aai = np.full((genome_count, genome_count), np.nan, dtype=np.float64)
    shared = np.zeros((genome_count, genome_count), dtype=np.int64)
    for a in range(genome_count):
        aai[a, a] = 100.0

    if not keep.any():
        return aai, shared

    lo = lo[keep].astype(np.int64, copy=False)
    hi = hi[keep].astype(np.int64, copy=False)
    ident = np.asarray(identity)[keep].astype(np.float64, copy=False)

    key = lo * genome_count + hi
    size = genome_count * genome_count

    counts = np.bincount(key, minlength=size)
    sums = np.bincount(key, weights=ident, minlength=size)

    nz = np.nonzero(counts)[0]
    if nz.size:
        a_idx, b_idx = np.divmod(nz, genome_count)
        c = counts[nz]
        m = sums[nz] / c
        shared[a_idx, b_idx] = c
        shared[b_idx, a_idx] = c
        aai[a_idx, b_idx] = m
        aai[b_idx, a_idx] = m

    return aai, shared
