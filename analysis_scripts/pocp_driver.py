#!/usr/bin/env python3
"""
POCP as published (Qin, Xie et al. 2014, J Bacteriol), computed across all 1,160 genera.

WHY THIS EXISTS
  The manuscript reports RaGCAn's screening performance against GTDB. The obvious reviewer
  question is how that compares with POCP, the established genome-based genus boundary. An
  earlier internal comparison used RaGCAn's own reciprocal best hits at 60% identity, which is
  STRICTER than POCP's criteria and therefore biased the comparison in RaGCAn's favour. That
  result was held back for exactly this reason. This run applies POCP's OWN criteria so the
  comparison is fair whichever way it falls.

PUBLISHED CRITERIA, applied here
  A protein in genome A counts as conserved with respect to genome B if BLASTP finds a hit with
    - sequence identity      > 40%
    - alignment covering     > 50% of the QUERY protein
    - e-value                < 1e-5
  POCP = (C1 + C2) / (T1 + T2) * 100, where C is the count of DISTINCT conserved proteins and T
  the total protein count. Note this is ONE-WAY and ASYMMETRIC: C1 and C2 are counted separately.
  Genus boundary: every pairwise POCP > 50%.

DEVIATION FROM THE PAPER, AND WHY IT IS ACCEPTABLE
  Qin et al. used BLASTP; we use DIAMOND in --very-sensitive mode with the same thresholds.
  This follows established practice - the reference POCP implementations (hoelzer/pocp, POCP-nf)
  substitute DIAMOND for BLASTP for the same reason. It must still be stated in any write-up.

  Everything else follows the publication. No RaGCAn code is used for the POCP arm; the only
  shared input is the proteome files.

Resumable: a genus already marked done is never recomputed. Hit tables are deleted per genus
once its POCP matrix is written, so peak disk stays bounded.
"""
import argparse, csv, fcntl, os, re, shutil, subprocess, sys, time
from collections import defaultdict
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.join(os.path.dirname(ROOT), 'LPSN_Genus_Survey')
PROTEOMES = os.path.join(SURVEY, 'proteomes')
SELECTION = os.path.join(SURVEY, 'reports', 'genome_selection.tsv')
RESULTS = os.path.join(ROOT, 'results')
REPORTS = os.path.join(ROOT, 'reports')
LOGS = os.path.join(ROOT, 'logs')
STATUS = os.path.join(REPORTS, 'pocp_status.tsv')
LOCK = os.path.join(REPORTS, 'pocp.lock')
PIDFILE = os.path.join(REPORTS, 'pocp.pid')
# DIAMOND is found at start-up: --diamond, then $RAGCAN_DIAMOND, then PATH. See main().
DIAMOND = None

MIN_ID, MIN_QCOV, MAX_E, POCP_THRESH = 40.0, 50.0, 1e-5, 50.0

def log(msg):
    print('%s  %s' % (time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), msg), flush=True)

def read_status():
    st = {}
    if os.path.exists(STATUS):
        for r in csv.DictReader(open(STATUS), delimiter='\t'):
            st[r['genus']] = r
    return st

def write_status(st):
    tmp = STATUS + '.tmp'
    cols = ['genus', 'n_genomes', 'state', 'bins', 'min_pocp', 'mean_pocp', 'wall_seconds', 'note']
    with open(tmp, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for g in sorted(st, key=lambda x: (int(st[x].get('n_genomes') or 0), x)):
            w.writerow({c: st[g].get(c, '') for c in cols})
        fh.flush(); os.fsync(fh.fileno())
    os.replace(tmp, STATUS)

def genus_sizes():
    per = defaultdict(int)
    for r in csv.DictReader(open(SELECTION), delimiter='\t'):
        per[r['genus']] += 1
    return per

def bins_complete_linkage(M, thresh):
    n = M.shape[0]
    bins = [[i] for i in range(n)]
    link = M.astype(float).copy()
    link[np.isnan(link)] = -np.inf
    np.fill_diagonal(link, -np.inf)
    alive = np.ones(n, dtype=bool)
    while alive.sum() > 1:
        masked = np.where(alive[:, None] & alive[None, :], link, -np.inf)
        f = int(np.argmax(masked)); i, j = divmod(f, n)
        if masked[i, j] < thresh: break
        if i > j: i, j = j, i
        bins[i] = bins[i] + bins[j]; bins[j] = []
        alive[j] = False
        link[i, :] = np.minimum(link[i, :], link[j, :]); link[:, i] = link[i, :]
        link[i, i] = -np.inf
    return [b for b in bins if b]

def do_genus(genus, threads):
    t0 = time.time()
    src = os.path.join(PROTEOMES, genus)
    faas = sorted(f for f in os.listdir(src) if f.endswith('.faa'))
    if len(faas) < 2:
        return dict(state='failed', note='fewer than 2 proteomes')
    out = os.path.join(RESULTS, genus); os.makedirs(out, exist_ok=True)
    work = os.path.join(out, 'working'); os.makedirs(work, exist_ok=True)

    # pool, prefixing every sequence id with its genome index so membership is recoverable
    pooled = os.path.join(work, 'pooled.faa')
    totals = []
    with open(pooled, 'w') as o:
        for gi, f in enumerate(faas):
            n = 0
            for line in open(os.path.join(src, f)):
                if line.startswith('>'):
                    n += 1
                    o.write('>%d|%d\n' % (gi, n))
                else:
                    o.write(line)
            totals.append(n)
    T = np.array(totals, dtype=float)

    db = os.path.join(work, 'pooled')
    subprocess.run([DIAMOND, 'makedb', '--in', pooled, '-d', db,
                    '--threads', str(threads), '--quiet'], check=True)
    hits = os.path.join(work, 'hits.tsv')
    subprocess.run([DIAMOND, 'blastp', '-q', pooled, '-d', db, '-o', hits,
                    '--outfmt', '6', 'qseqid', 'sseqid',
                    '--threads', str(threads), '--very-sensitive',
                    '--id', str(MIN_ID), '--query-cover', str(MIN_QCOV),
                    '--evalue', str(MAX_E), '--max-target-seqs', '500',
                    '--block-size', '0.4', '--index-chunks', '4', '--quiet'], check=True)

    # C[a][b] = number of DISTINCT proteins of genome a with a qualifying hit in genome b
    n = len(faas)
    seen = defaultdict(set)
    with open(hits) as fh:
        for line in fh:
            q, s = line.rstrip('\n').split('\t', 1)
            qa, qp = q.split('|', 1)
            sa = s.split('|', 1)[0]
            if qa != sa:
                seen[(int(qa), int(sa))].add(qp)
    C = np.zeros((n, n))
    for (a, b), st in seen.items():
        C[a, b] = len(st)

    P = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i == j: P[i, j] = 100.0
            elif T[i] + T[j] > 0: P[i, j] = (C[i, j] + C[j, i]) / (T[i] + T[j]) * 100.0
    ids = [f[:-4] for f in faas]

    os.makedirs(os.path.join(out, 'reports'), exist_ok=True)
    with open(os.path.join(out, 'reports', 'pocp_matrix.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['genome'] + ids)
        for i, g in enumerate(ids):
            w.writerow([g] + ['%.2f' % P[i, k] if P[i, k] == P[i, k] else 'NA' for k in range(n)])

    bins = bins_complete_linkage(P, POCP_THRESH)
    with open(os.path.join(out, 'reports', 'pocp_bins.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['bin', 'genomes', 'members'])
        for bi, b in enumerate(sorted(bins, key=lambda x: (-len(x), x[0])), 1):
            w.writerow([bi, len(b), ' '.join(ids[i] for i in b)])

    off = P[~np.eye(n, dtype=bool)]

    # RETAIN the hit table, gzipped, as the raw evidence behind this genus's POCP values.
    # 2026-08-30, at Manish's instruction: "we can provide the raw files. You can keep it
    # compressed if memory is the criteria." The searches are the expensive artefact; the
    # POCP arithmetic downstream is trivial. Deleting them once already cost a 30-hour re-run.
    #
    # NOTE ON WHAT THESE FILES ARE: DIAMOND applied POCP's criteria AT SEARCH TIME
    # (--id 40 --query-cover 50 --evalue 1e-5), so these tables are FILTERED, not raw hits.
    # They are sufficient to reproduce every POCP value reported here, and they are NOT
    # sufficient to recompute a metric at a looser threshold. That would need a new search.
    import gzip, shutil as _sh
    hitdir = os.path.join(out, 'hits'); os.makedirs(hitdir, exist_ok=True)
    gz = os.path.join(hitdir, 'pocp_hits.tsv.gz')
    try:
        with open(hits, 'rb') as fi, gzip.open(gz, 'wb', compresslevel=6) as fo:
            _sh.copyfileobj(fi, fo, 1 << 22)
        # The protein totals MUST travel with the hit table. The driver counts sequences in the
        # .faa files, and those differ from NCBI's annotation counts by a few proteins per genome
        # (verified 2026-08-30: 4601 counted vs 4611 annotated for one Alkalihalobacterium
        # genome). Without these numbers a third party cannot reproduce POCP from the shipped
        # files, because POCP's denominator is (T_i + T_j).
        with open(os.path.join(hitdir, 'genome_totals.tsv'), 'w') as fh:
            fh.write('genome_index\tgenome\tproteins_counted\n')
            for gi, f in enumerate(faas):
                fh.write('%d\t%s\t%d\n' % (gi, f[:-4], totals[gi]))
        with open(os.path.join(hitdir, 'README.txt'), 'w') as fh:
            fh.write('pocp_hits.tsv.gz - DIAMOND blastp output for %s\n' % genus)
            fh.write('columns: qseqid sseqid   (ids are <genome_index>|<protein_number>)\n')
            fh.write('genome_index maps to the sorted .faa list in pocp_matrix.tsv column order\n')
            fh.write('FILTERED AT SEARCH TIME to POCP criteria: --id 40 --query-cover 50 --evalue 1e-5\n')
            fh.write('--very-sensitive. Sufficient to reproduce POCP; NOT usable at looser thresholds.\n')
            fh.write('\ngenome_totals.tsv gives the protein count per genome, counted from the\n')
            fh.write('.faa files. POCP = (C_i + C_j)/(T_i + T_j)*100. These two files reproduce\n')
            fh.write('pocp_matrix.tsv exactly, with no other input.\n')
    except Exception as e:                                   # noqa: BLE001
        log('%s: WARNING could not retain hit table: %s' % (genus, e))

    for f in ('pooled.faa', 'pooled.dmnd', 'hits.tsv'):
        try: os.unlink(os.path.join(work, f))
        except OSError: pass
    try: os.rmdir(work)
    except OSError: pass
    return dict(state='done', bins=len(bins),
                min_pocp='%.2f' % np.nanmin(off) if off.size else '',
                mean_pocp='%.2f' % np.nanmean(off) if off.size else '',
                wall_seconds=int(time.time() - t0), note='')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-t', '--threads', type=int, default=96)
    ap.add_argument('--retry-failed', action='store_true')
    ap.add_argument('--diamond', metavar='PATH',
                    help='DIAMOND executable (default: $RAGCAN_DIAMOND, then PATH)')
    a = ap.parse_args()

    # Resolve DIAMOND before anything is written. If this were left to each genus, a missing
    # binary would mark every genus 'failed' one by one instead of stopping once.
    global DIAMOND
    DIAMOND = a.diamond or os.environ.get('RAGCAN_DIAMOND') or shutil.which('diamond')
    if not DIAMOND or not os.access(DIAMOND, os.X_OK):
        sys.exit('DIAMOND not found. Put it on PATH, set RAGCAN_DIAMOND, or pass --diamond.')
    for d in (RESULTS, REPORTS, LOGS): os.makedirs(d, exist_ok=True)
    fh = open(LOCK, 'w')
    try: fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        log('another POCP driver holds the lock - exiting'); return 0
    open(PIDFILE, 'w').write('%d\n' % os.getpid())

    sizes = genus_sizes()
    st = read_status()
    for g, n in sizes.items():
        st.setdefault(g, dict(genus=g, n_genomes=str(n), state='pending'))
        st[g]['n_genomes'] = str(n)
    write_status(st)

    # 'running' means a previous driver was killed mid-genus. Treat it as resumable, exactly as
    # the survey driver does - otherwise a kill silently drops that genus from the run forever.
    stale = [g for g in sizes if st[g]['state'] == 'running']
    for g in stale:
        log('RESUME: %s was left "running" by a killed driver - re-entering it' % g)
    todo = [g for g, _ in sorted(sizes.items(), key=lambda kv: kv[1])
            if st[g]['state'] in ('pending', 'running')
            or (a.retry_failed and st[g]['state'] == 'failed')]
    log('POCP driver start, pid %d, threads %d, %d genera to do' % (os.getpid(), a.threads, len(todo)))
    for g in todo:
        st[g]['state'] = 'running'; write_status(st)
        try:
            r = do_genus(g, a.threads)
        except Exception as e:                       # noqa: BLE001
            r = dict(state='failed', note='%s: %s' % (type(e).__name__, e))
        st[g].update(r); write_status(st)
        log('%s: %s%s' % (g, r['state'],
                          '  bins=%s min=%s wall=%ss' % (r.get('bins'), r.get('min_pocp'), r.get('wall_seconds'))
                          if r['state'] == 'done' else '  ' + r.get('note', '')))
    log('POCP driver finished')
    return 0

if __name__ == '__main__':
    sys.exit(main())
