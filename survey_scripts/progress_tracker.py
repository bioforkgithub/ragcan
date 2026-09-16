#!/usr/bin/env python3
"""
Progress tracker for the LPSN genus survey.

A plain detached process. It has NO dependency on any AI agent, and it holds no state
of its own: everything it prints is rebuilt from files on disk each cycle, so killing
it and restarting it is a no-op. An agent may read what this writes; an agent must
never be what produces it. That is the whole point - the session that was supposed to
be watching the last long run died overnight, which is why this file exists.

Every cycle (default 120 s) it:
  1. rebuilds reports/PROGRESS.md, atomically, from reports/survey_status.tsv and
     reports/survey_results.tsv - safe to `cat` at any instant;
  2. appends any newly observed genus state change to reports/events.log as an
     independent OBS_ observation (the driver writes its own lines to the same
     append-only file; two independent observers of the same run is deliberate);
  3. flags every completed genus MATCH or MISMATCH against the expectation registered
     in PLAN.md Phase 1. A MISMATCH is the go/no-go signal the pilot exists to
     produce. It is recorded prominently and never suppressed.

The ETA it prints is EXTRAPOLATED, not measured, and is labelled as such wherever it
appears. It is derived from wall times actually observed on genera already finished,
scaled by (total proteins)^2 - the cost of an all-against-all search - using NCBI's
per-assembly protein counts for the genera not yet run.
"""

import argparse
import csv
import datetime
import fcntl
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REPORTS = os.path.join(ROOT, 'reports')

STATUS = os.path.join(REPORTS, 'survey_status.tsv')
RESULTS = os.path.join(REPORTS, 'survey_results.tsv')
SELECTION = os.path.join(REPORTS, 'genome_selection.tsv')
EVENTS = os.path.join(REPORTS, 'events.log')
PROGRESS = os.path.join(REPORTS, 'PROGRESS.md')
FAILURES = os.path.join(REPORTS, 'failures.tsv')
LOCK = os.path.join(REPORTS, 'tracker.lock')
PIDFILE = os.path.join(REPORTS, 'tracker.pid')
DRIVER_PID = os.path.join(REPORTS, 'driver.pid')

EXPECTED = {
    'Chryseobacterium': 'coherent', 'Stutzerimonas': 'coherent',
    'Halopseudomonas': 'coherent', 'Bordetella': 'coherent',
    'Brucella': 'coherent', 'Yersinia': 'coherent',
    'Pseudomonas': 'split', 'Clostridium': 'split', 'Bacillus': 'split',
    'Lactobacillus': 'split', 'Streptococcus': 'split', 'Mycobacterium': 'split',
    'Aeromonas': 'split', 'Ralstonia': 'split',
    'Paraburkholderia': 'split', 'Burkholderia': 'split',
}
EXPECT_NOTE = {
    'Aeromonas': 'recently split, documented',
    'Ralstonia': 'recently split, documented',
    'Paraburkholderia': 'recently split, documented',
    'Burkholderia': 'recently split, documented',
}


def utc():
    return datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')


def event(kind, genus, detail=''):
    with open(EVENTS, 'a') as fh:
        fh.write('%s\t%s\t%s\t%s\n' % (utc(), kind, genus, detail))
        fh.flush()
        os.fsync(fh.fileno())


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def alive(pidfile):
    try:
        pid = int(open(pidfile).read().strip())
    except (OSError, ValueError):
        return None, False
    try:
        os.kill(pid, 0)
        return pid, True
    except OSError:
        return pid, False


def genus_proteins():
    """Total protein-coding genes per genus, from NCBI's own per-assembly counts.
    Used only to extrapolate an ETA; never used to filter anything."""
    tot = {}
    for r in read_tsv(SELECTION):
        try:
            tot[r['genus']] = tot.get(r['genus'], 0) + int(r['protein_coding_gene_count'])
        except (ValueError, KeyError):
            tot.setdefault(r['genus'], 0)
    return tot


def fmt_dur(seconds):
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return '%dh %02dm' % (h, m)
    if m:
        return '%dm %02ds' % (m, s)
    return '%ds' % s


def parse_utc(s):
    try:
        return datetime.datetime.strptime(s, '%Y-%m-%dT%H:%M:%SZ')
    except (ValueError, TypeError):
        return None


def build(status_rows, result_rows, proteins):
    res = {r['genus']: r for r in result_rows}
    by_state = {}
    for r in status_rows:
        by_state.setdefault(r['state'], []).append(r)

    n = len(status_rows)
    done = by_state.get('done', [])
    running = by_state.get('running', []) + by_state.get('fetching', [])
    failed = by_state.get('failed', [])
    pending = by_state.get('pending', [])

    L = []
    A = L.append
    A('# LPSN genus-coherence survey - Phase 1 progress')
    A('')
    A('Regenerated automatically by `scripts/progress_tracker.py` at **%s**.' % utc())
    A('This file is rebuilt from `survey_status.tsv` and `survey_results.tsv` every')
    A('cycle; it holds no state and is safe to read at any moment.')
    A('')

    dpid, dalive = alive(DRIVER_PID)
    tpid, _ = alive(PIDFILE)
    A('| process | pid | alive |')
    A('|---|---|---|')
    A('| driver (`run_survey.py`) | %s | %s |' % (dpid, 'YES' if dalive else '**NO**'))
    A('| tracker (`progress_tracker.py`) | %s | YES |' % tpid)
    A('')
    if not dalive and pending:
        A('> **The driver is not running and work remains.** Restart it with the one')
        A('> command in `HOW_TO_RESUME.md`. Nothing is lost; it resumes where it stopped.')
        A('')

    A('## Where it is')
    A('')
    A('| state | genera |')
    A('|---|---:|')
    A('| done | %d |' % len(done))
    A('| running / fetching | %d |' % len(running))
    A('| pending | %d |' % len(pending))
    A('| failed | %d |' % len(failed))
    A('| **total** | **%d** |' % n)
    A('')

    if running:
        A('### Running now')
        A('')
        A('| genus | species | state | since (UTC) | elapsed |')
        A('|---|---:|---|---|---:|')
        for r in running:
            t0 = parse_utc(r.get('started_utc', ''))
            el = fmt_dur((datetime.datetime.utcnow() - t0).total_seconds()) if t0 else '?'
            A('| %s | %s | %s | %s | %s |'
              % (r['genus'], r['n_species'], r['state'], r.get('started_utc', '?'), el))
        A('')

    # ---- results so far -------------------------------------------------
    A('## Results so far')
    A('')
    if not res:
        A('_No genus has completed yet._')
        A('')
    else:
        A('| genus | genomes | core genes | bins @65 | min AAI | mean AAI | pairs >=95 | '
          'first split | expected | verdict | wall |')
        A('|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|')
        order = [r['genus'] for r in status_rows if r['genus'] in res]
        for g in order:
            r = res[g]
            mark = r.get('agrees_with_expectation', '')
            mark = '**MISMATCH**' if mark == 'MISMATCH' else mark
            A('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |'
              % (g, r.get('genomes_used', ''), r.get('core_genes', ''),
                 r.get('bins_at_65', ''), r.get('min_aai', ''), r.get('mean_aai', ''),
                 r.get('pairs_ge_95', ''), r.get('first_split_threshold', '') or 'none <=85',
                 r.get('expectation', ''), mark,
                 fmt_dur(r['wall_seconds']) if r.get('wall_seconds') else ''))
        A('')

        A('### Agreement with the registered expectation (the go/no-go gate)')
        A('')
        mism = [r for r in res.values() if r.get('agrees_with_expectation') == 'MISMATCH']
        match = [r for r in res.values() if r.get('agrees_with_expectation') == 'MATCH']
        A('- MATCH: **%d** of %d completed' % (len(match), len(res)))
        A('- MISMATCH: **%d**' % len(mism))
        A('')
        if mism:
            A('> A mismatch is not a bug. It is the signal this pilot exists to produce.')
            A('')
            A('| genus | expected | observed | bins @65 | first split | where the expectation comes from |')
            A('|---|---|---|---:|---:|---|')
            for r in sorted(mism, key=lambda x: x['genus']):
                A('| **%s** | %s | %s | %s | %s | %s |'
                  % (r['genus'], r.get('expectation', ''), r.get('verdict', ''),
                     r.get('bins_at_65', ''),
                     r.get('first_split_threshold', '') or 'none <=85',
                     r.get('expectation_basis', '') or EXPECT_NOTE.get(r['genus'], '')))
            A('')
            if any(r['genus'] in EXPECT_NOTE for r in mism):
                A('> Note on the four "recently split, documented" genera '
                  '(*Aeromonas*, *Ralstonia*, *Paraburkholderia*, *Burkholderia*): '
                  'PLAN.md names them but does **not** state a bin count for them. '
                  'Scoring them as "expect >1 bin" is an interpretation made by this '
                  'pipeline. A mismatch there may mean the expectation was '
                  'under-specified rather than that the screen disagreed with a '
                  'documented answer - check the genus before reporting it as a '
                  'disagreement.')
                A('')

        A('### Quality covariates (recorded, NOT filtered on)')
        A('')
        A('No genome was removed for proteome size or assembly level. '
          'PROJECT_STATE.md 4.1/12.2/12.3: a protein-count filter deletes '
          'reduced-genome lineages (it removed all 7 *Halopseudomonas* from the '
          '389-genome *Pseudomonas* set) and an assembly-level filter does not '
          'control gene count. Quality is a covariate here, analysed afterwards.')
        A('')
        A('| genus | complete | chromosome | scaffold | contig | proteins min/median/max | '
          'reused | downloaded | no annotation |')
        A('|---|---:|---:|---:|---:|---|---:|---:|---:|')
        for g in [r['genus'] for r in status_rows if r['genus'] in res]:
            r = res[g]
            A('| %s | %s | %s | %s | %s | %s / %s / %s | %s | %s | %s |'
              % (g, r.get('asm_complete', ''), r.get('asm_chromosome', ''),
                 r.get('asm_scaffold', ''), r.get('asm_contig', ''),
                 r.get('proteins_min', ''), r.get('proteins_median', ''),
                 r.get('proteins_max', ''), r.get('proteomes_reused', ''),
                 r.get('proteomes_downloaded', ''), r.get('species_no_annotation', '')))
        A('')

        A('### Fast-AAI boundary assertion (condition of use, PROJECT_STATE 6.2a)')
        A('')
        b65 = sum(int(r.get('boundary_flags_65') or 0) for r in res.values())
        b95 = sum(int(r.get('boundary_flags_95') or 0) for r in res.values())
        A('Pairs whose stored AAI sits exactly on a threshold, and could therefore in '
          'principle bin differently between the original and fast AAI paths: '
          '**%d at 65.00**, **%d at 95.00** across %d completed genera.'
          % (b65, b95, len(res)))
        if b65 or b95:
            A('')
            A('> Flagged pairs are listed in `reports/boundary_flags.tsv` and need a '
              'full-precision check before their genus verdict is reported.')
        A('')

    # ---- failures --------------------------------------------------------
    fails = read_tsv(FAILURES)
    if fails:
        A('## Failures')
        A('')
        A('| when (UTC) | genus | reason |')
        A('|---|---|---|')
        for r in fails[-20:]:
            A('| %s | %s | %s |' % (r.get('utc', ''), r.get('genus', ''),
                                    r.get('reason', '')))
        A('')

    # ---- ETA -------------------------------------------------------------
    A('## Remaining time - EXTRAPOLATED, NOT MEASURED')
    A('')
    walls = [(r['genus'], float(r['wall_seconds']))
             for r in result_rows if r.get('wall_seconds')]
    if not walls:
        A('_No genus has finished yet, so there is nothing to extrapolate from._ '
          'No estimate is given rather than a guessed one - this project does not '
          'report numbers it has not observed.')
    else:
        num = sum(w for _g, w in walls)
        den = sum((proteins.get(g, 0) / 1e6) ** 2 for g, _w in walls)
        remaining = [r['genus'] for r in status_rows
                     if r['state'] in ('pending', 'running', 'fetching')]
        A('Observed wall times on the %d genus/genera already finished:' % len(walls))
        A('')
        A('| genus | proteins (NCBI counts) | observed wall |')
        A('|---|---:|---:|')
        for g, w in walls:
            A('| %s | %s | %s |' % (g, '{:,}'.format(proteins.get(g, 0)), fmt_dur(w)))
        A('')
        if den > 0 and remaining:
            k = num / den
            est = sum(k * (proteins.get(g, 0) / 1e6) ** 2 for g in remaining)
            A('Extrapolating those to the %d genera not yet finished, on the assumption '
              'that an all-against-all search costs ~(total proteins)^2:'
              % len(remaining))
            A('')
            A('- **extrapolated remaining wall time: ~%s**' % fmt_dur(est))
            A('')
            A('This is an EXTRAPOLATION from %d observed run(s). It is not a '
              'measurement, and this project has abandoned finish-time predictions '
              'before (PROJECT_STATE 12.3). Treat it as an order of magnitude only.'
              % len(walls))
        else:
            A('Nothing left to extrapolate to.')
    A('')
    A('---')
    A('')
    A('Audit trail: `reports/events.log` (append-only). Restart instructions: '
      '`HOW_TO_RESUME.md`.')
    A('')
    return '\n'.join(L)


def cycle(seen):
    status_rows = read_tsv(STATUS)
    result_rows = read_tsv(RESULTS)
    proteins = genus_proteins()

    # independent observation of state changes
    res = {r['genus']: r for r in result_rows}
    for r in status_rows:
        key = r['genus']
        cur = r['state']
        if seen.get(key) != cur:
            if seen.get(key) is not None or cur != 'pending':
                detail = r.get('note', '')
                event('OBS_STATE', key, '%s -> %s  %s' % (seen.get(key, '?'), cur, detail))
                if cur == 'done' and key in res:
                    v = res[key]
                    event('OBS_%s' % v.get('agrees_with_expectation', 'NA'), key,
                          'expected=%s observed=%s bins=%s core=%s'
                          % (v.get('expectation'), v.get('verdict'),
                             v.get('bins_at_65'), v.get('core_genes')))
            seen[key] = cur

    text = build(status_rows, result_rows, proteins)
    tmp = PROGRESS + '.tmp'
    with open(tmp, 'w') as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, PROGRESS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--interval', type=int, default=120)
    ap.add_argument('--once', action='store_true')
    args = ap.parse_args()

    os.makedirs(REPORTS, exist_ok=True)
    lock = open(LOCK, 'w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print('another tracker is already running - exiting')
        return 0
    with open(PIDFILE, 'w') as fh:
        fh.write('%d\n' % os.getpid())

    # State is rebuilt from disk, so a restart re-announces nothing it already knows.
    seen = {r['genus']: r['state'] for r in read_tsv(STATUS)}
    event('TRACKER_START', '-', 'pid=%d interval=%ds' % (os.getpid(), args.interval))
    cycle(seen)
    if args.once:
        return 0
    while True:
        time.sleep(args.interval)
        try:
            cycle(seen)
        except Exception as exc:                # noqa: BLE001
            event('TRACKER_ERROR', '-', str(exc))


if __name__ == '__main__':
    raise SystemExit(main())
