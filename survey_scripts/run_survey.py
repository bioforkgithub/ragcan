#!/usr/bin/env python3
"""
Phase 1 driver for the LPSN genus-coherence survey.

Runs one RaGCAn screen per genus, smallest genus first, resumably, in the background.
Safe to run repeatedly: a genus already marked `done` is never recomputed.

RESUME DESIGN - this driver assumes it WILL be killed mid-run, repeatedly.

  Genus level   reports/survey_status.tsv is the resume point. The state transition is
                written to disk BEFORE the work starts, so a genus killed mid-run is on
                disk as `running`, not silently back at `pending`. On restart a stale
                `fetching`/`running` is treated as resumable and said so in the log.
                The file is written temp-then-os.replace(), so a kill mid-write cannot
                corrupt the one file that must never be lost.

  Mid-genus     pygemini checkpoints its DIAMOND search with
                results/<Genus>/working/all_vs_all.tsv.complete and reuses it
                (pygemini.py:225). working/ is therefore NOT deleted until the genus is
                confirmed done. On a resume the output directory is non-empty, so
                pygemini needs --force to proceed - and --force deliberately KEEPS
                working/ (prepare_output_directory(..., keep_working=True),
                pygemini.py:687), so it discards partial results, never the finished
                search. --rerun-search is NEVER passed; passing it would throw away
                hours of completed DIAMOND work.
                Honest limit: a search killed *while running* has no .complete marker
                and is redone from the start. That is pygemini's own documented
                behaviour ("DIAMOND writes progressively and cannot resume"), not
                something this driver can fix.

  Downloads     see fetch_proteomes.py - atomic rename, gzip verified in full.

  Concurrency   an flock on reports/driver.lock means a second driver cannot start on
                top of a live one.

FIXED PARAMETERS, identical for every genus so the results are comparable:
  --pident 60 --bin-aai 65, DIAMOND very-sensitive (pygemini's default), -t THREADS.
  Threads default to 140, not 220: a SqueezeMeta coassembly for an unrelated project
  holds ~100 of this machine's 256 cores and must not be starved.

NO QUALITY FILTER IS APPLIED - see select_genomes.py and PROJECT_STATE.md 4.1/12.2/12.3.
Assembly level and protein count are recorded as covariates, not used to drop genomes.
"""

import argparse
import csv
import datetime
import fcntl
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROJ = os.path.dirname(ROOT)

REPORTS = os.path.join(ROOT, 'reports')
RESULTS = os.path.join(ROOT, 'results')
LOGS = os.path.join(ROOT, 'logs')
PROTEOMES = os.path.join(ROOT, 'proteomes')

STATUS = os.path.join(REPORTS, 'survey_status.tsv')
EVENTS = os.path.join(REPORTS, 'events.log')
FAILURES = os.path.join(REPORTS, 'failures.tsv')
LOCK = os.path.join(REPORTS, 'driver.lock')
PIDFILE = os.path.join(REPORTS, 'driver.pid')
SELECTION = os.path.join(REPORTS, 'genome_selection.tsv')

RUN_FAST = os.path.join(PROJ, 'AAI_Optimisation', 'run_fast.py')

STATUS_COLS = ['genus', 'n_species', 'state', 'attempts', 'started_utc',
               'updated_utc', 'wall_seconds', 'note']

sys.path.insert(0, HERE)


def utc():
    return datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')


def event(kind, genus, detail=''):
    """Append one line to the audit trail. Append-only: a crash cannot truncate it."""
    with open(EVENTS, 'a') as fh:
        fh.write('%s\t%s\t%s\t%s\n' % (utc(), kind, genus, detail))
        fh.flush()
        os.fsync(fh.fileno())


def log(msg):
    print('%s  %s' % (utc(), msg), flush=True)


# ---------------------------------------------------------------- status file

def read_status():
    if not os.path.exists(STATUS):
        return {}
    with open(STATUS, newline='') as fh:
        return {r['genus']: r for r in csv.DictReader(fh, delimiter='\t')}


def write_status(rows):
    """Atomic: temp file, fsync, rename. The resume point is never half-written."""
    tmp = STATUS + '.tmp'
    with open(tmp, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=STATUS_COLS, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, '') for c in STATUS_COLS})
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, STATUS)


def set_state(order, status, genus, state, **kw):
    r = status[genus]
    r['state'] = state
    r['updated_utc'] = utc()
    r.update({k: v for k, v in kw.items()})
    write_status([status[g] for g in order])


# ---------------------------------------------------------------- genus list

def genus_order():
    """Smallest genus first, so signal arrives early. Pseudomonas (391) is last."""
    counts = {}
    with open(SELECTION, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            counts[r['genus']] = counts.get(r['genus'], 0) + 1
    return sorted(counts.items(), key=lambda kv: (kv[1], kv[0]))


def init_status(order):
    status = read_status()
    changed = False
    for genus, n in order:
        if genus not in status:
            status[genus] = {'genus': genus, 'n_species': n, 'state': 'pending',
                             'attempts': '0', 'started_utc': '', 'updated_utc': utc(),
                             'wall_seconds': '', 'note': ''}
            changed = True
        else:
            if status[genus].get('n_species') != str(n):
                status[genus]['n_species'] = n
                changed = True
    if changed or not os.path.exists(STATUS):
        write_status([status[g] for g, _ in order])
    return status


# ---------------------------------------------------------------- run one genus

def run_ok(genus):
    """A pygemini run counts as successful ONLY on both of these, per the project rule:
    reports/summary.txt exists AND the log reached '[9] writing results'.
    `rc=` in driver logs is meaningless in this project (PROJECT_STATE 5.1); the
    subprocess return code IS checked here because it is read directly, not through
    a shell `$?` that $(date) has already reset."""
    summary = os.path.join(RESULTS, genus, 'reports', 'summary.txt')
    if not os.path.exists(summary):
        return False, 'reports/summary.txt missing'
    logpath = os.path.join(LOGS, '%s.log' % genus)
    try:
        with open(logpath, errors='replace') as fh:
            if '[9] writing results' not in fh.read():
                return False, "log never reached '[9] writing results'"
    except OSError as exc:
        return False, 'log unreadable: %s' % exc
    return True, ''


def do_genus(genus, order, status, threads, keep_working):
    from fetch_proteomes import fetch_genus
    import summarise_genus

    rec = status[genus]
    rec['attempts'] = str(int(rec.get('attempts') or 0) + 1)

    # --- fetch -------------------------------------------------------------
    set_state([g for g, _ in order], status, genus, 'fetching', started_utc=utc(), note='')
    event('FETCH_START', genus)
    indir, counts = fetch_genus(genus, log=log)
    event('FETCH_DONE', genus, str(counts))

    if counts['usable'] < 2:
        set_state([g for g, _ in order], status, genus, 'failed',
                  note='only %d usable proteomes' % counts['usable'])
        event('GENUS_FAILED', genus, 'only %d usable proteomes' % counts['usable'])
        return False

    # --- run ---------------------------------------------------------------
    outdir = os.path.join(RESULTS, genus)
    logpath = os.path.join(LOGS, '%s.log' % genus)
    resuming = os.path.isdir(outdir) and bool(os.listdir(outdir))
    marker = os.path.join(outdir, 'working', 'all_vs_all.tsv.complete')
    if resuming:
        log('[run] %s: output directory is not empty - RESUMING. --force is passed so '
            'partial results are cleared, but pygemini keeps working/ (keep_working=True), '
            'so the finished DIAMOND search is reused. --rerun-search is NOT passed.'
            % genus)
        log('[run] %s: finished-search checkpoint present: %s'
            % (genus, os.path.exists(marker)))

    cmd = [sys.executable, RUN_FAST, '-i', indir, '-o', outdir, '-t', str(threads),
           '--pident', '60', '--bin-aai', '65', '--log-file', logpath]
    if resuming:
        cmd.append('--force')

    set_state([g for g, _ in order], status, genus, 'running',
              started_utc=utc(), note='attempt %s' % rec['attempts'])
    event('RUN_START', genus, '%d genomes, threads=%d, resume=%s'
          % (counts['usable'], threads, resuming))
    log('[run] %s: %s' % (genus, ' '.join(cmd)))

    t0 = time.time()
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    wall = int(time.time() - t0)

    ok, why = run_ok(genus)
    if proc.returncode != 0 and not ok:
        why = why or 'exit code %d' % proc.returncode
    if not ok:
        set_state([g for g, _ in order], status, genus, 'failed',
                  wall_seconds=wall, note=why)
        record_failure(genus, why, wall)
        event('GENUS_FAILED', genus, why)
        log('[run] %s: FAILED - %s' % (genus, why))
        return False

    # --- summarise ---------------------------------------------------------
    try:
        row = summarise_genus.summarise(genus, wall_seconds=wall)
    except Exception as exc:                    # noqa: BLE001
        set_state([g for g, _ in order], status, genus, 'failed',
                  wall_seconds=wall, note='summarise failed: %s' % exc)
        record_failure(genus, 'summarise failed: %s' % exc, wall)
        event('GENUS_FAILED', genus, 'summarise failed: %s' % exc)
        return False

    set_state([g for g, _ in order], status, genus, 'done',
              wall_seconds=wall,
              note='%s bins=%s core=%s %s' % (row['verdict'], row['bins_at_65'],
                                              row['core_genes'],
                                              row['agrees_with_expectation']))
    event('GENUS_DONE', genus,
          'genomes=%s core=%s bins=%s min_aai=%s mean_aai=%s pairs>=95=%s '
          'first_split=%s expected=%s %s wall=%ds'
          % (row['genomes_used'], row['core_genes'], row['bins_at_65'],
             row['min_aai'], row['mean_aai'], row['pairs_ge_95'],
             row['first_split_threshold'], row['expectation'],
             row['agrees_with_expectation'], wall))
    log('[run] %s: DONE in %ds - %s bins, %s core genes, %s'
        % (genus, wall, row['bins_at_65'], row['core_genes'],
           row['agrees_with_expectation']))

    # --- cleanup, only AFTER done is on disk -------------------------------
    if not keep_working:
        work = os.path.join(outdir, 'working')
        if os.path.isdir(work):
            try:
                size = sum(os.path.getsize(os.path.join(work, f))
                           for f in os.listdir(work)
                           if os.path.isfile(os.path.join(work, f)))
                shutil.rmtree(work)
                log('[run] %s: removed working/ (%.2f GB of transient hit files)'
                    % (genus, size / 1e9))
            except OSError as exc:
                log('[run] %s: could not remove working/: %s' % (genus, exc))
    return True


def record_failure(genus, why, wall):
    new = not os.path.exists(FAILURES)
    with open(FAILURES, 'a', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        if new:
            w.writerow(['utc', 'genus', 'wall_seconds', 'reason'])
        w.writerow([utc(), genus, wall, why])


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-t', '--threads', type=int, default=140)
    ap.add_argument('--only', nargs='*', help='restrict to these genera')
    ap.add_argument('--keep-working', action='store_true',
                    help='do not delete working/ after a genus is done')
    ap.add_argument('--retry-failed', action='store_true',
                    help='also re-attempt genera currently marked failed')
    args = ap.parse_args()

    for d in (REPORTS, RESULTS, LOGS, PROTEOMES):
        os.makedirs(d, exist_ok=True)

    lock = open(LOCK, 'w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        log('another driver already holds %s - exiting. This is the safe outcome; '
            'the running driver will continue.' % LOCK)
        return 0
    with open(PIDFILE, 'w') as fh:
        fh.write('%d\n' % os.getpid())

    order = genus_order()
    status = init_status(order)
    event('DRIVER_START', '-', 'pid=%d threads=%d' % (os.getpid(), args.threads))
    log('driver start, pid %d, threads %d' % (os.getpid(), args.threads))
    log('genus order (smallest first): %s'
        % ', '.join('%s(%d)' % (g, n) for g, n in order))

    # report what resume found before doing anything
    for genus, _n in order:
        st = status[genus]['state']
        if st in ('running', 'fetching'):
            log('RESUME: %s was left in state "%s" by a previous run (killed mid-work). '
                'Treating as resumable and re-entering it.' % (genus, st))
            event('RESUME_STALE', genus, 'was %s' % st)
        elif st == 'done':
            log('RESUME: %s already done - will not be recomputed.' % genus)

    todo = []
    for genus, _n in order:
        if args.only and genus not in args.only:
            continue
        st = status[genus]['state']
        if st == 'done':
            continue
        if st == 'failed' and not args.retry_failed:
            log('SKIP: %s is marked failed (%s). Re-run with --retry-failed to retry.'
                % (genus, status[genus].get('note', '')))
            continue
        todo.append(genus)

    log('%d genera to process: %s' % (len(todo), ', '.join(todo)))

    for genus in todo:
        try:
            do_genus(genus, order, status, args.threads, args.keep_working)
        except KeyboardInterrupt:
            log('interrupted during %s - state on disk is "%s", re-run to resume'
                % (genus, status[genus]['state']))
            event('DRIVER_INTERRUPTED', genus)
            raise
        except Exception as exc:                # noqa: BLE001
            set_state([g for g, _ in order], status, genus, 'failed',
                      note='driver exception: %s' % exc)
            record_failure(genus, 'driver exception: %s' % exc, '')
            event('GENUS_FAILED', genus, 'driver exception: %s' % exc)
            log('[run] %s: driver exception: %s' % (genus, exc))

    done = sum(1 for g, _ in order if status[g]['state'] == 'done')
    log('driver finished: %d/%d genera done' % (done, len(order)))
    event('DRIVER_FINISHED', '-', '%d/%d done' % (done, len(order)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
