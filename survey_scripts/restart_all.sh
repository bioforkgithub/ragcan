#!/bin/bash
# Bring back ALL FOUR survey processes: driver, tracker, watchdog, load guard.
#
# start_survey.sh starts only the driver and tracker. This wrapper adds the watchdog
# (restarts driver/tracker/guard if they die) and the load guard (steps threads down if
# the machine strains). All four take flocks, so running this while they are already
# running is a safe no-op - each duplicate exits immediately.
#
# WHY IT EXISTS: the machine rebooted twice on 2026-08-26 and nothing brought the survey
# back. Cause never established (journalctl is root-only here). This is invoked from
# @reboot in crontab so a reboot costs at most one genus, not the whole run.
#
#   bash scripts/restart_all.sh [threads]     # threads default: reports/THREADS, else 96

set -u

# --- PATH FIX, 2026-08-28 -------------------------------------------------------
# cron runs with a minimal PATH where python3 is /usr/bin/python3, which has NO numpy.
# On 2026-08-28 that silently failed all 65 archaeal genera with
# "driver exception: No module named 'numpy'" - the bacterial genera had succeeded only
# because the driver happened to be launched from an interactive shell. Pin the
# interpreter explicitly so it cannot depend on who started the process.
export PATH="/localsoftware/anaconda3/bin:$PATH"
# --------------------------------------------------------------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 1
T="${1:-$(cat reports/THREADS 2>/dev/null)}"
case "$T" in ''|*[!0-9]*) T=96 ;; esac
echo "$T" > reports/THREADS

echo "$(date -u +%FT%TZ) restart_all: threads=$T"
bash "$ROOT/scripts/start_survey.sh" "$T" 2>&1 | sed 's/^/  /'

for s in watchdog load_guard; do
    pid=$(cat "$ROOT/reports/$s.pid" 2>/dev/null)
    if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then
        echo "  $s already running (pid $pid)"
    else
        setsid nohup bash "$ROOT/scripts/$s.sh" >/dev/null 2>&1 < /dev/null &
        disown 2>/dev/null
        sleep 3
        np=$(cat "$ROOT/reports/$s.pid" 2>/dev/null)
        [ -n "$np" ] && kill -0 "$np" 2>/dev/null \
            && echo "  $s started (pid $np)" || echo "  $s FAILED to start"
    fi
done
echo "$(date -u +%FT%TZ) restart_all: done"
