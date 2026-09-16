#!/bin/bash
# Start (or restart) the LPSN genus survey, fully detached.
#
# Safe to run repeatedly and safe to run while it is already running: both the driver
# and the tracker take an flock, so a second copy exits immediately rather than
# doubling up. A genus already marked `done` is never recomputed.
#
# setsid + nohup means neither process is a child of this shell, so both survive the
# terminal closing, the ssh session dropping, and the interactive session ending. That is
# the whole point: an AI session cannot be the durable layer.
#
#   bash scripts/start_survey.sh            # threads default to 140
#   bash scripts/start_survey.sh 100        # override thread count

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
THREADS="${1:-140}"
cd "$ROOT" || exit 1

mkdir -p logs reports results proteomes

# The manifest is cheap and deterministic; regenerate it if it is missing.
if [ ! -s reports/genome_selection.tsv ]; then
    python3 scripts/select_genomes.py || exit 1
fi

start_one () {           # name script logfile extra-args...
    local name="$1"; shift
    local script="$1"; shift
    local logfile="$1"; shift
    local pidfile="$ROOT/reports/${name}.pid"
    if [ -f "$pidfile" ] && kill -0 "$(cat "$pidfile" 2>/dev/null)" 2>/dev/null; then
        echo "$name already running (pid $(cat "$pidfile")) - leaving it alone"
        return 0
    fi
    setsid nohup python3 -u "$script" "$@" >> "$logfile" 2>&1 < /dev/null &
    disown 2>/dev/null
    echo "$name started, logging to $logfile"
}

start_one driver  "$ROOT/scripts/run_survey.py"       "$ROOT/logs/driver.log"  -t "$THREADS"
sleep 2
start_one tracker "$ROOT/scripts/progress_tracker.py" "$ROOT/logs/tracker.log" --interval 120
sleep 3

echo
echo "--- reports/survey_status.tsv ---"
cat "$ROOT/reports/survey_status.tsv" 2>/dev/null
echo
echo "driver pid : $(cat "$ROOT/reports/driver.pid" 2>/dev/null)"
echo "tracker pid: $(cat "$ROOT/reports/tracker.pid" 2>/dev/null)"
echo
echo "watch progress with:  cat $ROOT/reports/PROGRESS.md"
