#!/bin/bash
# Watchdog for the LPSN genus-coherence survey.
#
# WHY THIS EXISTS: an interactive session died overnight on 2026-08-25/26 and work stalled
# unnoticed. An AI agent cannot be the surveillance layer — it dies with its session.
# This is a plain shell loop with no dependency on any interactive session. It outlives every session.
#
# WHAT IT DOES, every 60s:
#   - is the driver alive?  if not, and the survey is not finished, restart it
#   - is the tracker alive? if not, restart it
#   - record every check and every action in watchdog.log
#
# It is idempotent and safe to run twice (it takes a lock). Restarting the driver is
# safe by design: the driver is resumable and never recomputes a genus already `done`.
#
# START:  setsid nohup bash scripts/watchdog.sh >/dev/null 2>&1 &
# STOP:   kill $(cat reports/watchdog.pid)

S=/data/prosjekt/15719-Res-Marine/Project_Work/Project_Pygemini/LPSN_Genus_Survey
R="$S/reports"
LOG="$R/watchdog.log"
cd "$S" || exit 1

# single-instance lock
exec 9>"$R/watchdog.lock"
flock -n 9 || { echo "$(date -u +%FT%TZ) another watchdog holds the lock; exiting" >>"$LOG"; exit 0; }
echo $$ > "$R/watchdog.pid"

log(){ echo "$(date -u +%FT%TZ) $*" >> "$LOG"; }

alive(){ # $1 = pidfile
  local p; p=$(cat "$1" 2>/dev/null) || return 1
  [ -n "$p" ] && kill -0 "$p" 2>/dev/null
}

survey_finished(){
  # finished when no genus is left in pending/running/fetching
  [ -f "$R/survey_status.tsv" ] || return 1
  ! awk -F'\t' 'NR>1 && ($0 ~ /pending|running|fetching/){found=1} END{exit !found}' "$R/survey_status.tsv"
}

log "watchdog START pid=$$"

while true; do
  if survey_finished; then
    log "survey finished — all genera terminal. watchdog exiting."
    # leave the tracker running one last cycle so PROGRESS.md reflects the final state
    sleep 150 9>&-
    exit 0
  fi

  if ! alive "$R/tracker.pid"; then
    log "TRACKER DEAD — restarting"
    setsid nohup python3 -u "$S/scripts/progress_tracker.py" >> "$S/logs/tracker.out" 2>&1 9>&- &
    sleep 3 9>&-
    alive "$R/tracker.pid" && log "tracker restarted pid=$(cat "$R/tracker.pid")" \
                           || log "tracker restart FAILED — see logs/tracker.out"
  fi

  if ! alive "$R/load_guard.pid"; then
    log "LOAD GUARD DEAD — restarting (it is the only thing capping core usage)"
    setsid nohup bash "$S/scripts/load_guard.sh" >/dev/null 2>&1 9>&- &
    sleep 3 9>&-
    alive "$R/load_guard.pid" && log "load_guard restarted pid=$(cat "$R/load_guard.pid")" \
                              || log "load_guard restart FAILED — see reports/load_guard.log"
  fi

  if ! alive "$R/driver.pid"; then
    log "DRIVER DEAD and survey unfinished — restarting (resumable; completed genera are skipped)"
    # No --resume flag exists: the driver resumes automatically from survey_status.tsv.
    # --retry-failed is the documented restart form (HOW_TO_RESUME.md) and also re-attempts
    # any genus left `failed` by a kill. Threads must match the original run.
    T=$(cat "$R/THREADS" 2>/dev/null); case "$T" in ''|*[!0-9]*) T=96 ;; esac
    log "restarting driver with -t $T (from reports/THREADS)"
    setsid nohup python3 -u "$S/scripts/run_survey.py" --retry-failed -t "$T" >> "$S/logs/driver.out" 2>&1 9>&- &
    sleep 5 9>&-
    alive "$R/driver.pid" && log "driver restarted pid=$(cat "$R/driver.pid")" \
                          || log "driver restart FAILED — see logs/driver.out"
  fi

  sleep 60 9>&-
done
