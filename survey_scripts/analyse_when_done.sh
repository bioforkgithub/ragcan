#!/bin/bash
# Run the full Phase 3 analysis automatically, ONCE, when the survey finishes.
#
# WHY THIS EXISTS: Manish asked for "all the analysis done by Tuesday". An AI session cannot do
# that - it ends when the conversation ends. This is a plain shell script on cron, so the
# analysis runs whether or not anyone is at a keyboard. Same reasoning as the watchdog.
#
# It does nothing until every genus is terminal (done/failed) AND the driver has exited.
# It writes a stamp file so it runs once, not hourly forever.
#
# CRON: 32 * * * * bash scripts/analyse_when_done.sh >> logs/analysis.log 2>&1

set -u
S=/data/prosjekt/15719-Res-Marine/Project_Work/Project_Pygemini/LPSN_Genus_Survey
cd "$S" || exit 1
STAMP="$S/phase3/ANALYSIS_COMPLETE.stamp"
OUT="$S/phase3/FINAL_ANALYSIS.txt"
log(){ echo "$(date -u +%FT%TZ) $*"; }

[ -f "$STAMP" ] && exit 0                       # already done

# 1. any genus still pending/running/fetching? then not finished.
if awk -F'\t' 'NR>1 && ($3=="pending"||$3=="running"||$3=="fetching"){f=1} END{exit !f}' \
        reports/survey_status.tsv; then
    exit 0
fi
# 2. driver still alive? give it time to write its last results.
p=$(cat reports/driver.pid 2>/dev/null)
[ -n "$p" ] && kill -0 "$p" 2>/dev/null && { log "all genera terminal but driver still alive; waiting"; exit 0; }
# 3. are the archaea even queued yet? if the manifest has genera the status file
#    has never seen, the survey is NOT finished - a later driver start will run them.
mg=$(awk -F'\t' 'NR>1{g[$1]=1} END{print length(g)}' reports/genome_selection.tsv)
sg=$(awk -F'\t' 'NR>1{n++} END{print n+0}' reports/survey_status.tsv)
if [ "$mg" -gt "$sg" ]; then
    log "manifest has $mg genera but status has $sg - archaea not queued yet; waiting"
    exit 0
fi

log "SURVEY FINISHED - running Phase 3 analysis"
{
  echo "FINAL PHASE 3 ANALYSIS"
  echo "generated $(date -u +%FT%TZ) by scripts/analyse_when_done.sh, unattended"
  echo "=================================================================="
  echo
  echo "Genera scored: $sg   Genomes: $(awk -F'\t' 'NR>1{n++} END{print n}' reports/genome_selection.tsv)"
  awk -F'\t' 'NR>1{c[$3]++} END{for(s in c) printf "  %-9s %d\n", s, c[s]}' reports/survey_status.tsv
  echo
  echo "### 1. AUTHORITATIVE: partition agreement with GTDB (adjusted Rand index)"
  echo "###    PLAN.md line 115 specifies this method. Report THESE numbers."
  python3 phase3/agreement.py 2>&1
  echo
  echo "### 2. Screening metrics - the 'first line of identification' framing"
  python3 phase3/screen_metrics.py 2>&1
  echo
  echo "### 3. Threshold-margin triage of genera RaGCAn split alone"
  python3 phase3/triage_oversplits.py 2>&1
  echo
  echo "### 4. Do singleton bins match GTDB reassignments?"
  python3 phase3/singleton_check.py 2>&1
  echo
  echo "### 5. SUPERSEDED 2x2 - kept only for comparison. Its >=2-members rule"
  echo "###    inverted the headline once. See MANUSCRIPT_CHANGES.md section P."
  python3 phase3/gtdb_compare.py 2>&1
} > "$OUT" 2>&1

date -u +%FT%TZ > "$STAMP"
log "analysis written to $OUT"
