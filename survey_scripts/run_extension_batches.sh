#!/bin/bash
# Batch wrapper for the two-to-three-species extension. Runs run_survey.py over the
# extension genera in batches of 25, ONE driver at a time, at nice 19 with 16 threads.
# Before each batch: waits while the 1-minute load average is above 200.
# Does not touch any process it did not start. All paths are relative to the run folder.
set -u
umask 077
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 1
echo $$ > reports/batches.pid
BATCH=25
# same order as the driver: smallest genus first, then name
mapfile -t ORDER < <(python3 -c "
import csv,collections
c=collections.Counter(r['genus'] for r in csv.DictReader(open('reports/genome_selection.tsv'),delimiter='\t'))
[print(g) for g,n in sorted(c.items(), key=lambda kv:(kv[1],kv[0]))]")
echo "$(date -u +%FT%TZ) batches start: ${#ORDER[@]} genera"
for ((i=0; i<${#ORDER[@]}; i+=BATCH)); do
  while :; do
    l=$(cut -d' ' -f1 /proc/loadavg)
    if awk -v l="$l" 'BEGIN{exit !(l>200)}'; then
      echo "$(date -u +%FT%TZ) load1=$l > 200, waiting 60s"; sleep 60
    else break; fi
  done
  batch=("${ORDER[@]:i:BATCH}")
  echo "$(date -u +%FT%TZ) batch $((i/BATCH+1)) load1=$l: ${batch[0]} .. ${batch[-1]}"
  nice -n 19 python3 -u scripts/run_survey.py -t 16 --only "${batch[@]}" >> logs/driver.log 2>&1
done
echo "$(date -u +%FT%TZ) batches finished"
