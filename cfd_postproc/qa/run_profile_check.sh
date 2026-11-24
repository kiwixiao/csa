#!/bin/bash
# Prescribed-profile sanity check for one subject: per-step boundary series (parallel shards) + plot.
#   run_profile_check.sh <case> <flow_profile.csv> <pressure_profile.csv> <outdir> <SUBJECT> [workers] [extra plot args...]
set -e
CASE=$1; FLOW=$2; PRES=$3; OUT=$4; SUBJ=$5; NW=${6:-3}; shift 6 || true
HERE=$(cd "$(dirname "$0")" && pwd); PVB=$HERE/../pvb.sh
DIR=$(dirname "$CASE")
BUSY=$(find "$DIR" -maxdepth 1 -type f -mmin -3 | wc -l)
if [ "$BUSY" -gt 0 ]; then
  echo "ERROR: $BUSY files in $DIR changed in the last 3 min - an export/copy is still in progress. Re-run when it has finished."; exit 2
fi
rm -rf "$OUT/shards"; mkdir -p "$OUT/shards"
N=$(grep -A2 'number of steps' "$CASE" | head -1 | awk '{print $NF}')
CH=$(( (N + NW - 1) / NW ))
for w in $(seq 0 $((NW - 1))); do
  $PVB "$HERE/cfd_boundary_series.py" --case "$CASE" --out "$OUT/shards/part$w.csv" --steps $((w * CH)):$(((w + 1) * CH)) > "$OUT/shards/part$w.log" 2>&1 &
done
wait
{ head -1 "$OUT/shards/part0.csv"; for f in "$OUT"/shards/part*.csv; do tail -n +2 "$f"; done | sort -t, -k1,1n; } > "$OUT/${SUBJ}_cfd_boundary_series.csv"
python3 "$HERE/plot_profiles_vs_cfd.py" --subject "$SUBJ" --flow "$FLOW" --pressure "$PRES" \
    --cfd "$OUT/${SUBJ}_cfd_boundary_series.csv" --out "$OUT/${SUBJ}_profiles_vs_cfd.png" "$@"
