#!/bin/bash
# Same-plane CSA comparison of two csa-bifurcation runs (baseline vs edited t0 surface, same registration).
#   run_samplane.sh <baseline_subject_dir> <edited_subject_dir> <out_dir> [workers]
# -> <out_dir>/samplane_base_vs_edit.csv   (then plot_samplane.py)
set -e
B=$1; E=$2; OUT=$3; NW=${4:-12}
HERE=$(cd "$(dirname "$0")" && pwd); PVB=$HERE/../pvb.sh
mkdir -p "$OUT/logs"; rm -f "$OUT"/recut_*.csv
for w in $(seq 0 $((NW - 1))); do
  $PVB $HERE/samplane_recut.py --base-dir "$B" --edit-dir "$E" --out-dir "$OUT" --shard $w:$NW > "$OUT/logs/recut_$w.log" 2>&1 &
done
wait
{ head -1 "$OUT/recut_00.csv"; for f in "$OUT"/recut_*.csv; do tail -n +2 "$f"; done; } > "$OUT/samplane_base_vs_edit.csv"
python3 -c "
import pandas as pd; d = pd.read_csv('$OUT/samplane_base_vs_edit.csv'); r = d.area_base / d.area_table
print('%d frames, %d planes; baseline re-cut / slicer table: median %.4f, %.1f%% of planes within 1%%' % (d.ms.nunique(), len(d), r.median(), 100 * ((r - 1).abs() < 0.01).mean()))"
echo "$OUT/samplane_base_vs_edit.csv"
