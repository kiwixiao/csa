#!/bin/bash
# Per-subject CFD post-processing (generalised from /storage4/ENT003/run_postproc_all.sh).
#
#   run_subject.sh <subjects/SUBJ.env> [stage ...]
#
# Upstream stages (run only when named, or with 'all'):
#   csa        CSA from the registration: csa-bifurcation on SUBJ_DIR (conda env mirtk; step 0 calls env vmtk)
#   facecheck  guard: face labels == faces of every deformed frame (else the slicer silently froze CSA at t0)
#   motioncheck does the CFD wall follow the registration? (REG_CSV control points vs CFD wall, every step)
#   velscan    per-step velocity percentiles -> VELSCAN (pick VMAX from it)
#   clean      nostril end of the CSA tables: off-centerline planes rejected, cap -> first plane fill planes
#              (after extract; writes POST/csa_clean = CSA_PLOT / PE_IN / PLANES_DIR)
#   all        = csa facecheck motioncheck velscan extract clean velgrid plots regions pdf vellight reportvideo
#
# Stages (default, in this order):
#   extract    CSA-plane total pressure / energy flux per region      (pvbatch, ~1 h, parallel)
#   velgrid    velocity sampling grid only (velocity_breath_final/manifest.json: bounds + dims shared by
#              the white renders); no frames
#   plots      anatomy / pressure / WOB plots per path (wob_summary.csv sets the peak-WOB frame)
#   regions    work of breathing by region and per segment, local loss at peak power
#   pdf        CFD report, one page per instant (locator / loss-coloured centerline / velocity, plots, WOB per segment)
#              (<S>_CFD_report.pdf, .docx)
#   vellight   velocity volume render on white, every CFD step (for the report video; pvbatch, 6 shards)
#   reportvideo the report page animated over the breath (<S>_report_breath.mp4)
#
# Old-format stages (superseded by the report + report video; run only when named):
#   velrender  velocity volume render, dark, every CFD step (also writes the velgrid manifest)
#   stlplanes  moving wall + CSA planes render, dark, every frame
#   wobgradient standalone WOB gradient plots (<S>_WOB_gradient.png, <S>_WOB_gradient_all.png)
#   refs       centerline locators at peak inhalation
#   video      breath multipanel video (needs velrender stlplanes refs)
#   csavideo   CSA breath video (needs stlplanes refs)
#   snapshot   hi-res light snapshot at peak dissipated power (needs refs)
#   onepager   text one-pager (<S>_CFD_report_text.pdf; needs snapshot, after regions)
#
# Fixes vs the ENT003 driver: plane-STL folder passed explicitly (normals), density per
# subject (CFD_RHO), locator frame + peak velocity computed, dark STL-plane render included.
set -e
CFG=$1; shift
[ -f "$CFG" ] || { echo "usage: $0 <subject.env> [stage ...]"; exit 1; }
source "$CFG"
STAGES=${*:-extract velgrid plots regions pdf vellight reportvideo}
[ "$STAGES" = all ] && STAGES="csa facecheck motioncheck velscan extract clean velgrid plots regions pdf vellight reportvideo"
HERE=$(cd "$(dirname "$0")" && pwd)
PVB=$HERE/pvb.sh
export CFD_RHO=$RHO SUBJECT=$SUBJ PYTHONPATH=$HERE
EXTRACT=$POST/csa_pe; PANELS=$POST/panels; VEL=$POST/velocity_breath_final
# stages after extraction read the CLEANED tables when the config provides them
# (clean_nostril_planes.py: off-centerline nostril planes dropped, nostril opening added as s = 0)
PE=${PE_IN:-$EXTRACT}; CSAP=${CSA_PLOT:-$CSA}; PLANES_ARG=${PLANES_DIR:+--planes-dir $PLANES_DIR}
# composite layout: columns = renders across the top, each over its plot (bigger geometry); stacked = ENT002/003 style
LAYOUT=${LAYOUT:-columns}; [ "$LAYOUT" = columns ] && REFSUF=_big || REFSUF=
mkdir -p "$POST" "$EXTRACT/logs" "$PANELS"
cd "$POST"
has() { [[ " $STAGES " == *" $1 "* ]]; }
log() { echo "=== [$SUBJ] $* $(date +%T) ==="; }

# derived numbers: frame of peak inspiration (locator), global peak velocity (narrative)
PEAK_INHALE_FRAME=$(python3 -c "
import numpy as np
F=np.array([[float(x) for x in l.split(',')[:2]] for l in open('$PROF').read().splitlines()[1:] if l.strip()])
t=F[np.argmin(F[:,1]),0]; print('out_%06d' % (25*round(1000*t/25)))")
VPK=$([ -f "$VELSCAN" ] && python3 -c "
import numpy as np; D=np.genfromtxt('$VELSCAN',delimiter=',',names=True); print('%.1f' % D['v_max'].max())" || echo nan)
# motion frame of peak dissipated power (snapshot / report), once the WOB plots exist
# (derived again after the plots stage, so a full run from scratch has it for the later stages)
WOBSUM=$PE/plots_wob/$MAIN_PATH/wob_summary.csv
peak_power_frame() {
  [ -f $WOBSUM ] || return 0
  K=$(python3 -c "
import csv,numpy as np
W=list(csv.DictReader(open('$WOBSUM')))
t=float(W[int(np.argmax([float(r['dissipated_power_mW']) for r in W]))]['time_s']); print(int(round(t/0.025)))")
  PEAK_POWER_FRAME=$(printf out_%06d $((K*25)))
}
peak_power_frame
echo "subject $SUBJ  rho=$RHO  peak-inhale frame=$PEAK_INHALE_FRAME  peak-power frame=$PEAK_POWER_FRAME  peak |v|=$VPK m/s  regions: $REGIONS"

CSA_ID=$(basename "${CSA%-Data.csv}")                       # csa-bifurcation --subject-id (ENT005-Data.csv -> ENT005)
NSTEPS=$(grep -A2 'number of steps' "$CASE" | head -1 | awk '{print $NF}')
NW=${WORKERS:-16}

if has csa; then
  log "0a. CSA from the registration: csa-bifurcation $SUBJ_DIR (stop $BREATH_MS ms)"
  [ "$(realpath -m "$(dirname "$CSA")")" = "$(realpath -m "$SUBJ_DIR/csa_bifurcation")" ] || { echo "CSA must be SUBJ_DIR/csa_bifurcation/<id>-Data.csv"; exit 1; }
  conda run --no-capture-output -n ${MIRTK_ENV:-mirtk} csa-bifurcation "$SUBJ_DIR" --subject-id $CSA_ID \
      --start 0 --stop $BREATH_MS --step ${CSA_STEP:-25} > "$SUBJ_DIR/csa_run.log" 2>&1
  tail -n 3 "$SUBJ_DIR/csa_run.log"
fi

if has facecheck; then
  log "0b. face-label guard (csa-bifurcation silent frame-0 fallback)"
  conda run --no-capture-output -n ${MIRTK_ENV:-mirtk} python $HERE/qa/facecheck.py "$SUBJ_DIR" --json $POST/facecheck.json
fi

if has motioncheck; then
  log "0c. CFD mesh-motion check ($NW shards)"
  MC=$POST/motion_check; mkdir -p $MC/logs
  python3 $HERE/qa/reg_traj.py --reg-csv "$REG_CSV" --out $MC/reg_traj.npy ${TABLES_DIR:+--tables-dir $TABLES_DIR}
  for w in $(seq 0 $((NW - 1))); do
    $PVB $HERE/qa/wall_track.py --case "$CASE" --traj $MC/reg_traj.npy --tag $SUBJ --out-dir $MC \
        --motion-dir "$SUBJ_DIR/motion/stl" --shard $w:$NW > $MC/logs/wall_track_$w.log 2>&1 &
  done
  wait
  python3 $HERE/qa/motion_check_fig.py --track-dir $MC --runs "$SUBJ=$SUBJ CFD" --out $MC/${SUBJ}_motion_check.png \
      || echo "WARNING: the CFD wall does not follow the registration - CFD results are not valid (see $MC)"
fi

if has velscan; then
  log "0d. velocity scan -> $VELSCAN"
  VS=$(dirname "$VELSCAN"); mkdir -p $VS; CH=$(( (NSTEPS + 3) / 4 ))
  for w in 0 1 2 3; do
    $PVB $HERE/qa/velocity_peak_scan.py --case "$CASE" --out $VS/part$w.csv --steps $((w * CH)):$(((w + 1) * CH)) > $VS/part$w.log 2>&1 &
  done
  wait
  { head -1 $VS/part0.csv; for f in $VS/part[0-3].csv; do tail -n +2 $f; done | sort -t, -k1,1n; } > "$VELSCAN"
  VPK=$(python3 -c "
import numpy as np; D=np.genfromtxt('$VELSCAN',delimiter=',',names=True); print('%.1f' % D['v_max'].max())")
  python3 -c "
import numpy as np; D=np.genfromtxt('$VELSCAN',delimiter=',',names=True)
print('peak |v| %.1f m/s at t=%.3f s; max over steps of the 99.9th pct %.1f m/s (-> VMAX)' % (D['v_max'].max(), D['t_s'][D['v_max'].argmax()], D['v_p999'].max()))"
fi

if has extract; then
  log "1. CSA-plane pressure / energy extraction"
  rm -f $EXTRACT/logs/_status
  for R in $REGIONS; do
    ( $PVB $HERE/csa_pressure_energy.py --case "$CASE" --csa "$CSA" --planes-dir "$(dirname $CSA)/frames" \
        --region $R --frame-stride 1 --out $EXTRACT/csa_pe_$R.csv > $EXTRACT/logs/$R.log 2>&1
      echo "$R exit=$? $(date +%T)" >> $EXTRACT/logs/_status ) &
  done
  wait; cat $EXTRACT/logs/_status
  grep -q "exit=[^0]" $EXTRACT/logs/_status && { echo "extraction failed, see $EXTRACT/logs"; exit 1; }
fi

if has clean; then
  log "1b. nostril-end cleaning -> $POST/csa_clean"
  CL=$POST/csa_clean; mkdir -p $CL
  python3 $HERE/clean_nostril_planes.py --subject-dir "$SUBJ_DIR" --csa "$CSA" --pe-dir $EXTRACT --out-dir $CL --cfd-steps $NSTEPS
  $PVB $HERE/fill_nostril_planes.py --subject-dir "$SUBJ_DIR" --csa "$CSA" --labels "${CSA%-Data.csv}_full_face_labels.npy" \
      --rejected $CL/rejected_planes.csv --case "$CASE" --out-dir $CL > $CL/fill.log 2>&1
  grep -v -i warn $CL/fill.log | tail -n 3
  python3 $HERE/clean_nostril_planes.py --subject-dir "$SUBJ_DIR" --csa "$CSA" --pe-dir $EXTRACT \
      --fill-csa $CL/fill_csa.csv --fill-pe $CL/fill_pe.csv --out-dir $CL --cfd-steps $NSTEPS
  python3 $HERE/clean_combined_planes.py --csa "$CSA" --frames-dir "$(dirname "$CSA")/frames" \
      --rejected $CL/rejected_planes.csv --out $CL/combined_planes
fi

if has velrender || has velgrid; then
  # velgrid alone: the same sampling pass, written to the manifest without rendering (--dry-run)
  log "2. velocity volume render (vmax $VMAX, max-dim $VEL_MAXDIM)$(has velrender || echo ', sampling grid only')"
  $PVB $HERE/render_velocity_volume.py --case-dir "$(dirname $(dirname $CASE))" --case-file "$CASE" --outdir $VEL \
      --max-dim $VEL_MAXDIM --max-voxels 60e6 --size 1280x960 --preset Turbo --view sagittal --parallel-projection \
      --surface-opacity 0.25 --surface-source part:1 --vmax $VMAX --clamp-above 1.3 $(has velrender || echo --dry-run) \
      > $POST/velocity_render.log 2>&1
fi

if has stlplanes; then
  log "3. moving wall + CSA planes render (dark, every frame)"
  $PVB $HERE/render_stl_planes_seq.py --subject-dir "$SUBJ_DIR" $PLANES_ARG --out $PANELS/stl_planes > $POST/stl_planes_render.log 2>&1
fi

if has plots; then
  log "4. anatomy / pressure / WOB plots"
  for P in $PATHS_ANAT; do
    python3 $HERE/plot_csa_anatomy.py --csa "$CSAP" --prof "$PROF" --path $P --out $POST/csa_anatomy/$P
    python3 $HERE/plot_csa_pe.py --pe-dir $PE --prof "$PROF" --path $P --out $PE/plots/$P
  done
  for P in $WOB_PATHS; do
    python3 $HERE/plot_csa_pe_wob.py --pe-dir $PE --prof "$PROF" --path $P --flow-share 1.0 --out $PE/plots_wob/$P
  done
  peak_power_frame
fi

if has wobgradient; then
  log "4b. standalone WOB gradient plots (old format)"
  python3 $HERE/plot_wob_gradient_all.py --pe-glob "$PE/csa_pe_*.csv" --out $POST/${SUBJ}_WOB_gradient_all.png ${WOB_EXCLUDE:+--exclude $WOB_EXCLUDE}
  python3 $HERE/plot_wob_gradient.py --pe-csv $PE/csa_pe_${MAIN_PATH}.csv --prof "$PROF" --region $MAIN_PATH --rho $RHO \
      --out $POST/${SUBJ}_WOB_gradient.png
fi

if has refs; then
  log "5. centerline locators (frame $PEAK_INHALE_FRAME)"
  for P in $PATHS_ANAT all; do
    python3 $HERE/make_centerline_ref.py --subject-dir "$SUBJ_DIR" --csa "$CSAP" --frame $PEAK_INHALE_FRAME --path $P --out $PANELS
  done
  # large-label versions for the columns layout (the panel title replaces the inner title)
  python3 $HERE/make_centerline_ref.py --subject-dir "$SUBJ_DIR" --csa "$CSAP" --frame $PEAK_INHALE_FRAME --path $MAIN_PATH \
      --out $PANELS --font-scale 1.8 --no-title --suffix _big --prox-label "${PROX_LABEL:-nostril}"
fi

if has video; then
  log "6. breath multipanel video"
  python3 $HERE/link_velocity_frames.py --vel-dir $VEL/frames --stl-dir "$SUBJ_DIR/motion/stl" --out $VEL/frames_by_motion
  python3 $HERE/compose_breath_video.py --csa-csv "$CSAP" --prof "$PROF" --wob-csv $PE --region $MAIN_PATH \
      --planes-dir $PANELS/stl_planes --vel-dir $VEL/frames_by_motion --ref $PANELS/centerline_ref_${MAIN_PATH}${REFSUF}.png \
      --vmax $VMAX --rho $RHO --layout $LAYOUT $([ "$LAYOUT" = columns ] && echo "--dpi 100") --out $PANELS/composite \
      --prox-label "${PROX_LABEL:-nostril}"
  ffmpeg -y -loglevel error -framerate $FPS -i $PANELS/composite/frame_%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 \
      -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
      $POST/${SUBJ}_breath_multipanel.mp4
fi

if has csavideo; then
  log "7. CSA breath video"
  python3 $HERE/compose_csa_video.py --csa-csv "$CSAP" --prof "$PROF" --planes-dir $PANELS/stl_planes \
      --ref $PANELS/centerline_ref_all.png --rho $RHO --regions "$(echo $REGIONS | tr ' ' ',')" --out $PANELS/csa_composite
  ffmpeg -y -loglevel error -framerate $FPS -i $PANELS/csa_composite/frame_%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 \
      $POST/${SUBJ}_CSA_breath.mp4
fi

if has snapshot; then
  log "8. hi-res light snapshot at peak dissipated power"
  I=$((K-1)); [ $I -lt 0 ] && I=0
  echo "    peak-power motion frame k=$K (CFD index $I)"
  $PVB $HERE/render_stl_planes_seq.py --subject-dir "$SUBJ_DIR" $PLANES_ARG --out $PANELS/stl_planes_white_hires \
      --theme light --size 1680x1560 --frames $K:$((K+1))
  # reuse the video render's bounds + grid so the snapshot is framed exactly like the video
  # (pass with --opt=value: the bounds string starts with '-')
  VB=$(python3 -c "import json; print(','.join('%.8f'%v for v in json.load(open('$VEL/manifest.json'))['union_bounds']))")
  VD=$(python3 -c "import json; print(','.join(str(int(v)) for v in json.load(open('$VEL/manifest.json'))['resample_dims']))")
  $PVB $HERE/render_velocity_volume.py --case-dir "$(dirname $(dirname $CASE))" --case-file "$CASE" \
      --outdir $POST/velocity_white_hires --frames $I:$((I+1)) --sampling-bounds="$VB" --dims="$VD" \
      --size 2560x1920 --preset Turbo --view sagittal --parallel-projection --surface-opacity 0.25 \
      --surface-source part:1 --vmax $VMAX --clamp-above 1.3 --background 1,1,1 > $POST/velocity_white_hires.log 2>&1
  mkdir -p $POST/velocity_white_hires/frames_by_motion
  ln -sf ../frames/frame_$(printf %04d $I).png $POST/velocity_white_hires/frames_by_motion/frame_$(printf %04d $K).png
  python3 $HERE/compose_breath_video.py --csa-csv "$CSAP" --prof "$PROF" --wob-csv $PE --region $MAIN_PATH \
      --planes-dir $PANELS/stl_planes_white_hires --vel-dir $POST/velocity_white_hires/frames_by_motion \
      --ref $PANELS/centerline_ref_${MAIN_PATH}_light${REFSUF}.png --vmax $VMAX --rho $RHO --theme light \
      --layout $LAYOUT $([ "$LAYOUT" = columns ] && echo "--figsize 14.5,8.2 --dpi 300" || echo "--dpi 200") \
      --frames $K:$((K+1)) --out $PANELS/composite_light
  cp $PANELS/composite_light/frame_$(printf %04d $K).png $PANELS/snapshot_light_hires.png
fi

if has regions; then
  log "9a. work of breathing by region (inhale vs exhale)"
  python3 $HERE/wob_regions.py --pe-dir $PE --prof "$PROF" --rho $RHO --main-path $MAIN_PATH --other-path ${OTHER_PATH:-RightNose} \
      --regions "$WOB_REGIONS" --subject $SUBJ --out-dir $POST/wob_regions \
      --bin-mm ${WOB_BIN_MM:-25} --peak-time $(python3 -c "print($K*0.025)") --series
fi

if has pdf; then
  # one page per instant, same layout: peak dissipated power (WOB), peak inhalation, peak exhalation
  K_IN=$(python3 -c "
import numpy as np
F=np.loadtxt('$PROF',delimiter=',',skiprows=1); print(int(round(F[np.argmin(F[:,1]),0]/0.025)))")
  K_EX=$(python3 -c "
import numpy as np
F=np.loadtxt('$PROF',delimiter=',',skiprows=1); print(int(round(F[np.argmax(F[:,1]),0]/0.025)))")
  log "9. report pages (no text): peak WOB k=$K, peak inhalation k=$K_IN, peak exhalation k=$K_EX"
  RP=$PANELS/report; mkdir -p $RP $POST/velocity_white_hires/frames_by_motion
  LD=$POST/wob_regions/loss_per_mm
  # WOB_CMAX (env) fixes the centerline colour scale, e.g. to share one scale between subjects being compared
  CMAX=${WOB_CMAX:-$(python3 -c "
import glob,numpy as np
print('%.3f' % max(np.loadtxt(f,delimiter=',',skiprows=1)[:,1].max() for f in glob.glob('$LD/t*.csv')))")}
  VB=$(python3 -c "import json; print(','.join('%.8f'%v for v in json.load(open('$VEL/manifest.json'))['union_bounds']))")
  VD=$(python3 -c "import json; print(','.join(str(int(v)) for v in json.load(open('$VEL/manifest.json'))['resample_dims']))")
  PAGES=
  for KK in $K $K_IN $K_EX; do
    FR=$(printf out_%06d $((KK*25))); II=$((KK-1)); [ $II -lt 0 ] && II=0
    case $KK in $K) TITLE="Peak work of breathing (WOB)";; $K_IN) TITLE="Peak inhalation";; *) TITLE="Peak exhalation";; esac
    if [ ! -f $POST/velocity_white_hires/frames/frame_$(printf %04d $II).png ]; then
      # one retry: the GPU volume render can segfault transiently (ENT005 reproduction run, step 59)
      for TRY in 1 2; do
        $PVB $HERE/render_velocity_volume.py --case-dir "$(dirname $(dirname $CASE))" --case-file "$CASE" \
            --outdir $POST/velocity_white_hires --frames $II:$((II+1)) --sampling-bounds="$VB" --dims="$VD" \
            --size 2560x1920 --preset Turbo --view sagittal --parallel-projection --surface-opacity 0.25 \
            --surface-source part:1 --vmax $VMAX --clamp-above 1.3 --background 1,1,1 > $POST/velocity_white_hires_$II.log 2>&1 && break
        [ $TRY = 2 ] && { echo "ERROR: hires velocity render of step $II failed twice (see $POST/velocity_white_hires_$II.log)"; exit 1; }
        echo "hires velocity render of step $II crashed, retrying once"; mv $POST/velocity_white_hires_$II.log $POST/velocity_white_hires_$II.try1.log
      done
    fi
    ln -sf ../frames/frame_$(printf %04d $II).png $POST/velocity_white_hires/frames_by_motion/frame_$(printf %04d $KK).png
    python3 $HERE/make_centerline_ref.py --subject-dir "$SUBJ_DIR" --csa "$CSAP" --frame $FR --path $MAIN_PATH \
        --out $RP --themes light --font-scale 1.5 --no-title --no-colorbar --suffix _report_k$KK \
        --mark-step ${WOB_BIN_MM:-25} --prox-label "${PROX_LABEL:-nostril}"
    python3 $HERE/make_centerline_ref.py --subject-dir "$SUBJ_DIR" --csa "$CSAP" --frame $FR --path $MAIN_PATH \
        --out $RP --themes light --font-scale 1.5 --no-title --suffix _report_k$KK --tag wob \
        --prox-label "${PROX_LABEL:-nostril}" --color-label "WOB rate (mW/mm)" --color-max $CMAX \
        --color-csv $LD/t$(python3 -c "print('%.3f' % (max($KK, 1) * 0.025))").csv
    python3 $HERE/compose_breath_video.py --csa-csv "$CSAP" --prof "$PROF" --wob-csv $PE --region $MAIN_PATH \
        --planes-dir $PANELS/stl_planes_white_hires --vel-dir $POST/velocity_white_hires/frames_by_motion \
        --ref $RP/centerline_ref_${MAIN_PATH}_light_report_k$KK.png --wob-ref $RP/centerline_wob_${MAIN_PATH}_light_report_k$KK.png \
        --bins-csv $POST/wob_regions/${SUBJ}_wob_by_${WOB_BIN_MM:-25}mm.csv --vmax $VMAX --rho $RHO --theme light \
        --layout report --dpi 300 --font-scale 0.68 --frames $KK:$((KK+1)) --out $RP --pdf-out $RP/page_k$KK.pdf \
        --page-title "$TITLE" \
        --mark-step ${WOB_BIN_MM:-25} --prox-label "${PROX_LABEL:-nostril}"
    PAGES="$PAGES $RP/page_k$KK.pdf"; PNGS="$PNGS $RP/frame_$(printf %04d $KK).png"
  done
  pdfunite $PAGES $POST/${SUBJ}_CFD_report.pdf
  # Word version: the same page images (npm package 'docx'; NODE_PATH may point to it)
  node $HERE/report_to_docx.js --title "$SUBJ CFD report" $POST/${SUBJ}_CFD_report.docx $PNGS \
      || echo "(Word version skipped: npm package 'docx' not found)"
  echo "wrote $POST/${SUBJ}_CFD_report.pdf (pages: k =$( echo $PAGES | grep -o 'k[0-9]*' | tr '\n' ' '))"
fi

if has onepager; then
  log "9b. text one-pager (old format)"
  python3 $HERE/build_onepager.py --subject $SUBJ --path $MAIN_PATH \
      --wob $PE/plots_wob/$MAIN_PATH/wob_summary.csv --anatomy $POST/csa_anatomy/$MAIN_PATH/csa_anatomy_summary.csv \
      --pe-dir $PE --prof "$PROF" --csa "$CSAP" --snapshot $PANELS/snapshot_light_hires.png \
      --peak-velocity $VPK ${WOB_REGIONS:+--regions-csv $POST/wob_regions/${SUBJ}_wob_by_region.csv --regions-png $POST/wob_regions/${SUBJ}_wob_by_region.png} \
      --out $POST/${SUBJ}_CFD_report_text.pdf
fi
if has vellight; then
  log "10a. velocity render on white, all steps (report video)"
  VL=$POST/velocity_breath_light; mkdir -p $VL/logs
  VB=$(python3 -c "import json; print(','.join('%.8f'%v for v in json.load(open('$VEL/manifest.json'))['union_bounds']))")
  VD=$(python3 -c "import json; print(','.join(str(int(v)) for v in json.load(open('$VEL/manifest.json'))['resample_dims']))")
  VLARGS=(--case-dir "$(dirname $(dirname $CASE))" --case-file "$CASE" --outdir $VL \
      --sampling-bounds="$VB" --dims="$VD" --size 1600x1200 --preset Turbo --view sagittal --parallel-projection \
      --surface-opacity 0.25 --surface-source part:1 --vmax $VMAX --clamp-above 1.3 --background 1,1,1)
  for w in 0 1 2 3 4 5; do
    $PVB $HERE/render_velocity_volume.py "${VLARGS[@]}" --shard $w:6 > $VL/logs/shard$w.log 2>&1 &
  done
  wait   # (returns 0 even when a shard crashed -> check every step below)
  # a shard that hit a transient GPU segfault leaves its remaining steps missing: render those once more, one by one
  ND=$(( ${#NSTEPS} > 4 ? ${#NSTEPS} : 4 ))
  for ((i = 0; i < NSTEPS; i++)); do
    [ -f $VL/frames/$(printf "frame_%0${ND}d.png" $i) ] && continue
    echo "velocity step $i missing after the shard pass, rendering it again"
    $PVB $HERE/render_velocity_volume.py "${VLARGS[@]}" --frames $i:$((i+1)) > $VL/logs/retry_$i.log 2>&1 \
        || { echo "ERROR: velocity render of step $i failed twice (see $VL/logs/retry_$i.log)"; exit 1; }
  done
  python3 $HERE/link_velocity_frames.py --vel-dir $VL/frames --stl-dir "$SUBJ_DIR/motion/stl" --out $VL/frames_by_motion
fi

if has reportvideo; then
  log "10. report video (the PDF page over the breath)"
  RV=$PANELS/report_video; LD=$POST/wob_regions/loss_per_mm
  # WOB_CMAX (env) fixes the centerline colour scale, e.g. to share one scale between subjects being compared
  CMAX=${WOB_CMAX:-$(python3 -c "
import glob,numpy as np
print('%.3f' % max(np.loadtxt(f,delimiter=',',skiprows=1)[:,1].max() for f in glob.glob('$LD/t*.csv')))")}
  python3 $HERE/make_centerline_ref.py --subject-dir "$SUBJ_DIR" --csa "$CSAP" --path $MAIN_PATH --out $RV/wob \
      --themes light --font-scale 1.5 --no-title --tag wob --prox-label "${PROX_LABEL:-nostril}" --color-label "WOB rate (mW/mm)" \
      --color-series $LD --color-max $CMAX --limits-from $PEAK_POWER_FRAME --no-tight --dpi 200
  python3 $HERE/compose_breath_video.py --csa-csv "$CSAP" --prof "$PROF" --wob-csv $PE --region $MAIN_PATH \
      --planes-dir $PANELS/stl_planes --vel-dir $POST/velocity_breath_light/frames_by_motion \
      --ref $PANELS/report/centerline_ref_${MAIN_PATH}_light_report_k$K.png --wob-ref-dir $RV/wob \
      --bins-csv $POST/wob_regions/${SUBJ}_wob_by_${WOB_BIN_MM:-25}mm.csv \
      --bins-power-csv $POST/wob_regions/${SUBJ}_power_by_${WOB_BIN_MM:-25}mm_vs_time.csv \
      --vmax $VMAX --rho $RHO --theme light --layout report --dpi 128 --font-scale 0.68 \
      --mark-step ${WOB_BIN_MM:-25} --prox-label "${PROX_LABEL:-nostril}" --page-title "Breathing cycle" --out $RV/frames
  ffmpeg -y -loglevel error -framerate $FPS -i $RV/frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 \
      -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" $POST/${SUBJ}_report_breath.mp4
fi
log "done: $STAGES"
