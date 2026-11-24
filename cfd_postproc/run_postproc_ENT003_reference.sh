#!/bin/bash
# ENT003 post-processing, rebuilt on the corrected 2026-08-20 EnSight export
# (87 steps, t=0.025..2.175 s, both MaskInlet and BronchiOut caps present).
#
# Prerequisites, both launched separately because they take ~1-2 h:
#   run_csa_pe_all.sh      -> postproc/csa_pe/csa_pe_<Region>.csv   (4 regions)
#   run_velocity_render.sh -> postproc/velocity_breath_final/frames (87 CFD frames)
#
# Everything below is cheap and re-runnable. rho = 1.2367 kg/m3 is the density the
# solver actually ran; see postproc/outlet_flow/inlet_outlet_massflow.csv.
set -e
cd /storage4/ENT003

RHO=1.2367
# smoothed CSA table, used for BOTH solvers so the CSA panels are identical geometry.
# RAW stays the canonical plane-STL ordering (smoothing drops planes, which shifts ranks).
CSA=/storage4/OF/OF_ENT003_LINUX/ENT003-Data-OFsmooth.csv
CSA_RAW=ENT003/csa_bifurcation/ENT003-Data.csv
PROF=ENT003FlowProfile.csv
PE=postproc/csa_pe
PANELS=postproc/panels
PVB=/tmp/claude-1000/-storage4-ENT003/ba64798d-019f-46f3-9853-bf4bf522e353/scratchpad/pvb.sh
MAIN_PATH=${MAIN_PATH:-LeftNose}          # airway route featured in the video + PDF
FPS=${FPS:-12}

echo "=== 1. velocity frame index bridge (CFD index -> motion index) ==="
python3 link_velocity_frames.py \
    --vel-dir postproc/velocity_breath_final/frames \
    --stl-dir ENT003/motion/stl \
    --out postproc/velocity_breath_final/frames_by_motion

echo "=== 2. anatomy + pressure/energy plots per path ==="
# Anatomy + pressure for every branch, but work-of-breathing ONLY for the routes that
# run all the way to the trachea. The Mouth branch stops at the oropharynx and carries
# ~1/3 of the flow (4.7 of 14.9 L/min, area-integrated at peak inhale), so its dP*Q is a
# partial-route number and is deliberately not produced.
for P in LeftNose RightNose Mouth; do
    python3 plot_csa_anatomy.py  --csa $CSA --prof $PROF --path $P --out postproc/csa_anatomy/$P
    python3 plot_csa_pe.py       --pe-dir $PE --prof $PROF --path $P --out $PE/plots/$P
done
for P in LeftNose RightNose; do
    # flow-share 1.0: the nasal and oral routes are parallel paths between the same two
    # ends (mask and trachea), so the measured opening->trachea drop is common to them
    # and total dissipation is that drop x the FULL prescribed flow.
    python3 plot_csa_pe_wob.py   --pe-dir $PE --prof $PROF --path $P --flow-share 1.0 \
                                 --out $PE/plots_wob/$P
done

echo "=== 3. WOB gradient plots (all regions overlaid) ==="
python3 plot_wob_gradient_all.py --pe-glob "$PE/csa_pe_*.csv" \
    --out postproc/ENT003_WOB_gradient_all.png --independent Mouth --exclude Mouth
python3 plot_wob_gradient.py --pe-csv $PE/csa_pe_${MAIN_PATH}.csv --prof $PROF \
    --region $MAIN_PATH --rho $RHO --out postproc/ENT003_WOB_gradient.png

echo "=== 4. anatomical locators ==="
for P in LeftNose RightNose Mouth all; do
    python3 make_centerline_ref.py --subject-dir ENT003 --csa $CSA \
        --frame out_001950 --path $P --out $PANELS
done

echo "=== 5. breath multipanel video (dark) ==="
python3 compose_breath_video.py \
    --csa-csv $CSA --prof $PROF --wob-csv $PE --region $MAIN_PATH \
    --planes-dir $PANELS/stl_planes \
    --vel-dir postproc/velocity_breath_final/frames_by_motion \
    --ref $PANELS/centerline_ref_${MAIN_PATH}.png \
    --vmax 16 --rho $RHO --out $PANELS/composite
ffmpeg -y -loglevel error -framerate $FPS -i $PANELS/composite/frame_%04d.png \
    -c:v libx264 -pix_fmt yuv420p -crf 18 postproc/ENT003_breath_multipanel.mp4

echo "=== 6. CSA breath video (all regions) ==="
python3 compose_csa_video.py \
    --csa-csv $CSA --prof $PROF --planes-dir $PANELS/stl_planes \
    --ref $PANELS/centerline_ref_all.png --rho $RHO \
    --out $PANELS/csa_composite
ffmpeg -y -loglevel error -framerate $FPS -i $PANELS/csa_composite/frame_%04d.png \
    -c:v libx264 -pix_fmt yuv420p -crf 18 postproc/ENT003_CSA_breath.mp4

echo "=== 7. hi-res light snapshot at peak dissipated power ==="
K=$(python3 -c "
import csv,numpy as np
W=list(csv.DictReader(open('$PE/plots_wob/$MAIN_PATH/wob_summary.csv')))
t=float(W[int(np.argmax([float(r['dissipated_power_mW']) for r in W]))]['time_s'])
print(int(round(t/0.025)))")
I=$((K-1)); [ $I -lt 0 ] && I=0
echo "    peak-power motion frame k=$K  (CFD index $I)"

$PVB render_stl_planes_seq.py --subject-dir ENT003 --out $PANELS/stl_planes_white_hires \
    --theme light --size 1680x1560 --frames $K:$((K+1))
# reuse the video render's union bounds + grid so the PDF snapshot is framed
# identically to the video (a 1-frame run would otherwise sample its own bounds)
# NOTE: pass these with --opt=value, not --opt value. The bounds string starts with a
# minus sign, and argparse treats a leading-dash value as the next option ("expected one
# argument") unless it is attached with '='.
VB=$(python3 -c "
import json; print(','.join('%.8f'%v for v in json.load(open('postproc/velocity_breath_final/manifest.json'))['union_bounds']))")
VD=$(python3 -c "
import json; print(','.join(str(int(v)) for v in json.load(open('postproc/velocity_breath_final/manifest.json'))['resample_dims']))")
[ -n "$VB" ] && [ -n "$VD" ] || { echo "ERROR: could not read render bounds from manifest"; exit 1; }
echo "    reusing video framing: bounds=$VB dims=$VD"
$PVB render_velocity_volume.py --case-dir /storage4/ENT003 --case-file case/ENT003.case \
    --outdir postproc/velocity_white_hires --frames $I:$((I+1)) \
    --sampling-bounds="$VB" --dims="$VD" \
    --size 2560x1920 --preset Turbo --view sagittal --parallel-projection \
    --surface-opacity 0.25 --surface-source part:1 --vmax 16 --clamp-above 1.3 \
    --background 1,1,1
mkdir -p postproc/velocity_white_hires/frames_by_motion
ln -sf ../frames/frame_$(printf %04d $I).png \
    postproc/velocity_white_hires/frames_by_motion/frame_$(printf %04d $K).png

python3 compose_breath_video.py \
    --csa-csv $CSA --prof $PROF --wob-csv $PE --region $MAIN_PATH \
    --planes-dir $PANELS/stl_planes_white_hires \
    --vel-dir postproc/velocity_white_hires/frames_by_motion \
    --ref $PANELS/centerline_ref_${MAIN_PATH}_light.png \
    --vmax 16 --rho $RHO --theme light --dpi 200 --frames $K:$((K+1)) \
    --out $PANELS/composite_light
cp $PANELS/composite_light/frame_$(printf %04d $K).png $PANELS/snapshot_light_hires.png

echo "=== 8. clinical one-pager ==="
VPK=19.9      # measured global max |v| over all 87 steps (at t=1.95 s, peak inspiration)
python3 build_onepager.py --subject ENT003 --path $MAIN_PATH \
    --wob $PE/plots_wob/$MAIN_PATH/wob_summary.csv \
    --anatomy postproc/csa_anatomy/$MAIN_PATH/csa_anatomy_summary.csv \
    --pe-dir $PE --prof $PROF \
    --csa $CSA --snapshot $PANELS/snapshot_light_hires.png \
    --peak-velocity $VPK \
    --out postproc/ENT003_CFD_onepager.pdf

echo "=== ALL DONE $(date +%T) ==="
