#!/usr/bin/env python3
"""Bridge the velocity render's CFD frame numbering to the motion/CSA frame numbering.

render_velocity_volume.py names its output by CFD timestep index (0..86, t=0.025..2.175 s).
compose_breath_video.py indexes by motion frame (0..88, out_000000..out_002180). The two
differ because the motion sequence carries a leading out_000000 and a trailing remainder
frame out_002180, both of which land on an existing CFD solve:

    cfd_idx = clamp(round(frame_ms / 25) - 1, 0, n_cfd - 1)

    out_000000 -> cfd 0     out_000025 -> cfd 0     out_000050 -> cfd 1   ...
    out_002175 -> cfd 86    out_002180 -> cfd 86

This writes a directory of symlinks named by motion index so the composer can stay
index-by-position. The old ENT003 run fed the composer the raw CFD numbering, which was
off by one for every frame after the first.

Usage:
  python3 link_velocity_frames.py --vel-dir postproc/velocity_breath_final/frames \
      --stl-dir ENT003/motion/stl --out postproc/velocity_breath_final/frames_by_motion
"""
import argparse
import glob
import os
import re

p = argparse.ArgumentParser()
p.add_argument("--vel-dir", required=True, help="dir of frame_XXXX.png named by CFD index")
p.add_argument("--stl-dir", required=True, help="motion/stl dir, to recover the frame ms list")
p.add_argument("--out", required=True)
p.add_argument("--dt-ms", type=float, default=25.0)
p.add_argument("--offset", type=int, default=-1,
               help="cfd_idx = round(ms/dt) + offset. -1 for the STAR export, whose first "
                    "step is t=0.025 s; 0 for the OpenFOAM export, whose first step is "
                    "t=0.0 s.")
a = p.parse_args()

os.makedirs(a.out, exist_ok=True)
for old in glob.glob(os.path.join(a.out, "frame_*.png")):
    os.unlink(old)

src = sorted(glob.glob(os.path.join(a.vel_dir, "frame_*.png")))
if not src:
    raise SystemExit("no velocity frames in %s" % a.vel_dir)
have = {int(re.search(r"frame_(\d+)", os.path.basename(f)).group(1)): f for f in src}
ncfd = max(have) + 1

stems = sorted(glob.glob(os.path.join(a.stl_dir, "out_*.stl")))
if not stems:
    raise SystemExit("no motion frames in %s" % a.stl_dir)

made, missing = 0, []
for k, sp in enumerate(stems):
    ms = int(re.search(r"out_(\d+)", os.path.basename(sp)).group(1))
    idx = int(round(ms / a.dt_ms)) + a.offset
    idx = max(0, min(ncfd - 1, idx))
    if idx not in have:
        missing.append((k, ms, idx))
        continue
    os.symlink(os.path.relpath(have[idx], a.out), os.path.join(a.out, "frame_%04d.png" % k))
    made += 1

print("linked %d/%d motion frames -> %d CFD velocity frames in %s"
      % (made, len(stems), len(have), a.out))
if missing:
    print("MISSING velocity frames for motion indices:",
          ", ".join("%d(ms=%d->cfd%d)" % m for m in missing[:20]),
          "..." if len(missing) > 20 else "")
