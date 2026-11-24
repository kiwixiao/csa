#!/usr/bin/env python3
"""Guard against csa-bifurcation's SILENT frame-0 fallback.

slice_bifurcation.run_single_frame slices the deformed frame with the t0 face labels only if the label count equals
the frame's face count as trimesh loads it; otherwise it quietly slices the t0 region meshes for every frame (CSA
frozen at t0, no warning). Cause seen on ENT005: sub-micron sliver triangles in the t0 surface turn degenerate after
the FFD and trimesh drops them. Fix: geometry/drop_slivers.py on the t0 surface, re-run csa-bifurcation.

Run in the env csa-bifurcation runs in (same trimesh loader):
    conda run -n mirtk python qa/facecheck.py <subject_dir> [<subject_dir> ...] [--json out.json]
Exit status 1 if any frame of any subject mismatches.
"""
import argparse, glob, json, os, sys
import numpy as np, trimesh

ap = argparse.ArgumentParser(); ap.add_argument("subject_dirs", nargs="+"); ap.add_argument("--json", default=None)
a = ap.parse_args()
out, ok = {}, True
for sd in a.subject_dirs:
    lab = glob.glob(os.path.join(sd, "csa_bifurcation", "*_full_face_labels.npy"))
    if len(lab) != 1:
        print("%s: expected one csa_bifurcation/*_full_face_labels.npy, found %d" % (sd, len(lab))); ok = False; continue
    n = len(np.load(lab[0]))
    frames = sorted(glob.glob(os.path.join(sd, "motion", "stl", "out_*.stl")))
    bad = [os.path.basename(f) for f in frames if len(trimesh.load_mesh(f).faces) != n]
    out[sd] = dict(labels=n, frames=len(frames), frames_mismatch=len(bad), first_bad=bad[:5])
    ok &= not bad and len(frames) > 0
    print("%s: %d labels, %d frames, %d mismatched %s -> %s" % (sd, n, len(frames), len(bad), bad[:5], "OK" if not bad and frames else "FAIL"), flush=True)
if a.json:
    json.dump(out, open(a.json, "w"), indent=1)
sys.exit(0 if ok else 1)
