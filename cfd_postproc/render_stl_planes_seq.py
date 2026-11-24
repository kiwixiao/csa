#!/usr/bin/env pvbatch
"""Render the airway wall + CSA planes panel for EVERY motion frame.

Ported from /storage4/ENT002/render_stl_white.py, which rendered a single frame
for the one-pager. Differences:
  * iterates the whole motion/stl sequence instead of one index k
  * --subject-dir / --out are arguments
  * frame index comes from the sorted file list, so it does not assume a fixed
    25*(k+1) naming (ENT003 starts at out_000000 and ends on an appended
    remainder frame at out_002180)
  * dark background to match the video theme (ENT002's white variant was for PDF)

Usage:
  pvbatch render_stl_planes_seq.py --subject-dir /storage4/ENT003/ENT003 \
      --out /storage4/ENT003/postproc/panels/stl_planes
"""
import argparse, glob, math, os, sys
from paraview.simple import *  # noqa

p = argparse.ArgumentParser()
p.add_argument("--subject-dir", required=True)
p.add_argument("--planes-dir", default=None, help="folder of <frame>-Planes-All.stl (default <subject>/csa_bifurcation/combined_planes)")
p.add_argument("--out", required=True)
p.add_argument("--size", default="820x760")
p.add_argument("--theme", default="dark", choices=["dark", "light"])
p.add_argument("--frames", default=None,
               help="index range a:b (0-based, b exclusive). Output names keep the absolute "
                    "index, so a single-frame hi-res render for the PDF still matches the "
                    "video's frame numbering.")
a = p.parse_args()

SZ = [int(x) for x in a.size.split("x")]
os.makedirs(a.out, exist_ok=True)
stl_dir = os.path.join(a.subject_dir, "motion", "stl")
plane_dir = a.planes_dir or os.path.join(a.subject_dir, "csa_bifurcation", "combined_planes")
frames = sorted(glob.glob(os.path.join(stl_dir, "out_*.stl")))
if not frames:
    print("ERROR: no frames in %s" % stl_dir); sys.exit(1)
print("frames:", len(frames), flush=True)

paraview.simple._DisableFirstRenderCameraReset()
view = CreateView("RenderView"); view.ViewSize = SZ
try:
    view.UseColorPaletteForBackground = 0
except Exception:
    pass
view.Background = [1.0, 1.0, 1.0] if a.theme == "light" else [0.059, 0.059, 0.078]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1

# fixed camera from the first frame so framing is identical across the sequence
ref = OpenDataFile(frames[0])
Show(ref, view)
b = ref.GetDataInformation().GetBounds()
center = [(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, (b[4] + b[5]) / 2]
diag = math.sqrt(sum((b[2 * i + 1] - b[2 * i]) ** 2 for i in range(3)))
view.CameraFocalPoint = center
view.CameraPosition = [center[0] - diag * 2.0, center[1], center[2]]
view.CameraViewUp = [0, 0, 1]
ResetCamera(view)
CAM = dict(fp=list(view.CameraFocalPoint), pos=list(view.CameraPosition),
           up=list(view.CameraViewUp), scale=view.CameraParallelScale)
Hide(ref, view); Delete(ref); del ref

wall_rgb = [0.42, 0.47, 0.55] if a.theme == "light" else [0.55, 0.62, 0.72]
wall_op = 0.30 if a.theme == "light" else 0.22

KLO, KHI = 0, len(frames)
if a.frames:
    _a, _b = a.frames.split(":")
    KLO, KHI = (int(_a) if _a else 0), (int(_b) if _b else len(frames))
    print("rendering frames %d:%d only" % (KLO, KHI), flush=True)

for k in range(KLO, KHI):
    fpath = frames[k]
    stem = os.path.splitext(os.path.basename(fpath))[0]        # out_000025
    wall = OpenDataFile(fpath)
    dw = Show(wall, view)
    dw.DiffuseColor = wall_rgb; dw.Opacity = wall_op
    dw.Specular = 0.15; dw.Ambient = 0.35

    ppath = os.path.join(plane_dir, "%s-Planes-All.stl" % stem)
    dp = pl = None
    if os.path.exists(ppath):
        pl = OpenDataFile(ppath)
        dp = Show(pl, view)
        dp.DiffuseColor = [0.90, 0.42, 0.06]; dp.Opacity = 0.85

    view.CameraFocalPoint = CAM["fp"]; view.CameraPosition = CAM["pos"]
    view.CameraViewUp = CAM["up"]; view.CameraParallelScale = CAM["scale"]
    Render(view)
    SaveScreenshot("%s/frame_%04d.png" % (a.out, k), view, ImageResolution=SZ)

    Hide(wall, view); Delete(wall); del wall
    if pl is not None:
        Hide(pl, view); Delete(pl); del pl
    if k % 10 == 0:
        print("  frame %d/%d (%s)" % (k, len(frames), stem), flush=True)

print("DONE -> %s (%d frames)" % (a.out, KHI - KLO), flush=True)
