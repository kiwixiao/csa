"""CFD mesh-motion check, one shard of the saved steps (pvbatch).

Per saved CFD step: (1) distance of every registration control point (its position at that ms, qa/reg_traj.py)
to the CFD wall (Wall + Nose blocks, mm) -- if the mesh follows the prescribed motion these sit on the wall;
(2) signed distance of every 4th CFD wall point to the registration surface out_<ms>.stl (> 0 outside), if
--motion-dir is given. Writes <out-dir>/track_<tag>_<k>.npz.
Usage: pvb.sh wall_track.py --case X.case --traj reg_traj.npy --tag NAME --out-dir DIR [--motion-dir SUBJ/motion/stl] [--shard w:N]
"""
import argparse, os, numpy as np, vtk
from paraview.simple import OpenDataFile
from paraview import servermanager as sm
from vtkmodules.util.numpy_support import vtk_to_numpy

ap = argparse.ArgumentParser(); ap.add_argument("--case", required=True); ap.add_argument("--traj", required=True)
ap.add_argument("--tag", required=True); ap.add_argument("--out-dir", required=True); ap.add_argument("--motion-dir", default=None)
ap.add_argument("--wall-blocks", default="Wall,Nose"); ap.add_argument("--shard", default="0:1")
a = ap.parse_args()
w, nsh = (int(x) for x in a.shard.split(":")); want = a.wall_blocks.split(",")
CT = np.load(a.traj)                                     # (points, ms, 3) mm
r = OpenDataFile(a.case); r.PointArrays = []; tv = list(r.TimestepValues)


def wall(t):
    r.UpdatePipeline(t); d = sm.Fetch(r); app = vtk.vtkAppendPolyData(); it = d.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(d.NAME()).strip().split(".")[-1] in want:
            g = vtk.vtkGeometryFilter(); g.SetInputData(it.GetCurrentDataObject()); g.Update(); app.AddInputData(g.GetOutput())
        it.GoToNextItem()
    app.Update(); tr = vtk.vtkTransform(); tr.Scale(1000, 1000, 1000)
    tf = vtk.vtkTransformPolyDataFilter(); tf.SetInputData(app.GetOutput()); tf.SetTransform(tr); tf.Update()
    tri = vtk.vtkTriangleFilter(); tri.SetInputData(tf.GetOutput()); tri.Update()
    o = vtk.vtkPolyData(); o.DeepCopy(tri.GetOutput()); return o


def stl(p):
    rd = vtk.vtkSTLReader(); rd.SetFileName(p); rd.Update()
    nm = vtk.vtkPolyDataNormals(); nm.SetInputData(rd.GetOutput()); nm.ConsistencyOn(); nm.AutoOrientNormalsOn(); nm.SplittingOff(); nm.Update()
    return nm.GetOutput()


def closest(poly, P):
    loc = vtk.vtkStaticCellLocator(); loc.SetDataSet(poly); loc.BuildLocator()
    cp = [0., 0., 0.]; ci = vtk.reference(0); si = vtk.reference(0); d2 = vtk.reference(0.)
    out = np.empty(len(P))
    for i, p in enumerate(P):
        loc.FindClosestPoint(p, cp, ci, si, d2); out[i] = np.sqrt(d2.get())
    return out


os.makedirs(a.out_dir, exist_ok=True)
for k in range(w, len(tv), nsh):
    ms = int(round(tv[k] * 1000))
    if ms >= CT.shape[1]:
        print("step %d t=%g s is beyond the trajectories (%d ms), skipped" % (k, tv[k], CT.shape[1])); continue
    W = wall(tv[k]); P = vtk_to_numpy(W.GetPoints().GetData()).astype(float)
    dc = closest(W, CT[:, ms].astype(float))                       # control points -> CFD wall
    res = dict(ms=ms, dc=dc, ctrl=CT[:, ms])
    msg = ""
    if a.motion_dir:
        f = os.path.join(a.motion_dir, "out_%06d.stl" % ms)
        if os.path.exists(f):
            ipd = vtk.vtkImplicitPolyDataDistance(); ipd.SetInput(stl(f))
            sub = P[::4]; sd = np.array([ipd.EvaluateFunction(p) for p in sub]); res.update(pts=sub, sd=sd)
            msg = " | wall->motion sd p5 %.2f med %.2f p95 %.2f" % (np.percentile(sd, 5), np.median(sd), np.percentile(sd, 95))
    np.savez_compressed(os.path.join(a.out_dir, "track_%s_%03d.npz" % (a.tag, k)), **res)
    print(a.tag, k, ms, "ctrl->wall med %.3f p95 %.3f max %.3f%s" % (np.median(dc), np.percentile(dc, 95), dc.max(), msg), flush=True)
