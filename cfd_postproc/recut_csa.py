#!/usr/bin/env pvbatch
# -*- coding: utf-8 -*-
"""Re-cut CSA planes on the CFD WALL of a (virtually edited) case, for a span of one region.

The CSA table (<S>-Data.csv) was cut from the baseline surface. For an edited geometry run on the
same breath (ENT005 'Aryt': larynx widened), the planes (origin = baseline centroid, normal = the
slicer's disk PCA normal, as csa_pressure_energy.py) are kept so distance along the centerline stays
the baseline's, and only the AREA is re-measured: the plane cuts the exported ENT002.Wall surface of
the matched CFD step; closed outlines whose centre lies within the slicer disk's reach (+ margin) are
summed, an outline nested inside another (the epiglottis inside the laryngeal inlet) subtracted. Run it on the baseline case too: the ratio
edited / baseline re-cut carries the edit onto the slicer's own areas (method differences cancel).

Usage: pvb.sh recut_csa.py --case X.case --csa <S>-Data.csv --region DescendingAirway
           --arc-min 80 --arc-max 130 --out recut.csv [--shard w:N]
"""
import builtins as _bi
import argparse, os, re, sys, time
import numpy as np
from paraview.simple import OpenDataFile
from paraview import servermanager as sm
import vtk
from vtkmodules.util.numpy_support import vtk_to_numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from csa_pressure_energy import REGION_TAG, CFD_DT, stl_normal, rank_map, read_planes  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--case", required=True); ap.add_argument("--csa", required=True, help="RAW <S>-Data.csv (plane-STL ranks)")
ap.add_argument("--region", default="DescendingAirway", choices=sorted(REGION_TAG))
ap.add_argument("--arc-min", type=float, required=True); ap.add_argument("--arc-max", type=float, required=True)
ap.add_argument("--wall-block", default="Wall,Nose", help="export block(s) holding the airway wall, comma-separated")
ap.add_argument("--lobe-radius-factor", type=float, default=2.0); ap.add_argument("--lobe-radius-min", type=float, default=6.0)
ap.add_argument("--reach-margin", type=float, default=5.0,
                help="mm beyond the slicer disk's reach in which outlines count (the edit widens the lumen)")
ap.add_argument("--shard", default=None, help="w:N - every Nth frame starting at w")
ap.add_argument("--out", required=True)
a = ap.parse_args()

planes_dir = os.path.join(os.path.dirname(a.csa), "frames")
tag = REGION_TAG[a.region]
frames = read_planes(a.csa, a.region)
RANK = rank_map(a.csa, a.region)
fnames = sorted(frames, key=lambda s: int(re.search(r"(\d+)", s).group(1)))
if a.shard:
    w, nsh = (int(x) for x in a.shard.split(":")); fnames = fnames[w::nsh]

r = OpenDataFile(a.case); r.PointArrays = []
tv = list(r.TimestepValues)
want = [s.strip() for s in a.wall_block.split(",")]


def wall_poly(t):
    """The wall block(s) at time t as vtkPolyData in mm."""
    r.UpdatePipeline(t)
    data = sm.Fetch(r)
    app = vtk.vtkAppendPolyData()
    it = data.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        nm = it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()).strip().split(".")[-1]
        if nm in want:
            gf = vtk.vtkGeometryFilter(); gf.SetInputData(it.GetCurrentDataObject()); gf.Update()
            app.AddInputData(gf.GetOutput())
        it.GoToNextItem()
    app.Update()
    tr = vtk.vtkTransform(); tr.Scale(1000, 1000, 1000)
    tf = vtk.vtkTransformPolyDataFilter(); tf.SetInputData(app.GetOutput()); tf.SetTransform(tr); tf.Update()
    out = vtk.vtkPolyData(); out.DeepCopy(tf.GetOutput())
    return out


def loops(poly, o, n):
    """Closed cut loops of the wall with plane (o, n): list of (centre, area mm2, 2-D outline in plane coords)."""
    pl = vtk.vtkPlane(); pl.SetOrigin(*o); pl.SetNormal(*n)
    cut = vtk.vtkCutter(); cut.SetInputData(poly); cut.SetCutFunction(pl); cut.Update()
    cl = vtk.vtkCleanPolyData(); cl.SetInputConnection(cut.GetOutputPort())
    cl.ToleranceIsAbsoluteOn(); cl.SetAbsoluteTolerance(1e-4); cl.Update()           # joins the cut across block seams
    st = vtk.vtkStripper(); st.SetInputConnection(cl.GetOutputPort()); st.JoinContiguousSegmentsOn(); st.Update()
    pd = st.GetOutput()
    if pd.GetNumberOfPoints() == 0:
        return []
    P = vtk_to_numpy(pd.GetPoints().GetData())
    u = np.cross(n, [1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0]); u /= np.linalg.norm(u); v = np.cross(n, u)
    res = []
    lines = pd.GetLines(); lines.InitTraversal(); ids = vtk.vtkIdList()
    while lines.GetNextCell(ids):
        idx = [ids.GetId(k) for k in range(ids.GetNumberOfIds())]
        if len(idx) < 4 or idx[0] != idx[-1] and np.linalg.norm(P[idx[0]] - P[idx[-1]]) > 0.05:
            continue                                   # open polyline: not a lumen outline
        Q = P[idx] - o; x, y = Q @ u, Q @ v
        A = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
        res.append((P[idx].mean(axis=0), A, np.c_[x, y]))
    return res


def inside(pt, poly2):
    """Even-odd ray test: is 2-D point pt inside the closed outline poly2?"""
    x, y = poly2[:, 0], poly2[:, 1]; x2, y2 = np.roll(x, -1), np.roll(y, -1)
    cross = (y > pt[1]) != (y2 > pt[1])
    xs = x + (pt[1] - y) * (x2 - x) / np.where(y2 != y, y2 - y, 1e-30)
    return bool(np.sum(cross & (pt[0] < xs)) % 2)


def lumen_area(L):
    """Lumen area of the selected outlines: an outline inside an odd number of others is a hole
    (the epiglottis free edge sits inside the laryngeal inlet) and is subtracted."""
    tot = 0.0
    for i, (_, A, xy) in enumerate(L):
        depth = sum(inside(xy[0], L[j][2]) for j in range(len(L)) if j != i)
        tot += -A if depth % 2 else A
    return tot


def disk_reach(path, c):
    """Largest distance from the centroid to the slicer's own cross-section (its plane-STL disk)."""
    rd = vtk.vtkSTLReader(); rd.SetFileName(path); rd.Update()
    pts = rd.GetOutput().GetPoints()
    return float(np.linalg.norm(vtk_to_numpy(pts.GetData()) - c, axis=1).max()) if pts else None


os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
out = open(a.out, "w")
out.write("frame,cfd_idx,region,plane_index,arc_mm,area_table_mm2,area_recut_mm2,n_loops\n")
t0 = time.time()
for fn in fnames:
    fms = int(re.search(r"(\d+)", fn).group(1))
    ci = _bi.max(0, _bi.min(len(tv) - 1, int(round(fms / (CFD_DT * 1000.0))) - 1))
    poly = wall_poly(tv[ci])
    planes = frames[fn]
    cen = np.array([p[2] for p in planes])
    nrm = np.zeros_like(cen); nrm[1:-1] = cen[2:] - cen[:-2]; nrm[0] = cen[1] - cen[0]; nrm[-1] = cen[-1] - cen[-2]
    for i, (arc, area, c, pidx) in enumerate(planes):
        if not (a.arc_min <= arc <= a.arc_max):
            continue
        prank = RANK.get(fn, {}).get(pidx); n = None; rad = None
        if prank is not None:
            f = os.path.join(planes_dir, fn, "%s_Planes" % a.region, "%s-%s-%03d.stl" % (fn, tag, prank))
            got = stl_normal(f)
            if got is not None and np.linalg.norm(got[1] - c) <= 2.0:
                n = got[0]; rad = disk_reach(f, c) + a.reach_margin
        if n is None:
            n = nrm[i] / np.linalg.norm(nrm[i])
        if rad is None:                       # no disk: the extraction's lobe rule
            rad = _bi.max(a.lobe_radius_min, a.lobe_radius_factor * np.sqrt(_bi.max(area, 1.0) / np.pi))
        L = [lp for lp in loops(poly, c, n) if np.linalg.norm(lp[0] - c) <= rad]
        out.write("%s,%d,%s,%d,%.4f,%.4f,%.4f,%d\n" % (fn, ci, a.region, pidx, arc, area, lumen_area(L), len(L)))
    out.flush()
    print("  %s (cfd %d) done, %.0fs" % (fn, ci, time.time() - t0), flush=True)
out.close()
print("DONE", flush=True)
