"""Same planes, two geometries (pvbatch): every valid BASELINE plane on the nostril -> trachea path (LeftNose +
DescendingAirway; origin = the baseline table centroid, normal = the baseline slicer's own plane disk) is cut through the
registration-deformed surface of the baseline AND of the edited geometry in the same frame. Distance along the
centerline is the baseline's for both, so the only difference is the geometry (the slicer's own septum split can move
between two surfaces - ENT005 Aryt: left nose 19k vs 45k faces - so the edited faces get the BASELINE region labels by
nearest face at t0). LeftNose planes cut the left-labelled faces, DescendingAirway planes the whole surface (as the
slicer: left + right + desc). Closed outlines whose centre lies within the baseline disk's reach + --margin are
summed, a nested outline (epiglottis inside the inlet) subtracted. The baseline re-cut reproduces the slicer table
(ENT005: ratio 1.000 on 99% of planes); use area_edit / area_base as the edit factor.
Usage: pvb.sh samplane_recut.py --base-dir SUBJ_BASE --edit-dir SUBJ_EDIT --out-dir DIR [--shard w:N]   -> DIR/recut_<w>.csv
(both subject dirs = csa-bifurcation runs with the SAME registration; run_samplane.sh runs all shards + merges)
"""
import argparse, csv, glob, os, re, sys, numpy as np, vtk
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
ap = argparse.ArgumentParser(); ap.add_argument("--base-dir", required=True); ap.add_argument("--edit-dir", required=True)
ap.add_argument("--out-dir", required=True); ap.add_argument("--shard", default="0:1"); ap.add_argument("--margin", type=float, default=5.0)
a = ap.parse_args()
MARGIN = a.margin; OUTD = a.out_dir; os.makedirs(OUTD, exist_ok=True)
TAG = {"LeftNose": "Left", "DescendingAirway": "Desc"}
w, nsh = (int(x) for x in a.shard.split(":"))
CSA = [f for f in glob.glob(a.base_dir + "/csa_bifurcation/*-Data.csv")][0]
LABELS = CSA.replace("-Data.csv", "_full_face_labels.npy")


def read_stl(p):
    r = vtk.vtkSTLReader(); r.SetFileName(p); r.Update(); o = vtk.vtkPolyData(); o.DeepCopy(r.GetOutput()); return o


def centres(poly):
    f = vtk.vtkCellCenters(); f.SetInputData(poly); f.Update(); return vtk_to_numpy(f.GetOutput().GetPoints().GetData()).astype(float)



def submesh(poly, mask):
    p = vtk.vtkPolyData(); p.DeepCopy(poly)
    a = numpy_to_vtk(mask.astype(np.uint8), deep=True); a.SetName("m"); p.GetCellData().AddArray(a)
    th = vtk.vtkThreshold(); th.SetInputData(p); th.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "m")
    th.SetLowerThreshold(0.5); th.SetUpperThreshold(1.5); th.Update()
    g = vtk.vtkGeometryFilter(); g.SetInputData(th.GetOutput()); g.Update(); o = vtk.vtkPolyData(); o.DeepCopy(g.GetOutput()); return o


def loops(poly, o, n):
    pl = vtk.vtkPlane(); pl.SetOrigin(*o); pl.SetNormal(*n)
    cut = vtk.vtkCutter(); cut.SetInputData(poly); cut.SetCutFunction(pl); cut.Update()
    cl = vtk.vtkCleanPolyData(); cl.SetInputConnection(cut.GetOutputPort()); cl.ToleranceIsAbsoluteOn(); cl.SetAbsoluteTolerance(1e-4); cl.Update()
    st = vtk.vtkStripper(); st.SetInputConnection(cl.GetOutputPort()); st.JoinContiguousSegmentsOn(); st.Update(); pd = st.GetOutput()
    if pd.GetNumberOfPoints() == 0: return []
    P = vtk_to_numpy(pd.GetPoints().GetData()).astype(float)
    u = np.cross(n, [1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0]); u /= np.linalg.norm(u); v = np.cross(n, u)
    res = []; lines = pd.GetLines(); lines.InitTraversal(); ids = vtk.vtkIdList()
    while lines.GetNextCell(ids):
        idx = [ids.GetId(k) for k in range(ids.GetNumberOfIds())]
        if len(idx) < 4 or (idx[0] != idx[-1] and np.linalg.norm(P[idx[0]] - P[idx[-1]]) > 0.05): continue
        Q = P[idx] - o; x, y = Q @ u, Q @ v
        res.append((P[idx].mean(0), 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))), np.c_[x, y]))
    return res


def inside(pt, xy):
    x, y = xy[:, 0], xy[:, 1]; x2, y2 = np.roll(x, -1), np.roll(y, -1); c = (y > pt[1]) != (y2 > pt[1])
    xs = x + (pt[1] - y) * (x2 - x) / np.where(y2 != y, y2 - y, 1e-30); return bool(np.sum(c & (pt[0] < xs)) % 2)


def lumen(L):
    return sum((-A if sum(inside(xy[0], L[j][2]) for j in range(len(L)) if j != i) % 2 else A) for i, (_, A, xy) in enumerate(L))


def disk(path):
    """(normal, reach-point set) of a slicer plane STL"""
    r = read_stl(path); V = vtk_to_numpy(r.GetPoints().GetData()).astype(float) if r.GetNumberOfPoints() else None
    if V is None or len(V) < 3: return None, None
    c = V.mean(0); return np.linalg.svd(V - c, full_matrices=False)[2][2], V


# plane table + plane-STL ranks (rank over -Data-all.csv, per frame and region)
rows = [d for d in csv.DictReader(open(CSA)) if d["is_valid"].lower() == "true" and d["region"] in TAG]
rank = {}
for d in csv.DictReader(open(CSA.replace("-Data.csv", "-Data-all.csv"))):
    if d["region"] in TAG: rank.setdefault((d["frame_name"], d["region"]), []).append(int(d["plane_index"]))
rank = {k: {p: i for i, p in enumerate(sorted(v))} for k, v in rank.items()}
frames = sorted({d["frame_name"] for d in rows}, key=lambda s: int(s[4:]))

# baseline labels -> Aryt faces by nearest face centre at t0
LB = np.load(LABELS).astype(int)
b0, a0 = read_stl(a.base_dir + "/motion/stl/out_000000.stl"), read_stl(a.edit_dir + "/motion/stl/out_000000.stl")
assert b0.GetNumberOfCells() == len(LB), (b0.GetNumberOfCells(), len(LB))
pts = vtk.vtkPoints(); pts.SetData(numpy_to_vtk(centres(b0), deep=True)); pp = vtk.vtkPolyData(); pp.SetPoints(pts)
loc = vtk.vtkStaticPointLocator(); loc.SetDataSet(pp); loc.BuildLocator()
LA = LB[[loc.FindClosestPoint(c) for c in centres(a0)]]
print("labels base", np.bincount(LB), "edit (transferred)", np.bincount(LA), flush=True)

out = open("%s/recut_%02d.csv" % (OUTD, w), "w"); wr = csv.writer(out)
wr.writerow(["frame", "ms", "region", "plane_index", "arc_mm", "area_table", "area_base", "area_edit", "n_base", "n_edit"])
for fn in frames[w::nsh]:
    ms = int(fn[4:]); S = {}
    for run, lab in (("base", LB), ("edit", LA)):
        p = read_stl("%s/motion/stl/%s.stl" % (a.base_dir if run == "base" else a.edit_dir, fn))
        assert p.GetNumberOfCells() == len(lab), (run, fn, p.GetNumberOfCells(), len(lab))
        S[run] = {"LeftNose": submesh(p, lab == 1), "DescendingAirway": p}
    for d in (x for x in rows if x["frame_name"] == fn):
        reg, pi = d["region"], int(d["plane_index"]); c = np.array([float(d["centroid_" + k]) for k in "xyz"])
        rk = rank.get((fn, reg), {}).get(pi)
        if rk is None: continue
        n, V = disk("%s/csa_bifurcation/frames/%s/%s_Planes/%s-%s-%03d.stl" % (a.base_dir, fn, reg, fn, TAG[reg], rk))
        if n is None: continue
        reach = np.linalg.norm(V - c, axis=1).max() + MARGIN
        A = {}
        for run in ("base", "edit"):
            L = [lp for lp in loops(S[run][reg], c, n) if np.linalg.norm(lp[0] - c) <= reach]; A[run] = (lumen(L), len(L))
        wr.writerow([fn, ms, reg, pi, d["arc_length_mm"], d["area_mm2"], "%.4f" % A["base"][0], "%.4f" % A["edit"][0], A["base"][1], A["edit"][1]])
    out.flush(); print(fn, "done", flush=True)
out.close(); print("DONE", flush=True)
