"""Nostril cut plane for build_capped_surface.py (pvbatch).

Fits a plane (PCA) to the face/nose opening loop of the airway surface (the highest open boundary loop; the others
are bronchi), orients its normal INTO the nose, and slides it inward in 0.5 mm steps, printing the cut curves at each
depth. The first depth at which it cuts exactly two separate closed nostril rings is recorded; use that (or a little
more: ENT005 found 4.0, used 4.5) as build_capped_surface.py --nostril-depth.
Inward = the side where the plane cuts the airway (checked at +-3 mm); --inward-point x,y,z overrides
(any point inside the nasal cavity).
Usage: pvb.sh nostril_sweep.py airway_noface_mm.vtp [--out plane.json] [--inward-point x,y,z]
Writes <input>_nostril_plane.json by default: centre, normal (inward), rms_radius, first_two_ring_depth.
"""
import argparse, json
import numpy as np, vtk
from vtkmodules.util.numpy_support import vtk_to_numpy

ap = argparse.ArgumentParser(); ap.add_argument("surface"); ap.add_argument("--out", default=None)
ap.add_argument("--inward-point", default=None); ap.add_argument("--max-depth", type=float, default=14.0)
a = ap.parse_args()
rd = vtk.vtkXMLPolyDataReader(); rd.SetFileName(a.surface); rd.Update(); A = rd.GetOutput()
fe = vtk.vtkFeatureEdges(); fe.SetInputData(A); fe.BoundaryEdgesOn(); fe.FeatureEdgesOff(); fe.ManifoldEdgesOff(); fe.NonManifoldEdgesOff(); fe.Update()
cf = vtk.vtkPolyDataConnectivityFilter(); cf.SetInputConnection(fe.GetOutputPort()); cf.SetExtractionModeToAllRegions(); cf.ColorRegionsOn(); cf.Update()
E = cf.GetOutput(); P = vtk_to_numpy(E.GetPoints().GetData()); rid = vtk_to_numpy(E.GetPointData().GetArray('RegionId'))
loops = [P[rid == k] for k in range(cf.GetNumberOfExtractedRegions())]
nose = max(loops, key=lambda q: q[:, 2].mean())             # the opening loop is the highest one (others are bronchi)
c = nose.mean(0); u, s, vt = np.linalg.svd(nose - c); n = vt[2]
rad = np.sqrt(s[0] ** 2 / len(nose)) * 1.0


def rings(d, R, n):
    o = c + d * n; pl = vtk.vtkPlane(); pl.SetOrigin(*o); pl.SetNormal(*n)
    ct = vtk.vtkCutter(); ct.SetInputData(A); ct.SetCutFunction(pl); ct.Update()
    st = vtk.vtkStripper(); st.SetInputConnection(ct.GetOutputPort()); st.JoinContiguousSegmentsOn(); st.Update(); out = st.GetOutput()
    if out.GetNumberOfPoints() == 0: return []
    Q = vtk_to_numpy(out.GetPoints().GetData()); res = []
    for k in range(out.GetNumberOfCells()):
        cell = out.GetCell(k); q = Q[[cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]]
        m = q.mean(0); rr = np.linalg.norm((m - o) - np.dot(m - o, n) * n)
        if rr > R: continue                                  # only curves within the nostril footprint
        res.append(dict(closed=bool(np.linalg.norm(q[0] - q[-1]) < 1e-6), perim=float(np.linalg.norm(np.diff(q, axis=0), axis=1).sum()), centre=m.tolist(), pts=q))
    return res


if a.inward_point:
    cav = np.array([float(x) for x in a.inward_point.split(",")])
    if np.dot(n, cav - c) < 0: n = -n
else:
    L = lambda nn: sum(q["perim"] for q in rings(3.0, 1.3 * rad, nn))
    if L(-n) > L(n): n = -n
print('opening loop: %d pts, centre (%.1f, %.1f, %.1f), inward normal (%.3f, %.3f, %.3f), rms radius %.1f mm' % (len(nose), *c, *n, rad))
best = None
for d in np.arange(0.0, a.max_depth + 0.01, 0.5):
    r = rings(d, 1.3 * rad, n)
    desc = '  '.join('%s %.0fmm x=%.1f' % ('ring' if q['closed'] else 'OPEN', q['perim'], q['centre'][0]) for q in r)
    print('d=%4.1f mm: %d curves  %s' % (d, len(r), desc))
    big = [q for q in r if q['perim'] > 8]
    if best is None and len(big) == 2 and all(q['closed'] for q in big) and abs(big[0]['centre'][0] - big[1]['centre'][0]) > 3:
        best = d
out = a.out or a.surface.replace('.vtp', '_nostril_plane.json')
json.dump(dict(centre=c.tolist(), normal=n.tolist(), rms_radius=rad, first_two_ring_depth=best), open(out, 'w'), indent=1)
print('first depth with exactly two closed nostril rings: %s mm -> %s' % (best, out))
