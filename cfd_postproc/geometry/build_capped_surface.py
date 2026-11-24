"""Build the capped airway surface the CSA slicer needs (named caps), from an airway surface in mm.

Input : airway surface in mm with the face/mask already removed (.vtp), i.e. nose + nasal cavity +
        pharynx + trachea + bronchi, open at the face/nose opening and at the bronchi.
Cuts  : * nostrils - ONE plane parallel to the face/nose opening loop, slid --nostril-depth mm into
                     the nose (the depth at which it cuts two separate nostril rings),
        * trachea  - plane --trachea-above mm above the carina, perpendicular to the trachea axis.
Output: ASCII STL with solids Leftnose, Rightnose, outlet, wall (ENT003's seg_0_boundary.stl layout),
        all normals pointing out of the airway; plus a JSON record of the cut planes.
Run with pvbatch (needs vtk):  pvb.sh build_capped_surface.py in.vtp out.stl --nostril-plane plane.json
"""
import argparse
import json

import numpy as np
import vtk
from vtkmodules.util.numpy_support import vtk_to_numpy

ap = argparse.ArgumentParser()
ap.add_argument('surface'); ap.add_argument('out_stl')
ap.add_argument('--nostril-plane', required=True, help='JSON with centre + inward normal of the face/nose opening loop')
ap.add_argument('--nostril-depth', type=float, default=4.5, help='mm into the nose from the opening plane')
ap.add_argument('--nostril-radius', type=float, default=16.0, help='mm; only cut within this radius of the nostril axis')
ap.add_argument('--carina-z', type=float, required=True, help='z of the carina (mm)')
ap.add_argument('--trachea-above', type=float, default=10.0, help='mm above the carina for the outlet cut')
ap.add_argument('--fill-pinholes-mm', type=float, default=0.0,
                help='close open rings with perimeter below this (mm) into the wall: sub-mm pinholes a CFD-wall export can have (0 = off)')
ap.add_argument('--decimate', type=float, default=0.9, help='fraction of wall triangles to remove (0 = keep all)')
a = ap.parse_args()

rd = vtk.vtkXMLPolyDataReader(); rd.SetFileName(a.surface); rd.Update(); A = rd.GetOutput()
print('input surface: %d triangles' % A.GetNumberOfCells())

# --- nostril cut: keep the inward side of the plane, but only near the nostril axis
P = json.load(open(a.nostril_plane)); c0 = np.array(P['centre']); n = np.array(P['normal']) / np.linalg.norm(P['normal'])
o_n = c0 + a.nostril_depth * n
keep_plane = vtk.vtkPlane(); keep_plane.SetOrigin(*o_n); keep_plane.SetNormal(*n)           # > 0 on the nasal-cavity side
cyl = vtk.vtkCylinder(); cyl.SetCenter(*o_n); cyl.SetAxis(*n); cyl.SetRadius(a.nostril_radius)  # < 0 inside the cylinder
f_n = vtk.vtkImplicitBoolean(); f_n.SetOperationTypeToIntersection(); f_n.AddFunction(keep_plane); f_n.AddFunction(cyl)

# --- trachea cut: axis from cross-section centroids 5..30 mm above the carina
X = vtk_to_numpy(A.GetPoints().GetData())
cents = []
for z in np.arange(a.carina_z + 5, a.carina_z + 30, 2.5):
    pl = vtk.vtkPlane(); pl.SetOrigin(0, 0, z); pl.SetNormal(0, 0, 1)
    ct = vtk.vtkCutter(); ct.SetInputData(A); ct.SetCutFunction(pl); ct.Update()
    Q = vtk_to_numpy(ct.GetOutput().GetPoints().GetData()); cents.append(Q.mean(0))
cents = np.array(cents); cc = cents.mean(0); ax_t = np.linalg.svd(cents - cc)[2][0]
if ax_t[2] < 0: ax_t = -ax_t                                                                     # pointing up the airway
o_t = cc + ax_t * ((a.carina_z + a.trachea_above - cc[2]) / ax_t[2])
f_t = vtk.vtkPlane(); f_t.SetOrigin(*o_t); f_t.SetNormal(*ax_t)                                  # > 0 above the cut = keep
print('trachea axis (%.3f, %.3f, %.3f), outlet cut through (%.1f, %.1f, %.1f)' % (*ax_t, *o_t))

clip1 = vtk.vtkClipPolyData(); clip1.SetInputData(A); clip1.SetClipFunction(f_n); clip1.Update()
clip2 = vtk.vtkClipPolyData(); clip2.SetInputConnection(clip1.GetOutputPort()); clip2.SetClipFunction(f_t); clip2.Update()
tri = vtk.vtkTriangleFilter(); tri.SetInputConnection(clip2.GetOutputPort()); tri.Update()
cl = vtk.vtkCleanPolyData(); cl.SetInputConnection(tri.GetOutputPort()); cl.Update()

# --- keep the connected piece that holds the trachea
conn = vtk.vtkPolyDataConnectivityFilter(); conn.SetInputConnection(cl.GetOutputPort()); conn.SetExtractionModeToAllRegions(); conn.ColorRegionsOn(); conn.Update()
sizes = vtk_to_numpy(conn.GetRegionSizes())
cs = vtk.vtkPolyDataConnectivityFilter(); cs.SetInputConnection(cl.GetOutputPort()); cs.SetExtractionModeToClosestPointRegion(); cs.SetClosestPoint(*(o_t + 3 * ax_t)); cs.Update()
W = cs.GetOutput()
print('after cuts: %d pieces (sizes %s); kept the airway piece: %d triangles' % (len(sizes), sorted(sizes.tolist(), reverse=True)[:6], W.GetNumberOfCells()))

if a.decimate > 0:
    # DecimatePro with boundary vertices locked: the cut rings (future caps) stay exactly on their planes
    dc = vtk.vtkDecimatePro(); dc.SetInputData(W); dc.SetTargetReduction(a.decimate); dc.PreserveTopologyOn()
    dc.BoundaryVertexDeletionOff(); dc.SplittingOff(); dc.Update()
    cl2 = vtk.vtkCleanPolyData(); cl2.SetInputConnection(dc.GetOutputPort()); cl2.Update(); Wd = cl2.GetOutput()
    loc = vtk.vtkCellLocator(); loc.SetDataSet(W); loc.BuildLocator()
    Y = vtk_to_numpy(Wd.GetPoints().GetData()); cp = [0.0] * 3; ci = vtk.reference(0); si = vtk.reference(0); d2 = vtk.reference(0.0)
    dev = np.empty(min(len(Y), 30000)); idx = np.random.default_rng(0).choice(len(Y), len(dev), replace=False)
    for j, i in enumerate(idx):
        loc.FindClosestPoint(Y[i], cp, ci, si, d2); dev[j] = d2.get() ** 0.5
    print('decimated wall: %d -> %d triangles; deviation from the full surface: median %.3f, p99 %.3f, max %.3f mm'
          % (W.GetNumberOfCells(), Wd.GetNumberOfCells(), np.median(dev), np.percentile(dev, 99), dev.max()))
    W = Wd

# --- the open rings: 2 nostrils + 1 outlet
fe = vtk.vtkFeatureEdges(); fe.SetInputData(W); fe.BoundaryEdgesOn(); fe.FeatureEdgesOff(); fe.ManifoldEdgesOff(); fe.NonManifoldEdgesOff(); fe.Update()
st = vtk.vtkStripper(); st.SetInputConnection(fe.GetOutputPort()); st.JoinContiguousSegmentsOn(); st.Update(); R = st.GetOutput()
Rp = vtk_to_numpy(R.GetPoints().GetData())
rings = []
for k in range(R.GetNumberOfCells()):
    cell = R.GetCell(k); q = Rp[[cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]]
    if np.linalg.norm(q[0] - q[-1]) < 1e-6: q = q[:-1]
    rings.append(q)
pin = []
if a.fill_pinholes_mm > 0:
    per = [float(np.linalg.norm(np.diff(np.vstack([q, q[:1]]), axis=0), axis=1).sum()) for q in rings]
    pin = [q for q, p in zip(rings, per) if p < a.fill_pinholes_mm]; rings = [q for q, p in zip(rings, per) if p >= a.fill_pinholes_mm]
    print('filled %d pinholes (perimeter < %g mm) into the wall' % (len(pin), a.fill_pinholes_mm))
print('open rings: %d  (perimeters %s mm)' % (len(rings), [round(float(np.linalg.norm(np.diff(np.vstack([q, q[:1]]), axis=0), axis=1).sum()), 1) for q in rings]))
assert len(rings) == 3, 'expected exactly 3 open rings (2 nostrils + outlet), got %d' % len(rings)
outlet = min(rings, key=lambda q: q[:, 2].mean())
nost = sorted([q for q in rings if q is not outlet], key=lambda q: q[:, 0].mean())              # smaller x first = Left (RAS, as ENT003)

def fan(q, out_normal):
    """Cap a ring with a fan from its centroid; winding gives a normal along out_normal."""
    cen = q.mean(0); tris = []
    for i in range(len(q)):
        t = np.array([cen, q[i], q[(i + 1) % len(q)]])
        if np.dot(np.cross(t[1] - t[0], t[2] - t[0]), out_normal) < 0: t = t[[0, 2, 1]]
        tris.append(t)
    return np.array(tris)

caps = {'Leftnose': fan(nost[0], -n), 'Rightnose': fan(nost[1], -n), 'outlet': fan(outlet, -ax_t)}
Wt = vtk_to_numpy(W.GetPoints().GetData())
Wc = vtk_to_numpy(W.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
wall = Wt[Wc]
# wall normals must point out of the airway too: test against the inward nostril-plane direction near the nose
nrm = np.cross(wall[:, 1] - wall[:, 0], wall[:, 2] - wall[:, 0])
centre_air = (o_t + 25 * ax_t)                                     # a point inside the trachea lumen
near = np.linalg.norm(wall.mean(1) - centre_air, axis=1) < 12
if np.mean(np.einsum('ij,ij->i', nrm[near], wall[near].mean(1) - centre_air)) < 0: wall = wall[:, [0, 2, 1]]
for q in pin:                                                      # pinhole fans, wound like their neighbours
    cen = q.mean(0); T = np.array([[cen, q[i], q[(i + 1) % len(q)]] for i in range(len(q))])
    jj = np.argmin(np.linalg.norm(wall.reshape(-1, 3) - q[0], axis=1)) // 3
    nb = np.cross(wall[jj, 1] - wall[jj, 0], wall[jj, 2] - wall[jj, 0])
    if np.dot(np.cross(T[0, 1] - T[0, 0], T[0, 2] - T[0, 0]), nb) < 0: T = T[:, [0, 2, 1]]
    wall = np.vstack([wall, T])

def write_solid(fh, name, T):
    fh.write('solid %s\n' % name)
    N = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); N /= np.linalg.norm(N, axis=1, keepdims=True) + 1e-30
    for t, nn in zip(T, N):
        fh.write('  facet normal %e %e %e\n    outer loop\n' % tuple(nn))
        for v in t: fh.write('      vertex %e %e %e\n' % tuple(v))
        fh.write('    endloop\n  endfacet\n')
    fh.write('endsolid %s\n' % name)

with open(a.out_stl, 'w') as fh:
    for name in ('Leftnose', 'Rightnose', 'outlet'): write_solid(fh, name, caps[name])
    write_solid(fh, 'wall', wall)
rec = dict(nostril_plane=dict(origin=o_n.tolist(), inward_normal=n.tolist(), depth_mm=a.nostril_depth, radius_mm=a.nostril_radius),
           outlet_plane=dict(origin=o_t.tolist(), up_normal=ax_t.tolist(), carina_z=a.carina_z, above_mm=a.trachea_above),
           caps={k: dict(centre=v[:, 0].mean(0).tolist(), n_tris=len(v)) for k, v in caps.items()}, wall_tris=len(wall))
json.dump(rec, open(a.out_stl.replace('.stl', '_cuts.json'), 'w'), indent=1)
for k, v in caps.items():
    print('cap %-9s centre (%.1f, %.1f, %.1f) mm, %d triangles' % (k, *v[:, 0].mean(0), len(v)))
print('wrote %s (wall %d triangles)' % (a.out_stl, len(wall)))
