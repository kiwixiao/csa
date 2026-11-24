"""Fill the nostril end of each nasal path with CSA planes that turn from the cap to the centerline.

Near the nostril the CSA slicer's planes (normal to the centerline) cross the tilted cap and fail
(ENT005: no usable plane for the first ~6-8 mm). Here, for every CSA frame and each nostril:
  * plane at s = 0 is the CAP itself (its ring in motion/stl/<frame>.stl, centre + PCA normal),
  * the first good slicer plane (its STL disk, PCA normal) sits at centerline distance s1,
  * planes in between sit on the per-side centerline every --step mm, normal blended linearly
    (renormalised) from the cap's to the first good plane's: they start parallel to the cap and
    turn gradually until they are normal to the centerline.
Each plane cuts that nostril's OWN passage surface (slicer face labels, closed with its cap) and,
like the slicer, the CSA is the SUM of every lumen piece of that passage near the plane (the
vestibule is cut into several channels here). The same plane cuts the CFD; the CFD pieces whose
centres fall inside those outlines are integrated (TP, Q, energy flux; outward normal, as in
csa_pressure_energy.py). Check printed: the first good slicer plane re-cut this way vs the
slicer's own area.
Writes fill_csa.csv (CSA-table schema + s_from_cap_mm), fill_pe.csv (csa_pe schema), and per-frame
disk STLs. Run with pvbatch.
"""
import argparse, csv, os, time
import numpy as np
import vtk
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk, numpy_to_vtkIdTypeArray

ap = argparse.ArgumentParser()
ap.add_argument('--subject-dir', required=True); ap.add_argument('--csa', required=True, help='ORIGINAL <S>-Data.csv (plane-STL ranks)')
ap.add_argument('--labels', required=True, help='<S>_full_face_labels.npy (1 = left nose, 2 = right nose)')
ap.add_argument('--rejected', required=True, help='rejected_planes.csv from clean_nostril_planes.py')
ap.add_argument('--case', required=True); ap.add_argument('--out-dir', required=True)
ap.add_argument('--step', type=float, default=1.0); ap.add_argument('--near-mm', type=float, default=15.0)
ap.add_argument('--cfd-dt-ms', type=float, default=25.0); ap.add_argument('--geometry-only', action='store_true')
a = ap.parse_args()
SIDES = {'LeftNose': ('leftnose', 'Left', 1), 'RightNose': ('rightnose', 'Right', 2)}
LAB = np.load(a.labels)
os.makedirs(os.path.join(a.out_dir, 'fill_disks'), exist_ok=True)

def vtk_points(path):
    tok = open(path, errors='ignore').read().split(); i = tok.index('POINTS'); n = int(tok[i + 1])
    return np.array(tok[i + 3:i + 3 + 3 * n], float).reshape(n, 3)

def stl_tris(path):
    raw = open(path, 'rb').read(); n = int(np.frombuffer(raw[80:84], '<u4')[0])
    dt = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
    return np.frombuffer(raw[84:], dtype=dt, count=n)['v'].astype(float)

def poly_from_tris(T):
    V, inv = np.unique(np.round(T.reshape(-1, 3), 6), axis=0, return_inverse=True); F = inv.reshape(-1, 3)
    pts = vtk.vtkPoints(); pts.SetData(numpy_to_vtk(V.copy(), deep=1))
    ca = vtk.vtkCellArray(); ca.SetCells(len(F), numpy_to_vtkIdTypeArray(np.c_[np.full(len(F), 3), F].ravel().astype(np.int64), deep=1))
    p = vtk.vtkPolyData(); p.SetPoints(pts); p.SetPolys(ca); return p

def open_rings(T):
    """Ordered boundary rings (closed polylines) of a triangle soup."""
    p = poly_from_tris(T)
    fe = vtk.vtkFeatureEdges(); fe.SetInputData(p); fe.BoundaryEdgesOn(); fe.FeatureEdgesOff(); fe.ManifoldEdgesOff(); fe.NonManifoldEdgesOff(); fe.Update()
    st = vtk.vtkStripper(); st.SetInputConnection(fe.GetOutputPort()); st.JoinContiguousSegmentsOn(); st.Update(); o = st.GetOutput()
    P = vtk_to_numpy(o.GetPoints().GetData())
    return [P[[o.GetCell(k).GetPointId(i) for i in range(o.GetCell(k).GetNumberOfPoints())]] for k in range(o.GetNumberOfCells())]

def pca_plane(P):
    c = P.mean(0); n = np.linalg.svd(P - c)[2][2]; return c, n / np.linalg.norm(n)

def basis(n):
    u = np.cross(n, [1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0]); u /= np.linalg.norm(u); return u, np.cross(n, u)

def area2d(q2):
    x, y = q2[:, 0], q2[:, 1]; return 0.5 * abs(np.dot(x[:-1], y[1:]) - np.dot(x[1:], y[:-1]))

def inside2d(pt, q2):
    x, y = q2[:, 0], q2[:, 1]; px, py = pt
    c = ((y[:-1] > py) != (y[1:] > py)) & (px < x[:-1] + (py - y[:-1]) * (x[1:] - x[:-1]) / (y[1:] - y[:-1] + 1e-30))
    return c.sum() % 2 == 1

def front_piece(Ts, p1, n1):
    """The passage surface between the cap and the first good slicer plane, closed at both ends:
    clip at plane 1 (keep the nostril side; n1 points OUT of the nose), then triangulate every
    open loop (the cap ring and the plane-1 outline, which may be several pieces)."""
    pl = vtk.vtkPlane(); pl.SetOrigin(*p1); pl.SetNormal(*n1)
    cp = vtk.vtkClipPolyData(); cp.SetInputData(poly_from_tris(Ts)); cp.SetClipFunction(pl); cp.Update()
    tf = vtk.vtkTriangleFilter(); tf.SetInputConnection(cp.GetOutputPort()); tf.Update()
    cl = vtk.vtkCleanPolyData(); cl.SetInputConnection(tf.GetOutputPort()); cl.Update(); S = cl.GetOutput()
    F = vtk_to_numpy(S.GetPolys().GetData()).reshape(-1, 4)[:, 1:]; T2 = vtk_to_numpy(S.GetPoints().GetData())[F]
    caps = []
    for q in open_rings(T2):
        q_ = q[:-1] if np.linalg.norm(q[0] - q[-1]) < 1e-9 else q
        pts = vtk.vtkPoints(); [pts.InsertNextPoint(*x) for x in q_]; ln = vtk.vtkCellArray(); ln.InsertNextCell(len(q_) + 1)
        for i in range(len(q_)): ln.InsertCellPoint(i)
        ln.InsertCellPoint(0)
        pd = vtk.vtkPolyData(); pd.SetPoints(pts); pd.SetLines(ln)
        ctri = vtk.vtkContourTriangulator(); ctri.SetInputData(pd); ctri.Update(); o = ctri.GetOutput()
        if o.GetNumberOfCells():
            Fc = vtk_to_numpy(o.GetPolys().GetData()).reshape(-1, 4)[:, 1:]; caps.append(vtk_to_numpy(o.GetPoints().GetData())[Fc])
        else:                                   # fallback: fan
            c = q_.mean(0); caps.append(np.array([[c, q_[i], q_[(i + 1) % len(q_)]] for i in range(len(q_))]))
    return poly_from_tris(np.vstack([T2] + caps))

def section(surf, origin, normal, near):
    """All closed loops of surf cut by the plane whose centroid is within `near` mm of origin:
    summed area / perimeter, area-weighted centroid, and the loops."""
    pl = vtk.vtkPlane(); pl.SetOrigin(*origin); pl.SetNormal(*normal)
    ct = vtk.vtkCutter(); ct.SetInputData(surf); ct.SetCutFunction(pl); ct.Update()
    st = vtk.vtkStripper(); st.SetInputConnection(ct.GetOutputPort()); st.JoinContiguousSegmentsOn(); st.Update(); o = st.GetOutput()
    if o.GetNumberOfPoints() == 0: return None
    P = vtk_to_numpy(o.GetPoints().GetData()); u, v = basis(normal); loops = []
    for k in range(o.GetNumberOfCells()):
        q = P[[o.GetCell(k).GetPointId(i) for i in range(o.GetCell(k).GetNumberOfPoints())]]
        if len(q) < 4 or np.linalg.norm(q[0] - q[-1]) > 1e-6 or np.linalg.norm(q[:-1].mean(0) - origin) > near: continue
        q2 = np.c_[(q - origin) @ u, (q - origin) @ v]; A = area2d(q2)
        if A > 0.05: loops.append((q, q2, A))
    if not loops: return None
    A = sum(l[2] for l in loops)
    return dict(loops=loops, area=A, perim=sum(float(np.linalg.norm(np.diff(l[0], axis=0), axis=1).sum()) for l in loops),
                centroid=sum(l[0][:-1].mean(0) * l[2] for l in loops) / A, u=u, v=v)

# ---- slicer planes that survive cleaning, and their STL rank -------------------------------
valid = {}; arcs = {}; sarea = {}
for d in csv.DictReader(open(a.csa)):
    if d['is_valid'].lower() == 'true' and d['region'] in SIDES:
        k = (d['frame_name'], d['region'], int(d['plane_index'])); valid.setdefault(k[:2], []).append(k[2])
        arcs[k] = float(d['arc_length_mm']); sarea[k] = float(d['area_mm2'])
rej = {(d['frame_name'], d['region'], int(d['plane_index'])) for d in csv.DictReader(open(a.rejected))}
allp = {}
for d in csv.DictReader(open(a.csa.replace('-Data.csv', '-Data-all.csv'))):      # disks are numbered over ALL planes
    allp.setdefault((d['frame_name'], d['region']), []).append(int(d['plane_index']))
frames = sorted({k[0] for k in valid}, key=lambda f: int(f.split('_')[1]))

geo = {}; check = []
t0 = time.time()
for fn in frames:
    T = stl_tris(os.path.join(a.subject_dir, 'motion', 'stl', fn + '.stl'))
    assert len(T) == len(LAB), 'face labels do not match %s' % fn
    disks = []
    for reg, (side, tag, lab) in SIDES.items():
        Ts = T[LAB == lab]
        P = vtk_points(os.path.join(a.subject_dir, 'motion', 'centerlines', side, fn + '.vtk'))
        if P[0, 2] < P[-1, 2]: P = P[::-1]                                      # cap end first
        sc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
        t_out = P[0] - P[min(4, len(P) - 1)]; t_out /= np.linalg.norm(t_out)  # pointing OUT of the nose
        ring = min(open_rings(Ts), key=lambda q: np.linalg.norm(q.mean(0) - P[0]))   # this nostril's cap ring
        rq = ring if np.linalg.norm(ring[0] - ring[-1]) < 1e-9 else np.vstack([ring, ring[:1]])
        c0, n0 = pca_plane(rq[:-1]); n0 = n0 if n0 @ t_out > 0 else -n0
        # first good slicer plane
        pl = sorted(valid[(fn, reg)]); good = [p for p in pl if (fn, reg, p) not in rej]
        p1 = max(good, key=lambda p: arcs[(fn, reg, p)]); rank = sorted(allp[(fn, reg)]).index(p1)
        D = stl_tris(os.path.join(a.subject_dir, 'csa_bifurcation', 'frames', fn, '%s_Planes' % reg, '%s-%s-%03d.stl' % (fn, tag, rank))).reshape(-1, 3)
        c1, n1 = pca_plane(D)
        sd = (P - c1) @ n1; cross = np.flatnonzero(np.sign(sd[:-1]) != np.sign(sd[1:]))
        j = cross[np.argmin(np.linalg.norm(P[cross] - c1, axis=1))]; f = sd[j] / (sd[j] - sd[j + 1])
        s1 = sc[j] + f * (sc[j + 1] - sc[j]); p_1 = P[j] + f * (P[j + 1] - P[j])
        n1 = n1 if n1 @ (P[j] - P[j + 1]) > 0 else -n1
        chk = section(poly_from_tris(Ts), p_1, n1, 1e9)                       # whole passage at plane 1 = the slicer's sum
        surf = front_piece(Ts, p_1, n1)                                       # cap .. plane 1, closed
        check.append((fn, reg, sarea[(fn, reg, p1)], chk['area'] if chk else np.nan))
        arc_cap = arcs[(fn, reg, p1)] + s1
        u0, v0 = basis(n0); q2c = np.c_[(rq - c0) @ u0, (rq - c0) @ v0]
        cap = dict(loops=[(rq, q2c, area2d(q2c))], area=area2d(q2c), perim=float(np.linalg.norm(np.diff(rq, axis=0), axis=1).sum()),
                   centroid=c0, u=u0, v=v0)
        out = [dict(s=0.0, origin=c0, normal=n0, arc=arc_cap, **cap)]
        for s in np.arange(a.step, s1 - 0.4 * a.step, a.step):
            w = s / s1; n = (1 - w) * n0 + w * n1; n /= np.linalg.norm(n)
            o = np.array([np.interp(s, sc, P[:, i]) for i in range(3)])
            sec = section(surf, o, n, 1e9)                                     # every piece between cap and plane 1
            if sec is None: continue
            out.append(dict(s=float(s), origin=o, normal=n, arc=arc_cap - s, **sec))
        geo[(fn, reg)] = out
        disks += [l[0] for p in out for l in p['loops']]
    Tq = []
    for q in disks:
        c = q[:-1].mean(0); Tq += [[c, q[i], q[i + 1]] for i in range(len(q) - 1)]
    Tq = np.array(Tq, '<f4'); dt = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')]); rec = np.zeros(len(Tq), dt); rec['v'] = Tq
    with open(os.path.join(a.out_dir, 'fill_disks', fn + '-fill.stl'), 'wb') as fh:
        fh.write(b'nostril transition planes'.ljust(80, b' ')); fh.write(np.uint32(len(rec)).tobytes()); fh.write(rec.tobytes())
C = np.array([(c[2], c[3]) for c in check])
print('geometry: %d frames x 2 nostrils in %.0fs' % (len(frames), time.time() - t0), flush=True)
print('check - first good slicer plane re-cut here vs the slicer: median |diff| %.1f%%, max %.1f%%'
      % (100 * np.nanmedian(abs(C[:, 1] / C[:, 0] - 1)), 100 * np.nanmax(abs(C[:, 1] / C[:, 0] - 1))), flush=True)

with open(os.path.join(a.out_dir, 'fill_csa.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['plane_index', 'region', 'arc_length_mm', 'arc_length_branch_mm', 'area_mm2', 'perimeter_mm', 'hydraulic_diameter_mm',
                                    'equivalent_diameter_mm', 'centroid_x', 'centroid_y', 'centroid_z', 'centerline_source', 'is_valid', 'frame_name', 'frame_index',
                                    's_from_cap_mm', 'n_pieces'])
    for fn in frames:
        for reg in SIDES:
            for k, p in enumerate(geo[(fn, reg)]):
                w.writerow([1000 + k, reg, '%.6f' % p['arc'], '%.6f' % p['arc'], '%.6f' % p['area'], '%.6f' % p['perim'], '%.6f' % (4 * p['area'] / p['perim']),
                            '%.6f' % (2 * np.sqrt(p['area'] / np.pi)), *('%.6f' % x for x in p['centroid']),
                            'nostril_cap' if k == 0 else 'nostril_transition', 'True', fn, frames.index(fn), '%.3f' % p['s'], len(p['loops'])])
if a.geometry_only:
    raise SystemExit(0)

# ---- CFD on the same planes ------------------------------------------------------------------
from paraview.simple import *
from paraview import servermanager as sm
from paraview.vtk.numpy_interface import dataset_adapter as dsa
r = OpenDataFile(a.case); r.PointArrays = ['TotalPressure', 'EnergyFlux', 'Velocity']; tv = list(r.TimestepValues)
mb = MergeBlocks(Input=r); slc = Slice(Input=mb); slc.SliceType = 'Plane'; conn = Connectivity(Input=slc)
thr = Threshold(Input=conn, Scalars=['POINTS', 'RegionId']); calc = Calculator(Input=thr); calc.AttributeType = 'Point Data'; calc.ResultArrayName = 'TPvn'
iv = IntegrateVariables(Input=calc)
cfd_of = lambda fn: min(max(int(round(int(fn.split('_')[1]) / a.cfd_dt_ms)) - 1, 0), len(tv) - 1)
fh = open(os.path.join(a.out_dir, 'fill_pe.csv'), 'w')
fh.write('frame,frame_ms,cfd_idx,cfd_t_s,region,arc_mm,area_mm2,TP_area_avg_Pa,TP_mdot_avg_Pa,energy_flux_W,Q_Lmin,area_int_mm2,n_lobes\n')
t0 = time.time(); nomatch = 0
for i, fn in enumerate(frames):
    k = cfd_of(fn); ct = tv[k]
    for reg in SIDES:
        for p in geo[(fn, reg)]:
            o, n = p['origin'], p['normal']
            slc.SliceType.Origin = list(o * 1e-3); slc.SliceType.Normal = list(n); conn.UpdatePipeline(ct)
            d = dsa.WrapDataObject(sm.Fetch(conn)); pts = np.asarray(d.Points) * 1e3
            if pts.size == 0: continue
            rid = np.asarray(d.PointData['RegionId']); u, v = p['u'], p['v']
            pieces = []
            for R in np.unique(rid):
                # the CSA outline comes from the moving CSA surface, the CFD nose barely moves: accept a CFD
                # piece whose centre is inside an outline OR within that outline's equivalent radius of it
                cR = pts[rid == R].mean(0); c2 = ((cR - o) @ u, (cR - o) @ v)
                for l in p['loops']:
                    q2 = np.c_[(l[0] - o) @ u, (l[0] - o) @ v]
                    if inside2d(c2, q2) or np.linalg.norm(cR - l[0][:-1].mean(0)) < max(2.0, np.sqrt(l[2] / np.pi)):
                        pieces.append(R); break
            if not pieces:                                      # fall back: the piece nearest the CSA centroid
                nomatch += 1; pieces = [min(np.unique(rid), key=lambda R: np.linalg.norm(pts[rid == R].mean(0) - p['centroid']))]
            A = Q = ef = iTP = iTPvn = 0.0
            calc.Function = 'TotalPressure*(Velocity_X*(%.12g)+Velocity_Y*(%.12g)+Velocity_Z*(%.12g))' % tuple(n)
            for R in pieces:
                thr.LowerThreshold = float(R); thr.UpperThreshold = float(R); iv.UpdatePipeline(ct); g = dsa.WrapDataObject(sm.Fetch(iv))
                Ai = float(np.asarray(g.CellData['Area'])[0])
                if Ai <= 0: continue
                A += Ai; Q += float(np.asarray(g.PointData['Velocity'])[0] @ n); ef += float(np.asarray(g.PointData['EnergyFlux'])[0] @ n)
                iTP += float(np.asarray(g.PointData['TotalPressure'])[0]); iTPvn += float(np.asarray(g.PointData['TPvn'])[0])
            if A <= 0: continue
            fh.write('%s,%d,%d,%.4f,%s,%.6f,%.6f,%.6f,%.6f,%.8g,%.6f,%.6f,%d\n' % (fn, int(fn.split('_')[1]), k, ct, reg, p['arc'], p['area'],
                     iTP / A, iTPvn / Q if abs(Q) > 1e-12 else iTP / A, ef, Q * 6e4, A * 1e6, len(pieces)))
    fh.flush()
    if i % 10 == 0: print('CFD: frame %d/%d  %.0fs' % (i, len(frames), time.time() - t0), flush=True)
print('CFD planes with no CFD piece inside the CSA outline (nearest used): %d' % nomatch, flush=True)
print('DONE', flush=True)
