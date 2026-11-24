"""Total pressure / flow / energy flux on the nostril OPENING plane (the CSA surface's cut plane).

The CSA slicer's first planes cross the cap opening, so their outline is not closed there;
the nostril end of the path then has no reliable plane for the first few mm. This slices the
CFD directly at the nostril cut plane every step and integrates over each nostril's piece of
the slice, exactly like csa_pressure_energy.py (IntegrateVariables over the region nearest each
cap centre; normal pointing OUT of the nose, i.e. toward increasing centerline arc, so the sign
convention matches the csa_pe tables: Q < 0 on inhale).

Usage: pvb.sh opening_slice.py --case X.case --cuts cut_planes.json --out opening.csv
"""
import argparse, json, time
import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from paraview.vtk.numpy_interface import dataset_adapter as dsa

ap = argparse.ArgumentParser(); ap.add_argument('--case', required=True); ap.add_argument('--cuts', required=True)
ap.add_argument('--out', required=True); a = ap.parse_args()
J = json.load(open(a.cuts))
o = np.array(J['nostril_plane']['origin']) * 1e-3
n = -np.array(J['nostril_plane']['inward_normal'], float); n /= np.linalg.norm(n)      # out of the nose
caps = {'LeftNose': np.array(J['caps']['Leftnose']['centre']) * 1e-3, 'RightNose': np.array(J['caps']['Rightnose']['centre']) * 1e-3}

r = OpenDataFile(a.case); r.PointArrays = ['TotalPressure', 'EnergyFlux', 'Velocity']; tv = list(r.TimestepValues)
mb = MergeBlocks(Input=r)
sl = Slice(Input=mb); sl.SliceType = 'Plane'; sl.SliceType.Origin = list(o); sl.SliceType.Normal = list(n)
conn = Connectivity(Input=sl)
thr = Threshold(Input=conn, Scalars=['POINTS', 'RegionId'])
calc = Calculator(Input=thr); calc.AttributeType = 'Point Data'; calc.ResultArrayName = 'TPvn'
calc.Function = 'TotalPressure*(Velocity_X*(%.12g)+Velocity_Y*(%.12g)+Velocity_Z*(%.12g))' % tuple(n)
iv = IntegrateVariables(Input=calc)
fh = open(a.out, 'w'); fh.write('cfd_idx,cfd_t_s,region,area_mm2,TP_area_avg_Pa,TP_mdot_avg_Pa,energy_flux_W,Q_Lmin,cx_mm,cy_mm,cz_mm,dist_to_cap_mm\n')
t0 = time.time()
for k, t in enumerate(tv):
    conn.UpdatePipeline(t)
    d = dsa.WrapDataObject(sm.Fetch(conn)); pts = np.asarray(d.Points); rid = np.asarray(d.PointData['RegionId'])
    cents = {R: pts[rid == R].mean(0) for R in np.unique(rid)}
    for reg, cc in caps.items():
        R = min(cents, key=lambda q: np.linalg.norm(cents[q] - cc))
        thr.LowerThreshold = float(R); thr.UpperThreshold = float(R); iv.UpdatePipeline(t)
        g = dsa.WrapDataObject(sm.Fetch(iv)); A = float(np.asarray(g.CellData['Area'])[0])
        Q = float(np.asarray(g.PointData['Velocity'])[0] @ n); ef = float(np.asarray(g.PointData['EnergyFlux'])[0] @ n)
        iTP = float(np.asarray(g.PointData['TotalPressure'])[0]); iTPvn = float(np.asarray(g.PointData['TPvn'])[0])
        tp_a = iTP / A; tp_m = iTPvn / Q if abs(Q) > 1e-12 else tp_a
        fh.write('%d,%.4f,%s,%.4f,%.5f,%.5f,%.8g,%.5f,%.3f,%.3f,%.3f,%.2f\n' % (k, t, reg, A * 1e6, tp_a, tp_m, ef, Q * 6e4,
                 *(cents[R] * 1e3), np.linalg.norm(cents[R] - cc) * 1e3))
    fh.flush()
    if k % 10 == 0: print('step %d/%d  %.0fs' % (k, len(tv), time.time() - t0), flush=True)
print('DONE', a.out, flush=True)
