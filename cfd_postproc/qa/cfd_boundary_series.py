"""Per-step boundary series from a STAR-CCM+ EnSight export (run with pvbatch).

For each step: volumetric flow through BronchiOut and MaskInlet (cap normals oriented OUT of the
domain, so + = air leaving), cap areas, and area-averaged static/total pressure on both caps.
Rows are written and flushed one step at a time (crash-safe; resumes by skipping done steps).

Usage: pvb.sh cfd_boundary_series.py --case X.case --out series.csv [--steps 0:40 | --steps 3,7,9]
"""
import argparse, os
import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

ap = argparse.ArgumentParser()
ap.add_argument('--case', required=True); ap.add_argument('--out', required=True); ap.add_argument('--steps', default=None)
a = ap.parse_args()
r = OpenDataFile(a.case); r.PointArrays = ['Velocity', 'Pressure', 'TotalPressure']; tv = list(r.TimestepValues)
if a.steps is None: steps = list(range(len(tv)))
elif ':' in a.steps: lo, hi = a.steps.split(':'); steps = list(range(int(lo), min(int(hi), len(tv))))
else: steps = [int(s) for s in a.steps.split(',')]

def blocks(mb):
    return {mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()).replace('Block 0_', '').split('.')[-1].strip(): mb.GetBlock(i)
            for i in range(mb.GetNumberOfBlocks())}

def cap(blk, dom_c):
    """(flow L/min out of domain, area mm2, area-avg p, area-avg TP) for one boundary cap."""
    X = vtk_to_numpy(blk.GetPoints().GetData()); pd = blk.GetPointData()
    V = vtk_to_numpy(pd.GetArray('Velocity')); P = vtk_to_numpy(pd.GetArray('Pressure')); T = vtk_to_numpy(pd.GetArray('TotalPressure'))
    flux = 0.0; Avec = np.zeros(3); A = []; pc = []; tc = []
    for c in range(blk.GetNumberOfCells()):
        cell = blk.GetCell(c); ids = [cell.GetPointId(k) for k in range(cell.GetNumberOfPoints())]; p = X[ids]; av = np.zeros(3)
        for k in range(len(ids)):
            u, w = p[k], p[(k + 1) % len(ids)]
            av += [(u[1] - w[1]) * (u[2] + w[2]), (u[2] - w[2]) * (u[0] + w[0]), (u[0] - w[0]) * (u[1] + w[1])]
        av *= 0.5; flux += V[ids].mean(0).dot(av); Avec += av; A.append(np.linalg.norm(av)); pc.append(P[ids].mean()); tc.append(T[ids].mean())
    A = np.array(A); sgn = np.sign(np.dot(Avec, X.mean(0) - dom_c))      # orient the area vector outward
    return flux * sgn * 6e4, A.sum() * 1e6, float((np.array(pc) * A).sum() / A.sum()), float((np.array(tc) * A).sum() / A.sum())

done = set()
if os.path.exists(a.out):
    done = {int(l.split(',')[0]) for l in open(a.out).read().splitlines()[1:] if l}
fh = open(a.out, 'a')
if not done: fh.write('step,label_t_s,Q_bronchi_out_Lmin,A_bronchi_mm2,p_bronchi_Pa,tp_bronchi_Pa,Q_mask_out_Lmin,A_mask_mm2,p_mask_Pa,tp_mask_Pa\n'); fh.flush()
for k in steps:
    if k in done: continue
    r.UpdatePipeline(tv[k]); B = blocks(sm.Fetch(r))
    dom_c = vtk_to_numpy(B['Region'].GetPoints().GetData()).mean(0)
    qb, ab, pb, tb = cap(B['BronchiOut'], dom_c); qm, am, pm, tm = cap(B['MaskInlet'], dom_c)
    fh.write('%d,%.6f,%.5f,%.5f,%.4f,%.4f,%.5f,%.5f,%.4f,%.4f\n' % (k, tv[k], qb, ab, pb, tb, qm, am, pm, tm)); fh.flush()
print('DONE', a.out, flush=True)
