"""Per-step velocity magnitude percentiles over the whole fluid domain (pvbatch). Picks the render colour scale.
Usage: pvb.sh velocity_peak_scan.py --case X.case --out scan.csv [--steps a:b]"""
import argparse, os
import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from vtkmodules.util.numpy_support import vtk_to_numpy
ap = argparse.ArgumentParser(); ap.add_argument('--case', required=True); ap.add_argument('--out', required=True); ap.add_argument('--steps', default=None)
a = ap.parse_args()
r = OpenDataFile(a.case); r.PointArrays = ['VelocityMagnitude']; tv = list(r.TimestepValues)
lo, hi = (int(x) for x in a.steps.split(':')) if a.steps else (0, len(tv))
fh = open(a.out, 'w'); fh.write('step,t_s,v_p99,v_p999,v_max,x_max_mm,y_max_mm,z_max_mm\n')
for k in range(lo, min(hi, len(tv))):
    r.UpdatePipeline(tv[k]); reg = sm.Fetch(r).GetBlock(0)
    v = vtk_to_numpy(reg.GetPointData().GetArray('VelocityMagnitude')); X = vtk_to_numpy(reg.GetPoints().GetData())
    i = int(np.argmax(v)); p99, p999 = np.percentile(v, [99, 99.9])
    fh.write('%d,%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%.2f\n' % (k, tv[k], p99, p999, v[i], *(X[i] * 1000))); fh.flush()
