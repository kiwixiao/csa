#!/usr/bin/env python3
"""Registration control-point trajectories for the CFD mesh-motion check (qa/wall_track.py).

The registration CSV (e.g. ENT5_noalign_be0.001.csv) holds the morphing control points, one row per point and
x,y,z per millisecond as columns: these are exactly the positions the CFD motion tables prescribe.
Writes them as a (points, ms, 3) float32 array. --tables-dir also checks a few of the CFD motion tables
(table_<k>_*.csv: x,y,z,dx,dy,dz) against the CSV.
Usage: python3 reg_traj.py --reg-csv REG.csv --out reg_traj.npy [--tables-dir CFD/Tables]
"""
import argparse, glob, os
import numpy as np, pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument("--reg-csv", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--tables-dir", default=None)
a = ap.parse_args()
df = pd.read_csv(a.reg_csv)
A = df.to_numpy().reshape(len(df), -1, 3)
print("trajectories: %d points x %d ms (columns %s .. %s)" % (A.shape[0], A.shape[1], df.columns[0], df.columns[-1]))
np.save(a.out, A.astype(np.float32))
if a.tables_dir:
    for k in sorted({0, 1, A.shape[1] // 2, A.shape[1] - 2}):
        f = glob.glob("%s/table_%05d_*.csv" % (a.tables_dir, k))
        if not f:
            print("table %d: not found" % k); continue
        t = np.loadtxt(f[0], delimiter=",", skiprows=1)
        print("table %5d %s: position diff %.3g mm, increment diff %.3g mm" % (
            k, os.path.basename(f[0]), np.abs(t[:, :3] - A[:, k]).max(), np.abs(t[:, 3:] - (A[:, k + 1] - A[:, k])).max()))
d = np.linalg.norm(A[:, -1] - A[:, 0], axis=1)
print("end - start displacement: mean %.2f, max %.2f mm (cycle closure)" % (d.mean(), d.max()))
