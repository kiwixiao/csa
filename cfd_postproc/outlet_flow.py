#!/usr/bin/env pvbatch
"""Integrate volumetric/mass flow through a named boundary cap at every CFD frame,
so the true time of each (mislabeled) frame can be recovered by matching to the
subject's flow profile. Masks garbage cells (|v|>100 m/s) so corrupt frames don't
poison the integral. Writes one CSV row per frame (incremental / crash-safe).

Ported from /storage4/ENT002/outlet_flow.py. Differences:
  * --case / --out / --part are arguments instead of module constants.
  * ENT003's export omits the BronchiOut cap and instead carries a Nose part, so
    the flux is taken through MaskInlet. By mass conservation that equals the
    tracheal outflow with opposite sign; --flip normalises the sign convention.
  * If the requested part is absent the available parts are listed and the run
    aborts, rather than silently writing empty rows.

Usage:
  pvbatch outlet_flow.py --case /storage4/ENT003/case/ENT003.case \
      --part MaskInlet --flip --out /storage4/ENT003/postproc/outlet_flow/outlet_Q.csv
"""
import argparse
import os
import sys
import time

import numpy as np
from paraview.simple import *  # noqa
from paraview import servermanager as sm
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

RHO = 1.2        # kg/m3, for mass-flow comparison to the profile
VGARB = 100.0    # m/s: above this a cell is corrupt -> excluded


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", required=True, help="EnSight .case file")
    p.add_argument("--out", required=True, help="output CSV path")
    p.add_argument("--part", default="BronchiOut",
                   help="substring of the boundary part to integrate (default: BronchiOut)")
    p.add_argument("--flip", action="store_true",
                   help="negate the flux (use when integrating an inlet instead of the outlet)")
    p.add_argument("--rho", type=float, default=RHO, help="density kg/m3")
    p.add_argument("--vgarb", type=float, default=VGARB,
                   help="reject cells with |v| above this (m/s)")
    return p.parse_args()


def block_names(mb):
    """Return [(index, name)] for every block in the multiblock dataset."""
    out = []
    for i in range(mb.GetNumberOfBlocks()):
        md = mb.GetMetaData(i)
        out.append((i, (md.Get(vtkMultiBlockDataSet.NAME()) if md else "").strip()))
    return out


def find_block(mb, needle):
    """Return the first block whose name contains *needle* (case-insensitive)."""
    for i, name in block_names(mb):
        if name and needle.lower() in name.lower():
            return mb.GetBlock(i)
    return None


def flux_through(blk, vgarb):
    """return (flux_all, flux_clean, area, nbad, ncell, vmax) via per-cell area-vector."""
    P = vtk_to_numpy(blk.GetPoints().GetData())
    V = vtk_to_numpy(blk.GetPointData().GetArray('Velocity'))
    vmax = float(np.abs(V).max())
    nc = blk.GetNumberOfCells()
    fa = fc = area = 0.0
    nbad = 0
    for c in range(nc):
        cell = blk.GetCell(c)
        npn = cell.GetNumberOfPoints()
        ids = [cell.GetPointId(k) for k in range(npn)]
        pts = P[ids]
        av = np.zeros(3)                      # Newell area-vector (|av| = area, dir = normal)
        for k in range(npn):
            a = pts[k]; b = pts[(k + 1) % npn]
            av[0] += (a[1] - b[1]) * (a[2] + b[2])
            av[1] += (a[2] - b[2]) * (a[0] + b[0])
            av[2] += (a[0] - b[0]) * (a[1] + b[1])
        av *= 0.5
        vc = V[ids].mean(axis=0)
        f = float(vc.dot(av))
        fa += f
        area += float(np.linalg.norm(av))
        if np.abs(V[ids]).max() > vgarb:
            nbad += 1
        else:
            fc += f
    return fa, fc, area, nbad, nc, vmax


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    r = OpenDataFile(args.case)
    r.PointArrays = ['Velocity']
    tv = list(r.TimestepValues)
    print("case:", args.case, flush=True)
    print("frames:", len(tv), "label range: %.4f -> %.4f s" % (tv[0], tv[-1]), flush=True)

    # Fail loudly if the requested part is not in the export.
    r.UpdatePipeline(tv[0])
    mb0 = sm.Fetch(r)
    if find_block(mb0, args.part) is None:
        print("ERROR: no block matching %r. Available parts:" % args.part, flush=True)
        for i, name in block_names(mb0):
            print("   [%d] %s" % (i, name), flush=True)
        sys.exit(1)
    print("integrating through part matching %r (flip=%s)" % (args.part, args.flip), flush=True)

    sign = -1.0 if args.flip else 1.0
    with open(args.out, "w") as out:
        out.write("idx,label_t_s,flux_all_m3s,flux_clean_m3s,Q_clean_Lmin,"
                  "mdot_clean_kgs,area_m2,nbad,ncell,vmax\n")
        t0 = time.time()
        for idx in range(len(tv)):
            r.UpdatePipeline(tv[idx])
            blk = find_block(sm.Fetch(r), args.part)
            if blk is None:
                out.write("%d,%.4f,,,,,,,,\n" % (idx, tv[idx])); out.flush(); continue
            fa, fc, area, nbad, nc, vmax = flux_through(blk, args.vgarb)
            fa *= sign
            fc *= sign
            out.write("%d,%.4f,%.6e,%.6e,%.4f,%.6e,%.6e,%d,%d,%.4g\n"
                      % (idx, tv[idx], fa, fc, fc * 60000.0, fc * args.rho,
                         area, nbad, nc, vmax))
            out.flush()
            if idx % 5 == 0 or nbad:
                print("idx=%3d lbl=%7.3fs Qclean=%8.2f L/min bad=%d/%d vmax=%.3g  (%.0fs)"
                      % (idx, tv[idx], fc * 60000.0, nbad, nc, vmax, time.time() - t0),
                      flush=True)
    print("DONE ->", args.out, flush=True)


main()
