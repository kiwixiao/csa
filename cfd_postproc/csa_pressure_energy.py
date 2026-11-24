#!/usr/bin/env pvbatch
# -*- coding: utf-8 -*-
"""CSA-plane total-pressure + energy-flux vs arc-length, per matched CFD time.

For each CSA frame (out_XXXXXX ms) -> matched CFD step; for each plane: slice the CFD
volume (origin=CSA centroid, normal=centerline tangent), isolate the lumen cross-section
nearest the centroid, and compute area-averaged & mass-flow-averaged TotalPressure, the
energy flux n.integral(EnergyFlux) and volume flux Q=n.integral(Velocity).

Ported from /storage4/ENT002/csa_pressure_energy.py. Differences:
  * --case / --csa / --out are arguments instead of module constants.
  * --region selects which CSA region to walk. ENT002 was single-region (Mouth);
    ENT003's septum path yields DescendingAirway / LeftNose / RightNose / Mouth,
    so rows are filtered by region and the plane-STL tag is looked up per region.
  * every CFD quantity is an AREA-WEIGHTED SURFACE INTEGRAL (ParaView
    IntegrateVariables), not an unweighted mean over the slice nodes. Q = int v.n dA,
    energy = int E.n dA, TP_mdot = int TP(v.n)dA / int(v.n)dA. The nodal mean is a
    quadrature error on a boundary-layer mesh -- nodes crowd the wall where v->0, so it
    read ~40% low and, fatally, did not conserve flux along the trunk (p90/p10 = 1.7x
    where the physics says flat).
  * all disconnected lumen pieces near the plane centroid are summed, so a multi-lobed
    cross-section (the airway is 2-3 separate channels around the epiglottis) is not
    reduced to whichever lobe happens to sit nearest the centreline.
  * --rank-csa names the canonical CSA table that defines the STL numbering. Needed when
    --csa is a filtered/smoothed table: dropping a plane shifts every later rank, which
    would silently fetch a neighbouring plane's disk.
  * the plane STL is addressed by the plane's RANK within its region, not by its
    plane_index. ENT002's plane_index happened to equal the file number; ENT003's
    does not (LeftNose is 16..73 over files 000..055, Mouth is 128..201 over
    000..074), so indexing by plane_index silently fetched a plane ~15 mm away --
    or nothing at all, falling back to the coarse central-difference normal.
"""
import builtins as _bi
import argparse, csv, os, re, sys, time
from paraview.simple import *  # noqa
from paraview import servermanager as sm
from paraview.vtk.numpy_interface import dataset_adapter as dsa
import numpy as np

CORRUPT_TP = 1.0e5            # |TotalPressure| above this (Pa) => corrupted step, skip
CFD_DT = 0.025                # s per CFD save; CFD idx = round(frame_ms/25)-1

# csa-slicer writes plane STLs as <frame>-<TAG>-<idx>.stl under <frame>/<Region>_Planes/
REGION_TAG = {
    "DescendingAirway": "Desc",
    "LeftNose": "Left",
    "RightNose": "Right",
    "Mouth": "Mouth",
}


def stl_normal(path):
    """(normal, centroid) of a plane-STL disk via PCA of its vertices. None if unreadable.
    This is the ACTUAL csa-slicer plane orientation (smooth centerline tangent) — far cleaner at
    bends than central-differencing the wobbly cross-section centroids."""
    try:
        sz = os.path.getsize(path)
        with open(path, "rb") as f:
            f.seek(80); n = int(np.frombuffer(f.read(4), dtype="<u4")[0])
            if 84 + n * 50 == sz:                      # binary STL
                dt = np.dtype([("nrm", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")])
                V = np.frombuffer(f.read(n * 50), dtype=dt, count=n)["v"].reshape(-1, 3).astype(float)
            else:                                      # ASCII STL
                V = np.array([[float(x) for x in ln.split()[1:4]]
                              for ln in open(path) if ln.strip().startswith("vertex")])
        if len(V) < 3:
            return None
        c = V.mean(axis=0)
        return np.linalg.svd(V - c, full_matrices=False)[2][2], c
    except Exception:
        return None


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--case", required=True, help="EnSight .case file")
    p.add_argument("--csa", required=True, help="<subject>-Data.csv from csa-bifurcation")
    p.add_argument("--region", default="Mouth", choices=sorted(REGION_TAG),
                   help="CSA region to walk (default: Mouth)")
    p.add_argument("--velocity-array", default="Velocity",
                   help="name of the velocity vector in the export (OpenFOAM writes 'U')")
    p.add_argument("--energy-array", default="EnergyFlux",
                   help="name of the energy-flux vector, or 'none' if the export has "
                        "no such field (OpenFOAM). The WOB panels use dP*Q and do not "
                        "need it. NB: pass 'none', not an empty string -- pvbatch cannot "
                        "decode an empty command-line argument.")
    p.add_argument("--cfd-index-offset", type=int, default=-1,
                   help="cfd_idx = round(frame_ms/25) + offset. -1 for the STAR export "
                        "whose first step is t=0.025 s; 0 for OpenFOAM's, which starts at 0.")
    p.add_argument("--planes-dir", default=None,
                   help="folder holding <frame>/<Region>_Planes/*.stl (plane normals); default: <dir of --rank-csa or --csa>/frames")
    p.add_argument("--rank-csa", default=None,
                   help="canonical CSA table defining the plane-STL numbering "
                        "(default: --csa). Pass the RAW table when --csa is filtered.")
    p.add_argument("--lobe-radius-factor", type=float, default=2.0,
                   help="sum lumen pieces whose centroid lies within this factor times "
                        "the plane's equivalent radius (never less than --lobe-radius-min)")
    p.add_argument("--lobe-radius-min", type=float, default=6.0, help="mm")
    p.add_argument("--frame-stride", type=int, default=4, help="use every Nth CSA frame")
    p.add_argument("--frames", default=None,
                   help="a:b - only frames a..b-1 of the (strided) frame list, to split a region over workers")
    p.add_argument("--plane-stride", type=int, default=1)
    p.add_argument("--out", required=True)
    return p.parse_args(argv)


def rank_map(csa_path, region):
    """{plane_index: rank} -- the order the STL disks were written in. The slicer writes a disk
    for EVERY plane, invalid ones included, so rank over <S>-Data-all.csv when it exists (ranking
    over the valid rows only shifted every disk after an invalid plane by one: ENT005 40% of
    planes took a neighbour's normal, median 1.1 deg off). Computed per frame."""
    all_path = csa_path.replace("-Data.csv", "-Data-all.csv")
    use_all = all_path != csa_path and os.path.exists(all_path)
    per = {}
    with open(all_path if use_all else csa_path) as f:
        for d in csv.DictReader(f):
            if d.get("region") != region or (not use_all and d["is_valid"].lower() != "true"):
                continue
            per.setdefault(d["frame_name"], []).append(int(d["plane_index"]))
    return {fn: {p: k for k, p in enumerate(sorted(v))} for fn, v in per.items()}


def read_planes(csa_path, region):
    """{frame_name: [(arc_mm, area_mm2, centroid_xyz_mm, plane_index), ...]} sorted by arc."""
    frames = {}
    with open(csa_path) as f:
        for d in csv.DictReader(f):
            if d["is_valid"].lower() != "true":
                continue
            if d.get("region") != region:
                continue
            fn = d["frame_name"]
            frames.setdefault(fn, []).append(
                [float(d["arc_length_mm"]), float(d["area_mm2"]),
                 np.array([float(d["centroid_x"]), float(d["centroid_y"]), float(d["centroid_z"])]),
                 int(d["plane_index"])])
    for fn in frames:
        frames[fn].sort(key=lambda x: x[0])
    return frames


def main():
    a = parse_args(sys.argv[1:])
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    planes_dir = a.planes_dir or os.path.join(os.path.dirname(a.rank_csa or a.csa), "frames")
    if not os.path.isdir(planes_dir):
        print("ERROR: plane-STL folder not found: %s (pass --planes-dir)" % planes_dir, flush=True); sys.exit(1)
    print("plane STLs (for plane normals): %s" % planes_dir, flush=True)
    tag = REGION_TAG[a.region]

    frames = read_planes(a.csa, a.region)
    if not frames:
        print("ERROR: no valid rows for region %r in %s" % (a.region, a.csa), flush=True)
        sys.exit(1)
    fnames = sorted(frames, key=lambda s: int(re.search(r"(\d+)", s).group(1)))
    fnames = fnames[::a.frame_stride]
    if a.frames:
        lo, hi = a.frames.split(":")
        fnames = fnames[int(lo or 0):int(hi) if hi else None]
    print("region %s: %d CSA frames to process (stride %d), ~%d planes/frame"
          % (a.region, len(fnames), a.frame_stride, len(frames[fnames[0]])), flush=True)

    VEL = a.velocity_array
    EFA = "" if a.energy_array.lower() in ("", "none") else a.energy_array
    r = OpenDataFile(a.case)
    r.PointArrays = [x for x in ['TotalPressure', EFA, VEL] if x]
    mb = MergeBlocks(Input=r)
    tv = list(r.TimestepValues)
    sl = Slice(Input=mb); sl.SliceType = 'Plane'
    conn = Connectivity(Input=sl)
    thr = Threshold(Input=conn, Scalars=['POINTS', 'RegionId'])
    calc = Calculator(Input=thr)
    calc.AttributeType = 'Point Data'
    calc.ResultArrayName = 'TPvn'
    iv = IntegrateVariables(Input=calc)

    RANK = rank_map(a.rank_csa or a.csa, a.region)
    bad_rank = [0]

    HDR = ("frame,frame_ms,cfd_idx,cfd_t_s,region,arc_mm,area_mm2,"
           "TP_area_avg_Pa,TP_mdot_avg_Pa,energy_flux_W,Q_Lmin,area_int_mm2,n_lobes\n")
    done_frames = set()
    if os.path.exists(a.out) and os.path.getsize(a.out) > 0:      # resume: skip frames already written
        with open(a.out) as f:
            next(f, None)
            for line in f:
                done_frames.add(line.split(",", 1)[0])
        print("resume: %d frames already in CSV, skipping" % len(done_frames), flush=True)
        out = open(a.out, "a")
    else:
        out = open(a.out, "w"); out.write(HDR); out.flush()
    nrows = 0
    skipped = []
    t_start = time.time()
    for fn in fnames:
        if fn in done_frames:
            continue
        fms = int(re.search(r"(\d+)", fn).group(1))
        cfd_idx = int(round(fms / (CFD_DT * 1000.0))) + a.cfd_index_offset
        cfd_idx = _bi.max(0, _bi.min(len(tv) - 1, cfd_idx))
        ct = tv[cfd_idx]
        mb.UpdatePipeline(ct)
        ai = mb.GetDataInformation().GetPointDataInformation().GetArrayInformation('TotalPressure')
        tpmax = _bi.max(abs(ai.GetComponentRange(0)[0]), abs(ai.GetComponentRange(0)[1])) if ai else 0
        if tpmax > CORRUPT_TP:
            skipped.append((fn, cfd_idx, tpmax))
            print("  SKIP %s (cfd idx %d, |TP|max=%.3g > %.0g corrupt)" % (fn, cfd_idx, tpmax, CORRUPT_TP), flush=True)
            continue
        planes = frames[fn]
        cen = np.array([p[2] for p in planes])
        nrm = np.zeros_like(cen)                         # central-diff normal = SIGN reference only
        nrm[1:-1] = cen[2:] - cen[:-2]; nrm[0] = cen[1] - cen[0]; nrm[-1] = cen[-1] - cen[-2]
        nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
        frow = []
        for i in range(0, len(planes), a.plane_stride):
            arc, area_mm2, c_mm, pidx = planes[i]
            o = c_mm * 0.001
            prank = RANK.get(fn, {}).get(pidx)
            sn = None
            if prank is not None:
                got = stl_normal(os.path.join(planes_dir, fn, "%s_Planes" % a.region,
                                              "%s-%s-%03d.stl" % (fn, tag, prank)))
                if got is not None:
                    sn, scen = got
                    # free self-check: the disk we fetched must sit on this plane
                    if np.linalg.norm(scen - c_mm) > 2.0:
                        bad_rank[0] += 1
                        sn = None
            n = nrm[i] if sn is None else (-sn if np.dot(sn, nrm[i]) < 0 else sn)
            sl.SliceType.Origin = list(o); sl.SliceType.Normal = list(n)
            conn.UpdatePipeline(ct)
            d = dsa.WrapDataObject(sm.Fetch(conn))
            pts = np.asarray(d.Points)
            if pts.size == 0:
                continue
            rid = np.asarray(d.PointData['RegionId'])

            # Sum every lumen piece belonging to this cross-section. Around the
            # epiglottis the airway is genuinely 2-3 separate channels; taking only the
            # nearest one would drop part of the flux and break conservation.
            r_eq = np.sqrt(max(area_mm2, 1.0) / np.pi)
            rad = _bi.max(a.lobe_radius_min, a.lobe_radius_factor * r_eq) * 1e-3
            calc.Function = ('TotalPressure*(%s_X*(%.12g)+%s_Y*(%.12g)+%s_Z*(%.12g))'
                             % (VEL, n[0], VEL, n[1], VEL, n[2]))
            A_tot = 0.0; Q = 0.0; eflux = 0.0; iTP = 0.0; iTPvn = 0.0; npc = 0
            for R in np.unique(rid):
                if np.linalg.norm(pts[rid == R].mean(axis=0) - o) > rad:
                    continue
                thr.LowerThreshold = float(R); thr.UpperThreshold = float(R)
                iv.UpdatePipeline(ct)
                g = dsa.WrapDataObject(sm.Fetch(iv))
                A = float(np.asarray(g.CellData['Area'])[0])
                if A <= 0:
                    continue
                # IntegrateVariables returns the INTEGRAL of each point array over the
                # surface, i.e. exactly int f dA -- the area-weighted quantity we want.
                Q += float(np.asarray(g.PointData[VEL])[0] @ n)
                if EFA:
                    eflux += float(np.asarray(g.PointData[EFA])[0] @ n)
                iTP += float(np.asarray(g.PointData['TotalPressure'])[0])
                iTPvn += float(np.asarray(g.PointData['TPvn'])[0])
                A_tot += A; npc += 1
            if A_tot <= 0:
                continue
            tp_area = iTP / A_tot                                  # int TP dA / A
            tp_mdot = (iTPvn / Q) if abs(Q) > 1e-12 else tp_area   # int TP(v.n)dA / int(v.n)dA
            frow.append((fn, fms, cfd_idx, round(ct, 4), a.region, arc, area_mm2,
                         tp_area, tp_mdot, eflux, Q * 60000.0, A_tot * 1e6, npc))
        for row in frow:                      # write this frame's rows immediately (crash-safe)
            out.write(",".join(str(x) for x in row) + "\n")
        out.flush(); nrows += len(frow)
        el = time.time() - t_start
        print("  %s (t=%.3fs cfd%d): %d planes done, %.0fs elapsed" %
              (fn, ct, cfd_idx, len(planes), el), flush=True)

    out.close()
    if bad_rank[0]:
        print("WARNING: %d planes had an STL disk further than 2 mm from the CSA centroid "
              "(rank/plane_index mismatch?) and fell back to the central-difference normal"
              % bad_rank[0], flush=True)
    print("wrote %d rows to %s; skipped %d corrupt frames" % (nrows, a.out, len(skipped)), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
