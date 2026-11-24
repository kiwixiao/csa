#!/usr/bin/env python3
"""Shared loaders for the ENT003 CSA / pressure-energy post-processing.

ENT002 was a single-region case: one centerline, mouth at the high-arc end, so every
plot could just do ``s = arc.max() - arc`` on the whole CSV. ENT003's septum-aware
csa-slicer emits four regions instead:

    DescendingAirway   arc   0.0 -> 160.4 mm   trachea  -> nasopharynx
    LeftNose           arc 164.4 -> 221.9 mm   nasopharynx -> left nostril
    RightNose          arc 165.4 -> 220.0 mm   nasopharynx -> right nostril
    Mouth              arc   0.0 ->  73.5 mm   oropharynx -> lips  (separate branch)

The nose regions continue DescendingAirway's arc numbering, so a nostril-to-trachea
route is just the concatenation of the two. A *path* here is that ordered region list;
distance is measured from the airway opening exactly as ENT002 measured it from the
mouth, so the downstream plots keep the same shape and axis convention:

    s = 0   at the airway opening (nostril / lips)
    s = max at the distal end (trachea / oropharynx)

Sign convention for the breath, established by integrating the exported caps against
ENT003FlowProfile.csv (see postproc/outlet_flow/inlet_outlet_massflow.csv):
positive prescribed mass flow = EXHALATION, negative = INHALATION.
"""
import csv
import os

import numpy as np

# Ordered opening -> distal. Named paths the plot scripts accept for --path.
PATHS = {
    "LeftNose":         ["LeftNose", "DescendingAirway"],
    "RightNose":        ["RightNose", "DescendingAirway"],
    "Mouth":            ["Mouth"],
    "DescendingAirway": ["DescendingAirway"],
}
PATH_LABEL = {
    "LeftNose":         ("distance from left nostril (mm)", "left nostril", "trachea"),
    "RightNose":        ("distance from right nostril (mm)", "right nostril", "trachea"),
    "Mouth":            ("distance from lips (mm)", "lips", "oropharynx"),
    "DescendingAirway": ("distance from nasopharynx (mm)", "nasopharynx", "trachea"),
}
RHO = float(os.environ.get('CFD_RHO', '1.2367'))   # kg/m3: the density the solver ran (fitted from the outlet flux); set per subject by the driver


def xlabel(path):
    """Compact axis label. The long ENT002 form ('... - mouth 0 -> trachea') overflows
    the composite figure for ENT003's longer region names and gets clipped."""
    lab, _prox, dist = PATH_LABEL[path]
    return u"%s   —   0 → %s" % (lab, dist)


def regions_of(path):
    if path in PATHS:
        return PATHS[path]
    return [r.strip() for r in path.split(",") if r.strip()]


# --------------------------------------------------------------------------- #
# geometry (csa-bifurcation <subject>-Data.csv)
# --------------------------------------------------------------------------- #
def load_anatomy(csv_path, path):
    """-> (times_s, {t: dict(s, area, dh, pidx, cxyz)}) with s measured from the opening."""
    regs = set(regions_of(path))
    frames = {}
    for r in csv.DictReader(open(csv_path)):
        if r["is_valid"].lower() != "true" or r["region"] not in regs:
            continue
        t = int(r["frame_name"].split("_")[1]) / 1000.0
        frames.setdefault(t, []).append(r)
    times = sorted(frames)
    per = {}
    for t in times:
        rs = sorted(frames[t], key=lambda r: -float(r["arc_length_mm"]))   # opening first
        arc = np.array([float(r["arc_length_mm"]) for r in rs])
        per[t] = dict(
            s=arc.max() - arc,
            area=np.array([float(r["area_mm2"]) for r in rs]),
            dh=np.array([float(r["hydraulic_diameter_mm"]) for r in rs]),
            pidx=np.array([int(r["plane_index"]) for r in rs]),
            region=np.array([r["region"] for r in rs]),
            cxyz=np.array([[float(r["centroid_x"]), float(r["centroid_y"]),
                            float(r["centroid_z"])] for r in rs]))
    return times, per


# --------------------------------------------------------------------------- #
# CFD (csa_pressure_energy.py output, one csa_pe_<Region>.csv per region)
# --------------------------------------------------------------------------- #
def load_pe(pe_dir, path, qmin=0.5):
    """Concatenate the per-region csa_pe CSVs of *path*, keyed by CFD time.

    Rows are keyed on frame_ms (1:1 with the motion frame) and then collapsed to
    cfd_t_s, because out_000000/out_000025 both map to CFD index 0 and would
    otherwise double up. Returns (times, {t: dict(...)}) with s from the opening.
    tp_m is NaN where |Q| < qmin L/min, since a mass-flow average is meaningless
    at a flow reversal.
    """
    byms = {}
    for reg in regions_of(path):
        p = os.path.join(pe_dir, "csa_pe_%s.csv" % reg)
        if not os.path.exists(p):
            raise SystemExit("missing %s -- run csa_pressure_energy.py --region %s" % (p, reg))
        for d in csv.DictReader(open(p)):
            byms.setdefault(int(d["frame_ms"]), []).append(d)

    # frame_ms -> cfd time; keep the LAST frame_ms for each cfd index (out_000025 over
    # out_000000, out_002180 over out_002175) so every plotted frame is a distinct solve
    bycfd = {}
    for ms in sorted(byms):
        bycfd[float(byms[ms][0]["cfd_t_s"])] = byms[ms]

    times = sorted(bycfd)
    data = {}
    for t in times:
        rs = sorted(bycfd[t], key=lambda r: -float(r["arc_mm"]))          # opening first
        arc = np.array([float(r["arc_mm"]) for r in rs])
        Q = np.array([float(r["Q_Lmin"]) for r in rs])
        tp_m = np.array([float(r["TP_mdot_avg_Pa"]) for r in rs])
        tp_m[np.abs(Q) < qmin] = np.nan
        data[t] = dict(s=arc.max() - arc,
                       area=np.array([float(r["area_mm2"]) for r in rs]),
                       tp_a=np.array([float(r["TP_area_avg_Pa"]) for r in rs]),
                       tp_m=tp_m,
                       ef=np.array([float(r["energy_flux_W"]) for r in rs]),
                       Q=Q,
                       region=np.array([r["region"] for r in rs]))
    return times, data


def cum_from_opening(a):
    """Cumulative loss referenced to the airway opening (a[0]), NaN-tolerant."""
    fin = a[np.isfinite(a)]
    ref = a[0] if np.isfinite(a[0]) else (fin[0] if fin.size else np.nan)
    return ref - a


# --------------------------------------------------------------------------- #
# prescribed breath
# --------------------------------------------------------------------------- #
def load_profile(prof_csv):
    """-> (t, mdot_kg_s). Positive = exhalation (verified against both caps)."""
    T, M = [], []
    for r in csv.reader(open(prof_csv)):
        try:
            T.append(float(r[0])); M.append(float(r[1]))
        except (ValueError, IndexError):
            continue
    return np.array(T), np.array(M)


def phase_of(prof_csv, times, dead=0.3):
    """Per-time breath phase from the prescribed profile.

    -> (q_Lmin, is_inhale, is_exhale). q is signed so that INHALATION is negative,
    matching the profile's own convention. dead = L/min band treated as transition.
    """
    T, M = load_profile(prof_csv)
    q = np.interp(np.asarray(times, float), T, M) / RHO * 60000.0     # kg/s -> L/min
    return q, q < -dead, q > dead


def peak_frames(prof_csv, times):
    """Indices of peak inhalation and peak exhalation within *times*."""
    q, _, _ = phase_of(prof_csv, times)
    return int(np.argmin(q)), int(np.argmax(q))


# --------------------------------------------------------------------------- #

def side_centerline(subject_dir, path, frame, s_max=None):
    """The real per-side centerline of a nostril path for *frame* (motion/centerlines/<side>/),
    ordered from the OPENING, with its arc length from the opening. If *s_max* (the path's plane
    distance range) is given the arc is scaled to it, so markers line up with the plotted s.
    Returns (P (N,3) mm, s) or None when no per-side centerline exists for this path.
    Drawing this instead of the cross-section centroids avoids their wobble where the section
    changes shape (nasal slits, epiglottis); the distances themselves are unaffected."""
    side = {"LeftNose": "leftnose", "RightNose": "rightnose"}.get(path)
    if side is None:
        return None
    p = os.path.join(subject_dir, "motion", "centerlines", side, frame + ".vtk")
    if not os.path.exists(p):
        return None
    tok = open(p, errors="ignore").read().split()
    i = tok.index("POINTS"); n = int(tok[i + 1])
    P = np.array(tok[i + 3:i + 3 + 3 * n], float).reshape(n, 3)
    if P[0, 2] < P[-1, 2]:                      # nostril end (high z) first
        P = P[::-1]
    sc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    if s_max:
        sc *= s_max / sc[-1]
    return P, sc

def stl_vertices(path):
    """(N,3) vertex array from a binary or ASCII STL."""
    sz = os.path.getsize(path)
    with open(path, "rb") as f:
        f.seek(80)
        n = int(np.frombuffer(f.read(4), dtype="<u4")[0])
        if 84 + n * 50 == sz:
            dt = np.dtype([("nrm", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")])
            return np.frombuffer(f.read(n * 50), dtype=dt, count=n)["v"].reshape(-1, 3).astype(float)
    return np.array([[float(x) for x in ln.split()[1:4]]
                     for ln in open(path) if ln.strip().startswith("vertex")])
