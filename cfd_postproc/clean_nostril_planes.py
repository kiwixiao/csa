#!/usr/bin/env python3
"""Clean the nostril end of the CSA / csa_pe tables and add the nostril OPENING as s = 0.

1. Off-centerline planes: near the nostril the slicer's plane crosses the cap opening, its outline
   is not closed there, and it falls back to an unrelated closed loop (ENT005: a pocket ~20 mm
   above the vestibule, 8-24 mm2). A plane whose area centroid is more than --ratio equivalent
   diameters from its own side's centerline is flagged invalid (ENT005: pocket planes 3.9-6.6,
   every genuine plane <= 2.4 - the nasal cavity is a tall slit with the centerline on its floor).
2. Opening row: one row per frame per nasal region at the nostril cut plane, from
   opening_slice.py (CFD sliced exactly there). Its arc = the region's outermost plane arc +
   --delta-mm (cap -> first plane spacing), so csa_paths' s = arc.max() - arc is 0 at the cap.
Writes <out-dir>/<S>-Data-clean.csv and <out-dir>/csa_pe/csa_pe_<Region>.csv (originals untouched).
"""
import sys
import argparse, csv, glob, os
import numpy as np
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument('--subject-dir', required=True); ap.add_argument('--csa', required=True); ap.add_argument('--pe-dir', required=True)
ap.add_argument('--opening', default=None, help='opening_slice.py output (single opening row per frame)')
ap.add_argument('--fill-csa', default=None); ap.add_argument('--fill-pe', default=None,
               help='fill_nostril_planes.py output: cap + transition planes (replaces --opening)')
ap.add_argument('--out-dir', required=True)
ap.add_argument('--ratio', type=float, default=3.0); ap.add_argument('--delta-mm', type=float, default=1.3)
ap.add_argument('--cfd-steps', type=int, required=True, help='number of CFD steps (for frame -> cfd index clamping)')
a = ap.parse_args()
SIDE = {'LeftNose': 'leftnose', 'RightNose': 'rightnose'}

def read_vtk_points(path):
    """Points of a legacy ASCII/binary-free VTK polydata written by the slicer (ASCII POINTS block)."""
    txt = open(path, errors='ignore').read().split()
    i = txt.index('POINTS'); n = int(txt[i + 1])
    return np.array(txt[i + 3:i + 3 + 3 * n], float).reshape(n, 3)

rows = list(csv.DictReader(open(a.csa))); cols = list(rows[0].keys())
bad = set(); trees = {}
for r in rows:
    if r['region'] not in SIDE: continue
    key = (SIDE[r['region']], r['frame_name'])
    if key not in trees:
        trees[key] = cKDTree(read_vtk_points(os.path.join(a.subject_dir, 'motion', 'centerlines', key[0], r['frame_name'] + '.vtk')))
    c = np.array([float(r['centroid_x']), float(r['centroid_y']), float(r['centroid_z'])])
    if trees[key].query(c)[0] / max(float(r['equivalent_diameter_mm']), 1e-6) > a.ratio:
        bad.add((r['frame_name'], r['region'], r['arc_length_mm']))
print('off-centerline planes flagged invalid: %d (%s)' % (len(bad), {g: sum(1 for b in bad if b[1] == g) for g in SIDE}))
with open(os.path.join(a.out_dir, 'rejected_planes.csv'), 'w') as fh:
    fh.write('frame_name,region,plane_index,arc_length_mm\n')
    for r in rows:
        if (r['frame_name'], r['region'], r['arc_length_mm']) in bad: fh.write('%s,%s,%s,%s\n' % (r['frame_name'], r['region'], r['plane_index'], r['arc_length_mm']))

if not (a.fill_csa and a.fill_pe) and not a.opening:     # first pass: only the rejected list (fill_nostril_planes.py needs it)
    print('no --opening / --fill-csa: wrote rejected_planes.csv only'); sys.exit(0)

# opening arc per frame/region = outermost plane arc (valid or not) + delta
arc_open = {}
for r in rows:
    if r['region'] in SIDE:
        k = (r['frame_name'], r['region']); arc_open[k] = max(arc_open.get(k, -1e9), float(r['arc_length_mm']) + a.delta_mm)
FILL = bool(a.fill_csa and a.fill_pe)
O = {} if FILL else {(int(r['cfd_idx']), r['region']): r for r in csv.DictReader(open(a.opening))}
FC = list(csv.DictReader(open(a.fill_csa))) if FILL else []
FP = list(csv.DictReader(open(a.fill_pe))) if FILL else []
frames = sorted({r['frame_name'] for r in rows}, key=lambda f: int(f.split('_')[1]))
cfd_of = lambda fn: min(max(int(round(int(fn.split('_')[1]) / 25.0)) - 1, 0), a.cfd_steps - 1)

os.makedirs(os.path.join(a.out_dir, 'csa_pe'), exist_ok=True)
subj = os.path.basename(a.csa).split('-Data')[0]
with open(os.path.join(a.out_dir, '%s-Data-clean.csv' % subj), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
    for r in rows:
        if (r['frame_name'], r['region'], r['arc_length_mm']) in bad: r = dict(r, is_valid='False')
        w.writerow(r)
    if FILL:
        for r in FC:
            w.writerow({c: r.get(c, '') for c in cols})
    else:
        for fn in frames:
            for reg in SIDE:
                o = O[(cfd_of(fn), reg)]; area = float(o['area_mm2']); deq = 2 * np.sqrt(area / np.pi)
                w.writerow({**{c: '' for c in cols}, 'plane_index': 0, 'region': reg, 'arc_length_mm': '%.6f' % arc_open[(fn, reg)],
                            'arc_length_branch_mm': '%.6f' % arc_open[(fn, reg)], 'area_mm2': '%.6f' % area, 'perimeter_mm': '%.6f' % (np.pi * deq), 'equivalent_diameter_mm': '%.6f' % deq,
                            'hydraulic_diameter_mm': '%.6f' % deq, 'centroid_x': o['cx_mm'], 'centroid_y': o['cy_mm'], 'centroid_z': o['cz_mm'],
                            'centerline_source': 'nostril_opening_cut_plane', 'is_valid': 'True', 'frame_name': fn,
                            'frame_index': frames.index(fn)})
for p in sorted(glob.glob(os.path.join(a.pe_dir, 'csa_pe_*.csv'))):
    reg = os.path.basename(p)[7:-4]; R = list(csv.DictReader(open(p))); pc = list(R[0].keys()); drop = 0
    with open(os.path.join(a.out_dir, 'csa_pe', os.path.basename(p)), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=pc); w.writeheader()
        for r in R:
            if (r['frame'], reg, r['arc_mm']) in bad: drop += 1; continue
            w.writerow(r)
        if reg in SIDE and FILL:
            for r in FP:
                if r['region'] == reg: w.writerow({c: r[c] for c in pc})
        elif reg in SIDE:
            for fn in frames:
                k = cfd_of(fn); o = O[(k, reg)]
                w.writerow(dict(frame=fn, frame_ms=int(fn.split('_')[1]), cfd_idx=k, cfd_t_s=o['cfd_t_s'], region=reg,
                                arc_mm='%.6f' % arc_open[(fn, reg)], area_mm2=o['area_mm2'], TP_area_avg_Pa=o['TP_area_avg_Pa'],
                                TP_mdot_avg_Pa=o['TP_mdot_avg_Pa'], energy_flux_W=o['energy_flux_W'], Q_Lmin=o['Q_Lmin'],
                                area_int_mm2=o['area_mm2'], n_lobes=1))
    print('%s: dropped %d off-centerline rows%s' % (os.path.basename(p), drop, ', added %d opening rows' % len(frames) if reg in SIDE else ''))
