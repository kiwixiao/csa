#!/usr/bin/env python3
"""Rebuild <frame>-Planes-All.stl per frame from the per-region plane STLs, leaving out rejected
planes (clean_nostril_planes.py's rejected_planes.csv). Plane STL files are numbered by rank among
the frame/region's VALID planes in the ORIGINAL CSA table (same rule as csa_pressure_energy.rank_map).
Usage: clean_combined_planes.py --csa ORIG-Data.csv --frames-dir .../csa_bifurcation/frames --rejected rejected_planes.csv --out DIR"""
import argparse, csv, os
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--csa', required=True); ap.add_argument('--frames-dir', required=True)
ap.add_argument('--rejected', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--extra-dir', default=None, help='per-frame <frame>-fill.stl disks to add (fill_nostril_planes.py)'); a = ap.parse_args()
TAG = {'DescendingAirway': 'Desc', 'LeftNose': 'Left', 'RightNose': 'Right', 'Mouth': 'Mouth'}
valid = {}; allp = {}
for d in csv.DictReader(open(a.csa)):
    if d['is_valid'].lower() == 'true': valid.setdefault((d['frame_name'], d['region']), []).append(int(d['plane_index']))
for d in csv.DictReader(open(a.csa.replace('-Data.csv', '-Data-all.csv'))):          # disks are numbered over ALL planes
    allp.setdefault((d['frame_name'], d['region']), []).append(int(d['plane_index']))
rej = {(d['frame_name'], d['region'], int(d['plane_index'])) for d in csv.DictReader(open(a.rejected))}
dt = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
def tris(path):
    raw = open(path, 'rb').read(); n = int(np.frombuffer(raw[80:84], '<u4')[0])
    if 84 + 50 * n == len(raw): return np.frombuffer(raw[84:], dtype=dt, count=n)['v']
    V = np.array([l.split()[1:4] for l in raw.decode().splitlines() if l.strip().startswith('vertex')], float)
    return V.reshape(-1, 3, 3)
os.makedirs(a.out, exist_ok=True); frames = sorted({k[0] for k in valid}); kept = dropped = 0
for fn in frames:
    T = []
    for reg, tag in TAG.items():
        pl = sorted(valid.get((fn, reg), [])); al = sorted(allp.get((fn, reg), []))
        for pidx in pl:
            rank = al.index(pidx)
            if (fn, reg, pidx) in rej: dropped += 1; continue
            p = os.path.join(a.frames_dir, fn, '%s_Planes' % reg, '%s-%s-%03d.stl' % (fn, tag, rank))
            if os.path.exists(p): T.append(tris(p)); kept += 1
    if a.extra_dir and os.path.exists(os.path.join(a.extra_dir, fn + '-fill.stl')):
        T.append(tris(os.path.join(a.extra_dir, fn + '-fill.stl')))
    T = np.concatenate(T).astype('<f4'); rec = np.zeros(len(T), dtype=dt); rec['v'] = T
    N = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); rec['n'] = N / (np.linalg.norm(N, axis=1, keepdims=True) + 1e-30)
    with open(os.path.join(a.out, '%s-Planes-All.stl' % fn), 'wb') as fh:
        fh.write(b'cleaned combined CSA planes'.ljust(80, b' ')); fh.write(np.uint32(len(rec)).tobytes()); fh.write(rec.tobytes())
print('frames %d: kept %d plane disks, left out %d rejected ones -> %s' % (len(frames), kept, dropped, a.out))
