"""Close holes in the nasal septum of a capped airway surface (CSA geometry only - not the CFD mesh).

A septal hole makes the left and right nasal passages one connected surface, which defeats the
CSA slicer's left/right split. Each hole is a short tunnel through the septum. We find it as the
place where surface walks from the two nostrils meet (find_septal_holes.py logic), cut the tunnel
along that meeting ring (removing the faces that touch it), and cap the two cut ends with flat
fans. Re-checks that the fronts no longer meet anywhere in front of the septum.

Usage: pvb.sh patch_septal_holes.py in.stl out.stl [--front-of-y 38]
"""
import argparse, re
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser(); ap.add_argument('stl_in'); ap.add_argument('stl_out')
ap.add_argument('--front-of-y', type=float, default=38.0); ap.add_argument('--min-z', type=float, default=1845.0)
a = ap.parse_args()

txt = open(a.stl_in).read(); S = {}
for m in re.finditer(r'solid (\w+)\n(.*?)endsolid \1', txt, re.S):
    S[m.group(1)] = np.array(re.findall(r'vertex (\S+) (\S+) (\S+)', m.group(2)), float).reshape(-1, 3, 3)

def mesh(T):
    V, inv = np.unique(np.round(T.reshape(-1, 3), 6), axis=0, return_inverse=True); return V, inv.reshape(-1, 3)

def meeting(V, F):
    E = np.unique(np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1), axis=0)
    w = np.linalg.norm(V[E[:, 0]] - V[E[:, 1]], axis=1)
    G = csr_matrix((np.r_[w, w], (np.r_[E[:, 0], E[:, 1]], np.r_[E[:, 1], E[:, 0]])), shape=(len(V), len(V)))
    tree = cKDTree(V)
    seeds = {}
    for cap in ('Leftnose', 'Rightnose'):
        d, i = tree.query(S[cap].reshape(-1, 3), distance_upper_bound=1e-4); seeds[cap] = np.unique(i[np.isfinite(d)])
    dL = dijkstra(G, directed=False, indices=seeds['Leftnose'], min_only=True)
    dR = dijkstra(G, directed=False, indices=seeds['Rightnose'], min_only=True)
    s = np.sign(dL - dR); cross = E[s[E[:, 0]] != s[E[:, 1]]]
    mid = 0.5 * (V[cross[:, 0]] + V[cross[:, 1]])
    keep = (mid[:, 1] > a.front_of_y) & (mid[:, 2] > a.min_z)
    return cross[keep]

V, F = mesh(S['wall'])
cross = meeting(V, F)
bad = np.zeros(len(V), bool); bad[np.unique(cross.ravel())] = True
print('meeting ring(s) in front of y=%.0f mm: %d edges, %d vertices' % (a.front_of_y, len(cross), bad.sum()))
if not bad.any():
    print('septum already intact - nothing to patch'); raise SystemExit(0)
drop = bad[F].any(axis=1)                              # faces touching the meeting ring = the cut
F2 = F[~drop]
# clean the cut edge: remove faces at "pinch" vertices (more than two boundary edges) and
# dangling faces (two or three open edges) until every open ring is a simple loop
for it in range(50):
    E = np.sort(np.vstack([F2[:, [0, 1]], F2[:, [1, 2]], F2[:, [2, 0]]]), axis=1)
    u, inv, cnt = np.unique(E, axis=0, return_inverse=True, return_counts=True)
    open_e = cnt[inv] == 1
    deg = np.bincount(u[cnt == 1].ravel(), minlength=len(V))
    pinch = deg > 2
    open_per_face = open_e.reshape(3, -1).sum(0)
    rm = pinch[F2].any(1) | (open_per_face >= 2)
    if not rm.any(): break
    F2 = F2[~rm]
print('cut: removed %d wall faces (%d cleanup passes)' % (len(F) - len(F2), it))

# open rings of the cut wall = 3 cap rings + the new tunnel ends
E = np.vstack([F2[:, [0, 1]], F2[:, [1, 2]], F2[:, [2, 0]]])
key = np.sort(E, axis=1); u, cnt = np.unique(key, axis=0, return_counts=True)
B = u[cnt == 1]
nb = np.unique(B.ravel()); idx = {v: i for i, v in enumerate(nb)}
A = csr_matrix((np.ones(len(B)), ([idx[x] for x in B[:, 0]], [idx[x] for x in B[:, 1]])), shape=(len(nb), len(nb)))
nc, lab = connected_components(A, directed=False)
cap_pts = cKDTree(np.vstack([S[c].reshape(-1, 3) for c in ('Leftnose', 'Rightnose', 'outlet')]))
patches = []
for k in range(nc):
    ring_v = nb[lab == k]
    if cap_pts.query(V[ring_v[0]])[0] < 1e-4:          # one of the original cap rings
        continue
    # order the ring by walking its boundary edges
    redges = B[np.isin(B[:, 0], ring_v)]
    nbr = {}
    for x, y in redges: nbr.setdefault(x, []).append(y); nbr.setdefault(y, []).append(x)
    if any(len(v) != 2 for v in nbr.values()):
        raise SystemExit('tunnel end ring is not a simple loop - inspect manually')
    order = [ring_v[0]]; prev = None
    while True:
        nx = [y for y in nbr[order[-1]] if y != prev][0]
        if nx == order[0]: break
        prev = order[-1]; order.append(nx)
    q = V[order]; c = q.mean(0)
    # orient the fan like the neighbouring wall: outward = consistent winding with the adjacent face
    tris = np.array([[c, q[i], q[(i + 1) % len(q)]] for i in range(len(q))])
    e0, e1 = order[0], order[1]
    adj = F2[(np.isin(F2, [e0]).any(1)) & (np.isin(F2, [e1]).any(1))][0]
    pos = list(adj); fwd = (pos[(pos.index(e0) + 1) % 3] == e1)   # wall face traverses e0->e1 ?
    if fwd: tris = tris[:, [0, 2, 1]]                              # the patch must traverse it e1->e0
    patches.append(tris)
    print('  patched tunnel end: ring of %d vertices, centre (%.1f, %.1f, %.1f) mm, max radius %.1f mm'
          % (len(q), *c, np.linalg.norm(q - c, axis=1).max()))
wall_new = np.vstack([V[F2]] + patches)

# re-check: rebuild and test the fronts again
Vn, Fn = mesh(wall_new)
left = meeting(Vn, Fn)
print('re-check: %d meeting edges in front of y=%.0f mm (0 = septum now separates left/right)' % (len(left), a.front_of_y))
Et = np.sort(np.vstack([Fn[:, [0, 1]], Fn[:, [1, 2]], Fn[:, [2, 0]]]), axis=1); _, ct = np.unique(Et, axis=0, return_counts=True)
print('wall boundary edges after patch: %d (should equal the 3 cap rings)' % int((ct == 1).sum()))

def write_solid(fh, name, T):
    fh.write('solid %s\n' % name)
    N = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); N /= np.linalg.norm(N, axis=1, keepdims=True) + 1e-30
    for t, n in zip(T, N):
        fh.write('  facet normal %e %e %e\n    outer loop\n' % tuple(n))
        for v in t: fh.write('      vertex %e %e %e\n' % tuple(v))
        fh.write('    endloop\n  endfacet\n')
    fh.write('endsolid %s\n' % name)
with open(a.stl_out, 'w') as fh:
    for nm in ('Leftnose', 'Rightnose', 'outlet'): write_solid(fh, nm, S[nm])
    write_solid(fh, 'wall', wall_new)
print('wrote', a.stl_out)
