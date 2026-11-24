"""Locate holes in the nasal septum of a capped airway surface (seg_0_boundary.stl layout).

Walk the surface from each nostril cap ring (Dijkstra on the vertex graph). In a nose with an
intact septum the left and right fronts only meet behind the septum, where the passages merge;
any other place where they meet is a tunnel through the septum. Reports each hole's centre/size.
Usage: pvb.sh find_septal_holes.py seg_0_boundary.stl --front-of-y 40 [--json out.json]
"""
import argparse, json, re
import numpy as np, vtk
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser(); ap.add_argument('stl'); ap.add_argument('--front-of-y', type=float, default=40.0,
    help='only report meetings anterior of this y (mm) - behind it the passages merge naturally')
ap.add_argument('--json', default=None); a = ap.parse_args()
txt = open(a.stl).read(); S = {}
for m in re.finditer(r'solid (\w+)\n(.*?)endsolid \1', txt, re.S):
    S[m.group(1)] = np.array(re.findall(r'vertex (\S+) (\S+) (\S+)', m.group(2)), float).reshape(-1, 3, 3)
W = S['wall'].reshape(-1, 3)
V, inv = np.unique(np.round(W, 6), axis=0, return_inverse=True); F = inv.reshape(-1, 3)
E = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]); E = np.unique(np.sort(E, axis=1), axis=0)
w = np.linalg.norm(V[E[:, 0]] - V[E[:, 1]], axis=1)
G = csr_matrix((np.r_[w, w], (np.r_[E[:, 0], E[:, 1]], np.r_[E[:, 1], E[:, 0]])), shape=(len(V), len(V)))
tree = cKDTree(V)
def ring_vertices(cap):
    return np.unique(tree.query(cap.reshape(-1, 3), distance_upper_bound=1e-4)[1])
seedL = ring_vertices(S['Leftnose']); seedR = ring_vertices(S['Rightnose'])
seedL = seedL[seedL < len(V)]; seedR = seedR[seedR < len(V)]
dL = dijkstra(G, directed=False, indices=seedL, min_only=True); dR = dijkstra(G, directed=False, indices=seedR, min_only=True)
side = np.sign(dL - dR)
cross = E[(side[E[:, 0]] != side[E[:, 1]])]
mid = 0.5 * (V[cross[:, 0]] + V[cross[:, 1]])
nasal = (mid[:, 1] > a.front_of_y) & (mid[:, 2] > 1845)
cross, mid = cross[nasal], mid[nasal]
# group the meeting edges into separate places (edges sharing a vertex = same place)
nodes = np.unique(cross.ravel()); idx = {v: i for i, v in enumerate(nodes)}
A = csr_matrix((np.ones(len(cross)), ([idx[u] for u in cross[:, 0]], [idx[v] for v in cross[:, 1]])), shape=(len(nodes), len(nodes)))
nc, lab = connected_components(A, directed=False)
holes = []
for k in range(nc):
    P = V[nodes[lab == k]]
    if len(P) < 6: continue
    c = P.mean(0); r = float(np.percentile(np.linalg.norm(P - c, axis=1), 95))
    holes.append(dict(centre=c.round(2).tolist(), radius_mm=round(r, 2), n_vertices=int(len(P)), geodesic_from_left=float(dL[nodes[lab == k]].mean()), geodesic_from_right=float(dR[nodes[lab == k]].mean())))
holes.sort(key=lambda h: -h['n_vertices'])
print('left/right fronts meet at %d places anterior of y=%.0f mm:' % (len(holes), a.front_of_y))
for h in holes: print('   centre (%.1f, %.1f, %.1f) mm, radius ~%.1f mm, %d vertices' % (*h['centre'], h['radius_mm'], h['n_vertices']))
if a.json: json.dump(holes, open(a.json, 'w'), indent=1)
