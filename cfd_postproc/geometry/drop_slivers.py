"""Drop sub-micron sliver triangles from the wall solid of a named multi-solid ASCII STL.

They turn degenerate after the FFD; trimesh then drops them, the face-label count no longer matches and
csa-bifurcation silently slices the frame-0 meshes for every frame (see qa/facecheck.py).
Usage: python3 drop_slivers.py in.stl out.stl   (wall triangles with min edge < 1e-3 mm or area < 1e-6 mm2)
"""
import re, sys, numpy as np
src, dst = sys.argv[1], sys.argv[2]
txt = open(src).read(); out = []
for m in re.finditer(r"solid (\S+)\n(.*?)endsolid \1", txt, re.S):
    nm = m.group(1); T = np.array(re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", m.group(2)), float).reshape(-1, 3, 3)
    if nm == "wall":
        e = np.stack([np.linalg.norm(T[:, i] - T[:, (i + 1) % 3], axis=1) for i in range(3)], 1)
        a = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)
        keep = (e.min(1) >= 1e-3) & (a >= 1e-6); print("wall: dropped %d of %d triangles" % ((~keep).sum(), len(T))); T = T[keep]
    out.append("solid %s\n" % nm)
    for t in T:
        n = np.cross(t[1] - t[0], t[2] - t[0]); l = np.linalg.norm(n); n = n / l if l > 0 else n
        out.append(" facet normal %e %e %e\n  outer loop\n" % tuple(n) + "".join("   vertex %.6f %.6f %.6f\n" % tuple(x) for x in t) + "  endloop\n endfacet\n")
    out.append("endsolid %s\n" % nm)
open(dst, "w").write("".join(out))
