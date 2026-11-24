#!/usr/bin/env python3
"""Arc-length landmark figure for marking anatomical regions along a nostril -> trachea path.

Left : sagittal airway silhouette (face to the left, as in the renders) with the path's
       centerline, a dot every 10 mm and a label every 20 mm of distance from the opening.
Right: CSA vs the same distance (median over the breath + min-max envelope), gridlines every
       10 mm, the exact nose -> shared-airway junction and the dissipation bottleneck marked.
Usage: python3 make_landmark_figure.py --subject-dir S --csa S-Data.csv --path LeftNose
                                      --frame out_000000 --out fig.png [--bottleneck-mm 157]
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import csa_paths as CP

p = argparse.ArgumentParser()
p.add_argument("--subject-dir", required=True); p.add_argument("--csa", required=True)
p.add_argument("--path", default="LeftNose"); p.add_argument("--frame", default="out_000000")
p.add_argument("--out", required=True); p.add_argument("--bottleneck-mm", type=float, default=None)
a = p.parse_args()

times, per = CP.load_anatomy(a.csa, a.path)
t0 = int(a.frame.split("_")[1]) / 1000.0
d0 = per[min(times, key=lambda t: abs(t - t0))]
o = np.argsort(d0["s"]); S, C, R = d0["s"][o], d0["cxyz"][o], d0["region"][o]
CL = CP.side_centerline(a.subject_dir, a.path, a.frame, s_max=float(S.max()))
LINE_P, LINE_S = (CL[0], CL[1]) if CL is not None else (C, S)      # draw the real centerline
nose_regs = [r for r in CP.regions_of(a.path) if r != "DescendingAirway"]
junction = float(S[np.isin(R, nose_regs)].max()) if np.isin(R, nose_regs).any() else None
_, prox, dist = CP.PATH_LABEL[a.path]

# CSA over the breath on a common distance grid (each frame's own distance axis)
grid = np.arange(0, np.floor(S.max()) + 1, 1.0)
A = np.array([np.interp(grid, np.sort(per[t]["s"]), per[t]["area"][np.argsort(per[t]["s"])], left=np.nan, right=np.nan) for t in times])
med, lo, hi = np.nanmedian(A, 0), np.nanmin(A, 0), np.nanmax(A, 0)

V = CP.stl_vertices(os.path.join(a.subject_dir, "motion", "stl", "%s.stl" % a.frame))
fig = plt.figure(figsize=(17, 11.5), facecolor="white")
ax = fig.add_axes([0.02, 0.05, 0.40, 0.88])
ax.scatter(V[:, 1], V[:, 2], s=0.8, c="#33404f", alpha=0.05, edgecolors="none", rasterized=True)
pts = np.column_stack([LINE_P[:, 1], LINE_P[:, 2]]).reshape(-1, 1, 2)
lc = LineCollection(np.concatenate([pts[:-1], pts[1:]], axis=1), cmap="viridis", array=LINE_S[:-1], linewidths=4)
lc.set_clim(0, S.max()); ax.add_collection(lc)
cmap = plt.get_cmap("viridis")
for d in range(0, int(S.max()) + 1, 10):
    i = int(np.argmin(np.abs(LINE_S - d))); col = cmap(d / S.max())
    big = d % 20 == 0
    ax.plot(LINE_P[i, 1], LINE_P[i, 2], "o", ms=8 if big else 5, mfc=col, mec="#222", mew=1.0, zorder=6)
    if big:
        # label on the side away from the face (face is at large y = left of the plot)
        ax.annotate("%d" % d, (LINE_P[i, 1], LINE_P[i, 2]), xytext=(-12, 0), textcoords="offset points", ha="left" if LINE_P[i, 1] > np.median(LINE_P[:, 1]) + 30 else "right",
                    va="center", fontsize=11, fontweight="bold", color="#15181d",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=col, lw=1.2, alpha=0.9))
ax.set_xlim(V[:, 1].max() + 12, V[:, 1].min() - 25)     # +y (face) on the LEFT, as in the renders
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("%s - %s path, frame %s\ncenterline distance from the %s (mm): dot every 10 mm, label every 20 mm"
             % (os.path.basename(a.subject_dir.rstrip("/")), a.path, a.frame, prox), fontsize=12.5, fontweight="bold", color="#12305c", loc="left")

bx = fig.add_axes([0.50, 0.08, 0.46, 0.80])
bx.fill_between(grid, lo, hi, color="#2a78d6", alpha=0.18, lw=0, label="min-max over the breath")
bx.plot(grid, med, color="#2a78d6", lw=2, label="median over the breath")
for d in range(0, int(S.max()) + 1, 10):
    bx.axvline(d, color="#b9b8b2" if d % 20 else "#7f7e78", lw=0.6 if d % 20 else 1.0, zorder=0)
if junction is not None:
    bx.axvline(junction, color="#2e7d32", lw=2.2, ls="--")
    bx.text(junction + 1, np.nanmax(hi) * 0.97, "nose ends here\n(%s joins the\nshared airway, %.0f mm)\nEXACT" % (a.path, junction),
            color="#2e7d32", fontsize=10.5, fontweight="bold", va="top")
if a.bottleneck_mm is not None:
    bx.axvline(a.bottleneck_mm, color="#c0392b", lw=2.2, ls="--")
    bx.text(a.bottleneck_mm + 1, np.nanmax(hi) * 0.60, "pressure-loss\nbottleneck\n(%.0f mm)" % a.bottleneck_mm, color="#c0392b", fontsize=10.5, fontweight="bold", va="top")
bx.set_xticks(range(0, int(S.max()) + 1, 20)); bx.set_xticks(range(0, int(S.max()) + 1, 10), minor=True)
bx.set_xlim(0, S.max()); bx.set_ylim(0, np.nanmax(hi) * 1.05)
bx.set_xlabel("distance from the %s along the centerline (mm)   ->   %s" % (prox, dist), fontsize=12)
bx.set_ylabel(u"cross-sectional area (mm²)", fontsize=12)
bx.set_title("Cross-section along the same path (0 = nostril opening = the cap)", fontsize=12.5, fontweight="bold", color="#12305c", loc="left")
bx.legend(loc="upper right", frameon=False); bx.grid(False)
for sp in ("top", "right"): bx.spines[sp].set_visible(False)
fig.savefig(a.out, dpi=110, facecolor="white")
print("wrote", a.out, "| path length %.0f mm, nose junction %s mm" % (S.max(), junction))
