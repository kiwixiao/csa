#!/usr/bin/env python3
"""Airway anatomy: cross-sectional area along a centerline path, through the breath.

Pure geometry (no CFD) from ENT003-Data.csv. Produces a distance x time CSA heatmap,
CSA profiles at several breath phases, and the narrowest-cross-section dynamics.

Ported from /storage4/ENT002/plot_csa_anatomy.py. Differences:
  * --csa / --out / --path are arguments instead of module constants
  * ENT002 was one region; ENT003 needs a *path* (see csa_paths.py), so distance is
    measured from the nostril/lips rather than the mouth
  * the breath-phase selection for the profile plot is driven by the prescribed flow
    profile instead of evenly spaced times, so peak inhale/exhale are always shown

Usage:
  python3 plot_csa_anatomy.py --csa ENT003/csa_bifurcation/ENT003-Data.csv \
      --prof ENT003FlowProfile.csv --path LeftNose --out postproc/csa_anatomy
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

import csa_paths as CP

p = argparse.ArgumentParser()
p.add_argument("--csa", required=True)
p.add_argument("--prof", required=True)
p.add_argument("--path", default="LeftNose")
p.add_argument("--out", required=True)
a = p.parse_args()

os.makedirs(a.out, exist_ok=True)
times, per = CP.load_anatomy(a.csa, a.path)
if not times:
    raise SystemExit("no valid CSA rows for path %r in %s" % (a.path, a.csa))
tt = np.array(times)
NS = max(len(per[t]["s"]) for t in times)
smax = min(per[t]["s"].max() for t in times)      # grid covered by every frame
grid = np.linspace(0, smax, NS)
XL = CP.xlabel(a.path)

A = np.full((len(times), NS), np.nan)
D = np.full((len(times), NS), np.nan)
for i, t in enumerate(times):
    s = per[t]["s"]; o = np.argsort(s)
    A[i] = np.interp(grid, s[o], per[t]["area"][o], left=np.nan, right=np.nan)
    D[i] = np.interp(grid, s[o], per[t]["dh"][o], left=np.nan, right=np.nan)
print("path %s: %d frames  t=%.2f-%.2fs  distance 0-%.0f mm  %d planes"
      % (a.path, len(times), tt[0], tt[-1], smax, NS), flush=True)

ki, ke = CP.peak_frames(a.prof, times)
q, _, _ = CP.phase_of(a.prof, times)

# 1. heatmap: distance x time coloured by CSA
fig, ax = plt.subplots(figsize=(11, 5.5))
im = ax.pcolormesh(grid, tt, A, cmap="turbo", shading="auto")
imin = np.nanargmin(A, axis=1)
ax.plot(grid[imin], tt, "w-", lw=1.0, alpha=0.7, label="narrowest CSA")
for k, lab, c in ((ki, "peak inhale", "w"), (ke, "peak exhale", "0.75")):
    ax.axhline(tt[k], color=c, ls=":", lw=1.0)
    ax.annotate(lab, (smax * 0.99, tt[k]), ha="right", va="bottom", fontsize=8, color=c)
ax.set_xlabel(XL); ax.set_ylabel("time (s)")
ax.set_title(u"%s airway cross-sectional area through the breathing cycle (mm²)  —  %s path" % (os.environ.get("SUBJECT", ""), a.path))
fig.colorbar(im, ax=ax, label=u"CSA (mm²)")
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(a.out + "/csa_anatomy_heatmap.png", dpi=300); plt.close(fig)

# 2. CSA profiles at several breath phases (always including both peaks)
sel = sorted(set(list(np.linspace(0, len(times) - 1, 6).round().astype(int)) + [ki, ke]))
colors = cm.viridis(np.linspace(0, 0.92, len(sel)))
fig, ax = plt.subplots(figsize=(10, 5))
for i, c in zip(sel, colors):
    t = times[i]
    lab = "t=%.2fs (%.1f L/min)" % (t, q[i])
    if i == ki: lab += "  PEAK INHALE"
    if i == ke: lab += "  PEAK EXHALE"
    ax.plot(per[t]["s"], per[t]["area"], color=c,
            lw=2.4 if i in (ki, ke) else 1.4, label=lab)
ax.plot(grid, np.nanmean(A, axis=0), "k--", lw=1.2, alpha=0.6, label="cycle mean")
ax.set_xlabel(XL); ax.set_ylabel(u"cross-section area (mm²)")
ax.set_title("Airway CSA profile along the centerline at breath phases")
ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(a.out + "/csa_anatomy_profiles.png", dpi=300); plt.close(fig)

# 3. narrowest-cross-section dynamics
minA = np.nanmin(A, axis=1)
minloc = grid[np.nanargmin(A, axis=1)]
fig, ax = plt.subplots(figsize=(10, 4.5))
ax.plot(tt, minA, "-o", ms=3, color="#c0392b", label="minimum CSA")
ax.set_xlabel("time (s)"); ax.set_ylabel(u"minimum CSA (mm²)", color="#c0392b")
ax.grid(alpha=0.3); ax.legend(loc="upper left", fontsize=8)
axr = ax.twinx(); axr.plot(tt, minloc, ":", color="#2980b9", lw=1.3)
axr.set_ylabel("location of narrowest point (mm from opening)", color="#2980b9")
ax.set_title("Narrowest airway cross-section through the breath  —  %s path\n"
             "(dynamic constriction: does the min CSA shrink during inspiration?)" % a.path)
fig.tight_layout(); fig.savefig(a.out + "/csa_anatomy_minCSA.png", dpi=300); plt.close(fig)

with open(a.out + "/csa_anatomy_summary.csv", "w") as f:
    f.write("time_s,flow_Lmin,min_CSA_mm2,min_loc_mm_from_opening,mean_CSA_mm2\n")
    for i, t in enumerate(times):
        f.write("%.4f,%.3f,%.3f,%.2f,%.3f\n"
                % (t, q[i], minA[i], minloc[i], np.nanmean(A[i])))

print("wrote 3 figures + summary to", a.out)
print("min CSA over breath: %.2f mm2 (t=%.2fs) .. %.2f mm2 (widest min)"
      % (minA.min(), tt[int(np.argmin(minA))], minA.max()))
print("dynamic narrowing: %.0f%%   narrowest point ~%.0f mm from the opening (median)"
      % (100 * (1 - minA.min() / minA.max()), np.median(minloc)))
