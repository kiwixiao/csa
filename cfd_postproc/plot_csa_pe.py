#!/usr/bin/env python3
"""Plots for the CSA total-pressure / energy-flux drop analysis along one airway path.

Reads the per-region csa_pe_<Region>.csv written by csa_pressure_energy.py and
concatenates the regions of --path (see csa_paths.py). x = distance from the airway
opening (0 = nostril or lips, increasing toward the trachea). Cumulative loss is
referenced from the opening.

Ported from /storage4/ENT002/plot_csa_pe.py. Differences:
  * reads a *path* of region CSVs rather than one whole-case CSV
  * --prof drives the inhale/exhale annotation instead of ENT002's hardcoded
    "positive = inhale" (ENT003's prescribed profile is positive on EXHALATION,
    and the record starts mid-exhalation)
  * no T_MAX cut: the corrected export is a clean 0.025-2.175 s already

Usage:
  python3 plot_csa_pe.py --pe-dir postproc/csa_pe --prof ENT003FlowProfile.csv \
      --path LeftNose --out postproc/csa_pe/plots/LeftNose
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
p.add_argument("--pe-dir", required=True)
p.add_argument("--prof", required=True)
p.add_argument("--path", default="LeftNose")
p.add_argument("--out", required=True)
p.add_argument("--qmin", type=float, default=0.5,
               help="L/min below which the mass-flow-averaged TP is masked")
a = p.parse_args()

os.makedirs(a.out, exist_ok=True)
times, data = CP.load_pe(a.pe_dir, a.path, qmin=a.qmin)
if not times:
    raise SystemExit("no rows for path %r under %s" % (a.path, a.pe_dir))
tt = np.array(times)
XL = CP.xlabel(a.path)
q, _, _ = CP.phase_of(a.prof, times)
ki, ke = CP.peak_frames(a.prof, times)
print("path %s: %d time frames, %d planes each" % (a.path, len(times), len(data[times[0]]["s"])),
      flush=True)

sel = sorted(set(list(np.linspace(0, len(times) - 1, 5).round().astype(int)) + [ki, ke]))
colors = cm.viridis(np.linspace(0, 0.92, len(sel)))


def lab(i):
    s = "t=%.2fs (%+.1f L/min)" % (times[i], q[i])
    return s + ("  PEAK INHALE" if i == ki else ("  PEAK EXHALE" if i == ke else ""))


def lw(i):
    return 2.4 if i in (ki, ke) else 1.4


# 1. CSA area
fig, ax = plt.subplots(figsize=(9.5, 4.8))
for i, c in zip(sel, colors):
    ax.plot(data[times[i]]["s"], data[times[i]]["area"], color=c, lw=lw(i), label=lab(i))
ax.set_xlabel(XL); ax.set_ylabel(u"cross-section area (mm²)")
ax.set_title("%s airway cross-sectional area  —  %s path" % (os.environ.get("SUBJECT", ""), a.path))
ax.legend(fontsize=7.5, ncol=2); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(a.out + "/csa_area_profile.png", dpi=300); plt.close(fig)

# 2. total pressure (area-avg solid, mass-flow-avg dashed)
fig, ax = plt.subplots(figsize=(9.5, 4.8))
for i, c in zip(sel, colors):
    d = data[times[i]]
    ax.plot(d["s"], d["tp_a"], color=c, lw=lw(i), label=lab(i))
    ax.plot(d["s"], d["tp_m"], color=c, lw=1.0, ls="--", alpha=0.7)
ax.set_xlabel(XL); ax.set_ylabel("total pressure (Pa)")
ax.set_title("Total pressure  (solid = area-avg, dashed = mass-flow-avg)")
ax.legend(fontsize=7.5, ncol=2); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(a.out + "/csa_pressure_profile.png", dpi=300); plt.close(fig)

# 3. cumulative total-pressure loss from the opening
fig, ax = plt.subplots(figsize=(9.5, 4.8))
for i, c in zip(sel, colors):
    d = data[times[i]]
    ax.plot(d["s"], CP.cum_from_opening(d["tp_a"]), color=c, lw=lw(i) + 0.2, label=lab(i))
    ax.plot(d["s"], CP.cum_from_opening(d["tp_m"]), color=c, lw=1.0, ls="--", alpha=0.7)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel(XL); ax.set_ylabel("cumulative total-pressure loss from opening (Pa)")
ax.set_title("Cumulative total-pressure loss  (solid = area-avg, dashed = mass-flow-avg)\n"
             "positive = inhalation loss; negative = exhalation (flow reversed)")
ax.legend(fontsize=7.5, ncol=2); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(a.out + "/csa_pressure_drop.png", dpi=300); plt.close(fig)

# 4. energy flux + cumulative energy-flux loss
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 4.8))
for i, c in zip(sel, colors):
    d = data[times[i]]
    a1.plot(d["s"], d["ef"], color=c, lw=lw(i), label=lab(i))
    a2.plot(d["s"], CP.cum_from_opening(d["ef"]), color=c, lw=lw(i))
a1.set_xlabel(XL); a1.set_ylabel("energy flux through plane (W)")
a1.set_title("Energy flux"); a1.legend(fontsize=7.5); a1.grid(alpha=0.3)
a2.set_xlabel(XL); a2.set_ylabel("cumulative energy-flux loss from opening (W)")
a2.set_title("Cumulative energy-flux loss"); a2.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(a.out + "/csa_energy_flux.png", dpi=300); plt.close(fig)

# 5. time series: opening->trachea loss + bulk flow
dP_a = np.array([data[t]["tp_a"][0] - data[t]["tp_a"][-1] for t in times])


def endpoints(v):
    f = v[np.isfinite(v)]
    return (f[0] - f[-1]) if f.size >= 2 else np.nan


dP_m = np.array([endpoints(data[t]["tp_m"]) for t in times])
Qb = np.array([np.median(data[t]["Q"][:max(1, len(data[t]["Q"]) // 3)]) for t in times])
dP_m[np.abs(Qb) < a.qmin] = np.nan
dEF = np.array([data[t]["ef"][0] - data[t]["ef"][-1] for t in times])

fig, ax = plt.subplots(figsize=(10.5, 4.8))
ax.plot(tt, dP_a, "-o", ms=3, color="#c0392b", label=u"ΔP opening→trachea (area-avg)")
ax.plot(tt, dP_m, "--s", ms=3, color="#e67e22", label=u"ΔP opening→trachea (mass-flow-avg)")
ax.axhline(0, color="k", lw=0.5)
ax.axvline(tt[ki], color="#2471a3", ls=":", lw=1.0)
ax.annotate("peak inhale", (tt[ki], ax.get_ylim()[1]), fontsize=8, color="#2471a3",
            ha="center", va="top")
ax.set_xlabel("time (s)"); ax.set_ylabel("total-pressure loss (Pa)", color="#c0392b")
ax.grid(alpha=0.3); ax.legend(loc="upper left", fontsize=8)
axr = ax.twinx(); axr.plot(tt, q, ":", color="#2980b9", lw=1.4)
axr.set_ylabel("prescribed flow (L/min, + = exhale)", color="#2980b9")
ax.set_title(u"Opening→trachea total-pressure loss vs time  —  %s path" % a.path)
fig.tight_layout(); fig.savefig(a.out + "/csa_drop_timeseries.png", dpi=300); plt.close(fig)

# 6. heatmap distance x time (cumulative pressure loss, area-avg)
s0 = data[times[0]]["s"]
M = np.full((len(times), len(s0)), np.nan)
for i, t in enumerate(times):
    dp = CP.cum_from_opening(data[t]["tp_a"])
    M[i] = dp if len(data[t]["s"]) == len(s0) else np.interp(s0, data[t]["s"][::-1], dp[::-1])
fig, ax = plt.subplots(figsize=(10.5, 5.2))
vmax = np.nanpercentile(np.abs(M), 99)
im = ax.pcolormesh(s0, tt, M, cmap="RdBu_r", shading="auto", vmin=-vmax, vmax=vmax)
ax.set_xlabel(XL); ax.set_ylabel("time (s)")
ax.set_title("Cumulative total-pressure loss from the opening (Pa) — area-avg\n"
             "red = inhalation loss, blue = exhalation")
fig.colorbar(im, ax=ax, label=u"ΔP from opening (Pa)")
fig.tight_layout(); fig.savefig(a.out + "/csa_pressure_heatmap.png", dpi=300); plt.close(fig)

with open(a.out + "/csa_drop_summary.csv", "w") as f:
    f.write("time_s,flow_Lmin,dP_opening_trachea_area_Pa,dP_opening_trachea_mdot_Pa,energy_flux_drop_W\n")
    for i, t in enumerate(times):
        f.write("%.4f,%.3f,%.3f,%.3f,%.6f\n" % (t, q[i], dP_a[i], dP_m[i], dEF[i]))

print("wrote 6 figures + summary to", a.out)
k = int(np.nanargmax(np.abs(dP_a)))
print("peak |dP| area-avg: %.1f Pa at t=%.2fs (prescribed flow %.1f L/min)" % (abs(dP_a[k]), tt[k], q[k]))
