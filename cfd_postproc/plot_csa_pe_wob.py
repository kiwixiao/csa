#!/usr/bin/env python3
"""Pressure-loss and energy-loss (work of breathing) maps along one airway path.

Distance x time envelopes of cumulative total-pressure loss and cumulative energy
dissipation, plus the instantaneous dissipated power and cumulative work-of-breathing
time series. Writes wob_summary.csv, which build_onepager.py consumes.

Ported from /storage4/ENT002/plot_csa_pe_wob.py. Differences:
  * reads a *path* of per-region csa_pe CSVs (see csa_paths.py) instead of one CSV
  * the inspiratory fraction comes from the prescribed profile rather than ENT002's
    hardcoded ``tt < 1.5``; ENT003 exhales first and inhales after t~1.36 s
  * bulk flow Q is taken from the prescribed profile (exact, and immune to the
    per-plane flux noise at a flow reversal) rather than the per-plane median

Usage:
  python3 plot_csa_pe_wob.py --pe-dir postproc/csa_pe --prof ENT003FlowProfile.csv \
      --path LeftNose --out postproc/csa_pe/plots_wob/LeftNose
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import csa_paths as CP

p = argparse.ArgumentParser()
p.add_argument("--pe-dir", required=True)
p.add_argument("--prof", required=True)
p.add_argument("--path", default="LeftNose")
p.add_argument("--out", required=True)
p.add_argument("--flow-share", type=float, default=1.0,
               help="fraction of the prescribed flow carried by this path "
                    "(0.5 for one nostril of a symmetric pair; 1.0 = whole airway)")
a = p.parse_args()

os.makedirs(a.out, exist_ok=True)
times, data = CP.load_pe(a.pe_dir, a.path)
if not times:
    raise SystemExit("no rows for path %r under %s" % (a.path, a.pe_dir))
tt = np.array(times)
XL = CP.xlabel(a.path)
NS = max(len(data[t]["s"]) for t in times)
smax = min(data[t]["s"].max() for t in times)
grid = np.linspace(0, smax, NS)

qL, inh, exh = CP.phase_of(a.prof, times)              # L/min, + = exhale
Qmag = np.abs(qL) / 60000.0 * a.flow_share             # m3/s carried by this path

# cumulative total-pressure loss from the opening (time x distance)
CPm = np.full((len(times), NS), np.nan)
for i, t in enumerate(times):
    s = data[t]["s"]; o = np.argsort(s)
    cp = CP.cum_from_opening(data[t]["tp_a"])
    CPm[i] = np.interp(grid, s[o], cp[o], left=np.nan, right=np.nan)

# energy dissipation = pressure loss x bulk flow (the raw per-plane EnergyFlux field
# integrates too noisily to differentiate along the airway)
CE = CPm * Qmag[:, None] * 1000.0                      # mW

tp_tot = np.array([data[t]["tp_a"][0] - data[t]["tp_a"][-1] for t in times])   # opening->trachea
powerP = np.abs(tp_tot) * Qmag                          # W
dt = np.diff(tt, prepend=tt[0])
wob_cum = np.cumsum(powerP * dt)                        # J
e_tot = tp_tot * Qmag

ki, ke = CP.peak_frames(a.prof, times)


def envelope(M, fname, title, ylabel):
    Mn, Mx, Mmean = np.nanmin(M, axis=0), np.nanmax(M, axis=0), np.nanmean(M, axis=0)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.fill_between(grid, Mn, Mx, color="#3498db", alpha=0.18, label=u"min–max over breath")
    ax.plot(grid, Mmean, "--", color="0.4", lw=1.3, label="cycle mean")
    ax.plot(grid, M[ki], color="#c0392b", lw=2.2, label="peak inhale (t=%.2fs)" % tt[ki])
    ax.plot(grid, M[ke], color="#2471a3", lw=2.2, label="peak exhale (t=%.2fs)" % tt[ke])
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel(XL); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.3); ax.set_xlim(0, smax)
    fig.tight_layout(); fig.savefig(a.out + "/" + fname, dpi=300); plt.close(fig)


envelope(CPm, "wob_pressure_loss_envelope.png",
         "Cumulative total-pressure loss along the airway, through the breath  —  %s path\n"
         u"(min–max band = full breath range)" % a.path,
         "cumulative total-pressure loss from opening (Pa)")
envelope(CE, "wob_energy_loss_envelope.png",
         u"Cumulative energy loss = ΔP·Q (work-of-breathing rate) along the airway  —  %s path\n"
         u"(min–max band = full breath range)" % a.path,
         "cumulative energy loss from opening (mW)")


def heat(M, fname, title, label, cmap="RdBu_r"):
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    v = np.nanpercentile(np.abs(M), 99)
    im = ax.pcolormesh(grid, tt, M, cmap=cmap, shading="auto", vmin=-v, vmax=v)
    ax.set_xlabel(XL); ax.set_ylabel("time (s)"); ax.set_title(title)
    fig.colorbar(im, ax=ax, label=label)
    fig.tight_layout(); fig.savefig(a.out + "/" + fname, dpi=300); plt.close(fig)


heat(CPm, "wob_pressure_loss_heatmap.png",
     "Cumulative total-pressure loss from the opening (Pa)  —  %s path" % a.path,
     u"ΔP (Pa)")
heat(CE, "wob_energy_loss_heatmap.png",
     u"Cumulative energy loss ΔP·Q (mW)  —  %s path" % a.path, "mW")

# LOCAL per-mm dissipation: where along the airway the loss actually concentrates
dEdx = np.gradient(np.abs(CE), grid, axis=1)
fig, ax = plt.subplots(figsize=(10.5, 5.5))
v = np.nanpercentile(np.abs(dEdx), 99)
im = ax.pcolormesh(grid, tt, dEdx, cmap="inferno", shading="auto", vmin=0, vmax=v)
ax.set_xlabel(XL); ax.set_ylabel("time (s)")
ax.set_title("LOCAL dissipation per mm of airway (mW/mm)  —  %s path\n"
             "bright bands = where the breathing effort is spent" % a.path)
fig.colorbar(im, ax=ax, label="mW/mm")
fig.tight_layout(); fig.savefig(a.out + "/wob_local_dissipation.png", dpi=300); plt.close(fig)

# time series: power + cumulative WOB
fig, ax = plt.subplots(figsize=(10.5, 4.6))
ax.plot(tt, powerP * 1000, "-o", ms=3, color="#c0392b", label=u"dissipated power ΔP·Q (mW)")
if inh.any():        # a partial record (resume in progress) can be all-exhalation
    ax.axvspan(tt[inh].min(), tt[inh].max(), color="#2471a3", alpha=0.08, label="inhalation")
ax.set_xlabel("time (s)"); ax.set_ylabel("resistive power (mW)", color="#c0392b")
ax.grid(alpha=0.3); ax.legend(loc="upper left", fontsize=8)
axr = ax.twinx(); axr.plot(tt, wob_cum * 1000, color="#2980b9", lw=1.6)
axr.set_ylabel("cumulative work of breathing (mJ)", color="#2980b9")
ax.set_title("Airway resistive power and cumulative work of breathing  —  %s path" % a.path)
fig.tight_layout(); fig.savefig(a.out + "/wob_timeseries.png", dpi=300); plt.close(fig)

wob_total = float(np.sum(powerP * dt))
wob_insp = float(np.sum((powerP * dt)[inh]))
with open(a.out + "/wob_summary.csv", "w") as f:
    f.write("time_s,flow_Lmin,is_inhale,dissipated_power_mW,cum_WOB_mJ,dP_opening_trachea_Pa\n")
    for i, t in enumerate(times):
        f.write("%.4f,%.3f,%d,%.4f,%.4f,%.3f\n"
                % (t, qL[i], int(inh[i]), powerP[i] * 1000, wob_cum[i] * 1000, tp_tot[i]))

print("wrote 6 figures + wob_summary.csv to", a.out)
print("peak dissipated power: %.2f mW at t=%.2fs" % (powerP.max() * 1000, tt[int(np.argmax(powerP))]))
print("resistive work of breathing over the record: %.3f mJ  (inspiratory %.3f mJ = %.0f%%)"
      % (wob_total * 1000, wob_insp * 1000, 100 * wob_insp / wob_total if wob_total else 0))
print("peak opening->trachea dP: %.1f Pa at t=%.2fs"
      % (np.abs(tp_tot).max(), tt[int(np.argmax(np.abs(tp_tot)))]))
