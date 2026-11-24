#!/usr/bin/env python3
"""Cumulative WOB E(x) and its spatial derivative dE/dx along the airway.

E(x) is the running energy loss from the distal end (mW), computed the same way
compose_breath_video.py does it: a distal-referenced total-pressure drop times
the bulk volumetric flow. dE/dx is that curve differentiated with respect to
centerline distance, so it reads as local dissipated power per unit length
(mW/mm) -- i.e. WHERE the airway is resistive, not just how much has been lost
by a given station.

Reads csa_pe.csv (from csa_pressure_energy.py).
"""
import argparse, csv, os, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--pe-csv", required=True)
p.add_argument("--prof", default="")
p.add_argument("--region", default="Mouth")
p.add_argument("--out", required=True)
p.add_argument("--rho", type=float, default=1.2)
p.add_argument("--smooth", type=int, default=5,
               help="odd window for smoothing dE/dx (plane spacing is ~1 mm and "
                    "the raw derivative is noisy); 0 disables")
p.add_argument("--theme", default="dark", choices=["dark", "light"])
a = p.parse_args()

# ---- read per-plane pressure/energy, group by CSA FRAME ----
# Not by cfd_t_s: csa_pressure_energy maps frame -> cfd_idx = round(ms/25) - 1
# clamped at 0, so out_000000 and out_000025 share cfd index 0. Grouping on the
# CFD time therefore merges two frames, duplicating every plane and giving
# dx = 0 in the derivative. frame_ms is 1:1 with the CSA frame.
wframes = {}
for r in csv.DictReader(open(a.pe_csv)):
    if r.get("region") and r["region"] != a.region:
        continue
    wframes.setdefault(float(r["frame_ms"]) / 1000.0, []).append(r)
if not wframes:
    raise SystemExit("no rows for region %r in %s" % (a.region, a.pe_csv))

times = sorted(wframes)
curves = {}
for tw in times:
    rs = sorted(wframes[tw], key=lambda r: -float(r["arc_mm"]))   # distal first
    arc = np.array([float(r["arc_mm"]) for r in rs])
    x = arc.max() - arc                                           # 0 at distal end
    Q = np.array([float(r["Q_Lmin"]) for r in rs])
    tpm = np.array([float(r["TP_mdot_avg_Pa"]) for r in rs])
    tpa = np.array([float(r["TP_area_avg_Pa"]) for r in rs])
    tpp = np.where(np.abs(Q) < 2.0, tpa, tpm)                     # mass-flow avg unless Q ~ 0
    dP = tpp[0] - tpp                                             # distal-referenced drop (Pa)
    q = abs(np.median(Q)) / 60000.0                               # bulk flow (m3/s)
    E = dP * q * 1000.0                                           # mW
    o = np.argsort(x)
    x, E = x[o], E[o]
    # dE/dx via central differences on a possibly non-uniform grid
    dEdx = np.gradient(E, x)
    if a.smooth and a.smooth > 1 and len(dEdx) > a.smooth:
        w = a.smooth | 1                                          # force odd
        ker = np.ones(w) / w
        dEdx = np.convolve(dEdx, ker, mode="same")
    curves[tw] = dict(x=x, E=E, dEdx=dEdx)

# common grid for the min-max envelope
xmax = max(c["x"].max() for c in curves.values())
grid = np.linspace(0, xmax, 300)
Es, Ds = [], []
for tw in times:
    c = curves[tw]
    Es.append(np.interp(grid, c["x"], c["E"], left=np.nan, right=np.nan))
    Ds.append(np.interp(grid, c["x"], c["dEdx"], left=np.nan, right=np.nan))
Es, Ds = np.array(Es), np.array(Ds)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Elo, Ehi = np.nanmin(Es, axis=0), np.nanmax(Es, axis=0)
    Dlo, Dhi = np.nanmin(Ds, axis=0), np.nanmax(Ds, axis=0)

# representative instant: the frame with the largest peak |dE/dx|
peak_t = max(times, key=lambda t: np.nanmax(np.abs(curves[t]["dEdx"])))
pc = curves[peak_t]
pk_i = int(np.nanargmax(np.abs(pc["dEdx"])))

if a.theme == "light":
    TH = dict(bg="white", fg="#15181d", tick="#4a4a4a", title="#12305c", spine="0.5",
              legfc="#eef2f7", legec="0.6", grid="0.78", zero="0.55")
else:
    TH = dict(bg="#0f0f14", fg="0.92", tick="0.82", title="0.95", spine="0.35",
              legfc="#1b1b24", legec="0.45", grid=None, zero="0.45")
plt.rcParams.update({
    "font.weight": "bold", "font.size": 14,
    "axes.titleweight": "bold", "axes.labelweight": "bold",
    "axes.titlesize": 15.5, "axes.labelsize": 15,
    "xtick.labelsize": 13, "ytick.labelsize": 13,
    "text.color": TH["fg"], "axes.labelcolor": TH["fg"],
    "xtick.color": TH["tick"], "ytick.color": TH["tick"],
})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 10), facecolor=TH["bg"], sharex=True,
                               gridspec_kw=dict(hspace=0.18, left=0.135, right=0.97,
                                                top=0.93, bottom=0.08))
for ax in (ax1, ax2):
    ax.set_facecolor(TH["bg"])
    for sp in ax.spines.values():
        sp.set_color(TH["spine"])
    ax.grid(alpha=0.35, color=TH["grid"]) if TH["grid"] else ax.grid(alpha=0.16)

ax1.fill_between(grid, Elo, Ehi, color="#16a085", alpha=0.22, linewidth=0,
                 label="min–max over breath")
ax1.plot(pc["x"], pc["E"], color="#f1c40f", lw=3.0, label="t = %.3f s (peak dE/dx)" % peak_t)
ax1.axhline(0, color=TH["zero"], lw=0.8)
ax1.set_ylabel("cumulative WOB  E  (mW)")
ax1.set_title("Cumulative work of breathing along the airway — %s" % a.region,
              color=TH["title"])
ax1.legend(fontsize=13, loc="upper left", facecolor=TH["legfc"], edgecolor=TH["legec"],
           labelcolor=TH["fg"], framealpha=0.95)

ax2.fill_between(grid, Dlo, Dhi, color="#e67e22", alpha=0.22, linewidth=0,
                 label="min–max over breath")
ax2.plot(pc["x"], pc["dEdx"], color="#e74c3c", lw=3.0, label="t = %.3f s" % peak_t)
ax2.axhline(0, color=TH["zero"], lw=0.8)
ax2.plot(pc["x"][pk_i], pc["dEdx"][pk_i], "o", color="#ffffff", ms=11, zorder=6,
         markeredgecolor="#e74c3c", markeredgewidth=2.5)
ax2.annotate("peak loss\n%.1f mm, %.3g mW/mm" % (pc["x"][pk_i], pc["dEdx"][pk_i]),
             xy=(pc["x"][pk_i], pc["dEdx"][pk_i]), xytext=(-150, 20),
             textcoords="offset points", fontsize=13, color=TH["fg"],
             arrowprops=dict(arrowstyle="->", color=TH["fg"], lw=1.4))
ax2.set_ylabel("local dissipation  dE/dx  (mW/mm)")
ax2.set_xlabel("centerline length (mm)   —   0 = distal opening")
ax2.set_title("Where the energy is lost — spatial derivative of E", color=TH["title"])
ax2.legend(fontsize=13, loc="upper left", facecolor=TH["legfc"], edgecolor=TH["legec"],
           labelcolor=TH["fg"], framealpha=0.95)

os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
fig.savefig(a.out, dpi=300, facecolor=TH["bg"])
print("wrote", a.out)
print("frames=%d  planes/frame=%d  x: 0..%.1f mm" % (len(times), len(pc["x"]), xmax))
print("E     range: %.4g .. %.4g mW" % (np.nanmin(Es), np.nanmax(Es)))
print("dE/dx range: %.4g .. %.4g mW/mm" % (np.nanmin(Ds), np.nanmax(Ds)))
print("peak |dE/dx| at t=%.3f s, x=%.1f mm -> %.4g mW/mm"
      % (peak_t, pc["x"][pk_i], pc["dEdx"][pk_i]))
