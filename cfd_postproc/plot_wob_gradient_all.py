#!/usr/bin/env python3
"""Cumulative WOB E(x) and local dissipation dE/dx for ALL regions overlaid.

Axis convention matches compose_csa_video.py:
  * connected group (nose + descending airway) shares one flip reference on the
    global arc, so the branches stay joined to the trunk they feed
  * independent regions (mouth) are flipped against their own maximum, since
    they sit on a separate centerline

dE/dx is local dissipated power per unit length (mW/mm) -- where the airway is
resistive, as opposed to the running total E(x).
"""
import argparse, csv, glob, os, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--pe-glob", required=True, help="glob for csa_pe*.csv files")
p.add_argument("--out", required=True)
p.add_argument("--independent", default="Mouth")
p.add_argument("--exclude", default="",
               help="comma-separated regions to leave out. Use for branches whose "
                    "energy loss is not a whole-airway number -- e.g. ENT003's Mouth "
                    "branch, which ends at the oropharynx and carries only ~1/3 of the "
                    "flow, so its dP*Q is a partial-route figure not comparable to the "
                    "nostril->trachea routes.")
p.add_argument("--smooth", type=int, default=5)
p.add_argument("--theme", default="dark", choices=["dark", "light"])
a = p.parse_args()
INDEP = {r.strip() for r in a.independent.split(",") if r.strip()}

RCOL = {"DescendingAirway": "#4aa3df", "Mouth": "#f39c12",
        "LeftNose": "#9b59b6", "RightNose": "#2ecc71"}

# ---- load every region, grouped by CSA frame (frame_ms is 1:1 with the frame;
# cfd_t_s is not -- out_000000 and out_000025 share cfd index 0) ----
EXCLUDE = {x.strip() for x in a.exclude.split(",") if x.strip()}
if EXCLUDE:
    print("excluding regions:", ", ".join(sorted(EXCLUDE)))
byreg = {}
for path in sorted(glob.glob(a.pe_glob)):
    for r in csv.DictReader(open(path)):
        reg = r.get("region")
        if not reg or reg in EXCLUDE:
            continue
        byreg.setdefault(reg, {}).setdefault(float(r["frame_ms"]) / 1000.0, []).append(r)
if not byreg:
    raise SystemExit("no rows matched %s" % a.pe_glob)

# shared flip reference across the connected group
conn_max = 0.0
indep_max = {}
for reg, frames in byreg.items():
    m = max(float(x["arc_mm"]) for rs in frames.values() for x in rs)
    if reg in INDEP:
        indep_max[reg] = m
    else:
        conn_max = max(conn_max, m)
REF = {reg: (indep_max[reg] if reg in INDEP else conn_max) for reg in byreg}

curves = {}
for reg, frames in byreg.items():
    curves[reg] = {}
    for tw, rs in frames.items():
        rs = sorted(rs, key=lambda r: -float(r["arc_mm"]))
        arc = np.array([float(r["arc_mm"]) for r in rs])
        x = REF[reg] - arc
        Q = np.array([float(r["Q_Lmin"]) for r in rs])
        tpm = np.array([float(r["TP_mdot_avg_Pa"]) for r in rs])
        tpa = np.array([float(r["TP_area_avg_Pa"]) for r in rs])
        tpp = np.where(np.abs(Q) < 2.0, tpa, tpm)
        E = (tpp[0] - tpp) * (abs(np.median(Q)) / 60000.0) * 1000.0     # mW
        o = np.argsort(x)
        x, E = x[o], E[o]
        keep = np.concatenate([[True], np.diff(x) > 1e-9])              # drop dup stations
        x, E = x[keep], E[keep]
        if len(x) < 3:
            continue
        d = np.gradient(E, x)
        if a.smooth > 1 and len(d) > a.smooth:
            w = a.smooth | 1
            d = np.convolve(d, np.ones(w) / w, mode="same")
        curves[reg][tw] = dict(x=x, E=E, d=d)

# envelopes per region
env = {}
for reg, cs in curves.items():
    lo_x = min(c["x"].min() for c in cs.values())
    hi_x = max(c["x"].max() for c in cs.values())
    g = np.linspace(lo_x, hi_x, 300)
    Es = np.array([np.interp(g, c["x"], c["E"], left=np.nan, right=np.nan) for c in cs.values()])
    Ds = np.array([np.interp(g, c["x"], c["d"], left=np.nan, right=np.nan) for c in cs.values()])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        env[reg] = dict(g=g, Elo=np.nanmin(Es, 0), Ehi=np.nanmax(Es, 0),
                        Dlo=np.nanmin(Ds, 0), Dhi=np.nanmax(Ds, 0))

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
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10.5), facecolor=TH["bg"], sharex=True,
                               gridspec_kw=dict(hspace=0.20, left=0.125, right=0.975,
                                                top=0.93, bottom=0.085))
for ax in (ax1, ax2):
    ax.set_facecolor(TH["bg"])
    for sp in ax.spines.values():
        sp.set_color(TH["spine"])
    ax.grid(alpha=0.35, color=TH["grid"]) if TH["grid"] else ax.grid(alpha=0.16)
    ax.axhline(0, color=TH["zero"], lw=0.8)

order = [r for r in ("DescendingAirway", "Mouth", "LeftNose", "RightNose") if r in curves]
peaks = []
for reg in order:
    c = RCOL.get(reg, "#aaa")
    e = env[reg]
    dash = "--" if reg in INDEP else "-"
    lbl = reg + ("  (own centerline)" if reg in INDEP else "")
    ax1.fill_between(e["g"], e["Elo"], e["Ehi"], color=c, alpha=0.16, linewidth=0)
    ax1.plot(e["g"], e["Ehi"], color=c, lw=2.6, ls=dash, label=lbl)
    ax2.fill_between(e["g"], e["Dlo"], e["Dhi"], color=c, alpha=0.16, linewidth=0)
    ax2.plot(e["g"], e["Dhi"], color=c, lw=2.6, ls=dash, label=lbl)
    i = int(np.nanargmax(np.abs(e["Dhi"])))
    peaks.append((reg, e["g"][i], e["Dhi"][i]))
    ax2.plot(e["g"][i], e["Dhi"][i], "o", color="#fff", ms=9, zorder=6,
             markeredgecolor=c, markeredgewidth=2.4)

ax1.set_ylabel("cumulative WOB  E  (mW)")
ax1.set_title("Cumulative work of breathing along the airway — all regions", color=TH["title"])
ax1.legend(fontsize=13, loc="upper left", ncol=2, facecolor=TH["legfc"],
           edgecolor=TH["legec"], labelcolor=TH["fg"], framealpha=0.95)
ax2.set_ylabel("local dissipation  dE/dx  (mW/mm)")
ax2.set_xlabel("centerline length (mm)   —   0 = nostril; mouth from its own opening")
ax2.set_title("Where the energy is lost — dE/dx  (envelope max over the breath; dots = peaks)",
              color=TH["title"])
ax2.legend(fontsize=13, loc="upper left", ncol=2, facecolor=TH["legfc"],
           edgecolor=TH["legec"], labelcolor=TH["fg"], framealpha=0.95)

os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
fig.savefig(a.out, dpi=300, facecolor=TH["bg"])
print("wrote", a.out)
for reg in order:
    cs = curves[reg]
    print("  %-18s frames=%3d  x: %.1f..%.1f mm" %
          (reg, len(cs), env[reg]["g"].min(), env[reg]["g"].max()))
print("\npeak dE/dx per region:")
for reg, x, v in sorted(peaks, key=lambda t: -abs(t[2])):
    print("  %-18s %8.1f mm   %.4g mW/mm" % (reg, x, v))
