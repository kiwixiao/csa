"""CSA along the centerline, baseline vs edited geometry, both from the registration (same FFDs), same planes
(compare/run_samplane.sh output). Baseline = the slicer's CSA table. Edited = baseline x (edited / baseline re-cut on the
same plane): the same cutting method on both deformed surfaces, so method differences cancel.
Usage: plot_samplane.py --csv samplane_base_vs_edit.csv --prof FlowProfile.csv --out fig.png
           [--edit-zone 143:170] [--times 0,450,1500] [--labels baseline,edited] [--title "..."]
--times default: t0, peak inhale, peak exhale (from the flow profile; inhale < 0). --edit-zone (mm from the nostril)
shades the edit and sets the 'narrowest in the edit' numbers; default = where the ratio exceeds 1.05.
"""
import argparse, csv, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
ap = argparse.ArgumentParser(); ap.add_argument("--csv", required=True); ap.add_argument("--prof", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--edit-zone", default=None); ap.add_argument("--times", default=None)
ap.add_argument("--labels", default="baseline,edited geometry")
ap.add_argument("--title", default="CSA along the centerline, baseline vs edited geometry, both from the registration")
ap.add_argument("--rho", type=float, default=1.2090, help="kg/m3: flow profile kg/s -> L/min")
args = ap.parse_args()
LB, LA = args.labels.split(","); CB, CA = "0.15", "#2471a3"
d = pd.read_csv(args.csv).rename(columns={"area_aryt": "area_edit", "n_aryt": "n_edit"})   # older files: aryt
d["s"] = d.groupby("ms").arc_mm.transform("max") - d.arc_mm
d["ratio"] = d.area_edit / d.area_base; d["A_b"] = d.area_table; d["A_a"] = d.area_table * d.ratio
d = d.sort_values(["ms", "s"])
if args.edit_zone:
    EDIT = tuple(float(x) for x in args.edit_zone.split(":"))
else:
    hi = d[d.ratio > 1.05].s; EDIT = (hi.quantile(0.02), hi.quantile(0.98)) if len(hi) else (0, 0)
tq, qq = [], []
for r in csv.reader(open(args.prof)):
    try: tq.append(float(r[0])); qq.append(float(r[1]) / args.rho * 60000)
    except ValueError: pass
tq, qq = np.array(tq), np.array(qq)
MSS = np.array(sorted(d.ms.unique()))
snap = lambda t: int(MSS[np.argmin(np.abs(MSS - 1000 * t))])
kin, kex = snap(tq[np.argmin(qq)]), snap(tq[np.argmax(qq)])
TMS = [int(x) for x in args.times.split(",")] if args.times else [int(MSS[0]), kin, kex]
TP = [(m, "t = %.2f s%s" % (m / 1000, "  peak inhale" if m == kin else "  peak exhale" if m == kex else "  (frame 1)" if m == MSS[0] else "")) for m in TMS]

fig = plt.figure(figsize=(18, 11.5))
fig.text(0.035, 0.975, args.title, fontsize=17, weight="bold", va="top")
fig.text(0.035, 0.94, "Frame-1 (t0) surface of each geometry moved through the breath with the same registration (MIRTK FFDs), then cut with the SAME "
         "centerline planes (the baseline's), so the only difference is the geometry.\nBaseline = the slicer's CSA.  Edited = baseline x (edited / baseline cut on that plane).  "
         "s = distance from the left nostril along the left nose -> trachea path.", fontsize=10, va="top", color="0.25")
for j, (ms, ttl) in enumerate(TP):
    q = d[d.ms == ms]
    a = fig.add_axes([0.045 + j * 0.322, 0.47, 0.29, 0.39]); b = fig.add_axes([0.045 + j * 0.322, 0.33, 0.29, 0.12], sharex=a)
    a.plot(q.s, q.A_b, color=CB, lw=2.4, label=LB)
    a.plot(q.s, q.A_a, color=CA, lw=1.6, label=LA)
    a.fill_between(q.s, q.A_b, q.A_a, where=q.A_a > q.A_b, color=CA, alpha=0.25, lw=0)
    b.plot(q.s, q.ratio, color=CA, lw=1.5); b.axhline(1, color="k", lw=0.8)
    for x in (a, b): x.axvspan(*EDIT, color="#fff0c8", zorder=0); x.grid(alpha=0.3); x.set_xlim(0, 265)
    e = q[(q.s >= EDIT[0]) & (q.s <= EDIT[1])]; k = e.A_b.idxmin()
    a.set_title(ttl, fontsize=12.5, weight="bold"); a.set_ylim(0, 420); plt.setp(a.get_xticklabels(), visible=False)
    a.text(np.mean(EDIT), 0.02, "edit\nzone", transform=a.get_xaxis_transform(), ha="center", va="bottom", fontsize=9)
    a.text(0.98, 0.97, "narrowest in the edit:\nbaseline %.0f mm²  ->  edited %.0f mm²  (x%.2f)" % (e.A_b.min(), e.loc[k, "A_a"], e.loc[k, "A_a"] / e.A_b.min()),
           transform=a.transAxes, va="top", ha="right", fontsize=9.5, bbox=dict(fc="w", ec="0.7", lw=0.6))
    b.set_ylim(0.85, 1.75); b.set_xlabel("distance from nostril s (mm)")
    if j == 0:
        a.set_ylabel("CSA (mm²)"); b.set_ylabel("edited / baseline"); a.legend(loc="upper left", fontsize=10)
# whole breath: ratio map + narrowest CSA in the edit vs time
h = fig.add_axes([0.045, 0.06, 0.45, 0.2])
d["sb"] = (d.s // 2) * 2
H = d.pivot_table(index="ms", columns="sb", values="ratio", aggfunc="median")
im = h.pcolormesh(H.columns.values, H.index.values / 1000, H.values, cmap="RdBu_r", norm=LogNorm(0.67, 1.5), shading="nearest")
h.axvline(EDIT[0], color="k", ls="--", lw=0.8); h.axvline(EDIT[1], color="k", ls="--", lw=0.8); h.set_xlim(0, 265)
h.set_xlabel("distance from nostril s (mm)"); h.set_ylabel("time (s)"); h.set_title("edited / baseline over the whole breath (white = same)", fontsize=11.5, weight="bold")
cb = fig.colorbar(im, ax=h, pad=0.01, fraction=0.03); cb.ax.minorticks_off(); cb.set_ticks([0.67, 0.8, 1, 1.25, 1.5]); cb.set_ticklabels(["0.67", "0.8", "1", "1.25", "1.5"])
g = fig.add_axes([0.58, 0.06, 0.39, 0.2])
e = d[(d.s >= EDIT[0]) & (d.s <= EDIT[1])].groupby("ms")
t = np.array(sorted(d.ms.unique())) / 1000
g.plot(t, e.A_b.min(), color=CB, lw=2.2, label=LB); g.plot(t, e.A_a.min(), color=CA, lw=1.8, label=LA)
g.set_ylabel("narrowest CSA in the\nedit zone (mm²)"); g.set_xlabel("time (s)"); g.grid(alpha=0.3); g.set_xlim(0, 2.07); g.set_ylim(0, None)
g.axvspan(0, 0.98, color="#e8f1fa", zorder=0); g.axvspan(0.98, 2.07, color="#faece8", zorder=0)
g.text(0.49, 0.04, "INHALE", transform=g.get_xaxis_transform(), ha="center", fontsize=9, color="#1b4f72"); g.text(1.52, 0.04, "EXHALE", transform=g.get_xaxis_transform(), ha="center", fontsize=9, color="#7b241c")
g.legend(fontsize=9, loc="lower left", bbox_to_anchor=(0, 0.1)); g.set_title("Edit zone: narrowest cross-section through the breath", fontsize=11.5, weight="bold")
out = args.out; fig.savefig(out, dpi=120); print(out, "edit zone s = %.0f-%.0f mm" % EDIT)
for lo, hi, nm in ((0, 140, "nose+pharynx"), (175, 270, "below edit"), EDIT + ("edit",)):
    q = d[(d.s >= lo) & (d.s <= hi)].ratio; print("%-12s ratio median %.3f  p5 %.3f  p95 %.3f  max %.2f" % (nm, q.median(), q.quantile(.05), q.quantile(.95), q.max()))
m = pd.DataFrame(dict(b=e.A_b.min(), a=e.A_a.min())); m["r"] = m.a / m.b; print("narrowest in edit: baseline %.0f-%.0f, edited %.0f-%.0f mm2, ratio %.2f-%.2f" % (m.b.min(), m.b.max(), m.a.min(), m.a.max(), m.r.min(), m.r.max()))
