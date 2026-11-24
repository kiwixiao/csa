#!/usr/bin/env python3
"""CSA breath video with a per-region centerline axis.

Each region is plotted against its OWN branch centerline length, flipped so 0 is
the anatomical opening (lips / nostril / pharyngeal end). VMTK builds centerlines
source->target starting at the tracheal outlet, so raw arc length grows inward;
flipping (x = arc.max() - arc) puts 0 at the opening for every branch.

That makes the mouth overlayable with the nasal passages: previously the mouth was
drawn on a nose-referenced "distance from mouth" axis, which aligned it against the
trachea and meant nothing. The axis is therefore labelled "centerline length",
not "distance to nose".

Panels (all driven by CSA + the measured flow profile; no CFD):
  LEFT  : airway wall + moving CSA planes / plane locator
  RIGHT : flow(t) with phase dot / CSA vs centerline length, all regions overlaid
"""
import argparse, csv, os, sys, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

p = argparse.ArgumentParser()
p.add_argument("--csa-csv", required=True)
p.add_argument("--prof", required=True)
p.add_argument("--planes-dir", required=True)
p.add_argument("--ref", default="")
p.add_argument("--out", required=True)
p.add_argument("--regions", default="DescendingAirway,Mouth,LeftNose,RightNose")
p.add_argument("--independent", default="Mouth",
               help="comma-separated regions measured on their OWN centerline. "
                    "Everything else shares the global arc so anatomically "
                    "continuous branches (nose -> nasopharynx -> trachea) line up.")
p.add_argument("--xaxis", default="arc", choices=["arc", "fixed", "index"],
               help="'arc' = per-frame centerline distance (x slides with bulk motion); "
                    "'fixed' = each plane pinned at its MEAN mm position (no drift, true spacing); "
                    "'index' = signed INTEGER plane position (no drift, evenly spaced; 0=nasopharynx)")
p.add_argument("--rho", type=float, default=1.2)
p.add_argument("--dpi", type=int, default=100)
p.add_argument("--theme", default="dark", choices=["dark", "light"])
a = p.parse_args()
os.makedirs(a.out, exist_ok=True)
REGIONS = [r.strip() for r in a.regions.split(",") if r.strip()]

RCOL = {"DescendingAirway": "#4aa3df", "Mouth": "#f39c12",
        "LeftNose": "#9b59b6", "RightNose": "#2ecc71"}

# ---- flow profile -> L/min ----
tp, qp = [], []
for r in csv.reader(open(a.prof)):
    try:
        tp.append(float(r[0])); qp.append(float(r[1]) / a.rho * 60000.0)
    except Exception:
        pass
tp, qp = np.array(tp), np.array(qp)
# light smooth so the prescribed waveform isn't staircased (display only)
try:
    from scipy.signal import savgol_filter
    if len(qp) >= 9:
        qp = savgol_filter(qp, 9, 2)
except Exception:
    pass

# ---- CSA keyed by ms then region ----
raw = {}
for r in csv.DictReader(open(a.csa_csv)):
    if r["is_valid"].lower() != "true":
        continue
    reg = r.get("region")
    if reg not in REGIONS:
        continue
    ms = int(r["frame_name"].split("_")[1])
    raw.setdefault(ms, {}).setdefault(reg, []).append(r)
if not raw:
    print("ERROR: no rows for regions %s" % REGIONS); sys.exit(1)

ms_sorted = sorted(raw)
NF = len(ms_sorted)
INDEP = {r.strip() for r in a.independent.split(",") if r.strip()}

# Which arc column each region uses:
#   connected group (nose + descending airway) -> global arc (one continuous line)
#   independent group (mouth) -> its own branch arc (separate centerline)
ARC_COL = {reg: ("arc_length_branch_mm" if reg in INDEP else "arc_length_mm")
           for reg in REGIONS}

# CONFLUENCE axis: 0 = nasopharynx (= max arc of the DescendingAirway trunk, where
# the nasal path merges). Computed ONCE so the axis does not drift frame-to-frame.
#   connected (nose + DA): s = conf - arc  -> nasopharynx 0, trachea +, nose - (inlets)
#   independent (mouth): kept on its OWN arc, anchored at its TRUE junction on the
#     trunk (the oropharynx, ~+38 mm downstream of the nasopharynx) and extending
#     upstream (negative) toward the lips -- NOT fake-aligned to the nasal path.
CONF_REG = "DescendingAirway"
_conf = 0.0
_ref_ms = None
_indep_min = {}
for ms in ms_sorted:
    if CONF_REG in raw[ms]:
        _conf = max(_conf, max(float(x["arc_length_mm"]) for x in raw[ms][CONF_REG]))
        if _ref_ms is None:
            _ref_ms = ms
    for reg, rs in raw[ms].items():
        if reg in INDEP:
            _indep_min[reg] = min(_indep_min.get(reg, 1e18),
                                  min(float(x[ARC_COL[reg]]) for x in rs))

def _trunk_dist_at_z(zq):     # map a z-level to distance along the DA trunk (conf=0)
    rs = raw[_ref_ms][CONF_REG]
    zs = np.array([float(x["centroid_z"]) for x in rs])
    ar = np.array([float(x["arc_length_mm"]) for x in rs])
    return _conf - ar[int(np.argmin(np.abs(zs - zq)))]

JUNC = {}
for reg in INDEP:
    if reg not in REGIONS:
        continue
    rs = raw[_ref_ms].get(reg) or next(v[reg] for v in raw.values() if reg in v)
    phar = min(rs, key=lambda x: float(x[ARC_COL[reg]]))      # pharyngeal end = min branch arc
    JUNC[reg] = _trunk_dist_at_z(float(phar["centroid_z"]))
print("confluence(nasopharynx)=0 ; trunk max=%.0f mm ; mouth junction(s)=%s" %
      (_conf, {k: round(v) for k, v in JUNC.items()}), flush=True)

# --xaxis fixed/index: freeze each plane's x so the centerline's bulk motion doesn't
# slide the curve; only CSA (height) moves.
#   fixed = each plane's MEAN centerline distance (mm, true anatomical spacing)
#   index = signed INTEGER plane position (nasopharynx 0, trachea +, nose/mouth -)
_MEAN_ARC = {}
IDX_POS = {}
if a.xaxis in ("fixed", "index"):
    _acc = {}
    for ms in ms_sorted:
        for reg, rs in raw[ms].items():
            for x in rs:
                _acc.setdefault((reg, int(x["plane_index"])), []).append(float(x[ARC_COL[reg]]))
    _MEAN_ARC = {k: float(np.mean(v)) for k, v in _acc.items()}
if a.xaxis == "index":
    # Each region gets its OWN ordinal (the two nostrils are PARALLEL, not in series).
    def _ranked(reg):
        return sorted([(v, pidx) for (r, pidx), v in _MEAN_ARC.items() if r == reg], key=lambda t: t[0])
    # DA trunk: nasopharynx (max arc) = 0, trachea (min arc) = +(N-1)
    da = _ranked(CONF_REG); n_da = len(da)
    for rank, (ma, pidx) in enumerate(da):
        IDX_POS[(CONF_REG, pidx)] = float((n_da - 1) - rank)
    # nasal branches: own ordinal, negative, each starting just above the nasopharynx (parallel)
    for reg in REGIONS:
        if reg == CONF_REG or reg in INDEP:
            continue
        for rank, (ma, pidx) in enumerate(_ranked(reg)):
            IDX_POS[(reg, pidx)] = float(-(rank + 1))
    # mouth: anchor its pharyngeal end at the junction index (nearest DA plane to the +38mm junction)
    for reg in INDEP:
        if reg not in REGIONS:
            continue
        junc_idx = min(((_conf - _MEAN_ARC[(CONF_REG, p)], IDX_POS[(CONF_REG, p)])
                        for (r, p) in _MEAN_ARC if r == CONF_REG),
                       key=lambda t: abs(t[0] - JUNC[reg]))[1]
        for rank, (ma, pidx) in enumerate(_ranked(reg)):
            IDX_POS[(reg, pidx)] = junc_idx - float(rank)

per = {}          # per[ms][reg] = dict(s=..., area=...)
smax = {}
smin = {}
for ms in ms_sorted:
    per[ms] = {}
    for reg, rs in raw[ms].items():
        area = np.array([float(x["area_mm2"]) for x in rs])
        if a.xaxis == "index":
            s = np.array([IDX_POS[(reg, int(x["plane_index"]))] for x in rs])
        else:
            arc = np.array([(_MEAN_ARC[(reg, int(x["plane_index"]))] if a.xaxis == "fixed"
                             else float(x[ARC_COL[reg]])) for x in rs])
            s = (JUNC[reg] - (arc - _indep_min[reg])) if reg in INDEP else (_conf - arc)
        o = np.argsort(s)
        per[ms][reg] = dict(s=s[o], area=area[o])
        smax[reg] = max(smax.get(reg, -1e18), float(s.max()))
        smin[reg] = min(smin.get(reg, 1e18), float(s.min()))

# per-region min/max envelope over the breath, on a common per-region grid
NG = 200
grid, env = {}, {}
for reg in REGIONS:
    if reg not in smax:
        continue
    g = np.linspace(smin[reg], smax[reg], NG)
    stack = []
    for ms in ms_sorted:
        if reg not in per[ms]:
            continue
        d = per[ms][reg]
        stack.append(np.interp(g, d["s"], d["area"], left=np.nan, right=np.nan))
    stack = np.array(stack)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        env[reg] = (np.nanmin(stack, axis=0), np.nanmax(stack, axis=0))
    grid[reg] = g
XMAX = max(smax.values())
# headroom so the legend never sits on top of the curves
YMAX = max(float(np.nanmax(e[1])) for e in env.values()) * 1.42

if a.theme == "light":
    TH = dict(bg="white", fg="#15181d", tick="#4a4a4a", title="#12305c", spine="0.5",
              legfc="#eef2f7", legec="0.6", grid="0.78", flow="0.32", zero="0.55")
else:
    TH = dict(bg="#0f0f14", fg="0.92", tick="0.82", title="0.95", spine="0.35",
              legfc="#1b1b24", legec="0.45", grid=None, flow="0.6", zero="0.4")
plt.rcParams.update({
    "font.weight": "bold", "font.size": 16,
    "axes.titleweight": "bold", "axes.labelweight": "bold",
    "axes.titlesize": 17, "axes.labelsize": 16.5,
    "xtick.labelsize": 14, "ytick.labelsize": 14,
    "text.color": TH["fg"], "axes.labelcolor": TH["fg"],
    "xtick.color": TH["tick"], "ytick.color": TH["tick"],
})
BG = TH["bg"]
REF = mpimg.imread(a.ref) if a.ref and os.path.exists(a.ref) else None


def style_axes(ax):
    for sp in ax.spines.values():
        sp.set_color(TH["spine"])
    ax.grid(alpha=0.35, color=TH["grid"]) if TH["grid"] else ax.grid(alpha=0.16)


def compose(k):
    ms = ms_sorted[k]
    t = ms / 1000.0
    Bi = mpimg.imread("%s/frame_%04d.png" % (a.planes_dir, k)) \
        if os.path.exists("%s/frame_%04d.png" % (a.planes_dir, k)) else None
    qnow = float(np.interp(t, tp, qp))
    # Sign convention verified against the exported caps (postproc/outlet_flow/
    # inlet_outlet_massflow.csv): a POSITIVE prescribed mass flow is air entering the
    # domain at BronchiOut and leaving at MaskInlet, i.e. EXHALATION. ENT003 exhales
    # first and inhales after t~1.36 s.
    phase = "INHALE" if qnow < -0.3 else ("EXHALE" if qnow > 0.3 else "transition")

    fig = plt.figure(figsize=(15.0, 9.0), facecolor=BG)
    gs = fig.add_gridspec(2, 12, height_ratios=[1, 1.15],
                          hspace=0.30, wspace=0.55, left=0.045, right=0.985,
                          top=0.935, bottom=0.085)

    axB = fig.add_subplot(gs[0, 0:4])
    if Bi is not None: axB.imshow(Bi)
    axB.axis("off")
    axB.set_title("airway wall + CSA planes (moving)", fontsize=16, color=TH["title"])

    axL = fig.add_subplot(gs[1, 0:4])
    if REF is not None: axL.imshow(REF)
    axL.axis("off")
    axL.set_title("plane locator", fontsize=16, color=TH["title"])

    axf = fig.add_subplot(gs[0, 4:12]); axf.set_facecolor(BG)
    axf.plot(tp, qp, color=TH["flow"], lw=2.0)
    axf.axhline(0, color=TH["zero"], lw=0.7)
    axf.axvline(t, color="#e74c3c", lw=1.1, alpha=0.55)
    axf.plot(t, qnow, "o", color="#e74c3c", ms=13, zorder=6)
    axf.set_xlim(tp.min(), tp.max())
    axf.set_ylabel("flow (L/min)"); axf.set_xlabel("time (s)")
    axf.set_title("breathing cycle    t = %.2f s    —    %s" % (t, phase),
                  color=TH["title"], fontsize=17)
    style_axes(axf)

    axC = fig.add_subplot(gs[1, 4:12]); axC.set_facecolor(BG)
    for reg in REGIONS:
        if reg not in per[ms]:
            continue
        c = RCOL.get(reg, "#aaaaaa")
        lo, hi = env[reg]
        axC.fill_between(grid[reg], lo, hi, color=c, alpha=0.16, linewidth=0)
        d = per[ms][reg]
        lbl = reg + ("  (own branch, joins +%d mm)" % round(JUNC[reg]) if reg in INDEP else "")
        axC.plot(d["s"], d["area"], color=c, lw=3.0, label=lbl,
                 ls=("--" if reg in INDEP else "-"), solid_capstyle="round")
    axC.axvline(0, color=TH["zero"], lw=1.4, ls=(0, (5, 3)), alpha=0.9)
    axC.text(2.5, YMAX * 0.965, "nasopharynx\n(confluence)", ha="left", va="top",
             fontsize=11.5, color=TH["title"], fontweight="bold")
    axC.annotate("inlets", xy=(min(smin.values()) + 4, YMAX * 0.06), fontsize=12,
                 color=TH["fg"], ha="left", fontweight="bold", alpha=0.7)
    axC.annotate("→ pharynx · throat · trachea", xy=(XMAX * 0.42, YMAX * 0.06), fontsize=12,
                 color=TH["fg"], ha="left", fontweight="bold", alpha=0.7)
    axC.set_xlabel(("plane index from nasopharynx   —   0 = confluence · nose/mouth < 0 · pharynx→trachea > 0  (evenly spaced, no drift)"
                    if a.xaxis == "index" else
                    "distance from nasopharynx (mm)   —   PLANE-FIXED (each plane at its mean position; only CSA moves)"
                    if a.xaxis == "fixed" else
                    "distance from nasopharynx (mm)   —   0 = confluence · inlets (nose/mouth) < 0 · pharynx→trachea > 0"))
    axC.set_ylabel("CSA (mm²)")
    axC.set_xlim(min(smin.values()) - 5, XMAX + 5); axC.set_ylim(0, YMAX)
    leg = axC.legend(fontsize=14, loc="upper right", ncol=2,
                     facecolor=TH["legfc"], edgecolor=TH["legec"],
                     labelcolor=TH["fg"], framealpha=0.95, borderpad=0.7,
                     labelspacing=0.5, handlelength=2.4, title="region")
    leg.get_title().set_color(TH["fg"])
    leg.get_title().set_fontsize(14)
    leg.get_title().set_fontweight("bold")
    style_axes(axC)
    axC.set_title(("cross-sectional area vs PLANE INDEX (fixed integer position; pure area change)"
                   if a.xaxis == "index" else
                   "cross-sectional area — PLANE-FIXED axis (no bulk drift; pure area change)"
                   if a.xaxis == "fixed" else
                   "cross-sectional area vs centerline length   (shaded = min–max over breath)"),
                  color=TH["title"], fontsize=16.5)

    fig.savefig("%s/frame_%04d.png" % (a.out, k), dpi=a.dpi, facecolor=BG)
    plt.close(fig)


print("frames=%d  regions=%s" % (NF, [r for r in REGIONS if r in smax]), flush=True)
for reg in REGIONS:
    if reg in smax:
        print("   %-18s centerline length %.1f .. %.1f mm%s" % (reg, smin[reg], smax[reg], "  (own centerline)" if reg in INDEP else ""), flush=True)
for k in range(NF):
    compose(k)
    if k % 20 == 0:
        print("  composited %d/%d" % (k, NF), flush=True)
print("DONE -> %s" % a.out, flush=True)
