#!/usr/bin/env python3
"""Assemble the synchronized multi-panel breath video frames (3-row x 2-col):
  LEFT  col : airway STL+CSA planes (moving) / anatomical LOCATOR / velocity render
  RIGHT col : flow(t)+phase dot / CSA vs arc / cumulative ENERGY-loss(WOB, mW) vs arc

Ported from /storage4/ENT002/compose_breath_video.py. Differences:
  * all paths are arguments instead of module constants
  * NF (frame count) and NPL (planes per frame) are derived from the data rather
    than hardcoded to ENT002's 135/224
  * --region selects one CSA region, because ENT002 was single-region (Mouth)
    while ENT003's septum path produces four regions with overlapping arc ranges;
    mixing them would make the arc axis meaningless
  * the velocity colourbar range is taken from the actual rendered range
    (--vmax) instead of ENT002's fixed 0-8 m/s
"""
import argparse, csv, os, sys, warnings

import csa_paths as CP
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

p = argparse.ArgumentParser()
p.add_argument("--csa-csv", required=True)
p.add_argument("--prof", required=True, help="FlowProfile.csv")
p.add_argument("--wob-csv", required=True,
               help="directory holding the per-region csa_pe_<Region>.csv files")
p.add_argument("--planes-dir", required=True, help="stl+planes render frames")
p.add_argument("--vel-dir", required=True, help="velocity render frames")
p.add_argument("--ref", default="", help="static locator PNG (optional)")
p.add_argument("--out", required=True)
p.add_argument("--region", default="LeftNose",
               help="csa_paths path name: LeftNose | RightNose | Mouth | DescendingAirway")
p.add_argument("--frames", default=None,
               help="index range a:b (0-based, b exclusive) -- render a single frame for the PDF")
p.add_argument("--vmax", type=float, default=None, help="velocity colourbar max (default: auto)")
p.add_argument("--dt", type=float, default=0.025, help="s per CSA frame")
p.add_argument("--rho", type=float, default=1.2)
p.add_argument("--dpi", type=int, default=92)
p.add_argument("--theme", default="dark", choices=["dark", "light"])
p.add_argument("--layout", default="stacked", choices=["stacked", "columns", "report"],
               help="stacked: renders in a left column (ENT002/003); columns: renders across the top, each over its plot; "
                    "report: PDF page - locator / loss-coloured centerline / velocity, the plots, WOB per segment")
p.add_argument("--wob-ref", default="", help="report: centerline coloured by local loss (make_centerline_ref.py --color-csv)")
p.add_argument("--bins-csv", default="", help="report: WOB per segment (wob_regions.py --bin-mm)")
p.add_argument("--bins-power-csv", default="", help="report video: segment power vs time (wob_regions.py --series) -> "
               "bars fill up as the breath goes, the whole-breath total drawn as an outline")
p.add_argument("--wob-ref-dir", default="", help="report video: per-frame loss-coloured centerline images frame_%%04d.png")
p.add_argument("--pdf-out", default="", help="also save the (single) frame as this PDF")
p.add_argument("--page-title", default="", help="report: page title, ' — t = <t> s' is appended")
p.add_argument("--prox-label", default=None, help="name for the opening in axis labels instead of e.g. 'left nostril'")
p.add_argument("--mark-step", type=float, default=0,
               help="distance marks every N mm on the plots along the centerline (match the locator labels / WOB segments)")
p.add_argument("--figsize", default=None, help="W,H inches (default 12,13 stacked / 19.2,10.8 columns)")
p.add_argument("--font-scale", type=float, default=1.0)
p.add_argument("--csa-ylim", default=None,
               help="fixed CSA y-axis as 'lo,hi' in mm^2. Needed when two runs are shown "
                    "side by side: auto-scaling gives each its own limits and the curves "
                    "stop being visually comparable.")
p.add_argument("--energy-ylim", default=None,
               help="fixed energy-loss y-axis as 'lo,hi' in mW, same reasoning.")
a = p.parse_args()
os.makedirs(a.out, exist_ok=True)

# ---- flow profile -> L/min ----
tp, qp = [], []
for r in csv.reader(open(a.prof)):
    try:
        tp.append(float(r[0])); qp.append(float(r[1]) / a.rho * 60000.0)
    except Exception:
        pass
tp, qp = np.array(tp), np.array(qp)
# light smooth so the prescribed waveform isn't staircased (display only; broad
# inhale/exhale peaks are ~0.3 s wide, a w9 window leaves them intact)
try:
    from scipy.signal import savgol_filter
    if len(qp) >= 9:
        qp = savgol_filter(qp, 9, 2)
except Exception:
    pass

# ---- CSA anatomy keyed by ms (every region on the path) ----
PATH_REGS = set(CP.regions_of(a.region))
frames = {}
for r in csv.DictReader(open(a.csa_csv)):
    if r["is_valid"].lower() != "true" or r.get("region") not in PATH_REGS:
        continue
    frames.setdefault(int(r["frame_name"].split("_")[1]), []).append(r)
if not frames:
    print("ERROR: no rows for path %r" % a.region); sys.exit(1)

ms_sorted = sorted(frames)
NF = len(ms_sorted)
NPL = max(len(v) for v in frames.values())
print("path=%s (%s)  frames=%d  planes/frame<=%d"
      % (a.region, "+".join(CP.regions_of(a.region)), NF, NPL), flush=True)

per, smax = {}, 0.0
for ms in ms_sorted:
    rs = sorted(frames[ms], key=lambda r: -float(r["arc_length_mm"]))     # distal first
    arc = np.array([float(r["arc_length_mm"]) for r in rs])
    area = np.array([float(r["area_mm2"]) for r in rs])
    s = arc.max() - arc
    per[ms] = dict(s=s, area=area)
    smax = max(smax, s.max())
gridA = np.linspace(0, smax, NPL)
allA = []
for ms in ms_sorted:
    s, area = per[ms]["s"], per[ms]["area"]
    o = np.argsort(s)
    allA.append(np.interp(gridA, s[o], area[o], left=np.nan, right=np.nan))
allA = np.array(allA)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Amin, Amax = np.nanmin(allA, axis=0), np.nanmax(allA, axis=0)
# the band is a min-max over frames: drop it where few frames reach (the path length varies a little)
Amin[np.isfinite(allA).mean(0) < 0.8] = np.nan; Amax[np.isfinite(allA).mean(0) < 0.8] = np.nan

# ---- WOB cumulative loss keyed by cfd_t_s ; ENERGY = dP*Q (mW) ----
_petimes, _pe = CP.load_pe(a.wob_csv, a.region)
wframes = {t: None for t in _petimes}
wper, wsmax = {}, 0.0
allE = []
for tw in sorted(wframes):
    d = _pe[tw]
    sw = d["s"]
    Q = d["Q"]; tpm = d["tp_m"]; tpa = d["tp_a"]
    tpp = np.where(np.abs(Q) < 2.0, tpa, np.where(np.isfinite(tpm), tpm, tpa))
    cum = tpp[0] - tpp                                        # distal-referenced
    q = abs(float(np.interp(tw, tp, qp))) / 60000.0           # bulk Q (m3/s), prescribed
    wper[tw] = dict(sw=sw, cum=cum, q=q)
    wsmax = max(wsmax, sw.max())
wtimes = np.array(sorted(wper))
wgridA = np.linspace(0, wsmax, NPL)
for tw in wtimes:
    sw, cum, q = wper[tw]["sw"], wper[tw]["cum"], wper[tw]["q"]
    o = np.argsort(sw)
    allE.append(np.interp(wgridA, sw[o], np.abs(cum * q * 1000.0)[o], left=np.nan, right=np.nan))
allE = np.array(allE)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Emin, Emax = np.nanmin(allE, axis=0), np.nanmax(allE, axis=0)
Emin[np.isfinite(allE).mean(0) < 0.8] = np.nan; Emax[np.isfinite(allE).mean(0) < 0.8] = np.nan

# ---- velocity colourbar range: use the real data range unless overridden ----
VMAX = a.vmax
if VMAX is None:
    VMAX = 1.0
    print("  (no --vmax given, colourbar 0-%.3g m/s)" % VMAX, flush=True)

if a.theme == "light":
    TH = dict(bg="white", fg="#15181d", tick="#4a4a4a", title="#12305c", spine="0.5",
              legfc="#eef2f7", legec="0.6", grid="0.78", flow="0.32", zero="0.55", cbtxt="#222")
else:
    TH = dict(bg="#0f0f14", fg="0.92", tick="0.82", title="0.95", spine="0.35",
              legfc="0.12", legec="0.35", grid=None, flow="0.6", zero="0.4", cbtxt="0.92")
FS = a.font_scale
plt.rcParams.update({
    "font.weight": "bold", "font.size": 16 * FS,
    "axes.titleweight": "bold", "axes.labelweight": "bold",
    "axes.titlesize": 17 * FS, "axes.labelsize": 16.5 * FS,
    "xtick.labelsize": 14 * FS, "ytick.labelsize": 14 * FS,
    "text.color": TH["fg"], "axes.labelcolor": TH["fg"],
    "xtick.color": TH["tick"], "ytick.color": TH["tick"],
})
BG = TH["bg"]
VCROP_FRAC = (0.02, 0.97, 0.02, 0.80)   # y0,y1,x0,x1 fractions: drop the built-in colorbar
REF = mpimg.imread(a.ref) if a.ref and os.path.exists(a.ref) else None
WREF = mpimg.imread(a.wob_ref) if a.wob_ref and os.path.exists(a.wob_ref) else None
VIR = plt.get_cmap("viridis")
DMARKS = (list(np.arange(a.mark_step, smax, a.mark_step)) if a.mark_step > 0 else
          [d for d in (25, 50, 75, 100, 150, 200) if d < smax])
C_IN, C_EX = "#2a78d6", "#eb6834"                     # inhale / exhale


def style_axes(ax):
    for sp in ax.spines.values():
        sp.set_color(TH["spine"])
    ax.grid(alpha=0.35, color=TH["grid"]) if TH["grid"] else ax.grid(alpha=0.15)


def _img(path):
    return mpimg.imread(path) if os.path.exists(path) else None


def content_box(paths, mask_corners=False, pad=0.03):
    """Union bounding box (fractions y0,y1,x0,x1) of non-background pixels over sample frames.
    mask_corners ignores the renderer's bottom-left text/triad and bottom-right colourbar."""
    box = None
    for pth in paths:
        im = _img(pth)
        if im is None: continue
        im = im[..., :3]; h, w = im.shape[:2]
        diff = np.abs(im - im[2, 2]).max(axis=2) > 0.06
        if mask_corners:
            diff[int(0.80 * h):, :int(0.22 * w)] = False; diff[int(0.55 * h):, int(0.80 * w):] = False
        ys, xs = np.nonzero(diff)
        if ys.size == 0: continue
        b = (ys.min() / h, ys.max() / h, xs.min() / w, xs.max() / w)
        box = b if box is None else (min(box[0], b[0]), max(box[1], b[1]), min(box[2], b[2]), max(box[3], b[3]))
    if box is None: return (0.0, 1.0, 0.0, 1.0)
    return (max(0.0, box[0] - pad), min(1.0, box[1] + pad), max(0.0, box[2] - pad), min(1.0, box[3] + pad))


def crop(im, box, blank_corners=False):
    h, w = im.shape[:2]
    if blank_corners:
        im = im.copy(); bgc = im[2, 2]
        im[int(0.80 * h):, :int(0.22 * w)] = bgc; im[int(0.55 * h):, int(0.80 * w):] = bgc
    return im[int(box[0] * h):int(np.ceil(box[1] * h)), int(box[2] * w):int(np.ceil(box[3] * w))]


def pad_right(im, frac):
    """Widen with background on the right so the colourbar sits beside the geometry, not on it."""
    return np.concatenate([im, np.broadcast_to(im[2, 2], (im.shape[0], int(frac * im.shape[1])) + im.shape[2:])], axis=1)


KLO, KHI = 0, NF
if a.frames:
    _a, _b = a.frames.split(":")
    KLO, KHI = (int(_a) if _a else 0), (int(_b) if _b else NF)

if a.layout in ("columns", "report"):
    # union over frames spread through the breath plus the requested ones (a snapshot run renders only those)
    _samp = sorted(set(range(0, NF, max(1, NF // 12))) | set(range(KLO, KHI)))
    BOX_B = content_box(["%s/frame_%04d.png" % (a.planes_dir, k) for k in _samp])
    BOX_V = content_box(["%s/frame_%04d.png" % (a.vel_dir, k) for k in _samp], mask_corners=True)
    if REF is not None:
        REF = crop(REF, content_box([a.ref]))
    if WREF is not None:
        WREF = crop(WREF, content_box([a.wob_ref]))
    if a.wob_ref_dir:
        BOX_W = content_box(["%s/frame_%04d.png" % (a.wob_ref_dir, k) for k in _samp])
    print("crop boxes (y0,y1,x0,x1): planes %s  velocity %s" % (np.round(BOX_B, 3), np.round(BOX_V, 3)), flush=True)


def compose(k):
    ms = ms_sorted[k]
    t = ms / 1000.0
    Bi = _img("%s/frame_%04d.png" % (a.planes_dir, k))
    Vimg = _img("%s/frame_%04d.png" % (a.vel_dir, k))
    qnow = float(np.interp(t, tp, qp))
    # Sign convention verified against the exported caps (postproc/outlet_flow/
    # inlet_outlet_massflow.csv): a POSITIVE prescribed mass flow is air entering the
    # domain at BronchiOut and leaving at MaskInlet, i.e. EXHALATION. ENT003 exhales
    # first and inhales after t~1.36 s.
    phase = "INHALE" if qnow < -0.3 else ("EXHALE" if qnow > 0.3 else "transition")

    if a.layout == "columns":
        W, H = (float(x) for x in (a.figsize or "19.2,10.8").split(","))
        fig = plt.figure(figsize=(W, H), facecolor=BG)
        gs = fig.add_gridspec(2, 3, height_ratios=[1.9, 1], hspace=0.17, wspace=0.20,
                              left=0.062, right=0.985, top=0.955, bottom=0.085)
        axB, axL, axV = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
        axf, axC, axE = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])
        if Bi is not None: axB.imshow(crop(Bi, BOX_B))
        if REF is not None: axL.imshow(REF)
        if Vimg is not None: axV.imshow(pad_right(crop(Vimg, BOX_V, blank_corners=True), 0.22))
        cax = axV.inset_axes([0.845, 0.10, 0.035, 0.72])
        cb_orient, cb_lab_pos = "vertical", None
    elif a.layout == "report":
        W, H = (float(x) for x in (a.figsize or "8.5,11").split(","))
        fig = plt.figure(figsize=(W, H), facecolor=BG)
        # title, then flow over time (where in the breath this page is), the renders at full page width,
        # then WOB per segment, CSA and WOB rate on one shared distance axis (same marks as the locator labels)
        if a.page_title:
            fig.suptitle(u"%s  —  t = %.2f s" % (a.page_title, t), x=0.5, y=0.992, va="top",
                         fontsize=27 * FS, fontweight="bold", color=TH["title"])
        g0 = fig.add_gridspec(1, 1, left=0.085, right=0.985, top=0.935, bottom=0.875)
        g1 = fig.add_gridspec(1, 3, width_ratios=[0.92, 0.98, 1.1], wspace=0.02,
                              left=0.01, right=0.99, top=0.815, bottom=0.475)
        gb = fig.add_gridspec(3, 1, height_ratios=[1.25, 0.95, 0.95], hspace=0.55,
                              left=0.085, right=0.985, top=0.445, bottom=0.045)
        axf = fig.add_subplot(g0[0])
        axB, axW, axV = (fig.add_subplot(g1[0, i]) for i in range(3))
        axS, axC, axE = (fig.add_subplot(gb[i]) for i in range(3))
        axL = axB                                     # first panel = the arc-length locator
        if REF is not None: axB.imshow(REF)
        Wimg = _img("%s/frame_%04d.png" % (a.wob_ref_dir, k)) if a.wob_ref_dir else WREF
        if Wimg is not None: axW.imshow(crop(Wimg, BOX_W) if a.wob_ref_dir else Wimg)
        axW.axis("off"); axW.set_title("WOB at t = %.2f s" % t, fontsize=17 * FS * 1.12, color=TH["title"])
        if Vimg is not None: axV.imshow(pad_right(crop(Vimg, BOX_V, blank_corners=True), 0.30))
        cax = axV.inset_axes([0.80, 0.12, 0.05, 0.70])
        cb_orient, cb_lab_pos = "vertical", None
    else:
        W, H = (float(x) for x in (a.figsize or "12.0,13.0").split(","))
        fig = plt.figure(figsize=(W, H), facecolor=BG)
        gs = fig.add_gridspec(3, 12, height_ratios=[1, 1, 1],
                              hspace=0.34, wspace=0.10, left=0.055, right=0.975, top=0.955, bottom=0.05)
        axB, axL, axV = fig.add_subplot(gs[0, 0:4]), fig.add_subplot(gs[1, 0:4]), fig.add_subplot(gs[2, 0:4])
        axf, axC, axE = fig.add_subplot(gs[0, 5:12]), fig.add_subplot(gs[1, 5:12]), fig.add_subplot(gs[2, 5:12])
        if Bi is not None: axB.imshow(Bi)
        if REF is not None: axL.imshow(REF)
        if Vimg is not None:
            # drop the renderer's built-in colourbar/annotations so only the drawn bar is shown
            vh, vw = Vimg.shape[:2]
            fy0, fy1, fx0, fx1 = VCROP_FRAC
            axV.imshow(Vimg[int(fy0 * vh):int(fy1 * vh), int(fx0 * vw):int(fx1 * vw)])
        # horizontal colourbar in the empty bottom strip of the render
        cax = axV.inset_axes([0.13, 0.055, 0.60, 0.030])
        cb_orient, cb_lab_pos = "horizontal", "top"
    TTL = {"columns": ("airway wall + CSA planes", "plane locator", "velocity"),
           "report": ("distance along centerline", None, "velocity"),
           "stacked": ("airway wall + CSA planes (moving)", "plane locator", "velocity flow (m/s)")}[a.layout]
    for ax_, ttl in zip((axB, axL, axV), TTL):
        if ttl is None:
            continue
        ax_.axis("off"); ax_.set_title(ttl, fontsize=16 * FS * (1.12 if a.layout != "stacked" else 1.0), color=TH["title"])
    smap = plt.cm.ScalarMappable(cmap="turbo", norm=plt.Normalize(0, VMAX))
    cbv = fig.colorbar(smap, cax=cax, orientation=cb_orient)
    cbv.set_label("velocity (m/s)", color=TH["title"], fontsize=14 * FS, fontweight="bold")
    if cb_lab_pos:
        cbv.ax.xaxis.set_label_position(cb_lab_pos)
    cbv.ax.tick_params(labelsize=12.5 * FS, color=TH["tick"])
    plt.setp(plt.getp(cbv.ax, "xticklabels" if cb_orient == "horizontal" else "yticklabels"), color=TH["cbtxt"], fontweight="bold")
    cbv.outline.set_edgecolor(TH["spine"])

    # ---- RIGHT column: line plots ----
    axf.set_facecolor(BG)
    axf.plot(tp, qp, color=TH["flow"], lw=1.9); axf.axhline(0, color=TH["zero"], lw=0.6)
    axf.axvline(t, color="#e74c3c", lw=1.0, alpha=0.5)
    axf.plot(t, qnow, "o", color="#e74c3c", ms=12, zorder=6)
    axf.set_xlim(tp.min(), tp.max()); axf.set_ylabel("flow (L/min)"); axf.set_xlabel("time (s)")
    # columns/report: legends sit above the axes (never on the data), titles above the legends
    # report: full-width plots, title on the left and legend on the right of the same line
    WIDE, REP = a.layout != "stacked", a.layout == "report"
    LEG = (dict(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=13 * FS,
                handlelength=1.4, columnspacing=1.0, borderaxespad=0.2, labelcolor=TH["fg"]) if WIDE else
           dict(fontsize=14 * FS, facecolor=TH["legfc"], edgecolor=TH["legec"], labelcolor=TH["fg"], framealpha=0.92))
    if REP:
        LEG.update(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3)
    TPAD = 5.0 if REP else (13 * FS * 2.0 if WIDE else 6.0)
    TLOC = "left" if REP else "center"
    L_RANGE, L_NOW = {"report": ("min–max over breath", "t = %.2f s" % t), "columns": ("min–max over breath", "now"),
                      "stacked": ("min–max over breath", None)}[a.layout]
    if REP:
        tr = axf.get_xaxis_transform()
        axf.fill_between(tp, 0, 1, where=qp < 0, color=C_IN, alpha=0.13, lw=0, transform=tr, label="inhale")
        axf.fill_between(tp, 0, 1, where=qp > 0, color=C_EX, alpha=0.13, lw=0, transform=tr, label="exhale")
        axf.legend(**LEG)
    axf.set_title({"report": "flow    t = %.2f s    %s", "columns": "t = %.2f s    %s",
                   "stacked": "breathing cycle    t = %.2f s    —    %s"}[a.layout] % (t, phase), pad=TPAD, loc=TLOC,
                  color=TH["title"], fontsize=17 * FS)
    style_axes(axf)

    axC.set_facecolor(BG)
    for d in DMARKS:
        axC.axvline(d, color=VIR(d / max(smax, 1e-9)), lw=1.3, ls=(0, (2, 2)), alpha=0.55, zorder=0)
    axC.fill_between(gridA, Amin, Amax, color="#3498db", alpha=0.20, label=L_RANGE)
    axC.plot(per[ms]["s"], per[ms]["area"], color="#f39c12", lw=2.6, label=L_NOW or "CSA now")
    XL = CP.xlabel(a.region).split("   ")[0].split("  —")[0] if WIDE else CP.xlabel(a.region)
    XL = XL.replace(CP.PATH_LABEL[a.region][1], a.prox_label or CP.PATH_LABEL[a.region][1])
    XMAX = min(smax, wsmax)
    if REP:                                          # shared distance axis: segment labels above, mm below
        axC.tick_params(labelbottom=False)
        axE.set_xticks([0] + (DMARKS[1::2] if len(DMARKS) > 8 else DMARKS))
    else:
        axC.set_xlabel(XL)
    axC.set_ylabel("CSA (mm²)")
    axC.legend(**(LEG if WIDE else dict(LEG, loc="upper right")))
    axC.set_xlim(0, XMAX if REP else smax)
    if a.csa_ylim:
        axC.set_ylim(*[float(x) for x in a.csa_ylim.split(",")])
    style_axes(axC)
    axC.set_title("cross-sectional area" if WIDE else "cross-sectional area along the airway", pad=TPAD, loc=TLOC,
                  color=TH["title"], fontsize=16.5 * FS)

    tw = wtimes[int(np.argmin(np.abs(wtimes - t)))]
    e_now = np.abs(wper[tw]["cum"] * wper[tw]["q"]) * 1000.0
    axE.set_facecolor(BG)
    for d in DMARKS:
        axE.axvline(d, color=VIR(d / max(smax, 1e-9)), lw=1.3, ls=(0, (2, 2)), alpha=0.55, zorder=0)
    axE.fill_between(wgridA, Emin, Emax, color="#16a085", alpha=0.20, label=L_RANGE)
    axE.plot(wper[tw]["sw"], e_now, color="#f1c40f", lw=2.6, label=L_NOW or "energy loss now")
    WOBN = WIDE                                      # name it work of breathing (WOB); stacked keeps the ENT003 wording
    axE.axhline(0, color=TH["zero"], lw=0.7)
    axE.set_xlabel(XL)
    axE.set_ylabel("WOB rate (mW)" if WOBN else "energy loss (mW)")
    axE.legend(**(LEG if WIDE else dict(LEG, loc="upper left")))
    axE.set_xlim(0, XMAX if REP else wsmax)
    if a.energy_ylim:
        axE.set_ylim(*[float(x) for x in a.energy_ylim.split(",")])
    style_axes(axE)
    # "... along the airway" overflowed the figure and clipped; the x-axis already says it
    axE.set_title({"report": "WOB rate = ΔP·Q, accumulated from the nostril",
                   "columns": "WOB rate = ΔP·Q", "stacked": "energy loss = ΔP·Q   (work of breathing)"}[a.layout],
                  pad=TPAD, loc=TLOC,
                  color=TH["title"], fontsize=16.5 * FS)

    if a.layout == "report" and a.bins_csv:
        # work of breathing per centerline segment, inhale vs exhale (wob_regions.py --bin-mm)
        SB = list(csv.DictReader(open(a.bins_csv)))
        lo, hi = (np.array([float(r[c]) for r in SB]) for c in ("lo_mm", "hi_mm"))
        xc, bw = 0.5 * (lo + hi), 0.38 * (hi - lo)
        axS.set_facecolor(BG)
        for d in DMARKS:
            axS.axvline(d, color=VIR(d / max(smax, 1e-9)), lw=1.3, ls=(0, (2, 2)), alpha=0.55, zorder=0)
        Ein, Eex = (np.array([float(r[c]) for r in SB]) for c in ("inhale_mJ", "exhale_mJ"))
        if a.bins_power_csv:                         # video: accumulated up to t, whole breath as an outline
            PT = list(csv.DictReader(open(a.bins_power_csv)))
            pc = [c for c in PT[0] if c.endswith("_mW")]
            done = [r for r in PT if float(r["time_s"]) <= t + 1e-6]
            acc = {ph: np.array([sum(float(r[c]) * float(r["dt_s"]) for r in done if r["phase"] == ph) for c in pc])
                   for ph in ("inhale", "exhale")}
            axS.bar(xc - bw / 2 - 0.4, Ein, bw, facecolor="none", edgecolor=C_IN, lw=1.0, label="whole breath")
            axS.bar(xc + bw / 2 + 0.4, Eex, bw, facecolor="none", edgecolor=C_EX, lw=1.0)
            Ein, Eex = acc["inhale"], acc["exhale"]
        axS.bar(xc - bw / 2 - 0.4, Ein, bw, color=C_IN, label="inhale")
        axS.bar(xc + bw / 2 + 0.4, Eex, bw, color=C_EX, label="exhale")
        step = np.median(hi - lo)                    # a short last segment gets a two-line label so it fits
        axS.set_xticks(xc); axS.set_xticklabels([("%s–%s" if h - l > 0.75 * step else "%s–\n%s") % (r["lo_mm"], r["hi_mm"])
                                                 for r, l, h in zip(SB, lo, hi)])
        axS.tick_params(axis="x", length=0, labelsize=12.5 * FS)
        axS.set_xlim(0, XMAX)
        axS.set_ylabel("WOB (mJ)")
        axS.grid(axis="y", alpha=0.35, color=TH["grid"]); axS.set_axisbelow(True)
        for sp_ in ("top", "right"): axS.spines[sp_].set_visible(False)
        for sp_ in ("left", "bottom"): axS.spines[sp_].set_color(TH["spine"])
        axS.legend(**LEG)
        axS.set_title(("WOB per %d mm, accumulated to t = %.2f s" % (round(hi[0] - lo[0]), t)) if a.bins_power_csv
                      else ("work of breathing (WOB) per %d mm of centerline, whole breath" % round(hi[0] - lo[0])),
                      loc="left", pad=TPAD, color=TH["title"], fontsize=16.5 * FS)

    fig.savefig("%s/frame_%04d.png" % (a.out, k), dpi=a.dpi, facecolor=BG)
    if a.pdf_out:
        fig.savefig(a.pdf_out, dpi=a.dpi, facecolor=BG)
    plt.close(fig)


for k in range(KLO, KHI):
    compose(k)
    if k % 10 == 0:
        print("  composited %d/%d" % (k, NF), flush=True)
print("DONE -> %s (%d frames)" % (a.out, KHI - KLO), flush=True)
