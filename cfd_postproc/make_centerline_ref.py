#!/usr/bin/env python3
"""Static anatomical locator for the breath video / clinical PDF.

Airway silhouette (sagittal Y-Z projection, matching the 'sagittal' camera used by
render_velocity_volume.py and render_stl_planes_seq.py) with the CSA centerline drawn
on top, coloured and labelled by distance from the airway opening.

Ported from /storage4/ENT002/make_centerline_ref.py. Differences:
  * --subject-dir / --csa / --out are arguments instead of module constants
  * ENT002 had one centerline; ENT003 has four regions, so --path selects a route
    (see csa_paths.py) and --path all draws every region in its own colour, matching
    the region colours used by compose_csa_video.py
  * distance marks are chosen from the actual path length rather than ENT002's
    fixed 0/50/100/150/200 mm

Usage:
  python3 make_centerline_ref.py --subject-dir ENT003 \
      --csa ENT003/csa_bifurcation/ENT003-Data.csv --frame out_001950 \
      --path LeftNose --out postproc/panels
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import csa_paths as CP

RCOL = {"DescendingAirway": "#4aa3df", "Mouth": "#f39c12",
        "LeftNose": "#2ecc71", "RightNose": "#e74c3c"}

THEMES = {
    "dark":  dict(bg="#0f0f14", scat="#8a94a6", sa=0.035, txt="white", tag="#e2e6ee",
                  box="#1b1b24", title="0.95", cbtxt="0.9", cbtick="0.7", mec="white"),
    "light": dict(bg="white", scat="#33404f", sa=0.055, txt="#15181d", tag="#33404f",
                  box="#eef2f7", title="#12305c", cbtxt="#222", cbtick="#555", mec="#222"),
}

p = argparse.ArgumentParser()
p.add_argument("--subject-dir", required=True)
p.add_argument("--csa", required=True)
p.add_argument("--frame", default="out_001950")
p.add_argument("--path", default="LeftNose", help="a csa_paths path name, or 'all'")
p.add_argument("--out", required=True)
p.add_argument("--font-scale", type=float, default=1.0)
p.add_argument("--no-title", action="store_true", help="omit the inner title (the composer titles the panel)")
p.add_argument("--suffix", default="", help="appended to the output file name")
p.add_argument("--themes", default="dark,light")
p.add_argument("--no-colorbar", action="store_true", help="omit the colourbar (distance labels say the same)")
p.add_argument("--mark-step", type=float, default=0, help="distance label every N mm (default: 25 or 50 by path length)")
p.add_argument("--prox-label", default=None, help="name for the opening instead of e.g. 'left nostril'")
p.add_argument("--color-csv", default=None,
               help="colour the centerline by a quantity instead of distance: CSV with columns s_mm,value "
                    "(e.g. wob_regions.py --peak-time output); no distance labels are drawn")
p.add_argument("--color-label", default="local energy loss (mW/mm)")
p.add_argument("--tag", default="ref", help="output name: centerline_<tag>_<path>...")
p.add_argument("--color-max", type=float, default=None, help="fixed top of the colour scale (same across video frames)")
p.add_argument("--out-file", default=None, help="explicit output PNG (one theme only)")
p.add_argument("--dpi", type=int, default=300)
p.add_argument("--limits-from", default=None, help="take the view limits from this frame's STL (steady framing in a video)")
p.add_argument("--no-tight", action="store_true", help="keep the full canvas (same size every frame)")
p.add_argument("--color-series", default=None,
               help="video: dir of t<time>.csv files (wob_regions.py --series); writes <out>/frame_<k>.png for every motion "
                    "frame k, coloured by the nearest time step (pass --color-max, --limits-from, --no-tight)")
p.add_argument("--jobs", type=int, default=12)
a = p.parse_args()

if a.color_series:
    # one process per motion frame (this script handles one frame), run in parallel
    import glob, subprocess, sys
    from concurrent.futures import ThreadPoolExecutor
    T = sorted(float(os.path.basename(f)[1:-4]) for f in glob.glob(os.path.join(a.color_series, "t*.csv")))
    argv, base = sys.argv[1:], []
    i = 0
    while i < len(argv):                                  # drop the per-frame options, keep the rest
        if argv[i] in ("--color-series", "--frame", "--jobs"):
            i += 2; continue
        base.append(argv[i]); i += 1
    os.makedirs(a.out, exist_ok=True)
    cmds = []
    for f in sorted(glob.glob(os.path.join(a.subject_dir, "motion", "stl", "out_*.stl"))):
        ms = int(os.path.basename(f)[4:10]); t = min(T, key=lambda x: abs(x - ms / 1000.0))
        cmds.append([sys.executable, os.path.abspath(__file__)] + base + [
            "--frame", "out_%06d" % ms, "--color-csv", os.path.join(a.color_series, "t%.3f.csv" % t),
            "--out-file", os.path.join(a.out, "frame_%04d.png" % (ms // 25))])
    with ThreadPoolExecutor(a.jobs) as ex:
        bad = [c for c, r in zip(cmds, ex.map(lambda c: subprocess.run(c, capture_output=True), cmds)) if r.returncode]
    print("color series: %d frames -> %s%s" % (len(cmds), a.out, ("  (%d FAILED)" % len(bad)) if bad else ""))
    sys.exit(1 if bad else 0)

os.makedirs(a.out, exist_ok=True)
V = CP.stl_vertices(os.path.join(a.subject_dir, "motion", "stl", "%s.stl" % a.frame))
vy, vz = V[:, 1], V[:, 2]
VL = CP.stl_vertices(os.path.join(a.subject_dir, "motion", "stl", "%s.stl" % a.limits_from)) if a.limits_from else V

import csv as _csv
rows = [r for r in _csv.DictReader(open(a.csa))
        if r["frame_name"] == a.frame and r["is_valid"].lower() == "true"]
if not rows:
    raise SystemExit("no valid CSA rows for frame %s" % a.frame)


def marks_for(smax):
    step = int(a.mark_step) or (25 if smax <= 120 else 50)
    m = list(range(0, int(smax) + 1, step))
    if m[-1] < smax - step * 0.35 and not a.mark_step:      # fixed step: the end is tagged, not numbered
        m.append(int(smax))
    return m


def place_labels(fig, ax, items, avoid_xy, blocked_axes=(), gap_pt=14):
    """Put each label next to its point but OFF the airway: candidate positions around the point
    at growing distance, rejected if the box touches the projected surface (avoid_xy, data
    coords, grown by gap_pt so the text stands clear of it), another label or a blocked axes
    (colourbar). items: (x, y, text, kwargs) -> annotate() with an arrow to the point."""
    fig.canvas.draw()
    T = ax.transData.transform
    W, H = fig.canvas.get_width_height(); c = 4                          # occupancy grid, 4 px cells
    occ = np.zeros((H // c + 2, W // c + 2), bool)
    P = T(avoid_xy); ok = (P[:, 0] >= 0) & (P[:, 0] < W) & (P[:, 1] >= 0) & (P[:, 1] < H)
    occ[(P[ok, 1] // c).astype(int), (P[ok, 0] // c).astype(int)] = True
    for _ in range(2 + int(gap_pt * fig.dpi / 72.0 / c)):              # close gaps + keep a clear margin
        occ = occ | np.roll(occ, 1, 0) | np.roll(occ, -1, 0) | np.roll(occ, 1, 1) | np.roll(occ, -1, 1)
    boxes = [bx.get_window_extent().extents for bx in blocked_axes]
    px2pt = 72.0 / fig.dpi
    for x, y, text, kw in items:
        probe = ax.text(0, 0, text, **{k: v for k, v in kw.items() if k not in ("arrowprops",)})
        m = (1.0 if "bbox" in kw else 0.4) * kw.get("fontsize", 12) * fig.dpi / 72.0   # box padding + gap
        bb = probe.get_window_extent(); w, h = bb.width + m, bb.height + m; probe.remove()
        x0, y0 = T((x, y))
        best = None
        for r in np.arange(24, 400, 4) / px2pt:                         # distance from the point, px
            for ang, pen in ((180, 0), (0, 2), (160, 4), (200, 4), (20, 4), (340, 4), (135, 8), (225, 8),
                             (45, 8), (315, 8), (90, 12), (270, 12)):
                ca, sa = np.cos(np.radians(ang)), np.sin(np.radians(ang))
                cx = x0 + ca * r + (np.sign(ca) * w / 2 if abs(ca) > 0.3 else 0)
                cy = y0 + sa * r + (np.sign(sa) * h / 2 if abs(sa) > 0.3 else 0)
                e = (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
                if e[0] < 0 or e[1] < 0 or e[2] > W or e[3] > H:
                    continue
                if occ[int(e[1] // c):int(e[3] // c) + 1, int(e[0] // c):int(e[2] // c) + 1].any():
                    continue
                if any(e[0] < b[2] and e[2] > b[0] and e[1] < b[3] and e[3] > b[1] for b in boxes):
                    continue
                cost = r * px2pt + pen
                if best is None or cost < best[0]:
                    best = (cost, cx, cy, e)
            if best is not None and r * px2pt > best[0]:
                break
        if best is None:
            continue
        _, cx, cy, e = best
        boxes.append(e)
        ax.annotate(text, (x, y), xytext=((cx - x0) * px2pt, (cy - y0) * px2pt), textcoords="offset points",
                    ha="center", va="center", **kw)


def render(theme):
    C = THEMES[theme]
    fig = plt.figure(figsize=(5.2, 8.8), facecolor=C["bg"])
    # untrimmed canvas (--no-tight): leave room for the colourbar label inside the figure
    ax = fig.add_axes([0.02, 0.03, 0.74 if a.no_tight else 0.82, 0.92]); ax.set_facecolor(C["bg"])
    ax.scatter(vy, vz, s=1.1, c=C["scat"], alpha=C["sa"], edgecolors="none", rasterized=True)

    padY = (VL[:, 1].max() - VL[:, 1].min()) * 0.10 + 8
    xl_left, xl_right = VL[:, 1].max() + padY, VL[:, 1].min() - padY   # +Y to the LEFT, as the renders show it

    if a.path == "all":
        drawn = set()
        for side in ("LeftNose", "RightNose"):                # real centerlines, split at the junction
            rs = [r for r in rows if r["region"] in (side, "DescendingAirway")]
            if not any(r["region"] == side for r in rs):
                continue
            arc = np.array([float(r["arc_length_mm"]) for r in rs]); sP = arc.max() - arc
            junc = max(x for x, r in zip(sP, rs) if r["region"] == side)
            CL = CP.side_centerline(a.subject_dir, side, a.frame, s_max=sP.max())
            if CL is None:
                continue
            P, sc = CL
            k = sc <= junc
            ax.plot(P[k, 1], P[k, 2], "-", color=RCOL[side], lw=5.0, alpha=0.9, label=side, solid_capstyle="round")
            drawn.add(side)
            if "DescendingAirway" not in drawn:
                k = sc >= junc
                ax.plot(P[k, 1], P[k, 2], "-", color=RCOL["DescendingAirway"], lw=5.0, alpha=0.9,
                        label="DescendingAirway", solid_capstyle="round")
                drawn.add("DescendingAirway")
        for reg in ("DescendingAirway", "LeftNose", "RightNose", "Mouth"):
            if reg in drawn:
                continue
            rs = [r for r in rows if r["region"] == reg]
            if not rs:
                continue
            arc = np.array([float(r["arc_length_mm"]) for r in rs])
            o = np.argsort(arc)
            cy = np.array([float(r["centroid_y"]) for r in rs])[o]
            cz = np.array([float(r["centroid_z"]) for r in rs])[o]
            ax.plot(cy, cz, "-", color=RCOL[reg], lw=5.0, alpha=0.9, label=reg, solid_capstyle="round")
        ax.legend(fontsize=10, loc="lower left", facecolor=C["box"],
                  edgecolor="0.4", labelcolor=C["txt"])
        title = "airway locator\nCSA regions"
        out = os.path.join(a.out, "centerline_ref_all" + ("" if theme == "dark" else "_light") + ".png")
        cb = None
    else:
        regs = set(CP.regions_of(a.path))
        rs = [r for r in rows if r["region"] in regs]
        arc = np.array([float(r["arc_length_mm"]) for r in rs])
        cy = np.array([float(r["centroid_y"]) for r in rs])
        cz = np.array([float(r["centroid_z"]) for r in rs])
        s = arc.max() - arc                                  # 0 at the opening
        o = np.argsort(s)
        cyO, czO, sO = cy[o], cz[o], s[o]
        CL = CP.side_centerline(a.subject_dir, a.path, a.frame, s_max=s.max())
        if CL is not None:                                   # draw the real (smoothed) centerline
            cyO, czO, sO = CL[0][:, 1], CL[0][:, 2], CL[1]
            cy, cz, s = cyO, czO, sO

        pts = np.column_stack([cyO, czO]).reshape(-1, 1, 2)
        seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
        if a.color_csv:                                      # colour by a quantity along the centerline
            Cq = np.loadtxt(a.color_csv, delimiter=",", skiprows=1)
            val = np.clip(np.interp(0.5 * (sO[:-1] + sO[1:]), Cq[:, 0], Cq[:, 1]), 0, None)
            cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
                "or", plt.get_cmap("OrRd")(np.linspace(0.12, 1.0, 256)))
            norm = plt.Normalize(0.0, a.color_max or float(np.percentile(val, 99.5)))
            ax.plot(cyO, czO, "-", color="#3a2a22", lw=8.0, solid_capstyle="round", zorder=3)   # outline: low values stay visible
            lc = LineCollection(seg, cmap=cmap, norm=norm, array=val, linewidths=6.0, zorder=4)
            ax.add_collection(lc)
            cb_label, labels = a.color_label, []
        else:
            cmap = plt.get_cmap("viridis"); norm = plt.Normalize(0.0, s.max())
            lc = LineCollection(seg, cmap=cmap, norm=norm, array=sO[:-1], linewidths=5.5)
            ax.add_collection(lc)
            cb_label, labels = CP.PATH_LABEL[a.path][0].replace(CP.PATH_LABEL[a.path][1], a.prox_label or CP.PATH_LABEL[a.path][1]), []
            for d in marks_for(s.max()):
                i = int(np.argmin(np.abs(s - d)))
                col = cmap(d / s.max())
                ax.plot(cy[i], cz[i], "o", ms=10 * a.font_scale ** 0.5, mfc=col, mec=C["mec"], mew=1.6, zorder=6)
                labels.append((cy[i], cz[i], "%d mm" % d, dict(
                    fontsize=13 * a.font_scale, color=C["txt"], fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", fc=C["box"], ec=col, lw=1.6),
                    arrowprops=dict(arrowstyle="-|>", mutation_scale=11 * a.font_scale, color=col, lw=1.6,
                                    shrinkA=0, shrinkB=6))))
        _, prox, dist = CP.PATH_LABEL[a.path]
        prox = a.prox_label or prox
        for i_, txt in ((int(np.argmin(s)), prox), (int(np.argmax(s)), dist)):
            labels.append((cy[i_], cz[i_], txt, dict(fontsize=12 * a.font_scale, color=C["tag"], fontweight="bold",
                                                      style="italic", arrowprops=dict(arrowstyle="-|>", color=C["tag"],
                                                                                      mutation_scale=9 * a.font_scale,
                                                                                      lw=1.1, shrinkA=0, shrinkB=4))))
        title = "airway locator\ncenterline distance from %s" % prox
        out = os.path.join(a.out, "centerline_%s_%s%s%s.png"
                           % (a.tag, a.path, "" if theme == "dark" else "_light", a.suffix))
        cax = fig.add_axes([0.79 if a.no_tight else 0.87, 0.14, 0.04, 0.68])
        if a.no_colorbar:
            cax.set_visible(False)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label(cb_label, color=C["title"], fontsize=13 * a.font_scale, fontweight="bold")
        cb.ax.tick_params(labelsize=12 * a.font_scale, color=C["cbtick"])
        plt.setp(plt.getp(cb.ax, "yticklabels"), color=C["cbtxt"], fontweight="bold")
        cb.outline.set_edgecolor("0.4" if theme == "dark" else "0.6")

    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(xl_left, xl_right); ax.set_ylim(VL[:, 2].min() - 6, VL[:, 2].max() + 6)
    if a.path != "all":                                      # labels go where they do not cover the airway
        place_labels(fig, ax, labels, np.column_stack([np.r_[vy, cyO], np.r_[vz, czO]]),
                     blocked_axes=() if a.no_colorbar else (cax,))
    if not a.no_title:
        ax.set_title(title, color=C["title"], fontsize=15 * a.font_scale, fontweight="bold", pad=8)
    out = a.out_file or out
    fig.savefig(out, dpi=a.dpi, facecolor=C["bg"], **({} if a.no_tight else dict(bbox_inches="tight", pad_inches=0.12)))
    plt.close(fig)
    print("wrote", out)


for th in a.themes.split(","):
    render(th)
