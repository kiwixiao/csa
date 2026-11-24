#!/usr/bin/env python3
"""One-page clinical PDF for one subject (ENT003 original, generalised 2026-10-21).

A single breath-video snapshot at the airway bottleneck (peak dissipated power) plus a
summary of what the CFD adds for this patient. Every number is recomputed here from the
post-processing CSVs so the page can never drift from the figures.

Ported from /storage4/ENT002/build_onepager.py. Differences:
  * all paths are arguments; the subject label and the airway path are parameters
  * ENT002 measured distance from the mouth along one region; ENT003 is a nasal-mask
    case, so the route is nostril -> trachea (see csa_paths.py) and the wording follows
  * the snapshot time is taken from the peak of the recomputed dissipated power rather
    than being hardcoded, and the inspiratory fraction comes from the prescribed profile
    (ENT003 exhales first, inhales after t~1.36 s)
  * the finding paragraph is assembled from the data, so re-running after a new export
    updates the prose as well as the figures

Usage:
  python3 build_onepager.py --subject ENT003 --path LeftNose \
      --wob postproc/csa_pe/plots_wob/LeftNose/wob_summary.csv \
      --anatomy postproc/csa_anatomy/LeftNose/csa_anatomy_summary.csv \
      --csa ENT003/csa_bifurcation/ENT003-Data.csv \
      --snapshot postproc/panels/snapshot_light_hires.png \
      --out postproc/ENT003_CFD_onepager.pdf
"""
import argparse
import csv
import os
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch

import csa_paths as CP

p = argparse.ArgumentParser()
p.add_argument("--subject", default="ENT003")
p.add_argument("--path", default="LeftNose")
p.add_argument("--wob", required=True)
p.add_argument("--anatomy", required=True)
p.add_argument("--pe-dir", required=True)
p.add_argument("--prof", required=True)
p.add_argument("--csa", required=True)
p.add_argument("--snapshot", required=True)
p.add_argument("--out", required=True)
p.add_argument("--peak-velocity", type=float, default=None,
               help="peak |v| in the domain (m/s) for the narrative; omit to leave it out")
p.add_argument("--regions-csv", default=None, help="wob_regions.py CSV -> adds page 2 (where/when)")
p.add_argument("--regions-png", default=None)
p.add_argument("--smooth-mm", type=float, default=5.0,
               help="window for smoothing the along-airway loss gradient")
a = p.parse_args()

NAVY, ACC, INK = "#12305c", "#c0392b", "#1a1a1a"
_, PROX, DIST = CP.PATH_LABEL[a.path]

# ---- breath-level metrics ---------------------------------------------------------
W = list(csv.DictReader(open(a.wob)))
T = np.array([float(r["time_s"]) for r in W])
P = np.array([float(r["dissipated_power_mW"]) for r in W])
C = np.array([float(r["cum_WOB_mJ"]) for r in W])
DP = np.array([abs(float(r["dP_opening_trachea_Pa"])) for r in W])
INH = np.array([r["is_inhale"] == "1" for r in W])
dt = np.diff(T, prepend=T[0])

wob_tot = float(C[-1])
peak_pow = float(P.max())
tpk = float(T[int(np.argmax(P))])
insp_frac = 100.0 * float(np.sum((P * dt)[INH])) / float(np.sum(P * dt))
dP_pk = float(DP[int(np.argmax(P))])
dP_max = float(DP.max())
# phase wording from the data: STAR convention, prescribed mass flow > 0 = exhalation
PHASE = "inspiration" if INH[int(np.argmax(P))] else "expiration"
prof = np.array([[float(x) for x in r[:2]] for r in csv.reader(open(a.prof)) if r and r[0].replace(".", "", 1).replace("-", "", 1).isdigit()])
first = prof[np.flatnonzero(np.abs(prof[:, 1]) > 1e-12)[0], 1]
ORDER = "one inhalation plus one exhalation" if first < 0 else "one exhalation plus one inhalation"

# ---- locate the bottleneck by DISSIPATION, not by minimum area --------------------
# The narrowest cross-section is not necessarily where the pressure is lost: for ENT003
# the nostril plane is as tight as the glottis but drops ~0 Pa, while the laryngeal
# segment carries most of the loss. Finding the bottleneck from the along-airway loss
# gradient at peak inspiration is the whole point of running the CFD.
petimes, pe = CP.load_pe(a.pe_dir, a.path)
tsnap = min(petimes, key=lambda t: abs(t - tpk))
d = pe[tsnap]
o = np.argsort(d["s"])
S = d["s"][o]
CUM = np.abs(CP.cum_from_opening(d["tp_a"])[o])
ds = float(np.median(np.diff(S)))
win = max(3, int(round(a.smooth_mm / max(ds, 1e-6))))
grad = np.convolve(np.gradient(CUM, S), np.ones(win) / win, mode="same")

i_star = int(np.nanargmax(grad))
s_star = float(S[i_star])
# Contiguous span where the loss gradient stays above a fraction of its peak.
# 0.5 (half-max) is too brittle: on the smoothed CSA the gradient peak sharpens and the
# span collapsed to 4 mm, reporting 37% of the drop where the anatomical segment view
# says 59%. 0.25 is stable (0.15 returns the identical span) and agrees with it.
half = grad[i_star] * 0.25
lo = i_star
while lo > 0 and grad[lo - 1] > half:
    lo -= 1
hi = i_star
while hi < len(grad) - 1 and grad[hi + 1] > half:
    hi += 1
seg_lo, seg_hi = float(S[lo]), float(S[hi])
seg_dP = float(CUM[hi] - CUM[lo])
seg_frac = 100.0 * seg_dP / max(CUM.max(), 1e-9)

# plane index nearest the bottleneck, so it can be found again in the video
atimes, aper = CP.load_anatomy(a.csa, a.path)
tref = min(atimes, key=lambda t: abs(t - tsnap))
j = int(np.argmin(np.abs(aper[tref]["s"] - s_star)))
constr_plane = int(aper[tref]["pidx"][j])
constr_reg = str(aper[tref]["region"][j])

# CSA at the bottleneck through the breath, tracked by PLANE INDEX. Interpolating at a
# fixed distance-from-opening instead would overstate the narrowing, because the
# centerline length changes frame to frame and a fixed s slides along the anatomy
# (67% vs the true 37% for this case).
# Track the THROAT, not one plane: per frame take the narrowest valid section among the
# bottleneck span's planes, then narrowing = 1 - min/max of that over the breath. A single
# plane at the edge of an abrupt widening slides between throat and funnel as the airway
# moves (ENT005 plane 53: 32 <-> 96 mm2 = a spurious 67%, neighbours 22-28%).
span_planes = {int(p) for p, sv, rg in zip(aper[tref]["pidx"], aper[tref]["s"], aper[tref]["region"])
               if seg_lo - 1.0 <= sv <= seg_hi + 1.0 and str(rg) == constr_reg}
span_planes |= {constr_plane}
print("bottleneck span planes (%s): %s" % (constr_reg, sorted(span_planes)))
per_frame = {}
for r in csv.DictReader(open(a.csa)):
    if r["is_valid"].lower() == "true" and r["region"] == constr_reg and int(r["plane_index"]) in span_planes:
        per_frame[r["frame_name"]] = min(per_frame.get(r["frame_name"], 1e18), float(r["area_mm2"]))
csa_star = np.array(list(per_frame.values()))
csa_min = float(csa_star.min())
dyn = 100.0 * (1.0 - csa_min / float(csa_star.max()))

# the anatomical minimum, reported only if it is somewhere else
Aro = list(csv.DictReader(open(a.anatomy)))
anat_min = float(min(float(r["min_CSA_mm2"]) for r in Aro))
anat_loc = float(np.median([float(r["min_loc_mm_from_opening"]) for r in Aro]))
anat_elsewhere = abs(anat_loc - s_star) > 15.0
anat_dP = float(np.interp(min(anat_loc + 10, S.max()), S, CUM)
                - np.interp(max(anat_loc - 10, 0.0), S, CUM))

print("dP@peak=%.0f Pa (max %.0f)  WOB=%.2f mJ  peak=%.2f mW at t=%.2fs  insp=%.0f%%"
      % (dP_pk, dP_max, wob_tot, peak_pow, tpk, insp_frac))
print("bottleneck by dissipation: s=%.0f mm (%s plane %d), span %.0f-%.0f mm carries "
      "%.0f Pa = %.0f%% of the drop; CSA %.1f mm2, dynamic narrowing %.0f%%"
      % (s_star, constr_reg, constr_plane, seg_lo, seg_hi, seg_dP, seg_frac, csa_min, dyn))
print("anatomical min CSA %.1f mm2 at s=%.0f mm -> drops only %.1f Pa over +-10 mm"
      % (anat_min, anat_loc, anat_dP))

# ---- page -------------------------------------------------------------------------
fig = plt.figure(figsize=(8.5, 11.0), facecolor="white")
LH = lambda fs, ls=1.4: fs * ls / 72.0 / 11.0        # one line height in figure fraction

fig.text(0.5, 0.973, u"Patient-Specific Airway CFD  —  Subject %s" % a.subject,
         ha="center", fontsize=19, fontweight="bold", color=NAVY)
fig.text(0.5, 0.951,
         "Dynamic airway bottleneck & work of breathing at peak %s (t = %.2f s)" % (PHASE, tpk),
         ha="center", fontsize=11.5, color="#444")

img = mpimg.imread(a.snapshot)
ih, iw = img.shape[:2]
if iw > ih:                                   # wide (columns layout): full page width
    W_ = 0.93; H_ = W_ * 8.5 / 11.0 * ih / iw; X_ = 0.035; Y_ = 0.935 - H_
else:                                         # tall (stacked layout, ENT002/003)
    W_, H_, X_ = 0.66, 0.455, 0.17; Y_ = 0.468
axi = fig.add_axes([X_, Y_, W_, H_]); axi.imshow(img); axi.axis("off")
axi.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0", transform=axi.transAxes,
                             fill=False, ec="#cccccc", lw=1.0))
cap = (("Snapshot at the bottleneck — Top: moving airway wall + section planes, centerline locator, "
        "velocity field.  Below each: flow phase, cross-sectional area, and work of breathing along the "
        "%s → %s route.") if iw > ih else
       ("Snapshot at the bottleneck — Left: moving airway wall + section planes, centerline "
        "locator, velocity field.  Right: flow phase, cross-sectional area, and work of "
        "breathing along the %s → %s route.")) % (PROX, DIST)
ycap = Y_ - 0.006
for ln in textwrap.fill(cap, 124).split("\n"):
    fig.text(0.5, ycap, ln, ha="center", va="top", fontsize=8.0, color="#555", style="italic")
    ycap -= LH(8.0, 1.25)

YB = ycap - 0.025                          # banner just under the caption
fig.patches.append(plt.Rectangle((0.06, YB), 0.88, 0.028, transform=fig.transFigure,
                                 facecolor=NAVY, edgecolor="none"))
fig.text(0.072, YB + 0.014, "The value of this analysis — what the CFD reveals for %s" % a.subject,
         va="center", fontsize=12.5, fontweight="bold", color="white")

jet = ((" Flow accelerates into a jet of ~%.0f m/s through it," % a.peak_velocity)
       if a.peak_velocity else " Flow accelerates sharply through it,")
_where = ("at the %s" % PROX) if anat_loc < 5 else ("about %.0f mm from the %s" % (anat_loc, PROX))
EQUAL = anat_min >= 0.85 * csa_min          # the tightest section is about as narrow as the bottleneck
if not anat_elsewhere:
    decoy = ""
elif EQUAL:
    decoy = (" Notably a second section %s is just as tight (~%.0f mm²) but drops only %.0f Pa: on imaging the "
             "two look alike, yet only one carries the airway's resistive load." % (_where, anat_min, anat_dP))
else:
    decoy = (" Notably this is NOT the narrowest cross-section: the tightest section (~%.0f mm², %s) drops "
             "only %.0f Pa, so ranking the airway by cross-sectional area alone would have pointed at the "
             "wrong segment." % (anat_min, _where, anat_dP))

find = (
    "FINDING (this patient).  Ranking the airway by where total pressure is actually lost, the "
    "dominant bottleneck on the %s → %s route sits about %.0f mm from the %s (%s, plane index "
    "%d). The %.0f-%.0f mm span across it carries %.0f Pa — %.0f%% of the entire %s-to-%s drop "
    "— while its cross-section falls to ~%.0f mm² and tightens a further ~%.0f%% through the "
    "breath, a dynamic narrowing static imaging cannot show.%s which concentrates the airway's "
    "resistive load: a total-pressure drop of ~%.0f Pa at peak %s and a resistive work "
    "of breathing of %.1f mJ across the recorded breath (%s), "
    "%.0f%% of it inspiratory and peaking at %.0f mW.%s"
) % (PROX, DIST, s_star, PROX, constr_reg, constr_plane, seg_lo, seg_hi, seg_dP, seg_frac,
     PROX, DIST, csa_min, dyn, jet, dP_pk, PHASE, wob_tot, ORDER, insp_frac, peak_pow, decoy)

why = [
    u"Localises the dominant resistance to one airway segment — a concrete surgical target, "
    u"not just “a narrowing somewhere.”",
    ("Ranks segments by pressure actually lost, not by how narrow they look: here two sections are "
     "about equally tight, but only one is the true bottleneck." if anat_elsewhere and EQUAL else
     "Ranks segments by pressure actually lost, not by how narrow they look: here the tightest "
     "cross-section and the true bottleneck are in different places." if anat_elsewhere else
     "Ranks segments by pressure actually lost, not by how narrow they look: here the tightest "
     "cross-section is also where the most pressure is lost, confirming it as the target."),
    "Quantifies the functional cost: work of breathing (mW, mJ/breath) and along-airway "
    "pressure loss that anatomy alone cannot give.",
    "Gives an objective pre-operative baseline and a virtual-surgery platform to predict the "
    "benefit of relieving the constriction.",
]

def layout(scale, draw):
    """Lay out the finding + bullets at *scale*; return the final y. The ENT003 finding
    runs longer than ENT002's, so a fixed font size ran the BOTTOM LINE box off the page
    -- pick the largest scale that still fits instead."""
    fs_f, fs_h, fs_b = 10.4 * scale, 11.2 * scale, 10.2 * scale
    wrap_f, wrap_b = int(104 / scale), int(96 / scale)
    y = YB - 0.017
    for ln in textwrap.fill(find, wrap_f).split("\n"):
        if draw:
            fig.text(0.072, y, ln, va="top", ha="left", fontsize=fs_f, color=INK)
        y -= LH(fs_f)
    y -= 0.014
    if draw:
        fig.text(0.072, y, "WHY IT MATTERS  (beyond CT / endoscopy)", va="top",
                 fontsize=fs_h, fontweight="bold", color=NAVY)
    y -= LH(fs_h) + 0.006
    for b in why:
        for i, ln in enumerate(textwrap.fill(b, wrap_b).split("\n")):
            if draw:
                if i == 0:
                    fig.text(0.082, y, u"•", va="top", fontsize=fs_h, color=ACC,
                             fontweight="bold")
                fig.text(0.104, y, ln, va="top", ha="left", fontsize=fs_b, color=INK)
            y -= LH(fs_b)
        y -= 0.006
    return y


FLOOR = 0.088                      # leave room for the BOTTOM LINE box
scale = 1.0
while scale > 0.80 and layout(scale, False) < FLOOR:
    scale -= 0.02
y = layout(scale, True)
if scale < 1.0:
    print("text block scaled to %.0f%% to fit the page" % (100 * scale))

y -= 0.006
fig.patches.append(plt.Rectangle((0.06, y - 0.052), 0.88, 0.050, transform=fig.transFigure,
                                 facecolor="#eef2f7", edgecolor="#cbd6e2", lw=1.0))
for ln in textwrap.fill(
        "BOTTOM LINE.  CFD turns moving-airway anatomy + airflow into quantitative, localized, "
        "breath-resolved function — pinpointing where and how much resistance there is, to guide "
        "and justify intervention.", 108).split("\n"):
    fig.text(0.5, y - 0.014, ln, ha="center", va="top", fontsize=9.4, color=NAVY,
             fontweight="bold")
    y -= LH(9.4)

os.makedirs(os.path.dirname(a.out), exist_ok=True)
from matplotlib.backends.backend_pdf import PdfPages
with PdfPages(a.out) as pdf:
    pdf.savefig(fig, facecolor="white", dpi=300)
    if a.regions_csv and a.regions_png and os.path.exists(a.regions_csv):
        # ---- page 2: where and when the work of breathing occurs ----
        R = list(csv.DictReader(open(a.regions_csv)))
        bkey = [k for k in R[0] if k.startswith("bounds")][0]
        nm = [r["region"] for r in R]; bd = [r[bkey] for r in R]
        Ei = np.array([float(r["inhale_mJ"]) for r in R]); Ee = np.array([float(r["exhale_mJ"]) for r in R])
        tot = Ei.sum() + Ee.sum(); ii, ie = int(np.argmax(Ei)), int(np.argmax(Ee))
        f2 = plt.figure(figsize=(8.5, 11.0), facecolor="white")
        f2.text(0.5, 0.962, u"Where and when the work of breathing occurs  —  %s" % a.subject, ha="center",
                fontsize=16, fontweight="bold", color=NAVY)
        f2.text(0.5, 0.940, "resistive work by airway region, inhalation vs exhalation (%s → %s route)" % (PROX, DIST),
                ha="center", fontsize=10.5, fontweight="bold", color="#444")
        im2 = mpimg.imread(a.regions_png); h2, w2 = im2.shape[:2]; W2 = 0.95; H2 = W2 * 8.5 / 11.0 * h2 / w2
        ax2 = f2.add_axes([0.025, 0.925 - H2, W2, H2]); ax2.imshow(im2); ax2.axis("off")
        y2 = 0.925 - H2 - 0.03
        f2.patches.append(plt.Rectangle((0.06, y2 - 0.005), 0.88, 0.028, transform=f2.transFigure,
                                        facecolor=NAVY, edgecolor="none"))
        f2.text(0.072, y2 + 0.009, "What the regional breakdown shows", va="center", fontsize=12.5,
                fontweight="bold", color="white")
        bullets = [
            "Exhalation carries %.0f%% of the resistive work (%.0f of %.0f mJ over the recorded breath); "
            "inhalation %.0f%% (%.0f mJ)." % (100 * Ee.sum() / tot, Ee.sum(), tot, 100 * Ei.sum() / tot, Ei.sum()),
            "Exhalation: the %s (%s mm from the %s) costs the most, %.0f mJ = %.0f%% of the exhalation work."
            % (nm[ie], bd[ie], PROX, Ee[ie], 100 * Ee[ie] / Ee.sum()),
            "Inhalation: the %s (%s mm) costs the most, %.0f mJ = %.0f%% of the inhalation work."
            % (nm[ii], bd[ii], Ei[ii], 100 * Ei[ii] / Ei.sum()),
        ]
        if "nose" in nm:
            k = nm.index("nose")
            bullets.append("The nose (both nostrils, opening to the top of the shared airway) costs %.0f mJ, "
                           "%.0f%% of the breath." % (Ei[k] + Ee[k], 100 * (Ei[k] + Ee[k]) / tot))
        bullets += [
            "Work is counted where the energy is dissipated: a jet leaving a narrowing loses most of its energy "
            "just downstream, so the region next to a constriction can carry more work than the constriction itself.",
            "Region boundaries (mm along the centerline from the %s opening): %s."
            % (PROX, "; ".join("%s %s" % (n, b) for n, b in zip(nm, bd))),
        ]
        yb = y2 - 0.04
        for b in bullets:
            for i, ln in enumerate(textwrap.fill(b, 92).split("\n")):
                if i == 0:
                    f2.text(0.072, yb, u"•", va="top", fontsize=11.5, color=ACC, fontweight="bold")
                f2.text(0.094, yb, ln, va="top", fontsize=11, color=INK)
                yb -= LH(11)
            yb -= 0.008
        pdf.savefig(f2, facecolor="white", dpi=300)
        print("page 2: regional work of breathing added")
print("wrote", a.out)
