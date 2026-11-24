#!/usr/bin/env python3
"""Summary + figure of the CFD mesh-motion check (qa/wall_track.py output).

Uses the control points that sit on the CFD wall at the first saved step (< --on-wall mm): for an edited geometry
the points inside the edit are off the wall from the start and are excluded. PASS when, over the whole breath,
the median stays <= --pass-median and the worst point <= --pass-max (mm). Exit status 1 if any run fails.
Usage: python3 motion_check_fig.py --track-dir DIR --runs TAG[=label] [TAG[=label] ...] --out fig.png
"""
import argparse, glob, os, sys
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(); ap.add_argument("--track-dir", required=True); ap.add_argument("--runs", nargs="+", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--on-wall", type=float, default=0.1)
ap.add_argument("--pass-median", type=float, default=0.1); ap.add_argument("--pass-max", type=float, default=1.0)
a = ap.parse_args()
fig, axs = plt.subplots(1, 2, figsize=(15, 5.6)); fig.subplots_adjust(left=0.06, right=0.98, top=0.82, bottom=0.12, wspace=0.18)
ok, lines = True, []
for i, spec in enumerate(a.runs):
    tag, lab = (spec.split("=", 1) + [spec])[:2]
    D = [np.load(f) for f in sorted(glob.glob(os.path.join(a.track_dir, "track_%s_*.npz" % tag)))]
    if not D:
        print("%s: no track files" % tag); ok = False; continue
    keep = D[0]["dc"] < a.on_wall
    t = np.array([int(z["ms"]) for z in D]) / 1000
    med = np.array([np.median(z["dc"][keep]) for z in D]); p95 = np.array([np.percentile(z["dc"][keep], 95) for z in D])
    mx = np.array([z["dc"][keep].max() for z in D])
    good = med.max() <= a.pass_median and mx.max() <= a.pass_max; ok &= good
    lines.append("%s: %d steps, %d/%d points on the wall at the first step; worst median %.3f, worst 95th pct %.3f, worst point %.3f mm -> %s"
                 % (lab, len(D), keep.sum(), len(keep), med.max(), p95.max(), mx.max(), "PASS" if good else "FAIL"))
    c = "C%d" % i
    axs[0].fill_between(t, med, p95, color=c, alpha=0.18, lw=0); axs[0].plot(t, med, color=c, lw=2.2, label=lab)
    axs[1].plot(t, mx, color=c, lw=2.2, label=lab)
for ax, ttl in zip(axs, ("median distance (band to 95th percentile)", "worst point (max distance)")):
    ax.set_yscale("log"); ax.set_ylim(1e-3, 30); ax.grid(alpha=0.3, which="both"); ax.set_xlabel("time (s)")
    ax.set_yticks([0.001, 0.01, 0.1, 1, 10]); ax.set_yticklabels(["0.001", "0.01", "0.1", "1", "10"])
    ax.set_ylabel("registration point to CFD wall (mm)"); ax.set_title(ttl, fontsize=11, weight="bold")
axs[0].axhline(a.pass_median, color="k", ls=":", lw=1); axs[1].axhline(a.pass_max, color="k", ls=":", lw=1)
axs[0].legend(fontsize=9, loc="upper left")
fig.text(0.06, 0.965, "CFD mesh-motion check: does the CFD wall follow the registration motion?", fontsize=15, weight="bold", va="top")
fig.text(0.06, 0.905, "\n".join(lines), fontsize=9, va="top", color="0.2")
fig.savefig(a.out, dpi=130); print("\n".join(lines)); print(a.out)
sys.exit(0 if ok else 1)
