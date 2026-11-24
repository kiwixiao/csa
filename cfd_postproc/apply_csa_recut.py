#!/usr/bin/env python3
"""Carry a virtual edit onto a subject's CSA table: area x (edited / native) from recut_csa.py.

recut_csa.py cuts the native and the edited CFD wall with the SAME planes; their ratio is the edit
(meshing and method differences cancel). Applied only over the edited span, faded to 1 over --taper-mm
at both ends, so planes where the two walls are the same keep the slicer's area exactly.
--reference table divides by the slicer's own area instead of the native re-cut: needed when the native
CFD wall does not follow the motion surfaces (ENT005: its mesh motion drifts, area down to ~40% by the end
of the breath, while the Aryt wall matches the motion surfaces within ~3%).

Usage: apply_csa_recut.py --csa <S>-Data-clean.csv --native recut_native.csv --edited recut_edited.csv
           --arc-lo 88 --arc-hi 123 --out <S>-Data-clean.csv(edited copy)
"""
import argparse, csv
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--csa", required=True); ap.add_argument("--native", required=True); ap.add_argument("--edited", required=True)
ap.add_argument("--arc-lo", type=float, required=True); ap.add_argument("--arc-hi", type=float, required=True)
ap.add_argument("--taper-mm", type=float, default=5.0); ap.add_argument("--out", required=True)
ap.add_argument("--min-factor", type=float, default=None,
                help="floor on the area factor: 1 when the edit only widens (a lower re-cut is a missed lumen piece "
                     "around the epiglottis, not the edit)")
ap.add_argument("--reference", choices=("native", "table"), default="native",
                help="denominator of the area factor: the native re-cut, or the CSA table's own area")
a = ap.parse_args()


def recut(path, col="area_recut_mm2"):
    return {(r["frame"], r["region"], r["plane_index"]): float(r[col]) for r in csv.DictReader(open(path))}


E = recut(a.edited)
N = recut(a.native) if a.reference == "native" else recut(a.edited, "area_table_mm2")
rows = list(csv.DictReader(open(a.csa))); cols = list(rows[0].keys())
done, ratios = 0, []
for r in rows:
    k = (r["frame_name"], r["region"], r["plane_index"])
    if k not in N or k not in E or N[k] <= 0 or E[k] <= 0:
        continue
    arc = float(r["arc_length_mm"])
    w = np.clip(min(arc - (a.arc_lo - a.taper_mm), (a.arc_hi + a.taper_mm) - arc) / a.taper_mm, 0.0, 1.0)
    if w <= 0:
        continue
    f = E[k] / N[k] if a.min_factor is None else max(E[k] / N[k], a.min_factor)
    f = 1.0 + w * (f - 1.0)
    A = float(r["area_mm2"]) * f
    r["area_mm2"] = "%.6f" % A
    r["equivalent_diameter_mm"] = "%.6f" % np.sqrt(4 * A / np.pi)
    done += 1; ratios.append(f)
with open(a.out, "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=cols); wr.writeheader(); wr.writerows(rows)
print("edited %d rows (arc %.0f-%.0f mm, taper %.0f mm, vs %s): area factor median %.3f, max %.3f -> %s"
      % (done, a.arc_lo, a.arc_hi, a.taper_mm, a.reference, np.median(ratios), np.max(ratios), a.out))
