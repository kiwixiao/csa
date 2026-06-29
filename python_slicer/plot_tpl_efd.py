#!/usr/bin/env python3
"""
Plot regional TPL and EFDiff as percentage of total across breathing cycle.

Standalone script — reads CFD output CSVs directly.

Usage:
    python plot_tpl_efd.py /path/to/subject/ --subject-id ENT001
    python plot_tpl_efd.py /path/to/subject/  # auto-detect subject ID
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path


def clean_column_name(col):
    """Extract region name from STAR-CCM+ column header.
    e.g. 'tpl_Larynx_Monitor_pts: tpl_Larynx_Monitor_pts' -> 'Larynx'
    """
    # Take the part before '_Monitor_pts'
    name = col.split("_Monitor_pts")[0]
    # Remove prefix (tpl_ or efd_)
    for prefix in ("tpl_", "efd_"):
        if name.startswith(prefix):
            name = name[len(prefix):]
    return name


def load_and_clean(csv_path, prefix):
    """Load CSV, clean column names, return DataFrame with Time + region columns."""
    df = pd.read_csv(csv_path)
    rename = {}
    for col in df.columns:
        if col.strip().lower() == "time":
            rename[col] = "Time"
        elif col.startswith(prefix):
            rename[col] = clean_column_name(col)
    df = df.rename(columns=rename)
    # Keep only renamed columns
    keep = [c for c in df.columns if c in rename.values()]
    return df[keep]


def main():
    parser = argparse.ArgumentParser(
        description="Plot regional TPL and EFDiff as % of total"
    )
    parser.add_argument("subject_dir", help="Subject directory")
    parser.add_argument("--subject-id", default="")
    parser.add_argument("--tpl-csv", default="", help="Override TPL CSV path")
    parser.add_argument("--efd-csv", default="", help="Override EFDiff CSV path")
    parser.add_argument("--output", default="", help="Output PDF path")
    parser.add_argument("--switch-time", type=float, required=True,
                        help="Time (s) when inhale switches to exhale")
    args = parser.parse_args()

    subject_dir = Path(args.subject_dir).resolve()
    subject_id = args.subject_id or subject_dir.name

    # Find CSV files
    tpl_path = Path(args.tpl_csv) if args.tpl_csv else None
    efd_path = Path(args.efd_csv) if args.efd_csv else None

    if tpl_path is None:
        candidates = list(subject_dir.glob("*TPL*.csv")) + list(subject_dir.glob("*tpl*.csv"))
        if candidates:
            tpl_path = candidates[0]
    if efd_path is None:
        candidates = list(subject_dir.glob("*EFDiff*.csv")) + list(subject_dir.glob("*efdiff*.csv"))
        if candidates:
            efd_path = candidates[0]

    if tpl_path is None or not tpl_path.exists():
        print(f"ERROR: TPL CSV not found in {subject_dir}")
        sys.exit(1)
    if efd_path is None or not efd_path.exists():
        print(f"ERROR: EFDiff CSV not found in {subject_dir}")
        sys.exit(1)

    print(f"TPL:  {tpl_path}")
    print(f"EFD:  {efd_path}")

    # Load data
    tpl = load_and_clean(tpl_path, "tpl_")
    efd = load_and_clean(efd_path, "efd_")

    time_tpl = tpl["Time"].values
    time_efd = efd["Time"].values

    # Identify regions
    tpl_regions = [c for c in tpl.columns if c not in ("Time", "Nose")]
    efd_regions = [c for c in efd.columns if c != "Time"]

    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
              "#ff7f00", "#a65628", "#f781bf", "#999999"]

    output_path = Path(args.output) if args.output else subject_dir / f"{subject_id}_WoB_report.pdf"
    t_switch = args.switch_time

    # --- Compute Work of Breathing by integrating |EFDiff| over time ---
    # WoB (Joules) = ∫ |EFDiff| dt, then average power (W) = WoB / period
    inhale_mask = time_efd <= t_switch
    exhale_mask = time_efd > t_switch
    t_total = time_efd[-1] - time_efd[0]
    t_inhale = t_switch - time_efd[0]
    t_exhale = time_efd[-1] - t_switch

    def integrate_region(values, time, mask):
        """Trapezoid integration of |values| over masked time range."""
        t_masked = time[mask]
        v_masked = np.abs(values[mask])
        if len(t_masked) < 2:
            return 0.0
        from scipy.integrate import trapezoid
        return trapezoid(v_masked, t_masked)

    # Compute per-region: total energy (J), inhale energy, exhale energy
    rows = []
    for r in efd_regions:
        vals = efd[r].values
        e_total = integrate_region(vals, time_efd, np.ones(len(time_efd), dtype=bool))
        e_inhale = integrate_region(vals, time_efd, inhale_mask)
        e_exhale = integrate_region(vals, time_efd, exhale_mask)
        rows.append({
            "region": r,
            "total_J": e_total,
            "inhale_J": e_inhale,
            "exhale_J": e_exhale,
        })

    wob = pd.DataFrame(rows).set_index("region")
    wob = wob.sort_values("total_J", ascending=True)
    total_e = wob["total_J"].sum()
    inhale_e = wob["inhale_J"].sum()
    exhale_e = wob["exhale_J"].sum()

    print(f"\nWork of Breathing (Energy per Cycle):")
    print(f"  Total:  {total_e:.4f} J")
    print(f"  Inhale: {inhale_e:.4f} J (0–{t_switch}s)")
    print(f"  Exhale: {exhale_e:.4f} J ({t_switch}–{time_efd[-1]:.1f}s)")

    # --- Bar chart: 3 groups (Total, Inhale, Exhale) ---
    fig, (ax_total, ax_phases) = plt.subplots(
        1, 2, figsize=(18, 7),
        gridspec_kw={"width_ratios": [1.0, 1.5]},
    )
    fig.suptitle(
        f"{subject_id} — Work of Breathing\n"
        f"Total = {total_e:.4f} J  |  "
        f"Inhale = {inhale_e:.4f} J (0–{t_switch}s)  |  "
        f"Exhale = {exhale_e:.4f} J ({t_switch}–{time_efd[-1]:.1f}s)",
        fontsize=13, fontweight="bold",
    )

    # Left: Total WoB
    bar_colors = [colors[efd_regions.index(r) % len(colors)]
                  if r in efd_regions else "gray"
                  for r in wob.index]
    ax_total.barh(wob.index, wob["total_J"], color=bar_colors, alpha=0.85)
    for i, (r, row) in enumerate(wob.iterrows()):
        pct = row["total_J"] / total_e * 100 if total_e > 0 else 0
        ax_total.text(row["total_J"] + wob["total_J"].max() * 0.02, i,
                      f"{row['total_J']:.4f} J ({pct:.0f}%)",
                      va="center", fontsize=9)
    ax_total.set_xlabel("Energy per Cycle (J)", fontsize=11, fontweight="bold")
    ax_total.set_title("Total (Full Cycle)", fontsize=11)
    ax_total.grid(True, alpha=0.3, axis="x")
    ax_total.set_xlim(0, wob["total_J"].max() * 1.4)
    ax_total.tick_params(axis="y", labelsize=11)

    # Right: Inhale vs Exhale grouped bars
    y_pos = np.arange(len(wob))
    bar_h = 0.35
    ax_phases.barh(y_pos + bar_h / 2, wob["inhale_J"], height=bar_h,
                   color="steelblue", alpha=0.8, label=f"Inhale (0–{t_switch}s)")
    ax_phases.barh(y_pos - bar_h / 2, wob["exhale_J"], height=bar_h,
                   color="coral", alpha=0.8, label=f"Exhale ({t_switch}–{time_efd[-1]:.1f}s)")

    # Labels: percentage of phase total + raw value in brackets
    for i, (r, row) in enumerate(wob.iterrows()):
        inh_pct = row["inhale_J"] / inhale_e * 100 if inhale_e > 0 else 0
        exh_pct = row["exhale_J"] / exhale_e * 100 if exhale_e > 0 else 0
        ax_phases.text(row["inhale_J"] + wob["inhale_J"].max() * 0.02,
                       i + bar_h / 2,
                       f"{inh_pct:.0f}% ({row['inhale_J']:.4f} J)",
                       va="center", fontsize=8, color="steelblue")
        ax_phases.text(row["exhale_J"] + wob["exhale_J"].max() * 0.02,
                       i - bar_h / 2,
                       f"{exh_pct:.0f}% ({row['exhale_J']:.4f} J)",
                       va="center", fontsize=8, color="coral")

    ax_phases.set_yticks(y_pos)
    ax_phases.set_yticklabels(wob.index, fontsize=11)
    ax_phases.set_xlabel("Energy per Phase (J)", fontsize=11, fontweight="bold")
    ax_phases.set_title("Inhale vs Exhale", fontsize=11)
    ax_phases.legend(fontsize=9, loc="lower right")
    ax_phases.grid(True, alpha=0.3, axis="x")

    # Summary sentence for clinicians
    top_total = wob["total_J"].idxmax()
    top_total_pct = (wob.loc[top_total, "total_J"] / total_e * 100) if total_e > 0 else 0.0

    top_inhale = wob["inhale_J"].idxmax()
    top_inhale_pct = (wob.loc[top_inhale, "inhale_J"] / inhale_e * 100) if inhale_e > 0 else 0.0

    top_exhale = wob["exhale_J"].idxmax()
    top_exhale_pct = (wob.loc[top_exhale, "exhale_J"] / exhale_e * 100) if exhale_e > 0 else 0.0

    inhale_ratio = inhale_e / exhale_e if exhale_e > 0 else 0

    summary = (
        f"The {top_total} accounts for {top_total_pct:.0f}% of total energy dissipation ({total_e:.3f} J/cycle). "
        f"Inhalation requires {inhale_ratio:.1f}× more energy than exhalation. "
        f"During inhale, {top_inhale} dominates ({top_inhale_pct:.0f}%); "
        f"during exhale, {top_exhale} dominates ({top_exhale_pct:.0f}%)."
    )

    fig.text(
        0.5, 0.035,
        summary,
        ha="center", fontsize=9, fontweight="bold", color="black",
        wrap=True,
    )

    fig.text(
        0.5, 0.005,
        f"WoB = ∫|EFDiff|dt (Joules per cycle).  Switch time = {t_switch}s.  "
        f"Inhale + Exhale = Total ({inhale_e:.4f} + {exhale_e:.4f} = {inhale_e+exhale_e:.4f} J).",
        ha="center", fontsize=8, color="gray", style="italic",
    )

    plt.tight_layout(rect=[0, 0.03, 1, 0.90])

    with PdfPages(str(output_path)) as pdf:
        pdf.savefig(fig, bbox_inches="tight", dpi=150)
    plt.close()

    print(f"Report: {output_path}")


if __name__ == "__main__":
    main()
