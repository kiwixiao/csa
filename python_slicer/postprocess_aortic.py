#!/usr/bin/env python3
"""
Post-processing for aortic CSA pipeline (single time point).

Generates per-branch analysis: CSA plots, circularity, plane reference PNG,
PDF report. No temporal analysis (single frame).

Usage:
    python python_slicer/postprocess_aortic.py AorticSubject/csa/ --subject-id AORT001
"""

import sys
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

log = logging.getLogger(__name__)


def plot_csa_per_branch(df, output_dir, subject_id):
    """CSA by plane index, one plot per branch."""
    regions = sorted(df["region"].unique())

    for region in regions:
        rdf = df[df["region"] == region].sort_values("plane_index")
        if len(rdf) == 0:
            continue

        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        axes[0].plot(rdf["plane_index"], rdf["area_mm2"],
                    "-o", color="steelblue", markersize=3, linewidth=1.5)
        axes[0].set_ylabel("CSA (mm²)")
        axes[0].set_title(f"{subject_id} — {region} — CSA")
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(rdf["plane_index"], rdf["hydraulic_diameter_mm"],
                    "-o", color="coral", markersize=3, linewidth=1.5)
        axes[1].set_ylabel("Hydraulic Diameter (mm)")
        axes[1].set_xlabel("Plane Index")
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = output_dir / f"{subject_id}_{region}_CSA.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info(f"  {out_path.name}")


def compute_circularity(df, output_dir, subject_id):
    """Add circularity metric and plot per branch."""
    df = df.copy()
    df["circularity"] = (4 * np.pi * df["area_mm2"]) / (df["perimeter_mm"] ** 2)
    df.loc[df["perimeter_mm"] <= 0, "circularity"] = 0

    csv_path = output_dir / f"{subject_id}_enhanced_metrics.csv"
    df.to_csv(csv_path, index=False)
    log.info(f"  Enhanced metrics: {csv_path.name}")

    for region in sorted(df["region"].unique()):
        rdf = df[df["region"] == region].sort_values("plane_index")
        if len(rdf) == 0:
            continue

        fig, ax = plt.subplots(figsize=(12, 5))
        ax.plot(rdf["plane_index"], rdf["circularity"],
               "-o", color="steelblue", markersize=3, linewidth=1.5)
        ax.set_xlabel("Plane Index")
        ax.set_ylabel("Circularity (1.0 = circle)")
        ax.set_title(f"{subject_id} — {region} — Circularity")
        ax.set_ylim(0, 1.1)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = output_dir / f"{subject_id}_{region}_circularity.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info(f"  {out_path.name}")

    return df


def generate_plane_reference_aortic(df, output_dir, subject_id):
    """Generate plane reference PNG with labeled indices per branch."""
    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv
    from PIL import Image
    subject_dir = output_dir.parent
    surface_stls = list((subject_dir / "surface").glob("*.stl"))
    if not surface_stls:
        return

    regions = sorted(df["region"].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(regions), 1)))
    color_map = {r: colors[i][:3] for i, r in enumerate(regions)}

    panels = []
    for region in regions:
        rdf = df[df["region"] == region].sort_values("plane_index")
        if len(rdf) == 0:
            continue

        # Find plane STLs
        planes_dir = output_dir / f"{region}_Planes"
        color_rgb = color_map[region]
        color_hex = "#{:02x}{:02x}{:02x}".format(
            int(color_rgb[0]*255), int(color_rgb[1]*255), int(color_rgb[2]*255))

        plotter = pv.Plotter(off_screen=True, window_size=(1000, 1600))
        plotter.set_background('white')

        # Full mesh (transparent)
        plotter.add_mesh(pv.read(str(surface_stls[0])),
                        color='lightgray', opacity=0.1, smooth_shading=True)

        # Branch planes
        if planes_dir.exists():
            for stl_file in sorted(planes_dir.glob('*.stl')):
                plotter.add_mesh(pv.read(str(stl_file)),
                                color=color_hex, opacity=0.5)

        # Labels every 2nd plane
        for _, row in rdf.iterrows():
            pi = row['plane_index']
            if pi % 2 != 0:
                continue
            centroid = [row['centroid_x'], row['centroid_y'], row['centroid_z']]
            plotter.add_point_labels(
                np.array([centroid]), [f"{pi}"],
                font_size=14, text_color=color_hex,
                point_size=0, shape=None, always_visible=True)

        plotter.add_text(region, position='upper_edge', font_size=12,
                        color=color_hex)

        plotter.view_yz()
        plotter.camera.up = (0, 0, 1)
        plotter.enable_parallel_projection()
        plotter.reset_camera()
        pos = list(plotter.camera_position)
        pos[0] = (-pos[0][0], pos[0][1], pos[0][2])
        plotter.camera_position = pos
        plotter.camera.zoom(0.85)

        img = plotter.screenshot(return_img=True)
        plotter.close()
        panels.append(img[:, :, :3])

    if not panels:
        return

    max_h = max(p.shape[0] for p in panels)
    resized = []
    for p in panels:
        pil = Image.fromarray(p)
        new_w = int(pil.width * max_h / pil.height)
        pil = pil.resize((new_w, max_h), Image.LANCZOS)
        resized.append(np.array(pil))

    composite = np.hstack(resized)
    out_path = output_dir / f"{subject_id}_plane_reference.png"
    Image.fromarray(composite).save(str(out_path))
    log.info(f"  Plane reference: {out_path.name}")


def generate_combined_pdf(df, output_dir, subject_id):
    """Generate one master PDF with all analysis: sanity check, CSA, circularity, summary."""
    from PIL import Image as PILImage

    pdf_path = output_dir / f"{subject_id}_report.pdf"
    regions = sorted(df["region"].unique())

    with PdfPages(str(pdf_path)) as pdf:
        # Page 1: Sanity check PNG
        sanity_png = output_dir / f"{subject_id}_sanity_check.png"
        if sanity_png.exists():
            fig, ax = plt.subplots(figsize=(10, 14))
            ax.imshow(PILImage.open(str(sanity_png)))
            ax.set_title(f"{subject_id} — Geometry + Planes Sanity Check", fontsize=14)
            ax.axis('off')
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

        # Page 2: CSA plot PNG
        csa_png = output_dir / f"{subject_id}-CSA_plot.png"
        if csa_png.exists():
            fig, ax = plt.subplots(figsize=(16, 10))
            ax.imshow(PILImage.open(str(csa_png)))
            ax.axis('off')
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

        # Page 3: Plane reference PNG
        ref_png = output_dir / f"{subject_id}_plane_reference.png"
        if ref_png.exists():
            fig, ax = plt.subplots(figsize=(14, 10))
            ax.imshow(PILImage.open(str(ref_png)))
            ax.set_title(f"{subject_id} — Plane Reference", fontsize=14)
            ax.axis('off')
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

        # Per-region pages
        for region in regions:
            rdf = df[df["region"] == region].sort_values("plane_index")
            if len(rdf) == 0:
                continue

            # CSA + diameter page
            fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
            axes[0].plot(rdf["arc_length_mm"], rdf["area_mm2"],
                        "-o", color="steelblue", markersize=3, linewidth=1.5)
            axes[0].set_ylabel("CSA (mm²)")
            axes[0].set_title(f"{subject_id} — {region} — Cross-Sectional Area")
            axes[0].grid(True, alpha=0.3)
            axes[1].plot(rdf["arc_length_mm"], rdf["hydraulic_diameter_mm"],
                        "-o", color="coral", markersize=3, linewidth=1.5)
            axes[1].set_ylabel("Hydraulic Diameter (mm)")
            axes[1].set_xlabel("Arc Length (mm)")
            axes[1].grid(True, alpha=0.3)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

            # Circularity page
            if "circularity" in rdf.columns:
                fig, ax = plt.subplots(figsize=(12, 5))
                ax.plot(rdf["arc_length_mm"], rdf["circularity"],
                       "-o", color="steelblue", markersize=3, linewidth=1.5)
                ax.set_xlabel("Arc Length (mm)")
                ax.set_ylabel("Circularity (1.0 = circle)")
                ax.set_title(f"{region} — Cross-Section Circularity")
                ax.set_ylim(0, 1.1)
                ax.grid(True, alpha=0.3)
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches="tight")
                plt.close()

            # Summary page
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.axis('off')
            stats_text = (
                f"Region: {region}\n"
                f"Planes: {len(rdf)}\n"
                f"Arc length: {rdf['arc_length_mm'].min():.1f} - "
                f"{rdf['arc_length_mm'].max():.1f} mm\n"
                f"Area: min={rdf['area_mm2'].min():.2f}, "
                f"max={rdf['area_mm2'].max():.2f}, "
                f"mean={rdf['area_mm2'].mean():.2f} mm²\n"
                f"Hyd. diameter: min={rdf['hydraulic_diameter_mm'].min():.2f}, "
                f"max={rdf['hydraulic_diameter_mm'].max():.2f}, "
                f"mean={rdf['hydraulic_diameter_mm'].mean():.2f} mm\n"
            )
            if "circularity" in rdf.columns:
                stats_text += f"Circularity: mean={rdf['circularity'].mean():.3f}\n"
            ax.text(0.1, 0.9, stats_text, transform=ax.transAxes,
                   fontsize=12, verticalalignment='top', fontfamily='monospace')
            ax.set_title(f"{subject_id} — {region} — Summary", fontsize=14)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

    log.info(f"  Combined report: {pdf_path.name}")


def generate_sanity_check(df, output_dir, subject_id):
    """STL + centerlines + planes colored by branch — sanity check PNG."""
    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv

    subject_dir = output_dir.parent
    surface_stls = list((subject_dir / "surface").glob("*.stl"))
    branches_dir = subject_dir / "branches"
    if not surface_stls or not branches_dir.exists():
        return

    regions = sorted(df["region"].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(regions), 1)))
    color_map = {r: "#{:02x}{:02x}{:02x}".format(
        int(colors[i][0]*255), int(colors[i][1]*255), int(colors[i][2]*255))
        for i, r in enumerate(regions)}

    # Branch label to folder name mapping for centerlines
    branch_cl_colors = {}
    for bd in sorted(branches_dir.iterdir()):
        if bd.is_dir():
            branch_cl_colors[bd.name] = 'gray'

    pl = pv.Plotter(off_screen=True, window_size=(1200, 1600))
    pl.set_background('white')
    pl.add_mesh(pv.read(str(surface_stls[0])), color='lightblue', opacity=0.1,
                smooth_shading=True)

    # Centerlines
    for bd in sorted(branches_dir.iterdir()):
        if not bd.is_dir():
            continue
        cl_files = list(bd.glob('*_centerline.vtk'))
        if cl_files:
            pl.add_mesh(pv.read(str(cl_files[0])), color='darkgray', line_width=3)

    # Planes colored by region
    for region in regions:
        pd2 = output_dir / f"{region}_Planes"
        if not pd2.exists():
            continue
        c = color_map.get(region, 'gray')
        for sf in sorted(pd2.glob('*.stl')):
            pl.add_mesh(pv.read(str(sf)), color=c, opacity=0.5, show_edges=False)

    pl.view_yz()
    pl.camera.up = (0, 0, 1)
    pl.enable_parallel_projection()
    pl.reset_camera()
    pos = list(pl.camera_position)
    pos[0] = (-pos[0][0], pos[0][1], pos[0][2])
    pl.camera_position = pos
    pl.camera.zoom(0.85)
    out_path = output_dir / f"{subject_id}_sanity_check.png"
    pl.screenshot(str(out_path))
    pl.close()
    log.info(f"  Sanity check: {out_path.name}")


def run_aortic_postprocessing(output_dir, subject_id):
    """Run all aortic post-processing."""
    output_dir = Path(output_dir)
    csv_path = output_dir / f"{subject_id}-Data.csv"
    if not csv_path.exists():
        log.error(f"CSV not found: {csv_path}")
        return

    df = pd.read_csv(csv_path)
    log.info(f"Loaded {len(df)} rows from {csv_path.name}")
    log.info(f"  Regions: {df['region'].value_counts().to_dict()}")

    # Per-branch CSA plots
    log.info("\nGenerating CSA plots...")
    plot_csa_per_branch(df, output_dir, subject_id)

    # Circularity
    log.info("\nComputing circularity...")
    df = compute_circularity(df, output_dir, subject_id)

    # CoA visualization + clinical metrics
    log.info("\nGenerating CoA visualizations...")
    try:
        from visualize_coa import render_all_views
        subject_dir = output_dir.parent
        surface_stls = list((subject_dir / "surface").glob("*.stl"))
        if surface_stls:
            stl_path = [s for s in surface_stls if not s.name.startswith("._")][0]
            branches_dir = subject_dir / "branches"
            coa_metrics = render_all_views(df, str(stl_path), str(branches_dir),
                                           output_dir, subject_id)
    except Exception as e:
        log.warning(f"  CoA visualization failed: {e}")

    # Sanity check PNG (STL + centerlines + planes)
    log.info("\nGenerating sanity check...")
    generate_sanity_check(df, output_dir, subject_id)

    # Plane reference PNG
    log.info("\nGenerating plane reference...")
    generate_plane_reference_aortic(df, output_dir, subject_id)

    # Combined PDF report (all plots + per-region analysis)
    log.info("\nGenerating combined PDF report...")
    generate_combined_pdf(df, output_dir, subject_id)

    log.info(f"\nAortic post-processing complete.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Aortic CSA post-processing")
    parser.add_argument("output_dir", help="CSA output directory")
    parser.add_argument("--subject-id", default="")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    subject_id = args.subject_id or output_dir.parent.name

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    run_aortic_postprocessing(output_dir, subject_id)
