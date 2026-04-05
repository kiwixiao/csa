#!/usr/bin/env python3
"""
Coarctation of the Aorta (CoA) visualization and clinical metrics.

Three visualization options + clinical quantification:
1. Side-by-side: STL geometry + CSA plot
2. Stacked: STL top, CSA bottom (x-axis aligned)
3. Color-mapped STL: surface colored by local CSA

Clinical metrics:
- Coarctation Index (CI = D_min / D_ref)
- Area stenosis (%)
- Severity classification (mild/moderate/severe)

Usage:
    from visualize_coa import compute_coa_metrics, render_all_views
    metrics = compute_coa_metrics(df)
    render_all_views(df, stl_path, branches_dir, output_dir, subject_id)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from pathlib import Path
import logging

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Clinical metrics
# ---------------------------------------------------------------------------

def compute_coa_metrics(df):
    """Compute CoA clinical metrics from CSA data.

    Anatomically correct locations:
    - Isthmus: region just AFTER Branch_3 (last arch branch = left subclavian)
      This is where coarctation occurs. Measure min CSA in this region.
    - Reference: descending aorta at diaphragm (last 10 planes of MainAorta)
    - CI = D_isthmus / D_reference

    Returns dict with clinical metrics and landmark locations.
    """
    main = df[df["region"] == "MainAorta"].sort_values("arc_length_mm")
    if len(main) == 0:
        return None

    # Find Branch_3 (last arch branch = left subclavian) arc length
    branch_arcs = _get_branch_arc_lengths(df)
    last_branch = max(branch_arcs.values()) if branch_arcs else 0

    # Reference: descending aorta at diaphragm (last 10 planes)
    n_ref = min(10, len(main) // 3)
    ref_planes = main.tail(n_ref)
    ref_csa = ref_planes["area_mm2"].mean()
    ref_diam = ref_planes["hydraulic_diameter_mm"].mean()
    ref_arc = ref_planes["arc_length_mm"].mean()

    # Isthmus region: planes just after last arch branch
    # Take the first ~15 planes after the last branch divergence
    isthmus_start = last_branch
    isthmus_region = main[main["arc_length_mm"] >= isthmus_start].head(15)

    if len(isthmus_region) > 0:
        # CoA measurement: minimum CSA in the isthmus region
        isthmus_min_idx = isthmus_region["area_mm2"].idxmin()
        isthmus_row = isthmus_region.loc[isthmus_min_idx]
        coa_csa = isthmus_row["area_mm2"]
        coa_diam = isthmus_row["hydraulic_diameter_mm"]
        coa_arc = isthmus_row["arc_length_mm"]
    else:
        # Fallback: global min
        isthmus_min_idx = main["area_mm2"].idxmin()
        isthmus_row = main.loc[isthmus_min_idx]
        coa_csa = isthmus_row["area_mm2"]
        coa_diam = isthmus_row["hydraulic_diameter_mm"]
        coa_arc = isthmus_row["arc_length_mm"]

    # Global min CSA (may differ from isthmus min)
    global_min_idx = main["area_mm2"].idxmin()
    global_min_row = main.loc[global_min_idx]

    # Max CSA (ascending)
    max_idx = main["area_mm2"].idxmax()
    max_row = main.loc[max_idx]

    # Coarctation Index (CI = D_isthmus / D_ref)
    ci = coa_diam / ref_diam if ref_diam > 0 else 0

    # Area stenosis at isthmus (%)
    area_stenosis = (1 - coa_csa / ref_csa) * 100 if ref_csa > 0 else 0

    # Severity classification
    if ci > 0.7:
        severity = "Mild"
        severity_color = "green"
    elif ci > 0.5:
        severity = "Moderate"
        severity_color = "orange"
    else:
        severity = "Severe"
        severity_color = "red"

    metrics = {
        "coarctation_index": ci,
        "area_stenosis_pct": area_stenosis,
        "severity": severity,
        "severity_color": severity_color,
        # Isthmus (CoA site)
        "coa_csa_mm2": coa_csa,
        "coa_diameter_mm": coa_diam,
        "coa_arc_length_mm": coa_arc,
        # Global min (may differ from isthmus)
        "min_csa_mm2": global_min_row["area_mm2"],
        "min_arc_length_mm": global_min_row["arc_length_mm"],
        # Max (ascending)
        "max_csa_mm2": max_row["area_mm2"],
        "max_arc_length_mm": max_row["arc_length_mm"],
        # Reference (diaphragm)
        "ref_csa_mm2": ref_csa,
        "ref_diameter_mm": ref_diam,
        "ref_arc_length_mm": ref_arc,
        # Ratios
        "max_min_ratio": max_row["area_mm2"] / coa_csa if coa_csa > 0 else 0,
        # Landmarks
        "isthmus_arc_start": isthmus_start,
        "last_branch_arc": last_branch,
    }
    return metrics


def _get_branch_arc_lengths(df):
    """Get arc lengths where each branch diverges from MainAorta.

    Maps branch first-plane centroid to nearest MainAorta plane
    to get the divergence point on the MainAorta arc length.
    """
    main = df[df["region"] == "MainAorta"].sort_values("arc_length_mm")
    if len(main) == 0:
        return {}

    main_centroids = main[["centroid_x", "centroid_y", "centroid_z"]].values
    main_arcs = main["arc_length_mm"].values

    from scipy.spatial import cKDTree
    tree = cKDTree(main_centroids)

    branches = {}
    for region in sorted(df["region"].unique()):
        if region == "MainAorta":
            continue
        rdf = df[df["region"] == region].sort_values("arc_length_mm")
        if len(rdf) == 0:
            continue
        # First plane of branch — find nearest MainAorta plane
        first_centroid = rdf.iloc[0][["centroid_x", "centroid_y", "centroid_z"]].values.astype(float)
        dist, idx = tree.query(first_centroid)
        branches[region] = float(main_arcs[idx])

    return branches


# ---------------------------------------------------------------------------
# Option 1: Side-by-side (STL left, CSA right)
# ---------------------------------------------------------------------------

def render_side_by_side(df, stl_path, branches_dir, output_path, subject_id,
                        metrics=None):
    """STL sagittal view (left) + CSA plot with landmarks (right)."""
    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv
    from PIL import Image

    if metrics is None:
        metrics = compute_coa_metrics(df)

    main = df[df["region"] == "MainAorta"].sort_values("arc_length_mm")
    branch_arcs = _get_branch_arc_lengths(df)

    # Left panel: STL with planes
    pl = pv.Plotter(off_screen=True, window_size=(800, 1200))
    pl.set_background('white')
    mesh = pv.read(str(stl_path))
    pl.add_mesh(mesh, color='lightblue', opacity=0.15, smooth_shading=True)

    # Add MainAorta planes
    csa_dir = Path(output_path).parent
    main_planes = csa_dir / "MainAorta_Planes"
    if main_planes.exists():
        for sf in sorted(main_planes.glob("*.stl")):
            pl.add_mesh(pv.read(str(sf)), color='steelblue', opacity=0.4)

    # Branch planes
    for region in branch_arcs:
        bp = csa_dir / f"{region}_Planes"
        if bp.exists():
            for sf in sorted(bp.glob("*.stl")):
                pl.add_mesh(pv.read(str(sf)), color='gray', opacity=0.5)

    pl.view_yz()
    pl.camera.up = (0, 0, 1)
    pl.enable_parallel_projection()
    pl.reset_camera()
    pos = list(pl.camera_position)
    pos[0] = (-pos[0][0], pos[0][1], pos[0][2])
    pl.camera_position = pos
    pl.camera.zoom(0.85)
    stl_img = pl.screenshot(return_img=True)
    pl.close()

    # Right panel: CSA plot with annotations
    fig, ax = plt.subplots(figsize=(8, 12))
    ax.plot(main["arc_length_mm"], main["area_mm2"],
            "-", color="steelblue", linewidth=2)

    # Branch divergence lines
    colors_br = plt.cm.Set1(np.linspace(0, 1, max(len(branch_arcs), 1)))
    for i, (name, arc) in enumerate(branch_arcs.items()):
        ax.axvline(x=arc, color=colors_br[i], linewidth=1.5, linestyle="--", alpha=0.7)
        ax.text(arc, ax.get_ylim()[1] * 0.98, f" {name}",
                color=colors_br[i], fontsize=9, va='top', rotation=90)

    # Arch zone shading
    if len(branch_arcs) >= 2:
        arc_vals = sorted(branch_arcs.values())
        ax.axvspan(arc_vals[0], arc_vals[-1], alpha=0.08, color='orange', label='Arch zone')

    # Min/max markers
    if metrics:
        ax.plot(metrics["coa_arc_length_mm"], metrics["coa_csa_mm2"],
                "v", color="red", markersize=12, zorder=5,
                label=f"CoA (isthmus): {metrics['coa_csa_mm2']:.1f} mm²")
        ax.plot(metrics["max_arc_length_mm"], metrics["max_csa_mm2"],
                "^", color="green", markersize=10, zorder=5,
                label=f"Max: {metrics['max_csa_mm2']:.1f} mm²")

        # Metrics box
        box_text = (
            f"CI (isthmus/ref): {metrics['coarctation_index']:.2f}\n"
            f"Isthmus CSA: {metrics['coa_csa_mm2']:.1f} mm²\n"
            f"Ref CSA (diaphragm): {metrics['ref_csa_mm2']:.1f} mm²\n"
            f"Stenosis: {metrics['area_stenosis_pct']:.0f}%\n"
            f"Severity: {metrics['severity']}"
        )
        ax.text(0.98, 0.02, box_text, transform=ax.transAxes,
                fontsize=10, va='bottom', ha='right',
                bbox=dict(boxstyle='round', fc='lightyellow', ec=metrics['severity_color'],
                         linewidth=2, alpha=0.9))

    ax.set_xlabel("Arc Length (mm)", fontsize=11)
    ax.set_ylabel("CSA (mm²)", fontsize=11)
    ax.set_title(f"{subject_id} — CSA along Main Aorta", fontsize=12, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    fig.canvas.draw()
    plot_img = np.array(fig.canvas.buffer_rgba())[:, :, :3]
    plt.close()

    # Composite side-by-side
    stl_pil = Image.fromarray(stl_img[:, :, :3])
    plot_pil = Image.fromarray(plot_img)
    # Match heights
    h = max(stl_pil.height, plot_pil.height)
    stl_pil = stl_pil.resize((int(stl_pil.width * h / stl_pil.height), h), Image.LANCZOS)
    plot_pil = plot_pil.resize((int(plot_pil.width * h / plot_pil.height), h), Image.LANCZOS)

    total_w = stl_pil.width + plot_pil.width
    composite = Image.new('RGB', (total_w, h), (255, 255, 255))
    composite.paste(stl_pil, (0, 0))
    composite.paste(plot_pil, (stl_pil.width, 0))
    composite.save(str(output_path))
    log.info(f"  Side-by-side: {Path(output_path).name}")


# ---------------------------------------------------------------------------
# Option 2: Stacked (STL top, CSA bottom aligned)
# ---------------------------------------------------------------------------

def render_stacked_aligned(df, stl_path, output_path, subject_id, metrics=None):
    """STL on top, CSA below, vertical lines spanning both."""
    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv
    from PIL import Image

    if metrics is None:
        metrics = compute_coa_metrics(df)

    main = df[df["region"] == "MainAorta"].sort_values("arc_length_mm")
    branch_arcs = _get_branch_arc_lengths(df)
    csa_dir = Path(output_path).parent

    # Top: STL
    pl = pv.Plotter(off_screen=True, window_size=(1600, 800))
    pl.set_background('white')
    mesh = pv.read(str(stl_path))
    pl.add_mesh(mesh, color='lightblue', opacity=0.15, smooth_shading=True)
    main_planes = csa_dir / "MainAorta_Planes"
    if main_planes.exists():
        for sf in sorted(main_planes.glob("*.stl")):
            pl.add_mesh(pv.read(str(sf)), color='steelblue', opacity=0.4)
    for region in branch_arcs:
        bp = csa_dir / f"{region}_Planes"
        if bp.exists():
            for sf in sorted(bp.glob("*.stl")):
                pl.add_mesh(pv.read(str(sf)), color='gray', opacity=0.5)

    pl.view_yz(); pl.camera.up = (0, 0, 1)
    pl.enable_parallel_projection(); pl.reset_camera()
    pos = list(pl.camera_position)
    pos[0] = (-pos[0][0], pos[0][1], pos[0][2])
    pl.camera_position = pos
    pl.camera.zoom(0.85)
    stl_img = pl.screenshot(return_img=True)[:, :, :3]
    pl.close()

    # Bottom: CSA plot
    dpi = 100
    fig, ax = plt.subplots(figsize=(16, 6), dpi=dpi)
    ax.plot(main["arc_length_mm"], main["area_mm2"],
            "-", color="steelblue", linewidth=2)

    # Branch lines + arch zone
    colors_br = plt.cm.Set1(np.linspace(0, 1, max(len(branch_arcs), 1)))
    for i, (name, arc) in enumerate(branch_arcs.items()):
        ax.axvline(x=arc, color=colors_br[i], linewidth=1.5, linestyle="--", alpha=0.7)
        ax.text(arc, ax.get_ylim()[1] * 0.95, f" {name}", color=colors_br[i],
                fontsize=9, va='top')
    if len(branch_arcs) >= 2:
        arc_vals = sorted(branch_arcs.values())
        ax.axvspan(arc_vals[0], arc_vals[-1], alpha=0.08, color='orange')

    if metrics:
        ax.plot(metrics["coa_arc_length_mm"], metrics["coa_csa_mm2"],
                "v", color="red", markersize=12, zorder=5)
        ax.annotate(f"CI={metrics['coarctation_index']:.2f} ({metrics['severity']})",
                    xy=(0.5, 0.05), xycoords='axes fraction',
                    fontsize=12, fontweight='bold', color=metrics['severity_color'],
                    ha='center', bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.9))

    ax.set_xlabel("Arc Length (mm)", fontsize=11)
    ax.set_ylabel("CSA (mm²)", fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.canvas.draw()
    plot_img = np.array(fig.canvas.buffer_rgba())[:, :, :3]
    plt.close()

    # Stack
    stl_pil = Image.fromarray(stl_img).resize((1600, 800), Image.LANCZOS)
    plot_pil = Image.fromarray(plot_img).resize((1600, 600), Image.LANCZOS)
    composite = Image.new('RGB', (1600, 1400), (255, 255, 255))
    composite.paste(stl_pil, (0, 0))
    composite.paste(plot_pil, (0, 800))
    composite.save(str(output_path))
    log.info(f"  Stacked aligned: {Path(output_path).name}")


# ---------------------------------------------------------------------------
# Option 3: CSA color-mapped STL
# ---------------------------------------------------------------------------

def render_csa_colormap(df, stl_path, output_path, subject_id, metrics=None):
    """STL surface colored by local CSA value."""
    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv
    from scipy.spatial import cKDTree

    if metrics is None:
        metrics = compute_coa_metrics(df)

    main = df[df["region"] == "MainAorta"].sort_values("arc_length_mm")
    if len(main) == 0:
        return

    mesh = pv.read(str(stl_path))

    # Map CSA to surface: for each face center, find nearest plane centroid
    # and assign that plane's CSA value
    face_centers = mesh.cell_centers().points
    plane_centroids = main[["centroid_x", "centroid_y", "centroid_z"]].values
    plane_csa = main["area_mm2"].values

    tree = cKDTree(plane_centroids)
    dists, idxs = tree.query(face_centers)
    face_csa = plane_csa[idxs]

    # Faces too far from any plane get NaN (transparent)
    max_dist = 5.0  # mm
    face_csa[dists > max_dist] = np.nan

    mesh.cell_data["CSA"] = face_csa

    pl = pv.Plotter(off_screen=True, window_size=(1200, 1600))
    pl.set_background('white')

    # Color-mapped mesh
    pl.add_mesh(mesh, scalars="CSA", cmap="RdYlBu",
                clim=[main["area_mm2"].min() * 0.8, main["area_mm2"].max() * 1.1],
                nan_color='lightgray', nan_opacity=0.1,
                scalar_bar_args={
                    'title': 'CSA (mm²)',
                    'title_font_size': 14,
                    'label_font_size': 12,
                    'position_x': 0.05, 'position_y': 0.05,
                    'width': 0.3, 'height': 0.04,
                },
                smooth_shading=True)

    # Annotate min CSA location
    if metrics:
        # Find isthmus location by arc length
        arc_diff = (main["arc_length_mm"] - metrics["coa_arc_length_mm"]).abs()
        isthmus_idx = arc_diff.idxmin()
        coa_pt = np.array([
            float(main.loc[isthmus_idx, "centroid_x"]),
            float(main.loc[isthmus_idx, "centroid_y"]),
            float(main.loc[isthmus_idx, "centroid_z"]),
        ])
        pl.add_point_labels(
            coa_pt.reshape(1, 3),
            [f"Isthmus: {metrics['coa_csa_mm2']:.0f}mm² (CI={metrics['coarctation_index']:.2f})"],
            font_size=16, text_color='red', point_size=12,
            shape=None, always_visible=True)

    pl.view_yz(); pl.camera.up = (0, 0, 1)
    pl.enable_parallel_projection(); pl.reset_camera()
    pos = list(pl.camera_position)
    pos[0] = (-pos[0][0], pos[0][1], pos[0][2])
    pl.camera_position = pos
    pl.camera.zoom(0.85)
    pl.screenshot(str(output_path))
    pl.close()
    log.info(f"  CSA colormap: {Path(output_path).name}")


# ---------------------------------------------------------------------------
# Option 4: Centerline tube colored by CSA + threshold bar
# ---------------------------------------------------------------------------

def render_centerline_csa(df, stl_path, output_path, subject_id, metrics=None):
    """STL with thick centerline tube colored by CSA + clinical threshold bar."""
    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv
    from PIL import Image

    if metrics is None:
        metrics = compute_coa_metrics(df)

    main = df[df["region"] == "MainAorta"].sort_values("arc_length_mm")
    if len(main) == 0:
        return

    # Build centerline polyline from plane centroids
    centroids = main[["centroid_x", "centroid_y", "centroid_z"]].values.astype(float)
    csa_values = main["area_mm2"].values.astype(float)

    # Create tube along centerline, radius scaled by sqrt(CSA/pi) for visual
    points = pv.PolyData(centroids)
    n = len(centroids)
    lines = np.zeros(n + 1, dtype=int)
    lines[0] = n
    lines[1:] = np.arange(n)
    poly = pv.PolyData(centroids)
    poly.lines = lines
    poly.point_data["CSA"] = csa_values

    # Tube with constant radius (visual, not scaled)
    tube = poly.tube(radius=1.5, n_sides=20)
    tube.point_data["CSA"] = np.interp(
        np.arange(tube.n_points),
        np.linspace(0, tube.n_points - 1, n),
        csa_values,
    )

    mesh = pv.read(str(stl_path))

    # Compute CI-relative values (CSA / CSA_ref)
    ref_csa = metrics["ref_csa_mm2"] if metrics else csa_values.mean()
    ci_values = csa_values / ref_csa

    # Also create CI-colored tube
    tube_ci = poly.tube(radius=1.5, n_sides=20)
    tube_ci.point_data["CI"] = np.interp(
        np.arange(tube_ci.n_points),
        np.linspace(0, tube_ci.n_points - 1, n),
        ci_values,
    )

    # Helper to set up camera
    def _setup_cam(pl):
        pl.view_yz()
        pl.camera.up = (0, 0, 1)
        pl.enable_parallel_projection()
        pl.reset_camera()
        pos = list(pl.camera_position)
        pos[0] = (-pos[0][0], pos[0][1], pos[0][2])
        pl.camera_position = pos
        pl.camera.zoom(0.85)

    # Helper to add labels
    def _add_labels(pl):
        if not metrics:
            return
        offset = np.array([15.0, 0.0, 0.0])

        arc_diff = (main["arc_length_mm"] - metrics["coa_arc_length_mm"]).abs()
        isthmus_idx = arc_diff.idxmin()
        coa_pt = np.array([float(main.loc[isthmus_idx, c])
                           for c in ["centroid_x", "centroid_y", "centroid_z"]])
        pl.add_point_labels((coa_pt + offset).reshape(1, 3),
            [f"Isthmus: {metrics['coa_csa_mm2']:.0f} mm²"],
            font_size=18, text_color='red', point_size=0, shape=None, always_visible=True)
        pl.add_mesh(pv.Line(coa_pt, coa_pt + offset * 0.8), color='red', line_width=2)
        pl.add_mesh(pv.Sphere(radius=0.8, center=coa_pt), color='red')

        arc_diff_ref = (main["arc_length_mm"] - metrics["ref_arc_length_mm"]).abs()
        ref_idx = arc_diff_ref.idxmin()
        ref_pt = np.array([float(main.loc[ref_idx, c])
                           for c in ["centroid_x", "centroid_y", "centroid_z"]])
        pl.add_point_labels((ref_pt + offset).reshape(1, 3),
            [f"Ref: {metrics['ref_csa_mm2']:.0f} mm²"],
            font_size=16, text_color='darkblue', point_size=0, shape=None, always_visible=True)
        pl.add_mesh(pv.Line(ref_pt, ref_pt + offset * 0.8), color='darkblue', line_width=2)
        pl.add_mesh(pv.Sphere(radius=0.8, center=ref_pt), color='darkblue')

    # ── Left panel: Absolute CSA ──
    pl1 = pv.Plotter(off_screen=True, window_size=(900, 1200))
    pl1.set_background('white')
    pl1.add_mesh(mesh, color='lightblue', opacity=0.1, smooth_shading=True)
    pl1.add_mesh(tube, scalars="CSA", cmap="RdYlBu",
                clim=[csa_values.min() * 0.9, csa_values.max() * 1.05],
                scalar_bar_args={
                    'title': 'CSA (mm²)', 'title_font_size': 24,
                    'label_font_size': 18, 'bold': True,
                    'position_x': 0.08, 'position_y': 0.15,
                    'width': 0.08, 'height': 0.55, 'vertical': True,
                })
    _add_labels(pl1)
    _setup_cam(pl1)
    img_left = pl1.screenshot(return_img=True)[:, :, :3]
    pl1.close()

    # ── Right panel: CI relative ──
    pl2 = pv.Plotter(off_screen=True, window_size=(900, 1200))
    pl2.set_background('white')
    pl2.add_mesh(mesh, color='lightblue', opacity=0.1, smooth_shading=True)
    pl2.add_mesh(tube_ci, scalars="CI", cmap="RdYlGn",
                clim=[0, 2.0],
                scalar_bar_args={
                    'title': 'CI (CSA/CSA_ref)', 'title_font_size': 24,
                    'label_font_size': 18, 'bold': True,
                    'position_x': 0.10, 'position_y': 0.15,
                    'width': 0.08, 'height': 0.55, 'vertical': True,
                })
    # Add CI-specific labels
    if metrics:
        offset = np.array([15.0, 0.0, 0.0])
        arc_diff = (main["arc_length_mm"] - metrics["coa_arc_length_mm"]).abs()
        isthmus_idx = arc_diff.idxmin()
        coa_pt = np.array([float(main.loc[isthmus_idx, c])
                           for c in ["centroid_x", "centroid_y", "centroid_z"]])
        pl2.add_point_labels((coa_pt + offset).reshape(1, 3),
            [f"CI = {metrics['coarctation_index']:.2f}"],
            font_size=20, text_color='red', point_size=0, shape=None, always_visible=True)
        pl2.add_mesh(pv.Line(coa_pt, coa_pt + offset * 0.8), color='red', line_width=2)
        pl2.add_mesh(pv.Sphere(radius=0.8, center=coa_pt), color='red')

        arc_diff_ref = (main["arc_length_mm"] - metrics["ref_arc_length_mm"]).abs()
        ref_idx = arc_diff_ref.idxmin()
        ref_pt = np.array([float(main.loc[ref_idx, c])
                           for c in ["centroid_x", "centroid_y", "centroid_z"]])
        pl2.add_point_labels((ref_pt + offset).reshape(1, 3),
            [f"CI = 1.0 (ref)"],
            font_size=16, text_color='darkgreen', point_size=0, shape=None, always_visible=True)
        pl2.add_mesh(pv.Line(ref_pt, ref_pt + offset * 0.8), color='darkgreen', line_width=2)
        pl2.add_mesh(pv.Sphere(radius=0.8, center=ref_pt), color='darkgreen')

    _setup_cam(pl2)
    img_right = pl2.screenshot(return_img=True)[:, :, :3]
    pl2.close()

    # Composite side by side
    from PIL import Image as PILImage
    h = max(img_left.shape[0], img_right.shape[0])
    pw = 900  # panel width
    left_pil = PILImage.fromarray(img_left).resize((pw, h), PILImage.LANCZOS)
    right_pil = PILImage.fromarray(img_right).resize((pw, h), PILImage.LANCZOS)
    stl_composite = PILImage.new('RGB', (pw * 2, h), (255, 255, 255))
    stl_composite.paste(left_pil, (0, 0))
    stl_composite.paste(right_pil, (pw, 0))
    stl_img = np.array(stl_composite)

    # ── Bottom: Clinical threshold bar ──
    dpi = 100
    fig, axes = plt.subplots(2, 1, figsize=(12, 3), dpi=dpi,
                              gridspec_kw={'height_ratios': [2, 1]})

    # CI gauge bar
    ax = axes[0]
    # Draw colored zones
    ax.barh(0, 0.5, left=0, height=0.6, color='red', alpha=0.7)
    ax.barh(0, 0.2, left=0.5, height=0.6, color='orange', alpha=0.7)
    ax.barh(0, 0.3, left=0.7, height=0.6, color='green', alpha=0.7)

    # Zone labels
    ax.text(0.25, -0.15, 'Severe\n(<0.5)', ha='center', va='top', fontsize=10, color='red', fontweight='bold')
    ax.text(0.6, -0.15, 'Moderate\n(0.5-0.7)', ha='center', va='top', fontsize=9, color='darkorange', fontweight='bold')
    ax.text(0.85, -0.15, 'Mild\n(>0.7)', ha='center', va='top', fontsize=10, color='green', fontweight='bold')

    # Patient marker
    if metrics:
        ci_val = metrics['coarctation_index']
        if ci_val > 1.0:
            # CI beyond normal — show arrow at right edge with annotation
            ax.annotate(f"CI = {ci_val:.2f}\n(No stenosis)",
                       xy=(0.98, 0), xytext=(0.85, 0.55),
                       fontsize=13, fontweight='bold', color='green',
                       ha='center', va='bottom',
                       arrowprops=dict(arrowstyle='->', color='green', lw=2))
        else:
            ci_plot = max(ci_val, 0.02)
            ax.plot(ci_plot, 0, 'v', color='black', markersize=20, zorder=5)
            ax.text(ci_plot, 0.45, f"CI = {ci_val:.2f}",
                    ha='center', va='bottom', fontsize=14, fontweight='bold')

    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.5, 0.8)
    ax.set_title(f"{subject_id} — Coarctation Index", fontsize=14, fontweight='bold')
    ax.set_xlabel("Coarctation Index (D_isthmus / D_reference)", fontsize=11)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Summary text
    ax2 = axes[1]
    ax2.axis('off')
    if metrics:
        stenosis_str = (f"{metrics['area_stenosis_pct']:.0f}%"
                        if metrics['area_stenosis_pct'] > 0 else "None")
        summary = (
            f"Isthmus CSA: {metrics['coa_csa_mm2']:.1f} mm²  |  "
            f"Reference CSA: {metrics['ref_csa_mm2']:.1f} mm²  |  "
            f"Stenosis: {stenosis_str}  |  "
            f"Max/Min: {metrics['max_min_ratio']:.1f}x  |  "
            f"Severity: {metrics['severity']}"
        )
        ax2.text(0.5, 0.5, summary, ha='center', va='center', fontsize=12,
                fontweight='bold', color=metrics['severity_color'],
                bbox=dict(boxstyle='round,pad=0.5', fc='lightyellow', ec=metrics['severity_color'],
                         linewidth=2))

    plt.tight_layout()
    fig.canvas.draw()
    bar_img = np.array(fig.canvas.buffer_rgba())[:, :, :3]
    plt.close()

    # Composite: STL top, threshold bar bottom
    stl_pil = Image.fromarray(stl_img)
    bar_pil = Image.fromarray(bar_img)
    w = max(stl_pil.width, bar_pil.width)
    stl_pil = stl_pil.resize((w, int(stl_pil.height * w / stl_pil.width)), Image.LANCZOS)
    bar_pil = bar_pil.resize((w, int(bar_pil.height * w / bar_pil.width)), Image.LANCZOS)

    composite = Image.new('RGB', (w, stl_pil.height + bar_pil.height), (255, 255, 255))
    composite.paste(stl_pil, (0, 0))
    composite.paste(bar_pil, (0, stl_pil.height))
    composite.save(str(output_path))
    log.info(f"  Centerline CSA: {Path(output_path).name}")


# ---------------------------------------------------------------------------
# Option 5: Spinning GIF of centerline tube
# ---------------------------------------------------------------------------

def render_centerline_gif(df, stl_path, output_path, subject_id, metrics=None,
                          n_frames=36, duration_ms=100):
    """Spinning GIF: STL + colored centerline tube rotating 360 degrees."""
    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv
    from PIL import Image

    if metrics is None:
        metrics = compute_coa_metrics(df)

    main = df[df["region"] == "MainAorta"].sort_values("arc_length_mm")
    if len(main) == 0:
        return

    centroids = main[["centroid_x", "centroid_y", "centroid_z"]].values.astype(float)
    csa_values = main["area_mm2"].values.astype(float)

    # Build tube
    n = len(centroids)
    poly = pv.PolyData(centroids)
    lines = np.zeros(n + 1, dtype=int)
    lines[0] = n
    lines[1:] = np.arange(n)
    poly.lines = lines
    poly.point_data["CSA"] = csa_values

    tube = poly.tube(radius=1.5, n_sides=20)
    tube.point_data["CSA"] = np.interp(
        np.arange(tube.n_points),
        np.linspace(0, tube.n_points - 1, n),
        csa_values,
    )

    mesh = pv.read(str(stl_path))

    frames = []
    for i in range(n_frames):
        azimuth = i * (360 / n_frames)

        pl = pv.Plotter(off_screen=True, window_size=(800, 1000))
        pl.set_background('white')
        pl.add_mesh(mesh, color='lightblue', opacity=0.1, smooth_shading=True)
        pl.add_mesh(tube, scalars="CSA", cmap="RdYlBu",
                    clim=[csa_values.min() * 0.9, csa_values.max() * 1.05],
                    show_scalar_bar=False)

        # Labels (only on front-facing frames to avoid clutter)
        if metrics and (i < 5 or i > n_frames - 5 or abs(i - n_frames // 2) < 3):
            arc_diff = (main["arc_length_mm"] - metrics["coa_arc_length_mm"]).abs()
            isthmus_idx = arc_diff.idxmin()
            coa_pt = np.array([
                float(main.loc[isthmus_idx, "centroid_x"]),
                float(main.loc[isthmus_idx, "centroid_y"]),
                float(main.loc[isthmus_idx, "centroid_z"]),
            ])
            pl.add_point_labels(
                coa_pt.reshape(1, 3),
                [f"Isthmus: {metrics['coa_csa_mm2']:.0f} mm²"],
                font_size=16, text_color='red',
                point_size=8, shape=None, always_visible=True)

        # Camera: rotate around Y axis (vertical)
        pl.camera.up = (0, 0, 1)
        center = mesh.center
        diag = mesh.length
        rad = np.radians(azimuth)
        cam_pos = (
            center[0] + diag * 0.9 * np.cos(rad),
            center[1] + diag * 0.9 * np.sin(rad),
            center[2],
        )
        pl.camera_position = [cam_pos, center, (0, 0, 1)]
        pl.enable_parallel_projection()

        img = pl.screenshot(return_img=True)[:, :, :3]
        pl.close()
        frames.append(Image.fromarray(img))

        if (i + 1) % 10 == 0:
            log.info(f"    GIF frame {i+1}/{n_frames}")

    # Save GIF
    frames[0].save(
        str(output_path),
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
    )
    log.info(f"  Spinning GIF: {Path(output_path).name} ({n_frames} frames)")


# ---------------------------------------------------------------------------
# Render all views
# ---------------------------------------------------------------------------

def render_all_views(df, stl_path, branches_dir, output_dir, subject_id):
    """Generate all 3 CoA visualization options + metrics."""
    output_dir = Path(output_dir)
    metrics = compute_coa_metrics(df)

    if metrics:
        log.info(f"\n  CoA Metrics:")
        log.info(f"    Isthmus CSA: {metrics['coa_csa_mm2']:.1f} mm² "
                 f"(at {metrics['coa_arc_length_mm']:.1f}mm, after Branch_3)")
        log.info(f"    Reference CSA: {metrics['ref_csa_mm2']:.1f} mm² "
                 f"(descending at diaphragm)")
        log.info(f"    Coarctation Index (CI): {metrics['coarctation_index']:.2f}")
        log.info(f"    Area Stenosis: {metrics['area_stenosis_pct']:.0f}%")
        log.info(f"    Severity: {metrics['severity']}")
        log.info(f"    Max/Min Ratio: {metrics['max_min_ratio']:.1f}x")

    # Centerline tube + threshold bar (doctor-friendly)
    render_centerline_csa(df, stl_path,
                          output_dir / f"{subject_id}_coa_centerline.png",
                          subject_id, metrics)


    return metrics
