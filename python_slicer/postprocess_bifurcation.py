#!/usr/bin/env python3
"""
Post-processing for bifurcation CSA pipeline.

Runs per-region analysis (heatmaps, compliance, dynamics) and
whole-airway visualization (video, interactive HTML).

Works dynamically with any number of frames.

Usage:
    python python_slicer/postprocess_bifurcation.py ENT001/csa_bifurcation/
"""

import sys
import argparse
import logging
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
import trimesh

sys.path.insert(0, str(Path(__file__).parent))

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-region CSA plots
# ---------------------------------------------------------------------------

def plot_csa_by_index_per_region(df, output_dir, subject_id):
    """CSA by plane index, one plot per region. One line per frame."""
    regions = [r for r in df["region"].unique() if r != "NasalCombined"]

    for region in regions:
        rdf = df[df["region"] == region].copy()
        frames = sorted(rdf["frame_name"].unique()) if "frame_name" in rdf.columns else [None]
        n_frames = len(frames)

        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        cmap = plt.cm.viridis(np.linspace(0, 1, max(n_frames, 1)))

        for fi, frame in enumerate(frames):
            fdf = rdf[rdf["frame_name"] == frame] if frame else rdf
            fdf = fdf.sort_values("plane_index")
            c = cmap[fi]
            label = frame if fi == 0 or fi == n_frames - 1 else None
            axes[0].plot(fdf["plane_index"], fdf["area_mm2"],
                        "-", color=c, linewidth=0.8, alpha=0.6, label=label)
            axes[1].plot(fdf["plane_index"], fdf["hydraulic_diameter_mm"],
                        "-", color=c, linewidth=0.8, alpha=0.6)

        axes[0].set_ylabel("CSA (mm²)")
        axes[0].set_title(f"{subject_id} — {region} — CSA by Plane Index")
        axes[0].legend(fontsize=7)
        axes[0].grid(True, alpha=0.3)
        axes[1].set_ylabel("Hydraulic Diameter (mm)")
        axes[1].set_xlabel("Plane Index")
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = output_dir / f"{subject_id}_{region}_CSA_by_index.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info(f"  {out_path.name}")


def plot_csa_band_per_region(df, output_dir, subject_id):
    """CSA band plot: mean line + min/max shaded region through breathing cycle."""
    regions = [r for r in df["region"].unique() if r != "NasalCombined"]
    if "frame_name" not in df.columns:
        return

    for region in regions:
        rdf = df[df["region"] == region]
        planes = sorted(rdf["plane_index"].unique())

        # Compute per-plane stats across all frames
        stats = rdf.groupby("plane_index").agg(
            area_mean=("area_mm2", "mean"),
            area_min=("area_mm2", "min"),
            area_max=("area_mm2", "max"),
            area_std=("area_mm2", "std"),
            arc=("arc_length_mm", "mean"),
            hyd_mean=("hydraulic_diameter_mm", "mean"),
            hyd_min=("hydraulic_diameter_mm", "min"),
            hyd_max=("hydraulic_diameter_mm", "max"),
        ).reindex(planes)

        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        # CSA band
        axes[0].fill_between(stats.index, stats["area_min"], stats["area_max"],
                             alpha=0.3, color="steelblue", label="Min-Max range")
        axes[0].plot(stats.index, stats["area_mean"], "-", color="steelblue",
                     linewidth=1.5, label="Mean")
        axes[0].set_ylabel("CSA (mm²)")
        axes[0].set_title(f"{subject_id} — {region} — CSA Band (min/mean/max)")
        axes[0].legend(fontsize=9)
        axes[0].grid(True, alpha=0.3)

        # Hydraulic diameter band
        axes[1].fill_between(stats.index, stats["hyd_min"], stats["hyd_max"],
                             alpha=0.3, color="coral", label="Min-Max range")
        axes[1].plot(stats.index, stats["hyd_mean"], "-", color="coral",
                     linewidth=1.5, label="Mean")
        axes[1].set_ylabel("Hydraulic Diameter (mm)")
        axes[1].set_xlabel("Plane Index")
        axes[1].legend(fontsize=9)
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = output_dir / f"{subject_id}_{region}_CSA_band.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info(f"  {out_path.name}")


def plot_csa_heatmap_per_region(df, output_dir, subject_id):
    """CSA heatmap: plane_index × breathing phase, one per region."""
    regions = [r for r in df["region"].unique() if r != "NasalCombined"]

    for region in regions:
        rdf = df[df["region"] == region].copy()
        if "frame_index" not in rdf.columns:
            continue

        planes = sorted(rdf["plane_index"].unique())
        frames = sorted(rdf["frame_index"].unique())

        if len(frames) < 2:
            continue

        # Build matrix
        matrix = np.full((len(planes), len(frames)), np.nan)
        plane_map = {p: i for i, p in enumerate(planes)}
        frame_map = {f: i for i, f in enumerate(frames)}

        for _, row in rdf.iterrows():
            pi = plane_map.get(row["plane_index"])
            fi = frame_map.get(row["frame_index"])
            if pi is not None and fi is not None:
                matrix[pi, fi] = row["area_mm2"]

        fig, ax = plt.subplots(figsize=(12, 8))
        im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd",
                       interpolation="nearest", origin="lower")
        ax.set_xlabel("Frame Index")
        ax.set_ylabel("Plane Index")
        ax.set_title(f"{subject_id} — {region} — CSA Heatmap")
        plt.colorbar(im, ax=ax, label="CSA (mm²)")

        plt.tight_layout()
        out_path = output_dir / f"{subject_id}_{region}_CSA_heatmap.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info(f"  {out_path.name}")


# ---------------------------------------------------------------------------
# Compliance analysis
# ---------------------------------------------------------------------------

def compute_compliance_per_region(df, output_dir, subject_id):
    """Compute compliance % per plane per region.

    Compliance = (max_CSA - min_CSA) / mean_CSA × 100%
    >30% = high collapse risk, 15-30% = moderate, <15% = stable
    """
    regions = [r for r in df["region"].unique() if r != "NasalCombined"]
    if "frame_name" not in df.columns:
        return

    all_compliance = []

    for region in regions:
        rdf = df[df["region"] == region]
        planes = sorted(rdf["plane_index"].unique())

        for pi in planes:
            pdata = rdf[rdf["plane_index"] == pi]
            if len(pdata) < 2:
                continue
            area_min = pdata["area_mm2"].min()
            area_max = pdata["area_mm2"].max()
            area_mean = pdata["area_mm2"].mean()
            area_std = pdata["area_mm2"].std()

            compliance = (area_max - area_min) / area_mean * 100 if area_mean > 0 else 0

            all_compliance.append({
                "region": region,
                "plane_index": pi,
                "arc_length_mm": pdata["arc_length_mm"].mean(),
                "area_mean_mm2": area_mean,
                "area_std_mm2": area_std,
                "area_min_mm2": area_min,
                "area_max_mm2": area_max,
                "area_range_mm2": area_max - area_min,
                "area_cv_percent": (area_std / area_mean * 100) if area_mean > 0 else 0,
                "compliance_percent": compliance,
            })

    if not all_compliance:
        return

    comp_df = pd.DataFrame(all_compliance)
    csv_path = output_dir / f"{subject_id}_compliance.csv"
    comp_df.to_csv(csv_path, index=False)
    log.info(f"  Compliance CSV: {csv_path.name}")

    # Plot compliance per region
    for region in regions:
        rcomp = comp_df[comp_df["region"] == region]
        if len(rcomp) == 0:
            continue

        fig, ax = plt.subplots(figsize=(12, 5))
        colors = []
        for _, row in rcomp.iterrows():
            c = row["compliance_percent"]
            if c > 30:
                colors.append("red")
            elif c > 15:
                colors.append("orange")
            else:
                colors.append("green")

        ax.bar(rcomp["plane_index"], rcomp["compliance_percent"], color=colors, alpha=0.8)
        ax.axhline(y=30, color="red", linestyle="--", alpha=0.5, label=">30% High risk")
        ax.axhline(y=15, color="orange", linestyle="--", alpha=0.5, label=">15% Moderate")
        ax.set_xlabel("Plane Index")
        ax.set_ylabel("Compliance (%)")
        ax.set_title(f"{subject_id} — {region} — Compliance (CSA variation)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = output_dir / f"{subject_id}_{region}_compliance.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info(f"  {out_path.name}")

    return comp_df


# ---------------------------------------------------------------------------
# Enhanced metrics (circularity, eccentricity, shape)
# ---------------------------------------------------------------------------

def compute_enhanced_metrics(df, output_dir, subject_id):
    """Add circularity and eccentricity to the dataframe."""
    df = df.copy()

    # Circularity = 4π × Area / Perimeter²
    df["circularity"] = (4 * np.pi * df["area_mm2"]) / (df["perimeter_mm"] ** 2)
    df.loc[df["perimeter_mm"] <= 0, "circularity"] = 0

    # Eccentricity from major/minor if available
    if "major_axis_mm" in df.columns and "minor_axis_mm" in df.columns:
        ratio = df["minor_axis_mm"] / df["major_axis_mm"].replace(0, np.nan)
        df["eccentricity"] = np.sqrt(1 - ratio**2)
    else:
        df["eccentricity"] = np.nan

    csv_path = output_dir / f"{subject_id}_enhanced_metrics.csv"
    df.to_csv(csv_path, index=False)
    log.info(f"  Enhanced metrics: {csv_path.name}")

    # Plot circularity per region
    regions = [r for r in df["region"].unique() if r != "NasalCombined"]
    for region in regions:
        rdf = df[df["region"] == region]
        if "frame_name" not in rdf.columns:
            continue

        # Average circularity per plane across frames
        plane_stats = rdf.groupby("plane_index").agg(
            circ_mean=("circularity", "mean"),
            circ_std=("circularity", "std"),
            arc=("arc_length_mm", "mean"),
        ).reset_index()

        fig, ax = plt.subplots(figsize=(12, 5))
        ax.fill_between(plane_stats["plane_index"],
                        plane_stats["circ_mean"] - plane_stats["circ_std"],
                        plane_stats["circ_mean"] + plane_stats["circ_std"],
                        alpha=0.3, color="steelblue")
        ax.plot(plane_stats["plane_index"], plane_stats["circ_mean"],
                "-", color="steelblue", linewidth=1.5)
        ax.set_xlabel("Plane Index")
        ax.set_ylabel("Circularity")
        ax.set_title(f"{subject_id} — {region} — Circularity (1.0=circle)")
        ax.set_ylim(0, 1.1)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = output_dir / f"{subject_id}_{region}_circularity.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info(f"  {out_path.name}")

    return df


# ---------------------------------------------------------------------------
# Centroid movement tracking
# ---------------------------------------------------------------------------

def compute_centroid_movement(df, output_dir, subject_id):
    """Track centroid displacement per plane across frames."""
    if "frame_name" not in df.columns:
        return

    regions = [r for r in df["region"].unique() if r != "NasalCombined"]
    all_movement = []

    for region in regions:
        rdf = df[df["region"] == region]
        planes = sorted(rdf["plane_index"].unique())

        for pi in planes:
            pdata = rdf[rdf["plane_index"] == pi].sort_values("frame_index")
            if len(pdata) < 2:
                continue

            centroids = pdata[["centroid_x", "centroid_y", "centroid_z"]].values
            mean_pos = centroids.mean(axis=0)
            displacements = np.linalg.norm(centroids - mean_pos, axis=1)

            # Path length (sum of frame-to-frame displacements)
            frame_diffs = np.linalg.norm(np.diff(centroids, axis=0), axis=1)

            all_movement.append({
                "region": region,
                "plane_index": pi,
                "arc_length_mm": pdata["arc_length_mm"].mean(),
                "max_displacement_mm": displacements.max(),
                "mean_displacement_mm": displacements.mean(),
                "std_displacement_mm": displacements.std(),
                "total_path_length_mm": frame_diffs.sum(),
            })

    if not all_movement:
        return

    mov_df = pd.DataFrame(all_movement)
    csv_path = output_dir / f"{subject_id}_centroid_movement.csv"
    mov_df.to_csv(csv_path, index=False)
    log.info(f"  Centroid movement: {csv_path.name}")

    # Plot per region
    for region in regions:
        rmov = mov_df[mov_df["region"] == region]
        if len(rmov) == 0:
            continue

        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        axes[0].bar(rmov["plane_index"], rmov["max_displacement_mm"],
                    color="steelblue", alpha=0.7)
        axes[0].set_ylabel("Max Displacement (mm)")
        axes[0].set_title(f"{subject_id} — {region} — Centroid Movement")
        axes[0].grid(True, alpha=0.3)

        axes[1].bar(rmov["plane_index"], rmov["mean_displacement_mm"],
                    color="coral", alpha=0.7, label="Mean")
        axes[1].errorbar(rmov["plane_index"], rmov["mean_displacement_mm"],
                        yerr=rmov["std_displacement_mm"], fmt="none",
                        color="gray", alpha=0.5)
        axes[1].set_ylabel("Mean Displacement (mm)")
        axes[1].set_xlabel("Plane Index")
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = output_dir / f"{subject_id}_{region}_centroid_movement.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info(f"  {out_path.name}")


# ---------------------------------------------------------------------------
# Breathing cycle statistics
# ---------------------------------------------------------------------------

def compute_breathing_stats(df, output_dir, subject_id):
    """Per-frame statistics across all planes."""
    if "frame_name" not in df.columns:
        return

    regions = [r for r in df["region"].unique() if r != "NasalCombined"]
    all_stats = []

    for region in regions:
        rdf = df[df["region"] == region]
        for frame, fdf in rdf.groupby("frame_name"):
            fi = fdf["frame_index"].iloc[0] if "frame_index" in fdf.columns else 0
            all_stats.append({
                "region": region,
                "frame_name": frame,
                "frame_index": fi,
                "n_planes": len(fdf),
                "area_min_mm2": fdf["area_mm2"].min(),
                "area_max_mm2": fdf["area_mm2"].max(),
                "area_mean_mm2": fdf["area_mm2"].mean(),
                "area_std_mm2": fdf["area_mm2"].std(),
                "area_median_mm2": fdf["area_mm2"].median(),
                "hyd_diam_mean_mm": fdf["hydraulic_diameter_mm"].mean(),
                "hyd_diam_min_mm": fdf["hydraulic_diameter_mm"].min(),
                "hyd_diam_max_mm": fdf["hydraulic_diameter_mm"].max(),
            })

    if all_stats:
        stats_df = pd.DataFrame(all_stats)
        csv_path = output_dir / f"{subject_id}_breathing_stats.csv"
        stats_df.to_csv(csv_path, index=False)
        log.info(f"  Breathing stats: {csv_path.name}")


# ---------------------------------------------------------------------------
# Whole-airway combined plane STLs (for video)
# ---------------------------------------------------------------------------

def generate_combined_plane_stls(frames_dir, output_dir, subject_id):
    """Combine regional plane STLs into one per frame (for video/visualization)."""
    frames_dir = Path(frames_dir)
    combined_dir = output_dir / "combined_planes"
    combined_dir.mkdir(parents=True, exist_ok=True)

    frame_dirs = sorted([d for d in frames_dir.iterdir() if d.is_dir()])

    for frame_dir in frame_dirs:
        frame_name = frame_dir.name
        all_stls = []

        for region_folder in ["DescendingAirway_Planes", "LeftNose_Planes", "RightNose_Planes"]:
            region_dir = frame_dir / region_folder
            if region_dir.exists():
                stl_files = sorted(region_dir.glob("*.stl"))
                for stl_file in stl_files:
                    try:
                        mesh = trimesh.load_mesh(str(stl_file))
                        all_stls.append(mesh)
                    except Exception:
                        continue

        if all_stls:
            combined = trimesh.util.concatenate(all_stls)
            out_path = combined_dir / f"{frame_name}-Planes-All.stl"
            combined.export(str(out_path))

    n_combined = len(list(combined_dir.glob("*.stl")))
    log.info(f"  Combined plane STLs: {n_combined} frames in {combined_dir.name}/")
    return combined_dir


# ---------------------------------------------------------------------------
# Breathing cycle video
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 4-panel breathing cycle video
# ---------------------------------------------------------------------------

def _render_pyvista_panel(airway_region_paths=None, airway_mesh_path=None,
                          plane_mesh_path=None, region_plane_paths=None,
                          window_size=(800, 800), clip_bounds=None,
                          view="sagittal", zoom_factor=0.85):
    """Render a single 3D panel with PyVista. Returns image as numpy array.

    airway_region_paths: dict {"DescendingAirway": path, "LeftNose": path, "RightNose": path}
        Colors airway surface by region (gray/blue/red).
    airway_mesh_path: single STL (fallback, rendered as lightblue).
    plane_mesh_path: single combined plane STL.
    region_plane_paths: dict of plane STL paths per region.
    """
    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv

    pl = pv.Plotter(off_screen=True, window_size=window_size)
    pl.set_background('white')

    region_color_map = {
        "DescendingAirway": "gray", "LeftNose": "blue",
        "RightNose": "red", "Mouth": "orange",
    }

    # Airway surface — colored by region or single color
    if airway_region_paths:
        for region, path in airway_region_paths.items():
            if path and Path(path).exists():
                aw = pv.read(str(path))
                if clip_bounds:
                    aw = aw.clip_box(clip_bounds, invert=False)
                c = region_color_map.get(region, "lightblue")
                pl.add_mesh(aw, color=c, opacity=0.15,
                            smooth_shading=True, show_edges=False)
    elif airway_mesh_path:
        aw = pv.read(str(airway_mesh_path))
        if clip_bounds:
            aw = aw.clip_box(clip_bounds, invert=False)
        pl.add_mesh(aw, color='lightblue', opacity=0.15,
                    smooth_shading=True, show_edges=False)

    # Plane meshes
    if region_plane_paths:
        for region, stl_paths in region_plane_paths.items():
            c = region_color_map.get(region, "gray")
            for sp in stl_paths:
                pm = pv.read(str(sp))
                if clip_bounds:
                    pm = pm.clip_box(clip_bounds, invert=False)
                pl.add_mesh(pm, color=c, opacity=0.6, show_edges=False)
    elif plane_mesh_path:
        pm = pv.read(str(plane_mesh_path))
        if clip_bounds:
            pm = pm.clip_box(clip_bounds, invert=False)
        centroids = pm.cell_centers().points
        if len(centroids) > 0:
            y_vals = centroids[:, 1]
            pm.cell_data['position'] = y_vals
            pl.add_mesh(pm, scalars='position', cmap='rainbow',
                        opacity=0.7, show_edges=False,
                        show_scalar_bar=False)

    if view == "sagittal":
        # Look along -X axis (from left side) so nose is on the left
        pl.view_yz()
        pl.camera.up = (0, 0, 1)
        pl.enable_parallel_projection()
        pl.reset_camera()
        # Flip camera to other side of X axis
        pos = list(pl.camera_position)
        pos[0] = (-pos[0][0], pos[0][1], pos[0][2])
        pl.camera_position = pos
        pl.camera.zoom(zoom_factor)
    elif view == "perspective_45":
        pl.reset_camera()
        center = pl.center
        diag = pl.length
        pl.camera_position = [
            (center[0] + diag * 0.5, center[1] - diag * 0.3, center[2] + diag * 0.4),
            center,
            (0, 0, 1),
        ]
        pl.camera.zoom(zoom_factor)

    img = pl.screenshot(return_img=True)
    pl.close()
    return img


def generate_4panel_video(combined_dir, motion_stl_dir, df, output_dir,
                          subject_id, fps=5, flow_profile_path=None):
    """Create professional 4-panel MP4 using PyVista rendering:
      Top-left:     Full airway STL (45° perspective)
      Top-right:    STL + CSA planes overlay
      Bottom-left:  CSA band plot with moving frame line
      Bottom-right: Zoomed nasopharynx/larynx with STL + planes
    """
    import subprocess
    from PIL import Image
    import io

    combined_dir = Path(combined_dir)
    motion_stl_dir = Path(motion_stl_dir)

    plane_files = sorted(combined_dir.glob("*-Planes-All.stl"))
    if not plane_files:
        log.warning("  No combined plane STLs, skipping 4-panel video")
        return

    frames_4p_dir = output_dir / "video_frames_4panel"
    frames_4p_dir.mkdir(parents=True, exist_ok=True)

    # Pre-compute CSA band per region
    region_stats = {}
    region_colors = {"DescendingAirway": "gray", "LeftNose": "blue", "RightNose": "red"}
    for region in ["DescendingAirway", "LeftNose", "RightNose"]:
        rdf = df[df["region"] == region]
        if len(rdf) == 0:
            continue
        stats = rdf.groupby("plane_index").agg(
            area_mean=("area_mm2", "mean"),
            area_min=("area_mm2", "min"),
            area_max=("area_mm2", "max"),
            arc_mean=("arc_length_mm", "mean"),
        ).reset_index().sort_values("arc_mean")
        region_stats[region] = stats

    frame_names = sorted(df["frame_name"].unique()) if "frame_name" in df.columns else []

    # Renormalize arc length for display: nose=0, trachea=max (nose-to-trachea left-to-right)
    arc_max = df["arc_length_mm"].max()
    for region in region_stats:
        region_stats[region]["arc_display"] = arc_max - region_stats[region]["arc_mean"]

    # Load flow profile if available
    flow_df = None
    if flow_profile_path is None:
        # Auto-detect in subject dir
        subject_dir = output_dir.parent
        for pattern in ["*FlowProfile.csv", "*flowprofile.csv", "*flow_profile.csv"]:
            matches = list(subject_dir.glob(pattern))
            if matches:
                flow_profile_path = matches[0]
                break
    if flow_profile_path and Path(flow_profile_path).exists():
        flow_df = pd.read_csv(flow_profile_path)
        flow_df.columns = ["time_s", "flow_kgs"]
        flow_df["time_ms"] = flow_df["time_s"] * 1000
        log.info(f"  Flow profile: {Path(flow_profile_path).name} "
                 f"({len(flow_df)} pts, 0-{flow_df['time_ms'].max():.0f}ms)")

    # Determine camera positions from first frame
    airway_stls = sorted(motion_stl_dir.glob("*.stl"))
    if not airway_stls:
        log.warning("  No airway STLs, skipping 4-panel video")
        return

    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv

    ref = pv.read(str(airway_stls[0]))
    center = ref.center
    bounds = ref.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)
    y_range = bounds[3] - bounds[2]
    diag = np.linalg.norm(np.array([bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]]))

    # Clip box for zoom: nasopharynx region (upper portion)
    clip = [bounds[0] - 10, bounds[1] + 10,
            bounds[2] + y_range * 0.4, bounds[3] + 10,
            bounds[4] - 10, bounds[5] + 10]

    panel_size = (800, 800)

    # Load face labels for splitting deformed mesh by region
    face_labels_path = output_dir / f"{subject_id}_full_face_labels.npy"
    face_labels = np.load(str(face_labels_path)) if face_labels_path.exists() else None
    # Labels: 0=mouth, 1=left, 2=right, 3=descending

    def _save_region_stls_for_frame(deformed_stl_path, tmp_dir):
        """Split deformed mesh by face labels, save as temp region STLs."""
        mesh = trimesh.load_mesh(str(deformed_stl_path))
        paths = {}
        if face_labels is not None and len(face_labels) == len(mesh.faces):
            from slicer.septum_refine import extract_submesh
            region_map = {0: "Mouth", 3: "DescendingAirway", 1: "LeftNose", 2: "RightNose"}
            for label_val, region in region_map.items():
                sub = extract_submesh(mesh, face_labels == label_val)
                if len(sub.faces) > 0:
                    p = tmp_dir / f"{region}.stl"
                    sub.export(str(p))
                    paths[region] = str(p)
        else:
            paths["DescendingAirway"] = str(deformed_stl_path)
        return paths

    import tempfile
    tmp_region_dir = Path(tempfile.mkdtemp())

    log.info(f"  Rendering {len(plane_files)} 4-panel frames (PyVista)...")

    for fi, pfile in enumerate(plane_files):
        frame_name = pfile.stem.replace("-Planes-All", "")
        airway_file = motion_stl_dir / f"{frame_name}.stl"
        aw_stl = str(airway_file) if airway_file.exists() else str(airway_stls[0])

        # Split deformed mesh into colored regions
        aw_region_paths = _save_region_stls_for_frame(aw_stl, tmp_region_dir)

        # Top-left: Airway only (colored by region), sagittal
        img_tl = _render_pyvista_panel(
            airway_region_paths=aw_region_paths,
            window_size=panel_size, view="sagittal", zoom_factor=0.85)

        # Top-right: Airway (colored) + planes colored by region
        frame_planes_dir = output_dir / "frames" / frame_name
        region_plane_paths = {}
        for region, folder in [("DescendingAirway", "DescendingAirway_Planes"),
                               ("LeftNose", "LeftNose_Planes"),
                               ("RightNose", "RightNose_Planes")]:
            rd = frame_planes_dir / folder
            if rd.exists():
                region_plane_paths[region] = sorted(rd.glob("*.stl"))
        img_tr = _render_pyvista_panel(
            airway_region_paths=aw_region_paths,
            region_plane_paths=region_plane_paths,
            window_size=panel_size, view="sagittal", zoom_factor=0.85)

        # Bottom: CSA band + flow profile (dual y-axis)
        w, h = panel_size
        phase_ms = frame_name.split("_")[-1] if "_" in frame_name else str(fi)
        dpi = 100
        fig_plot, ax = plt.subplots(figsize=(w * 2 / dpi, h / dpi), dpi=dpi)

        for region, stats in region_stats.items():
            c = region_colors.get(region, "black")
            ax.fill_between(stats["arc_display"], stats["area_min"], stats["area_max"],
                            alpha=0.2, color=c)
            ax.plot(stats["arc_display"], stats["area_mean"],
                    "-", color=c, linewidth=1, alpha=0.6, label=f"{region} mean")

            if frame_name in frame_names:
                frame_rdf = df[(df["region"] == region) & (df["frame_name"] == frame_name)]
                frame_rdf = frame_rdf.sort_values("arc_length_mm")
                if len(frame_rdf) > 0:
                    ax.plot(arc_max - frame_rdf["arc_length_mm"].values,
                            frame_rdf["area_mm2"].values,
                            "-", color=c, linewidth=2.5, alpha=0.9)

        ax.set_xlabel("Distance from Nose (mm)", fontsize=11)
        ax.set_ylabel("CSA (mm²)", fontsize=11)
        ax.set_title("CSA — All Regions (band = min/max, thick = current frame)",
                      fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='upper left')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()

        fig_plot.canvas.draw()
        plot_img = np.frombuffer(fig_plot.canvas.tostring_rgb(), dtype=np.uint8)
        plot_img = plot_img.reshape(fig_plot.canvas.get_width_height()[::-1] + (3,))
        plt.close(fig_plot)

        # Composite
        def crop_panel(img, tw, th):
            img = img[:th, :tw, :3]
            if img.shape[0] < th or img.shape[1] < tw:
                padded = np.full((th, tw, 3), 255, dtype=np.uint8)
                padded[:img.shape[0], :img.shape[1]] = img
                return padded
            return img

        top_row = np.hstack([crop_panel(img_tl, w, h), crop_panel(img_tr, w, h)])

        # Overlay flow profile on top row
        if flow_df is not None:
            fig_flow, ax_flow = plt.subplots(figsize=(w * 2 / dpi, h / dpi), dpi=dpi)
            fig_flow.patch.set_alpha(0)
            ax_flow.patch.set_alpha(0)
            ax_flow.plot(flow_df["time_ms"], flow_df["flow_kgs"],
                        "-", color="darkgreen", linewidth=2.5, alpha=0.8)
            ax_flow.axhline(y=0, color="gray", linewidth=0.5, linestyle="--", alpha=0.3)
            current_ms = float(phase_ms)
            current_flow = np.interp(current_ms, flow_df["time_ms"], flow_df["flow_kgs"])
            ax_flow.plot(current_ms, current_flow, "o", color="red",
                        markersize=14, zorder=5, markeredgecolor="white", markeredgewidth=2)
            ax_flow.set_xlabel("Time (ms)", fontsize=10, color="darkgreen")
            ax_flow.set_ylabel("Flow (kg/s)", fontsize=10, color="darkgreen")
            ax_flow.tick_params(axis='both', labelcolor='darkgreen', labelsize=8)
            for spine in ['top', 'right']:
                ax_flow.spines[spine].set_visible(False)
            for spine in ['bottom', 'left']:
                ax_flow.spines[spine].set_color('darkgreen')
            plt.tight_layout()
            fig_flow.canvas.draw()
            fw, fh = fig_flow.canvas.get_width_height()
            flow_rgba = np.frombuffer(fig_flow.canvas.buffer_rgba(),
                                      dtype=np.uint8).reshape(fh, fw, 4)
            plt.close(fig_flow)
            from PIL import Image as PILImage
            flow_pil = PILImage.fromarray(flow_rgba).resize((w * 2, h), PILImage.LANCZOS)
            flow_arr = np.array(flow_pil)
            alpha = flow_arr[:, :, 3:4].astype(float) / 255.0
            blended = top_row.astype(float) * (1 - alpha * 0.8) + \
                      flow_arr[:, :, :3].astype(float) * (alpha * 0.8)
            top_row = blended.clip(0, 255).astype(np.uint8)

        bot_row = crop_panel(plot_img, w * 2, h)
        composite = np.vstack([top_row, bot_row])

        # Add title bar
        title_fig, title_ax = plt.subplots(figsize=(16, 0.6))
        title_ax.text(0.5, 0.5, f"{subject_id} — Breathing Phase {phase_ms}ms",
                      ha='center', va='center', fontsize=16, fontweight='bold')
        title_ax.set_axis_off()
        title_fig.canvas.draw()
        title_img = np.frombuffer(title_fig.canvas.tostring_rgb(), dtype=np.uint8)
        title_img = title_img.reshape(title_fig.canvas.get_width_height()[::-1] + (3,))
        plt.close(title_fig)
        title_pil = Image.fromarray(title_img).resize((w * 2, 40), Image.LANCZOS)
        title_resized = np.array(title_pil)

        final = np.vstack([title_resized, composite])

        # Save frame
        out_frame = frames_4p_dir / f"frame_{fi:04d}.png"
        Image.fromarray(final).save(str(out_frame))

        if (fi + 1) % 5 == 0 or fi == len(plane_files) - 1:
            log.info(f"    Frame {fi+1}/{len(plane_files)}")

    # Encode to MP4
    video_path = output_dir / f"{subject_id}_4panel_breathing.mp4"
    frame_pattern = str(frames_4p_dir / "frame_%04d.png")
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(fps),
        "-i", frame_pattern,
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        "-c:v", "libx264", "-pix_fmt", "yuv420p",
        "-crf", "20", "-preset", "medium",
        str(video_path),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            log.info(f"  4-panel video: {video_path.name}")
        else:
            log.warning(f"  ffmpeg failed: {result.stderr[:200]}")
    except FileNotFoundError:
        log.warning("  ffmpeg not found, skipping video encoding")


# ---------------------------------------------------------------------------
# Highlighted-planes 4-panel video (custom plane selection)
# ---------------------------------------------------------------------------

def generate_highlighted_video(combined_dir, motion_stl_dir, df, output_dir,
                               subject_id, highlight_planes, fps=5):
    """Create 4-panel video with specific planes highlighted.

    Args:
        highlight_planes: dict like {"DescendingAirway": [46, 66, 76], "LeftNose": [20]}
            Planes to highlight in red/yellow on the STL + vertical lines on CSA plot.
    """
    import subprocess
    from PIL import Image

    combined_dir = Path(combined_dir)
    motion_stl_dir = Path(motion_stl_dir)

    plane_files = sorted(combined_dir.glob("*-Planes-All.stl"))
    if not plane_files:
        log.warning("  No combined plane STLs, skipping highlighted video")
        return

    frames_hl_dir = output_dir / "video_frames_highlighted"
    frames_hl_dir.mkdir(parents=True, exist_ok=True)

    # Pre-compute CSA band per region
    region_stats = {}
    region_colors_plot = {"DescendingAirway": "gray", "LeftNose": "blue", "RightNose": "red"}
    for region in ["DescendingAirway", "LeftNose", "RightNose"]:
        rdf = df[df["region"] == region]
        if len(rdf) == 0:
            continue
        stats = rdf.groupby("plane_index").agg(
            area_mean=("area_mm2", "mean"),
            area_min=("area_mm2", "min"),
            area_max=("area_mm2", "max"),
            arc_mean=("arc_length_mm", "mean"),
        ).reset_index().sort_values("arc_mean")
        region_stats[region] = stats

    frame_names = sorted(df["frame_name"].unique()) if "frame_name" in df.columns else []

    # Renormalize arc length for display
    arc_max = df["arc_length_mm"].max()
    for region in region_stats:
        region_stats[region]["arc_display"] = arc_max - region_stats[region]["arc_mean"]

    # Build highlight arc lengths for vertical lines (in display coords)
    highlight_arcs = {}
    for region, indices in highlight_planes.items():
        rdf = df[df["region"] == region]
        for pi in indices:
            arc_raw = rdf[rdf["plane_index"] == pi]["arc_length_mm"].mean()
            arc_disp = arc_max - arc_raw
            label = {"DescendingAirway": "D", "LeftNose": "L", "RightNose": "R"}.get(region, "")
            highlight_arcs[f"{label}{pi}"] = arc_disp

    # Collect highlight plane STL paths per frame for PyVista rendering
    frames_base = output_dir / "frames"

    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv

    airway_stls = sorted(motion_stl_dir.glob("*.stl"))
    if not airway_stls:
        return

    w, h = 800, 800
    region_color_map = {"DescendingAirway": "gray", "LeftNose": "blue", "RightNose": "red", "Mouth": "orange"}

    # Load face labels for splitting deformed mesh
    face_labels_path = output_dir / f"{subject_id}_full_face_labels.npy"
    face_labels_hl = np.load(str(face_labels_path)) if face_labels_path.exists() else None

    def _split_deformed(stl_path, tmp_dir):
        mesh = trimesh.load_mesh(str(stl_path))
        paths = {}
        if face_labels_hl is not None and len(face_labels_hl) == len(mesh.faces):
            from slicer.septum_refine import extract_submesh
            for lv, region in {0: "Mouth", 3: "DescendingAirway", 1: "LeftNose", 2: "RightNose"}.items():
                sub = extract_submesh(mesh, face_labels_hl == lv)
                if len(sub.faces) > 0:
                    p = tmp_dir / f"{region}.stl"
                    sub.export(str(p))
                    paths[region] = str(p)
        else:
            paths["DescendingAirway"] = str(stl_path)
        return paths

    import tempfile
    tmp_hl_dir = Path(tempfile.mkdtemp())

    log.info(f"  Rendering {len(plane_files)} highlighted frames...")

    for fi, pfile in enumerate(plane_files):
        frame_name = pfile.stem.replace("-Planes-All", "")
        frame_dir = frames_base / frame_name
        airway_file = motion_stl_dir / f"{frame_name}.stl"
        aw_stl = str(airway_file) if airway_file.exists() else str(airway_stls[0])

        # Split deformed mesh by region
        aw_region_paths = _split_deformed(aw_stl, tmp_hl_dir)

        # --- Top-left: Airway (colored, moving) + highlighted planes only ---
        pl1 = pv.Plotter(off_screen=True, window_size=(w, h))
        pl1.set_background('white')
        for region, path in aw_region_paths.items():
            pl1.add_mesh(pv.read(str(path)), color=region_color_map[region],
                        opacity=0.12, smooth_shading=True)

        # Add highlighted planes in bright colors with labels
        hl_colors = ['red', 'yellow', 'lime', 'cyan', 'magenta', 'orange']
        short_prefix = {"DescendingAirway": "D", "LeftNose": "L", "RightNose": "R"}
        ci = 0
        for region, indices in highlight_planes.items():
            prefix = {"DescendingAirway": "Desc", "LeftNose": "Left", "RightNose": "Right"}[region]
            region_plane_dir = frame_dir / f"{region}_Planes"
            if not region_plane_dir.exists():
                continue
            rdf = df[(df["region"] == region) & (df["frame_name"] == frame_name)].sort_values("plane_index")
            plane_to_seq = {pi: si for si, pi in enumerate(rdf["plane_index"])}

            for pi in indices:
                seq = plane_to_seq.get(pi)
                if seq is None:
                    continue
                stl_file = region_plane_dir / f"{frame_name}-{prefix}-{seq:03d}.stl"
                if stl_file.exists():
                    mesh = pv.read(str(stl_file))
                    color = hl_colors[ci % len(hl_colors)]
                    pl1.add_mesh(mesh, color=color,
                                opacity=0.9, show_edges=True, edge_color='black', line_width=1)
                    # Label at centroid
                    row = rdf[rdf["plane_index"] == pi]
                    if len(row) > 0:
                        centroid = row[["centroid_x", "centroid_y", "centroid_z"]].values[0]
                        label = f"{short_prefix.get(region, '')}{pi}"
                        pl1.add_point_labels(
                            np.array([centroid]), [label],
                            font_size=20, text_color=color,
                            point_size=0, shape=None,
                            always_visible=True)
                ci += 1

        pl1.view_yz(); pl1.camera.up = (0, 0, 1)
        pl1.enable_parallel_projection(); pl1.reset_camera()
        pos = list(pl1.camera_position)
        pos[0] = (-pos[0][0], pos[0][1], pos[0][2])
        pl1.camera_position = pos
        pl1.camera.zoom(0.85)
        img_tl = pl1.screenshot(return_img=True); pl1.close()

        # --- Top-right: Airway (colored, moving) + all planes + highlighted ---
        pl2 = pv.Plotter(off_screen=True, window_size=(w, h))
        pl2.set_background('white')
        for region, path in aw_region_paths.items():
            pl2.add_mesh(pv.read(str(path)), color=region_color_map[region],
                        opacity=0.08, smooth_shading=True)

        # All planes (subdued)
        all_pm = pv.read(str(pfile))
        pl2.add_mesh(all_pm, color='gray', opacity=0.2, show_edges=False)

        # Highlighted planes (bright) with labels
        ci = 0
        for region, indices in highlight_planes.items():
            prefix = {"DescendingAirway": "Desc", "LeftNose": "Left", "RightNose": "Right"}[region]
            region_plane_dir = frame_dir / f"{region}_Planes"
            if not region_plane_dir.exists():
                continue
            rdf = df[(df["region"] == region) & (df["frame_name"] == frame_name)].sort_values("plane_index")
            plane_to_seq = {pi: si for si, pi in enumerate(rdf["plane_index"])}

            for pi in indices:
                seq = plane_to_seq.get(pi)
                if seq is None:
                    continue
                stl_file = region_plane_dir / f"{frame_name}-{prefix}-{seq:03d}.stl"
                if stl_file.exists():
                    mesh = pv.read(str(stl_file))
                    color = hl_colors[ci % len(hl_colors)]
                    pl2.add_mesh(mesh, color=color,
                                opacity=0.9, show_edges=True, edge_color='black', line_width=1)
                    row = rdf[rdf["plane_index"] == pi]
                    if len(row) > 0:
                        centroid = row[["centroid_x", "centroid_y", "centroid_z"]].values[0]
                        label = f"{short_prefix.get(region, '')}{pi}"
                        pl2.add_point_labels(
                            np.array([centroid]), [label],
                            font_size=20, text_color=color,
                            point_size=0, shape=None,
                            always_visible=True)
                ci += 1

        pl2.view_yz(); pl2.camera.up = (0, 0, 1)
        pl2.enable_parallel_projection(); pl2.reset_camera()
        pos = list(pl2.camera_position)
        pos[0] = (-pos[0][0], pos[0][1], pos[0][2])
        pl2.camera_position = pos
        pl2.camera.zoom(0.85)
        img_tr = pl2.screenshot(return_img=True); pl2.close()

        # --- Bottom: CSA plot with highlight vertical lines ---
        dpi = 100
        fig_plot, ax = plt.subplots(figsize=(w * 2 / dpi, h / dpi), dpi=dpi)

        for region, stats in region_stats.items():
            c = region_colors_plot.get(region, "black")
            ax.fill_between(stats["arc_display"], stats["area_min"], stats["area_max"],
                            alpha=0.2, color=c)
            ax.plot(stats["arc_display"], stats["area_mean"],
                    "-", color=c, linewidth=1, alpha=0.6)

            if frame_name in frame_names:
                frame_rdf = df[(df["region"] == region) & (df["frame_name"] == frame_name)]
                frame_rdf = frame_rdf.sort_values("arc_length_mm")
                if len(frame_rdf) > 0:
                    ax.plot(arc_max - frame_rdf["arc_length_mm"].values,
                            frame_rdf["area_mm2"].values,
                            "-", color=c, linewidth=2.5, alpha=0.9)

        # Vertical lines for highlighted planes (already in display coords)
        ci = 0
        for label, arc_disp in highlight_arcs.items():
            ax.axvline(x=arc_disp, color=hl_colors[ci % len(hl_colors)],
                       linewidth=2, linestyle="--", alpha=0.8)
            ax.text(arc_disp, ax.get_ylim()[1] * 0.95, f" {label}",
                    color=hl_colors[ci % len(hl_colors)],
                    fontsize=11, fontweight='bold', va='top')
            ci += 1

        ax.set_xlabel("Distance from Nose (mm)", fontsize=11)
        ax.set_ylabel("CSA (mm²)", fontsize=11)
        ax.set_title("CSA — Highlighted Planes", fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()

        fig_plot.canvas.draw()
        plot_img = np.frombuffer(fig_plot.canvas.tostring_rgb(), dtype=np.uint8)
        plot_img = plot_img.reshape(fig_plot.canvas.get_width_height()[::-1] + (3,))
        plt.close(fig_plot)

        # Composite
        def crop_panel(img, tw, th):
            img = img[:th, :tw, :3]
            if img.shape[0] < th or img.shape[1] < tw:
                padded = np.full((th, tw, 3), 255, dtype=np.uint8)
                padded[:img.shape[0], :img.shape[1]] = img
                return padded
            return img

        top_row = np.hstack([crop_panel(img_tl, w, h), crop_panel(img_tr, w, h)])
        bot_row = crop_panel(plot_img, w * 2, h)
        composite = np.vstack([top_row, bot_row])

        # Title bar
        phase_ms = frame_name.split("_")[-1] if "_" in frame_name else str(fi)
        hl_str = ", ".join(highlight_arcs.keys())
        title_fig, title_ax = plt.subplots(figsize=(16, 0.6))
        title_ax.text(0.5, 0.5, f"{subject_id} — Phase {phase_ms}ms — Planes: {hl_str}",
                      ha='center', va='center', fontsize=16, fontweight='bold')
        title_ax.set_axis_off()
        title_fig.canvas.draw()
        title_img = np.frombuffer(title_fig.canvas.tostring_rgb(), dtype=np.uint8)
        title_img = title_img.reshape(title_fig.canvas.get_width_height()[::-1] + (3,))
        plt.close(title_fig)
        title_pil = Image.fromarray(title_img).resize((w * 2, 40), Image.LANCZOS)
        final = np.vstack([np.array(title_pil)[:, :, :3], composite])

        out_frame = frames_hl_dir / f"frame_{fi:04d}.png"
        Image.fromarray(final).save(str(out_frame))

        if (fi + 1) % 5 == 0 or fi == len(plane_files) - 1:
            log.info(f"    Frame {fi+1}/{len(plane_files)}")

    # Encode
    hl_tag = "_".join(f"{k}" for k in highlight_arcs.keys())
    video_path = output_dir / f"{subject_id}_highlighted_{hl_tag}.mp4"
    cmd = [
        "ffmpeg", "-y", "-framerate", str(fps),
        "-i", str(frames_hl_dir / "frame_%04d.png"),
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        "-c:v", "libx264", "-pix_fmt", "yuv420p",
        "-crf", "20", "-preset", "medium", str(video_path),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            log.info(f"  Highlighted video: {video_path.name}")
        else:
            log.warning(f"  ffmpeg failed: {result.stderr[:200]}")
    except FileNotFoundError:
        log.warning("  ffmpeg not found")


# ---------------------------------------------------------------------------
# Plane reference PNG (labeled plane indices on STL)
# ---------------------------------------------------------------------------

def generate_plane_reference(df, output_dir, subject_id):
    """Generate 3-panel side-by-side plane reference with labeled indices."""
    import os
    os.environ['PYVISTA_OFF_SCREEN'] = 'true'
    import pyvista as pv
    from PIL import Image

    f0 = df[df['frame_name'] == df['frame_name'].iloc[0]].copy() if 'frame_name' in df.columns else df
    frames_dir = output_dir / "frames"
    first_frame = sorted(frames_dir.iterdir())[0] if frames_dir.exists() else None
    if first_frame is None:
        log.warning("  No frames found, skipping plane reference")
        return

    # Airway mesh (refined, colored by region + mouth from face labels)
    aw_paths = {
        "DescendingAirway": output_dir / f"{subject_id}_DescendingAirway_refined.stl",
        "LeftNose": output_dir / f"{subject_id}_LeftNose_refined.stl",
        "RightNose": output_dir / f"{subject_id}_RightNose_refined.stl",
    }
    # Add mouth if face labels exist
    face_labels_path = output_dir / f"{subject_id}_full_face_labels.npy"
    surface_dir = output_dir.parent / "surface"
    if face_labels_path.exists() and (surface_dir / "frame0.stl").exists():
        fl = np.load(str(face_labels_path))
        full = trimesh.load_mesh(str(surface_dir / "frame0.stl"))
        if len(fl) == len(full.faces):
            from slicer.septum_refine import extract_submesh
            mouth = extract_submesh(full, fl == 0)
            mouth_path = output_dir / f"{subject_id}_Mouth_refined.stl"
            mouth.export(str(mouth_path))
            aw_paths["Mouth"] = mouth_path

    region_colors = {
        "DescendingAirway": "gray", "LeftNose": "blue",
        "RightNose": "red", "Mouth": "orange",
    }
    region_configs = [
        ("DescendingAirway", "Desc", "D", "gray", 3),
        ("LeftNose", "Left", "L", "blue", 2),
        ("RightNose", "Right", "R", "red", 2),
    ]

    panels = []
    for region, prefix, short, color, interval in region_configs:
        region_dir = first_frame / f"{region}_Planes"
        rdf = f0[f0['region'] == region].sort_values('plane_index')
        if len(rdf) == 0:
            continue

        plotter = pv.Plotter(off_screen=True, window_size=(1200, 2000))
        plotter.set_background('white')

        # Airway (all regions, transparent)
        for r, p in aw_paths.items():
            if p.exists():
                plotter.add_mesh(pv.read(str(p)), color=region_colors[r],
                                opacity=0.08, smooth_shading=True)

        # Region planes
        if region_dir.exists():
            for stl_file in sorted(region_dir.glob('*.stl')):
                plotter.add_mesh(pv.read(str(stl_file)), color=color,
                                opacity=0.5, show_edges=False)

        # Labels
        for _, row in rdf.iterrows():
            pi = row['plane_index']
            if pi % interval != 0:
                continue
            centroid = [row['centroid_x'], row['centroid_y'], row['centroid_z']]
            plotter.add_point_labels(
                np.array([centroid]), [f'{short}{pi}'],
                font_size=16, text_color=color,
                point_size=0, shape=None, always_visible=True)

        plotter.add_text(region, position='upper_edge', font_size=14, color=color)

        # Sagittal view with camera flip
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

    # Side by side, same height
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


# ---------------------------------------------------------------------------
# PDF report (per region)
# ---------------------------------------------------------------------------

def generate_pdf_report(df, comp_df, output_dir, subject_id):
    """Generate per-region PDF reports with key plots."""
    regions = [r for r in df["region"].unique() if r != "NasalCombined"]
    has_frames = "frame_name" in df.columns

    for region in regions:
        rdf = df[df["region"] == region]
        pdf_path = output_dir / f"{subject_id}_{region}_report.pdf"

        with PdfPages(str(pdf_path)) as pdf:
            # Page 1: CSA dynamics
            frames = sorted(rdf["frame_name"].unique()) if has_frames else [None]
            n_frames = len(frames)
            cmap = plt.cm.viridis(np.linspace(0, 1, max(n_frames, 1)))

            fig, ax = plt.subplots(figsize=(12, 6))
            for fi, frame in enumerate(frames):
                fdf = rdf[rdf["frame_name"] == frame] if frame else rdf
                fdf = fdf.sort_values("plane_index")
                ax.plot(fdf["plane_index"], fdf["area_mm2"],
                       "-", color=cmap[fi], linewidth=0.8, alpha=0.6)
            ax.set_xlabel("Plane Index")
            ax.set_ylabel("CSA (mm²)")
            ax.set_title(f"{subject_id} — {region} — CSA Dynamics ({n_frames} frames)")
            ax.grid(True, alpha=0.3)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

            # Page 2: Compliance
            if comp_df is not None:
                rcomp = comp_df[comp_df["region"] == region]
                if len(rcomp) > 0:
                    fig, ax = plt.subplots(figsize=(12, 5))
                    colors = ["red" if c > 30 else "orange" if c > 15 else "green"
                              for c in rcomp["compliance_percent"]]
                    ax.bar(rcomp["plane_index"], rcomp["compliance_percent"],
                           color=colors, alpha=0.8)
                    ax.axhline(y=30, color="red", linestyle="--", alpha=0.5)
                    ax.axhline(y=15, color="orange", linestyle="--", alpha=0.5)
                    ax.set_xlabel("Plane Index")
                    ax.set_ylabel("Compliance (%)")
                    ax.set_title(f"{region} — Compliance / Collapse Risk")
                    ax.grid(True, alpha=0.3)
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close()

            # Page 3: Circularity
            if "circularity" in rdf.columns and has_frames:
                plane_stats = rdf.groupby("plane_index").agg(
                    circ_mean=("circularity", "mean"),
                    circ_std=("circularity", "std"),
                ).reset_index()

                fig, ax = plt.subplots(figsize=(12, 5))
                ax.fill_between(plane_stats["plane_index"],
                               plane_stats["circ_mean"] - plane_stats["circ_std"],
                               plane_stats["circ_mean"] + plane_stats["circ_std"],
                               alpha=0.3, color="steelblue")
                ax.plot(plane_stats["plane_index"], plane_stats["circ_mean"],
                       "-", color="steelblue", linewidth=1.5)
                ax.set_xlabel("Plane Index")
                ax.set_ylabel("Circularity (1.0 = circle)")
                ax.set_title(f"{region} — Cross-Section Circularity")
                ax.set_ylim(0, 1.1)
                ax.grid(True, alpha=0.3)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close()

            # Page 4: Heatmap
            if has_frames:
                planes = sorted(rdf["plane_index"].unique())
                frame_indices = sorted(rdf["frame_index"].unique())
                if len(frame_indices) >= 2:
                    matrix = np.full((len(planes), len(frame_indices)), np.nan)
                    pmap = {p: i for i, p in enumerate(planes)}
                    fmap = {f: i for i, f in enumerate(frame_indices)}
                    for _, row in rdf.iterrows():
                        pi = pmap.get(row["plane_index"])
                        fi = fmap.get(row["frame_index"])
                        if pi is not None and fi is not None:
                            matrix[pi, fi] = row["area_mm2"]

                    fig, ax = plt.subplots(figsize=(12, 8))
                    im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd",
                                  interpolation="nearest", origin="lower")
                    ax.set_xlabel("Frame Index")
                    ax.set_ylabel("Plane Index")
                    ax.set_title(f"{region} — CSA Heatmap")
                    plt.colorbar(im, ax=ax, label="CSA (mm²)")
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close()

        log.info(f"  PDF report: {pdf_path.name}")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_postprocessing(output_dir, subject_id):
    """Run all post-processing on pipeline output."""
    output_dir = Path(output_dir)

    # Find the CSV
    csv_path = output_dir / f"{subject_id}-Data.csv"
    if not csv_path.exists():
        log.error(f"CSV not found: {csv_path}")
        return

    df = pd.read_csv(csv_path)
    log.info(f"Loaded {len(df)} rows from {csv_path.name}")
    log.info(f"  Regions: {df['region'].value_counts().to_dict()}")

    n_frames = df["frame_name"].nunique() if "frame_name" in df.columns else 1
    log.info(f"  Frames: {n_frames}")

    # Per-region CSA plots
    log.info("\nGenerating CSA plots...")
    plot_csa_by_index_per_region(df, output_dir, subject_id)

    # Band plots (min/mean/max through cycle)
    log.info("\nGenerating CSA band plots...")
    plot_csa_band_per_region(df, output_dir, subject_id)

    # Heatmaps
    log.info("\nGenerating heatmaps...")
    plot_csa_heatmap_per_region(df, output_dir, subject_id)

    # Enhanced metrics
    log.info("\nComputing enhanced metrics...")
    df = compute_enhanced_metrics(df, output_dir, subject_id)

    # Compliance
    log.info("\nComputing compliance...")
    comp_df = compute_compliance_per_region(df, output_dir, subject_id)

    # Centroid movement
    log.info("\nComputing centroid movement...")
    compute_centroid_movement(df, output_dir, subject_id)

    # Breathing stats
    log.info("\nComputing breathing statistics...")
    compute_breathing_stats(df, output_dir, subject_id)

    # Combined plane STLs + video
    frames_dir = output_dir / "frames"
    motion_stl_dir = output_dir.parent / "motion" / "stl"
    if frames_dir.exists():
        log.info("\nGenerating combined plane STLs...")
        combined_dir = generate_combined_plane_stls(frames_dir, output_dir, subject_id)

        log.info("\nGenerating 4-panel breathing video...")
        generate_4panel_video(combined_dir, motion_stl_dir, df, output_dir,
                              subject_id)

    # Plane reference PNG (labeled plane indices on STL)
    log.info("\nGenerating plane reference PNG...")
    generate_plane_reference(df, output_dir, subject_id)

    # PDF reports
    log.info("\nGenerating PDF reports...")
    generate_pdf_report(df, comp_df, output_dir, subject_id)

    log.info(f"\nPost-processing complete.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Post-process bifurcation CSA results",
    )
    parser.add_argument("output_dir", help="CSA output directory (e.g., ENT001/csa_bifurcation/)")
    parser.add_argument("--subject-id", default="")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    subject_id = args.subject_id or output_dir.parent.name

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    run_postprocessing(output_dir, subject_id)
