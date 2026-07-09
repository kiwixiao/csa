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
    """Compute airway collapse % per plane per region.

    Collapse = (CSAmax - CSAmin) / CSAmax × 100%
    This is the tracheomalacia criterion (Boiselle et al., Radiology 2006):
      >50% = tracheomalacia (abnormal collapse)
      >70% = severe tracheomalacia
      <50% = normal
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

            # Original compliance: (max-min)/mean
            compliance = (area_max - area_min) / area_mean * 100 if area_mean > 0 else 0
            # Tracheomalacia criterion: (max-min)/max (0-100%, clinically validated)
            collapse_pct = (area_max - area_min) / area_max * 100 if area_max > 0 else 0

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
                "collapse_percent": collapse_pct,
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
        df["eccentricity"] = np.sqrt(np.clip(1 - ratio**2, 0, None))
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

        for region_folder in ["DescendingAirway_Planes", "LeftNose_Planes",
                              "RightNose_Planes", "Mouth_Planes"]:
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
    region_colors = {"DescendingAirway": "gray", "LeftNose": "blue",
                     "RightNose": "red", "Mouth": "orange"}
    for region in ["DescendingAirway", "LeftNose", "RightNose", "Mouth"]:
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
            matches = [m for m in subject_dir.glob(pattern) if not m.name.startswith("._")]
            if matches:
                flow_profile_path = matches[0]
                break
    if flow_profile_path and Path(flow_profile_path).exists():
        _fdf = pd.read_csv(flow_profile_path)
        if _fdf.shape[1] >= 2:
            flow_df = _fdf.iloc[:, :2].copy()
            flow_df.columns = ["time_s", "flow_kgs"]
            flow_df["time_ms"] = flow_df["time_s"] * 1000
            log.info(f"  Flow profile: {Path(flow_profile_path).name} "
                     f"({len(flow_df)} pts, 0-{flow_df['time_ms'].max():.0f}ms)")
        else:
            flow_df = None
            log.warning(f"  Flow profile {Path(flow_profile_path).name} has "
                        f"{_fdf.shape[1]} column(s); expected >=2, skipping flow overlay")

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
    import shutil
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
                               ("RightNose", "RightNose_Planes"),
                               ("Mouth", "Mouth_Planes")]:
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
        plot_img = np.array(fig_plot.canvas.buffer_rgba())[:, :, :3]
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
        title_img = np.array(title_fig.canvas.buffer_rgba())[:, :, :3]
        plt.close(title_fig)
        title_pil = Image.fromarray(title_img).resize((w * 2, 40), Image.LANCZOS)
        title_resized = np.array(title_pil)

        final = np.vstack([title_resized, composite])

        # Save frame
        out_frame = frames_4p_dir / f"frame_{fi:04d}.png"
        Image.fromarray(final).save(str(out_frame))

        if (fi + 1) % 5 == 0 or fi == len(plane_files) - 1:
            log.info(f"    Frame {fi+1}/{len(plane_files)}")

    # Temp region STLs are no longer needed after all frames are rendered.
    shutil.rmtree(tmp_region_dir, ignore_errors=True)

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
    region_colors_plot = {"DescendingAirway": "gray", "LeftNose": "blue",
                          "RightNose": "red", "Mouth": "orange"}
    for region in ["DescendingAirway", "LeftNose", "RightNose", "Mouth"]:
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
    import shutil
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
        plot_img = np.array(fig_plot.canvas.buffer_rgba())[:, :, :3]
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
        title_img = np.array(title_fig.canvas.buffer_rgba())[:, :, :3]
        plt.close(title_fig)
        title_pil = Image.fromarray(title_img).resize((w * 2, 40), Image.LANCZOS)
        final = np.vstack([np.array(title_pil)[:, :, :3], composite])

        out_frame = frames_hl_dir / f"frame_{fi:04d}.png"
        Image.fromarray(final).save(str(out_frame))

        if (fi + 1) % 5 == 0 or fi == len(plane_files) - 1:
            log.info(f"    Frame {fi+1}/{len(plane_files)}")

    # Temp region STLs are no longer needed after all frames are rendered.
    shutil.rmtree(tmp_hl_dir, ignore_errors=True)

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
    surface_stls = list(surface_dir.glob("*.stl")) if surface_dir.exists() else []
    if face_labels_path.exists() and surface_stls:
        fl = np.load(str(face_labels_path))
        full = trimesh.load_mesh(str(surface_stls[0]))
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

def _render_mesh_with_planes(full_mesh_path, centroids, plane_indices, output_png, region):
    """Render full airway STL with region centerline and plane markers.

    Shows the complete airway surface (nose + descending) for context,
    with only the current region's centerline and plane labels drawn.
    Plane 0 = trachea (bottom), last plane = nasopharynx (top).
    """
    import pyvista as pv

    pv.OFF_SCREEN = True
    mesh = pv.read(str(full_mesh_path))

    plotter = pv.Plotter(off_screen=True, window_size=(800, 1600))
    plotter.set_background("white")

    # Full airway: semi-transparent for context
    plotter.add_mesh(mesh, color="lightsteelblue", opacity=0.3, show_edges=False)

    # Draw centerline through centroids (bold, visible)
    if len(centroids) > 1:
        cl_points = pv.PolyData(centroids)
        lines = np.column_stack([
            np.full(len(centroids) - 1, 2),
            np.arange(len(centroids) - 1),
            np.arange(1, len(centroids)),
        ]).ravel()
        cl_points.lines = lines
        plotter.add_mesh(cl_points, color="navy", line_width=3)

    # Add plane markers every 10 indices (less clutter)
    label_points = []
    label_texts = []
    for i, pi in enumerate(plane_indices):
        if int(pi) % 10 == 0:
            cx, cy, cz = centroids[i]
            sphere = pv.Sphere(radius=0.8, center=(cx, cy, cz))
            plotter.add_mesh(sphere, color="red", opacity=0.9)
            label_points.append((cx, cy, cz))
            label_texts.append(str(int(pi)))

    # Also mark first and last plane if not already included
    for idx in [0, len(plane_indices) - 1]:
        pi = plane_indices[idx]
        if int(pi) % 10 != 0:
            cx, cy, cz = centroids[idx]
            sphere = pv.Sphere(radius=0.8, center=(cx, cy, cz))
            plotter.add_mesh(sphere, color="darkred", opacity=0.9)
            label_points.append((cx, cy, cz))
            label_texts.append(str(int(pi)))

    if label_points:
        plotter.add_point_labels(
            label_points, label_texts,
            font_size=32, point_color="red", text_color="black",
            shape=None, show_points=False, always_visible=True,
            bold=True,
        )

    # Orient: sagittal view, full geometry visible
    # trachea (plane 0) at bottom, nasopharynx at top
    bounds = np.array(mesh.bounds)
    center = np.array([
        (bounds[0] + bounds[1]) / 2,
        (bounds[2] + bounds[3]) / 2,
        (bounds[4] + bounds[5]) / 2,
    ])
    extents = bounds[1::2] - bounds[::2]

    # Sagittal view: camera along -X axis (flipped side view)
    # Nose at top-left, trachea going down
    cam_dist = max(extents) * 3.0

    plotter.camera.position = (center[0] - cam_dist, center[1], center[2])
    plotter.camera.focal_point = center
    plotter.camera.up = (0, 0, 1)  # Z-up

    # Fit entire mesh in view, zoom to fill the panel height
    plotter.reset_camera()
    plotter.camera.zoom(1.8)

    plotter.screenshot(str(output_png))
    plotter.close()


def generate_pdf_report(df, comp_df, output_dir, subject_id, efd_csv=None):
    """Generate per-region single-page PDF reports.

    Layout: 3 columns
      Left:   STL mesh render with plane markers (pyvista)
      Middle: CSA band (min/mean/max), y = plane index
      Right:  Compliance bars (green/orange/red), y = plane index

    Plane 0 = trachea (bottom), last plane = nasopharynx (top).
    """
    regions = [r for r in df["region"].unique() if r != "NasalCombined"]
    has_frames = "frame_name" in df.columns

    # Find mesh files for rendering
    subject_dir = output_dir.parent
    branches_dir = subject_dir / "branches"

    for region in regions:
        rdf = df[df["region"] == region]
        pdf_path = output_dir / f"{subject_id}_{region}_report.pdf"

        # --- Compute per-plane stats across all frames ---
        plane_stats = rdf.groupby("plane_index").agg(
            area_mean=("area_mm2", "mean"),
            area_min=("area_mm2", "min"),
            area_max=("area_mm2", "max"),
            eq_diam_mean=("equivalent_diameter_mm", "mean"),
            cx=("centroid_x", "mean"),
            cy=("centroid_y", "mean"),
            cz=("centroid_z", "mean"),
        ).reset_index().sort_values("plane_index")

        planes = plane_stats["plane_index"].values
        n_planes = len(planes)
        if n_planes == 0:
            continue

        # --- Column 1: Render STL mesh ---
        mesh_png = output_dir / f"{subject_id}_{region}_mesh_render.png"

        # Find the full airway surface mesh (show complete anatomy for context)
        full_mesh_path = None
        surface_dir = subject_dir / "surface"
        if surface_dir.is_dir():
            for candidate in ["wall_mm.stl", "wall.stl", "seg_0.stl"]:
                if (surface_dir / candidate).exists():
                    full_mesh_path = surface_dir / candidate
                    break
            if full_mesh_path is None:
                surface_stls = list(surface_dir.glob("*.stl"))
                if surface_stls:
                    full_mesh_path = surface_stls[0]

        centroids = plane_stats[["cx", "cy", "cz"]].values
        mesh_image = None

        if full_mesh_path is not None:
            try:
                _render_mesh_with_planes(
                    full_mesh_path, centroids, planes, mesh_png, region,
                )
                mesh_image = plt.imread(str(mesh_png))
            except Exception as e:
                log.warning(f"  Mesh render failed for {region}: {e}")

        # --- Compute summary findings ---
        n_frames = len(rdf["frame_name"].unique()) if has_frames else 1
        min_csa_idx = plane_stats.loc[plane_stats["area_mean"].idxmin()]
        min_csa_plane = int(min_csa_idx["plane_index"])
        min_csa_val = min_csa_idx["area_mean"]

        # Compliance summary
        n_high_risk = 0
        max_comp_plane = None
        max_comp_val = 0
        if comp_df is not None:
            rcomp = comp_df[comp_df["region"] == region]
            if len(rcomp) > 0:
                n_high_risk = int((rcomp["collapse_percent"] > 50).sum())
                max_comp_row = rcomp.loc[rcomp["collapse_percent"].idxmax()]
                max_comp_plane = int(max_comp_row["plane_index"])
                max_comp_val = max_comp_row["collapse_percent"]

        summary = f"Narrowest: Plane {min_csa_plane} (mean CSA = {min_csa_val:.1f} mm²)"
        if max_comp_plane is not None:
            summary += f"  |  Max collapse: Plane {max_comp_plane} ({max_comp_val:.0f}%)"
            summary += f"  |  {n_high_risk}/{len(planes)} planes > 50% collapse"

        # --- Build figure: 3 columns ---
        fig = plt.figure(figsize=(20, 14))
        fig.suptitle(
            f"{subject_id} — {region} ({n_frames} frames)\n"
            f"Plane 0 = Trachea (bottom), Plane {int(planes[-1])} = Nasopharynx (top)\n"
            f"{summary}",
            fontsize=12, fontweight="bold", y=0.99,
        )

        gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.5, 1.0], wspace=0.12)

        # ── Column 1: Mesh image ──
        ax_geom = fig.add_subplot(gs[0])
        if mesh_image is not None:
            ax_geom.imshow(mesh_image, aspect="auto")
            ax_geom.set_title("Anatomy (STL)", fontsize=10)
        else:
            ax_geom.text(0.5, 0.5, "Mesh not\navailable",
                        transform=ax_geom.transAxes, ha="center", va="center",
                        fontsize=12, color="gray")
            ax_geom.set_title("Anatomy", fontsize=10)
        ax_geom.axis("off")

        # ── Column 2: CSA band (horizontal) ──
        # y = plane index, 0 at bottom (trachea), max at top (nasopharynx)
        ax_csa = fig.add_subplot(gs[1])
        ax_csa.fill_betweenx(
            planes, plane_stats["area_min"].values, plane_stats["area_max"].values,
            alpha=0.3, color="steelblue", label="Min–Max",
        )
        ax_csa.plot(
            plane_stats["area_mean"].values, planes,
            "-", color="steelblue", linewidth=1.5, label="Mean",
        )
        ax_csa.set_xlabel("CSA (mm²)", fontsize=11)
        ax_csa.set_ylabel("Plane Index\n(0 = Trachea → Nasopharynx)",
                          fontsize=12, fontweight="bold")
        ax_csa.set_title("Cross-Sectional Area", fontsize=11)
        ax_csa.set_ylim(0, None)  # start from 0
        ax_csa.set_xlim(0, None)  # start from 0
        ax_csa.grid(True, alpha=0.3)
        ax_csa.legend(fontsize=8, loc="upper right")
        ax_csa.tick_params(axis='y', labelsize=11, width=1.5)

        # ── Column 3: Tracheomalacia criterion ──
        ax_comp = fig.add_subplot(gs[2], sharey=ax_csa)
        if comp_df is not None:
            rcomp = comp_df[comp_df["region"] == region].sort_values("plane_index")
            if len(rcomp) > 0:
                comp_vals = rcomp["collapse_percent"].values
                comp_planes = rcomp["plane_index"].values
                colors = [
                    "red" if c > 70 else "orange" if c > 50 else "green"
                    for c in comp_vals
                ]
                ax_comp.barh(comp_planes, comp_vals, color=colors, alpha=0.8, height=0.8)

                # Color-banded background instead of just dashed lines
                ax_comp.axvspan(0, 50, color="green", alpha=0.06)
                ax_comp.axvspan(50, 70, color="orange", alpha=0.06)
                ax_comp.axvspan(70, 100, color="red", alpha=0.06)
                ax_comp.axvline(x=50, color="orange", linestyle="--", alpha=0.6, linewidth=1.0)
                ax_comp.axvline(x=70, color="red", linestyle="--", alpha=0.6, linewidth=1.0)

                # Threshold labels
                ax_comp.text(25, planes[-1] * 0.02, "Normal", ha="center",
                            fontsize=7, color="green", alpha=0.8)
                ax_comp.text(60, planes[-1] * 0.02, "Moderate", ha="center",
                            fontsize=7, color="orange", alpha=0.8)
                ax_comp.text(85, planes[-1] * 0.02, "Severe", ha="center",
                            fontsize=7, color="red", alpha=0.8)

                ax_comp.set_xlim(0, 100)
        ax_comp.set_xlabel("Airway Collapse (%)", fontsize=11)
        ax_comp.set_title("Airway Collapse (%)", fontsize=11)
        ax_comp.set_ylim(0, None)
        ax_comp.grid(True, alpha=0.3, axis="x")
        ax_comp.tick_params(axis='y', labelsize=11)

        # Footnote
        fig.text(
            0.5, 0.005,
            "Collapse = (CSAmax − CSAmin) / CSAmax × 100%.  "
            "Normal < 50%  |  Moderate 50–70%  |  Severe > 70%.  "
            "Collapse = (CSAmax − CSAmin) / CSAmax × 100%.",
            ha="center", fontsize=8, color="gray", style="italic",
        )


        with PdfPages(str(pdf_path)) as pdf:
            pdf.savefig(fig, bbox_inches="tight", dpi=150)
            plt.close()

            # Page 2: Work of Breathing (Joules per cycle, inhale vs exhale)
            if efd_csv is None:
                subject_dir = output_dir.parent
                for search_dir in [subject_dir, subject_dir.parent]:
                    for candidate in list(search_dir.glob(f"{subject_id}_EFDiff*.csv")) + \
                                     list(search_dir.glob("*EFDiff*.csv")):
                        if candidate.exists():
                            efd_csv = candidate
                            break
                    if efd_csv:
                        break

            if efd_csv is not None and Path(efd_csv).exists():
                try:
                    from plot_tpl_efd import load_and_clean
                    from scipy.integrate import trapezoid

                    efd_full = load_and_clean(efd_csv, "efd_")
                    time_efd = efd_full["Time"].values
                    efd_regions_full = [c for c in efd_full.columns if c != "Time"]

                    # Switch time (hardcoded for now — TODO: parameterize)
                    switch_time = 1.6
                    if not (time_efd[0] < switch_time < time_efd[-1]):
                        switch_time = None

                    if switch_time is not None:
                        inhale_mask = time_efd <= switch_time
                        exhale_mask = time_efd > switch_time

                        def _integrate(values, time, mask):
                            t_m, v_m = time[mask], np.abs(values[mask])
                            return trapezoid(v_m, t_m) if len(t_m) >= 2 else 0.0

                        wob_rows = []
                        for r in efd_regions_full:
                            vals = efd_full[r].values
                            e_total = _integrate(vals, time_efd, np.ones(len(time_efd), bool))
                            e_inh = _integrate(vals, time_efd, inhale_mask)
                            e_exh = _integrate(vals, time_efd, exhale_mask)
                            wob_rows.append({"region": r, "total_J": e_total,
                                             "inhale_J": e_inh, "exhale_J": e_exh})

                        wob_df = pd.DataFrame(wob_rows).set_index("region").sort_values("total_J", ascending=True)
                        total_e = wob_df["total_J"].sum()
                        inh_e = wob_df["inhale_J"].sum()
                        exh_e = wob_df["exhale_J"].sum()

                        colors_wob = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                                      "#ff7f00", "#a65628", "#f781bf", "#999999"]

                        fig_wob, (ax_t, ax_p) = plt.subplots(1, 2, figsize=(18, 7),
                            gridspec_kw={"width_ratios": [1.0, 1.5]})
                        fig_wob.suptitle(
                            f"{subject_id} — Work of Breathing (Energy per Cycle)\n"
                            f"Total = {total_e:.4f} J  |  "
                            f"Inhale = {inh_e:.4f} J (0–{switch_time}s)  |  "
                            f"Exhale = {exh_e:.4f} J ({switch_time}–{time_efd[-1]:.1f}s)",
                            fontsize=13, fontweight="bold")

                        # Left: total per cycle
                        bc = [colors_wob[efd_regions_full.index(r) % len(colors_wob)]
                              if r in efd_regions_full else "gray" for r in wob_df.index]
                        ax_t.barh(wob_df.index, wob_df["total_J"], color=bc, alpha=0.85)
                        for i, (r, row) in enumerate(wob_df.iterrows()):
                            pct = row["total_J"] / total_e * 100 if total_e > 0 else 0
                            ax_t.text(row["total_J"] + wob_df["total_J"].max() * 0.02, i,
                                      f"{row['total_J']:.4f} J ({pct:.0f}%)", va="center", fontsize=9)
                        ax_t.set_xlabel("Energy per Cycle (J)", fontsize=11, fontweight="bold")
                        ax_t.set_title("Total (Full Cycle)", fontsize=11)
                        ax_t.grid(True, alpha=0.3, axis="x")
                        ax_t.set_xlim(0, wob_df["total_J"].max() * 1.4)
                        ax_t.tick_params(axis="y", labelsize=11)

                        # Right: inhale vs exhale
                        y_pos = np.arange(len(wob_df))
                        bh = 0.35
                        ax_p.barh(y_pos + bh/2, wob_df["inhale_J"], height=bh,
                                  color="steelblue", alpha=0.8, label=f"Inhale (0–{switch_time}s)")
                        ax_p.barh(y_pos - bh/2, wob_df["exhale_J"], height=bh,
                                  color="coral", alpha=0.8, label=f"Exhale ({switch_time}–{time_efd[-1]:.1f}s)")
                        for i, (r, row) in enumerate(wob_df.iterrows()):
                            ip = row["inhale_J"] / inh_e * 100 if inh_e > 0 else 0
                            ep = row["exhale_J"] / exh_e * 100 if exh_e > 0 else 0
                            ax_p.text(row["inhale_J"] + wob_df["inhale_J"].max() * 0.02,
                                      i + bh/2, f"{ip:.0f}% ({row['inhale_J']:.4f} J)",
                                      va="center", fontsize=8, color="steelblue")
                            ax_p.text(row["exhale_J"] + wob_df["exhale_J"].max() * 0.02,
                                      i - bh/2, f"{ep:.0f}% ({row['exhale_J']:.4f} J)",
                                      va="center", fontsize=8, color="coral")
                        ax_p.set_yticks(y_pos)
                        ax_p.set_yticklabels(wob_df.index, fontsize=11)
                        ax_p.set_xlabel("Energy per Phase (J)", fontsize=11, fontweight="bold")
                        ax_p.set_title("Inhale vs Exhale", fontsize=11)
                        ax_p.legend(fontsize=9, loc="lower right")
                        ax_p.grid(True, alpha=0.3, axis="x")

                        # Summary
                        top_t = wob_df["total_J"].idxmax()
                        top_tp = wob_df.loc[top_t, "total_J"] / total_e * 100
                        top_i = wob_df["inhale_J"].idxmax()
                        top_ip = wob_df.loc[top_i, "inhale_J"] / inh_e * 100
                        top_e = wob_df["exhale_J"].idxmax()
                        top_ep = wob_df.loc[top_e, "exhale_J"] / exh_e * 100
                        ir = inh_e / exh_e if exh_e > 0 else 0
                        summary = (
                            f"The {top_t} accounts for {top_tp:.0f}% of total energy dissipation ({total_e:.3f} J/cycle). "
                            f"Inhalation requires {ir:.1f}× more energy than exhalation. "
                            f"During inhale, {top_i} dominates ({top_ip:.0f}%); "
                            f"during exhale, {top_e} dominates ({top_ep:.0f}%).")
                        fig_wob.text(0.5, 0.035, summary, ha="center", fontsize=9,
                                     fontweight="bold", color="black")
                        fig_wob.text(0.5, 0.005,
                            f"WoB = ∫|EFDiff|dt (Joules per cycle).  "
                            f"Inhale + Exhale = Total ({inh_e:.4f} + {exh_e:.4f} = {inh_e+exh_e:.4f} J).",
                            ha="center", fontsize=8, color="gray", style="italic")
                        plt.tight_layout(rect=[0, 0.05, 1, 0.90])
                        pdf.savefig(fig_wob, bbox_inches="tight", dpi=150)
                        plt.close()
                        log.info(f"  Added WoB page (J/cycle, switch={switch_time}s)")
                except Exception as e:
                    log.warning(f"  WoB page failed: {e}")

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
