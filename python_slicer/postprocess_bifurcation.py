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

def generate_plane_motion_video(combined_dir, motion_stl_dir, output_dir,
                                subject_id, fps=5):
    """Create MP4 video showing planes on the full airway through breathing cycle.

    Uses matplotlib 3D rendering. Shows deformed airway (transparent) +
    all cross-section planes (colored by position).
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import subprocess

    combined_dir = Path(combined_dir)
    motion_stl_dir = Path(motion_stl_dir)

    plane_files = sorted(combined_dir.glob("*-Planes-All.stl"))
    if not plane_files:
        log.warning("  No combined plane STLs found, skipping video")
        return

    video_frames_dir = output_dir / "video_frames"
    video_frames_dir.mkdir(parents=True, exist_ok=True)

    # Load first airway mesh for bounds
    first_airway = None
    airway_stls = sorted(motion_stl_dir.glob("*.stl"))
    if airway_stls:
        first_airway = trimesh.load_mesh(str(airway_stls[0]))

    # Determine global bounds from first frame
    if first_airway is not None:
        bounds_min = first_airway.bounds[0]
        bounds_max = first_airway.bounds[1]
    else:
        plane0 = trimesh.load_mesh(str(plane_files[0]))
        bounds_min = plane0.bounds[0] - 10
        bounds_max = plane0.bounds[1] + 10

    padding = 10
    xlim = (bounds_min[0] - padding, bounds_max[0] + padding)
    ylim = (bounds_min[1] - padding, bounds_max[1] + padding)
    zlim = (bounds_min[2] - padding, bounds_max[2] + padding)

    log.info(f"  Rendering {len(plane_files)} video frames...")

    for fi, pfile in enumerate(plane_files):
        frame_name = pfile.stem.replace("-Planes-All", "")

        # Load planes
        plane_mesh = trimesh.load_mesh(str(pfile))

        # Load corresponding deformed airway (subsample for speed)
        airway_file = motion_stl_dir / f"{frame_name}.stl"
        airway_mesh = None
        if airway_file.exists():
            airway_mesh = trimesh.load_mesh(str(airway_file))

        # Render
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Airway (subsampled, transparent)
        if airway_mesh is not None:
            step = max(1, len(airway_mesh.faces) // 5000)
            sub_faces = airway_mesh.faces[::step]
            airway_coll = Poly3DCollection(
                airway_mesh.vertices[sub_faces],
                alpha=0.08, facecolors='lightblue',
                edgecolors='gray', linewidths=0.05)
            ax.add_collection3d(airway_coll)

        # Planes (colored by Y position = along airway)
        if len(plane_mesh.faces) > 0:
            face_verts = plane_mesh.vertices[plane_mesh.faces]
            y_coords = face_verts.mean(axis=1)[:, 1]
            y_min, y_max = y_coords.min(), y_coords.max()
            colors = plt.cm.rainbow(
                (y_coords - y_min) / (y_max - y_min + 1e-6))
            plane_coll = Poly3DCollection(
                face_verts, alpha=0.6, facecolors=colors,
                edgecolors='black', linewidths=0.1)
            ax.add_collection3d(plane_coll)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)
        ax.view_init(elev=15, azim=45)
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
        ax.set_zlabel("Z (mm)")

        phase_ms = frame_name.split("_")[-1] if "_" in frame_name else str(fi)
        ax.set_title(f"{subject_id} — Breathing Phase {phase_ms}ms\n"
                     f"{len(plane_mesh.faces)} plane faces",
                     fontsize=12, fontweight="bold")

        out_frame = video_frames_dir / f"frame_{fi:04d}.png"
        plt.savefig(str(out_frame), dpi=100, bbox_inches="tight")
        plt.close()

        if (fi + 1) % 5 == 0 or fi == len(plane_files) - 1:
            log.info(f"    Frame {fi+1}/{len(plane_files)}")

    # Combine to MP4 with ffmpeg
    video_path = output_dir / f"{subject_id}_plane_motion.mp4"
    frame_pattern = str(video_frames_dir / "frame_%04d.png")
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(fps),
        "-i", frame_pattern,
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-crf", "23",
        "-preset", "medium",
        str(video_path),
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            log.info(f"  Video: {video_path.name}")
        else:
            log.warning(f"  ffmpeg failed: {result.stderr[:200]}")
    except FileNotFoundError:
        log.warning("  ffmpeg not found, skipping video encoding")


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

        log.info("\nGenerating breathing cycle video...")
        generate_plane_motion_video(combined_dir, motion_stl_dir, output_dir,
                                    subject_id)

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
