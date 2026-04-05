#!/usr/bin/env python3
"""
Enhanced Airway Dynamics Analysis
==================================

Analyzes anatomical changes of upper airway during breathing cycle.
Computes enhanced metrics beyond basic CSA measurements:

1. Shape Metrics: circularity, solidity, eccentricity
2. Motion Metrics: 3D trajectories, displacement, rotation
3. Deformation Metrics: compliance, collapse risk, variability

Generates comprehensive visualizations for each metric type.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import glob
import sys
import trimesh
from typing import List, Tuple, Dict
from scipy.spatial import ConvexHull
from matplotlib.gridspec import GridSpec

sys.path.insert(0, 'python_slicer')
from slicer.io_utils import read_vtk_centerline


def compute_circularity(area: float, perimeter: float) -> float:
    """
    Circularity = 4π × Area / Perimeter²
    1.0 = perfect circle, <1 = more elongated
    """
    if perimeter <= 0:
        return 0.0
    return (4.0 * np.pi * area) / (perimeter ** 2)


def compute_solidity(vertices_2d: np.ndarray, area: float) -> float:
    """
    Solidity = Area / ConvexHull_Area
    Measures concavity/irregularity (1.0 = convex, <1 = concave)
    """
    try:
        if len(vertices_2d) < 3:
            return 1.0

        hull = ConvexHull(vertices_2d)
        convex_area = hull.volume  # In 2D, 'volume' is actually area

        if convex_area <= 0:
            return 1.0

        return min(area / convex_area, 1.0)
    except Exception:
        return 1.0


def compute_eccentricity(major_axis: float, minor_axis: float) -> float:
    """
    Eccentricity = sqrt(1 - (b²/a²))
    where a = major axis, b = minor axis
    0 = circle, approaching 1 = very elongated
    """
    if major_axis <= 0 or minor_axis <= 0:
        return 0.0

    ratio = minor_axis / major_axis if major_axis > 0 else 1.0
    return np.sqrt(1 - ratio**2)


def load_plane_vertices_from_stl(stl_path: str) -> np.ndarray:
    """Load 2D vertices from plane STL file"""
    try:
        mesh = trimesh.load_mesh(stl_path)
        vertices = mesh.vertices

        # Project to 2D (use first two dimensions)
        vertices_2d = vertices[:, :2]

        return vertices_2d
    except Exception as e:
        print(f"  Warning: Could not load {stl_path}: {e}")
        return np.array([])


def compute_plane_normal(stl_path: str) -> np.ndarray:
    """Estimate plane normal from STL triangles"""
    try:
        mesh = trimesh.load_mesh(stl_path)
        # Average face normals (weighted by area)
        normals = mesh.face_normals
        areas = mesh.area_faces

        if len(areas) == 0:
            return np.array([0, 0, 1])

        weighted_normal = np.average(normals, axis=0, weights=areas)
        weighted_normal /= np.linalg.norm(weighted_normal)

        return weighted_normal
    except Exception:
        return np.array([0, 0, 1])


def analyze_enhanced_metrics(subject: str, partition: str, output_dir: Path) -> pd.DataFrame:
    """
    Compute enhanced metrics from existing STL files and CSVs

    Returns enhanced DataFrame with all metrics
    """
    print(f"\n{'='*80}")
    print(f"ENHANCED AIRWAY DYNAMICS ANALYSIS")
    print(f"{'='*80}")
    print(f"Subject: {subject}")
    print(f"Partition: {partition}")
    print()

    # Find all Data CSV files
    csv_pattern = f"{partition}SlicedSTLs/*-Data.csv"
    csv_files = sorted(glob.glob(csv_pattern))

    if not csv_files:
        print(f"Error: No CSV files found matching {csv_pattern}")
        return pd.DataFrame()

    print(f"Found {len(csv_files)} time points")

    # Load all CSVs
    all_data = []
    for csv_file in csv_files:
        time_point = Path(csv_file).stem.replace('-Data', '')
        df = pd.read_csv(csv_file)
        df['time_point'] = time_point
        all_data.append(df)

    # Combine all time points
    combined_df = pd.concat(all_data, ignore_index=True)

    print(f"Total rows: {len(combined_df)}")
    print(f"Unique planes: {combined_df['plane_index'].nunique()}")

    # Compute enhanced shape metrics
    print(f"\nComputing enhanced shape metrics...")

    enhanced_metrics = []

    for idx, row in combined_df.iterrows():
        if (idx + 1) % 100 == 0:
            print(f"  Processed {idx+1}/{len(combined_df)} rows...")

        # Circularity
        circularity = compute_circularity(row['area_mm2'], row['perimeter_mm'])

        # Eccentricity from existing major/minor axis
        if row['major_axis_mm'] > 0 and row['minor_axis_mm'] > 0:
            eccentricity = np.sqrt(1 - (row['minor_axis_mm'] / row['major_axis_mm'])**2)
        else:
            eccentricity = 0.0

        # Load STL to compute solidity
        stl_pattern = f"{partition}SlicedSTLs/{row['time_point']}-Planes-{row['plane_index']:03d}.stl"
        stl_files = glob.glob(stl_pattern)

        if stl_files:
            vertices_2d = load_plane_vertices_from_stl(stl_files[0])
            solidity = compute_solidity(vertices_2d, row['area_mm2'])

            # Also get plane normal
            normal = compute_plane_normal(stl_files[0])
        else:
            solidity = 1.0
            normal = np.array([0, 0, 1])

        enhanced_metrics.append({
            **row.to_dict(),
            'circularity': circularity,
            'eccentricity': eccentricity,
            'solidity': solidity,
            'normal_x': normal[0],
            'normal_y': normal[1],
            'normal_z': normal[2]
        })

    enhanced_df = pd.DataFrame(enhanced_metrics)

    print(f"✓ Enhanced metrics computed")

    # Save enhanced CSV
    output_csv = output_dir / f"{subject}_{partition}_enhanced_metrics.csv"
    enhanced_df.to_csv(output_csv, index=False)
    print(f"✓ Saved to: {output_csv}")

    return enhanced_df


def compute_summary_metrics(enhanced_df: pd.DataFrame, output_dir: Path, subject: str, partition: str):
    """
    Compute per-plane summary statistics across all time points
    """
    print(f"\nComputing per-plane summary statistics...")

    summary_metrics = []

    for plane_idx in sorted(enhanced_df['plane_index'].unique()):
        plane_data = enhanced_df[enhanced_df['plane_index'] == plane_idx]

        # CSA statistics
        area_mean = plane_data['area_mm2'].mean()
        area_std = plane_data['area_mm2'].std()
        area_min = plane_data['area_mm2'].min()
        area_max = plane_data['area_mm2'].max()
        area_range = area_max - area_min
        area_cv = (area_std / area_mean) * 100 if area_mean > 0 else 0  # Coefficient of variation (%)

        # Compliance (% change from mean)
        compliance = (area_range / area_mean) * 100 if area_mean > 0 else 0

        # Shape variability
        circularity_mean = plane_data['circularity'].mean()
        circularity_std = plane_data['circularity'].std()

        eccentricity_mean = plane_data['eccentricity'].mean()
        solidity_mean = plane_data['solidity'].mean()

        # 3D motion
        centroid_x = plane_data['centroid_x'].values
        centroid_y = plane_data['centroid_y'].values
        centroid_z = plane_data['centroid_z'].values

        # Total path length (trajectory)
        displacements = np.sqrt(
            np.diff(centroid_x)**2 +
            np.diff(centroid_y)**2 +
            np.diff(centroid_z)**2
        )
        total_path_length = np.sum(displacements)

        # Max displacement from mean position
        mean_pos = np.array([centroid_x.mean(), centroid_y.mean(), centroid_z.mean()])
        positions = np.column_stack([centroid_x, centroid_y, centroid_z])
        displacements_from_mean = np.linalg.norm(positions - mean_pos, axis=1)
        max_displacement = displacements_from_mean.max()

        # Normal vector rotation (angle change)
        normals = np.column_stack([
            plane_data['normal_x'].values,
            plane_data['normal_y'].values,
            plane_data['normal_z'].values
        ])

        # Compute max rotation angle between normals
        max_rotation = 0.0
        for i in range(len(normals)):
            for j in range(i+1, len(normals)):
                dot_product = np.dot(normals[i], normals[j])
                dot_product = np.clip(dot_product, -1.0, 1.0)
                angle = np.degrees(np.arccos(dot_product))
                max_rotation = max(max_rotation, angle)

        summary_metrics.append({
            'plane_index': plane_idx,
            'arc_length_mm': plane_data['arc_length_mm'].iloc[0],
            'area_mean_mm2': area_mean,
            'area_std_mm2': area_std,
            'area_min_mm2': area_min,
            'area_max_mm2': area_max,
            'area_range_mm2': area_range,
            'area_cv_percent': area_cv,
            'compliance_percent': compliance,
            'circularity_mean': circularity_mean,
            'circularity_std': circularity_std,
            'eccentricity_mean': eccentricity_mean,
            'solidity_mean': solidity_mean,
            'total_path_length_mm': total_path_length,
            'max_displacement_mm': max_displacement,
            'max_rotation_deg': max_rotation
        })

    summary_df = pd.DataFrame(summary_metrics)

    # Save summary
    output_csv = output_dir / f"{subject}_{partition}_summary_metrics.csv"
    summary_df.to_csv(output_csv, index=False)
    print(f"✓ Saved summary to: {output_csv}")

    return summary_df


def plot_enhanced_metrics(enhanced_df: pd.DataFrame, summary_df: pd.DataFrame,
                          output_dir: Path, subject: str, partition: str):
    """
    Create comprehensive visualization plots
    """
    print(f"\nGenerating visualization plots...")

    plot_dir = output_dir / f"{partition}_enhanced_plots"
    plot_dir.mkdir(exist_ok=True)

    # Get unique planes sorted by arc length
    planes = sorted(summary_df['plane_index'].unique())
    arc_lengths = summary_df.set_index('plane_index')['arc_length_mm']

    # ============================================================================
    # PLOT 1: CSA Dynamics (Area over time for each plane)
    # ============================================================================
    print(f"  Creating Plot 1: CSA Dynamics...")

    fig, ax = plt.subplots(figsize=(14, 8))

    for plane_idx in planes:
        plane_data = enhanced_df[enhanced_df['plane_index'] == plane_idx].sort_values('time_point')
        time_indices = range(len(plane_data))
        ax.plot(time_indices, plane_data['area_mm2'], marker='o', markersize=3,
                label=f"Plane {plane_idx} (@ {arc_lengths[plane_idx]:.1f}mm)", alpha=0.7)

    ax.set_xlabel('Time Point (Breathing Cycle)', fontsize=12)
    ax.set_ylabel('Cross-Sectional Area (mm²)', fontsize=12)
    ax.set_title(f'Airway CSA Dynamics During Breathing - {partition}', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

    plt.tight_layout()
    plt.savefig(plot_dir / "01_CSA_dynamics.png", dpi=150, bbox_inches='tight')
    plt.close()

    # ============================================================================
    # PLOT 2: Shape Metrics Along Airway
    # ============================================================================
    print(f"  Creating Plot 2: Shape Metrics...")

    fig, axes = plt.subplots(3, 1, figsize=(14, 12))

    # Circularity
    axes[0].errorbar(arc_lengths, summary_df['circularity_mean'],
                     yerr=summary_df['circularity_std'],
                     marker='o', capsize=5, linewidth=2, markersize=6)
    axes[0].set_ylabel('Circularity\n(1=circle)', fontsize=11)
    axes[0].set_title('Shape Metrics Along Airway', fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Perfect Circle')
    axes[0].legend()

    # Eccentricity
    axes[1].plot(arc_lengths, summary_df['eccentricity_mean'],
                 marker='s', linewidth=2, markersize=6, color='orange')
    axes[1].set_ylabel('Eccentricity\n(0=circle, 1=line)', fontsize=11)
    axes[1].grid(True, alpha=0.3)

    # Solidity
    axes[2].plot(arc_lengths, summary_df['solidity_mean'],
                 marker='^', linewidth=2, markersize=6, color='green')
    axes[2].set_xlabel('Arc Length Along Centerline (mm)', fontsize=12)
    axes[2].set_ylabel('Solidity\n(1=convex)', fontsize=11)
    axes[2].grid(True, alpha=0.3)
    axes[2].axhline(y=1.0, color='r', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(plot_dir / "02_shape_metrics.png", dpi=150)
    plt.close()

    # ============================================================================
    # PLOT 3: Motion Metrics
    # ============================================================================
    print(f"  Creating Plot 3: 3D Motion...")

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Displacement
    axes[0].bar(arc_lengths, summary_df['max_displacement_mm'], width=1.0, alpha=0.7, color='steelblue')
    axes[0].set_ylabel('Max Displacement (mm)', fontsize=11)
    axes[0].set_title('3D Motion Metrics During Breathing', fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3, axis='y')

    # Path length
    axes[1].bar(arc_lengths, summary_df['total_path_length_mm'], width=1.0, alpha=0.7, color='coral')
    axes[1].set_xlabel('Arc Length Along Centerline (mm)', fontsize=12)
    axes[1].set_ylabel('Total Path Length (mm)', fontsize=11)
    axes[1].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(plot_dir / "03_motion_metrics.png", dpi=150)
    plt.close()

    # ============================================================================
    # PLOT 4: Compliance and Variability
    # ============================================================================
    print(f"  Creating Plot 4: Compliance...")

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Compliance (% change in CSA)
    colors = ['red' if c > 30 else 'orange' if c > 15 else 'green'
              for c in summary_df['compliance_percent']]

    axes[0].bar(arc_lengths, summary_df['compliance_percent'], width=1.0, alpha=0.7, color=colors)
    axes[0].set_ylabel('Compliance (%)', fontsize=11)
    axes[0].set_title('Airway Compliance and Collapse Risk', fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3, axis='y')
    axes[0].axhline(y=30, color='red', linestyle='--', alpha=0.7, linewidth=2, label='High Collapse Risk (>30%)')
    axes[0].axhline(y=15, color='orange', linestyle='--', alpha=0.7, linewidth=2, label='Moderate Risk (>15%)')
    axes[0].legend()

    # CSA Coefficient of Variation
    axes[1].plot(arc_lengths, summary_df['area_cv_percent'],
                 marker='o', linewidth=2, markersize=6, color='purple')
    axes[1].set_xlabel('Arc Length Along Centerline (mm)', fontsize=12)
    axes[1].set_ylabel('CSA Variability (CV%)', fontsize=11)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_dir / "04_compliance.png", dpi=150)
    plt.close()

    # ============================================================================
    # PLOT 5: 3D Centroid Trajectories (selected planes)
    # ============================================================================
    print(f"  Creating Plot 5: 3D Trajectories...")

    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Select evenly spaced planes for visualization (max 10)
    n_planes_to_plot = min(10, len(planes))
    plot_indices = np.linspace(0, len(planes)-1, n_planes_to_plot, dtype=int)
    selected_planes = [planes[i] for i in plot_indices]

    for plane_idx in selected_planes:
        plane_data = enhanced_df[enhanced_df['plane_index'] == plane_idx].sort_values('time_point')

        x = plane_data['centroid_x'].values
        y = plane_data['centroid_y'].values
        z = plane_data['centroid_z'].values

        ax.plot(x, y, z, marker='o', markersize=4, linewidth=2,
                label=f"Plane {plane_idx} (@ {arc_lengths[plane_idx]:.1f}mm)", alpha=0.8)

    ax.set_xlabel('X (mm)', fontsize=11)
    ax.set_ylabel('Y (mm)', fontsize=11)
    ax.set_zlabel('Z (mm)', fontsize=11)
    ax.set_title('3D Centroid Trajectories During Breathing', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(plot_dir / "05_3d_trajectories.png", dpi=150)
    plt.close()

    print(f"✓ All plots saved to: {plot_dir}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Enhanced airway dynamics analysis')
    parser.add_argument('subject', help='Subject ID (e.g., OSAMRI037)')
    parser.add_argument('partition', help='Partition name (e.g., LeftNoseDecending)')
    parser.add_argument('-o', '--output', default='.', help='Output directory')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)

    # Step 1: Compute enhanced metrics
    enhanced_df = analyze_enhanced_metrics(args.subject, args.partition, output_dir)

    if enhanced_df.empty:
        print("Error: No data processed")
        return 1

    # Step 2: Compute summary statistics
    summary_df = compute_summary_metrics(enhanced_df, output_dir, args.subject, args.partition)

    # Step 3: Generate plots
    plot_enhanced_metrics(enhanced_df, summary_df, output_dir, args.subject, args.partition)

    print(f"\n{'='*80}")
    print(f"ANALYSIS COMPLETE!")
    print(f"{'='*80}")
    print(f"Enhanced metrics: {args.subject}_{args.partition}_enhanced_metrics.csv")
    print(f"Summary metrics: {args.subject}_{args.partition}_summary_metrics.csv")
    print(f"Plots directory: {args.partition}_enhanced_plots/")
    print(f"{'='*80}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
