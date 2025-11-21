#!/usr/bin/env python3
"""
Analyze Centerline Motion Across Breathing Cycle

Check if centerline points move through breathing phases due to FFD deformation.
"""

import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from slicer.io_utils import read_vtk_centerline


def analyze_centerline_motion(partition: str = "LeftNoseDecending"):
    """
    Analyze how centerline points move across breathing phases

    Args:
        partition: Partition name
    """
    # Find VTK files
    vtk_dir = Path(partition) / "FFD" / "vtk"
    vtk_files = sorted(glob.glob(str(vtk_dir / "*.vtk")))

    if not vtk_files:
        print(f"Error: No VTK files found in {vtk_dir}")
        return

    print(f"Found {len(vtk_files)} VTK files")

    # Load all centerlines
    centerlines = {}
    for vtk_path in vtk_files:
        file_id = Path(vtk_path).stem
        centerline = read_vtk_centerline(vtk_path)
        centerlines[file_id] = centerline
        print(f"{file_id}: {len(centerline)} points")

    # Check if all have same number of points
    n_points_list = [len(cl) for cl in centerlines.values()]
    if len(set(n_points_list)) > 1:
        print(f"\nWarning: Different centerline lengths: {set(n_points_list)}")
        print("Using minimum length for comparison")

    min_points = min(n_points_list)

    # Extract centerlines as arrays
    file_ids = sorted(centerlines.keys())
    n_phases = len(file_ids)

    # Build 3D array: [n_phases, n_points, 3]
    centerline_array = np.zeros((n_phases, min_points, 3))
    for i, file_id in enumerate(file_ids):
        centerline_array[i] = centerlines[file_id][:min_points]

    # Analyze motion at each centerline index
    print(f"\n{'='*80}")
    print("Centerline Point Motion Analysis")
    print(f"{'='*80}")

    motion_stats = []

    for idx in range(0, min_points, 10):  # Sample every 10th point
        # Get positions across all phases
        positions = centerline_array[:, idx, :]  # [n_phases, 3]

        # Compute range of motion
        x_range = positions[:, 0].max() - positions[:, 0].min()
        y_range = positions[:, 1].max() - positions[:, 1].min()
        z_range = positions[:, 2].max() - positions[:, 2].min()
        total_range = np.sqrt(x_range**2 + y_range**2 + z_range**2)

        # Compute std deviation
        std_3d = np.std(positions, axis=0)
        total_std = np.linalg.norm(std_3d)

        motion_stats.append({
            'centerline_index': idx,
            'x_range': x_range,
            'y_range': y_range,
            'z_range': z_range,
            'total_range': total_range,
            'total_std': total_std
        })

        if idx < 50:  # Print first few
            print(f"Index {idx:3d}: Range = {total_range:.2f} mm, Std = {total_std:.2f} mm")

    motion_df = pd.DataFrame(motion_stats)

    print(f"\nSummary:")
    print(f"Mean motion range: {motion_df['total_range'].mean():.2f} mm")
    print(f"Max motion range: {motion_df['total_range'].max():.2f} mm")
    print(f"Min motion range: {motion_df['total_range'].min():.2f} mm")

    # Save motion statistics
    motion_df.to_csv(f"{partition}_centerline_motion_stats.csv", index=False)
    print(f"\nSaved: {partition}_centerline_motion_stats.csv")

    # Plot motion vs centerline index
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Plot 1: Motion range
    ax = axes[0]
    ax.plot(motion_df['centerline_index'], motion_df['total_range'], 'b-', linewidth=2)
    ax.fill_between(motion_df['centerline_index'], 0, motion_df['total_range'], alpha=0.3)
    ax.set_xlabel('Centerline Index', fontsize=12, fontweight='bold')
    ax.set_ylabel('Motion Range (mm)', fontsize=12, fontweight='bold')
    ax.set_title(f'Centerline Point Motion Across Breathing Cycle\n{partition}',
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Plot 2: Component-wise ranges
    ax = axes[1]
    ax.plot(motion_df['centerline_index'], motion_df['x_range'], 'r-', label='X-range', linewidth=1.5)
    ax.plot(motion_df['centerline_index'], motion_df['y_range'], 'g-', label='Y-range', linewidth=1.5)
    ax.plot(motion_df['centerline_index'], motion_df['z_range'], 'b-', label='Z-range', linewidth=1.5)
    ax.set_xlabel('Centerline Index', fontsize=12, fontweight='bold')
    ax.set_ylabel('Component Range (mm)', fontsize=12, fontweight='bold')
    ax.set_title('Component-wise Motion', fontsize=12, fontweight='bold')
    ax.legend(loc='best', frameon=True)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{partition}_centerline_motion.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {partition}_centerline_motion.png")

    # Plot 3D trajectory for selected points
    plot_3d_trajectories(centerline_array, file_ids, partition, sample_indices=[0, 25, 50, 75, 100])

    return centerlines, centerline_array


def plot_3d_trajectories(centerline_array, file_ids, partition, sample_indices):
    """
    Plot 3D trajectories of selected centerline points

    Args:
        centerline_array: [n_phases, n_points, 3]
        file_ids: List of file IDs
        partition: Partition name
        sample_indices: Centerline indices to plot
    """
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    colors = plt.cm.viridis(np.linspace(0, 1, len(sample_indices)))

    for i, idx in enumerate(sample_indices):
        if idx >= centerline_array.shape[1]:
            continue

        # Get trajectory for this point
        trajectory = centerline_array[:, idx, :]  # [n_phases, 3]

        # Plot trajectory
        ax.plot(trajectory[:, 0], trajectory[:, 1], trajectory[:, 2],
                'o-', color=colors[i], linewidth=2, markersize=4,
                label=f'Index {idx}', alpha=0.7)

        # Mark start and end
        ax.scatter(trajectory[0, 0], trajectory[0, 1], trajectory[0, 2],
                   color=colors[i], s=100, marker='o', edgecolors='black', linewidths=2)
        ax.scatter(trajectory[-1, 0], trajectory[-1, 1], trajectory[-1, 2],
                   color=colors[i], s=100, marker='s', edgecolors='black', linewidths=2)

    ax.set_xlabel('X (mm)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Y (mm)', fontsize=11, fontweight='bold')
    ax.set_zlabel('Z (mm)', fontsize=11, fontweight='bold')
    ax.set_title(f'3D Centerline Point Trajectories\n{partition}\n(Circle=Start, Square=End)',
                 fontsize=13, fontweight='bold', pad=20)
    ax.legend(loc='best', frameon=True)

    plt.tight_layout()
    plt.savefig(f"{partition}_centerline_3d_trajectories.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {partition}_centerline_3d_trajectories.png")


if __name__ == "__main__":
    import sys
    partition = sys.argv[1] if len(sys.argv) > 1 else "LeftNoseDecending"
    analyze_centerline_motion(partition)
