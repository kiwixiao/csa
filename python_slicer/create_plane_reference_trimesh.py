#!/usr/bin/env python3
"""
Create plane reference figure using trimesh with isometric view
All planes labeled for anatomical identification
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import glob
import sys
import argparse
from pathlib import Path
import trimesh
import re
import pandas as pd


def format_partition_name(partition_name):
    """Format partition name for display"""
    formatted = re.sub(r'([a-z])([A-Z])', r'\1 \2', partition_name)
    return formatted


def load_plane_mesh(output_dir, time_point, plane_idx):
    """Load a single plane STL file"""
    plane_stl = output_dir / f"{time_point}-Planes-{plane_idx:03d}.stl"
    if plane_stl.exists():
        try:
            return trimesh.load_mesh(str(plane_stl))
        except Exception as e:
            return None
    return None


def get_plane_indices(output_dir, time_point):
    """Get list of plane indices from CSV file"""
    csv_file = output_dir / f"{time_point}-Data.csv"
    if not csv_file.exists():
        return []
    df = pd.read_csv(csv_file)
    return df['plane_index'].tolist()


def load_centerline_from_csv(output_dir, time_point):
    """Load centerline coordinates from CSV file"""
    csv_file = output_dir / f"{time_point}-Data.csv"
    if not csv_file.exists():
        return None
    df = pd.read_csv(csv_file)
    centerline_points = df[['centroid_x', 'centroid_y', 'centroid_z']].values
    return centerline_points


def main():
    parser = argparse.ArgumentParser(
        description='Create reference figure with trimesh - isometric view'
    )
    parser.add_argument(
        'partition', nargs='?', default=None,
        help='Partition name (e.g., LeftNoseDecending, RightNose)'
    )
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    if args.partition:
        partition_name = args.partition
    else:
        sliced_dirs = list(project_root.glob("*SlicedSTLs"))
        if not sliced_dirs:
            print("Error: No SlicedSTLs directories found.")
            sys.exit(1)
        sliced_dirs_sorted = sorted(
            sliced_dirs,
            key=lambda x: (not x.name.startswith("LeftNose"), x.name)
        )
        partition_name = sliced_dirs_sorted[0].name.replace("SlicedSTLs", "")

    partition_display_name = format_partition_name(partition_name)

    print("=" * 70)
    print(f"Creating Reference Figure - {partition_display_name}")
    print("=" * 70)

    # Setup paths
    output_dir = project_root / f"{partition_name}SlicedSTLs"
    partition_dir = project_root / partition_name / "FFD" / "stl"
    airway_stl_files = sorted(glob.glob(str(partition_dir / "*.stl")))

    print(f"\nLoading geometries...")

    # Find first time point
    plane_stl_files = sorted(glob.glob(str(output_dir / "*-Planes-All.stl")))
    if not plane_stl_files:
        print(f"Error: No plane files found")
        return

    time_point = Path(plane_stl_files[0]).stem.replace('-Planes-All', '')

    # Load data
    plane_indices = get_plane_indices(output_dir, time_point)
    centerline_points = load_centerline_from_csv(output_dir, time_point)

    print(f"  ✓ {len(plane_indices)} planes")

    # Load plane meshes
    plane_data_list = []
    for plane_idx in plane_indices:
        plane_mesh = load_plane_mesh(output_dir, time_point, plane_idx)
        if plane_mesh is not None:
            centroid = plane_mesh.vertices.mean(axis=0)
            plane_data_list.append((plane_mesh, plane_idx, centroid))

    # Load airway
    airway_mesh = None
    if airway_stl_files:
        airway_mesh = trimesh.load_mesh(airway_stl_files[0])

    print(f"\nCreating figure...")

    # Create large figure - VERY LARGE for readability
    fig = plt.figure(figsize=(40, 30), dpi=150)
    ax = fig.add_subplot(111, projection='3d')

    # Plot airway (MORE VISIBLE - 30% opacity instead of 2%)
    if airway_mesh is not None:
        airway_verts = airway_mesh.vertices
        airway_faces = airway_mesh.faces[::4]  # Less subsampling for better detail

        airway_collection = Poly3DCollection(
            airway_verts[airway_faces],
            alpha=0.3,
            facecolors='lightblue',
            edgecolors='gray',
            linewidths=0.1
        )
        ax.add_collection3d(airway_collection)

    # Plot centerline (THICKER AND MORE VISIBLE)
    if centerline_points is not None and len(centerline_points) > 0:
        ax.plot(
            centerline_points[:, 0],
            centerline_points[:, 1],
            centerline_points[:, 2],
            'r-', linewidth=4, alpha=1.0, label='Centerline', zorder=100
        )

    # Plot planes
    for plane_mesh, plane_idx, centroid in plane_data_list:
        plane_collection = Poly3DCollection(
            plane_mesh.vertices[plane_mesh.faces],
            alpha=0.6,
            facecolors='orange',
            edgecolors='darkgray',
            linewidths=0.5
        )
        ax.add_collection3d(plane_collection)

    # Calculate bounds
    if airway_mesh is not None:
        x_min, x_max = airway_mesh.vertices[:, 0].min(), airway_mesh.vertices[:, 0].max()
        y_min, y_max = airway_mesh.vertices[:, 1].min(), airway_mesh.vertices[:, 1].max()
        z_min, z_max = airway_mesh.vertices[:, 2].min(), airway_mesh.vertices[:, 2].max()
    else:
        all_verts = np.vstack([p[0].vertices for p in plane_data_list])
        x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
        y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
        z_min, z_max = all_verts[:, 2].min(), all_verts[:, 2].max()

    # Label placement - TWO COLUMNS on right side
    num_planes = len(plane_data_list)
    planes_per_column = (num_planes + 1) // 2  # Split into 2 columns

    # Sort planes by Z coordinate
    sorted_planes = sorted(plane_data_list, key=lambda p: p[2][2])

    # Right side position
    label_x = x_min - 5
    label_y_right = y_max + 60  # Far right
    label_y_middle = y_max + 30  # Middle column

    # Z spacing
    z_span = z_max - z_min
    z_label_start = z_min - z_span * 0.1
    z_label_end = z_max + z_span * 0.1
    z_spacing = (z_label_end - z_label_start) / max(1, planes_per_column - 1)

    for i, (plane_mesh, plane_idx, centroid) in enumerate(sorted_planes):
        # Determine which column (alternating)
        column = i % 2
        row = i // 2

        if column == 0:
            # Left column (closer to geometry)
            label_y = label_y_middle
        else:
            # Right column (farther from geometry)
            label_y = label_y_right

        label_z = z_label_start + row * z_spacing
        label_pos = np.array([label_x, label_y, label_z])

        # Leader line
        ax.plot(
            [centroid[0], label_pos[0]],
            [centroid[1], label_pos[1]],
            [centroid[2], label_pos[2]],
            'k-', linewidth=0.8, alpha=0.5
        )

        # Text label - LARGE
        ax.text(
            label_pos[0], label_pos[1], label_pos[2],
            str(plane_idx),
            fontsize=20,
            color='black',
            fontweight='bold',
            ha='left',
            va='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.8)
        )

    # Set axis limits with space for labels
    padding = 5
    label_space_y = 80

    ax.set_xlim(x_min - padding, x_max + padding)
    ax.set_ylim(y_min - padding, y_max + label_space_y)
    ax.set_zlim(z_min - padding, z_max + padding)

    # Labels
    ax.set_xlabel('X (mm)', fontsize=16, labelpad=15)
    ax.set_ylabel('Y (mm)', fontsize=16, labelpad=15)
    ax.set_zlabel('Z (mm)', fontsize=16, labelpad=15)

    # ORTHOGRAPHIC ISOMETRIC view (true isometric: 35.264° elevation, 45° azimuth)
    ax.view_init(elev=35.264, azim=45)
    ax.set_proj_type('ortho')  # Orthographic projection

    # Remove grid and panes
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('lightgray')
    ax.yaxis.pane.set_edgecolor('lightgray')
    ax.zaxis.pane.set_edgecolor('lightgray')
    ax.set_facecolor('white')

    # Title
    ax.set_title(
        f'Plane Index Reference - {partition_display_name}\n({len(plane_indices)} planes)',
        fontsize=22, fontweight='bold', pad=20
    )

    plt.tight_layout()

    # Output
    output_path = project_root / f"{partition_name}_planes_reference_iso.png"
    plt.savefig(str(output_path), dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    file_size = output_path.stat().st_size / (1024 * 1024)
    print(f"\n{'=' * 70}")
    print(f"✓ COMPLETE")
    print(f"  File: {output_path.name}")
    print(f"  Size: {file_size:.1f} MB")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
