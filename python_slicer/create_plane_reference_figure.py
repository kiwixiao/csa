#!/usr/bin/env python3
"""
Create a professional plane reference figure using PyVista
Clean, high-quality rendering with proper labels
"""

import numpy as np
import pyvista as pv
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
        description='Create professional reference figure with PyVista'
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
    print(f"Creating Professional Reference Figure - {partition_display_name}")
    print("=" * 70)

    # Setup paths
    output_dir = project_root / f"{partition_name}SlicedSTLs"

    # Load airway mesh (first time point)
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

    print(f"\nCreating figure with PyVista...")

    # Create PyVista plotter (off-screen)
    plotter = pv.Plotter(off_screen=True, window_size=(3000, 2400))

    # Add airway mesh (transparent)
    if airway_stl_files:
        airway_pv = pv.read(airway_stl_files[0])
        plotter.add_mesh(airway_pv, opacity=0.05, color='lightgray', smooth_shading=True)

    # Add centerline
    if centerline_points is not None and len(centerline_points) > 0:
        centerline_pv = pv.lines_from_points(centerline_points)
        plotter.add_mesh(centerline_pv, color='red', line_width=2, opacity=0.6)

    # Add planes (orange)
    for plane_mesh, plane_idx, centroid in plane_data_list:
        # Convert trimesh to PyVista
        faces_pv = np.column_stack((
            np.full(len(plane_mesh.faces), 3),
            plane_mesh.faces
        )).flatten()
        plane_pv = pv.PolyData(plane_mesh.vertices, faces_pv)
        plotter.add_mesh(plane_pv, color='orange', opacity=0.7, smooth_shading=True)

    # Calculate label positions (to the right, sorted by Z)
    if airway_stl_files:
        airway_mesh = trimesh.load_mesh(airway_stl_files[0])
        y_max = airway_mesh.vertices[:, 1].max()
    else:
        all_verts = np.vstack([p[0].vertices for p in plane_data_list])
        y_max = all_verts[:, 1].max()

    label_x = np.mean([p[2][0] for p in plane_data_list])
    label_y = y_max + 30

    # Sort by Z coordinate
    sorted_planes = sorted(plane_data_list, key=lambda p: p[2][2])

    # Add labels with leader lines
    min_spacing = 6.0
    prev_z = None

    for plane_mesh, plane_idx, centroid in sorted_planes:
        # Determine label Z position
        target_z = centroid[2]
        if prev_z is not None and target_z < prev_z + min_spacing:
            target_z = prev_z + min_spacing
        prev_z = target_z

        label_pos = np.array([label_x, label_y, target_z])

        # Leader line
        line_points = np.array([centroid, label_pos])
        line_pv = pv.lines_from_points(line_points)
        plotter.add_mesh(line_pv, color='black', line_width=1, opacity=0.6)

        # Label text
        plotter.add_point_labels(
            [label_pos],
            [str(plane_idx)],
            font_size=14,
            text_color='black',
            point_size=1,
            render_points_as_spheres=False,
            always_visible=True,
            shape_opacity=0
        )

    # Set camera to sagittal view (side view)
    plotter.camera_position = 'yz'
    plotter.camera.SetParallelProjection(True)

    # White background
    plotter.set_background('white')

    # Output file
    output_path = project_root / f"{partition_name}_planes_reference.png"

    # Render and save
    plotter.screenshot(str(output_path), return_img=False)
    plotter.close()

    file_size = output_path.stat().st_size / (1024 * 1024)
    print(f"\n{'=' * 70}")
    print(f"✓ COMPLETE")
    print(f"  File: {output_path.name}")
    print(f"  Size: {file_size:.1f} MB")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
