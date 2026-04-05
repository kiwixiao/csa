#!/usr/bin/env python3
"""
Create interactive 3D plane reference figure using PyVista
Produces both a high-quality static PNG and an interactive HTML viewer
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
        description='Create interactive PyVista plane reference figure'
    )
    parser.add_argument(
        'partition', nargs='?', default=None,
        help='Partition name (e.g., LeftNoseDecending, RightNose)'
    )
    parser.add_argument(
        '--interactive', action='store_true',
        help='Create interactive HTML output'
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
    print(f"Creating PyVista Reference Figure - {partition_display_name}")
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

    print(f"\nCreating PyVista scene...")

    # Create PyVista plotter with VERY HIGH RESOLUTION
    plotter = pv.Plotter(off_screen=True, window_size=(6000, 5000))

    # Add airway mesh (semi-transparent for context)
    if airway_mesh is not None:
        airway_pv = pv.read(airway_stl_files[0])
        plotter.add_mesh(
            airway_pv,
            opacity=0.15,
            color='lightblue',
            smooth_shading=True,
            show_edges=False
        )
        print(f"  ✓ Airway mesh added (15% opacity)")

    # Add centerline (bright and visible)
    if centerline_points is not None and len(centerline_points) > 0:
        centerline_pv = pv.lines_from_points(centerline_points)
        plotter.add_mesh(
            centerline_pv,
            color='red',
            line_width=8,
            opacity=1.0,
            label='Centerline'
        )
        print(f"  ✓ Centerline added")

    # Add planes with smart labeling
    print(f"\nAdding {len(plane_data_list)} planes with labels...")

    # Calculate bounds for label positioning
    if airway_mesh is not None:
        bounds = airway_mesh.bounds
        x_min, y_min, z_min = bounds[0]
        x_max, y_max, z_max = bounds[1]
    else:
        all_verts = np.vstack([p[0].vertices for p in plane_data_list])
        x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
        y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
        z_min, z_max = all_verts[:, 2].min(), all_verts[:, 2].max()

    # Sort planes by Z coordinate for organized labeling
    sorted_planes = sorted(plane_data_list, key=lambda p: p[2][2])

    # Label positioning strategy: THREE columns on the right side
    num_planes = len(sorted_planes)
    planes_per_column = (num_planes + 2) // 3  # Split into 3 columns

    # Y positions for three columns (right side of geometry)
    y_offset_base = y_max - y_min
    label_y_positions = [
        y_max + y_offset_base * 0.3,   # Near column
        y_max + y_offset_base * 0.5,   # Middle column
        y_max + y_offset_base * 0.7    # Far column
    ]

    label_x = x_min

    # Z spacing for labels
    z_span = z_max - z_min
    z_label_start = z_min - z_span * 0.15
    z_label_end = z_max + z_span * 0.15
    z_spacing = (z_label_end - z_label_start) / max(1, planes_per_column - 1)

    # Collect all label positions and texts
    label_points = []
    label_texts = []
    leader_lines = []

    for i, (plane_mesh, plane_idx, centroid) in enumerate(sorted_planes):
        # Add plane surface
        faces_pv = np.column_stack((
            np.full(len(plane_mesh.faces), 3),
            plane_mesh.faces
        )).flatten()
        plane_pv = pv.PolyData(plane_mesh.vertices, faces_pv)
        plotter.add_mesh(
            plane_pv,
            color='orange',
            opacity=0.7,
            smooth_shading=True,
            show_edges=False
        )

        # Determine which column (distribute evenly)
        column = i % 3
        row = i // 3

        label_y = label_y_positions[column]
        label_z = z_label_start + row * z_spacing
        label_pos = np.array([label_x, label_y, label_z])

        # Store label data
        label_points.append(label_pos)
        label_texts.append(str(plane_idx))

        # Create leader line
        line_points = np.array([centroid, label_pos])
        line_pv = pv.lines_from_points(line_points)
        plotter.add_mesh(
            line_pv,
            color='black',
            line_width=2,
            opacity=0.6
        )

    # Add all labels at once with PyVista's optimized label rendering
    if label_points:
        label_points_array = np.array(label_points)
        plotter.add_point_labels(
            label_points_array,
            label_texts,
            font_size=36,
            text_color='black',
            font_family='arial',
            point_size=0.1,
            render_points_as_spheres=False,
            always_visible=True,
            shape_opacity=0.9,
            shape_color='white',
            margin=3
        )

    print(f"  ✓ All planes and labels added")

    # Set isometric view with orthographic projection
    plotter.camera.SetParallelProjection(True)
    plotter.view_isometric()

    # Adjust camera position for better viewing angle
    plotter.camera.elevation = 35.264
    plotter.camera.azimuth = 45

    # White background for clarity
    plotter.set_background('white')

    # Add axes for reference
    plotter.add_axes(
        xlabel='X (mm)',
        ylabel='Y (mm)',
        zlabel='Z (mm)',
        line_width=3,
        labels_off=False
    )

    print(f"\nRendering outputs...")

    # Save high-resolution static image
    output_png = project_root / f"{partition_name}_planes_reference_pyvista.png"
    plotter.screenshot(str(output_png), return_img=False)
    png_size = output_png.stat().st_size / (1024 * 1024)
    print(f"  ✓ PNG saved: {output_png.name} ({png_size:.1f} MB)")

    # Save interactive HTML if requested
    if args.interactive:
        output_html = project_root / f"{partition_name}_planes_reference_interactive.html"
        plotter.export_html(str(output_html))
        html_size = output_html.stat().st_size / (1024 * 1024)
        print(f"  ✓ HTML saved: {output_html.name} ({html_size:.1f} MB)")
        print(f"    Open in browser to interact with 3D scene")

    plotter.close()

    print(f"\n{'=' * 70}")
    print(f"✓ COMPLETE")
    print(f"  Partition: {partition_display_name}")
    print(f"  Planes labeled: {len(plane_indices)}")
    print(f"  Output: {output_png.name}")
    if args.interactive:
        print(f"  Interactive: {output_html.name}")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
