#!/usr/bin/env python3
"""
Create interactive 3D plane reference using Plotly
Features:
- Hover over planes to see their index
- Rotate, zoom, pan freely
- Toggle visibility of different components
- Clean, professional rendering
"""

import numpy as np
import plotly.graph_objects as go
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


def mesh_to_plotly(mesh, opacity=1.0, color='orange', name='Mesh', showlegend=True, hovertext=None):
    """Convert trimesh to plotly mesh3d"""
    vertices = mesh.vertices
    faces = mesh.faces

    # Plotly uses i, j, k for triangle vertices
    i = faces[:, 0]
    j = faces[:, 1]
    k = faces[:, 2]

    return go.Mesh3d(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        i=i,
        j=j,
        k=k,
        opacity=opacity,
        color=color,
        name=name,
        showlegend=showlegend,
        hovertext=hovertext,
        hoverinfo='text' if hovertext else 'skip',
        flatshading=False,
        lighting=dict(
            ambient=0.5,
            diffuse=0.8,
            specular=0.2,
            roughness=0.5
        )
    )


def main():
    parser = argparse.ArgumentParser(
        description='Create interactive Plotly plane reference figure'
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
    print(f"Creating Interactive Plotly Reference - {partition_display_name}")
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

    print(f"\nCreating interactive Plotly figure...")

    # Create figure
    fig = go.Figure()

    # Add airway mesh (semi-transparent)
    if airway_mesh is not None:
        # Subsample airway for performance (every 10th face)
        airway_faces_sub = airway_mesh.faces[::10]
        airway_sub = trimesh.Trimesh(vertices=airway_mesh.vertices, faces=airway_faces_sub)
        airway_trace = mesh_to_plotly(
            airway_sub,
            opacity=0.15,
            color='lightblue',
            name='Airway',
            showlegend=True
        )
        fig.add_trace(airway_trace)
        print(f"  ✓ Airway mesh added")

    # Add centerline
    if centerline_points is not None and len(centerline_points) > 0:
        fig.add_trace(go.Scatter3d(
            x=centerline_points[:, 0],
            y=centerline_points[:, 1],
            z=centerline_points[:, 2],
            mode='lines',
            line=dict(color='red', width=6),
            name='Centerline',
            hoverinfo='skip'
        ))
        print(f"  ✓ Centerline added")

    # Add planes with hover labels
    print(f"  Adding {len(plane_data_list)} planes with hover labels...")
    for i, (plane_mesh, plane_idx, centroid) in enumerate(plane_data_list):
        # Create hover text with plane index
        hover_text = f"Plane {plane_idx}"

        plane_trace = mesh_to_plotly(
            plane_mesh,
            opacity=0.7,
            color='orange',
            name=f'Plane {plane_idx}',
            showlegend=False,  # Don't show individual planes in legend
            hovertext=hover_text
        )
        fig.add_trace(plane_trace)

        if (i + 1) % 25 == 0:
            print(f"    Progress: {i+1}/{len(plane_data_list)} planes")

    print(f"  ✓ All {len(plane_data_list)} planes added")

    # Calculate scene bounds
    if airway_mesh is not None:
        bounds = airway_mesh.bounds
        x_range = [bounds[0][0], bounds[1][0]]
        y_range = [bounds[0][1], bounds[1][1]]
        z_range = [bounds[0][2], bounds[1][2]]
    else:
        all_verts = np.vstack([p[0].vertices for p in plane_data_list])
        x_range = [all_verts[:, 0].min(), all_verts[:, 0].max()]
        y_range = [all_verts[:, 1].min(), all_verts[:, 1].max()]
        z_range = [all_verts[:, 2].min(), all_verts[:, 2].max()]

    # Update layout for isometric-like view
    fig.update_layout(
        title=dict(
            text=f'Interactive Plane Reference - {partition_display_name}<br><sub>{len(plane_indices)} planes | Hover over planes to see index</sub>',
            x=0.5,
            xanchor='center',
            font=dict(size=20)
        ),
        scene=dict(
            xaxis=dict(title='X (mm)', range=x_range),
            yaxis=dict(title='Y (mm)', range=y_range),
            zaxis=dict(title='Z (mm)', range=z_range),
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5),  # Isometric-like view
                projection=dict(type='orthographic')
            ),
            aspectmode='data',
            bgcolor='white'
        ),
        width=1600,
        height=1200,
        showlegend=True,
        legend=dict(
            x=0.02,
            y=0.98,
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='gray',
            borderwidth=1
        ),
        hovermode='closest',
        paper_bgcolor='white',
        plot_bgcolor='white'
    )

    # Save HTML
    output_html = project_root / f"{partition_name}_planes_reference_interactive.html"
    fig.write_html(str(output_html))
    html_size = output_html.stat().st_size / (1024 * 1024)

    print(f"\n{'=' * 70}")
    print(f"✓ COMPLETE")
    print(f"  Partition: {partition_display_name}")
    print(f"  Planes: {len(plane_indices)}")
    print(f"  Output: {output_html.name}")
    print(f"  Size: {html_size:.1f} MB")
    print(f"\nUsage:")
    print(f"  - Open {output_html.name} in web browser")
    print(f"  - Hover over any plane to see its index")
    print(f"  - Click and drag to rotate")
    print(f"  - Scroll to zoom")
    print(f"  - Click legend items to toggle visibility")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
