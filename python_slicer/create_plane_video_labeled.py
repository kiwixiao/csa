#!/usr/bin/env python3
"""
Create MP4 video showing plane motion across breathing cycle with clean annotation labels
Labels are positioned NEXT TO planes (not overlapping), like technical diagram annotations
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import glob
import os
import sys
import argparse
from pathlib import Path
import subprocess
import trimesh
import re
import pandas as pd


def format_partition_name(partition_name):
    """Format partition name for display (e.g., 'LeftNoseDecending' -> 'Left Nose Descending')"""
    formatted = re.sub(r'([a-z])([A-Z])', r'\1 \2', partition_name)
    return formatted


def load_plane_mesh(output_dir, time_point, plane_idx):
    """Load a single plane STL file"""
    plane_stl = output_dir / f"{time_point}-Planes-{plane_idx:03d}.stl"

    if plane_stl.exists():
        try:
            return trimesh.load_mesh(str(plane_stl))
        except Exception as e:
            print(f"Warning: Could not load plane {plane_idx}: {e}")
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


def calculate_label_positions_side_layout(plane_data_list, airway_bounds):
    """
    Calculate label positions arranged vertically on RIGHT side of geometry

    Args:
        plane_data_list: List of (plane_mesh, plane_idx, centroid) tuples
        airway_bounds: (x_min, x_max, y_min, y_max, z_min, z_max)

    Returns:
        List of (plane_idx, centroid, label_pos) tuples
    """
    x_min, x_max, y_min, y_max, z_min, z_max = airway_bounds

    # Place labels on the RIGHT side (positive Y direction for sagittal view)
    # Position: 30mm to the right of the geometry
    label_x = x_min  # Keep same X as geometry center
    label_y = y_max + 30  # 30mm to the right of airway

    # Spread labels vertically across Z range with extra padding
    num_labels = len(plane_data_list)
    z_range = z_max - z_min
    z_padding = z_range * 0.2  # 20% padding top and bottom

    z_start = z_min - z_padding
    z_end = z_max + z_padding
    z_spacing = (z_end - z_start) / max(1, num_labels - 1) if num_labels > 1 else 0

    label_positions = []
    for i, (plane_mesh, plane_idx, centroid) in enumerate(plane_data_list):
        label_z = z_start + i * z_spacing
        label_pos = np.array([label_x, label_y, label_z])
        label_positions.append((plane_idx, centroid, label_pos))

    return label_positions


def create_labeled_frame(airway_mesh, plane_data_list, centerline_points,
                         output_path, time_point, partition_display_name):
    """
    Create a single frame with planes and annotation labels

    Args:
        airway_mesh: Trimesh object for airway
        plane_data_list: List of (plane_mesh, plane_idx, centroid) tuples
        output_path: Output PNG path
        time_point: Time point name (e.g., 'out_000000')
        partition_display_name: Formatted partition name
    """
    fig = plt.figure(figsize=(24, 20), dpi=200)
    ax = fig.add_subplot(111, projection='3d')

    # Plot airway mesh (transparent blue background)
    if airway_mesh is not None:
        airway_verts = airway_mesh.vertices
        airway_faces = airway_mesh.faces

        # Subsample for faster rendering
        airway_faces_sub = airway_faces[::4]

        airway_collection = Poly3DCollection(
            airway_verts[airway_faces_sub],
            alpha=0.1,
            facecolors='lightblue',
            edgecolors='none',
            linewidths=0
        )
        ax.add_collection3d(airway_collection)

    # Draw centerline if available
    if centerline_points is not None and len(centerline_points) > 0:
        ax.plot(
            centerline_points[:, 0],
            centerline_points[:, 1],
            centerline_points[:, 2],
            'r-', linewidth=2, alpha=0.8, label='Centerline'
        )

    # Calculate cumulative arc length along centerline for natural spacing
    if centerline_points is not None and len(centerline_points) > 0:
        # Compute arc length at each centerline point
        arc_lengths = np.zeros(len(centerline_points))
        for i in range(1, len(centerline_points)):
            segment_length = np.linalg.norm(centerline_points[i] - centerline_points[i-1])
            arc_lengths[i] = arc_lengths[i-1] + segment_length

        total_arc_length = arc_lengths[-1]

    # Plot planes with labels distributed by arc length along centerline
    for i, (plane_mesh, plane_idx, centroid) in enumerate(plane_data_list):
        # Plot plane surface
        plane_collection = Poly3DCollection(
            plane_mesh.vertices[plane_mesh.faces],
            alpha=0.7,
            facecolors='orange',
            edgecolors='darkgray',
            linewidths=0.3
        )
        ax.add_collection3d(plane_collection)

        # Find closest centerline point to this plane's centroid
        if centerline_points is not None and len(centerline_points) > 0:
            distances = np.linalg.norm(centerline_points - centroid, axis=1)
            closest_idx = np.argmin(distances)
            arc_position = arc_lengths[closest_idx]

            # Calculate offset that increases with arc length position
            # This creates natural spacing following the centerline curve
            arc_ratio = arc_position / total_arc_length if total_arc_length > 0 else 0

            # CLOSER to geometry (reduced Y offset) but MORE SPREAD vertically (Z spacing)
            base_offset_y = 15  # Closer horizontal offset (reduced from 30)
            base_offset_z = -30 + arc_ratio * 100  # Vertical spread: -30 to +70 mm range

            label_offset = np.array([0, base_offset_y, base_offset_z])
        else:
            # Fallback if no centerline
            label_offset = np.array([0, 30, 5])

        label_pos = centroid + label_offset

        # Draw leader line from plane centroid to label
        ax.plot(
            [centroid[0], label_pos[0]],
            [centroid[1], label_pos[1]],
            [centroid[2], label_pos[2]],
            'k-', linewidth=1.0, alpha=0.7
        )

        # Draw text label with LARGER font
        ax.text(
            label_pos[0], label_pos[1], label_pos[2],
            str(plane_idx),
            fontsize=14,
            color='black',
            fontweight='bold',
            ha='left',
            va='center'
        )

    # Set axis limits based on airway mesh with MINIMAL padding to fill frame
    if airway_mesh is not None:
        airway_verts = airway_mesh.vertices
        x_min, x_max = airway_verts[:, 0].min(), airway_verts[:, 0].max()
        y_min, y_max = airway_verts[:, 1].min(), airway_verts[:, 1].max()
        z_min, z_max = airway_verts[:, 2].min(), airway_verts[:, 2].max()

        # Minimal padding to fit geometry to frame - extend Y to accommodate labels
        padding = 3  # Very small padding for tight fit
        label_space = 100  # Extra space on right side for labels

        ax.set_xlim(x_min - padding, x_max + padding)
        ax.set_ylim(y_min - padding, y_max + label_space)  # Extra space for labels
        ax.set_zlim(z_min - padding, z_max + padding)

    # Set labels
    ax.set_xlabel('X (mm)', fontsize=12)
    ax.set_ylabel('Y (mm)', fontsize=12)
    ax.set_zlabel('Z (mm)', fontsize=12)

    # Set SAGITTAL view (side view) with ORTHOGONAL projection
    ax.view_init(elev=0, azim=0)
    ax.set_proj_type('ortho')

    # REMOVE BACKGROUND GRID for cleaner look
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('none')
    ax.yaxis.pane.set_edgecolor('none')
    ax.zaxis.pane.set_edgecolor('none')

    # Add title
    phase = int(time_point.split('_')[-1])
    ax.set_title(
        f'Airway Cross-Sectional Planes (Labeled) - Breathing Phase {phase}\n{partition_display_name}',
        fontsize=13, fontweight='bold'
    )

    # Add annotation
    num_planes = len(plane_data_list)
    plane_text = f'{num_planes} planes' if num_planes else 'Multiple planes'
    ax.text2D(
        0.05, 0.95,
        f'Phase: {phase}/2000\n{plane_text}\nAnnotated',
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7)
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def create_video(frames_dir, output_video, fps=5):
    """Create MP4 video from frames using ffmpeg"""
    frame_pattern = os.path.join(frames_dir, 'frame_%04d.png')

    cmd = [
        'ffmpeg', '-y',
        '-framerate', str(fps),
        '-i', frame_pattern,
        '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        '-crf', '23',
        '-movflags', '+faststart',
        output_video
    ]

    print(f"\nCreating video: {output_video}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"✓ Video created successfully!")
        return True
    else:
        print(f"✗ Error creating video:")
        print(result.stderr)
        return False


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Create labeled plane motion video for breathing cycle'
    )
    parser.add_argument(
        'partition', nargs='?', default=None,
        help='Partition name (e.g., LeftNoseDecending, RightNose)'
    )
    args = parser.parse_args()

    # Auto-detect partition if not provided
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    if args.partition:
        partition_name = args.partition
    else:
        # Auto-detect: look for directories with SlicedSTLs suffix
        sliced_dirs = list(project_root.glob("*SlicedSTLs"))
        if not sliced_dirs:
            print("Error: No SlicedSTLs directories found. Please specify partition name.")
            print("Usage: python create_plane_video_labeled.py <partition_name>")
            sys.exit(1)

        # Use the first one found
        sliced_dirs_sorted = sorted(
            sliced_dirs,
            key=lambda x: (not x.name.startswith("LeftNose"), x.name)
        )
        partition_name = sliced_dirs_sorted[0].name.replace("SlicedSTLs", "")
        print(f"Auto-detected partition: {partition_name}")

    partition_display_name = format_partition_name(partition_name)

    print("=" * 80)
    print(f"Creating LABELED Plane Motion Video - {partition_display_name}")
    print("=" * 80)

    # Setup paths
    output_dir = project_root / f"{partition_name}SlicedSTLs"
    frames_dir = project_root / f"video_frames_{partition_name}_labeled"
    frames_dir.mkdir(exist_ok=True)

    # Load reference airway mesh (first time point)
    partition_dir = project_root / partition_name / "FFD" / "stl"
    airway_stl_files = sorted(glob.glob(str(partition_dir / "*.stl")))

    airway_mesh = None
    if airway_stl_files:
        print(f"\nLoading reference airway mesh: {Path(airway_stl_files[0]).name}")
        try:
            airway_mesh = trimesh.load_mesh(airway_stl_files[0])
            print(f"  Loaded: {len(airway_mesh.vertices)} vertices, {len(airway_mesh.faces)} faces")
        except Exception as e:
            print(f"  Warning: Could not load airway mesh: {e}")
            print("  Will only show planes")

    # Find all time points
    plane_stl_files = sorted(glob.glob(str(output_dir / "*-Planes-All.stl")))

    if not plane_stl_files:
        print(f"Error: No plane STL files found in {output_dir}")
        return

    print(f"\nFound {len(plane_stl_files)} time points")
    print(f"Creating frames with annotation labels...")

    # Create frames for each time point
    for i, plane_stl in enumerate(plane_stl_files):
        time_point = Path(plane_stl).stem.replace('-Planes-All', '')

        print(f"[{i+1}/{len(plane_stl_files)}] {time_point}...", end=' ')

        # Get plane indices for this time point
        plane_indices = get_plane_indices(output_dir, time_point)

        # Load centerline from CSV
        centerline_points = load_centerline_from_csv(output_dir, time_point)

        # Load individual planes and calculate centroids
        plane_data_list = []
        for plane_idx in plane_indices:
            plane_mesh = load_plane_mesh(output_dir, time_point, plane_idx)
            if plane_mesh is not None:
                centroid = plane_mesh.vertices.mean(axis=0)
                plane_data_list.append((plane_mesh, plane_idx, centroid))

        # Create frame
        frame_path = frames_dir / f"frame_{i:04d}.png"
        create_labeled_frame(
            airway_mesh,
            plane_data_list,
            centerline_points,
            str(frame_path),
            time_point,
            partition_display_name
        )

        print(f"✓")

    # Create video with "_labeled" suffix
    output_video = project_root / f"{partition_name}_plane_motion_breathing_cycle_labeled.mp4"
    success = create_video(str(frames_dir), str(output_video), fps=5)

    if success:
        video_size = output_video.stat().st_size / (1024 * 1024)
        print(f"\n{'=' * 80}")
        print(f"✓ LABELED VIDEO COMPLETE!")
        print(f"  File: {output_video.name}")
        print(f"  Size: {video_size:.1f} MB")
        print(f"  Duration: {len(plane_stl_files) / 5:.1f} seconds (5 fps)")
        print(f"  Frames: {len(plane_stl_files)}")
        print(f"{'=' * 80}")

        # Clean up frames
        print(f"\nCleaning up temporary frames...")
        import shutil
        shutil.rmtree(frames_dir)
        print(f"✓ Done!")


if __name__ == "__main__":
    main()
