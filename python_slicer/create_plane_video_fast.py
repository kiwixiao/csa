#!/usr/bin/env python3
"""
Create MP4 video showing plane motion across breathing cycle (Fast version using trimesh)
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


def format_partition_name(partition_name):
    """
    Format partition name for display (e.g., 'LeftNoseDecending' -> 'Left Nose Descending')
    """
    # Insert spaces before capital letters
    formatted = re.sub(r'([a-z])([A-Z])', r'\1 \2', partition_name)
    return formatted


def create_frame_fast(airway_mesh, plane_mesh, time_point, output_path, view_angle=0, num_planes=None, partition_display_name="Airway Partition"):
    """Create a single frame showing airway and planes using trimesh"""
    fig = plt.figure(figsize=(8, 8), dpi=100)  # 800x800 pixels (both divisible by 2)
    ax = fig.add_subplot(111, projection='3d')

    # Plot airway mesh (transparent)
    if airway_mesh is not None:
        airway_verts = airway_mesh.vertices
        airway_faces = airway_mesh.faces

        # Subsample for faster rendering (every 4th face)
        airway_faces_sub = airway_faces[::4]

        airway_collection = Poly3DCollection(
            airway_verts[airway_faces_sub],
            alpha=0.1,
            facecolors='lightblue',
            edgecolors='none',
            linewidths=0
        )
        ax.add_collection3d(airway_collection)

    # Plot planes (ALL GRAY COLOR)
    if plane_mesh is not None:
        plane_verts = plane_mesh.vertices
        plane_faces = plane_mesh.faces

        plane_collection = Poly3DCollection(
            plane_verts[plane_faces],
            alpha=0.7,
            facecolors='gray',
            edgecolors='darkgray',
            linewidths=0.3
        )
        ax.add_collection3d(plane_collection)

    # Set axis limits
    all_verts = []
    if airway_mesh is not None:
        all_verts.append(airway_mesh.vertices)
    if plane_mesh is not None:
        all_verts.append(plane_mesh.vertices)

    if all_verts:
        all_verts = np.vstack(all_verts)
        x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
        y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
        z_min, z_max = all_verts[:, 2].min(), all_verts[:, 2].max()

        padding = 10
        ax.set_xlim(x_min - padding, x_max + padding)
        ax.set_ylim(y_min - padding, y_max + padding)
        ax.set_zlim(z_min - padding, z_max + padding)

    # Set labels
    ax.set_xlabel('X (mm)', fontsize=11)
    ax.set_ylabel('Y (mm)', fontsize=11)
    ax.set_zlabel('Z (mm)', fontsize=11)

    # Set SAGITTAL view (side view) with ORTHOGONAL projection
    ax.view_init(elev=0, azim=0)
    ax.set_proj_type('ortho')

    # Add title
    phase = int(time_point.split('_')[-1])
    ax.set_title(f'Airway Cross-Sectional Planes - Breathing Phase {phase}\n{partition_display_name}',
                 fontsize=13, fontweight='bold')

    # Add annotation
    plane_text = f'{num_planes} planes' if num_planes else 'Multiple planes'
    ax.text2D(0.05, 0.95, f'Phase: {phase}/2000\n{plane_text}',
              transform=ax.transAxes,
              fontsize=11, verticalalignment='top',
              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    plt.tight_layout()
    plt.savefig(output_path, dpi=100)
    plt.close()


def create_video(frames_dir, output_video, fps=5):
    """Create MP4 video from frames using ffmpeg - simple and robust"""
    frame_pattern = os.path.join(frames_dir, 'frame_%04d.png')

    # Simple, robust ffmpeg settings - just works everywhere
    # The scale filter ensures dimensions are divisible by 2 (required for H.264)
    cmd = [
        'ffmpeg', '-y',
        '-framerate', str(fps),
        '-i', frame_pattern,
        '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',  # Ensure even dimensions
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
    parser = argparse.ArgumentParser(description='Create plane motion video for breathing cycle')
    parser.add_argument('partition', nargs='?', default=None,
                        help='Partition name (e.g., LeftNoseDecending, RightNose)')
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
            print("Usage: python create_plane_video_fast.py <partition_name>")
            sys.exit(1)

        # Use the first one found (prefer LeftNoseDecending if available)
        sliced_dirs_sorted = sorted(sliced_dirs, key=lambda x: (not x.name.startswith("LeftNose"), x.name))
        partition_name = sliced_dirs_sorted[0].name.replace("SlicedSTLs", "")
        print(f"Auto-detected partition: {partition_name}")

    partition_display_name = format_partition_name(partition_name)

    print("="*80)
    print(f"Creating Plane Motion Video - {partition_display_name}")
    print("="*80)

    # Setup paths dynamically based on partition
    output_dir = project_root / f"{partition_name}SlicedSTLs"
    frames_dir = project_root / f"video_frames_{partition_name}"
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

    # Find all plane STL files
    plane_stl_files = sorted(glob.glob(str(output_dir / "*-Planes-All.stl")))

    if not plane_stl_files:
        print(f"Error: No plane STL files found in {output_dir}")
        return

    # Count individual plane files from one time point to get plane count
    first_timepoint = Path(plane_stl_files[0]).stem.replace('-Planes-All', '')
    individual_plane_files = sorted(glob.glob(str(output_dir / f"{first_timepoint}-Planes-*.stl")))
    # Exclude the combined file
    individual_plane_files = [f for f in individual_plane_files if not f.endswith('-Planes-All.stl')]
    num_planes = len(individual_plane_files)

    print(f"\nFound {len(plane_stl_files)} time points")
    print(f"Each time point has {num_planes} planes")
    print(f"Creating frames...")

    # Create frames for each time point
    for i, plane_stl in enumerate(plane_stl_files):
        time_point = Path(plane_stl).stem.replace('-Planes-All', '')

        print(f"[{i+1}/{len(plane_stl_files)}] {time_point}...", end=' ')

        # Load planes for this time point
        try:
            plane_mesh = trimesh.load_mesh(plane_stl)
            print(f"{len(plane_mesh.faces)} faces", end=' ')
        except Exception as e:
            print(f"ERROR: {e}")
            continue

        # Create frame (no rotation - fixed sagittal view)
        frame_path = frames_dir / f"frame_{i:04d}.png"

        create_frame_fast(airway_mesh, plane_mesh, time_point, str(frame_path), 0, num_planes, partition_display_name)
        print(f"✓")

    # Create video
    output_video = project_root / f"{partition_name}_plane_motion_breathing_cycle.mp4"
    success = create_video(str(frames_dir), str(output_video), fps=5)

    if success:
        video_size = output_video.stat().st_size / (1024 * 1024)
        print(f"\n{'='*80}")
        print(f"✓ VIDEO COMPLETE!")
        print(f"  File: {output_video.name}")
        print(f"  Size: {video_size:.1f} MB")
        print(f"  Duration: {len(plane_stl_files) / 5:.1f} seconds (5 fps)")
        print(f"  Frames: {len(plane_stl_files)}")
        print(f"{'='*80}")

        # Clean up frames
        print(f"\nCleaning up temporary frames...")
        import shutil
        shutil.rmtree(frames_dir)
        print(f"✓ Done!")


if __name__ == "__main__":
    main()
