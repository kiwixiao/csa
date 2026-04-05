#!/usr/bin/env python3
"""
Create MP4 video showing plane motion across breathing cycle
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import glob
import os
from pathlib import Path
import subprocess

def read_stl_file(stl_path):
    """Read binary STL file and return vertices and faces"""
    with open(stl_path, 'rb') as f:
        # Skip header (80 bytes)
        f.read(80)

        # Read number of triangles (4 bytes)
        n_triangles = np.fromfile(f, dtype=np.uint32, count=1)[0]

        # Read triangles
        # Each triangle: normal (3 floats) + 3 vertices (3 floats each) + attribute (2 bytes)
        vertices = []
        faces = []

        for i in range(n_triangles):
            # Skip normal (12 bytes)
            f.read(12)

            # Read 3 vertices
            v1 = np.fromfile(f, dtype=np.float32, count=3)
            v2 = np.fromfile(f, dtype=np.float32, count=3)
            v3 = np.fromfile(f, dtype=np.float32, count=3)

            # Store vertices
            start_idx = len(vertices)
            vertices.extend([v1, v2, v3])
            faces.append([start_idx, start_idx+1, start_idx+2])

            # Skip attribute (2 bytes)
            f.read(2)

        return np.array(vertices), np.array(faces)


def load_mesh_stl(stl_path):
    """Load airway mesh STL"""
    try:
        vertices, faces = read_stl_file(stl_path)
        return vertices, faces
    except Exception as e:
        print(f"Error loading {stl_path}: {e}")
        return None, None


def create_frame(airway_verts, airway_faces, plane_verts, plane_faces,
                 time_point, output_path, view_angle=0):
    """Create a single frame showing airway and planes"""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot airway mesh (wireframe or transparent surface)
    if airway_verts is not None and len(airway_verts) > 0:
        airway_collection = Poly3DCollection(
            airway_verts[airway_faces],
            alpha=0.15,
            facecolors='lightblue',
            edgecolors='gray',
            linewidths=0.1
        )
        ax.add_collection3d(airway_collection)

    # Plot planes (colored by position along centerline)
    if plane_verts is not None and len(plane_verts) > 0:
        # Color planes by z-coordinate (position along airway)
        z_coords = plane_verts[plane_faces].mean(axis=1)[:, 2]
        z_min, z_max = z_coords.min(), z_coords.max()

        # Normalize colors
        colors = plt.cm.rainbow((z_coords - z_min) / (z_max - z_min + 1e-6))

        plane_collection = Poly3DCollection(
            plane_verts[plane_faces],
            alpha=0.6,
            facecolors=colors,
            edgecolors='black',
            linewidths=0.3
        )
        ax.add_collection3d(plane_collection)

    # Set axis limits
    all_verts = []
    if airway_verts is not None:
        all_verts.append(airway_verts)
    if plane_verts is not None:
        all_verts.append(plane_verts)

    if all_verts:
        all_verts = np.vstack(all_verts)
        x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
        y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
        z_min, z_max = all_verts[:, 2].min(), all_verts[:, 2].max()

        # Add padding
        padding = 10
        ax.set_xlim(x_min - padding, x_max + padding)
        ax.set_ylim(y_min - padding, y_max + padding)
        ax.set_zlim(z_min - padding, z_max + padding)

    # Set labels
    ax.set_xlabel('X (mm)', fontsize=12)
    ax.set_ylabel('Y (mm)', fontsize=12)
    ax.set_zlabel('Z (mm)', fontsize=12)

    # Set view angle (rotate around z-axis)
    ax.view_init(elev=15, azim=45 + view_angle)

    # Add title
    phase = int(time_point.split('_')[-1])
    ax.set_title(f'Airway Slicing Planes - Breathing Phase {phase}\nTime Point: {time_point}',
                 fontsize=14, fontweight='bold')

    # Add text annotation
    ax.text2D(0.05, 0.95, f'Phase: {phase}/2000', transform=ax.transAxes,
              fontsize=12, verticalalignment='top',
              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close()
    print(f"  Saved frame: {output_path}")


def create_video(frames_dir, output_video, fps=5):
    """Create MP4 video from frames using ffmpeg"""
    frame_pattern = os.path.join(frames_dir, 'frame_%04d.png')

    cmd = [
        'ffmpeg', '-y',
        '-framerate', str(fps),
        '-i', frame_pattern,
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        '-crf', '23',
        '-preset', 'medium',
        output_video
    ]

    print(f"\nCreating video: {output_video}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"✓ Video created successfully: {output_video}")
        return True
    else:
        print(f"✗ Error creating video:")
        print(result.stderr)
        return False


def main():
    print("="*80)
    print("Creating Plane Motion Video - Breathing Cycle")
    print("="*80)

    # Setup paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    output_dir = project_root / "LeftNoseDecendingSlicedSTLs"
    frames_dir = project_root / "video_frames"
    frames_dir.mkdir(exist_ok=True)

    # Find airway mesh STL (use first time point as reference)
    partition_dir = project_root / "LeftNoseDecending" / "FFD" / "stl"
    airway_stl_files = sorted(glob.glob(str(partition_dir / "*.stl")))

    if not airway_stl_files:
        print(f"Error: No airway STL files found in {partition_dir}")
        return

    print(f"\nLoading reference airway mesh: {Path(airway_stl_files[0]).name}")
    airway_verts, airway_faces = load_mesh_stl(airway_stl_files[0])

    if airway_verts is None:
        print("Warning: Could not load airway mesh, will only show planes")
    else:
        print(f"  Loaded: {len(airway_verts)} vertices, {len(airway_faces)} faces")

    # Find all plane STL files
    plane_stl_files = sorted(glob.glob(str(output_dir / "*-Planes-All.stl")))

    if not plane_stl_files:
        print(f"Error: No plane STL files found in {output_dir}")
        return

    print(f"\nFound {len(plane_stl_files)} time points")
    print(f"Creating frames with rotating view...")

    # Create frames for each time point
    for i, plane_stl in enumerate(plane_stl_files):
        time_point = Path(plane_stl).stem.replace('-Planes-All', '')

        print(f"\n[{i+1}/{len(plane_stl_files)}] Processing {time_point}...")

        # Load planes for this time point
        plane_verts, plane_faces = load_mesh_stl(plane_stl)

        if plane_verts is None:
            print(f"  Warning: Could not load planes from {plane_stl}")
            continue

        print(f"  Loaded: {len(plane_verts)} vertices, {len(plane_faces)} faces")

        # Create frame with slight rotation to add dynamics
        view_angle = i * (360 / len(plane_stl_files))  # Full rotation over cycle

        frame_path = frames_dir / f"frame_{i:04d}.png"
        create_frame(airway_verts, airway_faces, plane_verts, plane_faces,
                    time_point, str(frame_path), view_angle)

    # Create video
    output_video = project_root / "LeftNoseDecending_plane_motion_breathing_cycle.mp4"
    success = create_video(str(frames_dir), str(output_video), fps=5)

    if success:
        # Get video file size
        video_size = output_video.stat().st_size / (1024 * 1024)  # MB
        print(f"\n{'='*80}")
        print(f"✓ Video creation complete!")
        print(f"  Output: {output_video}")
        print(f"  Size: {video_size:.1f} MB")
        print(f"  Duration: {len(plane_stl_files) / 5:.1f} seconds at 5 fps")
        print(f"{'='*80}")

        # Clean up frames
        print(f"\nCleaning up temporary frames...")
        import shutil
        shutil.rmtree(frames_dir)
        print(f"✓ Cleaned up: {frames_dir}")


if __name__ == "__main__":
    main()
