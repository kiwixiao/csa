#!/usr/bin/env python3
"""
Diagnose Plane Intersection Issues

Load a combined planes STL and analyze:
1. Are planes properly oriented (perpendicular to centerline)?
2. Are planes intersecting each other?
3. What is the plane spacing?
"""

import sys
import numpy as np
import trimesh
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D


def analyze_combined_planes_stl(stl_path: str, vtk_path: str = None):
    """
    Analyze a combined planes STL for intersection issues

    Args:
        stl_path: Path to *-Planes-All.stl file
        vtk_path: Optional path to corresponding VTK centerline
    """
    print("="*80)
    print(f"Analyzing: {stl_path}")
    print("="*80)

    # Load STL
    mesh = trimesh.load_mesh(stl_path)
    print(f"\nMesh loaded:")
    print(f"  Vertices: {len(mesh.vertices)}")
    print(f"  Faces: {len(mesh.faces)}")
    print(f"  Components: {len(mesh.split())}")

    # Split into individual planes (separate components)
    components = mesh.split()
    print(f"\n{len(components)} individual planes detected")

    if len(components) == 0:
        print("Error: No components found")
        return

    # Analyze each plane
    plane_data = []

    for i, plane_mesh in enumerate(components):
        if len(plane_mesh.vertices) < 3:
            continue

        # Get plane centroid
        centroid = plane_mesh.centroid

        # Compute plane normal using PCA on vertices
        vertices = plane_mesh.vertices
        centered = vertices - vertices.mean(axis=0)

        # SVD to find plane normal (smallest singular vector)
        U, S, Vt = np.linalg.svd(centered, full_matrices=False)
        plane_normal = Vt[-1]  # Last row = smallest component = normal

        # Get plane bounds
        z_min = vertices[:, 2].min()
        z_max = vertices[:, 2].max()
        z_thickness = z_max - z_min

        plane_data.append({
            'index': i,
            'centroid': centroid,
            'normal': plane_normal,
            'z_min': z_min,
            'z_max': z_max,
            'z_thickness': z_thickness,
            'n_vertices': len(vertices),
            'area': plane_mesh.area
        })

    # Sort by Z coordinate of centroid
    plane_data.sort(key=lambda p: p['centroid'][2])

    # Check for Z-overlaps (intersection indicator)
    print(f"\n{'='*80}")
    print("Plane Overlap Analysis (sorted by Z)")
    print(f"{'='*80}")
    print(f"{'Index':<6} {'Z_centroid':<12} {'Z_range':<20} {'Z_thick':<10} {'Normal_Z':<10} {'Overlap?':<10}")
    print("-"*80)

    overlaps = []
    for i, plane in enumerate(plane_data):
        z_cent = plane['centroid'][2]
        z_range = f"[{plane['z_min']:.1f}, {plane['z_max']:.1f}]"
        z_thick = plane['z_thickness']
        normal_z = plane['normal'][2]

        # Check overlap with next plane
        overlap = ""
        if i < len(plane_data) - 1:
            next_plane = plane_data[i + 1]
            # If current plane's z_max > next plane's z_min, they overlap
            if plane['z_max'] > next_plane['z_min']:
                overlap = "YES"
                overlaps.append((i, i+1, plane['z_max'] - next_plane['z_min']))

        print(f"{plane['index']:<6} {z_cent:<12.2f} {z_range:<20} {z_thick:<10.2f} {normal_z:<10.3f} {overlap:<10}")

    # Report overlaps
    print(f"\n{'='*80}")
    print(f"OVERLAP SUMMARY: {len(overlaps)} pairs of overlapping planes")
    print(f"{'='*80}")

    if len(overlaps) > 0:
        print("\nOverlapping plane pairs:")
        for i, j, overlap_amount in overlaps[:10]:  # Show first 10
            print(f"  Plane {plane_data[i]['index']} and Plane {plane_data[j]['index']}: {overlap_amount:.2f} mm Z-overlap")
        if len(overlaps) > 10:
            print(f"  ... and {len(overlaps) - 10} more overlaps")

    # Analyze plane spacing
    print(f"\n{'='*80}")
    print("Plane Spacing Analysis")
    print(f"{'='*80}")

    spacings = []
    for i in range(len(plane_data) - 1):
        spacing = plane_data[i+1]['centroid'][2] - plane_data[i]['centroid'][2]
        spacings.append(spacing)

    if spacings:
        print(f"Mean spacing: {np.mean(spacings):.2f} mm")
        print(f"Std spacing: {np.std(spacings):.2f} mm")
        print(f"Min spacing: {np.min(spacings):.2f} mm")
        print(f"Max spacing: {np.max(spacings):.2f} mm")

    # Analyze plane normals
    print(f"\n{'='*80}")
    print("Plane Normal Analysis")
    print(f"{'='*80}")

    normals = np.array([p['normal'] for p in plane_data])
    z_components = normals[:, 2]

    print(f"Normal Z-component range: [{z_components.min():.3f}, {z_components.max():.3f}]")
    print(f"Mean |Normal_Z|: {np.abs(z_components).mean():.3f}")

    # Check if normals are roughly aligned with Z-axis
    # For "descending" planes, normals should have significant Z-component
    non_z_aligned = np.sum(np.abs(z_components) < 0.5)
    print(f"Planes with |Normal_Z| < 0.5: {non_z_aligned} / {len(plane_data)}")

    # Visualize
    visualize_plane_orientations(plane_data, stl_path)

    return plane_data, overlaps


def visualize_plane_orientations(plane_data, stl_path):
    """
    Visualize plane centroids and normals in 3D

    Args:
        plane_data: List of plane data dictionaries
        stl_path: Original STL path for output naming
    """
    fig = plt.figure(figsize=(16, 6))

    # 3D plot: plane positions and normals
    ax1 = fig.add_subplot(131, projection='3d')

    centroids = np.array([p['centroid'] for p in plane_data])
    normals = np.array([p['normal'] for p in plane_data])

    # Plot centroids
    ax1.scatter(centroids[:, 0], centroids[:, 1], centroids[:, 2],
                c=centroids[:, 2], cmap='viridis', s=30, alpha=0.6)

    # Plot normals as arrows (every 5th plane for clarity)
    for i in range(0, len(plane_data), 5):
        c = plane_data[i]['centroid']
        n = plane_data[i]['normal'] * 5  # Scale for visibility
        ax1.quiver(c[0], c[1], c[2], n[0], n[1], n[2],
                   color='red', alpha=0.5, arrow_length_ratio=0.3)

    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')
    ax1.set_zlabel('Z (mm)')
    ax1.set_title('Plane Positions & Normals', fontweight='bold')

    # 2D plot: Z vs Normal_Z
    ax2 = fig.add_subplot(132)
    z_coords = centroids[:, 2]
    normal_z = normals[:, 2]
    ax2.scatter(z_coords, normal_z, alpha=0.6, s=20)
    ax2.axhline(y=0, color='r', linestyle='--', label='Normal_Z = 0')
    ax2.set_xlabel('Plane Z Coordinate (mm)', fontweight='bold')
    ax2.set_ylabel('Normal Z Component', fontweight='bold')
    ax2.set_title('Normal Orientation vs Position', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 3D plot: Plane thickness visualization
    ax3 = fig.add_subplot(133)
    thicknesses = [p['z_thickness'] for p in plane_data]
    indices = [p['index'] for p in plane_data]
    ax3.scatter(indices, thicknesses, alpha=0.6, s=20, c=z_coords, cmap='viridis')
    ax3.set_xlabel('Plane Index', fontweight='bold')
    ax3.set_ylabel('Z Thickness (mm)', fontweight='bold')
    ax3.set_title('Plane Thickness (Z-extent)', fontweight='bold')
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save
    output_path = Path(stl_path).parent / f"{Path(stl_path).stem}_diagnosis.png"
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"\nSaved visualization: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python diagnose_plane_intersections.py <path-to-Planes-All.stl>")
        print("\nExample:")
        print("  python diagnose_plane_intersections.py LeftNoseDecendingSlicedSTLs/out_002000-Planes-All.stl")
        sys.exit(1)

    stl_path = sys.argv[1]

    if not Path(stl_path).exists():
        print(f"Error: File not found: {stl_path}")
        sys.exit(1)

    analyze_combined_planes_stl(stl_path)
