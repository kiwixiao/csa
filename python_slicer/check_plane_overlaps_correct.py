#!/usr/bin/env python3
"""
Correct Plane Overlap Detection

Check if adjacent planes actually intersect by testing if vertices
from one plane cross through the plane equation of another.
"""

import sys
import glob
import numpy as np
import trimesh
from pathlib import Path
import matplotlib.pyplot as plt


def signed_distance_to_plane(points, plane_point, plane_normal):
    """
    Compute signed distance from points to a plane

    Args:
        points: Nx3 array of points
        plane_point: 3D point on the plane
        plane_normal: 3D normal vector (should be unit length)

    Returns:
        N array of signed distances (positive = in direction of normal)
    """
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    # Distance = (point - plane_point) · normal
    return np.dot(points - plane_point, plane_normal)


def check_planes_for_real_overlaps(file_pattern: str):
    """
    Check for actual geometric intersections between adjacent planes

    Args:
        file_pattern: Glob pattern for plane STL files
    """
    plane_files = sorted(glob.glob(file_pattern))

    if not plane_files:
        print(f"Error: No files found matching pattern: {file_pattern}")
        return

    print("="*80)
    print(f"Checking {len(plane_files)} planes for GEOMETRIC overlaps")
    print("="*80)

    # Load all planes
    planes = []
    for plane_file in plane_files:
        try:
            mesh = trimesh.load_mesh(plane_file)
            filename = Path(plane_file).stem
            plane_idx = int(filename.split('-')[-1])

            # Fit plane to vertices
            vertices = mesh.vertices
            centroid = vertices.mean(axis=0)

            # Use PCA to find plane normal
            centered = vertices - centroid
            U, S, Vt = np.linalg.svd(centered, full_matrices=False)
            normal = Vt[-1]  # Last singular vector = plane normal

            planes.append({
                'index': plane_idx,
                'file': plane_file,
                'vertices': vertices,
                'centroid': centroid,
                'normal': normal,
                'mesh': mesh
            })
        except Exception as e:
            print(f"Warning: Failed to load {plane_file}: {e}")

    # Sort by index
    planes.sort(key=lambda p: p['index'])

    print(f"\nLoaded {len(planes)} planes\n")

    # Check for geometric overlaps
    overlaps = []

    for i in range(len(planes) - 1):
        plane_a = planes[i]
        plane_b = planes[i + 1]

        # Test 1: Do vertices of B cross through plane A?
        distances_b_to_a = signed_distance_to_plane(
            plane_b['vertices'],
            plane_a['centroid'],
            plane_a['normal']
        )

        # Test 2: Do vertices of A cross through plane B?
        distances_a_to_b = signed_distance_to_plane(
            plane_a['vertices'],
            plane_b['centroid'],
            plane_b['normal']
        )

        # Overlap if:
        # - Some vertices of B are on the "wrong side" of A
        # - Some vertices of A are on the "wrong side" of B

        # Since planes should be ordered along centerline:
        # All of B should be "ahead" of A in the normal direction

        # Check if B's vertices have mixed signs relative to A
        # (some before, some after the plane)
        b_crosses_a = (np.min(distances_b_to_a) < -0.1 and np.max(distances_b_to_a) > 0.1)

        # Check if A's vertices have mixed signs relative to B
        a_crosses_b = (np.min(distances_a_to_b) < -0.1 and np.max(distances_a_to_b) > 0.1)

        # More direct check: vertices from B that are "behind" plane A
        b_behind_a_count = np.sum(distances_b_to_a < -0.1)
        b_behind_a_pct = 100 * b_behind_a_count / len(distances_b_to_a)

        # Vertices from A that are "ahead" of plane B
        a_ahead_b_count = np.sum(distances_a_to_b > 0.1)
        a_ahead_b_pct = 100 * a_ahead_b_count / len(distances_a_to_b)

        has_overlap = (b_behind_a_pct > 1) or (a_ahead_b_pct > 1)

        if has_overlap:
            overlaps.append({
                'plane_a_idx': plane_a['index'],
                'plane_b_idx': plane_b['index'],
                'b_behind_a_pct': b_behind_a_pct,
                'a_ahead_b_pct': a_ahead_b_pct,
                'b_behind_a_count': b_behind_a_count,
                'a_ahead_b_count': a_ahead_b_count
            })

    # Report
    print(f"{'='*80}")
    print(f"GEOMETRIC OVERLAP DETECTION RESULTS")
    print(f"{'='*80}\n")

    if len(overlaps) == 0:
        print("✓ NO GEOMETRIC OVERLAPS DETECTED!")
        print("  All planes are properly ordered without intersection.")
    else:
        print(f"✗ FOUND {len(overlaps)} OVERLAPPING PLANE PAIRS:\n")

        for overlap in overlaps[:20]:  # Show first 20
            print(f"  Plane {overlap['plane_a_idx']:3d} ←→ Plane {overlap['plane_b_idx']:3d}:")
            if overlap['b_behind_a_pct'] > 1:
                print(f"    {overlap['b_behind_a_count']} vertices ({overlap['b_behind_a_pct']:.1f}%) of Plane {overlap['plane_b_idx']} are BEHIND Plane {overlap['plane_a_idx']}")
            if overlap['a_ahead_b_pct'] > 1:
                print(f"    {overlap['a_ahead_b_count']} vertices ({overlap['a_ahead_b_pct']:.1f}%) of Plane {overlap['plane_a_idx']} are AHEAD of Plane {overlap['plane_b_idx']}")
            print()

        if len(overlaps) > 20:
            print(f"  ... and {len(overlaps) - 20} more overlaps")

    print(f"\n{'='*80}")
    print(f"SUMMARY")
    print(f"{'='*80}")
    print(f"Total adjacent plane pairs: {len(planes) - 1}")
    print(f"Pairs with geometric overlap: {len(overlaps)}")
    print(f"Overlap percentage: {100 * len(overlaps) / (len(planes) - 1):.1f}%")

    return planes, overlaps


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_plane_overlaps_correct.py <file_pattern>")
        print("\nExample:")
        print("  python check_plane_overlaps_correct.py 'LeftNoseDecendingSlicedSTLs/out_002000-Planes-*.stl'")
        sys.exit(1)

    file_pattern = sys.argv[1]
    check_planes_for_real_overlaps(file_pattern)
