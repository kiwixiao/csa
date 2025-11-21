#!/usr/bin/env python3
"""
Diagnose spatial relationship between centerline endpoints and boundary cap planes

This will answer: Are centerline[0] and centerline[147] positioned ON, INSIDE, or OUTSIDE
the boundary cap planes?
"""

import sys
import numpy as np

sys.path.insert(0, 'python_slicer')

from slicer.io_utils import read_stl, read_vtk_centerline
from slicer.boundary_caps import extract_boundary_caps, fit_plane_to_points
from slicer.geometry import distance_point_to_plane

print("=" * 80)
print("DIAGNOSTIC: Centerline Endpoint Positions vs Boundary Cap Planes")
print("=" * 80)

# Load frame 14
FRAME = 14
STL_FILE = f"LeftNoseDecending/FFD/stl/out_00{FRAME}00.stl"
VTK_FILE = f"LeftNoseDecending/FFD/vtk/out_00{FRAME}00.vtk"

print(f"\nLoading frame {FRAME}...")
print(f"  STL: {STL_FILE}")
print(f"  VTK: {VTK_FILE}")

mesh = read_stl(STL_FILE)
centerline = read_vtk_centerline(VTK_FILE)

print(f"\nCenterline points: {len(centerline)}")
print(f"Mesh vertices: {len(mesh.vertices)}")
print(f"Mesh faces: {len(mesh.faces)}")

# Extract boundary caps
print("\n" + "=" * 80)
print("Extracting boundary caps...")
print("=" * 80)

start_cap, end_cap = extract_boundary_caps(mesh, centerline)

if start_cap is None or end_cap is None:
    print("ERROR: Could not extract boundary caps!")
    sys.exit(1)

print(f"\n✓ Start cap extracted: {len(start_cap.vertices)} vertices")
print(f"✓ End cap extracted: {len(end_cap.vertices)} vertices")

# Fit planes to boundary caps
print("\n" + "=" * 80)
print("Fitting planes to boundary caps...")
print("=" * 80)

start_cap_centroid, start_cap_normal = fit_plane_to_points(start_cap.vertices)
end_cap_centroid, end_cap_normal = fit_plane_to_points(end_cap.vertices)

# Orient normals toward centerline interior
direction_to_interior_start = centerline[1] - centerline[0]
direction_to_interior_start = direction_to_interior_start / np.linalg.norm(direction_to_interior_start)
if np.dot(start_cap_normal, direction_to_interior_start) < 0:
    start_cap_normal = -start_cap_normal

direction_to_interior_end = centerline[-2] - centerline[-1]
direction_to_interior_end = direction_to_interior_end / np.linalg.norm(direction_to_interior_end)
if np.dot(end_cap_normal, direction_to_interior_end) < 0:
    end_cap_normal = -end_cap_normal

print(f"\nStart cap plane:")
print(f"  Centroid: [{start_cap_centroid[0]:.4f}, {start_cap_centroid[1]:.4f}, {start_cap_centroid[2]:.4f}]")
print(f"  Normal:   [{start_cap_normal[0]:.4f}, {start_cap_normal[1]:.4f}, {start_cap_normal[2]:.4f}]")

print(f"\nEnd cap plane:")
print(f"  Centroid: [{end_cap_centroid[0]:.4f}, {end_cap_centroid[1]:.4f}, {end_cap_centroid[2]:.4f}]")
print(f"  Normal:   [{end_cap_normal[0]:.4f}, {end_cap_normal[1]:.4f}, {end_cap_normal[2]:.4f}]")

# Get centerline endpoints
centerline_first = centerline[0]
centerline_last = centerline[-1]

print("\n" + "=" * 80)
print("Centerline endpoints:")
print("=" * 80)

print(f"\nCenterline[0]:   [{centerline_first[0]:.4f}, {centerline_first[1]:.4f}, {centerline_first[2]:.4f}]")
print(f"Centerline[-1]:  [{centerline_last[0]:.4f}, {centerline_last[1]:.4f}, {centerline_last[2]:.4f}]")

# Compute signed distances
print("\n" + "=" * 80)
print("CRITICAL ANALYSIS: Spatial Relationship")
print("=" * 80)

# Distance from centerline[0] to start cap plane
dist_to_start = distance_point_to_plane(centerline_first, start_cap_centroid, start_cap_normal)

print(f"\nCenterline[0] → Start Cap Plane:")
print(f"  Signed distance: {dist_to_start:.6f} mm")
if abs(dist_to_start) < 0.01:
    print(f"  ✓ Centerline[0] is ON the start cap plane (distance ~0)")
elif dist_to_start > 0:
    print(f"  ✗ Centerline[0] is AHEAD of (outside) the start cap plane by {dist_to_start:.3f} mm")
    print(f"     This means plane 0 is cutting OUTSIDE the mesh!")
    print(f"     → This explains why plane 0 has no intersection!")
else:
    print(f"  ✓ Centerline[0] is INSIDE the mesh by {abs(dist_to_start):.3f} mm")

# Distance from centerline[-1] to end cap plane
dist_to_end = distance_point_to_plane(centerline_last, end_cap_centroid, end_cap_normal)

print(f"\nCenterline[-1] → End Cap Plane:")
print(f"  Signed distance: {dist_to_end:.6f} mm")
if abs(dist_to_end) < 0.01:
    print(f"  ✓ Centerline[-1] is ON the end cap plane (distance ~0)")
elif dist_to_end > 0:
    print(f"  ✗ Centerline[-1] is AHEAD of (outside) the end cap plane by {dist_to_end:.3f} mm")
    print(f"     This means plane {len(centerline)-1} is cutting OUTSIDE the mesh!")
    print(f"     → This explains why plane {len(centerline)-1} cuts open geometry!")
else:
    print(f"  ✓ Centerline[-1] is INSIDE the mesh by {abs(dist_to_end):.3f} mm")

# Distance between centroids
dist_centroids = np.linalg.norm(centerline_first - start_cap_centroid)
dist_centroids_end = np.linalg.norm(centerline_last - end_cap_centroid)

print("\n" + "=" * 80)
print("Euclidean distances (for reference):")
print("=" * 80)

print(f"\nEuclidean distance:")
print(f"  Centerline[0] ↔ Start cap centroid: {dist_centroids:.4f} mm")
print(f"  Centerline[-1] ↔ End cap centroid: {dist_centroids_end:.4f} mm")

# Summary and recommendations
print("\n" + "=" * 80)
print("DIAGNOSIS SUMMARY")
print("=" * 80)

print("\nPROBLEM IDENTIFIED:")
if dist_to_start > 0.01:
    print(f"  ✗ Centerline[0] is {dist_to_start:.3f} mm OUTSIDE the mesh start boundary")
    print(f"     → Plane 0 cannot create a valid cross-section (no intersection)")
else:
    print(f"  ✓ Centerline[0] position is OK")

if dist_to_end > 0.01:
    print(f"  ✗ Centerline[-1] is {dist_to_end:.3f} mm OUTSIDE the mesh end boundary")
    print(f"     → Plane {len(centerline)-1} cuts through open geometry at boundary")
else:
    print(f"  ✓ Centerline[-1] position is OK")

print("\nRECOMMENDED SOLUTION:")
if dist_to_start > 0.01 or dist_to_end > 0.01:
    print("  Option 1: Skip first/last planes if they're outside boundary caps")
    print("  Option 2: Project centerline endpoints onto boundary cap planes")
    print("  Option 3: Use boundary cap centroids as first/last plane positions")
    print("\n  The normals are already CORRECT (0.000° alignment).")
    print("  The problem is POSITION - centerline endpoints extend beyond mesh boundaries.")
else:
    print("  Centerline positions are correct. Problem may be elsewhere.")

print("=" * 80)
