#!/usr/bin/env python3
"""
Analyze frame 14 mesh boundaries to understand the plane 76 anomaly
"""

import trimesh
import numpy as np
import pandas as pd
from pathlib import Path

print("="*80)
print("Frame 14 Mesh Boundary Analysis")
print("="*80)

# Load frame 14 mesh
mesh_file = Path("LeftNoseDecending/FFD/stl/out_001400.stl")
print(f"\nLoading mesh: {mesh_file}")
mesh = trimesh.load_mesh(str(mesh_file))

print(f"  Vertices: {len(mesh.vertices)}")
print(f"  Faces: {len(mesh.faces)}")
print(f"  Is watertight: {mesh.is_watertight}")
print(f"  Is volume: {mesh.is_volume}")

# Get boundary edges
print("\n" + "="*80)
print("Boundary Edge Analysis")
print("="*80)

# In trimesh, boundary edges are edges that appear only once (not shared by 2 faces)
edges = mesh.edges_unique
edges_face_count = mesh.edges_unique_inverse

# Count how many faces each edge belongs to
unique_edges, counts = np.unique(edges_face_count, return_counts=True)
edge_face_counts = dict(zip(unique_edges, counts))

# Boundary edges appear in only 1 face
boundary_edge_mask = np.array([edge_face_counts.get(i, 0) == 1 for i in range(len(edges))])
boundary_edges = edges[boundary_edge_mask]

print(f"\nTotal unique edges: {len(edges)}")
print(f"Boundary edges (appear in only 1 face): {len(boundary_edges)}")

# Get vertices of boundary edges
boundary_vertices_indices = np.unique(boundary_edges.flatten())
boundary_vertices = mesh.vertices[boundary_vertices_indices]

print(f"Vertices on boundary: {len(boundary_vertices)}")

# Analyze Z-coordinates of boundary
print("\n" + "="*80)
print("Boundary Z-Coordinate Range")
print("="*80)

z_min = boundary_vertices[:, 2].min()
z_max = boundary_vertices[:, 2].max()
z_mean = boundary_vertices[:, 2].mean()
z_std = boundary_vertices[:, 2].std()

print(f"  Boundary Z min: {z_min:.2f} mm")
print(f"  Boundary Z max: {z_max:.2f} mm")
print(f"  Boundary Z mean: {z_mean:.2f} mm")
print(f"  Boundary Z std: {z_std:.2f} mm")

# Compare with overall mesh Z-range
mesh_z_min = mesh.vertices[:, 2].min()
mesh_z_max = mesh.vertices[:, 2].max()
print(f"\n  Overall mesh Z min: {mesh_z_min:.2f} mm")
print(f"  Overall mesh Z max: {mesh_z_max:.2f} mm")

# Load plane information from CSV to find plane 76's position
print("\n" + "="*80)
print("Plane 76 Position Analysis")
print("="*80)

csv_file = "LeftNoseDecendingSlicedSTLs/out_001400-Data.csv"
df = pd.read_csv(csv_file)

# Find plane 76
plane_76 = df[df['plane_index'] == 76].iloc[0]
plane_z = plane_76['plane_z']

print(f"Plane 76 Z-position: {plane_z:.2f} mm")
print(f"  CSA: {plane_76['area_mm2']:.2f} mm²")

# Check if plane 76 is near boundary region
distance_to_boundary_min = abs(plane_z - z_min)
distance_to_boundary_max = abs(plane_z - z_max)

print(f"\n  Distance to nearest boundary Z-min: {distance_to_boundary_min:.2f} mm")
print(f"  Distance to nearest boundary Z-max: {distance_to_boundary_max:.2f} mm")

# Find boundary vertices within ±5mm of plane 76
tolerance = 5.0
nearby_boundary_vertices = boundary_vertices[
    np.abs(boundary_vertices[:, 2] - plane_z) < tolerance
]

print(f"\n  Boundary vertices within ±{tolerance}mm of plane 76 Z: {len(nearby_boundary_vertices)}")

if len(nearby_boundary_vertices) > 0:
    print(f"  ⚠️  PLANE 76 IS NEAR MESH BOUNDARIES!")
    print(f"  This explains why the cross-section is incomplete.")

    # Show Z-distribution of nearby boundary points
    nearby_z = nearby_boundary_vertices[:, 2]
    print(f"\n  Nearby boundary Z-range: [{nearby_z.min():.2f}, {nearby_z.max():.2f}] mm")
else:
    print(f"  ✓ Plane 76 is NOT near mesh boundaries")
    print(f"  The incomplete cross-section may have another cause")

# Analyze plane distribution
print("\n" + "="*80)
print("All Planes Z-Position Distribution")
print("="*80)

plane_z_values = df['plane_z'].values
print(f"  Plane Z range: [{plane_z_values.min():.2f}, {plane_z_values.max():.2f}] mm")
print(f"  Number of planes: {len(plane_z_values)}")

# Check which planes might be affected by boundaries
planes_near_boundary_min = np.sum(np.abs(plane_z_values - z_min) < tolerance)
planes_near_boundary_max = np.sum(np.abs(plane_z_values - z_max) < tolerance)

print(f"\n  Planes within ±{tolerance}mm of Z-min boundary: {planes_near_boundary_min}")
print(f"  Planes within ±{tolerance}mm of Z-max boundary: {planes_near_boundary_max}")

# Check if other planes near plane 76 have issues
nearby_planes = df[(df['plane_index'] >= 70) & (df['plane_index'] <= 82)].copy()
nearby_planes = nearby_planes.sort_values('plane_index')

print("\n" + "="*80)
print("Planes Near Index 76 (70-82)")
print("="*80)
print(f"{'Index':>6} {'Z-pos':>10} {'CSA':>10} {'Near Boundary?':>20}")
print("-"*50)

for _, plane in nearby_planes.iterrows():
    idx = int(plane['plane_index'])
    z = plane['plane_z']
    area = plane['area_mm2']

    # Check if near boundary
    boundary_count = len(boundary_vertices[np.abs(boundary_vertices[:, 2] - z) < tolerance])
    near = "YES" if boundary_count > 0 else "no"

    marker = "  ⚠️" if idx == 76 else ""
    print(f"{idx:>6} {z:>10.2f} {area:>10.2f} {near:>20}{marker}")

print("\n" + "="*80)
print("Summary")
print("="*80)

if len(nearby_boundary_vertices) > 0:
    print("✓ EXPLANATION FOUND:")
    print(f"  Plane 76 at Z={plane_z:.2f}mm intersects the mesh boundary region")
    print(f"  {len(nearby_boundary_vertices)} boundary vertices found within ±{tolerance}mm")
    print(f"  When trimesh.section() slices at this location, it hits open edges")
    print(f"  This produces an incomplete/open contour instead of a closed loop")
    print(f"  Result: Incomplete STL with only {plane_76['area_mm2']:.2f} mm² instead of expected area")
else:
    print("⚠️  UNEXPECTED:")
    print(f"  Plane 76 is NOT near the boundary regions")
    print(f"  The mesh is non-watertight but boundary is elsewhere")
    print(f"  Further investigation needed to explain the incomplete cross-section")

print("="*80)
