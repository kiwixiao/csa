"""
Boundary cap extraction for mesh open ends

New feature: Use mesh boundary edges as natural first/last planes
This solves data loss at geometry ends and eliminates open plane issues.
"""

import numpy as np
import trimesh
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from typing import Tuple, Optional
from shapely.geometry import Polygon
from scipy.spatial import Delaunay

from .mesh_intersection import CrossSection


def extract_boundary_edges(mesh: trimesh.Trimesh) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract boundary edges from mesh (edges that appear in only one face)

    Args:
        mesh: Input mesh (can be non-watertight)

    Returns:
        Tuple of (boundary_edges, boundary_vertices)
        - boundary_edges: Mx2 array of edge indices
        - boundary_vertices: Nx3 array of unique boundary vertex positions
    """
    # Get all unique edges and their face counts
    edges = mesh.edges_unique
    edges_face = mesh.edges_unique_inverse

    # Count how many faces each edge belongs to
    unique_edges, counts = np.unique(edges_face, return_counts=True)
    edge_face_counts = dict(zip(unique_edges, counts))

    # Boundary edges appear in only 1 face
    boundary_mask = np.array([edge_face_counts.get(i, 0) == 1 for i in range(len(edges))])
    boundary_edges = edges[boundary_mask]

    if len(boundary_edges) == 0:
        return np.array([]), np.array([])

    # Get unique boundary vertices
    boundary_vertex_indices = np.unique(boundary_edges.flatten())
    boundary_vertices = mesh.vertices[boundary_vertex_indices]

    return boundary_edges, boundary_vertices


def separate_boundary_regions(boundary_edges: np.ndarray,
                               boundary_vertices: np.ndarray,
                               mesh: trimesh.Trimesh,
                               centerline: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Separate boundary vertices into start and end regions using
    connected components and centerline projection

    Strategy:
    1. Find connected components in boundary graph (each opening = 1 component)
    2. For each component, compute centroid
    3. Project centroid onto centerline to find closest centerline point
    4. Component closest to centerline[0] = start, closest to centerline[-1] = end

    Args:
        boundary_edges: Mx2 array of boundary edge indices (into mesh.vertices)
        boundary_vertices: Nx3 array of boundary vertex positions
        mesh: Original mesh (for vertex indexing)
        centerline: Kx3 array of centerline points

    Returns:
        Tuple of (start_boundary, end_boundary)
    """
    if len(boundary_vertices) == 0 or len(boundary_edges) == 0:
        return np.array([]), np.array([])

    # Build vertex index mapping: mesh vertex index -> boundary vertex index
    boundary_vertex_indices = np.unique(boundary_edges.flatten())
    vertex_map = {mesh_idx: i for i, mesh_idx in enumerate(boundary_vertex_indices)}

    # Build adjacency graph from boundary edges
    # Remap edge indices from mesh space to boundary space
    n_boundary = len(boundary_vertices)
    from scipy.sparse import lil_matrix
    from scipy.sparse.csgraph import connected_components

    adjacency = lil_matrix((n_boundary, n_boundary), dtype=int)
    for edge in boundary_edges:
        i = vertex_map[edge[0]]
        j = vertex_map[edge[1]]
        adjacency[i, j] = 1
        adjacency[j, i] = 1

    # Find connected components
    n_components, labels = connected_components(adjacency, directed=False)

    print(f"  Found {n_components} connected boundary loops")

    if n_components == 0:
        return np.array([]), np.array([])
    elif n_components == 1:
        # Only one loop - use Z-coordinate as fallback
        z_median = np.median(boundary_vertices[:, 2])
        closer_to_min = boundary_vertices[:, 2] < z_median
        return boundary_vertices[closer_to_min], boundary_vertices[~closer_to_min]

    # Compute centroid for each component
    component_centroids = []
    component_vertices = []

    for comp_id in range(n_components):
        mask = labels == comp_id
        comp_verts = boundary_vertices[mask]
        centroid = np.mean(comp_verts, axis=0)
        component_centroids.append(centroid)
        component_vertices.append(comp_verts)

    # Project each centroid onto centerline and compute arc length
    component_arc_lengths = []
    for centroid in component_centroids:
        # Find closest centerline point
        distances = np.linalg.norm(centerline - centroid, axis=1)
        closest_idx = np.argmin(distances)

        # Arc length = index along centerline (0 = start, N-1 = end)
        arc_length = closest_idx
        component_arc_lengths.append(arc_length)

    # Smallest arc length = start boundary
    # Largest arc length = end boundary
    start_comp_id = np.argmin(component_arc_lengths)
    end_comp_id = np.argmax(component_arc_lengths)

    print(f"  Start loop: component {start_comp_id} at centerline index {component_arc_lengths[start_comp_id]}")
    print(f"  End loop: component {end_comp_id} at centerline index {component_arc_lengths[end_comp_id]}")

    start_boundary = component_vertices[start_comp_id]
    end_boundary = component_vertices[end_comp_id]

    return start_boundary, end_boundary


def fit_plane_to_points(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit a plane to a set of 3D points using PCA

    The plane is defined by a position (centroid) and normal vector.

    Args:
        points: Nx3 array of points

    Returns:
        Tuple of (plane_position, plane_normal)
    """
    if len(points) < 3:
        raise ValueError("Need at least 3 points to fit a plane")

    # Compute centroid
    centroid = np.mean(points, axis=0)

    # Center points
    centered = points - centroid

    # Compute covariance matrix
    cov = np.cov(centered.T)

    # Get eigenvectors (principal components)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)

    # Normal is eigenvector with smallest eigenvalue
    # (perpendicular to plane of best fit)
    normal = eigenvectors[:, 0]

    # Normalize
    normal = normal / np.linalg.norm(normal)

    return centroid, normal


def order_boundary_loop(boundary_points: np.ndarray,
                        plane_normal: np.ndarray) -> np.ndarray:
    """
    Order boundary points to form a closed loop

    Uses nearest neighbor traversal to create ordered boundary.

    Args:
        boundary_points: Nx3 array of boundary points
        plane_normal: Normal vector of fitted plane

    Returns:
        Ordered Nx3 array forming closed loop
    """
    if len(boundary_points) < 3:
        return boundary_points

    # Start with first point
    ordered = [boundary_points[0]]
    remaining = list(range(1, len(boundary_points)))
    current_idx = 0

    # Greedily connect to nearest neighbor
    while remaining:
        current_point = boundary_points[current_idx]

        # Find nearest remaining point
        distances = np.linalg.norm(boundary_points[remaining] - current_point, axis=1)
        nearest_idx = remaining[np.argmin(distances)]

        ordered.append(boundary_points[nearest_idx])
        remaining.remove(nearest_idx)
        current_idx = nearest_idx

    return np.array(ordered)


def create_boundary_cap_cross_section(boundary_vertices: np.ndarray,
                                       plane_position: np.ndarray,
                                       plane_normal: np.ndarray,
                                       reference_point: np.ndarray) -> Optional[CrossSection]:
    """
    Create a CrossSection object from boundary vertices

    This "caps" the mesh boundary by creating a triangulated surface.

    Args:
        boundary_vertices: Nx3 array of boundary vertex positions
        plane_position: Fitted plane position (centroid)
        plane_normal: Fitted plane normal
        reference_point: Reference point (first/last centerline point) to orient normal

    Returns:
        CrossSection object or None if creation fails
    """
    if len(boundary_vertices) < 3:
        print(f"Warning: Too few boundary vertices ({len(boundary_vertices)})")
        return None

    try:
        # Orient normal toward reference point (interior)
        direction_to_ref = reference_point - plane_position
        if np.dot(plane_normal, direction_to_ref) < 0:
            plane_normal = -plane_normal

        # Project boundary vertices onto fitted plane
        # This creates a planar boundary loop
        from .geometry import rotate_to_xy_plane

        # Rotate to XY plane for 2D triangulation
        rotated_vertices, rotation = rotate_to_xy_plane(boundary_vertices, plane_normal)

        # Extract 2D coordinates (X, Y)
        boundary_2d = rotated_vertices[:, :2]

        # Order points to form closed loop
        # Reorder based on angle from centroid
        centroid_2d = np.mean(boundary_2d, axis=0)
        angles = np.arctan2(boundary_2d[:, 1] - centroid_2d[1],
                           boundary_2d[:, 0] - centroid_2d[0])
        sorted_indices = np.argsort(angles)
        boundary_2d_ordered = boundary_2d[sorted_indices]

        # Create Shapely polygon to compute area
        try:
            poly = Polygon(boundary_2d_ordered)
            if not poly.is_valid:
                poly = poly.buffer(0)  # Fix invalid geometry
            area = poly.area
            perimeter = poly.length
        except Exception as e:
            print(f"Warning: Shapely polygon creation failed: {e}")
            # Fallback: use shoelace formula
            area = 0.5 * abs(sum(boundary_2d_ordered[i,0] * boundary_2d_ordered[i+1,1] -
                                 boundary_2d_ordered[i+1,0] * boundary_2d_ordered[i,1]
                                 for i in range(len(boundary_2d_ordered)-1)))
            perimeter = sum(np.linalg.norm(boundary_2d_ordered[i+1] - boundary_2d_ordered[i])
                           for i in range(len(boundary_2d_ordered)-1))
            perimeter += np.linalg.norm(boundary_2d_ordered[0] - boundary_2d_ordered[-1])

        # Triangulate the boundary loop using Delaunay
        try:
            tri = Delaunay(boundary_2d_ordered)
            faces_2d = tri.simplices
        except Exception as e:
            print(f"Warning: Delaunay triangulation failed: {e}")
            return None

        # Use reordered original boundary vertices (preserves 3D positions)
        vertices_3d_final = boundary_vertices[sorted_indices]

        # Create CrossSection object
        cross_section = CrossSection(
            vertices=vertices_3d_final,
            faces=faces_2d,
            area=area,
            centroid=plane_position,
            perimeter=perimeter,
            boundary_2d=boundary_2d_ordered,
            side=3  # Mark as "joined" / special cap
        )

        return cross_section

    except Exception as e:
        print(f"Warning: Boundary cap creation failed: {e}")
        return None


def extract_boundary_caps(mesh: trimesh.Trimesh,
                          centerline: np.ndarray) -> Tuple[Optional[CrossSection],
                                                           Optional[CrossSection]]:
    """
    Extract start and end boundary caps from mesh

    Main function that combines all steps:
    1. Extract boundary edges
    2. Separate into start/end regions
    3. Fit planes to each region
    4. Create triangulated cross-sections

    Args:
        mesh: Input mesh (can be non-watertight)
        centerline: Nx3 array of centerline points

    Returns:
        Tuple of (start_cap, end_cap)
        - start_cap: CrossSection for mesh start boundary (plane 0)
        - end_cap: CrossSection for mesh end boundary (plane N)
        Both can be None if boundary extraction fails
    """
    print("\\nExtracting boundary caps...")

    # Extract boundary edges
    boundary_edges, boundary_vertices = extract_boundary_edges(mesh)

    if len(boundary_vertices) == 0:
        print("  No boundary edges found (mesh is watertight)")
        return None, None

    print(f"  Found {len(boundary_edges)} boundary edges")
    print(f"  Found {len(boundary_vertices)} boundary vertices")

    # Separate into start and end regions
    start_boundary, end_boundary = separate_boundary_regions(boundary_edges, boundary_vertices, mesh, centerline)

    print(f"  Start boundary: {len(start_boundary)} vertices")
    print(f"  End boundary: {len(end_boundary)} vertices")

    # Create start cap
    start_cap = None
    if len(start_boundary) >= 3:
        try:
            start_pos, start_normal = fit_plane_to_points(start_boundary)
            start_cap = create_boundary_cap_cross_section(
                start_boundary,
                start_pos,
                start_normal,
                centerline[0]  # Reference: first centerline point
            )
            if start_cap:
                print(f"  ✓ Start cap: area={start_cap.area:.2f} mm²")
        except Exception as e:
            print(f"  ✗ Start cap creation failed: {e}")

    # Create end cap
    end_cap = None
    if len(end_boundary) >= 3:
        try:
            end_pos, end_normal = fit_plane_to_points(end_boundary)
            end_cap = create_boundary_cap_cross_section(
                end_boundary,
                end_pos,
                end_normal,
                centerline[-1]  # Reference: last centerline point
            )
            if end_cap:
                print(f"  ✓ End cap: area={end_cap.area:.2f} mm²")
        except Exception as e:
            print(f"  ✗ End cap creation failed: {e}")

    return start_cap, end_cap
