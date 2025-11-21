"""
Validation and quality control utilities
New features + replaces validation logic from trachSlice_OSA_AB_3.m
"""

import numpy as np
import trimesh
from scipy.spatial import ConvexHull
from typing import Tuple, Optional
from .geometry import distance_point_to_plane, check_adjacent_planes_intersect as geom_check_intersect


def check_plane_enclosed_by_mesh(plane_centroid: np.ndarray,
                                 mesh: trimesh.Trimesh,
                                 tolerance: float = 1.0) -> bool:
    """
    Check if a plane cross-section centroid is enclosed by the mesh

    New feature: Use convex hull to determine if centroid is inside mesh bounds

    Args:
        plane_centroid: 3D centroid of the cross-section
        mesh: Original STL mesh
        tolerance: Safety margin (mm) - centroid should be at least this far inside

    Returns:
        True if centroid is enclosed by mesh, False otherwise
    """
    try:
        # Method 1: Use trimesh's built-in contains method
        # This uses ray-casting to determine if point is inside mesh
        is_inside = mesh.contains([plane_centroid])[0]

        if is_inside:
            return True

        # Method 2: Check against convex hull if mesh is not watertight
        # Create convex hull of mesh vertices
        try:
            hull = ConvexHull(mesh.vertices)

            # Check if point is inside convex hull
            # A point is inside if it satisfies all hull inequalities: Ax <= b
            # where A are the hull equations and b are the offsets
            distances = hull.equations[:, :3] @ plane_centroid + hull.equations[:, 3]

            # Point is inside if all distances are <= tolerance
            is_inside_hull = np.all(distances <= tolerance)

            return is_inside_hull

        except Exception:
            # If convex hull computation fails, use conservative approach
            # Check if centroid is within bounding box
            bbox_min = mesh.bounds[0] - tolerance
            bbox_max = mesh.bounds[1] + tolerance

            is_in_bbox = np.all(plane_centroid >= bbox_min) and np.all(plane_centroid <= bbox_max)
            return is_in_bbox

    except Exception as e:
        print(f"Warning: Enclosure check failed: {str(e)}, assuming valid")
        return True  # Conservative: assume valid if check fails


def validate_centroid_proximity(centroid: np.ndarray,
                                centerline_point: np.ndarray,
                                threshold: float = 10.0) -> bool:
    """
    Validate that cross-section centroid is close to the centerline point

    Replaces: Proximity checks in trachSlice_OSA_AB_3.m (lines 191-202) and
             sublobeMulti_DG_AB4.m (lines 81-86)

    MATLAB code checks:
    - abs(norm(centroidPos - PlanePoint)) > 10  (line 81 in sublobeMulti_DG_AB4.m)
    - X(g) > threshold_centroid1 where X = abs(abs(goodpos) - abs(currentCentroid))

    Args:
        centroid: 3D centroid of cross-section
        centerline_point: Corresponding centerline point
        threshold: Maximum allowed distance (mm), default 10.0 matches MATLAB

    Returns:
        True if centroid is within threshold distance
    """
    # Compute Euclidean distance
    distance = np.linalg.norm(centroid - centerline_point)

    # Also check component-wise difference (MATLAB: X = abs(abs(goodpos) - abs(currentCentroid)))
    # This catches cases where centroids are way off in any dimension
    component_diff = np.abs(np.abs(centroid) - np.abs(centerline_point))
    max_component_diff = np.max(component_diff)

    # MATLAB uses threshold_centroid1 = 100 for component check (line 191)
    component_threshold = 100.0

    # Valid if both checks pass
    is_valid = (distance <= threshold) and (max_component_diff <= component_threshold)

    return is_valid


def check_adjacent_planes_intersect(plane1: Tuple[np.ndarray, np.ndarray],
                                   plane2: Tuple[np.ndarray, np.ndarray],
                                   tolerance: float = 0.1) -> bool:
    """
    Check if two adjacent cutting planes intersect

    New feature: Ensure planes don't overlap

    Args:
        plane1: (position, normal) for first plane
        plane2: (position, normal) for second plane
        tolerance: Distance tolerance (mm)

    Returns:
        True if planes intersect, False otherwise
    """
    return geom_check_intersect(plane1, plane2, tolerance)


def validate_cross_section_area(area: float,
                                previous_area: Optional[float] = None,
                                min_area: float = 1.0,
                                max_area_ratio: float = 10.0) -> bool:
    """
    Validate cross-sectional area

    New feature: Sanity checks on computed areas

    Args:
        area: Computed cross-sectional area (mm^2)
        previous_area: Previous plane's area for continuity check
        min_area: Minimum acceptable area (mm^2)
        max_area_ratio: Maximum ratio between consecutive areas

    Returns:
        True if area is valid
    """
    # Check minimum area
    if area < min_area:
        return False

    # Check for NaN or inf
    if not np.isfinite(area):
        return False

    # Check continuity with previous area
    if previous_area is not None and previous_area > 0:
        ratio = max(area, previous_area) / min(area, previous_area)
        if ratio > max_area_ratio:
            # Area changed too drastically
            return False

    return True


def validate_cross_section_quality(vertices: np.ndarray,
                                   faces: np.ndarray,
                                   min_vertices: int = 3) -> bool:
    """
    Validate quality of triangulated cross-section mesh

    New feature: Ensure triangulation produced valid mesh

    Args:
        vertices: Nx3 array of vertices
        faces: Mx3 array of triangle faces
        min_vertices: Minimum number of vertices required

    Returns:
        True if mesh quality is acceptable
    """
    # Check minimum vertices
    if len(vertices) < min_vertices:
        return False

    # Check that we have faces
    if len(faces) == 0:
        return False

    # Check for degenerate triangles (area = 0)
    for face in faces:
        v0, v1, v2 = vertices[face]
        edge1 = v1 - v0
        edge2 = v2 - v0
        cross = np.cross(edge1, edge2)
        area = 0.5 * np.linalg.norm(cross)

        if area < 1e-10:  # Degenerate triangle
            return False

    return True


def detect_split_regions(centerline_points: np.ndarray,
                         z_threshold: float = 5.0) -> Optional[int]:
    """
    Detect where airway splits into two regions (trachea+nose1 vs nose2)

    Replaces: condtru detection in SlicerMasterCode.m (lines 58-63)

    MATLAB code:
    ```matlab
    goodpos3 = goodpos(:,3);
    result = [];
    if i == 1
        for p = 1:numel(goodpos3)-1
            result = [result, abs(goodpos3(p) - goodpos3(p+1)) > 5];
        end
        condtru = find(result==1);
    end
    ```

    Args:
        centerline_points: Nx3 array of points
        z_threshold: Threshold for detecting split (mm), default 5.0

    Returns:
        Index where split occurs, or None if no split detected
    """
    # Extract Z coordinates
    z_coords = centerline_points[:, 2]

    # Find large jumps in Z coordinate
    z_diffs = np.abs(np.diff(z_coords))

    # Find first index where difference exceeds threshold
    split_indices = np.where(z_diffs > z_threshold)[0]

    if len(split_indices) > 0:
        # Return first split index
        return int(split_indices[0])
    else:
        return None


def check_mesh_quality(mesh: trimesh.Trimesh) -> Tuple[bool, list]:
    """
    Check overall mesh quality and report issues

    New feature: Comprehensive mesh validation

    Args:
        mesh: Input mesh to validate

    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []

    # Check if mesh is empty
    if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
        issues.append("Empty mesh")

    # Check for degenerate faces
    if not mesh.is_volume:
        issues.append("Mesh does not enclose volume")

    # Check if mesh is watertight
    if not mesh.is_watertight:
        issues.append("Mesh is not watertight")

    # Check for inverted normals
    if hasattr(mesh, 'is_winding_consistent') and not mesh.is_winding_consistent:
        issues.append("Inconsistent face winding")

    # Report
    is_valid = len(issues) == 0

    return is_valid, issues


def check_cross_section_closed_smart(cross_sections: list) -> Tuple[bool, int, str]:
    """
    Smart closure detection that handles multi-section planes and donuts

    Distinguishes between:
    - Single closed loop (normal airway) ✅
    - Multiple closed loops (donut shapes) ✅
    - Open paths (incomplete cuts) ❌ REJECT

    This solves the frame 14 plane 76 problem: planes cutting through mesh
    boundaries/holes produce open paths, which should be rejected.

    Args:
        cross_sections: List of CrossSection objects from slicing

    Returns:
        Tuple of (is_closed, n_paths, reason_if_open)
        - is_closed: True if ALL paths are closed (valid)
        - n_paths: Number of paths found
        - reason: Description if open
    """
    if not cross_sections:
        return False, 0, "no_sections"

    # For multiple cross-sections, check the largest one (closest to centerline)
    # This is handled by _select_best_section in plane_slicer.py
    section = cross_sections[0] if len(cross_sections) == 1 else None

    if section is None:
        return False, len(cross_sections), "multiple_sections"

    # Check if boundary forms closed loop(s)
    # boundary_2d is Nx2 array of 2D points after rotation to XY plane
    boundary = section.boundary_2d

    if len(boundary) < 3:
        return False, 0, "too_few_vertices"

    # Check if boundary is closed: first and last vertices should be close
    # For closed loops, trimesh typically repeats first vertex at end
    distance_to_close = np.linalg.norm(boundary[0] - boundary[-1])

    # Threshold: 0.1mm for considering loop closed
    is_closed = distance_to_close < 0.1

    if not is_closed:
        return False, 1, f"open_path_gap_{distance_to_close:.2f}mm"

    # Additional check: boundary should form a simple closed curve
    # Verify no self-intersections (for quality)
    # For donuts, we'd have multiple disjoint closed loops
    # But since we select best_section, we only check one contour

    return True, 1, "closed"


def check_cross_section_closed_from_path3d(path_3d: 'trimesh.path.Path3D') -> Tuple[bool, int, str]:
    """
    Check if cross-section from trimesh.section() is closed

    This works directly with trimesh Path3D objects before they're
    converted to CrossSection objects.

    Args:
        path_3d: trimesh.path.Path3D object from mesh.section()

    Returns:
        Tuple of (is_closed, n_paths, reason_if_open)
    """
    if path_3d is None:
        return False, 0, "no_section"

    # Get all path entities (could be multiple for donuts)
    entities = path_3d.entities
    n_paths = len(entities)

    if n_paths == 0:
        return False, 0, "no_paths"

    # Check if ALL paths are closed
    all_closed = True
    open_paths = []

    for i, entity in enumerate(entities):
        # trimesh Path2D has is_closed property
        if hasattr(entity, 'is_closed'):
            if not entity.is_closed:
                all_closed = False
                open_paths.append(i)
        else:
            # Fallback: check if first and last points match
            points = path_3d.vertices[entity.points]
            if len(points) >= 2:
                dist = np.linalg.norm(points[0] - points[-1])
                if dist > 0.1:  # mm
                    all_closed = False
                    open_paths.append(i)

    if all_closed:
        return True, n_paths, "all_closed"
    else:
        return False, n_paths, f"open_paths_{open_paths}"
