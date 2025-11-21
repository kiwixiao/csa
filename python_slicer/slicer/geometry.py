"""
Geometry utilities for 3D operations
Replaces: AxelRot.m and normal calculation logic from trachSlice_OSA_AB_3.m
"""

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.signal import savgol_filter
from typing import Tuple, Optional


def compute_plane_normal(centerline_points: np.ndarray, index: int,
                        smooth: bool = True, znormal: bool = False) -> np.ndarray:
    """
    Compute normal vector for plane at given centerline point

    Replaces: normals calculation in trachSlice_OSA_AB_3.m (lines 81-96)

    The normal is computed as the tangent to the centerline (forward difference).
    For middle points, use central difference. For end points, use forward/backward difference.

    Args:
        centerline_points: Nx3 array of centerline points
        index: Index of current point
        smooth: If True, apply smoothing to normals (replaces smoothdata in MATLAB)
        znormal: If True, use fixed [0, 0, 1] normal (for special cases)

    Returns:
        3D unit normal vector
    """
    n_points = len(centerline_points)

    if znormal:
        # Fixed Z-axis normal (used for balloon experiments in MATLAB)
        return np.array([0.0, 0.0, 1.0])

    # Compute tangent vector based on position
    if index == 0:
        # First point: forward difference
        normal = centerline_points[1] - centerline_points[0]
    elif index == n_points - 1:
        # Last point: backward difference
        normal = centerline_points[-1] - centerline_points[-2]
    else:
        # Middle points: central difference (MATLAB: pos(i+1,:) - pos(i-1,:))
        normal = centerline_points[index + 1] - centerline_points[index - 1]

    # Normalize
    normal = normal / np.linalg.norm(normal)

    return normal


def compute_all_plane_normals(centerline_points: np.ndarray,
                              smooth: bool = True,
                              znormal: bool = False) -> np.ndarray:
    """
    Compute normal vectors for all centerline points

    Replaces: normals calculation loop in trachSlice_OSA_AB_3.m

    Args:
        centerline_points: Nx3 array of centerline points
        smooth: If True, apply smoothing to normals
        znormal: If True, use fixed [0, 0, 1] normal

    Returns:
        Nx3 array of unit normal vectors
    """
    n_points = len(centerline_points)
    normals = np.zeros_like(centerline_points)

    # Compute centerline normals for all points
    for i in range(n_points):
        normals[i] = compute_plane_normal(centerline_points, i, smooth=False, znormal=znormal)

    # Apply smoothing FIRST (MATLAB: normals = smoothdata(normals))
    if smooth and not znormal:
        # Use moving average to match MATLAB's smoothdata() default behavior
        # MATLAB smoothdata() uses automatic window size, we use 20 for smoother normals
        window = 20
        if n_points >= window:
            # Pad edges to handle boundaries
            normals_padded = np.pad(normals, ((window//2, window//2), (0, 0)), mode='edge')
            # Apply moving average
            from scipy.ndimage import uniform_filter1d
            normals_smooth = uniform_filter1d(normals_padded, size=window, axis=0)
            # Remove padding
            normals = normals_smooth[window//2:window//2+n_points]
            # Re-normalize after smoothing
            norms = np.linalg.norm(normals, axis=1, keepdims=True)
            norms[norms == 0] = 1  # Avoid division by zero
            normals = normals / norms

    # CRITICAL FIX: Ensure normals are consistently oriented along centerline
    # If consecutive normals point in opposite directions, flip them
    if not znormal and n_points > 1:
        for i in range(1, n_points):
            # Check if normal[i] points in similar direction to normal[i-1]
            dot_product = np.dot(normals[i], normals[i-1])
            if dot_product < 0:
                # Normals are pointing in opposite directions - flip current one
                normals[i] = -normals[i]

    return normals


def rotate_to_xy_plane(points_3d: np.ndarray,
                       plane_normal: np.ndarray) -> Tuple[np.ndarray, Rotation]:
    """
    Rotate 3D points so that the plane defined by plane_normal becomes the XY plane

    Replaces: AxelRot rotation logic in NewMeshPlane (sublobeMulti_DG_AB4.m, lines 974-986)

    Args:
        points_3d: Nx3 array of 3D points
        plane_normal: 3D unit normal vector of the plane

    Returns:
        Tuple of:
            - Nx3 array of rotated points (Z coordinate should be constant)
            - Rotation object for inverse transformation
    """
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    z_axis = np.array([0.0, 0.0, 1.0])

    # Check if already aligned with Z-axis
    if np.allclose(plane_normal, z_axis) or np.allclose(plane_normal, -z_axis):
        # No rotation needed
        return points_3d.copy(), Rotation.from_matrix(np.eye(3))

    # Compute rotation axis (cross product) and angle
    # MATLAB: trans = cross(PlaneNormal, zaxis)
    rotation_axis = np.cross(plane_normal, z_axis)
    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)

    # Compute angle
    # MATLAB: theta = acosd(dot(PlaneNormal, zaxis)/(norm(PlaneNormal)*norm(zaxis)))
    cos_theta = np.dot(plane_normal, z_axis)
    theta_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))

    # Create rotation using axis-angle representation
    # MATLAB uses AxelRot which implements Rodrigues' formula
    rotation = Rotation.from_rotvec(theta_rad * rotation_axis)

    # Apply rotation
    rotated_points = rotation.apply(points_3d)

    return rotated_points, rotation


def rotate_back_to_3d(points_2d: np.ndarray,
                     rotation: Rotation,
                     z_coord: float = 0.0) -> np.ndarray:
    """
    Rotate 2D points (in XY plane) back to original 3D orientation

    Replaces: Inverse AxelRot in NewMeshPlane (sublobeMulti_DG_AB4.m, lines 1063-1066)

    Args:
        points_2d: Nx2 array of 2D points
        rotation: Rotation object from rotate_to_xy_plane
        z_coord: Z coordinate to use for all points (default 0)

    Returns:
        Nx3 array of 3D points in original orientation
    """
    # Add Z coordinate
    points_3d = np.column_stack([points_2d, np.full(len(points_2d), z_coord)])

    # Apply inverse rotation
    rotated_back = rotation.inv().apply(points_3d)

    return rotated_back


def compute_arc_length(points: np.ndarray) -> np.ndarray:
    """
    Compute cumulative arc length along a path of points

    Replaces: arcLength calculation in trachSlice_OSA_AB_3.m (lines 233-237)

    Args:
        points: Nx3 array of points

    Returns:
        N-length array of cumulative distances from first point
    """
    n_points = len(points)
    arc_length = np.zeros(n_points)

    # MATLAB: arcLength(i) = arcLength(i-1) + norm(goodpos(i,:) - goodpos(i-1,:))
    for i in range(1, n_points):
        arc_length[i] = arc_length[i - 1] + np.linalg.norm(points[i] - points[i - 1])

    return arc_length


def validate_planes_no_intersection(plane_positions: np.ndarray,
                                    plane_normals: np.ndarray,
                                    min_distance: float = 0.5) -> np.ndarray:
    """
    Check if adjacent planes intersect (quality control)

    New feature: Ensure planes are well-separated and don't intersect

    Args:
        plane_positions: Nx3 array of plane origins
        plane_normals: Nx3 array of plane normal vectors
        min_distance: Minimum distance between adjacent plane positions (mm)

    Returns:
        Boolean array indicating which planes are valid (True = valid, no intersection)
    """
    n_planes = len(plane_positions)
    is_valid = np.ones(n_planes, dtype=bool)

    for i in range(n_planes - 1):
        # Check distance between adjacent planes
        distance = np.linalg.norm(plane_positions[i + 1] - plane_positions[i])

        if distance < min_distance:
            # Planes too close - mark second one as invalid
            is_valid[i + 1] = False

    return is_valid


def check_adjacent_planes_intersect(plane1: Tuple[np.ndarray, np.ndarray],
                                   plane2: Tuple[np.ndarray, np.ndarray],
                                   tolerance: float = 0.1) -> bool:
    """
    Check if two adjacent planes intersect each other (pre-check using origins only)

    This is a preliminary check before slicing. Full validation happens
    after slicing using check_plane_vertices_same_side().

    Args:
        plane1: Tuple of (position, normal) for first plane
        plane2: Tuple of (position, normal) for second plane
        tolerance: Minimum distance plane2 origin must be ahead of plane1 (mm)

    Returns:
        True if planes likely intersect, False otherwise
    """
    pos1, normal1 = plane1
    pos2, normal2 = plane2

    # Pre-check: Verify plane2 origin is ahead of plane1
    signed_distance = np.dot(pos2 - pos1, normal1)
    return signed_distance < tolerance


def check_plane_vertices_same_side(vertices: np.ndarray,
                                   previous_plane_pos: np.ndarray,
                                   previous_plane_normal: np.ndarray,
                                   tolerance: float = 0.01) -> bool:
    """
    Check if ALL vertices of a plane are on the same side of the previous plane

    This is the PROPER sequential validation: ensures plane N+1 doesn't
    fold back and intersect plane N by verifying all vertices are ahead.

    Args:
        vertices: Nx3 array of vertex positions for current plane
        previous_plane_pos: Position of previous plane
        previous_plane_normal: Normal of previous plane
        tolerance: Minimum distance vertices must be ahead (mm)

    Returns:
        True if all vertices on same side (valid), False if any vertex on wrong side
    """
    if len(vertices) == 0:
        return True

    # Compute signed distance for each vertex to previous plane
    # signed_dist = dot((vertex - prev_pos), prev_normal)
    # Positive = ahead, Negative = behind
    vectors = vertices - previous_plane_pos
    signed_distances = np.dot(vectors, previous_plane_normal)

    # All vertices should be on the positive side (ahead of previous plane)
    # with at least tolerance distance
    all_ahead = np.all(signed_distances >= tolerance)

    return all_ahead


def convert_units_if_needed(arc_length: np.ndarray, area: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert from meters to millimeters if needed

    Replaces: Unit conversion in trachSlice_OSA_AB_3.m (lines 239-243)

    Args:
        arc_length: Array of arc lengths
        area: Array of areas

    Returns:
        Tuple of (arc_length, area) in millimeters and mm^2
    """
    # MATLAB: if max(arcLength) < 1 then data is in meters
    if np.max(arc_length) < 1.0:
        arc_length = arc_length * 1000.0  # m to mm
        area = area * 1e6  # m^2 to mm^2

    return arc_length, area


def get_plane_equation(position: np.ndarray, normal: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Get plane equation in form: ax + by + cz = d

    Args:
        position: Point on the plane
        normal: Normal vector of the plane

    Returns:
        Tuple of (normal_vector, d_value) where plane equation is normal·x = d
    """
    normal = normal / np.linalg.norm(normal)
    d = np.dot(normal, position)
    return normal, d


def distance_point_to_plane(point: np.ndarray, plane_position: np.ndarray,
                           plane_normal: np.ndarray) -> float:
    """
    Compute perpendicular distance from point to plane

    Args:
        point: 3D point
        plane_position: Point on the plane
        plane_normal: Normal vector of the plane (must be unit vector)

    Returns:
        Signed distance (positive if point is on normal side)
    """
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    return np.dot(point - plane_position, plane_normal)
