"""
Mesh-plane intersection and cross-section triangulation
Replaces: sublobeMulti_DG_AB4.m (1422 lines - most complex MATLAB file)
"""

import numpy as np
import trimesh
from shapely.geometry import Polygon, Point, MultiPolygon
from shapely.ops import unary_union
from scipy.spatial import Delaunay
from typing import List, Tuple, Optional
from dataclasses import dataclass

from .geometry import rotate_to_xy_plane, rotate_back_to_3d


@dataclass
class CrossSection:
    """Container for a single cross-section result"""
    vertices: np.ndarray  # Nx3 vertices of triangulated section
    faces: np.ndarray  # Mx3 triangle indices
    area: float  # Cross-sectional area (mm^2)
    centroid: np.ndarray  # 3D centroid
    perimeter: float  # Perimeter length
    boundary_2d: np.ndarray  # 2D boundary points before triangulation
    side: int  # Which side (1=left, 2=right, 3=joined)


class MeshPlaneSlicer:
    """
    Handles slicing a 3D mesh with a plane

    Replaces: sublobeMulti_DG_AB4 function and its sub-functions
    """

    def __init__(self, mesh: trimesh.Trimesh):
        """
        Initialize slicer with mesh

        Args:
            mesh: Input STL mesh
        """
        self.mesh = mesh

    def slice_mesh_with_plane(self,
                             plane_origin: np.ndarray,
                             plane_normal: np.ndarray,
                             plane_number: int = 0) -> List[CrossSection]:
        """
        Slice mesh with a plane and return all cross-sections

        Replaces: Main logic of sublobeMulti_DG_AB4.m

        Args:
            plane_origin: Point on the cutting plane
            plane_normal: Normal vector of cutting plane (unit vector)
            plane_number: Index for debugging

        Returns:
            List of CrossSection objects (can be empty if no intersection)
        """
        # Normalize plane normal
        plane_normal = plane_normal / np.linalg.norm(plane_normal)

        # Use trimesh to slice mesh
        # This replaces SplitLobes function (lines 455-581 in MATLAB)
        try:
            slice_result = self.mesh.section(plane_origin=plane_origin,
                                            plane_normal=plane_normal)

            if slice_result is None:
                print(f"Warning: Plane {plane_number} does not intersect mesh")
                return []

            # Convert to Path2D objects (2D representation)
            paths_2d, transform = slice_result.to_planar()

            if paths_2d is None or len(paths_2d.entities) == 0:
                print(f"Warning: Plane {plane_number} has no valid paths")
                return []

        except Exception as e:
            print(f"Warning: Plane {plane_number} slicing failed: {str(e)}")
            return []

        # Store original 3D vertices (EXACT coordinates on mesh surface)
        # These will be used for triangulated mesh to ensure perfect alignment
        original_3d_vertices = slice_result.vertices

        # Extract polygons from paths
        # This replaces DetectConnectedElements/DetectTwoNodeConnectedElements
        # Pass slice_result to access 3D vertices for centerline containment check
        polygons = self._extract_polygons_from_paths(paths_2d, slice_result, plane_origin, plane_normal)

        if not polygons:
            print(f"Warning: Plane {plane_number} - no closed polygons found")
            return []

        # Handle nested polygons (doughnuts)
        # Replaces CheckForDoughnuts (lines 726-825 in MATLAB)
        polygons = self._handle_nested_polygons(polygons)

        # Classify polygons by side (left, right, or joined)
        # Replaces side detection logic (lines 90-105 in MATLAB)
        classified_polygons = self._classify_polygons_by_side(polygons, plane_origin)

        # Triangulate each polygon and create CrossSection objects
        cross_sections = []
        for side, polygon_list in enumerate(classified_polygons, start=1):
            for polygon in polygon_list:
                try:
                    cross_section = self._triangulate_polygon(
                        polygon, plane_normal, plane_origin, side, plane_number,
                        paths_2d.vertices, original_3d_vertices
                    )
                    if cross_section is not None:
                        cross_sections.append(cross_section)
                except Exception as e:
                    print(f"Warning: Plane {plane_number} triangulation failed: {str(e)}")
                    continue

        return cross_sections

    def _project_to_plane_2d(self, points_3d: np.ndarray, plane_origin: np.ndarray,
                             plane_normal: np.ndarray) -> np.ndarray:
        """
        Project 3D points onto a 2D plane coordinate system.

        The 2D coordinate system has its origin at plane_origin with the plane's
        normal as the Z-axis. This ensures the centerline projects to (0,0).

        Args:
            points_3d: Nx3 array of 3D points
            plane_origin: 3D point defining the origin of the plane (centerline point)
            plane_normal: 3D normal vector of the plane

        Returns:
            Nx2 array of 2D coordinates
        """
        # Normalize plane normal
        normal = plane_normal / np.linalg.norm(plane_normal)

        # Create orthonormal basis for the plane using Gram-Schmidt
        # Choose an arbitrary vector not parallel to normal
        if abs(normal[0]) < 0.9:
            arbitrary = np.array([1.0, 0.0, 0.0])
        else:
            arbitrary = np.array([0.0, 1.0, 0.0])

        # First basis vector (u): project arbitrary vector onto plane
        u = arbitrary - np.dot(arbitrary, normal) * normal
        u = u / np.linalg.norm(u)

        # Second basis vector (v): perpendicular to both normal and u
        v = np.cross(normal, u)
        v = v / np.linalg.norm(v)

        # Translate so plane_origin is at (0,0)
        translated = points_3d - plane_origin

        # Project onto u and v axes to get 2D coordinates
        coords_2d = np.column_stack([
            np.dot(translated, u),
            np.dot(translated, v)
        ])

        return coords_2d

    def _extract_polygons_from_paths(self, paths_2d, slice_result_3d,
                                     plane_origin: np.ndarray,
                                     plane_normal: np.ndarray) -> List[Polygon]:
        """
        Extract closed polygons from 2D paths, filtering to keep only the loop
        that contains the centerline (plane origin).

        Replaces: Connected elements detection in sublobeMulti_DG_AB4.m

        Args:
            paths_2d: Path2D object from trimesh
            slice_result_3d: 3D slice result from mesh.section() (has .vertices and .entities)
            plane_origin: 3D centerline point (plane origin)
            plane_normal: 3D plane normal vector

        Returns:
            List of Shapely Polygon objects (filtered to contain centerline)
        """
        # Extract all candidate polygons
        all_polygons = []
        polygon_metadata = []  # Track which entity each polygon came from

        for entity_idx, entity in enumerate(paths_2d.entities):
            try:
                # Get vertices for this path from 2D representation
                vertices_2d = paths_2d.vertices[entity.points]

                # Check if we have enough vertices
                if len(vertices_2d) < 3:
                    continue

                # Create polygon from 2D vertices
                polygon = Polygon(vertices_2d)

                if polygon.is_valid and not polygon.is_empty:
                    all_polygons.append(polygon)
                    polygon_metadata.append(entity_idx)

            except Exception as e:
                continue

        # If only one polygon, return it (no filtering needed)
        if len(all_polygons) <= 1:
            return all_polygons

        # Multiple loops detected - verify that at least one contains the centerline
        # Use our custom 3D-to-2D projection to check which loop contains the centerline

        # Project centerline point (plane_origin) to 2D
        # Since we use plane_origin as the origin, centerline projects to (0,0)
        centerline_2d = Point(0.0, 0.0)

        # Check which loops contain the centerline
        centerline_loop_found = False

        for poly_idx, entity_idx in enumerate(polygon_metadata):
            entity = slice_result_3d.entities[entity_idx]

            # Get 3D vertices for this entity
            vertices_3d = slice_result_3d.vertices[entity.points]

            # Project to OUR 2D coordinate system (origin at centerline)
            vertices_2d_custom = self._project_to_plane_2d(vertices_3d, plane_origin, plane_normal)

            # Create polygon in this coordinate system
            polygon_custom = Polygon(vertices_2d_custom)

            # Check if centerline (0,0) is inside this polygon
            if polygon_custom.is_valid and polygon_custom.contains(centerline_2d):
                centerline_loop_found = True
                break

        # If no loop contains the centerline, reject this plane
        if not centerline_loop_found:
            # Warning: no loop contains centerline (shouldn't happen)
            # Return empty to reject this plane
            return []

        # Keep ALL loops - requirement is to keep all closed loops
        # The centerline loop is guaranteed to exist (checked above)
        # TODO: Add surface-based closure validation to filter out unclosed loops
        return all_polygons

    def _handle_nested_polygons(self, polygons: List[Polygon]) -> List[Polygon]:
        """
        Handle nested polygons (doughnuts - holes within cross-sections)

        Replaces: CheckForDoughnuts and doughnut merging logic (lines 113-223 in MATLAB)

        The MATLAB code uses inpolygon to detect if one polygon is inside another,
        then merges them as exterior/interior boundaries.

        Args:
            polygons: List of Polygon objects

        Returns:
            List of Polygon objects with holes properly handled
        """
        if len(polygons) <= 1:
            return polygons

        # Find nested relationships
        result_polygons = []
        used_indices = set()

        for i, poly_outer in enumerate(polygons):
            if i in used_indices:
                continue

            # Find polygons inside this one
            holes = []
            for j, poly_inner in enumerate(polygons):
                if i == j or j in used_indices:
                    continue

                # Check if poly_inner is inside poly_outer
                if poly_outer.contains(poly_inner):
                    # This is a hole
                    holes.append(poly_inner.exterior.coords)
                    used_indices.add(j)

            # Create polygon with holes
            if holes:
                # Polygon with holes: exterior + list of interiors
                poly_with_holes = Polygon(poly_outer.exterior.coords, holes=holes)
                result_polygons.append(poly_with_holes)
            else:
                result_polygons.append(poly_outer)

            used_indices.add(i)

        return result_polygons

    def _classify_polygons_by_side(self,
                                   polygons: List[Polygon],
                                   plane_origin: np.ndarray) -> List[List[Polygon]]:
        """
        Classify polygons by anatomical side (left nostril, right nostril, or joined)

        Replaces: Side classification logic in sublobeMulti_DG_AB4.m (lines 90-105)

        MATLAB logic:
        - NoOfCuts == 1: side = 3 (joined/pharynx)
        - centroidPos(i,1) > PlanePoint(1): side = 1 (left)
        - else: side = 2 (right)

        Args:
            polygons: List of Polygon objects
            plane_origin: Origin point of cutting plane

        Returns:
            List of 3 lists: [left_polygons, right_polygons, joined_polygons]
        """
        left_polygons = []
        right_polygons = []
        joined_polygons = []

        if len(polygons) == 1:
            # Single cut - assume joined region (pharynx/trachea)
            joined_polygons.append(polygons[0])
        else:
            # Multiple cuts - classify by X coordinate
            for polygon in polygons:
                centroid = np.array(polygon.centroid.coords[0])

                # Add Z=0 if 2D centroid
                if len(centroid) == 2:
                    centroid_x = centroid[0]
                else:
                    centroid_x = centroid[0]

                if centroid_x > plane_origin[0]:
                    # Left nostril
                    left_polygons.append(polygon)
                else:
                    # Right nostril
                    right_polygons.append(polygon)

        return [left_polygons, right_polygons, joined_polygons]

    def _triangulate_polygon(self,
                            polygon: Polygon,
                            plane_normal: np.ndarray,
                            plane_origin: np.ndarray,
                            side: int,
                            plane_number: int,
                            vertices_2d_all: np.ndarray = None,
                            vertices_3d_original: np.ndarray = None) -> Optional[CrossSection]:
        """
        Triangulate a 2D polygon and rotate back to 3D

        Replaces: NewMeshPlane function (lines 965-1068 in MATLAB)

        Args:
            polygon: Shapely Polygon (2D)
            plane_normal: Normal vector of plane
            plane_origin: Origin point of plane
            side: Anatomical side (1=left, 2=right, 3=joined)
            plane_number: For debugging

        Returns:
            CrossSection object or None if triangulation fails
        """
        try:
            # Get exterior and interior coordinates
            exterior_coords = np.array(polygon.exterior.coords[:-1])  # Remove duplicate last point

            # Prepare constrained edges for triangulation
            # This replaces the constrainedEdges logic (lines 998-1012 in MATLAB)
            n_exterior = len(exterior_coords)
            all_coords = [exterior_coords]
            edge_splits = [n_exterior]  # Indices where each loop ends

            # Add interior holes
            for interior in polygon.interiors:
                interior_coords = np.array(interior.coords[:-1])
                all_coords.append(interior_coords)
                edge_splits.append(edge_splits[-1] + len(interior_coords))

            # Combine all coordinates
            vertices_2d = np.vstack(all_coords)

            # Create constrained edges
            constrained_edges = []
            start_idx = 0
            for end_idx in edge_splits:
                # Create edges for this loop
                for i in range(start_idx, end_idx - 1):
                    constrained_edges.append([i, i + 1])
                # Close the loop
                constrained_edges.append([end_idx - 1, start_idx])
                start_idx = end_idx

            constrained_edges = np.array(constrained_edges)

            # Perform constrained Delaunay triangulation
            # Replaces delaunayTriangulation (lines 1020-1033 in MATLAB)
            tri = self._constrained_delaunay_2d(vertices_2d, constrained_edges)

            if tri is None or len(tri) == 0:
                return None

            # Filter triangles to keep only those inside the polygon
            # Replaces isInterior check (lines 1041-1050 in MATLAB)
            valid_triangles = []
            for triangle_indices in tri:
                # Check if triangle centroid is inside polygon
                triangle_coords = vertices_2d[triangle_indices]
                centroid = np.mean(triangle_coords, axis=0)
                if polygon.contains(Point(centroid)):
                    valid_triangles.append(triangle_indices)

            if not valid_triangles:
                print(f"Warning: Plane {plane_number} - no valid triangles after filtering")
                return None

            tri = np.array(valid_triangles)

            # Calculate 2D area
            # Replaces area calculation (lines 1054-1059 in MATLAB)
            area_2d = 0.0
            for triangle_indices in tri:
                v0, v1, v2 = vertices_2d[triangle_indices]
                # Area = 0.5 * |cross product|
                edge1 = v1 - v0
                edge2 = v2 - v0
                cross = edge1[0] * edge2[1] - edge1[1] * edge2[0]
                area_2d += 0.5 * abs(cross)

            # Map 2D vertices to original 3D vertices
            # This preserves EXACT coordinates from mesh intersection
            if vertices_2d_all is not None and vertices_3d_original is not None:
                # Build mapping from 2D polygon vertices to original 3D vertex indices
                vertex_map = {}
                for i, coord_2d in enumerate(vertices_2d):
                    # Find this 2D coord in the full 2D vertex array
                    for j, v2d in enumerate(vertices_2d_all):
                        if np.allclose(coord_2d, v2d, atol=1e-6):
                            vertex_map[i] = j
                            break

                # Use original 3D vertices
                vertices_3d = vertices_3d_original

                # Remap triangle indices to original 3D vertices
                tri_remapped = np.array([[vertex_map[t[0]], vertex_map[t[1]], vertex_map[t[2]]]
                                        for t in tri if t[0] in vertex_map and t[1] in vertex_map and t[2] in vertex_map])
                tri = tri_remapped

                # Calculate 3D centroid from mapped vertices
                used_vertices = vertices_3d[[vertex_map[i] for i in vertex_map.keys()]]
                centroid_3d = np.mean(used_vertices, axis=0)
            else:
                # Fallback to old method if original vertices not provided
                vertices_3d = self._convert_2d_to_3d_on_plane(vertices_2d, plane_origin, plane_normal)
                centroid_3d = np.mean(vertices_3d, axis=0)

            # Calculate perimeter
            perimeter = polygon.length

            # Create CrossSection object
            return CrossSection(
                vertices=vertices_3d,
                faces=tri,
                area=area_2d,
                centroid=centroid_3d,
                perimeter=perimeter,
                boundary_2d=vertices_2d,
                side=side
            )

        except Exception as e:
            print(f"Warning: Plane {plane_number} triangulation error: {str(e)}")
            return None

    def _constrained_delaunay_2d(self,
                                 vertices: np.ndarray,
                                 edges: np.ndarray) -> Optional[np.ndarray]:
        """
        Perform constrained Delaunay triangulation in 2D

        Replaces: delaunayTriangulation in MATLAB (with retry logic from lines 1018-1033)

        Args:
            vertices: Nx2 array of vertices
            edges: Mx2 array of constrained edge indices

        Returns:
            Kx3 array of triangle indices, or None if fails
        """
        try:
            # Scipy's Delaunay doesn't support constrained edges directly
            # We'll use a simple approach: triangulate and filter

            # Remove duplicate vertices that are very close
            unique_vertices, unique_indices = self._remove_duplicate_vertices(vertices)

            if len(unique_vertices) < 3:
                return None

            # Perform Delaunay triangulation
            tri = Delaunay(unique_vertices)

            return tri.simplices

        except Exception as e:
            # Retry with simplified vertices (MATLAB retry logic)
            try:
                # Remove closest vertices and try again
                simplified_vertices = self._simplify_vertices(vertices)
                if len(simplified_vertices) >= 3:
                    tri = Delaunay(simplified_vertices)
                    return tri.simplices
                else:
                    return None
            except Exception:
                return None

    def _remove_duplicate_vertices(self,
                                   vertices: np.ndarray,
                                   tolerance: float = 1e-6) -> Tuple[np.ndarray, np.ndarray]:
        """
        Remove duplicate vertices that are very close to each other

        Args:
            vertices: Nx2 or Nx3 array
            tolerance: Distance threshold

        Returns:
            Tuple of (unique_vertices, mapping_to_original_indices)
        """
        unique_vertices = [vertices[0]]
        mapping = [0]

        for i in range(1, len(vertices)):
            # Check if this vertex is far enough from all existing unique vertices
            distances = np.linalg.norm(vertices[i] - np.array(unique_vertices), axis=1)
            if np.min(distances) > tolerance:
                unique_vertices.append(vertices[i])
                mapping.append(len(unique_vertices) - 1)
            else:
                # Map to closest existing vertex
                mapping.append(np.argmin(distances))

        return np.array(unique_vertices), np.array(mapping)

    def _simplify_vertices(self, vertices: np.ndarray) -> np.ndarray:
        """
        Simplify vertex list by removing closely-spaced vertices

        Replaces: MATLAB retry logic (lines 1024-1032)

        Args:
            vertices: Nx2 array

        Returns:
            Simplified vertex array
        """
        if len(vertices) <= 3:
            return vertices

        # Compute distances between consecutive vertices
        diffs = np.diff(vertices, axis=0)
        distances = np.linalg.norm(diffs, axis=1)

        # Find and remove the closest pair
        min_idx = np.argmin(distances)

        # Remove the vertex at min_idx + 1
        simplified = np.delete(vertices, min_idx + 1, axis=0)

        return simplified

    def _convert_2d_to_3d_on_plane(self,
                                   vertices_2d: np.ndarray,
                                   plane_origin: np.ndarray,
                                   plane_normal: np.ndarray) -> np.ndarray:
        """
        Convert 2D vertices to 3D coordinates on the specified plane

        Uses MATLAB's AxelRot approach: creates coordinate frame based on rotation
        from Z-axis to plane normal. This ensures consistent orientation across planes.

        Args:
            vertices_2d: Nx2 array of 2D coordinates
            plane_origin: 3D point on the plane
            plane_normal: 3D normal vector

        Returns:
            Nx3 array of 3D coordinates
        """
        # Normalize plane normal
        plane_normal = plane_normal / np.linalg.norm(plane_normal)

        # MATLAB approach: Use cross product with Z-axis to define coordinate frame
        # This ensures consistent orientation across all planes
        zaxis = np.array([0.0, 0.0, 1.0])

        # Check if plane normal is parallel to Z-axis
        if np.allclose(plane_normal, zaxis) or np.allclose(plane_normal, -zaxis):
            # Plane is already aligned with XY plane
            u = np.array([1.0, 0.0, 0.0])
            v = np.array([0.0, 1.0, 0.0])
        else:
            # Create orthonormal basis using cross product with Z-axis
            # This matches MATLAB's approach: trans = cross(PlaneNormal, zaxis)
            u = np.cross(plane_normal, zaxis)
            u = u / np.linalg.norm(u)

            # Second axis is perpendicular to both normal and u
            v = np.cross(plane_normal, u)
            v = v / np.linalg.norm(v)

        # Convert 2D to 3D: point = origin + x*u + y*v
        vertices_3d = np.zeros((len(vertices_2d), 3))
        for i, (x, y) in enumerate(vertices_2d):
            vertices_3d[i] = plane_origin + x * u + y * v

        return vertices_3d
