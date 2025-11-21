"""
Main airway slicer class
Replaces: trachSlice_OSA_AB_3.m (282 lines)
"""

import numpy as np
from typing import List, Tuple, Optional
from dataclasses import dataclass, field

from .io_utils import read_stl, read_vtk_centerline
from .geometry import (
    compute_all_plane_normals,
    compute_arc_length,
    convert_units_if_needed,
    check_adjacent_planes_intersect,
    check_plane_vertices_same_side
)
from .validation import (
    check_plane_enclosed_by_mesh,
    validate_centroid_proximity
)
from .mesh_intersection import MeshPlaneSlicer, CrossSection
from .measurements import DiameterProfile


@dataclass
class SlicingResults:
    """
    Container for all slicing results

    Replaces: Loose variables in trachSlice_OSA_AB_3.m return values
    """
    # Valid cross-sections
    valid_sections: List[CrossSection] = field(default_factory=list)
    valid_positions: List[np.ndarray] = field(default_factory=list)
    valid_plane_indices: List[int] = field(default_factory=list)

    # Skipped planes (with reasons)
    skipped_planes: List[Tuple[int, str]] = field(default_factory=list)

    # Arc length for each valid section
    arc_lengths: np.ndarray = field(default=None)

    # Diameter profile
    diameter_profile: Optional[DiameterProfile] = None

    def add_valid_plane(self,
                       index: int,
                       cross_section: CrossSection,
                       position: np.ndarray,
                       normal: np.ndarray) -> None:
        """Add a valid slicing plane"""
        self.valid_sections.append(cross_section)
        self.valid_positions.append(position)
        self.valid_plane_indices.append(index)

    def add_skipped_plane(self, index: int, reason: str) -> None:
        """Record a skipped plane with reason"""
        self.skipped_planes.append((index, reason))

    def compute_arc_lengths(self) -> None:
        """Compute cumulative arc lengths along centerline"""
        if len(self.valid_positions) > 0:
            positions_array = np.array(self.valid_positions)
            self.arc_lengths = compute_arc_length(positions_array)
        else:
            self.arc_lengths = np.array([])

    def compute_diameter_profile(self) -> DiameterProfile:
        """
        Compute diameter profile from valid sections

        Returns:
            DiameterProfile object with all measurements
        """
        if self.diameter_profile is None:
            self.diameter_profile = DiameterProfile()

            for i, (section, arc_len, idx) in enumerate(zip(
                self.valid_sections,
                self.arc_lengths,
                self.valid_plane_indices
            )):
                self.diameter_profile.add_cross_section(
                    arc_length=arc_len,
                    area=section.area,
                    perimeter=section.perimeter,
                    centroid=section.centroid,
                    boundary_2d=section.boundary_2d,
                    plane_index=idx,
                    is_valid=True
                )

        return self.diameter_profile

    def get_summary(self) -> dict:
        """Get summary statistics"""
        return {
            'n_valid_planes': len(self.valid_sections),
            'n_skipped_planes': len(self.skipped_planes),
            'n_total_planes': len(self.valid_sections) + len(self.skipped_planes),
            'skip_reasons': dict((reason, sum(1 for _, r in self.skipped_planes if r == reason))
                               for _, reason in self.skipped_planes)
        }


class AirwaySlicer:
    """
    Main class for slicing airway mesh along centerline

    Replaces: trachSlice_OSA_AB_3 function

    Usage:
        slicer = AirwaySlicer("mesh.stl", "centerline.vtk")
        results = slicer.slice_along_centerline()
    """

    def __init__(self, stl_path: str, vtk_centerline_path: str, znormal: bool = False):
        """
        Initialize slicer with mesh and centerline

        Args:
            stl_path: Path to STL mesh file
            vtk_centerline_path: Path to VTK centerline file
            znormal: If True, use fixed Z-axis normal (for special cases)
        """
        # Load mesh and centerline (replaces STL_Import and read_vtk)
        print(f"\nInitializing AirwaySlicer...")
        self.mesh = read_stl(stl_path)
        self.centerline = read_vtk_centerline(vtk_centerline_path)
        self.znormal = znormal

        # Compute convex hull for validation
        self.convex_hull = self.mesh.convex_hull

        # Compute plane normals for all centerline points
        print("Computing plane normals...")
        self.plane_normals = compute_all_plane_normals(
            self.centerline,
            smooth=True,
            znormal=znormal
        )

        # Initialize slicer
        self.mesh_slicer = MeshPlaneSlicer(self.mesh)

        print(f"Initialized with {len(self.centerline)} centerline points")

    def slice_along_centerline(self,
                              quality_checks: bool = True,
                              specific_indices: list = None) -> SlicingResults:
        """
        Slice mesh at each centerline point

        Replaces: Main loop in trachSlice_OSA_AB_3.m (lines 107-211)

        Args:
            quality_checks: Enable quality validation checks
            specific_indices: If provided, only slice at these indices

        Returns:
            SlicingResults object with all valid cross-sections
        """
        print(f"\nSlicing mesh along centerline...")
        print(f"Quality checks: {'enabled' if quality_checks else 'disabled'}")

        results = SlicingResults()
        n_points = len(self.centerline)

        # Determine which indices to process
        if specific_indices is not None:
            indices_to_process = specific_indices
            print(f"Using specific indices: {len(indices_to_process)} planes")
        else:
            indices_to_process = list(range(n_points))

        # Storage for tracking previous plane
        previous_position = None
        previous_normal = None

        for idx_pos, i in enumerate(indices_to_process):
            # Bounds check
            if i >= n_points:
                print(f"  Warning: Index {i} out of bounds, skipping")
                continue

            position = self.centerline[i]
            normal = self.plane_normals[i]  # Original normal from centerline

            print(f"  Plane {idx_pos+1}/{len(indices_to_process)} (index {i})...", end=" ")

            # ADAPTIVE NORMAL ADJUSTMENT: Prevent intersections dynamically
            if previous_position is not None and previous_normal is not None:
                # Compute forward direction from previous to current position
                forward_dir = position - previous_position
                forward_dist = np.linalg.norm(forward_dir)

                if forward_dist > 0.01:  # Avoid division by zero
                    forward_dir = forward_dir / forward_dist

                    # Check if planes would intersect with current normal
                    # by verifying current position is ahead of previous plane
                    signed_dist = np.dot(forward_dir, previous_normal)

                    # If angle is too shallow (< 30 degrees), adjust normal
                    if signed_dist < 0.5:  # cos(60°) = 0.5
                        # Blend normal with forward direction to ensure separation
                        # More blending for sharper curves (lower signed_dist)
                        blend_weight = max(0.3, signed_dist)  # 30-50% forward bias

                        adjusted_normal = blend_weight * normal + (1 - blend_weight) * forward_dir
                        adjusted_normal = adjusted_normal / np.linalg.norm(adjusted_normal)

                        normal = adjusted_normal
                        # Note: This adjustment happens silently - no skip message

            # Perform slicing
            cross_sections = self.mesh_slicer.slice_mesh_with_plane(
                    plane_origin=position,
                    plane_normal=normal,
                    plane_number=i
                )

            # Check if slicing failed
            if not cross_sections:
                print("SKIP (no intersection)")
                results.add_skipped_plane(i, "no_intersection")
                continue

            # BIFURCATION FIX: Merge all cross-sections (don't drop loops!)
            # At bifurcations, planes cut through multiple closed loops.
            # We must SUM all areas, not pick just one.
            # Replaces logic in trachSlice_OSA_AB_3.m lines 136-173
            merged_section = self._merge_cross_sections(cross_sections, position)

            if merged_section is None:
                print("SKIP (no valid section)")
                results.add_skipped_plane(i, "no_valid_section")
                continue

            # Quality check: Centroid proximity
            # Replaces lines 191-202 in trachSlice_OSA_AB_3.m
            if quality_checks:
                if not validate_centroid_proximity(
                    merged_section.centroid,
                    position,
                    threshold=10.0
                ):
                    print("SKIP (centroid too far)")
                    results.add_skipped_plane(i, "centroid_too_far")
                    continue

            # Quality check: Enclosed by mesh
            # New feature using convex hull
            if quality_checks:
                if not check_plane_enclosed_by_mesh(
                    merged_section.centroid,
                    self.convex_hull,
                    tolerance=1.0
                ):
                    print("SKIP (outside mesh)")
                    results.add_skipped_plane(i, "outside_mesh")
                    continue

            # Valid plane - add to results
            results.add_valid_plane(i, merged_section, position, normal)
            n_loops = len(cross_sections)
            if n_loops > 1:
                print(f"OK (area={merged_section.area:.2f} mm², {n_loops} loops merged)")
            else:
                print(f"OK (area={merged_section.area:.2f} mm²)")

            # Update tracking for next iteration
            previous_position = position
            previous_normal = normal

        # Compute arc lengths
        # Replaces lines 233-237 in trachSlice_OSA_AB_3.m
        results.compute_arc_lengths()

        # Convert units if needed (m to mm)
        # Replaces lines 239-243 in trachSlice_OSA_AB_3.m
        if len(results.arc_lengths) > 0:
            areas = np.array([s.area for s in results.valid_sections])
            arc_lengths_converted, areas_converted = convert_units_if_needed(
                results.arc_lengths,
                areas
            )
            results.arc_lengths = arc_lengths_converted

            # Update areas in sections
            for section, new_area in zip(results.valid_sections, areas_converted):
                try:
                    section.area = new_area
                except AttributeError:
                    # Trimesh area is read-only, skip for boundary cap meshes
                    pass

        # Print summary
        summary = results.get_summary()
        print(f"\nSlicing complete:")
        print(f"  Valid planes: {summary['n_valid_planes']}")
        print(f"  Skipped planes: {summary['n_skipped_planes']}")
        if summary['skip_reasons']:
            print(f"  Skip reasons:")
            for reason, count in summary['skip_reasons'].items():
                print(f"    - {reason}: {count}")

        return results

    def _merge_cross_sections(self,
                             cross_sections: List[CrossSection],
                             centerline_point: np.ndarray) -> Optional[CrossSection]:
        """
        Merge multiple cross-sections at bifurcations (SUM ALL AREAS!)

        BIFURCATION FIX: When a plane cuts through a bifurcation, it creates
        multiple closed loops. We must SUM the areas of ALL loops, not pick one.

        Args:
            cross_sections: List of cross-sections (each is a closed loop)
            centerline_point: Centerline point for centroid calculation

        Returns:
            Merged CrossSection with combined area, or None
        """
        if not cross_sections:
            return None

        if len(cross_sections) == 1:
            # Only one loop - return as-is
            return cross_sections[0]

        # BIFURCATION: Multiple closed loops detected!
        # Sum all areas and combine geometries

        total_area = sum(cs.area for cs in cross_sections)
        total_perimeter = sum(cs.perimeter for cs in cross_sections)

        # Combine all vertices and faces
        all_vertices = []
        all_faces = []
        vertex_offset = 0

        for cs in cross_sections:
            all_vertices.append(cs.vertices)
            # Offset face indices to account for combined vertex array
            offset_faces = cs.faces + vertex_offset
            all_faces.append(offset_faces)
            vertex_offset += len(cs.vertices)

        combined_vertices = np.vstack(all_vertices)
        combined_faces = np.vstack(all_faces)

        # Compute weighted centroid (weight by area)
        centroids = np.array([cs.centroid for cs in cross_sections])
        areas = np.array([cs.area for cs in cross_sections])
        weighted_centroid = np.average(centroids, axis=0, weights=areas)

        # Combine boundary 2D points (for visualization/measurements)
        all_boundary_2d = []
        for cs in cross_sections:
            if cs.boundary_2d is not None:
                all_boundary_2d.append(cs.boundary_2d)

        if all_boundary_2d:
            combined_boundary_2d = np.vstack(all_boundary_2d)
        else:
            combined_boundary_2d = cross_sections[0].boundary_2d

        # Use the side of the largest cross-section
        largest_cs = max(cross_sections, key=lambda cs: cs.area)

        # Create merged CrossSection
        return CrossSection(
            vertices=combined_vertices,
            faces=combined_faces,
            area=total_area,
            centroid=weighted_centroid,
            perimeter=total_perimeter,
            boundary_2d=combined_boundary_2d,
            side=largest_cs.side
        )

    def _select_best_section(self,
                            cross_sections: List[CrossSection],
                            centerline_point: np.ndarray) -> Optional[CrossSection]:
        """
        Select best cross-section from multiple candidates

        Replaces: Selection logic in trachSlice_OSA_AB_3.m (lines 136-173)

        MATLAB logic prioritizes:
        1. If only one section exists, use it
        2. If multiple sections, choose the one closest to centerline
        3. If centroids are equally close, choose larger area

        Args:
            cross_sections: List of candidate cross-sections
            centerline_point: Centerline point for comparison

        Returns:
            Best CrossSection or None
        """
        if not cross_sections:
            return None

        if len(cross_sections) == 1:
            return cross_sections[0]

        # Multiple sections - choose based on distance and area
        best_section = None
        best_distance = float('inf')
        best_area = 0.0

        for section in cross_sections:
            distance = np.linalg.norm(section.centroid - centerline_point)

            # If this is closer, it wins
            if distance < best_distance:
                best_section = section
                best_distance = distance
                best_area = section.area
            # If equally close (within 2mm), choose larger area
            # Replaces lines 150-155 in trachSlice_OSA_AB_3.m
            elif abs(distance - best_distance) < 2.0:
                if section.area > best_area:
                    best_section = section
                    best_distance = distance
                    best_area = section.area

        return best_section

    def export_results(self,
                      results: SlicingResults,
                      output_base_name: str,
                      output_dir: str = ".") -> None:
        """
        Export slicing results to files (ATOMIC: temp dir → final dir)

        Replaces: File writing in trachSlice_OSA_AB_3.m (lines 245-256)

        Ensures all-or-nothing export to prevent inconsistent file states if interrupted.

        Args:
            results: SlicingResults object
            output_base_name: Base name for output files (from STL filename)
            output_dir: Output directory
        """
        from pathlib import Path
        import shutil
        import tempfile
        from .io_utils import write_cross_section_stl, write_measurements_csv, write_all_planes_stl

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create temporary directory for atomic export
        # Files only appear in final location when ALL exports succeed
        temp_dir = Path(tempfile.mkdtemp(prefix=f"export_{output_base_name}_", dir=output_dir))

        try:
            # Export individual cross-section STLs to temp dir
            # Replaces lines 252-255 in MATLAB
            print(f"\nExporting {len(results.valid_sections)} cross-section STLs...")
            for i, section in enumerate(results.valid_sections):
                plane_idx = results.valid_plane_indices[i]
                stl_path = temp_dir / f"{output_base_name}-Planes-{plane_idx:03d}.stl"
                write_cross_section_stl(str(stl_path), section.vertices, section.faces)

            # Export combined STL with all planes to temp dir
            # Replaces line 251 in MATLAB
            all_vertices = [s.vertices for s in results.valid_sections]
            all_faces = [s.faces for s in results.valid_sections]
            combined_stl = temp_dir / f"{output_base_name}-Planes-All.stl"
            write_all_planes_stl(str(combined_stl), all_vertices, all_faces)

            # Export measurements CSV to temp dir
            # Replaces csvwrite in line 248 in MATLAB
            profile = results.compute_diameter_profile()
            df = profile.to_dataframe(flip_order=True)  # Flip to match MATLAB
            csv_path = temp_dir / f"{output_base_name}-Data.csv"
            df.to_csv(csv_path, index=False)

            # ATOMIC MOVE: All files written successfully, now move them
            # This prevents partial/inconsistent exports if process crashes mid-write
            for temp_file in temp_dir.iterdir():
                final_file = output_dir / temp_file.name
                # Remove old file if exists, then move new file
                if final_file.exists():
                    final_file.unlink()
                shutil.move(str(temp_file), str(final_file))

            print(f"Exported to: {output_dir}")
            print(f"  - {len(results.valid_sections)} individual plane STLs")
            print(f"  - 1 combined planes STL")
            print(f"  - 1 measurements CSV")

        finally:
            # Clean up temp directory
            if temp_dir.exists():
                shutil.rmtree(temp_dir)
