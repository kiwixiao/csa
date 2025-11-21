#!/usr/bin/env python
"""
Debug script to detect incomplete planes (boundary edges)

A plane with boundary edges means the cross-section hits a mesh opening
and is incomplete - NOT all contour points are on the surface.

This script:
1. Checks ALL sliced planes for boundary edges
2. Logs detailed statistics
3. Identifies problem indices
4. Validates against expected closed planes
"""

import sys
import logging
from pathlib import Path
import trimesh
import numpy as np

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('boundary_edges_debug.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def check_plane_boundary_edges(stl_path):
    """
    Check if a plane has boundary edges (incomplete cross-section)

    Returns:
        is_complete (bool): True if no boundary edges (fully enclosed)
        n_boundary_edges (int): Number of open boundary edges
        area (float): Cross-section area
        n_vertices (int): Number of vertices
        n_faces (int): Number of faces
    """
    try:
        mesh = trimesh.load_mesh(str(stl_path))

        # Find boundary edges (edges that appear only once, not shared by 2 faces)
        boundary_edges = mesh.edges[
            trimesh.grouping.group_rows(mesh.edges_sorted, require_count=1)
        ]

        is_complete = len(boundary_edges) == 0 and mesh.is_watertight

        return is_complete, len(boundary_edges), mesh.area, len(mesh.vertices), len(mesh.faces)

    except Exception as e:
        logger.error(f"Error checking {stl_path.name}: {e}")
        return None, 0, 0, 0, 0


def analyze_time_point(time_point_name, stl_dir, n_expected_planes=148):
    """
    Analyze all planes for a single time point

    Returns:
        incomplete_indices: List of plane indices with boundary edges
        stats: Dictionary of statistics
    """
    stl_dir = Path(stl_dir)

    logger.info(f"\n{'='*70}")
    logger.info(f"Analyzing: {time_point_name}")
    logger.info(f"{'='*70}")

    incomplete_indices = []
    complete_indices = []

    total_planes = 0
    total_area = 0

    # Check each plane
    for plane_idx in range(n_expected_planes):
        stl_file = stl_dir / f"{time_point_name}-Planes-{plane_idx:03d}.stl"

        if not stl_file.exists():
            logger.warning(f"  Plane {plane_idx:3d}: MISSING")
            continue

        is_complete, n_boundary, area, n_verts, n_faces = check_plane_boundary_edges(stl_file)

        if is_complete is None:
            continue

        total_planes += 1
        total_area += area

        if not is_complete:
            incomplete_indices.append(plane_idx)
            logger.warning(
                f"  Plane {plane_idx:3d}: INCOMPLETE - "
                f"{n_boundary} boundary edges, area={area:.2f} mm², "
                f"verts={n_verts}, faces={n_faces}"
            )
        else:
            complete_indices.append(plane_idx)
            logger.debug(
                f"  Plane {plane_idx:3d}: OK - "
                f"area={area:.2f} mm², verts={n_verts}, faces={n_faces}"
            )

    # Summary statistics
    stats = {
        'total_planes': total_planes,
        'complete_planes': len(complete_indices),
        'incomplete_planes': len(incomplete_indices),
        'completion_rate': len(complete_indices) / total_planes * 100 if total_planes > 0 else 0,
        'total_area': total_area,
        'incomplete_indices': incomplete_indices,
        'complete_indices': complete_indices
    }

    logger.info(f"\nSUMMARY for {time_point_name}:")
    logger.info(f"  Total planes checked: {stats['total_planes']}")
    logger.info(f"  Complete (no boundary edges): {stats['complete_planes']}")
    logger.info(f"  Incomplete (has boundary edges): {stats['incomplete_planes']}")
    logger.info(f"  Completion rate: {stats['completion_rate']:.1f}%")

    if incomplete_indices:
        logger.info(f"  Incomplete plane indices: {incomplete_indices[:20]}")
        if len(incomplete_indices) > 20:
            logger.info(f"    ... and {len(incomplete_indices) - 20} more")

    return incomplete_indices, stats


def find_consistently_incomplete_planes(all_time_points_data):
    """
    Find planes that are incomplete across ALL time points

    Args:
        all_time_points_data: Dict of {time_point_name: (incomplete_indices, stats)}

    Returns:
        consistently_incomplete: Set of plane indices incomplete in ALL time points
        sometimes_incomplete: Set of plane indices incomplete in SOME time points
    """
    logger.info(f"\n{'='*70}")
    logger.info("CROSS-TIME-POINT ANALYSIS")
    logger.info(f"{'='*70}")

    # Get all plane indices
    all_incomplete = set()
    for incomplete_indices, _ in all_time_points_data.values():
        all_incomplete.update(incomplete_indices)

    # Check which ones are incomplete in ALL time points
    consistently_incomplete = set(all_incomplete)
    for incomplete_indices, _ in all_time_points_data.values():
        consistently_incomplete &= set(incomplete_indices)

    sometimes_incomplete = all_incomplete - consistently_incomplete

    logger.info(f"\nPlanes incomplete in ALL {len(all_time_points_data)} time points: {len(consistently_incomplete)}")
    if consistently_incomplete:
        logger.info(f"  Indices: {sorted(consistently_incomplete)}")

    logger.info(f"\nPlanes incomplete in SOME time points: {len(sometimes_incomplete)}")
    if sometimes_incomplete:
        logger.info(f"  Indices: {sorted(sometimes_incomplete)[:30]}")
        if len(sometimes_incomplete) > 30:
            logger.info(f"    ... and {len(sometimes_incomplete) - 30} more")

    return consistently_incomplete, sometimes_incomplete


def main():
    """Main debug script"""
    logger.info("="*70)
    logger.info("BOUNDARY EDGE DETECTION DEBUG")
    logger.info("="*70)
    logger.info("Purpose: Identify incomplete planes (cross-sections with boundary edges)")
    logger.info("Problem: When planes cut through mesh openings, trimesh creates incomplete geometry")
    logger.info("="*70)

    # Configuration
    stl_dir = Path("LeftNoseDecendingSlicedSTLs_temp")

    if not stl_dir.exists():
        logger.error(f"Directory not found: {stl_dir}")
        logger.error("Run the slicer first to generate STL files")
        return

    # Get all time points
    time_point_files = sorted(stl_dir.glob("out_*-Planes-000.stl"))
    time_points = [f.stem.replace("-Planes-000", "") for f in time_point_files]

    logger.info(f"\nFound {len(time_points)} time points")

    # Analyze each time point
    all_time_points_data = {}

    for time_point in time_points:
        incomplete_indices, stats = analyze_time_point(time_point, stl_dir)
        all_time_points_data[time_point] = (incomplete_indices, stats)

    # Cross-time-point analysis
    if len(all_time_points_data) > 1:
        consistently_incomplete, sometimes_incomplete = find_consistently_incomplete_planes(all_time_points_data)

        # Write validated indices (complete in ALL time points)
        all_planes = set(range(148))
        all_incomplete = consistently_incomplete | sometimes_incomplete
        validated_complete = sorted(all_planes - all_incomplete)

        output_file = Path("LeftNoseDecending_boundary_validated_indices.txt")
        with open(output_file, 'w') as f:
            for idx in validated_complete:
                f.write(f"{idx}\n")

        logger.info(f"\n{'='*70}")
        logger.info("FINAL VALIDATED PLANES")
        logger.info(f"{'='*70}")
        logger.info(f"Complete in ALL time points: {len(validated_complete)}")
        logger.info(f"Saved to: {output_file}")
        logger.info(f"{'='*70}")

    logger.info("\nDebug complete! Check 'boundary_edges_debug.log' for full details")


if __name__ == "__main__":
    main()
