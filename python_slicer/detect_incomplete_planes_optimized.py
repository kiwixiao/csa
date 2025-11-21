#!/usr/bin/env python
"""
Optimized Two-Stage Incomplete Plane Detection

Stage 1 (Fast): Edge length analysis to identify suspicious planes
Stage 2 (Accurate): Surface validation to eliminate false positives

This approach saves computation time by only running expensive surface
validation on planes that are flagged as suspicious in the first stage.
"""

import numpy as np
import trimesh
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def check_edge_on_surface(edge_start, edge_end, source_mesh, n_samples=10, tolerance=0.5):
    """
    Check if a contour edge lies on the source mesh surface

    Args:
        edge_start, edge_end: 3D coordinates of edge endpoints
        source_mesh: The original mesh being sliced
        n_samples: Number of points to sample along edge
        tolerance: Maximum distance from surface (mm)

    Returns:
        is_on_surface: True if edge lies on surface
        max_distance: Maximum distance from surface
        off_surface_ratio: Fraction of samples off surface
    """
    # Sample points along the edge
    t_values = np.linspace(0, 1, n_samples)
    edge_samples = np.array([
        edge_start + t * (edge_end - edge_start)
        for t in t_values
    ])

    # Check distance to surface for each sample
    closest_points, distances, face_ids = source_mesh.nearest.on_surface(edge_samples)

    # Count how many samples are off surface
    off_surface = distances > tolerance
    off_surface_ratio = off_surface.sum() / len(distances)
    max_distance = distances.max()

    # Edge is "on surface" if most samples are close to surface
    is_on_surface = off_surface_ratio < 0.3  # Allow 30% tolerance

    return is_on_surface, max_distance, off_surface_ratio


def stage1_edge_length_check(cross_section_stl, long_edge_threshold=2.0):
    """
    Stage 1 (Fast): Check for unusually long edges

    Returns:
        has_long_edges: True if plane has suspiciously long edges
        long_edge_info: List of long edge details for stage 2
    """
    try:
        cross_section = trimesh.load_mesh(str(cross_section_stl))
    except Exception as e:
        return False, []

    # Get boundary edges
    boundary_edges = cross_section.edges[
        trimesh.grouping.group_rows(cross_section.edges_sorted, require_count=1)
    ]

    if len(boundary_edges) == 0:
        return False, []

    # Calculate edge lengths
    edge_data = []
    for edge in boundary_edges:
        v1 = cross_section.vertices[edge[0]]
        v2 = cross_section.vertices[edge[1]]
        length = np.linalg.norm(v2 - v1)
        edge_data.append({'v1': v1, 'v2': v2, 'length': length})

    edge_lengths = np.array([e['length'] for e in edge_data])
    mean_length = edge_lengths.mean()

    # Find long edges
    long_edges = []
    for data in edge_data:
        if data['length'] > long_edge_threshold * mean_length:
            long_edges.append({
                'v1': data['v1'],
                'v2': data['v2'],
                'length': data['length'],
                'length_ratio': data['length'] / mean_length
            })

    return len(long_edges) > 0, long_edges


def stage2_surface_validation(long_edges, source_mesh_stl):
    """
    Stage 2 (Accurate): Validate long edges against source mesh surface

    Returns:
        has_artificial_gaps: True if any edges are off surface
        artificial_gap_count: Number of artificial gaps found
    """
    try:
        source_mesh = trimesh.load_mesh(str(source_mesh_stl))
    except Exception as e:
        return False, 0

    artificial_gaps = 0

    for edge in long_edges:
        is_on_surface, max_dist, off_ratio = check_edge_on_surface(
            edge['v1'], edge['v2'], source_mesh
        )

        if not is_on_surface:
            artificial_gaps += 1

    return artificial_gaps > 0, artificial_gaps


def detect_incomplete_plane_optimized(cross_section_stl, source_mesh_stl,
                                      long_edge_threshold=2.0):
    """
    Two-stage optimized detection of incomplete planes

    Stage 1: Fast edge length check
    Stage 2: Surface validation (only if stage 1 flags suspicious edges)

    Returns:
        is_incomplete: True if plane is incomplete
        reason: Description of why plane is incomplete (or "complete")
    """
    # Stage 1: Fast edge length check
    has_long_edges, long_edges = stage1_edge_length_check(
        cross_section_stl,
        long_edge_threshold=long_edge_threshold
    )

    if not has_long_edges:
        # No long edges - plane is complete
        return False, "complete (no long edges)"

    # Stage 2: Surface validation (only for planes with long edges)
    has_artificial_gaps, gap_count = stage2_surface_validation(
        long_edges,
        source_mesh_stl
    )

    if has_artificial_gaps:
        return True, f"incomplete ({gap_count} artificial gap(s), {len(long_edges)} long edges)"
    else:
        return False, f"complete ({len(long_edges)} long edges on surface)"


def detect_incomplete_planes_for_timepoint(temp_dir, source_mesh_path,
                                           time_point, n_planes=148):
    """
    Detect incomplete planes for a single time point using optimized two-stage method

    Returns:
        incomplete_indices: List of plane indices that are incomplete
        stage1_count: Number of planes flagged in stage 1 (had long edges)
        stage2_count: Number of planes confirmed incomplete in stage 2
    """
    temp_dir = Path(temp_dir)
    incomplete_indices = []
    stage1_flagged = 0
    stage2_confirmed = 0

    for plane_idx in range(n_planes):
        cross_section_path = temp_dir / f"{time_point}-Planes-{plane_idx:03d}.stl"

        if not cross_section_path.exists():
            continue

        is_incomplete, reason = detect_incomplete_plane_optimized(
            cross_section_path,
            source_mesh_path,
            long_edge_threshold=2.0
        )

        # Track statistics
        if "long edges" in reason:
            stage1_flagged += 1
        if is_incomplete:
            stage2_confirmed += 1
            incomplete_indices.append(plane_idx)

    return incomplete_indices, stage1_flagged, stage2_confirmed


def main():
    """Test the optimized two-stage detection"""
    logger.info("="*70)
    logger.info("OPTIMIZED TWO-STAGE INCOMPLETE PLANE DETECTION")
    logger.info("="*70)
    logger.info("Stage 1: Fast edge length check")
    logger.info("Stage 2: Surface validation (only for suspicious planes)")
    logger.info("="*70)
    logger.info("")

    # Configuration
    temp_dir = Path("LeftNoseDecendingSlicedSTLs_temp")
    source_mesh_dir = Path("LeftNoseDecending/FFD/stl")

    if not temp_dir.exists():
        logger.error(f"Directory not found: {temp_dir}")
        return

    # Get all time points
    time_point_files = sorted(temp_dir.glob("out_*-Planes-000.stl"))
    time_points = [f.stem.replace("-Planes-000", "") for f in time_point_files]

    logger.info(f"Found {len(time_points)} time points\n")

    # Track results
    all_incomplete = {}
    total_stage1 = 0
    total_stage2 = 0

    for i, time_point in enumerate(time_points):
        logger.info(f"[{i+1}/{len(time_points)}] {time_point}...")

        source_mesh_path = source_mesh_dir / f"{time_point}.stl"
        if not source_mesh_path.exists():
            logger.warning(f"  Source mesh not found")
            continue

        incomplete, stage1, stage2 = detect_incomplete_planes_for_timepoint(
            temp_dir, source_mesh_path, time_point
        )

        all_incomplete[time_point] = incomplete
        total_stage1 += stage1
        total_stage2 += stage2

        logger.info(f"  Stage 1 flagged: {stage1} planes")
        logger.info(f"  Stage 2 confirmed: {stage2} planes")
        if incomplete:
            logger.info(f"  Incomplete: {incomplete}")
        logger.info("")

    # Find consistently incomplete planes
    logger.info("="*70)
    logger.info("CROSS-TIME-POINT SUMMARY")
    logger.info("="*70)
    logger.info(f"Total planes checked in stage 1: ~{total_stage1}")
    logger.info(f"Total planes confirmed incomplete: ~{total_stage2}")
    logger.info(f"Computation saved: {total_stage1 - total_stage2} surface validations skipped")
    logger.info("")

    # Count plane occurrences
    plane_counts = {}
    for incomplete_list in all_incomplete.values():
        for idx in incomplete_list:
            plane_counts[idx] = plane_counts.get(idx, 0) + 1

    # Consistently incomplete (>=70% of time points)
    n_timepoints = len(time_points)
    threshold = int(n_timepoints * 0.7)

    consistently_incomplete = {idx: count for idx, count in plane_counts.items()
                              if count >= threshold}

    logger.info(f"Planes incomplete in ≥70% of time points:")
    if consistently_incomplete:
        logger.info(f"  Count: {len(consistently_incomplete)}")
        for idx, count in sorted(consistently_incomplete.items()):
            logger.info(f"  Plane {idx:3d}: {count}/{n_timepoints} time points ({count/n_timepoints*100:.1f}%)")
    else:
        logger.info("  None!")

    logger.info("")
    logger.info("✓ Optimized detection complete!")


if __name__ == "__main__":
    main()
