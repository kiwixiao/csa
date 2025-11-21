#!/usr/bin/env python
"""
Analyze first and last 10 planes across all time points using gap detection

This script validates whether open profile detection is working correctly
by checking if edge planes have proper mesh intersection geometry.
"""

import sys
import numpy as np
import trimesh
from pathlib import Path
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s',
    handlers=[
        logging.FileHandler('edge_planes_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def analyze_contour_edges(stl_path):
    """
    Analyze the edge geometry of a cross-section contour

    Returns:
        edge_lengths: array of edge lengths
        mean_length: mean edge length
        max_length: maximum edge length
        n_vertices: number of vertices
        area: cross-section area
        gap_score: 0-100 indicating likelihood of incomplete intersection
    """
    try:
        mesh = trimesh.load_mesh(str(stl_path))
    except Exception as e:
        return None, 0, 0, 0, 0, 100

    # Get boundary edges
    boundary_edges = mesh.edges[
        trimesh.grouping.group_rows(mesh.edges_sorted, require_count=1)
    ]

    if len(boundary_edges) == 0:
        return None, 0, 0, 0, mesh.area, 100

    # Calculate edge lengths
    edge_lengths = []
    for edge in boundary_edges:
        v1 = mesh.vertices[edge[0]]
        v2 = mesh.vertices[edge[1]]
        length = np.linalg.norm(v2 - v1)
        edge_lengths.append(length)

    edge_lengths = np.array(edge_lengths)
    mean_length = edge_lengths.mean()
    max_length = edge_lengths.max()

    # Calculate gap score
    max_ratio = max_length / mean_length if mean_length > 0 else 0
    cv = np.std(edge_lengths) / mean_length if mean_length > 0 else 0

    threshold = mean_length + 3 * np.std(edge_lengths)
    outlier_edges = np.where(edge_lengths > threshold)[0]
    outlier_percentage = len(outlier_edges) / len(edge_lengths) * 100

    gap_score = min(100, (
        (max_ratio - 1) * 10 +
        cv * 30 +
        outlier_percentage * 2
    ))

    return edge_lengths, mean_length, max_length, len(mesh.vertices), mesh.area, gap_score


def classify_plane(gap_score, max_edge, mean_edge):
    """
    Classify plane as complete, incomplete, or suspicious

    Returns:
        status: "COMPLETE", "INCOMPLETE", or "SUSPICIOUS"
        reason: explanation
    """
    if gap_score > 80 or (max_edge > 10 * mean_edge):
        return "INCOMPLETE", f"gap_score={gap_score:.1f}, max_edge={max_edge:.2f}mm (>{10}×mean)"
    elif gap_score > 50:
        return "SUSPICIOUS", f"gap_score={gap_score:.1f}"
    else:
        return "COMPLETE", f"gap_score={gap_score:.1f}"


def analyze_time_point(stl_dir, time_point_name, first_n=10, last_n=10, total_planes=148):
    """
    Analyze first and last N planes for a time point

    Returns:
        results: dict with analysis results
    """
    stl_dir = Path(stl_dir)

    # Analyze first N planes
    first_planes = []
    for idx in range(first_n):
        stl_file = stl_dir / f"{time_point_name}-Planes-{idx:03d}.stl"

        if not stl_file.exists():
            continue

        edge_lengths, mean_length, max_length, n_vertices, area, gap_score = analyze_contour_edges(stl_file)

        if edge_lengths is None:
            continue

        status, reason = classify_plane(gap_score, max_length, mean_length)

        first_planes.append({
            'idx': idx,
            'area': area,
            'n_vertices': n_vertices,
            'mean_edge': mean_length,
            'max_edge': max_length,
            'gap_score': gap_score,
            'status': status,
            'reason': reason
        })

    # Analyze last N planes
    last_planes = []
    for idx in range(total_planes - last_n, total_planes):
        stl_file = stl_dir / f"{time_point_name}-Planes-{idx:03d}.stl"

        if not stl_file.exists():
            continue

        edge_lengths, mean_length, max_length, n_vertices, area, gap_score = analyze_contour_edges(stl_file)

        if edge_lengths is None:
            continue

        status, reason = classify_plane(gap_score, max_length, mean_length)

        last_planes.append({
            'idx': idx,
            'area': area,
            'n_vertices': n_vertices,
            'mean_edge': mean_length,
            'max_edge': max_length,
            'gap_score': gap_score,
            'status': status,
            'reason': reason
        })

    return {
        'time_point': time_point_name,
        'first_planes': first_planes,
        'last_planes': last_planes
    }


def print_plane_summary(planes, label):
    """Print summary of planes"""
    incomplete = [p for p in planes if p['status'] == 'INCOMPLETE']
    suspicious = [p for p in planes if p['status'] == 'SUSPICIOUS']
    complete = [p for p in planes if p['status'] == 'COMPLETE']

    logger.info(f"\n  {label}:")
    logger.info(f"    Total: {len(planes)}")
    logger.info(f"    Complete: {len(complete)}")
    logger.info(f"    Suspicious: {len(suspicious)}")
    logger.info(f"    Incomplete: {len(incomplete)}")

    if incomplete:
        logger.info(f"    Incomplete indices: {[p['idx'] for p in incomplete]}")
        for p in incomplete:
            logger.info(f"      Plane {p['idx']:3d}: {p['reason']}")

    if suspicious:
        logger.info(f"    Suspicious indices: {[p['idx'] for p in suspicious]}")


def main():
    """Main analysis"""
    logger.info("="*70)
    logger.info("EDGE PLANES ANALYSIS - ALL TIME POINTS")
    logger.info("="*70)
    logger.info("Analyzing first 10 and last 10 planes for all time points")
    logger.info("Using gap detection to identify incomplete intersections")
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

    logger.info(f"\nFound {len(time_points)} time points\n")

    # Analyze each time point
    all_results = []

    for i, time_point in enumerate(time_points):
        logger.info(f"[{i+1}/{len(time_points)}] {time_point}")

        results = analyze_time_point(stl_dir, time_point, first_n=10, last_n=10)
        all_results.append(results)

        # Print summary for this time point
        print_plane_summary(results['first_planes'], "First 10 planes")
        print_plane_summary(results['last_planes'], "Last 10 planes")

    # Cross-time-point summary
    logger.info("\n" + "="*70)
    logger.info("CROSS-TIME-POINT SUMMARY")
    logger.info("="*70)

    # Count time points with issues
    timepoints_with_incomplete_first = 0
    timepoints_with_incomplete_last = 0
    timepoints_with_any_incomplete = 0

    all_incomplete_first_indices = set()
    all_incomplete_last_indices = set()

    for results in all_results:
        first_incomplete = [p for p in results['first_planes'] if p['status'] == 'INCOMPLETE']
        last_incomplete = [p for p in results['last_planes'] if p['status'] == 'INCOMPLETE']

        if first_incomplete:
            timepoints_with_incomplete_first += 1
            all_incomplete_first_indices.update([p['idx'] for p in first_incomplete])

        if last_incomplete:
            timepoints_with_incomplete_last += 1
            all_incomplete_last_indices.update([p['idx'] for p in last_incomplete])

        if first_incomplete or last_incomplete:
            timepoints_with_any_incomplete += 1

    logger.info(f"\nTime points with incomplete planes:")
    logger.info(f"  - First 10 planes: {timepoints_with_incomplete_first}/{len(time_points)} time points")
    logger.info(f"  - Last 10 planes: {timepoints_with_incomplete_last}/{len(time_points)} time points")
    logger.info(f"  - Any incomplete: {timepoints_with_any_incomplete}/{len(time_points)} time points")

    logger.info(f"\nIncomplete plane indices (across all time points):")
    if all_incomplete_first_indices:
        logger.info(f"  - First 10: {sorted(all_incomplete_first_indices)}")
    else:
        logger.info(f"  - First 10: None")

    if all_incomplete_last_indices:
        logger.info(f"  - Last 10: {sorted(all_incomplete_last_indices)}")
    else:
        logger.info(f"  - Last 10: None")

    # Identify consistently incomplete planes
    logger.info("\n" + "="*70)
    logger.info("CONSISTENTLY INCOMPLETE PLANES")
    logger.info("="*70)

    # Count occurrences of each incomplete plane index
    first_incomplete_counts = {}
    last_incomplete_counts = {}

    for results in all_results:
        for p in results['first_planes']:
            if p['status'] == 'INCOMPLETE':
                idx = p['idx']
                first_incomplete_counts[idx] = first_incomplete_counts.get(idx, 0) + 1

        for p in results['last_planes']:
            if p['status'] == 'INCOMPLETE':
                idx = p['idx']
                last_incomplete_counts[idx] = last_incomplete_counts.get(idx, 0) + 1

    # Report planes incomplete in majority of time points
    threshold = len(time_points) * 0.7  # 70% threshold

    consistently_incomplete_first = {idx: count for idx, count in first_incomplete_counts.items()
                                     if count >= threshold}
    consistently_incomplete_last = {idx: count for idx, count in last_incomplete_counts.items()
                                    if count >= threshold}

    if consistently_incomplete_first:
        logger.info(f"\nFirst 10 planes incomplete in ≥70% of time points:")
        for idx, count in sorted(consistently_incomplete_first.items()):
            logger.info(f"  Plane {idx:3d}: {count}/{len(time_points)} time points ({count/len(time_points)*100:.1f}%)")
    else:
        logger.info(f"\n✓ No consistently incomplete planes in first 10")

    if consistently_incomplete_last:
        logger.info(f"\nLast 10 planes incomplete in ≥70% of time points:")
        for idx, count in sorted(consistently_incomplete_last.items()):
            logger.info(f"  Plane {idx:3d}: {count}/{len(time_points)} time points ({count/len(time_points)*100:.1f}%)")
    else:
        logger.info(f"\n✓ No consistently incomplete planes in last 10")

    # Final verdict
    logger.info("\n" + "="*70)
    logger.info("FINAL VERDICT")
    logger.info("="*70)

    if timepoints_with_any_incomplete == 0:
        logger.info("✓ ALL edge planes are complete across all time points!")
        logger.info("✓ Open profile detection appears to be working correctly")
    elif timepoints_with_any_incomplete < len(time_points) * 0.3:
        logger.info(f"⚠️ Some edge planes incomplete ({timepoints_with_any_incomplete}/{len(time_points)} time points)")
        logger.info("  - May need adjustment to detection thresholds")
    else:
        logger.info(f"❌ Many edge planes incomplete ({timepoints_with_any_incomplete}/{len(time_points)} time points)")
        logger.info("  - Open profile detection may need significant fixes")

    logger.info("\n✓ Analysis complete! Check 'edge_planes_analysis.log' for full details")


if __name__ == "__main__":
    main()
