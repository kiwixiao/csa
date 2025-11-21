#!/usr/bin/env python
"""
Analyze contour geometry to detect incomplete plane intersections

A proper plane cut should have:
1. All contour edges correspond to actual mesh face intersections
2. Contour point density consistent with mesh resolution
3. No unusually long "gap spanning" edges

This script analyzes edge lengths and point density to detect incomplete cuts.
"""

import sys
import numpy as np
import trimesh
from pathlib import Path
import matplotlib.pyplot as plt


def analyze_contour_edges(stl_path):
    """
    Analyze the edge geometry of a cross-section contour

    Returns:
        edge_lengths: array of edge lengths
        mean_length: mean edge length
        max_length: maximum edge length
        std_length: standard deviation of edge lengths
        n_long_edges: number of edges > 2*mean (potential gaps)
    """
    mesh = trimesh.load_mesh(str(stl_path))

    # Get boundary edges (the contour)
    boundary_edges = mesh.edges[
        trimesh.grouping.group_rows(mesh.edges_sorted, require_count=1)
    ]

    if len(boundary_edges) == 0:
        return None, 0, 0, 0, 0

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
    std_length = edge_lengths.std()

    # Count "long edges" (potential gaps)
    # These might be artificial closures rather than real mesh intersections
    n_long_edges = np.sum(edge_lengths > 2 * mean_length)

    return edge_lengths, mean_length, max_length, std_length, n_long_edges


def calculate_gap_score(edge_lengths, mean_length):
    """
    Calculate a "gap score" indicating likelihood of incomplete intersection

    High gap score = likely incomplete (has artificial closures)
    Low gap score = likely complete (all edges from mesh faces)

    Returns:
        gap_score: 0-100 (higher = more likely incomplete)
        outlier_edges: indices of edges that are suspiciously long
    """
    if edge_lengths is None or len(edge_lengths) == 0:
        return 100, []

    # Identify outlier edges (> 3 standard deviations)
    threshold = mean_length + 3 * np.std(edge_lengths)
    outlier_edges = np.where(edge_lengths > threshold)[0]

    # Calculate gap score based on:
    # 1. Presence of very long edges (>> mean)
    # 2. High variance in edge lengths
    # 3. Percentage of outlier edges

    max_ratio = edge_lengths.max() / mean_length if mean_length > 0 else 0
    cv = np.std(edge_lengths) / mean_length if mean_length > 0 else 0  # Coefficient of variation
    outlier_percentage = len(outlier_edges) / len(edge_lengths) * 100

    # Combine metrics into gap score
    # Max ratio > 5: definitely suspicious
    # CV > 1: high variance indicates gaps
    # Outlier percentage > 5%: multiple gaps

    gap_score = min(100, (
        (max_ratio - 1) * 10 +  # Penalize extreme length ratios
        cv * 30 +                # Penalize high variance
        outlier_percentage * 2    # Penalize many outliers
    ))

    return gap_score, outlier_edges.tolist()


def analyze_plane_comparison(stl_dir, time_point, plane_indices):
    """
    Compare multiple planes to identify anomalies

    Args:
        stl_dir: directory containing STL files
        time_point: time point name
        plane_indices: list of plane indices to compare

    Returns:
        results: dict of {plane_idx: analysis_results}
    """
    results = {}

    for idx in plane_indices:
        stl_file = Path(stl_dir) / f"{time_point}-Planes-{idx:03d}.stl"

        if not stl_file.exists():
            print(f"Plane {idx}: MISSING")
            continue

        mesh = trimesh.load_mesh(str(stl_file))
        edge_lengths, mean_length, max_length, std_length, n_long_edges = analyze_contour_edges(stl_file)

        if edge_lengths is None:
            print(f"Plane {idx}: NO EDGES")
            continue

        gap_score, outlier_edges = calculate_gap_score(edge_lengths, mean_length)

        results[idx] = {
            'area': mesh.area,
            'n_vertices': len(mesh.vertices),
            'n_edges': len(edge_lengths),
            'edge_lengths': edge_lengths,
            'mean_edge_length': mean_length,
            'max_edge_length': max_length,
            'std_edge_length': std_length,
            'n_long_edges': n_long_edges,
            'gap_score': gap_score,
            'outlier_edges': outlier_edges,
            'max_to_mean_ratio': max_length / mean_length if mean_length > 0 else 0
        }

        # Print summary
        status = "⚠️ SUSPICIOUS" if gap_score > 30 else "✓ OK"
        print(f"Plane {idx:3d}: {status}")
        print(f"  Area: {mesh.area:7.2f} mm²")
        print(f"  Vertices: {len(mesh.vertices):4d}")
        print(f"  Edges: {len(edge_lengths):4d}")
        print(f"  Mean edge length: {mean_length:.3f} mm")
        print(f"  Max edge length: {max_length:.3f} mm (ratio: {max_length/mean_length:.1f}x)")
        print(f"  Std dev: {std_length:.3f} mm")
        print(f"  Long edges (>2*mean): {n_long_edges}")
        print(f"  Gap score: {gap_score:.1f}")
        if outlier_edges:
            print(f"  Outlier edges: {len(outlier_edges)} ({len(outlier_edges)/len(edge_lengths)*100:.1f}%)")
        print()

    return results


def visualize_edge_distributions(results, output_file="edge_distribution_comparison.png"):
    """
    Create visualization comparing edge length distributions
    """
    n_planes = len(results)
    fig, axes = plt.subplots(1, n_planes, figsize=(6*n_planes, 5))

    if n_planes == 1:
        axes = [axes]

    for ax, (idx, data) in zip(axes, results.items()):
        edge_lengths = data['edge_lengths']
        mean_length = data['mean_edge_length']
        gap_score = data['gap_score']

        ax.hist(edge_lengths, bins=30, alpha=0.7, edgecolor='black')
        ax.axvline(mean_length, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_length:.3f}')
        ax.axvline(2*mean_length, color='orange', linestyle='--', linewidth=2, label=f'2×Mean: {2*mean_length:.3f}')

        status = "⚠️ INCOMPLETE" if gap_score > 30 else "✓ Complete"
        ax.set_title(f'Plane {idx} {status}\nGap Score: {gap_score:.1f}')
        ax.set_xlabel('Edge Length (mm)')
        ax.set_ylabel('Count')
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Saved edge distribution plot to: {output_file}")


def main():
    """Analyze plane 76 and neighbors to demonstrate incomplete intersection"""
    print("="*70)
    print("CONTOUR GEOMETRY ANALYSIS")
    print("="*70)
    print("Purpose: Detect incomplete plane intersections")
    print("Method: Analyze edge lengths and point density")
    print("="*70)
    print()

    # Configuration
    stl_dir = Path("LeftNoseDecendingSlicedSTLs_temp")
    time_point = "out_001400"

    # Analyze plane 76 and its neighbors
    print(f"Analyzing {time_point} (frame 15)")
    print(f"Comparing planes 75, 76, 77")
    print()

    results = analyze_plane_comparison(stl_dir, time_point, [75, 76, 77])

    # Compare results
    print("="*70)
    print("COMPARISON SUMMARY")
    print("="*70)

    if 76 in results and 75 in results and 77 in results:
        plane75 = results[75]
        plane76 = results[76]
        plane77 = results[77]

        print("\nArea comparison:")
        print(f"  Plane 75: {plane75['area']:.2f} mm²")
        print(f"  Plane 76: {plane76['area']:.2f} mm² ({plane76['area']/plane75['area']*100:.1f}% of plane 75)")
        print(f"  Plane 77: {plane77['area']:.2f} mm²")

        print("\nVertex count comparison:")
        print(f"  Plane 75: {plane75['n_vertices']} vertices")
        print(f"  Plane 76: {plane76['n_vertices']} vertices ({plane76['n_vertices']/plane75['n_vertices']*100:.1f}% of plane 75)")
        print(f"  Plane 77: {plane77['n_vertices']} vertices")

        print("\nEdge length comparison:")
        print(f"  Plane 75: mean={plane75['mean_edge_length']:.3f} mm, max={plane75['max_edge_length']:.3f} mm")
        print(f"  Plane 76: mean={plane76['mean_edge_length']:.3f} mm, max={plane76['max_edge_length']:.3f} mm")
        print(f"  Plane 77: mean={plane77['mean_edge_length']:.3f} mm, max={plane77['max_edge_length']:.3f} mm")

        print("\nGap score comparison:")
        print(f"  Plane 75: {plane75['gap_score']:.1f} ({'SUSPICIOUS' if plane75['gap_score'] > 30 else 'OK'})")
        print(f"  Plane 76: {plane76['gap_score']:.1f} ({'SUSPICIOUS' if plane76['gap_score'] > 30 else 'OK'})")
        print(f"  Plane 77: {plane77['gap_score']:.1f} ({'SUSPICIOUS' if plane77['gap_score'] > 30 else 'OK'})")

        # Visualize
        visualize_edge_distributions(results)

        print("\n" + "="*70)
        print("CONCLUSION")
        print("="*70)

        if plane76['gap_score'] > 30 and plane75['gap_score'] < 30 and plane77['gap_score'] < 30:
            print("✓ Plane 76 is CONFIRMED INCOMPLETE:")
            print("  - Gap score significantly higher than neighbors")
            print("  - Fewer vertices indicating missing mesh intersection data")
            print("  - Likely caused by plane intersecting mesh opening")
        elif plane76['gap_score'] > 30:
            print("⚠️ Plane 76 may be incomplete (high gap score)")
        else:
            print("? All planes appear geometrically valid")
            print("  (but visual inspection suggests plane 76 is incomplete)")

    print("\n✓ Analysis complete!")


if __name__ == "__main__":
    main()
