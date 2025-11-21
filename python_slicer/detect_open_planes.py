#!/usr/bin/env python
"""
Detect open vs closed planes (half-planes vs fully enclosed)

An open plane cuts through an airway opening (like nostril entrance)
A closed plane is fully enclosed within the airway geometry
"""

import trimesh
import numpy as np
import pandas as pd
from pathlib import Path
import sys

def check_plane_closure(mesh_path, centerline, plane_index, edge_threshold=10):
    """
    Check if a plane intersection is open or closed

    An "open" plane cuts through an airway opening (nostril entrance/exit).
    Planes in the middle with multiple contours (bifurcations, donuts) are still "closed".

    Strategy:
    - Only check for open planes near beginning/end of centerline (spatial filter)
    - Middle planes are assumed closed even if they have multiple contours

    Args:
        edge_threshold: Only check first/last N indices for openings

    Returns:
        is_closed (bool): True if fully enclosed, False if open
        n_boundary_edges (int): Number of open edges
        n_paths (int): Number of contour paths
    """
    # Load airway mesh
    mesh = trimesh.load_mesh(mesh_path)

    # Get plane from centerline
    plane_origin = centerline[plane_index]

    # Compute plane normal (simplified - just use direction to next point)
    if plane_index < len(centerline) - 1:
        direction = centerline[plane_index + 1] - centerline[plane_index]
    else:
        direction = centerline[plane_index] - centerline[plane_index - 1]

    plane_normal = direction / np.linalg.norm(direction)

    # Intersect mesh with plane
    try:
        slice_result = mesh.section(plane_origin=plane_origin,
                                    plane_normal=plane_normal)

        if slice_result is None:
            return None, 0, 0

        # Convert to 2D paths
        paths_2d, transform = slice_result.to_planar()

        if paths_2d is None:
            return None, 0, 0

        # Check each path entity
        n_paths = len(paths_2d.entities)
        total_boundary_edges = 0
        all_closed = True

        # SPATIAL FILTER: Only check for open contours at beginning/end
        # Middle planes with multiple contours (bifurcations, donuts) are OK
        is_near_edge = (plane_index < edge_threshold or
                       plane_index > len(centerline) - edge_threshold)

        if is_near_edge:
            # Near beginning/end - check for actual open contours
            for entity in paths_2d.entities:
                vertices = paths_2d.vertices[entity.points]

                # Check if path is closed (first == last point)
                if len(vertices) >= 2:
                    dist = np.linalg.norm(vertices[0] - vertices[-1])
                    is_path_closed = dist < 1e-6

                    if not is_path_closed:
                        all_closed = False
                        # Count open edges
                        total_boundary_edges += 1
        else:
            # Middle of airway - accept all planes
            # Multiple contours are OK (bifurcations, donuts)
            all_closed = True
            total_boundary_edges = 0

        return all_closed, total_boundary_edges, n_paths

    except Exception as e:
        print(f"Error checking plane {plane_index}: {e}")
        return None, 0, 0


def analyze_all_planes(stl_path, vtk_path, indices_file=None):
    """Analyze all planes for closure status"""
    from slicer.io_utils import read_vtk_centerline

    # Load centerline
    centerline = read_vtk_centerline(vtk_path)

    # Get indices to check
    if indices_file:
        with open(indices_file, 'r') as f:
            indices = [int(line.strip()) for line in f if line.strip()]
    else:
        indices = list(range(len(centerline)))

    print(f"Analyzing {len(indices)} planes for closure status...")
    print(f"Airway mesh: {Path(stl_path).name}")
    print()

    results = []
    open_planes = []
    closed_planes = []

    for idx in indices:
        is_closed, n_boundary, n_paths = check_plane_closure(stl_path, centerline, idx)

        if is_closed is None:
            status = "NO_INTERSECTION"
        elif is_closed:
            status = "CLOSED"
            closed_planes.append(idx)
        else:
            status = "OPEN"
            open_planes.append(idx)

        results.append({
            'plane_index': idx,
            'closure_status': status,
            'is_closed': is_closed if is_closed is not None else False,
            'n_boundary_edges': n_boundary,
            'n_contour_paths': n_paths
        })

        if idx % 20 == 0:
            print(f"  Checked {idx}/{len(indices)} planes...")

    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Total planes analyzed: {len(indices)}")
    print(f"Closed (fully enclosed): {len(closed_planes)}")
    print(f"Open (half-planes): {len(open_planes)}")

    if open_planes:
        print(f"\nOpen plane indices: {open_planes[:10]}{'...' if len(open_planes) > 10 else ''}")

    if closed_planes:
        print(f"Closed plane range: {min(closed_planes)} to {max(closed_planes)}")

    return pd.DataFrame(results)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Detect open vs closed planes")
    parser.add_argument('stl', help='STL mesh file')
    parser.add_argument('vtk', help='VTK centerline file')
    parser.add_argument('--indices', help='File with plane indices to check')
    parser.add_argument('-o', '--output', help='Output CSV file')

    args = parser.parse_args()

    df = analyze_all_planes(args.stl, args.vtk, args.indices)

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\nSaved results to: {args.output}")
    else:
        print(f"\n{df}")
