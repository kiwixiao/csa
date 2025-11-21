#!/usr/bin/env python3
"""
Check for plane overlaps in individual plane STL files
"""

import sys
import glob
import numpy as np
import trimesh
from pathlib import Path
import matplotlib.pyplot as plt


def check_planes_for_overlaps(file_pattern: str):
    """
    Load all individual plane STLs and check for Z-coordinate overlaps

    Args:
        file_pattern: Glob pattern like "LeftNoseDecendingSlicedSTLs/out_002000-Planes-*.stl"
    """
    # Get all plane files
    plane_files = sorted(glob.glob(file_pattern))

    if not plane_files:
        print(f"Error: No files found matching pattern: {file_pattern}")
        return

    print("="*80)
    print(f"Checking {len(plane_files)} plane files for overlaps")
    print("="*80)

    # Load each plane and extract geometry info
    plane_info = []

    for plane_file in plane_files:
        try:
            mesh = trimesh.load_mesh(plane_file)

            # Extract plane index from filename
            filename = Path(plane_file).stem
            # out_002000-Planes-145 -> 145
            plane_idx = int(filename.split('-')[-1])

            # Get Z bounds
            z_min = mesh.vertices[:, 2].min()
            z_max = mesh.vertices[:, 2].max()
            z_centroid = mesh.centroid[2]

            # Get plane normal by fitting plane to vertices
            vertices = mesh.vertices
            centered = vertices - vertices.mean(axis=0)
            U, S, Vt = np.linalg.svd(centered, full_matrices=False)
            plane_normal = Vt[-1]

            plane_info.append({
                'file': plane_file,
                'index': plane_idx,
                'z_min': z_min,
                'z_max': z_max,
                'z_centroid': z_centroid,
                'z_thickness': z_max - z_min,
                'centroid': mesh.centroid,
                'normal': plane_normal,
                'n_vertices': len(mesh.vertices)
            })

        except Exception as e:
            print(f"Warning: Failed to load {plane_file}: {e}")

    # Sort by plane index
    plane_info.sort(key=lambda p: p['index'])

    print(f"\nLoaded {len(plane_info)} planes successfully\n")

    # Print table
    print(f"{'Index':<8} {'Z_centroid':<12} {'Z_min':<12} {'Z_max':<12} {'Thickness':<10} {'Normal_Z':<10}")
    print("-"*80)

    for p in plane_info:
        print(f"{p['index']:<8} {p['z_centroid']:<12.2f} {p['z_min']:<12.2f} {p['z_max']:<12.2f} "
              f"{p['z_thickness']:<10.2f} {p['normal'][2]:<10.3f}")

    # Check for overlaps
    print(f"\n{'='*80}")
    print("OVERLAP DETECTION")
    print(f"{'='*80}\n")

    overlaps = []

    for i in range(len(plane_info) - 1):
        curr = plane_info[i]
        next_plane = plane_info[i + 1]

        # Check if Z ranges overlap
        if curr['z_max'] > next_plane['z_min']:
            overlap_amount = curr['z_max'] - next_plane['z_min']
            overlaps.append({
                'plane1': curr['index'],
                'plane2': next_plane['index'],
                'overlap_mm': overlap_amount,
                'plane1_range': (curr['z_min'], curr['z_max']),
                'plane2_range': (next_plane['z_min'], next_plane['z_max'])
            })

    if len(overlaps) == 0:
        print("✓ NO OVERLAPS DETECTED - All planes are properly separated!")
    else:
        print(f"✗ FOUND {len(overlaps)} OVERLAPPING PLANE PAIRS:\n")

        for overlap in overlaps:
            print(f"  Plane {overlap['plane1']:3d} [Z: {overlap['plane1_range'][0]:.2f} - {overlap['plane1_range'][1]:.2f}]")
            print(f"    ↓ OVERLAPS by {overlap['overlap_mm']:.2f} mm")
            print(f"  Plane {overlap['plane2']:3d} [Z: {overlap['plane2_range'][0]:.2f} - {overlap['plane2_range'][1]:.2f}]")
            print()

    # Statistics
    print(f"{'='*80}")
    print("STATISTICS")
    print(f"{'='*80}")

    z_centroids = [p['z_centroid'] for p in plane_info]
    spacings = np.diff(z_centroids)
    thicknesses = [p['z_thickness'] for p in plane_info]

    print(f"\nPlane spacing (centroid to centroid):")
    print(f"  Mean: {np.mean(spacings):.2f} mm")
    print(f"  Std:  {np.std(spacings):.2f} mm")
    print(f"  Min:  {np.min(spacings):.2f} mm")
    print(f"  Max:  {np.max(spacings):.2f} mm")

    print(f"\nPlane thickness (Z extent):")
    print(f"  Mean: {np.mean(thicknesses):.2f} mm")
    print(f"  Std:  {np.std(thicknesses):.2f} mm")
    print(f"  Min:  {np.min(thicknesses):.2f} mm")
    print(f"  Max:  {np.max(thicknesses):.2f} mm")

    print(f"\nOverlap summary:")
    print(f"  Total plane pairs: {len(plane_info) - 1}")
    print(f"  Overlapping pairs: {len(overlaps)}")
    print(f"  Overlap percentage: {100 * len(overlaps) / (len(plane_info) - 1):.1f}%")

    if overlaps:
        overlap_amounts = [o['overlap_mm'] for o in overlaps]
        print(f"  Mean overlap amount: {np.mean(overlap_amounts):.2f} mm")
        print(f"  Max overlap amount: {np.max(overlap_amounts):.2f} mm")

    # Visualize
    visualize_overlaps(plane_info, overlaps, file_pattern)

    return plane_info, overlaps


def visualize_overlaps(plane_info, overlaps, file_pattern):
    """Create visualization of plane positions and overlaps"""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    indices = [p['index'] for p in plane_info]
    z_cents = [p['z_centroid'] for p in plane_info]
    z_mins = [p['z_min'] for p in plane_info]
    z_maxs = [p['z_max'] for p in plane_info]
    thicknesses = [p['z_thickness'] for p in plane_info]

    # Plot 1: Z position vs index
    ax = axes[0, 0]
    ax.plot(indices, z_cents, 'b-o', label='Centroid', markersize=3)
    ax.fill_between(indices, z_mins, z_maxs, alpha=0.3, label='Z extent')

    # Highlight overlaps
    for overlap in overlaps:
        idx1 = overlap['plane1']
        idx2 = overlap['plane2']
        ax.axvspan(idx1, idx2, color='red', alpha=0.2)

    ax.set_xlabel('Plane Index', fontweight='bold')
    ax.set_ylabel('Z Coordinate (mm)', fontweight='bold')
    ax.set_title('Plane Positions (Red = Overlap)', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Plane spacing
    ax = axes[0, 1]
    spacings = np.diff(z_cents)
    spacing_indices = indices[:-1]
    ax.plot(spacing_indices, spacings, 'g-o', markersize=3)
    ax.axhline(y=0, color='r', linestyle='--', label='Zero spacing')
    ax.set_xlabel('Plane Index', fontweight='bold')
    ax.set_ylabel('Spacing to Next Plane (mm)', fontweight='bold')
    ax.set_title('Inter-Plane Spacing', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Plane thickness
    ax = axes[1, 0]
    ax.plot(indices, thicknesses, 'purple', marker='o', markersize=3)
    ax.set_xlabel('Plane Index', fontweight='bold')
    ax.set_ylabel('Plane Thickness (mm)', fontweight='bold')
    ax.set_title('Plane Z-Extent (Should be thin)', fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Plot 4: Overlap amounts
    ax = axes[1, 1]
    if overlaps:
        overlap_indices = [o['plane1'] for o in overlaps]
        overlap_amounts = [o['overlap_mm'] for o in overlaps]
        ax.bar(overlap_indices, overlap_amounts, color='red', alpha=0.6, width=0.8)
        ax.set_xlabel('Plane Index', fontweight='bold')
        ax.set_ylabel('Overlap Amount (mm)', fontweight='bold')
        ax.set_title(f'Overlap Amounts ({len(overlaps)} pairs)', fontweight='bold', color='red')
        ax.grid(True, alpha=0.3, axis='y')
    else:
        ax.text(0.5, 0.5, 'No Overlaps Detected!',
                ha='center', va='center', fontsize=16, color='green',
                fontweight='bold', transform=ax.transAxes)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')

    plt.tight_layout()

    # Save
    output_name = Path(file_pattern).parent / f"{Path(file_pattern).stem.split('-')[0]}_overlap_analysis.png"
    plt.savefig(output_name, dpi=200, bbox_inches='tight')
    print(f"\nVisualization saved: {output_name}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_plane_overlaps.py <file_pattern>")
        print("\nExample:")
        print("  python check_plane_overlaps.py 'LeftNoseDecendingSlicedSTLs/out_002000-Planes-*.stl'")
        sys.exit(1)

    file_pattern = sys.argv[1]
    check_planes_for_overlaps(file_pattern)
