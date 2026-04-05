#!/usr/bin/env python3
"""
Python Slicer - Consistent Planes Across Breathing Cycle

Uses pre-determined valid plane indices for all time points to ensure:
1. No overlapping planes
2. Consistent anatomical locations across breathing phases
3. Smooth, forward-progressing planes

Usage:
    python main_consistent.py <subject_name> <partition> --indices <indices_file>

Example:
    python main_consistent.py OSAMRI037 LeftNoseDecending --indices LeftNoseDecending_valid_plane_indices.txt
"""

import os
import sys
import glob
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.ndimage import uniform_filter1d

from slicer import AirwaySlicer
from slicer.io_utils import save_results_pickle, read_vtk_centerline
from slicer.validation import detect_split_regions


def smooth_centerline(centerline: np.ndarray, window: int = 5) -> np.ndarray:
    """Smooth centerline points using moving average"""
    if len(centerline) < window:
        return centerline

    padded = np.pad(centerline, ((window//2, window//2), (0, 0)), mode='edge')
    smoothed = uniform_filter1d(padded, size=window, axis=0)
    return smoothed[window//2:window//2+len(centerline)]


def load_valid_indices(indices_file: str) -> list:
    """Load pre-computed valid plane indices"""
    with open(indices_file, 'r') as f:
        indices = [int(line.strip()) for line in f if line.strip()]
    return indices


def process_subject_consistent(subject_name: str,
                               partition: str,
                               valid_indices: list,
                               centerline_smooth_window: int = 5,
                               output_base_dir: str = "."):
    """
    Process all time points using consistent plane indices

    Args:
        subject_name: Subject ID
        partition: Partition name
        valid_indices: List of valid plane indices to use
        centerline_smooth_window: Window for centerline smoothing
        output_base_dir: Output directory
    """
    print("="*80)
    print(f"Python Slicer - Consistent Planes Mode")
    print("="*80)
    print(f"Subject: {subject_name}")
    print(f"Partition: {partition}")
    print(f"Using {len(valid_indices)} pre-determined plane indices")
    print(f"Centerline smoothing window: {centerline_smooth_window}")
    print("="*80)

    # Setup paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    vtk_dir = project_root / partition / "FFD" / "vtk"
    stl_dir = project_root / partition / "FFD" / "stl"
    output_dir = Path(output_base_dir) / f"{partition}SlicedSTLs_Consistent"

    # Validate
    if not vtk_dir.exists() or not stl_dir.exists():
        print(f"Error: Input directories not found")
        sys.exit(1)

    # Clean output
    if output_dir.exists():
        print(f"\nCleaning output directory: {output_dir}")
        import shutil
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get files
    vtk_files = sorted(glob.glob(str(vtk_dir / "*.vtk")))
    stl_files = sorted(glob.glob(str(stl_dir / "*.stl")))

    if len(vtk_files) == 0 or len(stl_files) == 0:
        print(f"Error: No files found")
        sys.exit(1)

    print(f"\nFound {len(vtk_files)} VTK files and {len(stl_files)} STL files")

    # Storage
    all_results = {}
    all_measurements = []

    n_files = min(len(vtk_files), len(stl_files))

    for i in range(n_files):
        vtk_path = vtk_files[i]
        stl_path = stl_files[i]

        vtk_name = Path(vtk_path).name
        stl_name = Path(stl_path).name

        print(f"\n{'='*80}")
        print(f"Processing time point {i+1}/{n_files}")
        print(f"  VTK: {vtk_name}")
        print(f"  STL: {stl_name}")
        print(f"{'='*80}")

        # Load centerline
        centerline = read_vtk_centerline(vtk_path)
        centerline_smooth = smooth_centerline(centerline, window=centerline_smooth_window)

        # Create slicer with ORIGINAL centerline positions
        slicer = AirwaySlicer(stl_path, vtk_path)

        # Keep original centerline positions (don't override!)
        # Only compute normals from smoothed centerline for smooth orientation
        from slicer.geometry import compute_all_plane_normals
        slicer.plane_normals = compute_all_plane_normals(centerline_smooth, smooth=True)

        # Slice only at valid indices
        results = slicer.slice_along_centerline(
            quality_checks=False,  # Already validated
            specific_indices=valid_indices
        )

        print(f"\nSlicing complete:")
        print(f"  Valid planes: {len(results.valid_sections)} / {len(valid_indices)}")

        # Compute diameter profile
        diameter_profile = results.compute_diameter_profile()

        # Export
        base_name = Path(stl_path).stem
        slicer.export_results(results, base_name, str(output_dir))

        # Store
        all_results[base_name] = {
            'slicing_results': results,
            'diameter_profile': diameter_profile
        }

        # Collect measurements
        df = diameter_profile.to_dataframe()
        df['file_name'] = base_name
        df['file_index'] = i
        all_measurements.append(df)

    # Save combined results
    print(f"\n{'='*80}")
    print("Saving combined results...")
    results_file = Path(output_base_dir) / f"{subject_name}_{partition}_results_consistent.pkl"
    save_results_pickle(str(results_file), all_results)

    if all_measurements:
        combined_df = pd.concat(all_measurements, ignore_index=True)
        summary_csv = Path(output_base_dir) / f"{subject_name}_{partition}_all_measurements_consistent.csv"
        combined_df.to_csv(summary_csv, index=False)
        print(f"Saved combined measurements: {summary_csv}")

    print(f"\n{'='*80}")
    print("Processing complete!")
    print(f"Output: {output_dir}")
    print(f"{'='*80}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Python Slicer - Consistent Planes Mode"
    )

    parser.add_argument('subject', type=str, help='Subject ID')
    parser.add_argument('partition', type=str, help='Partition name')
    parser.add_argument('--indices', type=str, required=True,
                       help='File with valid plane indices (one per line)')
    parser.add_argument('--window', type=int, default=5,
                       help='Centerline smoothing window (default: 5)')
    parser.add_argument('-o', '--output', type=str, default='.',
                       help='Output directory (default: current)')

    args = parser.parse_args()

    # Load valid indices
    if not Path(args.indices).exists():
        print(f"Error: Indices file not found: {args.indices}")
        sys.exit(1)

    valid_indices = load_valid_indices(args.indices)

    # Process
    try:
        process_subject_consistent(
            args.subject,
            args.partition,
            valid_indices,
            centerline_smooth_window=args.window,
            output_base_dir=args.output
        )
    except KeyboardInterrupt:
        print("\n\nProcessing interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nError: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
