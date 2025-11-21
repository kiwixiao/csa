#!/usr/bin/env python3
"""
Plot CSA vs Centerline Distance Over Breathing Cycle

This script automatically detects result folders and plots cross-sectional area
vs arc length along the centerline for all breathing cycle phases.

Usage:
    python plot_csa_cycle.py                          # Auto-detect all partitions
    python plot_csa_cycle.py --partition LeftNoseDecending
    python plot_csa_cycle.py --subject OSAMRI037 --partition RightNose
"""

import os
import sys
import glob
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import List, Dict, Tuple


def find_result_folders(base_dir: str = ".") -> List[Path]:
    """
    Auto-detect result folders (e.g., LeftNoseDecendingSlicedSTLs, RightNoseSlicedSTLs)

    Args:
        base_dir: Base directory to search (default: current directory)

    Returns:
        List of Path objects for detected result folders
    """
    base_path = Path(base_dir)

    # Look for folders ending with "SlicedSTLs"
    result_folders = []
    for folder in base_path.glob("*SlicedSTLs"):
        if folder.is_dir():
            result_folders.append(folder)

    return sorted(result_folders)


def load_measurement_files(folder_path: Path) -> Dict[str, pd.DataFrame]:
    """
    Load all *-Data.csv files from a result folder

    Args:
        folder_path: Path to result folder

    Returns:
        Dictionary mapping file_id (e.g., 'out_000000') to DataFrame
    """
    csv_files = sorted(glob.glob(str(folder_path / "*-Data.csv")))

    if len(csv_files) == 0:
        print(f"Warning: No CSV files found in {folder_path}")
        return {}

    data_dict = {}
    for csv_path in csv_files:
        file_id = Path(csv_path).stem.replace("-Data", "")
        try:
            df = pd.read_csv(csv_path)
            # Filter to only valid planes
            if 'is_valid' in df.columns:
                df = df[df['is_valid'] == True].copy()
            data_dict[file_id] = df
        except Exception as e:
            print(f"Warning: Failed to load {csv_path}: {e}")

    return data_dict


def plot_csa_breathing_cycle(data_dict: Dict[str, pd.DataFrame],
                             partition_name: str,
                             output_dir: str = ".") -> str:
    """
    Plot CSA vs arc length for all breathing cycle phases

    Args:
        data_dict: Dictionary of file_id -> DataFrame
        partition_name: Name of partition (e.g., 'LeftNoseDecending')
        output_dir: Directory to save plot

    Returns:
        Path to saved plot file
    """
    if not data_dict:
        print("Error: No data to plot")
        return None

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))

    # Get colormap for breathing phases
    n_phases = len(data_dict)
    colors = cm.viridis(np.linspace(0, 1, n_phases))

    # Sort file IDs to ensure chronological order
    file_ids = sorted(data_dict.keys())

    print(f"\nPlotting {n_phases} breathing phases...")

    # Plot each breathing phase
    for i, file_id in enumerate(file_ids):
        df = data_dict[file_id]

        if 'arc_length_mm' not in df.columns or 'area_mm2' not in df.columns:
            print(f"Warning: Missing required columns in {file_id}")
            continue

        # Sort by arc length for proper line plotting
        df_sorted = df.sort_values('arc_length_mm')

        arc_length = df_sorted['arc_length_mm'].values
        area = df_sorted['area_mm2'].values

        # Extract phase number from file_id (e.g., out_000000 -> 0)
        phase_num = int(file_id.split('_')[-1])

        # Plot line with marker
        ax.plot(arc_length, area,
                color=colors[i],
                linewidth=1.5,
                alpha=0.7,
                label=f'Phase {phase_num}')

    # Formatting
    ax.set_xlabel('Arc Length Along Centerline (mm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cross-Sectional Area (mm²)', fontsize=12, fontweight='bold')
    ax.set_title(f'Cross-Sectional Area vs Centerline Distance\n{partition_name} - Breathing Cycle',
                 fontsize=14, fontweight='bold', pad=20)

    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Legend - put outside plot area
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
              frameon=True, fancybox=True, shadow=True,
              fontsize=8, ncol=1)

    # Tight layout to prevent legend cutoff
    plt.tight_layout()

    # Save figure
    output_path = Path(output_dir) / f"{partition_name}_CSA_breathing_cycle.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved plot: {output_path}")

    # Also save as PDF for high-quality vector graphics
    output_path_pdf = Path(output_dir) / f"{partition_name}_CSA_breathing_cycle.pdf"
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved plot (PDF): {output_path_pdf}")

    plt.close()

    return str(output_path)


def plot_csa_heatmap(data_dict: Dict[str, pd.DataFrame],
                     partition_name: str,
                     output_dir: str = ".") -> str:
    """
    Create heatmap showing CSA variation over breathing cycle

    Args:
        data_dict: Dictionary of file_id -> DataFrame
        partition_name: Name of partition
        output_dir: Directory to save plot

    Returns:
        Path to saved plot file
    """
    if not data_dict:
        return None

    # Sort file IDs
    file_ids = sorted(data_dict.keys())

    # Build 2D array: rows = arc length positions, cols = breathing phases
    # First, find common arc length grid
    all_arc_lengths = []
    for df in data_dict.values():
        if 'arc_length_mm' in df.columns:
            all_arc_lengths.extend(df['arc_length_mm'].values)

    if not all_arc_lengths:
        print("Warning: No arc length data found for heatmap")
        return None

    # Create uniform grid
    arc_min = min(all_arc_lengths)
    arc_max = max(all_arc_lengths)
    n_grid_points = 200
    arc_grid = np.linspace(arc_min, arc_max, n_grid_points)

    # Interpolate each phase onto common grid
    csa_matrix = np.zeros((n_grid_points, len(file_ids)))

    for i, file_id in enumerate(file_ids):
        df = data_dict[file_id]
        if 'arc_length_mm' not in df.columns or 'area_mm2' not in df.columns:
            continue

        df_sorted = df.sort_values('arc_length_mm')
        arc_length = df_sorted['arc_length_mm'].values
        area = df_sorted['area_mm2'].values

        # Interpolate onto grid
        csa_interp = np.interp(arc_grid, arc_length, area, left=np.nan, right=np.nan)
        csa_matrix[:, i] = csa_interp

    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 8))

    phase_labels = [f"{int(fid.split('_')[-1])}" for fid in file_ids]

    im = ax.imshow(csa_matrix, aspect='auto', cmap='hot', origin='lower',
                   extent=[0, len(file_ids), arc_min, arc_max],
                   interpolation='bilinear')

    ax.set_xlabel('Breathing Phase', fontsize=12, fontweight='bold')
    ax.set_ylabel('Arc Length Along Centerline (mm)', fontsize=12, fontweight='bold')
    ax.set_title(f'CSA Heatmap - {partition_name}\nBreathing Cycle Variation',
                 fontsize=14, fontweight='bold', pad=20)

    # Set x-ticks to show phase numbers
    tick_positions = np.arange(len(file_ids))
    ax.set_xticks(tick_positions + 0.5)
    ax.set_xticklabels(phase_labels, rotation=45, ha='right', fontsize=8)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Cross-Sectional Area (mm²)', fontsize=11, fontweight='bold')

    plt.tight_layout()

    # Save
    output_path = Path(output_dir) / f"{partition_name}_CSA_heatmap.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved heatmap: {output_path}")

    output_path_pdf = Path(output_dir) / f"{partition_name}_CSA_heatmap.pdf"
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved heatmap (PDF): {output_path_pdf}")

    plt.close()

    return str(output_path)


def generate_statistics_summary(data_dict: Dict[str, pd.DataFrame],
                                partition_name: str,
                                output_dir: str = ".") -> str:
    """
    Generate statistical summary across breathing cycle

    Args:
        data_dict: Dictionary of file_id -> DataFrame
        partition_name: Name of partition
        output_dir: Directory to save summary

    Returns:
        Path to saved summary CSV
    """
    summary_data = []

    file_ids = sorted(data_dict.keys())

    for file_id in file_ids:
        df = data_dict[file_id]

        if 'area_mm2' not in df.columns:
            continue

        phase_num = int(file_id.split('_')[-1])

        stats = {
            'phase': phase_num,
            'file_id': file_id,
            'n_planes': len(df),
            'area_min_mm2': df['area_mm2'].min(),
            'area_max_mm2': df['area_mm2'].max(),
            'area_mean_mm2': df['area_mm2'].mean(),
            'area_std_mm2': df['area_mm2'].std(),
            'area_median_mm2': df['area_mm2'].median(),
        }

        # Add hydraulic diameter if available
        if 'hydraulic_diameter_mm' in df.columns:
            stats['hydraulic_diameter_mean_mm'] = df['hydraulic_diameter_mm'].mean()
            stats['hydraulic_diameter_min_mm'] = df['hydraulic_diameter_mm'].min()
            stats['hydraulic_diameter_max_mm'] = df['hydraulic_diameter_mm'].max()

        summary_data.append(stats)

    summary_df = pd.DataFrame(summary_data)

    # Save
    output_path = Path(output_dir) / f"{partition_name}_breathing_cycle_statistics.csv"
    summary_df.to_csv(output_path, index=False)
    print(f"\nSaved statistics: {output_path}")

    # Print summary
    print(f"\n{'='*80}")
    print(f"Breathing Cycle Statistics - {partition_name}")
    print(f"{'='*80}")
    print(f"Total phases: {len(summary_df)}")
    print(f"Mean CSA across cycle: {summary_df['area_mean_mm2'].mean():.2f} ± {summary_df['area_mean_mm2'].std():.2f} mm²")
    print(f"CSA range: {summary_df['area_mean_mm2'].min():.2f} - {summary_df['area_mean_mm2'].max():.2f} mm²")
    print(f"Max CSA variation: {summary_df['area_max_mm2'].max():.2f} mm²")
    print(f"Min CSA variation: {summary_df['area_min_mm2'].min():.2f} mm²")

    return str(output_path)


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Plot CSA vs Centerline Distance Over Breathing Cycle",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot_csa_cycle.py                                    # Auto-detect all partitions
  python plot_csa_cycle.py --partition LeftNoseDecending     # Specific partition
  python plot_csa_cycle.py -o plots/                         # Custom output directory
        """
    )

    parser.add_argument(
        '--subject',
        type=str,
        default=None,
        help='Subject name (optional, for labeling only)'
    )

    parser.add_argument(
        '--partition',
        type=str,
        default=None,
        help='Specific partition to plot (e.g., LeftNoseDecending). If not specified, auto-detect all.'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default='.',
        help='Output directory for plots (default: current directory)'
    )

    parser.add_argument(
        '--no-heatmap',
        action='store_true',
        help='Skip heatmap generation (only create line plots)'
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*80)
    print("CSA vs Centerline Distance - Breathing Cycle Analysis")
    print("="*80)

    # Find result folders
    if args.partition:
        # Look for specific partition folder
        folder_name = f"{args.partition}SlicedSTLs"
        result_folders = [Path(folder_name)]
        if not result_folders[0].exists():
            print(f"Error: Folder not found: {folder_name}")
            sys.exit(1)
    else:
        # Auto-detect all result folders
        result_folders = find_result_folders()
        if not result_folders:
            print("Error: No result folders found (*SlicedSTLs)")
            print("Please run the slicer first or specify --partition")
            sys.exit(1)

    print(f"\nFound {len(result_folders)} result folder(s):")
    for folder in result_folders:
        print(f"  - {folder}")

    # Process each folder
    for folder in result_folders:
        partition_name = folder.name.replace("SlicedSTLs", "")

        print(f"\n{'='*80}")
        print(f"Processing: {partition_name}")
        print(f"{'='*80}")

        # Load measurement files
        data_dict = load_measurement_files(folder)

        if not data_dict:
            print(f"Warning: No data found in {folder}, skipping...")
            continue

        print(f"Loaded {len(data_dict)} breathing phase files")

        # Generate plots
        plot_csa_breathing_cycle(data_dict, partition_name, str(output_dir))

        if not args.no_heatmap:
            plot_csa_heatmap(data_dict, partition_name, str(output_dir))

        # Generate statistics
        generate_statistics_summary(data_dict, partition_name, str(output_dir))

    print(f"\n{'='*80}")
    print("All plots generated successfully!")
    print(f"Output directory: {output_dir.absolute()}")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
