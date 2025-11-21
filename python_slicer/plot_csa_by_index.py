#!/usr/bin/env python3
"""
Plot CSA vs Centerline Index (Not Arc Length!)

When centerlines are deformed by FFD registration, the same centerline INDEX
represents the same anatomical location across breathing phases, even though
the 3D coordinates and arc lengths change.

This avoids the "spiky" appearance caused by plotting vs arc length.

Usage:
    python plot_csa_by_index.py
    python plot_csa_by_index.py --partition RightNose
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
from typing import List, Dict


def load_measurement_files(folder_path: Path) -> Dict[str, pd.DataFrame]:
    """
    Load all *-Data.csv files from a result folder

    Args:
        folder_path: Path to result folder

    Returns:
        Dictionary mapping file_id to DataFrame
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


def plot_csa_by_index_cycle(data_dict: Dict[str, pd.DataFrame],
                            partition_name: str,
                            output_dir: str = ".") -> str:
    """
    Plot CSA vs centerline INDEX for all breathing cycle phases

    Args:
        data_dict: Dictionary of file_id -> DataFrame
        partition_name: Name of partition
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

    # Sort file IDs
    file_ids = sorted(data_dict.keys())

    print(f"\nPlotting {n_phases} breathing phases (by centerline index)...")

    # Plot each breathing phase
    for i, file_id in enumerate(file_ids):
        df = data_dict[file_id]

        if 'plane_index' not in df.columns or 'area_mm2' not in df.columns:
            print(f"Warning: Missing required columns in {file_id}")
            continue

        # Sort by plane_index (which corresponds to centerline index)
        df_sorted = df.sort_values('plane_index')

        plane_index = df_sorted['plane_index'].values
        area = df_sorted['area_mm2'].values

        # Extract phase number
        phase_num = int(file_id.split('_')[-1])

        # Plot line
        ax.plot(plane_index, area,
                color=colors[i],
                linewidth=1.5,
                alpha=0.7,
                label=f'Phase {phase_num}')

    # Formatting
    ax.set_xlabel('Centerline Index (Anatomical Position)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cross-Sectional Area (mm²)', fontsize=12, fontweight='bold')
    ax.set_title(f'CSA vs Centerline Index (Anatomical Registration)\n{partition_name} - Breathing Cycle',
                 fontsize=14, fontweight='bold', pad=20)

    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Legend
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
              frameon=True, fancybox=True, shadow=True,
              fontsize=8, ncol=1)

    plt.tight_layout()

    # Save figure
    output_path = Path(output_dir) / f"{partition_name}_CSA_by_index_breathing_cycle.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved plot: {output_path}")

    output_path_pdf = Path(output_dir) / f"{partition_name}_CSA_by_index_breathing_cycle.pdf"
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved plot (PDF): {output_path_pdf}")

    plt.close()

    return str(output_path)


def plot_csa_heatmap_by_index(data_dict: Dict[str, pd.DataFrame],
                              partition_name: str,
                              output_dir: str = ".") -> str:
    """
    Create heatmap: rows = centerline index, columns = breathing phase

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

    # Find common index range
    all_indices = []
    for df in data_dict.values():
        if 'plane_index' in df.columns:
            all_indices.extend(df['plane_index'].values)

    if not all_indices:
        print("Warning: No plane index data found")
        return None

    min_idx = int(min(all_indices))
    max_idx = int(max(all_indices))
    index_range = np.arange(min_idx, max_idx + 1)

    # Build 2D array: rows = centerline index, cols = breathing phase
    csa_matrix = np.full((len(index_range), len(file_ids)), np.nan)

    for phase_idx, file_id in enumerate(file_ids):
        df = data_dict[file_id]
        if 'plane_index' not in df.columns or 'area_mm2' not in df.columns:
            continue

        for _, row in df.iterrows():
            plane_idx = int(row['plane_index'])
            row_idx = plane_idx - min_idx
            if 0 <= row_idx < len(index_range):
                csa_matrix[row_idx, phase_idx] = row['area_mm2']

    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 10))

    phase_labels = [f"{int(fid.split('_')[-1])}" for fid in file_ids]

    im = ax.imshow(csa_matrix, aspect='auto', cmap='hot', origin='lower',
                   extent=[0, len(file_ids), min_idx, max_idx],
                   interpolation='nearest')

    ax.set_xlabel('Breathing Phase', fontsize=12, fontweight='bold')
    ax.set_ylabel('Centerline Index (Anatomical Position)', fontsize=12, fontweight='bold')
    ax.set_title(f'CSA Heatmap by Centerline Index\n{partition_name} - Breathing Cycle',
                 fontsize=14, fontweight='bold', pad=20)

    # Set x-ticks
    tick_positions = np.arange(len(file_ids))
    ax.set_xticks(tick_positions + 0.5)
    ax.set_xticklabels(phase_labels, rotation=45, ha='right', fontsize=8)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Cross-Sectional Area (mm²)', fontsize=11, fontweight='bold')

    plt.tight_layout()

    # Save
    output_path = Path(output_dir) / f"{partition_name}_CSA_heatmap_by_index.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved heatmap: {output_path}")

    output_path_pdf = Path(output_dir) / f"{partition_name}_CSA_heatmap_by_index.pdf"
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved heatmap (PDF): {output_path_pdf}")

    plt.close()

    return str(output_path)


def plot_csa_comparison_index_vs_arclength(data_dict: Dict[str, pd.DataFrame],
                                           partition_name: str,
                                           output_dir: str = ".") -> str:
    """
    Side-by-side comparison: CSA by index vs CSA by arc length
    Show why index-based plotting is better for deformed centerlines

    Args:
        data_dict: Dictionary of file_id -> DataFrame
        partition_name: Name of partition
        output_dir: Directory to save plot

    Returns:
        Path to saved plot file
    """
    if not data_dict:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    file_ids = sorted(data_dict.keys())
    n_phases = len(file_ids)
    colors = cm.viridis(np.linspace(0, 1, n_phases))

    # Left plot: By arc length (spiky)
    ax = axes[0]
    for i, file_id in enumerate(file_ids):
        df = data_dict[file_id]
        if 'arc_length_mm' not in df.columns or 'area_mm2' not in df.columns:
            continue

        df_sorted = df.sort_values('arc_length_mm')
        phase_num = int(file_id.split('_')[-1])

        ax.plot(df_sorted['arc_length_mm'], df_sorted['area_mm2'],
                color=colors[i], linewidth=1.5, alpha=0.6)

    ax.set_xlabel('Arc Length (mm)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Cross-Sectional Area (mm²)', fontsize=11, fontweight='bold')
    ax.set_title('CSA vs Arc Length\n(Spiky - Different Anatomical Locations)',
                 fontsize=12, fontweight='bold', color='red')
    ax.grid(True, alpha=0.3)

    # Right plot: By centerline index (smooth)
    ax = axes[1]
    for i, file_id in enumerate(file_ids):
        df = data_dict[file_id]
        if 'plane_index' not in df.columns or 'area_mm2' not in df.columns:
            continue

        df_sorted = df.sort_values('plane_index')
        phase_num = int(file_id.split('_')[-1])

        ax.plot(df_sorted['plane_index'], df_sorted['area_mm2'],
                color=colors[i], linewidth=1.5, alpha=0.6,
                label=f'Phase {phase_num}')

    ax.set_xlabel('Centerline Index', fontsize=11, fontweight='bold')
    ax.set_ylabel('Cross-Sectional Area (mm²)', fontsize=11, fontweight='bold')
    ax.set_title('CSA vs Centerline Index\n(Smooth - Same Anatomical Locations)',
                 fontsize=12, fontweight='bold', color='green')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7, ncol=1)

    plt.suptitle(f'{partition_name} - Comparison: Arc Length vs Index-Based Plotting',
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()

    output_path = Path(output_dir) / f"{partition_name}_comparison_index_vs_arclength.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved comparison plot: {output_path}")

    plt.close()

    return str(output_path)


def plot_csa_by_index_cycle_flipped(data_dict: Dict[str, pd.DataFrame],
                                    partition_name: str,
                                    output_dir: str = ".") -> str:
    """
    Plot CSA vs flipped centerline INDEX (Nose → Trachea direction)
    Clean version without plane index labels

    Args:
        data_dict: Dictionary of file_id -> DataFrame
        partition_name: Name of partition
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

    # Sort file IDs
    file_ids = sorted(data_dict.keys())

    print(f"\nPlotting {n_phases} breathing phases (flipped: nose → trachea)...")

    # Plot each breathing phase
    for i, file_id in enumerate(file_ids):
        df = data_dict[file_id]

        if 'plane_index' not in df.columns or 'area_mm2' not in df.columns:
            print(f"Warning: Missing required columns in {file_id}")
            continue

        # Sort by plane_index
        df_sorted = df.sort_values('plane_index')

        # FLIP the order: reverse both x and y arrays
        plane_index = df_sorted['plane_index'].values[::-1]
        area = df_sorted['area_mm2'].values[::-1]

        # Extract phase number
        phase_num = int(file_id.split('_')[-1])

        # Plot line
        ax.plot(plane_index, area,
                color=colors[i],
                linewidth=1.5,
                alpha=0.7,
                label=f'Phase {phase_num}')

    # Formatting
    ax.set_xlabel('Anatomical Position (Nose → Trachea)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cross-Sectional Area (mm²)', fontsize=12, fontweight='bold')
    ax.set_title(f'CSA Along Airway (Nose → Trachea Direction)\n{partition_name} - Breathing Cycle',
                 fontsize=14, fontweight='bold', pad=20)

    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Invert x-axis to show nose → trachea
    ax.invert_xaxis()

    # Legend
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
              frameon=True, fancybox=True, shadow=True,
              fontsize=8, ncol=1)

    plt.tight_layout()

    # Save figure
    output_path = Path(output_dir) / f"{partition_name}_CSA_flipped.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved flipped plot: {output_path}")

    output_path_pdf = Path(output_dir) / f"{partition_name}_CSA_flipped.pdf"
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved flipped plot (PDF): {output_path_pdf}")

    plt.close()

    return str(output_path)


def plot_csa_by_index_cycle_flipped_labeled(data_dict: Dict[str, pd.DataFrame],
                                            partition_name: str,
                                            output_dir: str = ".") -> str:
    """
    Plot CSA vs flipped centerline INDEX (Nose → Trachea direction)
    Labeled version with plane indices for traceability

    Args:
        data_dict: Dictionary of file_id -> DataFrame
        partition_name: Name of partition
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

    # Sort file IDs
    file_ids = sorted(data_dict.keys())

    print(f"\nPlotting {n_phases} breathing phases (flipped + labeled)...")

    # Use first phase to get plane indices for labeling
    first_df = data_dict[file_ids[0]]
    if 'plane_index' in first_df.columns:
        df_sorted = first_df.sort_values('plane_index')
        all_plane_indices = df_sorted['plane_index'].values[::-1]  # Flipped

    # Plot each breathing phase
    for i, file_id in enumerate(file_ids):
        df = data_dict[file_id]

        if 'plane_index' not in df.columns or 'area_mm2' not in df.columns:
            print(f"Warning: Missing required columns in {file_id}")
            continue

        # Sort by plane_index
        df_sorted = df.sort_values('plane_index')

        # FLIP the order: reverse both x and y arrays
        plane_index = df_sorted['plane_index'].values[::-1]
        area = df_sorted['area_mm2'].values[::-1]

        # Extract phase number
        phase_num = int(file_id.split('_')[-1])

        # Plot line
        ax.plot(plane_index, area,
                color=colors[i],
                linewidth=1.5,
                alpha=0.7,
                label=f'Phase {phase_num}')

    # Add plane index labels at reasonable intervals
    if len(all_plane_indices) > 0:
        # Label every 10th plane
        label_interval = max(1, len(all_plane_indices) // 15)
        for idx in range(0, len(all_plane_indices), label_interval):
            plane_idx = all_plane_indices[idx]
            # Add subtle vertical line
            ax.axvline(x=plane_idx, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            # Add text label
            ax.text(plane_idx, ax.get_ylim()[1] * 0.98, f'{plane_idx}',
                   ha='center', va='top', fontsize=7, color='gray', alpha=0.7)

    # Formatting
    ax.set_xlabel('Anatomical Position (Nose → Trachea) [Plane Index]', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cross-Sectional Area (mm²)', fontsize=12, fontweight='bold')
    ax.set_title(f'CSA Along Airway (Nose → Trachea Direction) - With Plane Indices\n{partition_name} - Breathing Cycle',
                 fontsize=14, fontweight='bold', pad=20)

    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Invert x-axis to show nose → trachea
    ax.invert_xaxis()

    # Legend
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
              frameon=True, fancybox=True, shadow=True,
              fontsize=8, ncol=1)

    plt.tight_layout()

    # Save figure
    output_path = Path(output_dir) / f"{partition_name}_CSA_flipped_labeled.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved flipped labeled plot: {output_path}")

    output_path_pdf = Path(output_dir) / f"{partition_name}_CSA_flipped_labeled.pdf"
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved flipped labeled plot (PDF): {output_path_pdf}")

    plt.close()

    return str(output_path)


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Plot CSA by Centerline Index (Anatomically Registered)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot_csa_by_index.py
  python plot_csa_by_index.py --partition LeftNoseDecending
  python plot_csa_by_index.py -o plots/
        """
    )

    parser.add_argument(
        '--partition',
        type=str,
        default='LeftNoseDecending',
        help='Partition to plot (default: LeftNoseDecending)'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default='.',
        help='Output directory (default: current directory)'
    )

    parser.add_argument(
        '--comparison',
        action='store_true',
        help='Generate comparison plot (index vs arc length)'
    )

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*80)
    print("CSA vs Centerline Index - Anatomically Registered Plotting")
    print("="*80)

    # Find result folder
    folder_name = f"{args.partition}SlicedSTLs"
    result_folder = Path(folder_name)

    if not result_folder.exists():
        print(f"Error: Folder not found: {folder_name}")
        sys.exit(1)

    print(f"\nProcessing: {args.partition}")

    # Load data
    data_dict = load_measurement_files(result_folder)

    if not data_dict:
        print(f"Error: No data found in {result_folder}")
        sys.exit(1)

    print(f"Loaded {len(data_dict)} breathing phase files")

    # Generate original plots (trachea → nose)
    plot_csa_by_index_cycle(data_dict, args.partition, str(output_dir))
    plot_csa_heatmap_by_index(data_dict, args.partition, str(output_dir))

    # Generate flipped plots (nose → trachea)
    plot_csa_by_index_cycle_flipped(data_dict, args.partition, str(output_dir))
    plot_csa_by_index_cycle_flipped_labeled(data_dict, args.partition, str(output_dir))

    if args.comparison:
        plot_csa_comparison_index_vs_arclength(data_dict, args.partition, str(output_dir))

    print(f"\n{'='*80}")
    print("Index-based plots generated successfully!")
    print(f"  - Original plots (trachea → nose)")
    print(f"  - Flipped plots (nose → trachea, clean + labeled)")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
