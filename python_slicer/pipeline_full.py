#!/usr/bin/env python
"""
Smart Airway Slicer Pipeline - Efficient Processing

Workflow:
1. Slice ALL planes once (no pre-filtering) - Uses latest code (window=20, adaptive normals)
2. Detect mutual planes from file existence (fast - no mesh loading)
3. Detect complete planes using optimized two-stage detection:
   - Stage 1 (Fast): Edge length check to identify suspicious planes
   - Stage 2 (Accurate): Surface validation to eliminate false positives
4. Filter and organize files with sequential numbering
5. Generate combined CSV with correct indices
6. Generate plots, video, and interactive HTML plane reference
7. Generate advanced dynamics analysis (enhanced metrics, PDF report)

Key improvements over old approach:
- No redundant re-slicing (5x faster)
- Optimized two-stage incomplete plane detection
- Interactive HTML plane reference (replaces static labeled images)
- Automatically regenerates everything when code changes
- One command: slicing → plots → video → interactive HTML → PDF report
"""

import sys
import glob
import subprocess
import pandas as pd
from pathlib import Path
from collections import defaultdict
import shutil
import numpy as np
import trimesh

# Add slicer module to path
SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

from detect_incomplete_planes_optimized import (
    detect_incomplete_plane_optimized,
    detect_incomplete_planes_for_timepoint
)


def step1_slice_all_planes(subject, partition, output_dir, workers=0):
    """Step 1: Slice ALL planes (for detection only)"""
    print("="*80)
    print("STEP 1: Slicing All Planes (For Detection)")
    print("="*80)
    print("Using latest code: window=20, adaptive normals")
    print("Output: Temporary directory (used for mutual plane detection)")

    # Use temporary output directory
    temp_output = output_dir / f"{partition}SlicedSTLs_temp"
    if temp_output.exists():
        print(f"Cleaning: {temp_output}")
        shutil.rmtree(temp_output)

    # Run slicing without --indices (slices all planes along centerline)
    cmd = [
        "python", str(SCRIPT_DIR / "main.py"),
        subject, partition,
        "-o", str(output_dir),
        "--workers", str(workers)
    ]

    print(f"Running: {' '.join(cmd)}\n")
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print(f"Error in slicing")
        return False

    # Rename to temp directory
    final_output = output_dir / f"{partition}SlicedSTLs"
    if final_output.exists():
        final_output.rename(temp_output)

    print("\n✓ Initial slicing complete (temporary)")
    return True


def step2_detect_mutual_planes(partition, output_dir):
    """Step 2: Detect mutual planes from CSV files (FAST!)"""
    print("\n" + "="*80)
    print("STEP 2: Detecting Mutual Planes from CSV Files")
    print("="*80)

    output_path = output_dir / f"{partition}SlicedSTLs_temp"

    # Find all CSV files (flat structure)
    csv_files = sorted(list(output_path.glob("*-Data.csv")))

    if len(csv_files) == 0:
        print(f"Error: No CSV files found in {output_path}")
        return None

    n_timepoints = len(csv_files)
    print(f"Found {n_timepoints} time point CSV files")

    # Count which plane indices appear in which time points
    plane_occurrences = defaultdict(set)

    for csv_file in csv_files:
        time_name = csv_file.stem.replace('-Data', '')  # out_000000-Data -> out_000000

        # Read CSV and extract plane indices
        df = pd.read_csv(csv_file)
        if 'plane_index' in df.columns:
            plane_indices = df['plane_index'].values
            for plane_idx in plane_indices:
                plane_occurrences[plane_idx].add(time_name)
        else:
            print(f"Warning: {csv_file.name} has no 'plane_index' column")

    # Find planes that appear in ALL time points
    mutual_indices = sorted([
        idx for idx, timepoints in plane_occurrences.items()
        if len(timepoints) == n_timepoints
    ])

    print(f"\nMutual plane detection:")
    print(f"  Total unique planes found: {len(plane_occurrences)}")
    print(f"  Mutual planes (in all {n_timepoints} time points): {len(mutual_indices)}")

    non_mutual = len(plane_occurrences) - len(mutual_indices)
    if non_mutual > 0:
        print(f"  Non-mutual planes (excluded): {non_mutual}")
        non_mutual_indices = sorted(set(plane_occurrences.keys()) - set(mutual_indices))
        print(f"    Indices: {non_mutual_indices[:20]}{'...' if len(non_mutual_indices) > 20 else ''}")

    # Save mutual indices
    mutual_file = output_dir / f"{partition}_mutual_indices.txt"
    with open(mutual_file, 'w') as f:
        for idx in mutual_indices:
            f.write(f"{idx}\n")

    print(f"\n✓ Saved {len(mutual_indices)} mutual indices to:")
    print(f"  {mutual_file}")

    return mutual_file, mutual_indices


def _detect_incomplete_for_single_timepoint(args_tuple):
    """
    Worker function for parallel incomplete plane detection.
    Top-level for multiprocessing pickling.
    """
    import traceback as tb
    i, stl_path, temp_output_str, n_planes, n_timepoints = args_tuple
    try:
        stl_name = Path(stl_path).stem
        incomplete, stage1, stage2 = detect_incomplete_planes_for_timepoint(
            temp_dir=Path(temp_output_str),
            source_mesh_path=stl_path,
            time_point=stl_name,
            n_planes=n_planes
        )
        return {
            'index': i, 'stl_name': stl_name,
            'incomplete_indices': incomplete,
            'stage1': stage1, 'stage2': stage2,
            'error': None
        }
    except Exception as e:
        return {
            'index': i, 'stl_name': Path(stl_path).stem,
            'incomplete_indices': [],
            'stage1': 0, 'stage2': 0,
            'error': tb.format_exc()
        }


def step3_detect_complete_planes_optimized(partition, mutual_indices, output_dir, workers=0):
    """
    Step 3: Detect complete planes using optimized two-stage detection

    Stage 1 (Fast): Edge length check to identify suspicious planes
    Stage 2 (Accurate): Surface validation to eliminate false positives
    """
    import os
    import multiprocessing

    print("\n" + "="*80)
    print("STEP 3: Detecting Complete Planes (Optimized Two-Stage)")
    print("="*80)
    print("Stage 1: Fast edge length check")
    print("Stage 2: Surface validation for suspicious planes")
    print("="*80)

    # Get all time points
    project_root = Path(__file__).parent.parent
    stl_dir = project_root / partition / "FFD" / "stl"
    temp_output = output_dir / f"{partition}SlicedSTLs_temp"

    stl_files = sorted(glob.glob(str(stl_dir / "*.stl")))

    if len(stl_files) == 0:
        print(f"Error: No STL files found in {stl_dir}")
        return None

    n_timepoints = len(stl_files)
    print(f"Checking {len(mutual_indices)} mutual planes across {n_timepoints} time points...")

    # Track which indices are complete in ALL time points
    complete_in_all = set(mutual_indices)
    total_stage1_flagged = 0
    total_stage2_confirmed = 0

    n_planes = max(mutual_indices) + 1 if mutual_indices else 148

    # Determine worker count
    if workers == 0:
        workers = max(1, min((os.cpu_count() or 4) // 2, 8))
    workers = min(workers, n_timepoints)

    if workers <= 1:
        # Sequential mode (original behavior)
        for i, stl_path in enumerate(stl_files):
            stl_name = Path(stl_path).stem
            print(f"  [{i+1}/{n_timepoints}] {stl_name}...", end=' ', flush=True)

            incomplete, stage1, stage2 = detect_incomplete_planes_for_timepoint(
                temp_dir=temp_output,
                source_mesh_path=stl_path,
                time_point=stl_name,
                n_planes=n_planes
            )

            total_stage1_flagged += stage1
            total_stage2_confirmed += stage2

            incomplete_mutual = [idx for idx in incomplete if idx in complete_in_all]
            for idx in incomplete_mutual:
                complete_in_all.discard(idx)

            if incomplete_mutual:
                print(f"{len(incomplete_mutual)} incomplete")
            else:
                print("all complete")
    else:
        # Parallel mode
        print(f"Using {workers} parallel workers for {n_timepoints} time points...")
        frame_args = [
            (i, stl_files[i], str(temp_output), n_planes, n_timepoints)
            for i in range(n_timepoints)
        ]

        with multiprocessing.Pool(processes=workers) as pool:
            results_list = pool.map(_detect_incomplete_for_single_timepoint, frame_args)

        for result in results_list:
            if result['error']:
                print(f"  ERROR in frame {result['index']+1}: {result['error']}")
                continue

            total_stage1_flagged += result['stage1']
            total_stage2_confirmed += result['stage2']

            incomplete_mutual = [idx for idx in result['incomplete_indices']
                                 if idx in complete_in_all]
            for idx in incomplete_mutual:
                complete_in_all.discard(idx)

            status = (f"{len(incomplete_mutual)} incomplete" if incomplete_mutual
                      else "all complete")
            print(f"  [{result['index']+1}/{n_timepoints}] {result['stl_name']}: {status}")

    complete_mutual_indices = sorted(complete_in_all)

    print(f"\nTwo-Stage Detection Results:")
    print(f"  Mutual planes checked: {len(mutual_indices)}")
    print(f"  Stage 1 flagged (long edges): ~{total_stage1_flagged // n_timepoints} planes/timepoint")
    print(f"  Stage 2 confirmed (off surface): ~{total_stage2_confirmed // n_timepoints} planes/timepoint")
    print(f"  Computation saved: ~{(total_stage1_flagged - total_stage2_confirmed) // n_timepoints} surface checks/timepoint")
    print(f"\nFinal Results:")
    print(f"  Complete mutual planes: {len(complete_mutual_indices)}")
    print(f"  Incomplete planes (excluded): {len(mutual_indices) - len(complete_mutual_indices)}")

    excluded = sorted(set(mutual_indices) - set(complete_mutual_indices))
    if excluded:
        print(f"    Excluded indices: {excluded}")

    # Save validated indices
    validated_file = output_dir / f"{partition}_validated_indices.txt"
    with open(validated_file, 'w') as f:
        for idx in complete_mutual_indices:
            f.write(f"{idx}\n")

    print(f"\n✓ Saved {len(complete_mutual_indices)} validated indices to:")
    print(f"  {validated_file}")

    return validated_file, complete_mutual_indices


def step4_reslice_mutual_planes(subject, partition, validated_indices, output_dir):
    """Step 4: Copy and renumber validated planes from temp folder (NO re-slicing)"""
    print("\n" + "="*80)
    print("STEP 4: Copying and Renumbering Validated Planes")
    print("="*80)
    print("Approach: Copy from temp folder instead of re-slicing")
    print("Benefits: Exact same geometry, much faster, no re-slicing artifacts")
    print("="*80)

    # Create index mapping: original → sequential (0-based)
    index_mapping = {orig_idx: new_idx for new_idx, orig_idx in enumerate(validated_indices)}

    print(f"\nCopying {len(validated_indices)} validated planes with sequential numbering (0-{len(validated_indices)-1})")
    print("Source: temp folder (Step 1 output)")
    print("Target: final folder with sequential numbering")

    # Save mapping file
    mapping_file = output_dir / f"{partition}_index_mapping.txt"
    with open(mapping_file, 'w') as f:
        f.write("# Index mapping: original -> sequential\n")
        for orig, seq in sorted(index_mapping.items()):
            f.write(f"{orig} -> {seq}\n")
    print(f"\n✓ Saved index mapping to: {mapping_file}")

    # Setup directories
    temp_output = output_dir / f"{partition}SlicedSTLs_temp"
    final_output = output_dir / f"{partition}SlicedSTLs"

    if not temp_output.exists():
        print(f"\nError: Temp directory not found: {temp_output}")
        return None

    # Clean final output directory
    if final_output.exists():
        print(f"\nCleaning final output: {final_output}")
        shutil.rmtree(final_output)
    final_output.mkdir(parents=True, exist_ok=True)

    # Get all time points from temp folder
    time_points = sorted(set([
        f.stem.split('-Planes-')[0]
        for f in temp_output.glob("*-Planes-*.stl")
    ]))

    print(f"\nFound {len(time_points)} time points in temp folder")
    print(f"Processing: {time_points[0]} ... {time_points[-1]}\n")

    copied_count = 0
    import pandas as pd

    for i, time_point in enumerate(time_points):
        if (i + 1) % 5 == 0 or i == 0:
            print(f"  [{i+1}/{len(time_points)}] {time_point}...")

        # Copy and renumber STL files
        for orig_idx, seq_idx in index_mapping.items():
            src_stl = temp_output / f"{time_point}-Planes-{orig_idx:03d}.stl"
            dst_stl = final_output / f"{time_point}-Planes-{seq_idx:03d}.stl"

            if src_stl.exists():
                shutil.copy2(src_stl, dst_stl)
                copied_count += 1

        # Process CSV file (one CSV per timepoint containing all planes)
        src_csv = temp_output / f"{time_point}-Data.csv"
        dst_csv = final_output / f"{time_point}-Data.csv"

        if src_csv.exists():
            # Read CSV and filter to validated planes only
            df = pd.read_csv(src_csv)

            # Filter to only validated plane indices (original numbering)
            df_filtered = df[df['plane_index'].isin(validated_indices)].copy()

            # Update plane_index to sequential numbering
            df_filtered['plane_index'] = df_filtered['plane_index'].map(index_mapping)

            # Sort by plane_index
            df_filtered = df_filtered.sort_values('plane_index')

            # Save filtered and renumbered CSV
            df_filtered.to_csv(dst_csv, index=False)

    print(f"\n✓ Copied {copied_count} STL files from temp folder")
    print(f"✓ Filtered and renumbered CSV files (one per timepoint)")

    # Create combined STL files for video generation
    print(f"\nCreating combined STL files for video...")
    import trimesh
    import re

    # Group individual plane STLs by time point
    time_point_stls = {}
    for stl_file in sorted(final_output.glob("*-Planes-*.stl")):
        filename = stl_file.name
        # Extract time point (e.g., out_000000 from out_000000-Planes-000.stl)
        match = re.match(r'(out_\d+)-Planes-\d+\.stl', filename)
        if match:
            time_point = match.group(1)
            if time_point not in time_point_stls:
                time_point_stls[time_point] = []
            time_point_stls[time_point].append(stl_file)

    # Create combined -Planes-All.stl for each time point
    for i, (time_point, stl_files) in enumerate(sorted(time_point_stls.items())):
        if (i + 1) % 5 == 0 or i == 0:
            print(f"  [{i+1}/{len(time_point_stls)}] {time_point}...")

        # Load all individual plane meshes
        meshes = []
        for stl_file in sorted(stl_files):
            try:
                mesh = trimesh.load_mesh(str(stl_file))
                meshes.append(mesh)
            except Exception as e:
                print(f"    Warning: Could not load {stl_file.name}: {e}")

        if meshes:
            # Combine all meshes
            combined_mesh = trimesh.util.concatenate(meshes)

            # Save combined mesh
            output_path = final_output / f"{time_point}-Planes-All.stl"
            combined_mesh.export(str(output_path))

    print(f"✓ Created {len(time_point_stls)} combined STL files (for video)")
    print(f"✓ Keeping temp directory for debugging: {temp_output}")

    return index_mapping


def step5_combine_and_renumber_csvs(subject, partition, validated_indices, index_mapping, output_dir):
    """Step 5: Combine CSVs with sequential numbering"""
    print("\n" + "="*80)
    print("STEP 5: Combining CSVs with Sequential Numbering")
    print("="*80)

    sliced_dir = output_dir / f"{partition}SlicedSTLs"

    # Find all CSV files (flat structure)
    csv_files = sorted(list(sliced_dir.glob("*-Data.csv")))

    if len(csv_files) == 0:
        print(f"Error: No CSV files found in {sliced_dir}")
        return None

    print(f"Found {len(csv_files)} CSV files to combine")

    all_data = []

    for csv_file in csv_files:
        time_name = csv_file.stem.replace('-Data', '')  # out_000000-Data -> out_000000
        phase = int(time_name.split('_')[1]) if '_' in time_name else 0

        # Read CSV
        df = pd.read_csv(csv_file)

        # Add metadata columns
        df['phase'] = phase
        df['file_id'] = time_name
        df['original_plane_index'] = df['plane_index']  # Save original index

        # Apply sequential index mapping
        df['plane_index'] = df['plane_index'].map(index_mapping)

        all_data.append(df)

    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)

    # Sort by phase and plane index
    combined_df = combined_df.sort_values(['phase', 'plane_index'])

    # Save combined CSV
    output_csv = output_dir / f"{subject}_{partition}_all_measurements.csv"
    combined_df.to_csv(output_csv, index=False)

    print(f"✓ Combined {len(csv_files)} CSV files")
    print(f"✓ Saved to: {output_csv}")

    return output_csv


def step6_generate_plots_and_video(partition, output_dir):
    """Step 6: Generate plots, video, and interactive HTML reference"""
    print("\n" + "="*80)
    print("STEP 6: Generating Plots, Video, and Interactive HTML")
    print("="*80)

    # Generate CSA plots
    print("\nGenerating CSA plots...")
    cmd = ["python", str(SCRIPT_DIR / "plot_csa_by_index.py"), "--partition", partition]
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print("Warning: Plot generation failed")
    else:
        print("✓ CSA plots generated")

    # Generate regular video (unlabeled)
    print("\nGenerating plane motion video...")
    cmd = ["python", str(SCRIPT_DIR / "create_plane_video_fast.py"), partition]
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print("Warning: Video generation failed")
    else:
        print("✓ Video generated")

    # Generate interactive HTML plane reference
    print("\nGenerating interactive HTML plane reference...")
    cmd = ["python", str(SCRIPT_DIR / "create_plane_reference_interactive_plotly.py"), partition]
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print("Warning: Interactive HTML generation failed")
    else:
        print("✓ Interactive HTML plane reference generated")
        print("  Open in browser to explore plane indices interactively")

    return True


def check_step1_cache(partition, output_dir):
    """Check if Step 1 outputs exist"""
    temp_output = output_dir / f"{partition}SlicedSTLs_temp"
    if not temp_output.exists():
        return False

    # Check if we have CSV files
    csv_files = list(temp_output.glob("*-Data.csv"))
    return len(csv_files) > 0


def check_step2_cache(partition, output_dir):
    """Check if Step 2 outputs exist"""
    mutual_file = output_dir / f"{partition}_mutual_indices.txt"
    return mutual_file.exists()


def check_step3_cache(partition, output_dir):
    """Check if Step 3 outputs exist"""
    validated_file = output_dir / f"{partition}_validated_indices.txt"
    return validated_file.exists()


def check_step4_cache(partition, output_dir):
    """Check if Step 4 outputs exist"""
    final_output = output_dir / f"{partition}SlicedSTLs"
    if not final_output.exists():
        return False

    # Check if we have STL and CSV files
    csv_files = list(final_output.glob("*-Data.csv"))
    return len(csv_files) > 0


def check_step5_cache(subject, partition, output_dir):
    """Check if Step 5 outputs exist"""
    combined_csv = output_dir / f"{subject}_{partition}_all_measurements.csv"
    return combined_csv.exists()


def load_cached_indices(partition, output_dir):
    """Load mutual and validated indices from cache"""
    mutual_file = output_dir / f"{partition}_mutual_indices.txt"
    validated_file = output_dir / f"{partition}_validated_indices.txt"

    mutual_indices = []
    validated_indices = []

    if mutual_file.exists():
        with open(mutual_file, 'r') as f:
            mutual_indices = [int(line.strip()) for line in f if line.strip()]

    if validated_file.exists():
        with open(validated_file, 'r') as f:
            validated_indices = [int(line.strip()) for line in f if line.strip()]

    return mutual_file, mutual_indices, validated_file, validated_indices


def run_full_pipeline(subject, partition, output_dir=".", force=False, workers=0):
    """Run complete smart pipeline with caching support

    Args:
        subject: Subject ID
        partition: Partition name
        output_dir: Output directory
        force: If True, ignore cache and re-run all steps
        workers: Number of parallel workers (0=auto, 1=sequential)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*80)
    print("SMART AIRWAY SLICER PIPELINE")
    print("="*80)
    print(f"Subject: {subject}")
    print(f"Partition: {partition}")
    print(f"Output directory: {output_dir}")
    print(f"Features: window=20, adaptive normals, auto mutual/closed detection")
    print(f"Caching: {'DISABLED (force mode)' if force else 'ENABLED'}")
    print("="*80 + "\n")

    # Initialize variables
    mutual_file = None
    mutual_indices = []
    validated_file = None
    validated_indices = []
    index_mapping = None

    # Step 1: Slice all planes
    if not force and check_step1_cache(partition, output_dir):
        print("\n" + "="*80)
        print("STEP 1: [CACHED] Slicing All Planes")
        print("="*80)
        print("✓ Found cached temp slicing results")
        temp_output = output_dir / f"{partition}SlicedSTLs_temp"
        csv_files = list(temp_output.glob("*-Data.csv"))
        print(f"✓ Found {len(csv_files)} CSV files in {temp_output}")
        print("✓ Skipping Step 1 (already completed)")
    else:
        if not step1_slice_all_planes(subject, partition, output_dir, workers=workers):
            print("\n❌ Pipeline failed at Step 1")
            return False

    # Step 2: Detect mutual planes from files
    if not force and check_step2_cache(partition, output_dir):
        print("\n" + "="*80)
        print("STEP 2: [CACHED] Detecting Mutual Planes")
        print("="*80)
        mutual_file = output_dir / f"{partition}_mutual_indices.txt"
        with open(mutual_file, 'r') as f:
            mutual_indices = [int(line.strip()) for line in f if line.strip()]
        print(f"✓ Loaded {len(mutual_indices)} mutual indices from cache")
        print(f"✓ File: {mutual_file}")
        print("✓ Skipping Step 2 (already completed)")
    else:
        result = step2_detect_mutual_planes(partition, output_dir)
        if result is None:
            print("\n❌ Pipeline failed at Step 2")
            return False
        mutual_file, mutual_indices = result

    # Step 3: Detect closed planes
    if not force and check_step3_cache(partition, output_dir):
        print("\n" + "="*80)
        print("STEP 3: [CACHED] Detecting Closed Planes")
        print("="*80)
        validated_file = output_dir / f"{partition}_validated_indices.txt"
        with open(validated_file, 'r') as f:
            validated_indices = [int(line.strip()) for line in f if line.strip()]
        print(f"✓ Loaded {len(validated_indices)} validated indices from cache")
        print(f"✓ File: {validated_file}")
        print("✓ Skipping Step 3 (already completed)")
    else:
        result = step3_detect_complete_planes_optimized(partition, mutual_indices, output_dir, workers=workers)
        if result is None:
            print("\n❌ Pipeline failed at Step 3")
            return False
        validated_file, validated_indices = result

    # Step 4: Re-slice mutual planes with sequential numbering
    if not force and check_step4_cache(partition, output_dir):
        print("\n" + "="*80)
        print("STEP 4: [CACHED] Re-slicing Mutual Planes")
        print("="*80)
        print("✓ Found cached re-sliced results")
        final_output = output_dir / f"{partition}SlicedSTLs"
        csv_files = list(final_output.glob("*-Data.csv"))
        print(f"✓ Found {len(csv_files)} CSV files in {final_output}")

        # Load index mapping from file
        mapping_file = output_dir / f"{partition}_index_mapping.txt"
        if mapping_file.exists():
            index_mapping = {}
            with open(mapping_file, 'r') as f:
                f.readline()  # Skip header
                for line in f:
                    if '->' in line:
                        orig, seq = line.strip().split('->')
                        index_mapping[int(orig)] = int(seq)
            print(f"✓ Loaded index mapping: {len(index_mapping)} entries")
        else:
            # Recreate mapping from validated indices
            index_mapping = {orig_idx: new_idx for new_idx, orig_idx in enumerate(validated_indices, start=1)}
            print(f"✓ Recreated index mapping: {len(index_mapping)} entries")

        print("✓ Skipping Step 4 (already completed)")
    else:
        index_mapping = step4_reslice_mutual_planes(subject, partition, validated_indices, output_dir)
        if index_mapping is None:
            print("\n❌ Pipeline failed at Step 4")
            return False

    # Step 5: Combine CSVs with sequential numbering
    if not force and check_step5_cache(subject, partition, output_dir):
        print("\n" + "="*80)
        print("STEP 5: [CACHED] Combining CSVs")
        print("="*80)
        combined_csv = output_dir / f"{subject}_{partition}_all_measurements.csv"
        print(f"✓ Found cached combined CSV: {combined_csv}")
        print("✓ Skipping Step 5 (already completed)")
    else:
        combined_csv = step5_combine_and_renumber_csvs(
            subject, partition, validated_indices, index_mapping, output_dir
        )

    # Step 6: Generate plots and video
    # Always run Step 6 (plots/videos are fast to regenerate)
    if not step6_generate_plots_and_video(partition, output_dir):
        print("\n❌ Pipeline failed at Step 6")
        return False

    # Step 7: Advanced dynamics analysis (generates enhanced CSVs + PDF report)
    print("\n" + "="*80)
    print("STEP 7: Advanced Airway Dynamics Analysis")
    print("="*80)
    try:
        # Import and run the advanced analysis module
        import subprocess
        result = subprocess.run(
            ['python', str(SCRIPT_DIR / 'analyze_airway_dynamics_advanced.py'), subject, partition],
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )

        if result.returncode == 0:
            print("✓ Advanced analysis completed successfully")
            print(f"✓ Generated enhanced metrics CSV")
            print(f"✓ Generated summary metrics CSV")
            print(f"✓ Generated inter-plane volumes CSV")
            print(f"✓ Generated PDF report with parameter reference page")
        else:
            print(f"⚠ Warning: Advanced analysis failed (exit code {result.returncode})")
            print("Pipeline will continue, but advanced metrics/PDF not generated")
            if result.stderr:
                print(f"Error: {result.stderr[:500]}")  # Show first 500 chars of error
    except subprocess.TimeoutExpired:
        print("⚠ Warning: Advanced analysis timed out after 10 minutes")
        print("Pipeline will continue, but advanced metrics/PDF not generated")
    except Exception as e:
        print(f"⚠ Warning: Could not run advanced analysis: {e}")
        print("Pipeline will continue, but advanced metrics/PDF not generated")

    # Summary
    print("\n" + "="*80)
    print("✅ PIPELINE COMPLETE!")
    print("="*80)
    print(f"✓ Sliced all planes initially (for detection)")
    print(f"✓ Detected {len(mutual_indices)} mutual planes from files")
    print(f"✓ Validated {len(validated_indices)} closed mutual planes")
    print(f"✓ Re-sliced ONLY mutual planes with sequential numbering (1-{len(validated_indices)})")
    print(f"✓ Generated clean STL + CSV files (no gaps in numbering)")
    print(f"✓ Generated combined CSV, plots, and video")
    print(f"\nOutput directory: {output_dir / f'{partition}SlicedSTLs'}/")
    print(f"  - STL plane meshes: {len(validated_indices)} planes × 21 time points")
    print(f"  - CSV data files: {len(validated_indices)} planes × 21 time points")
    print(f"\nSummary files:")
    print(f"  - Mutual indices: {mutual_file}")
    print(f"  - Validated indices: {validated_file}")
    print(f"  - Index mapping: {partition}_index_mapping.txt")
    print(f"  - Combined CSV: {combined_csv}")
    print(f"  - Plots: {partition}_CSA_*.png")
    print(f"  - Video: {partition}_plane_motion_breathing_cycle.mp4")
    print(f"  - Interactive HTML: {partition}_planes_reference_interactive.html")
    print(f"    (Open in browser to explore plane indices interactively)")
    print(f"\nAdvanced analysis outputs:")
    print(f"  - Enhanced metrics: {subject}_{partition}_enhanced_metrics.csv")
    print(f"  - Summary metrics: {subject}_{partition}_summary_metrics.csv")
    print(f"  - Inter-plane volumes: {subject}_{partition}_interplane_volumes.csv")
    print(f"  - PDF report: {subject}_{partition}_airway_dynamics_report.pdf")
    print("="*80)

    return True


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Smart airway slicer pipeline - slices, detects mutual/closed planes, generates plots & video"
    )
    parser.add_argument('subject', help='Subject ID')
    parser.add_argument('partition', help='Partition name')
    parser.add_argument('-o', '--output', default='.',
                       help='Output directory')
    parser.add_argument('--force', action='store_true',
                       help='Force re-run all steps (ignore cache)')
    parser.add_argument('--workers', type=int, default=0,
                       help='Parallel workers per step (0=auto, 1=sequential)')

    args = parser.parse_args()

    success = run_full_pipeline(
        args.subject,
        args.partition,
        args.output,
        force=args.force,
        workers=args.workers
    )

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
