#!/usr/bin/env python3
"""
Master orchestration script to run pipeline on all partitions
Validates inputs, then runs pipeline_full.py on all partitions in parallel
"""

import os
import sys
import glob
import subprocess
import time
from pathlib import Path
from datetime import datetime


def validate_partition(partition_name):
    """
    Validate that a partition has required input files

    Returns:
        (bool, list, dict): (is_valid, list_of_issues, file_counts)
    """
    issues = []
    file_counts = {'stl': 0, 'vtk': 0}

    partition_dir = Path(partition_name)

    # Check if partition directory exists
    if not partition_dir.exists():
        issues.append(f"Partition directory does not exist: {partition_name}")
        return False, issues, file_counts

    # Check FFD directory
    ffd_dir = partition_dir / "FFD"
    if not ffd_dir.exists():
        issues.append(f"FFD directory missing: {ffd_dir}")
        return False, issues, file_counts

    # Check STL directory
    stl_dir = ffd_dir / "stl"
    if not stl_dir.exists():
        issues.append(f"STL directory missing: {stl_dir}")
    else:
        stl_files = sorted(glob.glob(str(stl_dir / "*.stl")))
        file_counts['stl'] = len(stl_files)
        if len(stl_files) == 0:
            issues.append(f"No STL files found in {stl_dir}")

    # Check VTK directory
    vtk_dir = ffd_dir / "vtk"
    if not vtk_dir.exists():
        issues.append(f"VTK directory missing: {vtk_dir}")
    else:
        vtk_files = sorted(glob.glob(str(vtk_dir / "*.vtk")))
        file_counts['vtk'] = len(vtk_files)
        if len(vtk_files) == 0:
            issues.append(f"No VTK files found in {vtk_dir}")

    # Check if STL and VTK counts match (if both exist)
    if file_counts['stl'] > 0 and file_counts['vtk'] > 0:
        if file_counts['stl'] != file_counts['vtk']:
            issues.append(f"Mismatch: {file_counts['stl']} STL files but {file_counts['vtk']} VTK files")

    is_valid = len(issues) == 0
    return is_valid, issues, file_counts


def run_pipeline(patient_id, partition_name, log_file, force=False, workers=0):
    """
    Run pipeline_full.py for a single partition

    Args:
        force: If True, bypass all cache and regenerate everything
        workers: Number of parallel workers per step (0=auto, 1=sequential)

    Returns:
        subprocess.Popen: The process object
    """
    script_dir = Path(__file__).parent / "python_slicer"
    cmd = [
        "python",
        str(script_dir / "pipeline_full.py"),
        patient_id,
        partition_name
    ]

    if force:
        cmd.append("--force")
    if workers != 0:
        cmd.extend(["--workers", str(workers)])

    log_path = Path(log_file)
    with open(log_path, 'w') as f:
        f.write(f"Starting pipeline at {datetime.now()}\n")
        f.write(f"Command: {' '.join(cmd)}\n")
        f.write("="*80 + "\n\n")

    # Run in background, append to log file
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1
    )

    return process, log_path


def monitor_processes(processes):
    """
    Monitor multiple processes and write their output to log files

    Args:
        processes: list of (process, log_path, partition_name) tuples
    """
    while any(p[0].poll() is None for p in processes):
        for process, log_path, partition_name in processes:
            if process.poll() is None:
                # Process still running, read available output
                try:
                    line = process.stdout.readline()
                    if line:
                        with open(log_path, 'a') as f:
                            f.write(line)
                        # Print to console with partition prefix
                        print(f"[{partition_name}] {line.rstrip()}")
                except:
                    pass
        time.sleep(0.1)

    # Read remaining output
    for process, log_path, partition_name in processes:
        remaining = process.stdout.read()
        if remaining:
            with open(log_path, 'a') as f:
                f.write(remaining)
            for line in remaining.split('\n'):
                if line.strip():
                    print(f"[{partition_name}] {line}")


def main():
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Run airway slicing pipeline on all partitions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_all_partitions.py OSAMRI037
  python run_all_partitions.py OSAMRI042
        """
    )
    parser.add_argument('subject_id', type=str,
                       help='Subject ID (e.g., OSAMRI037)')
    parser.add_argument('--partitions', nargs='+',
                       default=["LeftNoseDecending", "RightNose"],
                       help='Partitions to process (default: LeftNoseDecending RightNose)')
    parser.add_argument('--force', action='store_true',
                       help='Force regeneration, bypass all cache')
    parser.add_argument('--workers', type=int, default=0,
                       help='Parallel workers per partition (0=auto, 1=sequential). '
                            'Note: partitions already run in parallel, so total CPU = partitions x workers')

    args = parser.parse_args()

    print("="*80)
    print("MASTER PIPELINE - Run All Partitions")
    print("="*80)
    print()

    # Configuration from arguments
    PATIENT_ID = args.subject_id
    PARTITIONS = args.partitions

    print(f"Patient ID: {PATIENT_ID}")
    print(f"Partitions to process: {', '.join(PARTITIONS)}")
    print()

    # Step 1: Validate all partitions
    print("="*80)
    print("STEP 1: VALIDATING INPUT FILES")
    print("="*80)
    print()

    validation_results = {}
    all_valid = True

    for partition in PARTITIONS:
        print(f"Checking partition: {partition}")
        is_valid, issues, file_counts = validate_partition(partition)
        validation_results[partition] = (is_valid, issues, file_counts)

        if is_valid:
            print(f"  ✓ VALID - {file_counts['stl']} STL files, {file_counts['vtk']} VTK files")
        else:
            print(f"  ✗ INVALID")
            for issue in issues:
                print(f"    - {issue}")
            # Show file counts even if invalid (helps with debugging)
            if file_counts['stl'] > 0 or file_counts['vtk'] > 0:
                print(f"    Found: {file_counts['stl']} STL files, {file_counts['vtk']} VTK files")
            all_valid = False
        print()

    # If any partition is invalid, stop
    if not all_valid:
        print("="*80)
        print("VALIDATION FAILED")
        print("="*80)
        print()
        print("Cannot proceed. Please fix the issues above and try again.")
        print()

        # Summary of what's missing
        print("ISSUES SUMMARY:")
        for partition, (is_valid, issues, file_counts) in validation_results.items():
            if not is_valid:
                print(f"\n{partition}:")
                for issue in issues:
                    print(f"  - {issue}")

        sys.exit(1)

    # Step 2: Run pipeline on all partitions in parallel
    print("="*80)
    print("STEP 2: RUNNING PIPELINE ON ALL PARTITIONS")
    print("="*80)
    print()
    print(f"Starting {len(PARTITIONS)} pipelines in parallel...")
    print("This will take approximately 5-10 minutes per partition.")
    print()

    # Start all processes
    processes = []
    for partition in PARTITIONS:
        log_file = f"pipeline_{partition}.log"
        print(f"Starting: {partition}")
        print(f"  Log file: {log_file}")

        process, log_path = run_pipeline(PATIENT_ID, partition, log_file, args.force, args.workers)
        processes.append((process, log_path, partition))

    print()
    print("All pipelines started. Monitoring progress...")
    print("="*80)
    print()

    start_time = time.time()

    # Monitor processes
    monitor_processes(processes)

    # Check results
    print()
    print("="*80)
    print("PIPELINE RESULTS")
    print("="*80)
    print()

    all_success = True
    for process, log_path, partition in processes:
        return_code = process.returncode
        if return_code == 0:
            print(f"✓ {partition}: SUCCESS")
        else:
            print(f"✗ {partition}: FAILED (exit code {return_code})")
            print(f"  Check log: {log_path}")
            all_success = False

    elapsed_time = time.time() - start_time
    print()
    print(f"Total time: {elapsed_time/60:.1f} minutes")
    print()

    # Summary of outputs
    if all_success:
        print("="*80)
        print("ALL PIPELINES COMPLETED SUCCESSFULLY")
        print("="*80)
        print()
        print("Generated outputs:")

        for partition in PARTITIONS:
            print(f"\n{partition}:")

            # Check for expected outputs
            outputs = {
                f"{partition}SlicedSTLs/": "Sliced plane STL files",
                f"{PATIENT_ID}_{partition}_all_measurements.csv": "Combined measurements CSV",
                f"{partition}_CSA_flipped.png": "Flipped CSA plot",
                f"{partition}_plane_motion_breathing_cycle.mp4": "Video",
                f"{partition}_index_mapping.txt": "Index mapping"
            }

            for output_path, description in outputs.items():
                if Path(output_path).exists():
                    print(f"  ✓ {description}")
                    print(f"    {output_path}")
                else:
                    print(f"  ✗ {description} (not found)")
                    print(f"    Expected: {output_path}")

        print()
        print("="*80)
        sys.exit(0)
    else:
        print("="*80)
        print("SOME PIPELINES FAILED")
        print("="*80)
        print()
        print("Check the log files above for details.")
        print()
        sys.exit(1)


if __name__ == "__main__":
    main()
