#!/usr/bin/env python3
"""
Step 0: Automatic Airway Branch Detection

Standalone entry point that runs branch_detector.py inside the 'vmtk' conda
environment via `conda run -n vmtk`. This way the user does not need to
manually activate the vmtk env.

Usage:
    python python_slicer/step0_branch_detect.py /path/to/full_airway.stl
    python python_slicer/step0_branch_detect.py /path/to/full_airway.stl --output-dir ./branches
    python python_slicer/step0_branch_detect.py /path/to/full_airway.stl --subject-id OSAMRI037
"""

import sys
import subprocess
import shutil
from pathlib import Path


VMTK_ENV_NAME = "vmtk"
BRANCH_DETECTOR_SCRIPT = Path(__file__).parent / "branch_detector.py"


def check_vmtk_env():
    """Check if the vmtk conda environment exists"""
    result = subprocess.run(
        ["conda", "env", "list", "--json"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("ERROR: conda not found. Please install conda first.")
        sys.exit(1)

    import json
    envs = json.loads(result.stdout).get("envs", [])
    for env_path in envs:
        if Path(env_path).name == VMTK_ENV_NAME:
            return True

    return False


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Step 0: Automatic airway branch detection using VMTK",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python python_slicer/step0_branch_detect.py /path/to/full_airway.stl
  python python_slicer/step0_branch_detect.py airway.stl --output-dir ./branches --subject-id OSAMRI037
  python python_slicer/step0_branch_detect.py airway.stl --verbose

This script runs inside the 'vmtk' conda environment automatically.
To create the vmtk environment:
  conda create -n vmtk -c conda-forge python=3.10 vmtk trimesh pyvista shapely scipy pandas networkx rtree matplotlib pytest
        """
    )
    parser.add_argument("stl_path", help="Path to full airway STL file")
    parser.add_argument("--output-dir", default="./branch_output",
                        help="Output directory (default: ./branch_output)")
    parser.add_argument("--subject-id", default="",
                        help="Subject ID for naming output files")
    parser.add_argument("--resampling-step", type=float, default=2.0,
                        help="Centerline resampling step length (default: 2.0)")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose output")
    args = parser.parse_args()

    # Validate input
    stl_path = Path(args.stl_path).resolve()
    if not stl_path.exists():
        print(f"ERROR: STL file not found: {stl_path}")
        sys.exit(1)

    # Check vmtk env
    print(f"Checking for '{VMTK_ENV_NAME}' conda environment...")
    if not check_vmtk_env():
        print(f"ERROR: Conda environment '{VMTK_ENV_NAME}' not found.")
        print(f"Create it with:")
        print(f"  conda create -n {VMTK_ENV_NAME} -c conda-forge python=3.10 vmtk "
              f"trimesh pyvista shapely scipy pandas networkx rtree matplotlib pytest")
        sys.exit(1)
    print(f"  Found '{VMTK_ENV_NAME}' environment.")

    # Build command to run branch_detector.py inside vmtk env
    cmd = [
        "conda", "run", "-n", VMTK_ENV_NAME,
        "python", str(BRANCH_DETECTOR_SCRIPT),
        str(stl_path),
        "--output-dir", str(Path(args.output_dir).resolve()),
    ]
    if args.subject_id:
        cmd.extend(["--subject-id", args.subject_id])
    if args.resampling_step != 2.0:
        cmd.extend(["--resampling-step", str(args.resampling_step)])
    if args.verbose:
        cmd.append("--verbose")

    print(f"\nRunning branch detection...")
    print(f"  Input: {stl_path}")
    print(f"  Output: {Path(args.output_dir).resolve()}")
    print(f"  Environment: {VMTK_ENV_NAME}")
    print()

    # Run with real-time output
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    for line in process.stdout:
        print(line, end="")

    process.wait()

    if process.returncode == 0:
        print(f"\nBranch detection completed successfully.")
        print(f"Results in: {Path(args.output_dir).resolve()}")
    else:
        print(f"\nERROR: Branch detection failed (exit code {process.returncode})")
        sys.exit(process.returncode)


if __name__ == "__main__":
    main()
