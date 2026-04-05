#!/usr/bin/env python3
"""
Interpolate deformed frames from FFD registration data.

Takes a reference mesh (frame 0) + FFD DOF files and generates N deformed
frames across the breathing cycle using MIRTK transform-points + cubic
spline interpolation.

Handles both STL surfaces and VTK/VTP centerlines.

Requires: MIRTK installed (mirtk command-line tool)

Usage:
    python python_slicer/interpolate_frames.py ENT001/
    python python_slicer/interpolate_frames.py ENT001/ --step 100 --start 0 --stop 2000
"""

import argparse
import csv
import logging
import os
import sys
import numpy as np
import scipy.interpolate
from pathlib import Path
from subprocess import check_call, CalledProcessError
from tempfile import mkstemp

import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Mesh I/O
# ---------------------------------------------------------------------------

def read_stl_vtk(path: str) -> vtk.vtkPolyData:
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(path))
    reader.Update()
    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(reader.GetOutput())
    return mesh


def write_stl_vtk(mesh: vtk.vtkPolyData, path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(mesh)
    writer.SetFileTypeToBinary()
    writer.Update()


def read_polydata(path: str) -> vtk.vtkPolyData:
    """Read VTK or VTP polydata."""
    path = str(path)
    if path.endswith(".vtp"):
        reader = vtk.vtkXMLPolyDataReader()
    else:
        reader = vtk.vtkPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(reader.GetOutput())
    return mesh


def write_polydata(mesh: vtk.vtkPolyData, path: str):
    """Write VTK polydata (legacy format)."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(mesh)
    writer.SetFileTypeToASCII()
    writer.Write()


def replace_mesh_points(mesh: vtk.vtkPolyData, points: np.ndarray) -> vtk.vtkPolyData:
    """Create new mesh with same topology but new point coordinates."""
    new_points = vtk.vtkPoints()
    new_points.SetData(numpy_to_vtk(points))
    new_mesh = vtk.vtkPolyData()
    new_mesh.ShallowCopy(mesh)
    new_mesh.SetPoints(new_points)
    return new_mesh


# ---------------------------------------------------------------------------
# DOF handling
# ---------------------------------------------------------------------------

def read_dofs(dofs_csv_path: Path) -> list:
    """Read DOF file paths and timepoints from CSV."""
    dofs = []
    with open(dofs_csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            dofs.append((str(Path(row["dof"])), float(row["t"])))
    return sorted(dofs, key=lambda x: x[1])


def compute_time_range(input_txt: Path) -> dict:
    """Read timing info from input.txt."""
    info = {}
    with open(input_txt) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            if "," in line:
                key, val = line.split(",", 1)
                info[key.strip()] = val.strip()
    return info


# ---------------------------------------------------------------------------
# Core interpolation
# ---------------------------------------------------------------------------

def interpolate_mesh(
    reference_mesh: vtk.vtkPolyData,
    reference_path: str,
    dofs: list,
    dofs_dir: Path,
    start: float,
    stop: float,
    step: float,
    output_template: str,
    mesh_type: str = "stl",
):
    """Deform a reference mesh at DOF timepoints, then interpolate.

    Args:
        reference_mesh: VTK polydata of the reference mesh
        reference_path: Path to reference mesh file (for mirtk)
        dofs: List of (dof_path, time) tuples
        dofs_dir: Directory containing DOF files (for relative paths)
        start, stop, step: Time range for output frames
        output_template: Path template like "out_{t:06.0f}.stl"
        mesh_type: "stl" or "vtk" (determines temp file suffix and I/O)
    """
    suffix = ".stl" if mesh_type == "stl" else ".vtk"

    time_points = []
    mesh_points = []

    for dof_path, t in dofs:
        log.info(f"  Deforming to t={t:.0f}ms using {dof_path}")
        time_points.append(t)

        if dof_path in ("id", "Id", "identity"):
            points = vtk_to_numpy(reference_mesh.GetPoints().GetData())
        else:
            # Resolve DOF path relative to dofs_dir
            full_dof_path = str(dofs_dir / dof_path) if not os.path.isabs(dof_path) else dof_path

            fp, tmp_path = mkstemp(suffix=suffix)
            os.close(fp)
            try:
                check_call([
                    "mirtk", "transform-points",
                    str(reference_path), tmp_path,
                    "-dofin", full_dof_path,
                ])
                if mesh_type == "stl":
                    deformed = read_stl_vtk(tmp_path)
                else:
                    deformed = read_polydata(tmp_path)
                points = vtk_to_numpy(deformed.GetPoints().GetData())
            finally:
                os.remove(tmp_path)

        mesh_points.append(points)

    # Stack and build cubic spline interpolator
    mesh_points = np.stack(mesh_points, axis=0)
    interpolator = scipy.interpolate.interp1d(
        time_points, mesh_points, axis=0,
        bounds_error=True, assume_sorted=True, kind="cubic",
    )

    # Generate output frames
    ts = np.arange(start, stop + step, step)
    output_paths = []
    for t in ts:
        points = interpolator(t)
        mesh = replace_mesh_points(reference_mesh, points)
        out_path = output_template.format(t=t)
        if mesh_type == "stl":
            write_stl_vtk(mesh, out_path)
        else:
            write_polydata(mesh, out_path)
        output_paths.append(out_path)

    log.info(f"  Generated {len(output_paths)} frames ({start}-{stop}ms, step={step}ms)")
    return output_paths


# ---------------------------------------------------------------------------
# High-level API
# ---------------------------------------------------------------------------

def generate_motion_frames(
    subject_dir: Path,
    start: float = 0,
    stop: float = 2000,
    step: float = 100,
):
    """Generate deformed STL + centerline frames for a subject.

    Expects:
        subject_dir/surface/frame0.stl
        subject_dir/registration/ffds.csv + ffd_*.dof.gz + img_0.nii.gz
        subject_dir/branches/LeftNose/*_centerline.vtp (or .vtk)
        subject_dir/branches/RightNose/*_centerline.vtp (or .vtk)

    Outputs:
        subject_dir/motion/stl/out_XXXXXX.stl (deformed surfaces)
        subject_dir/motion/centerlines/left_nose/out_XXXXXX.vtk
        subject_dir/motion/centerlines/right_nose/out_XXXXXX.vtk
        subject_dir/motion/centerlines/midline/out_XXXXXX.vtk
    """
    subject_dir = Path(subject_dir)
    reg_dir = subject_dir / "registration"
    branches_dir = subject_dir / "branches"
    motion_dir = subject_dir / "motion"

    # Read DOFs
    dofs_csv = reg_dir / "ffds.csv"
    if not dofs_csv.exists():
        raise FileNotFoundError(f"ffds.csv not found: {dofs_csv}")
    dofs = read_dofs(dofs_csv)
    log.info(f"Loaded {len(dofs)} DOF entries from {dofs_csv}")

    # Read timing info
    input_txt = reg_dir / "input.txt"
    if input_txt.exists():
        timing = compute_time_range(input_txt)
        bre_time = float(timing.get("breTime", stop))
        log.info(f"Breathing time: {bre_time}ms")
        # Use breTime as stop if provided
        if stop == 2000 and bre_time != 2000:
            stop = bre_time
            log.info(f"  Adjusted stop to {stop}ms from input.txt")

    # Check MIRTK available
    import subprocess
    try:
        subprocess.run(["mirtk", "help"], capture_output=True, check=False)
    except FileNotFoundError:
        raise RuntimeError("MIRTK not found. Install MIRTK and ensure 'mirtk' is in PATH.")

    # ── Deform surface STL ──────────────────────────────────
    # Find the single STL in surface/ (any name)
    surface_stls = list((subject_dir / "surface").glob("*.stl"))
    if not surface_stls:
        raise FileNotFoundError(f"No STL found in {subject_dir / 'surface'}")
    frame0_stl = surface_stls[0]

    stl_out_dir = motion_dir / "stl"
    stl_out_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"\nDeforming surface mesh ({frame0_stl.name})...")
    ref_mesh = read_stl_vtk(str(frame0_stl))
    log.info(f"  Reference: {ref_mesh.GetNumberOfPoints()} pts, {ref_mesh.GetNumberOfCells()} cells")

    interpolate_mesh(
        reference_mesh=ref_mesh,
        reference_path=str(frame0_stl),
        dofs=dofs,
        dofs_dir=reg_dir,
        start=start, stop=stop, step=step,
        output_template=str(stl_out_dir / "out_{t:06.0f}.stl"),
        mesh_type="stl",
    )

    # ── Deform centerlines ──────────────────────────────────
    # Find centerline files in branches/
    centerline_configs = []

    # Left nose
    for label in ["LeftNose", "RightNose"]:
        cl_dir = branches_dir / label
        if not cl_dir.exists():
            continue
        # Find the full-path centerline (Trachea_to_*)
        cl_files = list(cl_dir.glob("*Trachea_to_*_centerline.vtk")) + \
                   list(cl_dir.glob("*Trachea_to_*_centerline.vtp"))
        if cl_files:
            cl_path = cl_files[0]
            out_subdir = motion_dir / "centerlines" / label.lower()
            out_subdir.mkdir(parents=True, exist_ok=True)
            centerline_configs.append((label, cl_path, out_subdir))

    # Merged midline (if exists)
    midline_files = list(branches_dir.glob("*NoseMidline*centerline.vtk")) + \
                    list(branches_dir.glob("*NoseMidline*centerline.vtp"))
    if midline_files:
        ml_path = midline_files[0]
        out_subdir = motion_dir / "centerlines" / "midline"
        out_subdir.mkdir(parents=True, exist_ok=True)
        centerline_configs.append(("Midline", ml_path, out_subdir))

    for label, cl_path, out_subdir in centerline_configs:
        log.info(f"\nDeforming {label} centerline ({cl_path.name})...")
        ref_cl = read_polydata(str(cl_path))
        log.info(f"  Reference: {ref_cl.GetNumberOfPoints()} pts")

        # Determine file type from input
        suffix = cl_path.suffix  # .vtk or .vtp
        out_suffix = ".vtk"  # always output VTK for consistency

        interpolate_mesh(
            reference_mesh=ref_cl,
            reference_path=str(cl_path),
            dofs=dofs,
            dofs_dir=reg_dir,
            start=start, stop=stop, step=step,
            output_template=str(out_subdir / f"out_{{t:06.0f}}{out_suffix}"),
            mesh_type="vtk",
        )

    log.info(f"\n{'='*50}")
    log.info(f"Motion frames generated in: {motion_dir}")
    n_stl = len(list(stl_out_dir.glob("*.stl")))
    log.info(f"  STL frames: {n_stl}")
    for label, _, out_subdir in centerline_configs:
        n_cl = len(list(out_subdir.glob("*.vtk")))
        log.info(f"  {label} CL frames: {n_cl}")
    log.info(f"{'='*50}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate deformed frames from FFD registration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python python_slicer/interpolate_frames.py ENT001/
  python python_slicer/interpolate_frames.py ENT001/ --step 100 --start 0 --stop 2000
        """,
    )
    parser.add_argument("subject_dir", help="Subject directory (e.g., ENT001/)")
    parser.add_argument("--start", type=float, default=0, help="Start time (ms)")
    parser.add_argument("--stop", type=float, default=2000, help="Stop time (ms)")
    parser.add_argument("--step", type=float, default=100, help="Time step (ms)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    generate_motion_frames(
        subject_dir=Path(args.subject_dir),
        start=args.start,
        stop=args.stop,
        step=args.step,
    )
