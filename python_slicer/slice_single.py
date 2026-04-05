#!/usr/bin/env python3
"""
Single-frame CSA slicing with loop filtering for bifurcation geometries.

Slices a full airway STL along a centerline, selecting only the cross-section
loop that encloses the centerline point at each plane. This handles the case
where the plane cuts through multiple branches (nostrils, mouth, etc.).

Usage:
    python python_slicer/slice_single.py surface.stl centerline.vtk
    python python_slicer/slice_single.py surface.stl centerline.vtp -o ./results --name MyCase
"""

import sys
import argparse
import time
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path


def setup_logging(output_dir, name):
    """Setup logging to both console and file."""
    log_path = Path(output_dir) / f"{name}.log"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("slice_single")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    # File handler — detailed with timestamps
    fh = logging.FileHandler(str(log_path), mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s',
                                       datefmt='%H:%M:%S'))

    # Console handler — concise
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(message)s'))

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

sys.path.insert(0, str(Path(__file__).parent))
from slicer.io_utils import read_stl, read_vtk_centerline, write_cross_section_stl
from slicer.mesh_intersection import MeshPlaneSlicer
from slicer.geometry import compute_all_plane_normals, compute_arc_length, convert_units_if_needed
from slicer.loop_filter import filter_loops_by_surface_proximity, merge_loops


def slice_with_loop_filter(stl_path, centerline_path, output_dir, name,
                           window=20, branch_mesh_paths=None):
    """
    Slice mesh along centerline. At each plane, filter loops by branch
    convex hull — keep all loops inside the hull, discard loops from
    other branches, then merge kept loops.
    """
    import trimesh
    log = logging.getLogger("slice_single")

    # Load mesh
    t0 = time.time()
    mesh = read_stl(stl_path)
    log.debug(f"Loaded mesh: {time.time()-t0:.1f}s")

    # Load centerline
    t0 = time.time()
    centerline = read_vtk_centerline(centerline_path)
    n_points = len(centerline)
    log.debug(f"Loaded centerline ({n_points} pts): {time.time()-t0:.1f}s")

    # Load branch surface mesh(es) for proximity filtering
    branch_surface = None
    if branch_mesh_paths:
        t0 = time.time()
        combined_meshes = []
        for bm_path in branch_mesh_paths:
            bm = read_stl(bm_path)
            combined_meshes.append(bm)
            log.info(f"  Branch mesh loaded: {Path(bm_path).name}")
        branch_surface = trimesh.util.concatenate(combined_meshes)
        log.debug(f"Combined branch surface: {time.time()-t0:.1f}s")
        log.info(f"  Combined surface: {len(branch_surface.faces)} faces")

    # Compute normals
    t0 = time.time()
    log.info("Computing plane normals...")
    normals = compute_all_plane_normals(centerline, smooth=True)
    log.debug(f"Computed normals: {time.time()-t0:.1f}s")

    # Prepare slicer
    mesh_slicer = MeshPlaneSlicer(mesh)

    # Slice each plane
    valid_indices = []
    valid_sections = []
    valid_positions = []
    skipped = []

    total_start = time.time()
    log.info(f"\nSlicing {n_points} planes...")
    for i in range(n_points):
        position = centerline[i]
        normal = normals[i]
        plane_start = time.time()

        # Get individual loops
        t0 = time.time()
        loops = mesh_slicer.slice_mesh_with_plane(
            plane_origin=position,
            plane_normal=normal,
            plane_number=i
        )
        slice_time = time.time() - t0

        if not loops:
            log.info(f"  Plane {i+1}/{n_points} (index {i})... SKIP (no intersection)")
            log.debug(f"  Plane {i}: slice={slice_time:.2f}s, no intersection")
            skipped.append((i, "no_intersection"))
            continue

        # Filter loops by branch surface proximity, then merge kept loops
        t0 = time.time()
        if branch_surface is not None and len(loops) > 1:
            kept = filter_loops_by_surface_proximity(loops, branch_surface)
        else:
            kept = loops
        filter_time = time.time() - t0

        selected = merge_loops(kept)

        if selected is None:
            log.info(f"  Plane {i+1}/{n_points} (index {i})... SKIP (no valid loop)")
            log.debug(f"  Plane {i}: slice={slice_time:.2f}s, filter={filter_time:.2f}s, "
                      f"{len(loops)} loops, 0 kept")
            skipped.append((i, "no_valid_loop"))
            continue

        valid_indices.append(i)
        valid_sections.append(selected)
        valid_positions.append(position)

        n_total = len(loops)
        n_kept = len(kept)
        plane_time = time.time() - plane_start

        if n_total > 1:
            log.info(f"  Plane {i+1}/{n_points} (index {i})... "
                     f"OK (area={selected.area:.2f} mm², kept {n_kept}/{n_total} loops)")
        else:
            log.info(f"  Plane {i+1}/{n_points} (index {i})... "
                     f"OK (area={selected.area:.2f} mm²)")
        log.debug(f"  Plane {i}: slice={slice_time:.2f}s, filter={filter_time:.2f}s, "
                  f"total={plane_time:.2f}s, {n_total} loops, {n_kept} kept, "
                  f"area={selected.area:.2f}")

    total_time = time.time() - total_start
    avg_time = total_time / max(n_points, 1)

    # Compute arc lengths
    if valid_positions:
        positions_array = np.array(valid_positions)
        arc_lengths = compute_arc_length(positions_array)

        areas = np.array([s.area for s in valid_sections])
        arc_lengths, areas = convert_units_if_needed(arc_lengths, areas)

        for section, new_area in zip(valid_sections, areas):
            try:
                section.area = new_area
            except AttributeError:
                pass
    else:
        arc_lengths = np.array([])

    log.info(f"\nSlicing complete ({total_time:.1f}s, {avg_time:.2f}s/plane):")
    log.info(f"  Valid planes: {len(valid_sections)}")
    log.info(f"  Skipped: {len(skipped)}")
    if skipped:
        reasons = {}
        for _, reason in skipped:
            reasons[reason] = reasons.get(reason, 0) + 1
        for reason, count in reasons.items():
            log.info(f"    - {reason}: {count}")

    return valid_indices, valid_sections, valid_positions, arc_lengths, skipped


def export_results(valid_indices, valid_sections, arc_lengths, output_dir, name):
    """Export cross-section STLs and measurements CSV."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Individual plane STLs
    all_vertices = []
    all_faces = []
    for k, section in enumerate(valid_sections):
        stl_path = output_dir / f"{name}-Planes-{k:03d}.stl"
        write_cross_section_stl(str(stl_path), section.vertices, section.faces)
        all_vertices.append(section.vertices)
        all_faces.append(section.faces)

    # Combined all-planes STL
    if all_vertices:
        combined_verts = []
        combined_faces = []
        offset = 0
        for verts, faces in zip(all_vertices, all_faces):
            combined_verts.append(verts)
            combined_faces.append(faces + offset)
            offset += len(verts)
        from slicer.io_utils import write_all_planes_stl
        write_all_planes_stl(
            str(output_dir / f"{name}-Planes-All.stl"),
            all_vertices, all_faces)

    # Measurements CSV
    data = {
        "plane_index": valid_indices,
        "arc_length_mm": arc_lengths,
        "area_mm2": [s.area for s in valid_sections],
        "perimeter_mm": [s.perimeter for s in valid_sections],
        "centroid_x": [s.centroid[0] for s in valid_sections],
        "centroid_y": [s.centroid[1] for s in valid_sections],
        "centroid_z": [s.centroid[2] for s in valid_sections],
        "is_valid": [True] * len(valid_sections),
    }

    # Add diameter measurements where possible
    from slicer.measurements import DiameterCalculator
    hyd_diams = []
    equiv_diams = []
    for s in valid_sections:
        hyd_diams.append(DiameterCalculator.compute_hydraulic_diameter(s.area, s.perimeter))
        equiv_diams.append(DiameterCalculator.compute_equivalent_diameter(s.area))
    data["hydraulic_diameter_mm"] = hyd_diams
    data["equivalent_diameter_mm"] = equiv_diams

    df = pd.DataFrame(data)
    csv_path = output_dir / f"{name}-Data.csv"
    df.to_csv(csv_path, index=False)
    print(f"Wrote CSV: {csv_path.name} ({len(df)} rows)")

    return csv_path


def generate_csa_plot(csv_path, output_path, name):
    df = pd.read_csv(csv_path)
    df = df[df["is_valid"] == True] if "is_valid" in df.columns else df

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    axes[0].plot(df["arc_length_mm"], df["area_mm2"], "b-o", markersize=2, linewidth=1)
    axes[0].set_ylabel("Cross-Sectional Area (mm²)")
    axes[0].set_title(f"{name} — Cross-Sectional Area")
    axes[0].grid(True, alpha=0.3)

    if "hydraulic_diameter_mm" in df.columns:
        axes[1].plot(df["arc_length_mm"], df["hydraulic_diameter_mm"],
                     "r-o", markersize=2, linewidth=1)
    axes[1].set_xlabel("Arc Length (mm)")
    axes[1].set_ylabel("Hydraulic Diameter (mm)")
    axes[1].set_title(f"{name} — Hydraulic Diameter")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote plot: {Path(output_path).name}")


def generate_summary(csv_path, output_path, name, n_total, n_valid, n_skipped):
    df = pd.read_csv(csv_path)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"CSA Slicing Summary: {name}\n{'='*50}\n\n")
        f.write(f"Centerline points: {n_total}\n")
        f.write(f"Valid planes: {n_valid}\n")
        f.write(f"Skipped planes: {n_skipped}\n")
        f.write(f"Success rate: {n_valid/max(n_total,1)*100:.1f}%\n\n")
        if len(df) > 0:
            f.write(f"Arc length: {df['arc_length_mm'].min():.1f} - {df['arc_length_mm'].max():.1f} mm\n\n")
            f.write(f"Area (mm²): min={df['area_mm2'].min():.2f}, max={df['area_mm2'].max():.2f}, mean={df['area_mm2'].mean():.2f}\n")
            if "hydraulic_diameter_mm" in df.columns:
                f.write(f"Hyd. Diameter (mm): min={df['hydraulic_diameter_mm'].min():.2f}, max={df['hydraulic_diameter_mm'].max():.2f}, mean={df['hydraulic_diameter_mm'].mean():.2f}\n")
    print(f"Wrote summary: {Path(output_path).name}")


def main():
    parser = argparse.ArgumentParser(
        description="Single-frame CSA slicing with loop filtering",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python python_slicer/slice_single.py surface.stl centerline.vtk
  python python_slicer/slice_single.py surface.stl centerline.vtp -o ./results --name MyCase
        """)
    parser.add_argument("stl_path", help="Path to surface STL file")
    parser.add_argument("centerline_path", help="Path to centerline VTK or VTP file")
    parser.add_argument("-o", "--output", default="./slice_output")
    parser.add_argument("--name", default=None,
                        help="Name prefix (default: from STL filename)")
    parser.add_argument("--branch-mesh", nargs='+', default=None,
                        help="Path(s) to branch surface STL(s) for convex hull filtering. "
                             "Multiple meshes are combined into one hull. "
                             "Example: --branch-mesh LeftNose.stl DescendingAirway.stl")
    parser.add_argument("--window", type=int, default=20,
                        help="Centerline smoothing window (default: 20)")
    args = parser.parse_args()

    stl_path = Path(args.stl_path).resolve()
    cl_path = Path(args.centerline_path).resolve()
    output_dir = Path(args.output).resolve()

    branch_mesh_paths = None
    if args.branch_mesh:
        branch_mesh_paths = [str(Path(p).resolve()) for p in args.branch_mesh]

    if not stl_path.exists():
        print(f"ERROR: STL not found: {stl_path}"); sys.exit(1)
    if not cl_path.exists():
        print(f"ERROR: Centerline not found: {cl_path}"); sys.exit(1)
    if branch_mesh_paths:
        for bmp in branch_mesh_paths:
            if not Path(bmp).exists():
                print(f"ERROR: Branch mesh not found: {bmp}"); sys.exit(1)

    name = args.name or stl_path.stem
    output_dir.mkdir(parents=True, exist_ok=True)

    log = setup_logging(str(output_dir), name)

    log.info("=" * 60)
    log.info(f"SINGLE-FRAME CSA SLICING")
    log.info(f"  Surface: {stl_path.name}")
    log.info(f"  Centerline: {cl_path.name}")
    if branch_mesh_paths:
        for bmp in branch_mesh_paths:
            log.info(f"  Branch mesh: {Path(bmp).name}")
    log.info(f"  Output: {output_dir}")
    log.info(f"  Name: {name}")
    log.info(f"  Log: {output_dir / f'{name}.log'}")
    log.info("=" * 60)

    run_start = time.time()

    # Slice with loop filtering
    valid_indices, valid_sections, valid_positions, arc_lengths, skipped = \
        slice_with_loop_filter(
            str(stl_path), str(cl_path), str(output_dir), name,
            window=args.window,
            branch_mesh_paths=branch_mesh_paths)

    n_total = len(valid_sections) + len(skipped)
    n_valid = len(valid_sections)
    n_skipped = len(skipped)

    if n_valid == 0:
        log.error("No valid cross-sections found.")
        sys.exit(1)

    # Export
    t0 = time.time()
    log.info(f"\n--- Exporting ---")
    csv_path = export_results(valid_indices, valid_sections, arc_lengths,
                              output_dir, name)
    log.debug(f"Export: {time.time()-t0:.1f}s")

    # Plot + summary
    t0 = time.time()
    generate_csa_plot(csv_path, output_dir / f"{name}-CSA_plot.png", name)
    generate_summary(csv_path, output_dir / f"{name}-summary.txt",
                     name, n_total, n_valid, n_skipped)
    log.debug(f"Plot + summary: {time.time()-t0:.1f}s")

    total_time = time.time() - run_start
    log.info(f"\n{'='*60}")
    log.info(f"DONE in {total_time:.1f}s. Results in: {output_dir}")
    log.info(f"  Planes: {n_valid} valid, {n_skipped} skipped")
    log.info(f"  CSV:  {name}-Data.csv")
    log.info(f"  Plot: {name}-CSA_plot.png")
    log.info(f"  Log:  {name}.log")
    log.info(f"{'='*60}")


if __name__ == "__main__":
    main()
