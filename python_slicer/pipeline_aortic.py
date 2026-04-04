#!/usr/bin/env python3
"""
Aortic CSA Pipeline (Single Time Point)

Takes an aortic STL with 1 inlet + N outlets and produces:
- Branch detection + centerline extraction (VMTK)
- Cross-sectional slicing per branch
- CSA measurements, plots, PDF reports

Usage:
    conda run -n mirtk python python_slicer/pipeline_aortic.py AorticSubject/ --subject-id AORT001
    conda run -n mirtk python python_slicer/pipeline_aortic.py AorticSubject/ --subject-id AORT001 --skip-branch-detect
"""

import sys
import argparse
import logging
import time
import shutil
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from slicer.io_utils import read_stl, read_vtk_centerline, write_cross_section_stl
from slicer.geometry import compute_arc_length
from slicer.measurements import DiameterCalculator
from slicer.loop_filter import merge_loops
from slice_bifurcation import (
    slice_along_centerline, slice_along_midline, _save_loops_stl,
    setup_logging,
)

log = logging.getLogger(__name__)


def find_surface_stl(subject_dir):
    """Find the single STL file in surface/ folder."""
    surface_dir = Path(subject_dir) / "surface"
    stl_files = list(surface_dir.glob("*.stl"))
    if not stl_files:
        raise FileNotFoundError(f"No STL found in {surface_dir}")
    return stl_files[0]


def discover_branches_aortic(branches_dir, subject_id=""):
    """Auto-discover branches from aortic branch detector output.

    Identifies the main trunk (longest full-path CL) vs side branches.
    Side branches use segment CLs (branch-only portion).
    Trunk uses full-path CL (inlet→largest outlet).
    """
    prefix = f"{subject_id}_" if subject_id else ""
    branches = {}

    for d in sorted(branches_dir.iterdir()):
        if not d.is_dir():
            continue
        label = d.name

        # Find full-path centerline
        cl_files = list(d.glob(f"*_to_{label}_centerline.vtk"))
        cl_path = cl_files[0] if cl_files else None

        # Find segment centerline
        seg_files = list(d.glob(f"*{label}_segment_centerline.vtk"))
        seg_path = seg_files[0] if seg_files else None

        # Find mesh
        mesh_files = list(d.glob(f"*{label}_mesh.stl"))
        mesh_path = mesh_files[0] if mesh_files else None

        branches[label] = {
            "dir": d,
            "centerline": cl_path,
            "segment_centerline": seg_path,
            "mesh": mesh_path,
        }

    # Identify main trunk = branch with longest full-path centerline
    # Side branches use segment CLs, trunk uses full-path CL
    longest_label = None
    longest_len = 0
    for label, info in branches.items():
        if info["centerline"]:
            from slicer.io_utils import read_vtk_centerline
            cl = read_vtk_centerline(str(info["centerline"]))
            if len(cl) > longest_len:
                longest_len = len(cl)
                longest_label = label

    # Load main CL for ordering branches by divergence point
    main_info = branches[longest_label]
    main_cl = read_vtk_centerline(str(main_info["centerline"]))

    # Collect side branches with their divergence index along the main CL
    side_branches = []
    for label, info in branches.items():
        if label == longest_label or label == "Trunk":
            continue
        # Find where this branch's full-path CL diverges from the main CL
        if info["centerline"]:
            branch_cl = read_vtk_centerline(str(info["centerline"]))
            # Walk along both CLs, find last shared point
            from scipy.spatial import cKDTree
            main_tree = cKDTree(main_cl)
            dists, _ = main_tree.query(branch_cl)
            # Divergence = first point where distance exceeds threshold
            diverge_idx = 0
            for j in range(len(dists)):
                if dists[j] < 2.0:  # within 2mm = still on main CL
                    diverge_idx = j
                else:
                    break
            side_branches.append((diverge_idx, label, info))

    # Sort by divergence index (appearance order along main CL)
    side_branches.sort(key=lambda x: x[0])

    # Build result: MainAorta first, then branches in anatomical order
    result = {}
    result["MainAorta"] = {
        "dir": main_info["dir"],
        "centerline": main_info["centerline"],
        "segment_centerline": main_info["segment_centerline"],
        "mesh": main_info["mesh"],
        "is_trunk": True,
    }

    for i, (div_idx, old_label, info) in enumerate(side_branches):
        new_label = f"Branch_{i+1}"
        best_cl = info["segment_centerline"] or info["centerline"]
        result[new_label] = {
            "dir": info["dir"],
            "centerline": best_cl,
            "segment_centerline": info["segment_centerline"],
            "mesh": info["mesh"],
            "is_trunk": False,
            "_original_label": old_label,
        }

    return result


def _extend_cl_to_profile(centerline, mesh, offset_mm=1.0):
    """Extend centerline start/end to 1mm inside the open profile boundary.

    Detects open profile edges near each CL endpoint, computes the
    profile barycenter, and prepends/appends a point 1mm inward.
    """
    import trimesh

    extended = centerline.copy()

    for end in ["start", "end"]:
        if end == "start":
            cl_tip = centerline[0]
            cl_next = centerline[1]
        else:
            cl_tip = centerline[-1]
            cl_next = centerline[-2]

        # Direction from tip into the lumen
        inward = cl_next - cl_tip
        inward = inward / (np.linalg.norm(inward) + 1e-12)

        # Find boundary edges near the tip
        edges = mesh.edges_unique
        edge_midpoints = mesh.vertices[edges].mean(axis=1)
        dists = np.linalg.norm(edge_midpoints - cl_tip, axis=1)

        # Boundary edges: edges that appear in only 1 face
        from collections import Counter
        edge_face_count = Counter()
        for face in mesh.faces:
            for i in range(3):
                e = tuple(sorted([face[i], face[(i+1) % 3]]))
                edge_face_count[e] += 1

        boundary_mask = np.array([
            edge_face_count.get(tuple(sorted(e)), 0) == 1
            for e in edges
        ])

        # Boundary edges near the tip (within 2x the local radius)
        nearby = dists < np.percentile(dists, 10)
        profile_edges = boundary_mask & nearby

        if profile_edges.sum() < 3:
            continue

        # Profile barycenter
        profile_pts = edge_midpoints[profile_edges]
        barycenter = profile_pts.mean(axis=0)

        # New point: barycenter (right at the cap opening)
        new_pt = barycenter

        if end == "start":
            extended = np.vstack([new_pt.reshape(1, 3), extended])
        else:
            extended = np.vstack([extended, new_pt.reshape(1, 3)])

    return extended


def _slice_with_containment(mesh, centerline, label, output_dir, name, log):
    """Slice mesh along CL, keep only loops containing the CL point.

    Uses trimesh.section + shapely containment check. This correctly
    discards cross-lumen cuts where the CL point is outside the loop.
    """
    import trimesh
    from shapely.geometry import Point
    from slicer.geometry import compute_all_plane_normals
    from trimesh.creation import triangulate_polygon

    output_dir = Path(output_dir)
    planes_dir = output_dir / f"{label}_Planes"
    planes_dir.mkdir(parents=True, exist_ok=True)

    normals = compute_all_plane_normals(centerline, smooth=True)
    all_rows = []
    positions = []
    count = 0

    for i in range(len(centerline)):
        section = mesh.section(plane_origin=centerline[i], plane_normal=normals[i])
        if section is None:
            continue

        try:
            planar, transform = section.to_planar()
        except Exception:
            continue

        # Project CL point to 2D
        inv_transform = np.linalg.inv(transform)
        cl_hom = np.append(centerline[i], 1.0)
        cl_2d = (inv_transform @ cl_hom)[:2]
        cl_point = Point(cl_2d[0], cl_2d[1])

        # Check each polygon — keep only if CL point is inside
        for poly in planar.polygons_full:
            if poly.is_empty or poly.area < 1e-6:
                continue
            if not poly.contains(cl_point):
                continue

            # Triangulate and save
            try:
                tri_v2d, tri_f = triangulate_polygon(poly)
            except Exception:
                continue

            tri_v3d = trimesh.transformations.transform_points(
                np.column_stack([tri_v2d, np.zeros(len(tri_v2d))]), transform)
            centroid = tri_v3d.mean(axis=0)

            stl_path = planes_dir / f"{name}-{label}-{count:03d}.stl"
            write_cross_section_stl(str(stl_path), tri_v3d, tri_f)

            area = poly.area
            perimeter = poly.length
            hyd_diam = DiameterCalculator.compute_hydraulic_diameter(area, perimeter)
            equiv_diam = DiameterCalculator.compute_equivalent_diameter(area)

            positions.append(centerline[i])
            all_rows.append({
                "plane_index": i, "region": label,
                "arc_length_mm": 0.0,  # filled after
                "area_mm2": area, "perimeter_mm": perimeter,
                "hydraulic_diameter_mm": hyd_diam, "equivalent_diameter_mm": equiv_diam,
                "centroid_x": centroid[0], "centroid_y": centroid[1],
                "centroid_z": centroid[2], "is_valid": True,
            })
            count += 1
            break  # one loop per plane

        if (i + 1) % 20 == 0 or i == len(centerline) - 1:
            log.info(f"    {label} Plane {i+1}/{len(centerline)}... {count} kept")

    # Fill arc lengths from positions
    if positions:
        arc = compute_arc_length(np.array(positions))
        for k, row in enumerate(all_rows):
            row["arc_length_mm"] = float(arc[k]) if k < len(arc) else 0.0

    log.info(f"    {label}: {count} planes")
    return all_rows


def slice_branch(mesh, centerline, branch_mesh, log, label, output_dir, name):
    """Slice a single branch and return measurement rows.

    Uses slice_along_midline on the branch mesh (keeps all loops).
    Falls back to slice_along_centerline on full mesh if branch mesh unavailable.
    """
    output_dir = Path(output_dir)
    planes_dir = output_dir / f"{label}_Planes"
    planes_dir.mkdir(parents=True, exist_ok=True)

    if branch_mesh is not None:
        # Slice the BRANCH MESH directly — only produces loops on this surface
        # No cross-contamination from adjacent branches
        results = slice_along_midline(branch_mesh, centerline, log,
                                       label_prefix=f"{label} ")
    else:
        # Fallback: slice full mesh with containment filter
        results = slice_along_centerline(mesh, centerline, None, log,
                                          label_prefix=f"{label} ")

    if not results:
        log.warning(f"  {label}: no valid planes")
        return []

    # Compute arc lengths
    if branch_mesh is not None:
        # Results from slice_along_midline: (plane_idx, pos, normal, loop_data)
        # Arc length from actual result positions (not raw CL indices)
        result_positions = np.array([pos for _, pos, _, _ in results])
        result_arc = compute_arc_length(result_positions) if len(result_positions) > 0 else np.array([])
        rows = []
        count = 0
        for k, (plane_idx, pos, normal, loop_data) in enumerate(results):
            arc_len = float(result_arc[k]) if k < len(result_arc) else 0.0
            total_area = sum(L["area"] for L in loop_data)
            total_perimeter = sum(L["perimeter"] for L in loop_data)
            if len(loop_data) == 1:
                mc = loop_data[0]["centroid"]
            else:
                cs = np.array([L["centroid"] for L in loop_data])
                ar = np.array([L["area"] for L in loop_data])
                mc = np.average(cs, axis=0, weights=ar)

            stl_path = planes_dir / f"{name}-{label}-{count:03d}.stl"
            _save_loops_stl(loop_data, str(stl_path))
            count += 1

            hyd_diam = DiameterCalculator.compute_hydraulic_diameter(total_area, total_perimeter)
            equiv_diam = DiameterCalculator.compute_equivalent_diameter(total_area)
            rows.append({
                "plane_index": plane_idx,
                "region": label,
                "arc_length_mm": arc_len,
                "area_mm2": total_area,
                "perimeter_mm": total_perimeter,
                "hydraulic_diameter_mm": hyd_diam,
                "equivalent_diameter_mm": equiv_diam,
                "centroid_x": mc[0], "centroid_y": mc[1], "centroid_z": mc[2],
                "is_valid": True,
            })
    else:
        # Results from slice_along_centerline: (plane_idx, pos, loops)
        positions = np.array([pos for _, pos, _ in results])
        arc = compute_arc_length(positions) if len(positions) > 0 else np.array([])
        rows = []
        count = 0
        for k, (plane_idx, pos, loops) in enumerate(results):
            arc_len = float(arc[k]) if k < len(arc) else 0.0
            merged = merge_loops(loops)
            if merged is None:
                continue
            stl_path = planes_dir / f"{name}-{label}-{count:03d}.stl"
            write_cross_section_stl(str(stl_path), merged.vertices, merged.faces)
            count += 1

            hyd_diam = DiameterCalculator.compute_hydraulic_diameter(merged.area, merged.perimeter)
            equiv_diam = DiameterCalculator.compute_equivalent_diameter(merged.area)
            rows.append({
                "plane_index": plane_idx,
                "region": label,
                "arc_length_mm": arc_len,
                "area_mm2": merged.area,
                "perimeter_mm": merged.perimeter,
                "hydraulic_diameter_mm": hyd_diam,
                "equivalent_diameter_mm": equiv_diam,
                "centroid_x": merged.centroid[0],
                "centroid_y": merged.centroid[1],
                "centroid_z": merged.centroid[2],
                "is_valid": True,
            })

    log.info(f"  {label}: {count} planes")
    return rows


def generate_aortic_csa_plot(csv_path, output_path, name):
    """CSA plot: continuous main aorta (ascending+descending) + side branches.

    Clinical focus: CoA patients — show the narrowing around the arch.
    Main aorta plotted as one continuous line with ascending and descending
    joined by arc length. Arch branches shown as markers.
    """
    df = pd.read_csv(csv_path)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 2, width_ratios=[3, 1], hspace=0.3, wspace=0.3)

    # ── Left column: Main aorta continuous CSA ──
    # MainAorta = already merged in CSV (continuous arc length)
    main = df[df["region"] == "MainAorta"].sort_values("arc_length_mm").copy()

    # Top-left: CSA along main aorta
    ax1 = fig.add_subplot(gs[0, 0])
    if len(main) > 0:
        ax1.plot(main["arc_length_mm"], main["area_mm2"],
                 "-o", color="steelblue", markersize=3, linewidth=2,
                 label="MainAorta", alpha=0.9)

        # Mark minimum CSA (coarctation site)
        min_idx = main["area_mm2"].idxmin()
        min_row = main.loc[min_idx]
        ax1.plot(min_row["arc_length_mm"], min_row["area_mm2"],
                 "v", color="red", markersize=12, zorder=5,
                 label=f"Min CSA: {min_row['area_mm2']:.1f} mm²")

        # Mark max CSA
        max_idx = main["area_mm2"].idxmax()
        max_row = main.loc[max_idx]
        ax1.plot(max_row["arc_length_mm"], max_row["area_mm2"],
                 "^", color="green", markersize=10, zorder=5,
                 label=f"Max CSA: {max_row['area_mm2']:.1f} mm²")

        # Ratio annotation
        ratio = max_row["area_mm2"] / min_row["area_mm2"] if min_row["area_mm2"] > 0 else 0
        ax1.annotate(f"Max/Min ratio: {ratio:.1f}x",
                    xy=(0.5, 0.05), xycoords='axes fraction',
                    fontsize=11, fontweight='bold', color='darkred',
                    ha='center', bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

    ax1.set_ylabel("Cross-Sectional Area (mm²)", fontsize=11)
    ax1.set_title(f"{name} — Main Aorta CSA (Ascending → Descending)",
                  fontsize=12, fontweight='bold')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Bottom-left: Hydraulic diameter
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    if len(main) > 0:
        ax2.plot(main["arc_length_mm"], main["hydraulic_diameter_mm"],
                 "-o", color="steelblue", markersize=3, linewidth=2, alpha=0.9)
    ax2.set_xlabel("Arc Length along Main Aorta (mm)", fontsize=11)
    ax2.set_ylabel("Hydraulic Diameter (mm)", fontsize=11)
    ax2.set_title(f"{name} — Hydraulic Diameter", fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # ── Right column: Branch CSA ──
    branch_regions = [r for r in df["region"].unique()
                      if r != "MainAorta"]
    branch_colors = plt.cm.Set2(np.linspace(0, 1, max(len(branch_regions), 1)))

    ax3 = fig.add_subplot(gs[0, 1])
    for i, region in enumerate(sorted(branch_regions)):
        grp = df[df["region"] == region].sort_values("arc_length_mm")
        ax3.plot(grp["arc_length_mm"], grp["area_mm2"],
                 "-o", color=branch_colors[i], markersize=4, linewidth=1.5,
                 label=region, alpha=0.8)
    ax3.set_ylabel("CSA (mm²)", fontsize=10)
    ax3.set_title("Arch Branches", fontsize=11, fontweight='bold')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    # Bottom-right: summary stats table
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    stats_lines = [f"{'Region':<14} {'Planes':>6} {'Min CSA':>8} {'Max CSA':>8}"]
    stats_lines.append("-" * 40)
    for region in ["MainAorta"] + sorted(branch_regions):
        rdf = df[df["region"] == region]
        if len(rdf) == 0:
            continue
        stats_lines.append(
            f"{region:<14} {len(rdf):>6} {rdf['area_mm2'].min():>8.1f} "
            f"{rdf['area_mm2'].max():>8.1f}")
    ax4.text(0.05, 0.95, "\n".join(stats_lines), transform=ax4.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace')
    ax4.set_title("Summary", fontsize=11, fontweight='bold')

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Aortic CSA pipeline (single time point)",
    )
    parser.add_argument("subject_dir", help="Subject directory with surface/")
    parser.add_argument("--subject-id", default="")
    parser.add_argument("--skip-branch-detect", action="store_true",
                        help="Skip if branches/ already exists")
    parser.add_argument("--resampling-step", type=float, default=2.0)
    args = parser.parse_args()

    subject_dir = Path(args.subject_dir).resolve()
    subject_id = args.subject_id or subject_dir.name
    output_dir = subject_dir / "csa"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging(str(output_dir), subject_id)
    pipeline_start = time.time()

    logger.info("=" * 60)
    logger.info("AORTIC CSA PIPELINE (Single Time Point)")
    logger.info(f"  Subject: {subject_dir}")
    logger.info(f"  Subject ID: {subject_id}")
    logger.info("=" * 60)

    # ── Step 1: Branch detection ────────────────────────────
    branches_dir = subject_dir / "branches"

    if args.skip_branch_detect and branches_dir.exists():
        logger.info(f"\n--- Step 1: Skipping branch detection (exists) ---")
    else:
        logger.info(f"\n--- Step 1: Branch detection (vmtk env) ---")
        import subprocess
        script_dir = Path(__file__).parent
        stl_path = find_surface_stl(subject_dir)

        # Run aortic branch detector in vmtk env
        cmd = [
            "conda", "run", "-n", "vmtk",
            "python", str(script_dir / "aortic_branch_detector.py"),
            str(stl_path),
            "--output-dir", str(branches_dir),
            "--subject-id", subject_id,
            "--resampling-step", str(args.resampling_step),
        ]
        result = subprocess.run(cmd, capture_output=False)
        if result.returncode != 0:
            logger.error("Branch detection failed!")
            sys.exit(1)

    # Discover branches
    branches = discover_branches_aortic(branches_dir, subject_id)
    logger.info(f"  Found {len(branches)} branches: {list(branches.keys())}")

    # ── Step 2: Slice each branch ───────────────────────────
    logger.info(f"\n--- Step 2: Slicing ---")

    full_mesh = read_stl(str(find_surface_stl(subject_dir)))
    all_rows = []

    for label, info in branches.items():
        cl_path = info.get("centerline")
        if not cl_path:
            logger.warning(f"  {label}: no centerline found, skipping")
            continue

        centerline = read_vtk_centerline(str(cl_path))

        # Resample to uniform 2mm spacing
        from scipy.interpolate import interp1d
        arc = compute_arc_length(centerline)
        total_len = arc[-1] if len(arc) > 0 else 0
        if total_len > 0:
            n_pts = max(int(total_len / 2.0), 5)
            t_uniform = np.linspace(0, total_len, n_pts)
            centerline = interp1d(arc, centerline, axis=0, kind='cubic')(t_uniform)

        # Extend CL to open profile boundaries (1mm inside)
        centerline = _extend_cl_to_profile(centerline, full_mesh, offset_mm=1.0)

        logger.info(f"  {label}: CL {len(centerline)} pts "
                     f"({total_len:.1f}mm, resampled to 2mm, extended to profiles)")

        # Load branch mesh if available
        branch_mesh = None
        if info.get("mesh") and Path(info["mesh"]).exists():
            branch_mesh = read_stl(str(info["mesh"]))

        # For MainAorta: slice Ascending and Descending separately
        # Both use same full CL, slice their own mesh, filter by containment
        # (CL point must be inside the loop — discards cross-lumen cuts)
        if info.get("is_trunk"):
            prefix_str = f"{subject_id}_" if subject_id else ""
            trunk_mesh_path = branches_dir / "Trunk" / f"{prefix_str}Trunk_mesh.stl"
            trunk_mesh = read_stl(str(trunk_mesh_path)) if trunk_mesh_path.exists() else None
            desc_mesh = branch_mesh  # Branch_2 mesh = descending

            # Slice Ascending
            if trunk_mesh is not None:
                logger.info(f"    Slicing Ascending (Trunk: {len(trunk_mesh.faces)} faces)")
                rows_asc = _slice_with_containment(
                    trunk_mesh, centerline, "Ascending", output_dir, subject_id, logger)
                all_rows.extend(rows_asc)

            # Slice Descending
            if desc_mesh is not None:
                logger.info(f"    Slicing Descending (Branch: {len(desc_mesh.faces)} faces)")
                rows_desc = _slice_with_containment(
                    desc_mesh, centerline, "Descending", output_dir, subject_id, logger)
                all_rows.extend(rows_desc)
        else:
            rows = slice_branch(
                full_mesh, centerline, branch_mesh,
                logger, label, str(output_dir), subject_id,
            )
            all_rows.extend(rows)

    logger.info(f"\n  Total rows: {len(all_rows)}")

    # Save all planes combined STL per region
    import trimesh
    for region_name in set(r["region"] for r in all_rows):
        planes_dir = output_dir / f"{region_name}_Planes"
        if planes_dir.exists():
            stl_files = sorted(planes_dir.glob("*.stl"))
            if stl_files:
                meshes = [trimesh.load_mesh(str(f)) for f in stl_files]
                combined = trimesh.util.concatenate(meshes)
                combined.export(str(output_dir / f"{region_name}_All_Planes.stl"))
                logger.info(f"  {region_name}_All_Planes.stl ({len(stl_files)} planes)")

    # Create MainAorta_Planes folder (combine Ascending + Descending STLs)
    main_planes_dir = output_dir / "MainAorta_Planes"
    main_planes_dir.mkdir(parents=True, exist_ok=True)
    for src_folder in ["Ascending_Planes", "Descending_Planes"]:
        src_dir = output_dir / src_folder
        if src_dir.exists():
            for stl_file in sorted(src_dir.glob("*.stl")):
                shutil.copy2(str(stl_file), str(main_planes_dir / stl_file.name))
    logger.info(f"  MainAorta_Planes: {len(list(main_planes_dir.glob('*.stl')))} planes")

    # Combine Ascending + Descending into MainAorta_All_Planes
    asc_stl = output_dir / "Ascending_All_Planes.stl"
    desc_stl = output_dir / "Descending_All_Planes.stl"
    if asc_stl.exists() and desc_stl.exists():
        main_combined = trimesh.util.concatenate([
            trimesh.load_mesh(str(asc_stl)),
            trimesh.load_mesh(str(desc_stl)),
        ])
        main_combined.export(str(output_dir / "MainAorta_All_Planes.stl"))
        logger.info(f"  MainAorta_All_Planes.stl (Ascending + Descending)")

    if not all_rows:
        logger.error("No valid cross-sections found!")
        sys.exit(1)

    # ── Step 3: Save CSV ────────────────────────────────────
    logger.info(f"\n--- Step 3: Saving CSV ---")
    df = pd.DataFrame(all_rows)

    # Save raw CSV (internal regions: Ascending, Descending, Branch_X)
    raw_csv = output_dir / f"{subject_id}-Data-raw.csv"
    df.to_csv(raw_csv, index=False)

    # Merge Ascending + Descending → MainAorta with continuous arc length
    # Renumber arch branches sequentially (Branch_1, Branch_2, Branch_3)
    df_out = df.copy()

    # MainAorta: combine Ascending + Descending
    asc_mask = df_out["region"] == "Ascending"
    desc_mask = df_out["region"] == "Descending"
    if asc_mask.any() and desc_mask.any():
        asc_max_arc = df_out.loc[asc_mask, "arc_length_mm"].max()
        desc_min_arc = df_out.loc[desc_mask, "arc_length_mm"].min()
        df_out.loc[desc_mask, "arc_length_mm"] = (
            df_out.loc[desc_mask, "arc_length_mm"] - desc_min_arc + asc_max_arc + 2.0
        )
    df_out.loc[asc_mask | desc_mask, "region"] = "MainAorta"

    # Branches are already named Branch_1/2/3 in anatomical order from discovery

    csv_path = output_dir / f"{subject_id}-Data.csv"
    df_out.to_csv(csv_path, index=False)
    logger.info(f"  {csv_path.name}: {len(df_out)} rows")

    for region in df_out["region"].unique():
        n = len(df_out[df_out["region"] == region])
        logger.info(f"    {region}: {n} planes")

    # ── Step 4: Plots ───────────────────────────────────────
    logger.info(f"\n--- Step 4: Visualization ---")
    plot_path = output_dir / f"{subject_id}-CSA_plot.png"
    generate_aortic_csa_plot(csv_path, plot_path, subject_id)
    logger.info(f"  Plot: {plot_path.name}")

    # ── Step 5: Post-processing ─────────────────────────────
    logger.info(f"\n--- Step 5: Post-processing ---")
    try:
        from postprocess_aortic import run_aortic_postprocessing
        run_aortic_postprocessing(output_dir, subject_id)
    except ImportError:
        logger.info("  postprocess_aortic.py not found, skipping")
    except Exception as e:
        logger.warning(f"  Post-processing failed: {e}")
        logger.info("  CSV data saved successfully — rerun postprocessing manually")

    total_time = time.time() - pipeline_start
    logger.info(f"\n{'='*60}")
    logger.info(f"DONE in {total_time:.1f}s")
    logger.info(f"  Branches: {len(branches)}")
    logger.info(f"  Total planes: {len(all_rows)}")
    logger.info(f"  Output: {output_dir}")
    logger.info(f"{'='*60}")


if __name__ == "__main__":
    main()
