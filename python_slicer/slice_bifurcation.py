#!/usr/bin/env python3
"""
Bifurcation-aware CSA slicing pipeline.

Slices a full airway mesh along anatomically-correct centerlines, producing
cross-section measurements labeled by region (DescendingAirway, LeftNose,
RightNose, Mouth).

Requires branch_detector.py output (VMTK centerlines + surface partitions)
as input. Refines the nasal L/R split at the posterior septum using a
gap-centered midline.

Usage:
    python python_slicer/slice_bifurcation.py surface.stl branches_dir/
    python python_slicer/slice_bifurcation.py surface.stl branches_dir/ --name ENT001 -o ./results
"""

import sys
import argparse
import time
import logging
import numpy as np
import pandas as pd
import trimesh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
from slicer.io_utils import (
    read_stl, read_vtk_centerline, write_cross_section_stl,
)
from slicer.mesh_intersection import MeshPlaneSlicer
from slicer.geometry import compute_all_plane_normals, compute_arc_length, convert_units_if_needed
from slicer.loop_filter import filter_loops_by_surface_proximity, merge_loops
from slicer.measurements import DiameterCalculator
from slicer.septum_refine import refine_nasal_partition


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(output_dir, name):
    log_path = Path(output_dir) / f"{name}_bifurcation.log"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    fh = logging.FileHandler(str(log_path), mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                                       datefmt="%H:%M:%S"))
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(message)s"))

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# ---------------------------------------------------------------------------
# Branch discovery
# ---------------------------------------------------------------------------

def discover_branches(branches_dir: Path, subject_id: str = ""):
    """Auto-discover available branches from branch_detector output."""
    prefix = f"{subject_id}_" if subject_id else ""
    branches = {}

    for label in ["LeftNose", "RightNose", "Mouth", "DescendingAirway"]:
        d = branches_dir / label
        if not d.exists():
            continue

        # Find centerline (full-path: Trachea_to_{label})
        cl_patterns = [
            f"{prefix}Trachea_to_{label}_centerline.vtk",
            f"Trachea_to_{label}_centerline.vtk",
        ]
        cl_path = None
        for pat in cl_patterns:
            p = d / pat
            if p.exists():
                cl_path = p
                break

        # Find segment centerline (bifurcation → outlet)
        seg_patterns = [
            f"{prefix}{label}_segment_centerline.vtk",
            f"{label}_segment_centerline.vtk",
        ]
        seg_path = None
        for pat in seg_patterns:
            p = d / pat
            if p.exists():
                seg_path = p
                break

        # Find mesh
        mesh_patterns = [
            f"{prefix}{label}_mesh.stl",
            f"{label}_mesh.stl",
        ]
        mesh_path = None
        for pat in mesh_patterns:
            p = d / pat
            if p.exists():
                mesh_path = p
                break

        branches[label] = {
            "dir": d,
            "centerline": cl_path,
            "segment_centerline": seg_path,
            "mesh": mesh_path,
        }

    return branches


# ---------------------------------------------------------------------------
# Region classification
# ---------------------------------------------------------------------------

def build_vertex_classifier(refine_result):
    """Build a classifier that labels loops by vertex matching to refined meshes.

    Uses KDTree on vertices (not face centers) — instant and exact since the
    cross-section vertices come from the same mesh.
    """
    left_tree = cKDTree(refine_result["left_mesh"].triangles_center)
    right_tree = cKDTree(refine_result["right_mesh"].triangles_center)
    septum_pos = refine_result["septum_pos"]
    septum_normal = refine_result["septum_normal"]

    def classify_loop(loop):
        """Classify a CrossSection loop by majority vote of its vertices.

        For each vertex in the loop, find nearest face center on left_mesh
        vs right_mesh. Majority vote across ALL vertices determines the region.
        Robust for irregular/T-shaped cross-sections.
        """
        if np.dot(loop.centroid - septum_pos, septum_normal) < 0:
            return "DescendingAirway"

        # Above septum: check ALL vertices against L vs R face centers
        verts = loop.vertices
        if len(verts) == 0:
            return "DescendingAirway"

        d_left, _ = left_tree.query(verts)
        d_right, _ = right_tree.query(verts)

        n_left = (d_left < d_right).sum()
        n_right = (d_right <= d_left).sum()

        if n_left > n_right:
            return "LeftNose"
        else:
            return "RightNose"

    return classify_loop


# ---------------------------------------------------------------------------
# Slicing helpers
# ---------------------------------------------------------------------------

def slice_along_centerline(mesh, centerline, branch_surface, log,
                           label_prefix=""):
    """Slice mesh along centerline WITH containment filtering.

    Uses MeshPlaneSlicer which keeps only loops containing the centerline
    point. Good for simple tube geometries (descending airway).
    Returns list of (plane_idx, position, loops_kept) tuples.
    """
    normals = compute_all_plane_normals(centerline, smooth=True)
    slicer = MeshPlaneSlicer(mesh)
    results = []
    n = len(centerline)

    for i in range(n):
        loops = slicer.slice_mesh_with_plane(
            plane_origin=centerline[i],
            plane_normal=normals[i],
            plane_number=i,
        )
        if not loops:
            continue

        if branch_surface is not None and len(loops) > 1:
            kept = filter_loops_by_surface_proximity(loops, branch_surface)
        else:
            kept = loops

        if not kept:
            continue

        results.append((i, centerline[i], kept))

        if (i + 1) % 20 == 0 or i == n - 1:
            log.info(f"  {label_prefix}Plane {i+1}/{n}... "
                     f"{len(kept)} loop(s)")

    return results


def slice_along_midline(mesh, centerline, log, label_prefix=""):
    """Slice mesh along centerline using raw plane intersection.

    Uses trimesh.section() directly — keeps ALL loops, no containment
    filtering. The midline defines plane position + normal only.

    Returns list of (plane_idx, position, normal, section_paths) tuples.
    Each section_paths is a list of Nx3 arrays (one per loop).
    """
    normals = compute_all_plane_normals(centerline, smooth=True)
    results = []
    n = len(centerline)

    for i in range(n):
        try:
            section = mesh.section(
                plane_origin=centerline[i],
                plane_normal=normals[i],
            )
            if section is None:
                continue

            # Get 2D planar representation + transform for proper triangulation
            try:
                planar, transform = section.to_planar()
            except Exception:
                continue
            polygons = planar.polygons_full
            paths = section.discrete
            if not paths or not polygons:
                continue

            # Build loop data with proper triangulation
            from trimesh.creation import triangulate_polygon
            loop_data = []
            for poly in polygons:
                if poly.is_empty or poly.area < 1e-6:
                    continue
                area = poly.area
                perimeter = poly.length

                # Proper constrained triangulation (respects concavity)
                try:
                    tri_verts_2d, tri_faces = triangulate_polygon(poly)
                except Exception:
                    continue

                # Transform back to 3D
                tri_verts_3d = trimesh.transformations.transform_points(
                    np.column_stack([tri_verts_2d, np.zeros(len(tri_verts_2d))]),
                    transform,
                )
                centroid_3d = tri_verts_3d.mean(axis=0)

                loop_data.append({
                    "vertices": tri_verts_3d,
                    "faces": tri_faces,
                    "area": area,
                    "perimeter": perimeter,
                    "centroid": centroid_3d,
                })

            if not loop_data:
                continue

            results.append((i, centerline[i], normals[i], loop_data))

            if (i + 1) % 20 == 0 or i == n - 1:
                log.info(f"  {label_prefix}Plane {i+1}/{n}... "
                         f"{len(loop_data)} loop(s)")

        except Exception as e:
            log.debug(f"  {label_prefix}Plane {i+1}/{n}... error: {e}")
            continue

    return results


def generate_csa_plot(csv_path, output_path, name):
    """CSA plot colored by region."""
    df = pd.read_csv(csv_path)
    # Only plot individual regions (not NasalCombined to avoid double-counting)
    df_plot = df[~df["region"].isin(["NasalCombined"])].copy()

    region_colors = {
        "DescendingAirway": "gray",
        "LeftNose": "blue",
        "RightNose": "red",
        "NasalCombined": "purple",
        "Mouth": "orange",
    }

    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)

    for region, grp in df_plot.groupby("region"):
        c = region_colors.get(region, "black")
        axes[0].plot(grp["arc_length_mm"], grp["area_mm2"],
                     "-o", color=c, markersize=2, linewidth=1, label=region, alpha=0.8)
        if "hydraulic_diameter_mm" in grp.columns:
            axes[1].plot(grp["arc_length_mm"], grp["hydraulic_diameter_mm"],
                         "-o", color=c, markersize=2, linewidth=1, label=region, alpha=0.8)

    axes[0].set_ylabel("Cross-Sectional Area (mm²)")
    axes[0].set_title(f"{name} — CSA by Region")
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel("Arc Length from Trachea (mm)")
    axes[1].set_ylabel("Hydraulic Diameter (mm)")
    axes[1].set_title(f"{name} — Hydraulic Diameter by Region")
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def generate_summary(csv_path, output_path, name, refine_result):
    """Per-region summary statistics."""
    df = pd.read_csv(csv_path)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"Bifurcation CSA Summary: {name}\n{'='*50}\n\n")

        if refine_result:
            f.write(f"Septum position: {refine_result['septum_pos']}\n")
            f.write(f"Merged centerline: {len(refine_result['merged_centerline'])} pts\n\n")

        for region in ["DescendingAirway", "LeftNose", "RightNose",
                       "NasalCombined", "Mouth"]:
            sub = df[df["region"] == region]
            if len(sub) == 0:
                continue
            f.write(f"\n--- {region} ({len(sub)} planes) ---\n")
            f.write(f"  Arc length: {sub['arc_length_mm'].min():.1f} - "
                    f"{sub['arc_length_mm'].max():.1f} mm\n")
            f.write(f"  Area: min={sub['area_mm2'].min():.2f}, "
                    f"max={sub['area_mm2'].max():.2f}, "
                    f"mean={sub['area_mm2'].mean():.2f} mm²\n")
            if "hydraulic_diameter_mm" in sub.columns:
                f.write(f"  Hyd. diameter: min={sub['hydraulic_diameter_mm'].min():.2f}, "
                        f"max={sub['hydraulic_diameter_mm'].max():.2f}, "
                        f"mean={sub['hydraulic_diameter_mm'].mean():.2f} mm\n")


# ---------------------------------------------------------------------------
# Save centerline as VTK
# ---------------------------------------------------------------------------

def save_centerline_vtk(points, filepath):
    """Write a polyline centerline as legacy VTK."""
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Centerline\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        n = len(points)
        f.write(f"POINTS {n} float\n")
        for pt in points:
            f.write(f"{pt[0]} {pt[1]} {pt[2]}\n")
        f.write(f"LINES 1 {n + 1}\n")
        f.write(f"{n}")
        for i in range(n):
            f.write(f" {i}")
        f.write("\n")


# ---------------------------------------------------------------------------
# Single-frame slicing (callable from multi-frame pipeline)
# ---------------------------------------------------------------------------

def run_single_frame(
    stl_path,
    output_dir,
    name,
    refine_result,
    left_cl,
    right_cl,
    midline_cl=None,
    log=None,
):
    """Slice a single frame using pre-computed septum refinement.

    This is the core slicing function called by pipeline_bifurcation.py
    for each of the 21 deformed frames. The refine_result (face labels,
    septum, meshes) is computed once on frame 0 and reused here.

    For deformed frames: the STL has same topology as frame 0, so face
    labels from frame 0 apply directly. Centerlines are deformed versions.

    Args:
        stl_path: Path to deformed full airway STL
        output_dir: Where to save per-region plane STLs + CSV
        name: Prefix for output files
        refine_result: Dict from refine_nasal_partition() on frame 0
        left_cl: Deformed left nose centerline (Nx3)
        right_cl: Deformed right nose centerline (Nx3)
        midline_cl: Deformed merged midline centerline (Nx3), optional
        log: Logger instance

    Returns:
        dict with per-region results: {region: [(plane_idx, area, perimeter, ...)]}
    """
    import trimesh
    if log is None:
        log = logging.getLogger(__name__)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    full_mesh = read_stl(str(stl_path))

    from slicer.septum_refine import (
        extract_submesh, _smooth_centerline,
        find_divergence_index, build_matched_curves,
    )

    # Extract deformed region meshes using FACE INDICES from frame 0.
    # FFD preserves topology: same face count, same face order, just moved vertices.
    # No vertex matching needed — face indices transfer directly.
    full_face_labels = refine_result.get("_full_face_labels")
    if full_face_labels is not None and len(full_face_labels) == len(full_mesh.faces):
        # Labels: 0=mouth, 1=left, 2=right, 3=descending
        left_mesh = extract_submesh(full_mesh, full_face_labels == 1)
        right_mesh = extract_submesh(full_mesh, full_face_labels == 2)
        desc_mesh = extract_submesh(full_mesh, full_face_labels == 3)
        no_mouth_mask = full_face_labels > 0  # everything except mouth
        no_mouth_mesh = extract_submesh(full_mesh, no_mouth_mask)
    else:
        # Fallback: use refine_result meshes (frame 0 only)
        left_mesh = refine_result["left_mesh"]
        right_mesh = refine_result["right_mesh"]
        desc_mesh = refine_result["desc_mesh"]
        no_mouth_mesh = trimesh.util.concatenate([left_mesh, right_mesh, desc_mesh])

    septum_pos = refine_result["septum_pos"]
    septum_normal = refine_result["septum_normal"]

    # Build midline if not provided
    if midline_cl is None:
        diverge_idx = find_divergence_index(left_cl, right_cl)
        left_matched, right_matched, midline_avg, _ = build_matched_curves(
            left_cl, right_cl, diverge_idx)
        midline_cl_full = np.vstack([
            (left_cl[:diverge_idx] + right_cl[:diverge_idx]) / 2.0,
            midline_avg[1:]
        ])
        midline_cl = _smooth_centerline(midline_cl_full)

    # Truncate midline at septum + resample to 2mm
    dists_to_septum = np.linalg.norm(midline_cl - septum_pos, axis=1)
    septum_cl_idx = int(np.argmin(dists_to_septum))
    desc_cl_raw = midline_cl[:septum_cl_idx]

    from scipy.interpolate import interp1d as scipy_interp1d
    arc = compute_arc_length(desc_cl_raw)
    total_len = arc[-1] if len(arc) > 0 else 0
    if total_len > 0:
        n_resampled = max(int(total_len / 2.0), 10)
        t_uniform = np.linspace(0, total_len, n_resampled)
        desc_cl = scipy_interp1d(arc, desc_cl_raw, axis=0, kind='cubic')(t_uniform)
    else:
        desc_cl = desc_cl_raw

    # Create output dirs
    for d in ["DescendingAirway_Planes", "LeftNose_Planes", "RightNose_Planes"]:
        (output_dir / d).mkdir(parents=True, exist_ok=True)

    all_rows = []

    # ── Slice descending ──
    no_mouth_combined = trimesh.util.concatenate([left_mesh, right_mesh, desc_mesh])
    desc_results = slice_along_centerline(
        no_mouth_combined, desc_cl, None, log, label_prefix="Desc ")

    desc_positions = np.array([pos for _, pos, _ in desc_results])
    desc_arc = compute_arc_length(desc_positions) if len(desc_positions) > 0 else np.array([])
    desc_count = 0

    for k, (plane_idx, pos, loops) in enumerate(desc_results):
        arc_len = float(desc_arc[k]) if k < len(desc_arc) else 0.0
        merged = merge_loops(loops)
        if merged is None:
            continue
        stl_path_out = output_dir / "DescendingAirway_Planes" / f"{name}-Desc-{desc_count:03d}.stl"
        write_cross_section_stl(str(stl_path_out), merged.vertices, merged.faces)
        desc_count += 1

        hyd_diam = DiameterCalculator.compute_hydraulic_diameter(merged.area, merged.perimeter)
        equiv_diam = DiameterCalculator.compute_equivalent_diameter(merged.area)
        all_rows.append({
            "plane_index": plane_idx, "region": "DescendingAirway",
            "arc_length_mm": arc_len, "arc_length_branch_mm": arc_len,
            "area_mm2": merged.area, "perimeter_mm": merged.perimeter,
            "hydraulic_diameter_mm": hyd_diam, "equivalent_diameter_mm": equiv_diam,
            "centroid_x": merged.centroid[0], "centroid_y": merged.centroid[1],
            "centroid_z": merged.centroid[2],
            "centerline_source": "nasal_midline", "is_valid": True,
        })

    # Arc offset at divergence
    diverge_idx = refine_result["diverge_idx"]
    trunk_portion = midline_cl[:diverge_idx + 1] if diverge_idx < len(midline_cl) else midline_cl
    arc_offset = compute_arc_length(trunk_portion)[-1] if len(trunk_portion) > 1 else 0.0

    # ── Slice left nose ──
    left_seg = left_cl[diverge_idx:]
    left_seg_arc = compute_arc_length(left_seg)
    left_results = slice_along_midline(left_mesh, left_seg, log, label_prefix="Left ")
    left_count = 0

    for k, (plane_idx, pos, normal, loop_data) in enumerate(left_results):
        branch_arc = float(left_seg_arc[plane_idx]) if plane_idx < len(left_seg_arc) else 0.0
        arc_len = arc_offset + branch_arc
        total_area = sum(L["area"] for L in loop_data)
        total_perimeter = sum(L["perimeter"] for L in loop_data)
        if len(loop_data) == 1:
            mc = loop_data[0]["centroid"]
        else:
            cs = np.array([L["centroid"] for L in loop_data])
            ar = np.array([L["area"] for L in loop_data])
            mc = np.average(cs, axis=0, weights=ar)

        stl_out = output_dir / "LeftNose_Planes" / f"{name}-Left-{left_count:03d}.stl"
        _save_loops_stl(loop_data, str(stl_out))
        left_count += 1

        hyd_diam = DiameterCalculator.compute_hydraulic_diameter(total_area, total_perimeter)
        equiv_diam = DiameterCalculator.compute_equivalent_diameter(total_area)
        all_rows.append({
            "plane_index": plane_idx, "region": "LeftNose",
            "arc_length_mm": arc_len, "arc_length_branch_mm": branch_arc,
            "area_mm2": total_area, "perimeter_mm": total_perimeter,
            "hydraulic_diameter_mm": hyd_diam, "equivalent_diameter_mm": equiv_diam,
            "centroid_x": mc[0], "centroid_y": mc[1], "centroid_z": mc[2],
            "centerline_source": "left_nose", "is_valid": True,
        })

    # ── Slice right nose ──
    right_seg = right_cl[diverge_idx:]
    right_seg_arc = compute_arc_length(right_seg)
    right_results = slice_along_midline(right_mesh, right_seg, log, label_prefix="Right ")
    right_count = 0

    for k, (plane_idx, pos, normal, loop_data) in enumerate(right_results):
        branch_arc = float(right_seg_arc[plane_idx]) if plane_idx < len(right_seg_arc) else 0.0
        arc_len = arc_offset + branch_arc
        total_area = sum(L["area"] for L in loop_data)
        total_perimeter = sum(L["perimeter"] for L in loop_data)
        if len(loop_data) == 1:
            mc = loop_data[0]["centroid"]
        else:
            cs = np.array([L["centroid"] for L in loop_data])
            ar = np.array([L["area"] for L in loop_data])
            mc = np.average(cs, axis=0, weights=ar)

        stl_out = output_dir / "RightNose_Planes" / f"{name}-Right-{right_count:03d}.stl"
        _save_loops_stl(loop_data, str(stl_out))
        right_count += 1

        hyd_diam = DiameterCalculator.compute_hydraulic_diameter(total_area, total_perimeter)
        equiv_diam = DiameterCalculator.compute_equivalent_diameter(total_area)
        all_rows.append({
            "plane_index": plane_idx, "region": "RightNose",
            "arc_length_mm": arc_len, "arc_length_branch_mm": branch_arc,
            "area_mm2": total_area, "perimeter_mm": total_perimeter,
            "hydraulic_diameter_mm": hyd_diam, "equivalent_diameter_mm": equiv_diam,
            "centroid_x": mc[0], "centroid_y": mc[1], "centroid_z": mc[2],
            "centerline_source": "right_nose", "is_valid": True,
        })

    log.info(f"  Frame {name}: Desc={desc_count}, Left={left_count}, Right={right_count}")
    return all_rows


def _save_loops_stl(loop_data, stl_path):
    """Save loop data as properly triangulated STL."""
    all_verts, all_faces = [], []
    offset = 0
    for L in loop_data:
        if "faces" in L and len(L["vertices"]) >= 3:
            all_verts.append(L["vertices"])
            all_faces.append(L["faces"] + offset)
            offset += len(L["vertices"])
    if all_verts:
        write_cross_section_stl(stl_path,
            np.vstack(all_verts), np.vstack(all_faces))


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Bifurcation-aware CSA slicing with anatomical region labels",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python python_slicer/slice_bifurcation.py surface.stl branches/
  python python_slicer/slice_bifurcation.py surface.stl branches/ --name ENT001 -o ./results
        """,
    )
    parser.add_argument("stl_path", help="Path to full airway STL file")
    parser.add_argument("branches_dir", help="Directory with branch_detector output")
    parser.add_argument("-o", "--output", default="./slice_bifurcation_output")
    parser.add_argument("--name", default=None)
    parser.add_argument("--subject-id", default="")
    parser.add_argument("--window", type=int, default=20)
    args = parser.parse_args()

    stl_path = Path(args.stl_path).resolve()
    branches_dir = Path(args.branches_dir).resolve()
    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    name = args.name or stl_path.stem
    subject_id = args.subject_id

    if not stl_path.exists():
        print(f"ERROR: STL not found: {stl_path}")
        sys.exit(1)
    if not branches_dir.exists():
        print(f"ERROR: Branches dir not found: {branches_dir}")
        sys.exit(1)

    log = setup_logging(str(output_dir), name)
    run_start = time.time()

    log.info("=" * 60)
    log.info("BIFURCATION-AWARE CSA SLICING")
    log.info(f"  Surface: {stl_path.name}")
    log.info(f"  Branches: {branches_dir}")
    log.info(f"  Output: {output_dir}")
    log.info(f"  Name: {name}")
    log.info("=" * 60)

    # ── Phase 1: Discover branches ──────────────────────────
    log.info("\n--- Phase 1: Discovering branches ---")
    branches = discover_branches(branches_dir, subject_id)

    has_left = "LeftNose" in branches and branches["LeftNose"]["centerline"]
    has_right = "RightNose" in branches and branches["RightNose"]["centerline"]
    has_mouth = "Mouth" in branches and branches["Mouth"]["segment_centerline"]
    has_desc = "DescendingAirway" in branches and branches["DescendingAirway"]["mesh"]
    has_both_nostrils = has_left and has_right

    log.info(f"  LeftNose:  {'found' if has_left else 'MISSING'}")
    log.info(f"  RightNose: {'found' if has_right else 'MISSING'}")
    log.info(f"  Mouth:     {'found' if has_mouth else 'MISSING'}")
    log.info(f"  Descending: {'found' if has_desc else 'MISSING'}")

    if not has_left and not has_right:
        log.error("No nasal centerlines found. Cannot proceed.")
        sys.exit(1)

    # Load full mesh
    full_mesh = read_stl(str(stl_path))

    # ── Phase 2: Septum refinement ──────────────────────────
    refine_result = None
    merged_cl = None

    if has_both_nostrils:
        log.info("\n--- Phase 2: Septum refinement ---")
        left_cl = read_vtk_centerline(str(branches["LeftNose"]["centerline"]))
        right_cl = read_vtk_centerline(str(branches["RightNose"]["centerline"]))

        # Use ORIGINAL full mesh + mouth mesh for removal
        # Do NOT use VMTK's L+R+Desc — those have holes from fragment removal
        mouth_mesh_for_removal = None
        if has_mouth and branches["Mouth"]["mesh"]:
            mouth_mesh_for_removal = read_stl(str(branches["Mouth"]["mesh"]))

        refine_result = refine_nasal_partition(
            left_cl, right_cl, full_mesh,
            mouth_mesh=mouth_mesh_for_removal,
        )
        merged_cl = refine_result["merged_centerline"]

        # Save refined outputs
        refine_result["left_mesh"].export(str(output_dir / f"{name}_LeftNose_refined.stl"))
        refine_result["right_mesh"].export(str(output_dir / f"{name}_RightNose_refined.stl"))
        refine_result["desc_mesh"].export(str(output_dir / f"{name}_DescendingAirway_refined.stl"))
        save_centerline_vtk(merged_cl, output_dir / f"{name}_merged_centerline.vtk")
        log.info(f"  Saved refined meshes and merged centerline")
    else:
        # Fallback: use the single available nostril's full-path centerline
        label = "LeftNose" if has_left else "RightNose"
        log.info(f"\n--- Phase 2: Single nostril mode ({label}) ---")
        merged_cl = read_vtk_centerline(str(branches[label]["centerline"]))

    # ── Phase 3a: Slice DESCENDING AIRWAY ────────────────────
    # Use no-mouth mesh + midline centerline + containment filter (simple tube)
    log.info(f"\n--- Phase 3a: Slicing descending airway ---")

    if refine_result:
        no_mouth_mesh = trimesh.util.concatenate([
            refine_result["left_mesh"],
            refine_result["right_mesh"],
            refine_result["desc_mesh"],
        ])
        no_mouth_mesh.export(str(output_dir / f"{name}_no_mouth_mesh.stl"))
        log.info(f"  Saved no-mouth mesh ({len(no_mouth_mesh.faces)} faces)")
    else:
        no_mouth_mesh = full_mesh

    # Use MeshPlaneSlicer with containment check (midline is inside the lumen)
    # Truncate centerline at septum + resample to uniform spacing (~2mm)
    if refine_result:
        septum_pos = refine_result["septum_pos"]
        septum_normal = refine_result["septum_normal"]
        dists_to_septum = np.linalg.norm(merged_cl - septum_pos, axis=1)
        septum_cl_idx = int(np.argmin(dists_to_septum))
        desc_cl_raw = merged_cl[:septum_cl_idx]

        # Resample to uniform ~2mm spacing
        arc = compute_arc_length(desc_cl_raw)
        total_len = arc[-1]
        target_spacing = 2.0  # mm
        n_resampled = max(int(total_len / target_spacing), 10)
        from scipy.interpolate import interp1d
        t_uniform = np.linspace(0, total_len, n_resampled)
        desc_cl = interp1d(arc, desc_cl_raw, axis=0, kind='cubic')(t_uniform)
        log.info(f"  Descending CL: {len(desc_cl)} pts "
                 f"(truncated at septum, resampled to {target_spacing}mm spacing)")
    else:
        desc_cl = merged_cl

    desc_slicer_results = slice_along_centerline(
        no_mouth_mesh, desc_cl, None, log, label_prefix="Desc "
    )
    log.info(f"  Descending: {len(desc_slicer_results)} valid planes")

    # Create output directories
    for region_dir in ["DescendingAirway_Planes", "LeftNose_Planes", "RightNose_Planes"]:
        (output_dir / region_dir).mkdir(parents=True, exist_ok=True)

    # Save descending planes + CSV rows
    all_rows = []
    desc_positions = np.array([pos for _, pos, _ in desc_slicer_results])
    desc_arc = compute_arc_length(desc_positions) if len(desc_positions) > 0 else np.array([])
    desc_count = 0

    for k, (plane_idx, pos, loops) in enumerate(desc_slicer_results):
        arc_len = float(desc_arc[k]) if k < len(desc_arc) else 0.0
        merged = merge_loops(loops)
        if merged is None:
            continue

        # Save STL
        stl_path = output_dir / "DescendingAirway_Planes" / f"{name}-Desc-{desc_count:03d}.stl"
        write_cross_section_stl(str(stl_path), merged.vertices, merged.faces)
        desc_count += 1

        # CSV row
        hyd_diam = DiameterCalculator.compute_hydraulic_diameter(merged.area, merged.perimeter)
        equiv_diam = DiameterCalculator.compute_equivalent_diameter(merged.area)
        all_rows.append({
            "plane_index": plane_idx,
            "region": "DescendingAirway",
            "arc_length_mm": arc_len,
            "arc_length_branch_mm": arc_len,
            "area_mm2": merged.area,
            "perimeter_mm": merged.perimeter,
            "hydraulic_diameter_mm": hyd_diam,
            "equivalent_diameter_mm": equiv_diam,
            "centroid_x": merged.centroid[0],
            "centroid_y": merged.centroid[1],
            "centroid_z": merged.centroid[2],
            "centerline_source": "nasal_midline",
            "is_valid": True,
        })

    log.info(f"  DescendingAirway: {desc_count} planes saved")

    # Compute arc length offset: trachea to divergence point on the merged CL
    arc_offset_at_divergence = 0.0
    if refine_result:
        div_idx = refine_result["diverge_idx"]
        trunk_portion = merged_cl[:div_idx + 1] if div_idx < len(merged_cl) else merged_cl
        if len(trunk_portion) > 1:
            arc_offset_at_divergence = compute_arc_length(trunk_portion)[-1]
        log.info(f"  Arc length offset at divergence: {arc_offset_at_divergence:.1f}mm")

    # ── Phase 3b: Slice LEFT NOSE ─────────────────────────
    # Use refined left mesh + left centerline + raw slicing (keep all loops)
    if refine_result and has_left:
        log.info(f"\n--- Phase 3b: Slicing left nose ---")
        left_nose_seg = left_cl[refine_result["diverge_idx"]:]
        log.info(f"  Left nose CL: {len(left_nose_seg)} pts")

        left_results = slice_along_midline(
            refine_result["left_mesh"], left_nose_seg, log, label_prefix="LeftNose "
        )
        log.info(f"  LeftNose: {len(left_results)} valid planes")

        # Arc length along the FULL nose segment (not just valid results)
        # so plane_idx maps to the correct distance from divergence
        left_seg_arc = compute_arc_length(left_nose_seg)
        left_count = 0

        for k, (plane_idx, pos, normal, loop_data) in enumerate(left_results):
            branch_arc = float(left_seg_arc[plane_idx]) if plane_idx < len(left_seg_arc) else 0.0
            arc_len = arc_offset_at_divergence + branch_arc
            # Merge all loops at this plane (all belong to left nose)
            total_area = sum(L["area"] for L in loop_data)
            total_perimeter = sum(L["perimeter"] for L in loop_data)
            if len(loop_data) == 1:
                merged_centroid = loop_data[0]["centroid"]
            else:
                centroids = np.array([L["centroid"] for L in loop_data])
                areas = np.array([L["area"] for L in loop_data])
                merged_centroid = np.average(centroids, axis=0, weights=areas)

            # Save merged STL (proper triangulation, not fan)
            stl_path = output_dir / "LeftNose_Planes" / f"{name}-Left-{left_count:03d}.stl"
            all_verts, all_faces = [], []
            offset = 0
            for L in loop_data:
                if "faces" in L and len(L["vertices"]) >= 3:
                    all_verts.append(L["vertices"])
                    all_faces.append(L["faces"] + offset)
                    offset += len(L["vertices"])
            if all_verts:
                write_cross_section_stl(str(stl_path),
                    np.vstack(all_verts), np.vstack(all_faces))
            left_count += 1

            hyd_diam = DiameterCalculator.compute_hydraulic_diameter(total_area, total_perimeter)
            equiv_diam = DiameterCalculator.compute_equivalent_diameter(total_area)
            all_rows.append({
                "plane_index": plane_idx,
                "region": "LeftNose",
                "arc_length_mm": arc_len,
                "arc_length_branch_mm": branch_arc,
                "area_mm2": total_area,
                "perimeter_mm": total_perimeter,
                "hydraulic_diameter_mm": hyd_diam,
                "equivalent_diameter_mm": equiv_diam,
                "centroid_x": merged_centroid[0],
                "centroid_y": merged_centroid[1],
                "centroid_z": merged_centroid[2],
                "centerline_source": "left_nose",
                "is_valid": True,
            })

        log.info(f"  LeftNose: {left_count} planes saved")

    # ── Phase 3c: Slice RIGHT NOSE ────────────────────────
    if refine_result and has_right:
        log.info(f"\n--- Phase 3c: Slicing right nose ---")
        right_nose_seg = right_cl[refine_result["diverge_idx"]:]
        log.info(f"  Right nose CL: {len(right_nose_seg)} pts")

        right_results = slice_along_midline(
            refine_result["right_mesh"], right_nose_seg, log, label_prefix="RightNose "
        )
        log.info(f"  RightNose: {len(right_results)} valid planes")

        right_seg_arc = compute_arc_length(right_nose_seg)
        right_count = 0

        for k, (plane_idx, pos, normal, loop_data) in enumerate(right_results):
            branch_arc = float(right_seg_arc[plane_idx]) if plane_idx < len(right_seg_arc) else 0.0
            arc_len = arc_offset_at_divergence + branch_arc
            total_area = sum(L["area"] for L in loop_data)
            total_perimeter = sum(L["perimeter"] for L in loop_data)
            if len(loop_data) == 1:
                merged_centroid = loop_data[0]["centroid"]
            else:
                centroids = np.array([L["centroid"] for L in loop_data])
                areas = np.array([L["area"] for L in loop_data])
                merged_centroid = np.average(centroids, axis=0, weights=areas)

            stl_path = output_dir / "RightNose_Planes" / f"{name}-Right-{right_count:03d}.stl"
            all_verts, all_faces = [], []
            offset = 0
            for L in loop_data:
                if "faces" in L and len(L["vertices"]) >= 3:
                    all_verts.append(L["vertices"])
                    all_faces.append(L["faces"] + offset)
                    offset += len(L["vertices"])
            if all_verts:
                write_cross_section_stl(str(stl_path),
                    np.vstack(all_verts), np.vstack(all_faces))
            right_count += 1

            hyd_diam = DiameterCalculator.compute_hydraulic_diameter(total_area, total_perimeter)
            equiv_diam = DiameterCalculator.compute_equivalent_diameter(total_area)
            all_rows.append({
                "plane_index": plane_idx,
                "region": "RightNose",
                "arc_length_mm": arc_len,
                "arc_length_branch_mm": branch_arc,
                "area_mm2": total_area,
                "perimeter_mm": total_perimeter,
                "hydraulic_diameter_mm": hyd_diam,
                "equivalent_diameter_mm": equiv_diam,
                "centroid_x": merged_centroid[0],
                "centroid_y": merged_centroid[1],
                "centroid_z": merged_centroid[2],
                "centerline_source": "right_nose",
                "is_valid": True,
            })

        log.info(f"  RightNose: {right_count} planes saved")

    # ── Phase 4: Mouth slicing (not implemented yet) ────────
    # TODO: When needed, slice mouth using no-nose mesh + mouth centerline
    log.info("\n--- Phase 4: Mouth slicing (not yet implemented) ---")

    # ── Phase 5: Combine + export ───────────────────────────
    log.info(f"\n--- Phase 5: Exporting ---")

    if len(all_rows) == 0:
        log.error("No valid cross-sections found.")
        sys.exit(1)

    # Unit conversion
    df = pd.DataFrame(all_rows)
    arc_vals = df["arc_length_mm"].values
    area_vals = df["area_mm2"].values
    arc_conv, area_conv = convert_units_if_needed(arc_vals, area_vals)
    df["arc_length_mm"] = arc_conv
    df["area_mm2"] = area_conv

    # CSV
    csv_path = output_dir / f"{name}-Data.csv"
    df.to_csv(csv_path, index=False)
    log.info(f"  CSV: {csv_path.name} ({len(df)} rows)")

    # Region summary
    for region in df["region"].unique():
        n = len(df[df["region"] == region])
        log.info(f"    {region}: {n} rows")

    # Combined All_Planes STLs per region (for easy ParaView viewing)
    import glob
    for region, folder in [("DescendingAirway", "DescendingAirway_Planes"),
                           ("LeftNose", "LeftNose_Planes"),
                           ("RightNose", "RightNose_Planes")]:
        files = sorted(glob.glob(str(output_dir / folder / "*.stl")))
        if files:
            combined = trimesh.util.concatenate(
                [trimesh.load_mesh(f) for f in files])
            combined.export(str(output_dir / f"{region}_All_Planes.stl"))
            log.info(f"  {region}_All_Planes.stl ({len(files)} planes)")

    # Plot
    plot_path = output_dir / f"{name}-CSA_plot.png"
    generate_csa_plot(csv_path, plot_path, name)
    log.info(f"  Plot: {plot_path.name}")

    # Summary
    summary_path = output_dir / f"{name}-summary.txt"
    generate_summary(csv_path, summary_path, name, refine_result)
    log.info(f"  Summary: {summary_path.name}")

    total_time = time.time() - run_start
    log.info(f"\n{'='*60}")
    log.info(f"DONE in {total_time:.1f}s. Results in: {output_dir}")
    log.info(f"{'='*60}")


if __name__ == "__main__":
    main()
