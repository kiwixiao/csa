"""
Nasal septum refinement for bifurcation geometries.

Refines VMTK's Voronoi-based L/R nasal surface split by:
1. Building a gap-centered midline between L/R nose centerlines
2. Detecting the posterior septum via 1→2 loop counting
3. Re-partitioning the surface at the septum using the midline as L/R divider

The gap-centered midline is constructed by casting rays left and right from the
averaged midline to find cavity walls, then placing each point at the midpoint
between the two walls. This ensures the midline always lies in the gap between
left and right nasal passages.
"""

import logging
import numpy as np
import trimesh
from collections import Counter, defaultdict
from typing import Dict, Optional, Tuple

from scipy.interpolate import interp1d
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from .geometry import compute_all_plane_normals

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Centerline utilities
# ---------------------------------------------------------------------------

def find_divergence_index(left_cl: np.ndarray, right_cl: np.ndarray,
                          threshold: float = 1.0) -> int:
    """Find the index where left and right centerlines diverge.

    Walks both centerlines from the trachea end, returns the first index
    where the distance between left[i] and right[i] exceeds *threshold* mm.
    """
    max_check = min(len(left_cl), len(right_cl))
    for i in range(max_check):
        if np.linalg.norm(left_cl[i] - right_cl[i]) > threshold:
            return i
    return max_check - 1


def _arc_length_param(pts: np.ndarray) -> np.ndarray:
    """Cumulative arc-length along a polyline, starting at 0."""
    diffs = np.diff(pts, axis=0)
    return np.concatenate([[0], np.cumsum(np.linalg.norm(diffs, axis=1))])


def build_matched_curves(
    left_cl: np.ndarray,
    right_cl: np.ndarray,
    diverge_idx: int,
    n_samples: int = 200,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build arc-length-matched L/R curves and their averaged midline.

    Returns (left_matched, right_matched, midline_avg, t_common) where each
    is sampled at *n_samples* uniformly spaced parameter values in [0, 1].
    """
    left_nose = left_cl[diverge_idx:]
    right_nose = right_cl[diverge_idx:]

    s_left = _arc_length_param(left_nose)
    s_right = _arc_length_param(right_nose)

    t_left = s_left / s_left[-1]
    t_right = s_right / s_right[-1]

    t_common = np.linspace(0, 1, n_samples)

    left_matched = interp1d(t_left, left_nose, axis=0, kind="cubic")(t_common)
    right_matched = interp1d(t_right, right_nose, axis=0, kind="cubic")(t_common)

    midline_avg = (left_matched + right_matched) / 2.0
    return left_matched, right_matched, midline_avg, t_common


# ---------------------------------------------------------------------------
# Gap-centered midline
# ---------------------------------------------------------------------------

def build_gap_centered_midline(
    midline_avg: np.ndarray,
    lr_unit: np.ndarray,
    mesh: trimesh.Trimesh,
    smooth_window: int = 5,
) -> np.ndarray:
    """Construct a midline that lies in the gap between L/R cavity walls.

    For each averaged-midline point, rays are cast in +lr and -lr directions.
    The gap center is the midpoint between the two nearest wall hits.
    The result is smoothed to remove ray-casting jitter.
    """
    midline_gap = np.copy(midline_avg)
    OFFSET = 0.01  # mm, avoids self-intersection

    for i in range(len(midline_avg)):
        origin = midline_avg[i]
        direction = lr_unit[i]

        # Cast right
        r_locs, _, _ = mesh.ray.intersects_location(
            ray_origins=[origin + direction * OFFSET],
            ray_directions=[direction],
        )
        # Cast left
        l_locs, _, _ = mesh.ray.intersects_location(
            ray_origins=[origin - direction * OFFSET],
            ray_directions=[-direction],
        )

        r_pt, l_pt = None, None
        if len(r_locs) > 0:
            dists = np.linalg.norm(r_locs - origin, axis=1)
            r_pt = r_locs[np.argmin(dists)]
        if len(l_locs) > 0:
            dists = np.linalg.norm(l_locs - origin, axis=1)
            l_pt = l_locs[np.argmin(dists)]

        if r_pt is not None and l_pt is not None:
            midline_gap[i] = (l_pt + r_pt) / 2.0

    # Smooth
    smoothed = np.copy(midline_gap)
    for dim in range(3):
        smoothed[:, dim] = uniform_filter1d(midline_gap[:, dim], size=smooth_window)

    return smoothed


# ---------------------------------------------------------------------------
# Septum detection
# ---------------------------------------------------------------------------

def detect_posterior_septum(
    midline: np.ndarray,
    normals: np.ndarray,
    mesh: trimesh.Trimesh,
    min_area: float = 3.0,
) -> Tuple[int, np.ndarray]:
    """Detect the posterior septum by finding the 1→2 loop transition.

    Slices *mesh* at each midline point (walking from bifurcation toward nose)
    and counts significant loops (area > *min_area*). A median filter is applied
    to the loop-count sequence for robustness. Returns (septum_idx, septum_pos).

    If no clean 1→2 transition is found, returns index 0 (the bifurcation point)
    with a warning.
    """
    n = len(midline)
    loop_counts = np.zeros(n, dtype=int)

    for i in range(n):
        try:
            section = mesh.section(plane_origin=midline[i], plane_normal=normals[i])
            if section is None:
                continue
            paths = section.discrete
            n_sig = 0
            for path_pts in paths:
                xy = path_pts[:, :2]
                x, y = xy[:, 0], xy[:, 1]
                area = 0.5 * abs(
                    np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
                    + x[-1] * y[0] - x[0] * y[-1]
                )
                if area > min_area:
                    n_sig += 1
            loop_counts[i] = n_sig
        except Exception:
            continue

    # Median filter to remove single-plane noise
    filtered = median_filter(loop_counts, size=5)

    # Walk from bifurcation (idx 0) toward nose, find first 1→2 transition
    # Skip 0-count entries (partition boundary gaps)
    prev = 0
    for i in range(n):
        cur = int(filtered[i])
        if cur == 0:
            continue
        if prev == 1 and cur >= 2:
            log.info(f"Posterior septum detected at midline index {i} "
                     f"(position {midline[i]})")
            return i, midline[i]
        prev = cur

    # Fallback: no clean transition
    log.warning("No clean 1→2 loop transition found. "
                "Using bifurcation point as septum location.")
    return 0, midline[0]


# ---------------------------------------------------------------------------
# Surface re-partitioning
# ---------------------------------------------------------------------------

def _build_face_adjacency(mesh: trimesh.Trimesh) -> csr_matrix:
    """Build a sparse face-adjacency matrix from shared edges."""
    edge_to_faces = defaultdict(list)
    for fi, face in enumerate(mesh.faces):
        v0, v1, v2 = sorted(face)
        for edge in [(v0, v1), (v0, v2), (v1, v2)]:
            edge_to_faces[edge].append(fi)

    rows, cols = [], []
    for fids in edge_to_faces.values():
        for i in range(len(fids)):
            for j in range(i + 1, len(fids)):
                rows.extend([fids[i], fids[j]])
                cols.extend([fids[j], fids[i]])

    n = len(mesh.faces)
    return csr_matrix(
        (np.ones(len(rows), dtype=bool), (rows, cols)), shape=(n, n)
    )


def split_at_septum(
    mesh: trimesh.Trimesh,
    septum_pos: np.ndarray,
    septum_normal: np.ndarray,
    left_cl: np.ndarray,
    right_cl: np.ndarray,
    min_component_size: int = 500,
) -> np.ndarray:
    """Split mesh at the septum plane using connected components.

    Below the septum → DescendingAirway (label 0).
    Above the septum → connected components naturally separate L and R
    because the septum wall physically divides them.
    The two largest components are identified as LeftNose (1) or RightNose (2)
    by checking proximity to the L/R centerlines.
    Small fragments are assigned to the nearest large component.
    """
    face_centers = mesh.triangles_center

    # Split above/below septum
    signed_dist = np.dot(face_centers - septum_pos, septum_normal)
    above_septum = signed_dist > 0
    above_indices = np.where(above_septum)[0]

    labels = np.zeros(len(mesh.faces), dtype=int)  # 0 = descending

    if len(above_indices) == 0:
        log.warning("No faces above septum!")
        return labels

    # Build adjacency for above-septum faces and find connected components
    log.info(f"  Finding connected components above septum ({len(above_indices)} faces)...")
    above_mesh_adj = _build_face_adjacency(
        trimesh.Trimesh(vertices=mesh.vertices, faces=mesh.faces[above_indices], process=False)
    )
    # We need adjacency in the original mesh's face space for the above-septum subset
    sub_adj = _build_face_adjacency(mesh)
    sub_adj_above = sub_adj[above_indices][:, above_indices]
    n_comp, comp_labels = connected_components(sub_adj_above, directed=False)

    comp_sizes = Counter(comp_labels)
    log.info(f"  {n_comp} connected components above septum")
    for cid, sz in comp_sizes.most_common(5):
        log.info(f"    Component {cid}: {sz} faces")

    # Two largest = LeftNose and RightNose
    top2 = comp_sizes.most_common(2)
    if len(top2) < 2:
        log.warning("Only 1 component above septum — cannot separate L/R")
        labels[above_indices] = 1  # label everything as left
        return labels

    comp_a_id, comp_a_size = top2[0]
    comp_b_id, comp_b_size = top2[1]

    comp_a_faces = above_indices[comp_labels == comp_a_id]
    comp_b_faces = above_indices[comp_labels == comp_b_id]

    # Identify which is left, which is right by centerline proximity
    centroid_a = face_centers[comp_a_faces].mean(axis=0)
    centroid_b = face_centers[comp_b_faces].mean(axis=0)

    # Distance to left and right centerlines (use the nose portion)
    left_tree = cKDTree(left_cl)
    right_tree = cKDTree(right_cl)

    dist_a_left = left_tree.query(centroid_a)[0]
    dist_a_right = right_tree.query(centroid_a)[0]
    dist_b_left = left_tree.query(centroid_b)[0]
    dist_b_right = right_tree.query(centroid_b)[0]

    # Component closer to left CL → LeftNose (1), other → RightNose (2)
    if dist_a_left < dist_a_right:
        left_comp, right_comp = comp_a_id, comp_b_id
    else:
        left_comp, right_comp = comp_b_id, comp_a_id

    labels[above_indices[comp_labels == left_comp]] = 1
    labels[above_indices[comp_labels == right_comp]] = 2

    # Assign small components to nearest large component
    for cid, sz in comp_sizes.items():
        if cid in (left_comp, right_comp):
            continue
        frag_faces = above_indices[comp_labels == cid]
        frag_centroid = face_centers[frag_faces].mean(axis=0)
        d_left = left_tree.query(frag_centroid)[0]
        d_right = right_tree.query(frag_centroid)[0]
        assign_label = 1 if d_left < d_right else 2
        labels[frag_faces] = assign_label
        side_name = "LeftNose" if assign_label == 1 else "RightNose"
        log.info(f"    Small component ({sz} faces) → {side_name}")

    n_left = (labels == 1).sum()
    n_right = (labels == 2).sum()
    n_desc = (labels == 0).sum()
    log.info(f"  Split result: Desc={n_desc}, Left={n_left}, Right={n_right}")

    return labels


def extract_submesh(mesh: trimesh.Trimesh, face_mask: np.ndarray) -> trimesh.Trimesh:
    """Extract a submesh by boolean face mask, reindexing vertices."""
    faces = mesh.faces[face_mask]
    used = np.unique(faces.flatten())
    remap = np.full(len(mesh.vertices), -1, dtype=int)
    remap[used] = np.arange(len(used))
    return trimesh.Trimesh(
        vertices=mesh.vertices[used], faces=remap[faces], process=False
    )


# ---------------------------------------------------------------------------
# Merged centerline
# ---------------------------------------------------------------------------

def build_merged_centerline(
    left_cl: np.ndarray,
    right_cl: np.ndarray,
    midline_nose: np.ndarray,
    diverge_idx: int,
) -> np.ndarray:
    """Build a single centerline: trachea trunk + gap-centered nose midline.

    The trunk is the average of L/R before the divergence point.
    """
    trunk = (left_cl[:diverge_idx] + right_cl[:diverge_idx]) / 2.0
    merged = np.vstack([trunk, midline_nose[1:]])
    return merged


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def _remove_faces_by_vertex_matching(
    full_mesh: trimesh.Trimesh,
    region_mesh: trimesh.Trimesh,
    threshold: float = 0.01,
) -> trimesh.Trimesh:
    """Remove faces from full_mesh that belong to region_mesh using vertex matching.

    Matches vertices between the two meshes via KDTree (nearly instant).
    A face is removed if ALL 3 of its vertices match region_mesh vertices.
    """
    # Build KDTree on region mesh vertices
    region_tree = cKDTree(region_mesh.vertices)

    # For each full mesh vertex, find nearest region vertex
    dists, _ = region_tree.query(full_mesh.vertices)
    is_region_vertex = dists < threshold

    n_matched = is_region_vertex.sum()
    log.info(f"  Vertex matching: {n_matched}/{len(full_mesh.vertices)} vertices matched")

    # Remove faces where ALL 3 vertices are region vertices
    face_all_matched = np.all(is_region_vertex[full_mesh.faces], axis=1)
    keep_mask = ~face_all_matched

    n_removed = face_all_matched.sum()
    log.info(f"  Removing {n_removed} faces, keeping {keep_mask.sum()}")
    return extract_submesh(full_mesh, keep_mask)


def _smooth_centerline(points: np.ndarray, iterations: int = 30,
                       factor: float = 0.2) -> np.ndarray:
    """Smooth a centerline using Laplacian smoothing.

    Replaces each interior point with a weighted average of its neighbors,
    keeping endpoints fixed. This removes zig-zag from the gap-centering.
    """
    pts = points.copy()
    for _ in range(iterations):
        smoothed = pts.copy()
        for i in range(1, len(pts) - 1):
            mid = (pts[i - 1] + pts[i + 1]) / 2.0
            smoothed[i] = pts[i] + factor * (mid - pts[i])
        pts = smoothed
    return pts


def refine_nasal_partition(
    left_cl: np.ndarray,
    right_cl: np.ndarray,
    full_mesh: trimesh.Trimesh,
    mouth_mesh: Optional[trimesh.Trimesh] = None,
    n_samples: int = 200,
) -> Dict:
    """Refine the VMTK nasal partition at the posterior septum.

    Takes the **original full mesh** and removes mouth faces to build a
    clean no-mouth surface with no holes. VMTK's L+R+Desc meshes have
    small fragments removed during extraction, creating holes — so we
    work from the original instead.

    Returns a dict with refined meshes, merged centerline, and septum info.
    """
    log.info("=== Nasal Septum Refinement ===")

    # Step 1: Find divergence
    diverge_idx = find_divergence_index(left_cl, right_cl)
    bif_point = (left_cl[diverge_idx] + right_cl[diverge_idx]) / 2.0
    log.info(f"Bifurcation at CL index {diverge_idx}: {bif_point}")

    # Step 2: Build no-mouth mesh from ORIGINAL geometry
    # We do NOT combine VMTK's L+R+Desc because those have had small
    # fragments removed (holes). Instead, remove mouth from the original.
    if mouth_mesh is not None:
        log.info("Removing mouth faces from original mesh...")
        no_mouth = _remove_faces_by_vertex_matching(full_mesh, mouth_mesh)
    else:
        no_mouth = full_mesh
    log.info(f"No-mouth mesh: {len(no_mouth.faces)} faces "
             f"(from {len(full_mesh.faces)} original)")

    # Step 3: Build matched curves
    left_matched, right_matched, midline_avg, t_common = build_matched_curves(
        left_cl, right_cl, diverge_idx, n_samples
    )
    log.info(f"Matched curves: {n_samples} pts, "
             f"L-R separation {np.linalg.norm(left_matched - right_matched, axis=1).mean():.1f}mm mean")

    # L→R unit vectors
    lr_vectors = right_matched - left_matched
    lr_norms = np.linalg.norm(lr_vectors, axis=1, keepdims=True)
    lr_norms[lr_norms < 1e-6] = 1.0
    lr_unit = lr_vectors / lr_norms

    # Step 4: Build gap-centered midline
    midline_gap = build_gap_centered_midline(midline_avg, lr_unit, no_mouth)
    shift = np.linalg.norm(midline_gap - midline_avg, axis=1)
    log.info(f"Gap-centered midline: shift from avg "
             f"min={shift.min():.2f}, mean={shift.mean():.2f}, max={shift.max():.2f} mm")

    # Step 5: Detect septum
    normals = compute_all_plane_normals(midline_gap, smooth=True)
    septum_idx, septum_pos = detect_posterior_septum(
        midline_gap, normals, no_mouth
    )
    septum_normal = normals[septum_idx]
    dist_from_bif = np.linalg.norm(septum_pos - bif_point)
    log.info(f"Septum at midline index {septum_idx}, "
             f"{dist_from_bif:.1f}mm from bifurcation")

    # Step 6: Split at septum using connected components
    # The septum wall physically separates L and R — no midline projection needed
    face_labels = split_at_septum(
        no_mouth, septum_pos, septum_normal, left_cl, right_cl
    )

    # Step 7: Extract submeshes
    refined_left = extract_submesh(no_mouth, face_labels == 1)
    refined_right = extract_submesh(no_mouth, face_labels == 2)
    refined_desc = extract_submesh(no_mouth, face_labels == 0)

    # Step 8: Clean descending mesh — remove small disconnected fragments
    desc_components = refined_desc.split(only_watertight=False)
    if len(desc_components) > 1:
        desc_components.sort(key=lambda c: len(c.faces), reverse=True)
        removed = sum(len(c.faces) for c in desc_components[1:])
        refined_desc = desc_components[0]
        log.info(f"  Cleaned descending: removed {removed} fragment faces, "
                 f"kept {len(refined_desc.faces)}")

    # Step 9: Build merged centerline and smooth it
    merged_cl = build_merged_centerline(
        left_cl, right_cl, midline_gap, diverge_idx
    )
    merged_cl = _smooth_centerline(merged_cl, iterations=30, factor=0.2)
    log.info(f"Merged centerline: {len(merged_cl)} pts, "
             f"arc length {_arc_length_param(merged_cl)[-1]:.1f}mm (smoothed)")

    log.info("=== Septum refinement complete ===")

    return {
        "left_mesh": refined_left,
        "right_mesh": refined_right,
        "desc_mesh": refined_desc,
        "merged_centerline": merged_cl,
        "septum_pos": septum_pos,
        "septum_idx": septum_idx,
        "septum_normal": septum_normal,
        "midline": midline_gap,
        "lr_unit": lr_unit,
        "left_matched": left_matched,
        "right_matched": right_matched,
        "diverge_idx": diverge_idx,
    }
