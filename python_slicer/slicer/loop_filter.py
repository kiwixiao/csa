"""
Loop filter for bifurcation geometries.

Filters cross-section loops by checking if their vertices lie on a branch
surface. Keeps loops that touch the branch, discards loops from other branches.
"""

import numpy as np
from typing import List, Optional

from .mesh_intersection import CrossSection


def filter_loops_by_surface_proximity(cross_sections: List[CrossSection],
                                       branch_surface,
                                       threshold: float = 0.1
                                       ) -> List[CrossSection]:
    """
    Keep loops whose vertices are on the branch surface.
    Discard loops from other branches.

    For each loop, extracts the vertices actually used by its faces,
    computes distance to the branch surface, and keeps the loop if
    any vertex is within threshold distance.

    Args:
        cross_sections: list of CrossSection objects (individual loops)
        branch_surface: trimesh of the branch surface (original open mesh)
        threshold: max distance in mm to count as "on surface" (default 0.1)

    Returns:
        List of CrossSection objects that belong to this branch
    """
    if len(cross_sections) <= 1:
        return cross_sections

    kept = []
    for cs in cross_sections:
        used_ids = np.unique(cs.faces.flatten())
        used_verts = cs.vertices[used_ids]

        _, dists, _ = branch_surface.nearest.on_surface(used_verts)
        n_on_surface = (dists < threshold).sum()

        if n_on_surface > 0:
            kept.append(cs)

    if len(kept) > 0:
        return kept
    else:
        # No loop touches the surface — return all (single-branch region)
        return cross_sections


def merge_loops(cross_sections: List[CrossSection]) -> Optional[CrossSection]:
    """
    Merge multiple cross-section loops into one by summing areas and
    combining geometry.
    """
    if not cross_sections:
        return None
    if len(cross_sections) == 1:
        return cross_sections[0]

    total_area = sum(cs.area for cs in cross_sections)
    total_perimeter = sum(cs.perimeter for cs in cross_sections)

    all_verts = []
    all_faces = []
    offset = 0
    for cs in cross_sections:
        all_verts.append(cs.vertices)
        all_faces.append(cs.faces + offset)
        offset += len(cs.vertices)

    combined_verts = np.vstack(all_verts)
    combined_faces = np.vstack(all_faces)

    centroids = np.array([cs.centroid for cs in cross_sections])
    areas = np.array([cs.area for cs in cross_sections])
    weighted_centroid = np.average(centroids, axis=0, weights=areas)

    all_boundary = [cs.boundary_2d for cs in cross_sections if cs.boundary_2d is not None]
    combined_boundary = np.vstack(all_boundary) if all_boundary else cross_sections[0].boundary_2d

    return CrossSection(
        vertices=combined_verts,
        faces=combined_faces,
        area=total_area,
        centroid=weighted_centroid,
        perimeter=total_perimeter,
        boundary_2d=combined_boundary,
        side=cross_sections[0].side
    )
