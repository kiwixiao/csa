#!/usr/bin/env python3
"""
Automatic Airway Branch Detection using VMTK

Takes a single full airway STL (with open profiles at nostrils, mouth, trachea)
and automatically:
1. Classifies open profiles by anatomy (trachea, nostrils, mouth)
2. Extracts centerlines (1 source → N targets in single VMTK call)
3. Splits the surface mesh into labeled branches
4. Outputs separate STL + centerline VTP files per branch

Designed to run inside the 'vmtk' conda environment.
"""

import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple

import vtk
from vtk.util.numpy_support import vtk_to_numpy

from vmtk import vmtkscripts


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_stl_as_vtk(stl_path: str) -> vtk.vtkPolyData:
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(stl_path))
    reader.Update()
    polydata = reader.GetOutput()
    if polydata.GetNumberOfPoints() == 0:
        raise ValueError(f"Empty mesh loaded from {stl_path}")
    print(f"Loaded STL: {Path(stl_path).name} "
          f"({polydata.GetNumberOfPoints()} pts, {polydata.GetNumberOfCells()} cells)")
    return polydata


def write_stl_from_vtk(polydata: vtk.vtkPolyData, output_path: str):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(output_path))
    writer.SetInputData(polydata)
    writer.SetFileTypeToBinary()
    writer.Update()
    print(f"  Wrote STL: {Path(output_path).name}")


def write_vtp(polydata: vtk.vtkPolyData, output_path: str):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(output_path))
    writer.SetInputData(polydata)
    writer.Update()
    print(f"  Wrote VTP: {Path(output_path).name}")

    # Also save legacy VTK format (ParaView compatible)
    vtk_path = str(output_path).replace(".vtp", ".vtk")
    writer2 = vtk.vtkPolyDataWriter()
    writer2.SetFileName(vtk_path)
    writer2.SetInputData(polydata)
    writer2.Update()
    print(f"  Wrote VTK: {Path(vtk_path).name}")


def centerline_to_vtp(points: np.ndarray) -> vtk.vtkPolyData:
    """Convert Nx3 numpy array to a vtkPolyData polyline"""
    vtk_points = vtk.vtkPoints()
    for pt in points:
        vtk_points.InsertNextPoint(pt[0], pt[1], pt[2])

    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(len(points))
    for i in range(len(points)):
        polyline.GetPointIds().SetId(i, i)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetLines(cells)
    return polydata


# ---------------------------------------------------------------------------
# Open profile detection and classification
# ---------------------------------------------------------------------------

def detect_open_profiles(surface: vtk.vtkPolyData) -> List[Dict]:
    """Detect boundary loops and compute barycenters."""
    fe = vtk.vtkFeatureEdges()
    fe.SetInputData(surface)
    fe.BoundaryEdgesOn()
    fe.FeatureEdgesOff()
    fe.ManifoldEdgesOff()
    fe.NonManifoldEdgesOff()
    fe.Update()

    conn = vtk.vtkPolyDataConnectivityFilter()
    conn.SetInputData(fe.GetOutput())
    conn.SetExtractionModeToAllRegions()
    conn.ColorRegionsOn()
    conn.Update()

    n_profiles = conn.GetNumberOfExtractedRegions()
    region_ids = vtk_to_numpy(conn.GetOutput().GetPointData().GetArray("RegionId"))
    all_points = vtk_to_numpy(conn.GetOutput().GetPoints().GetData())

    profiles = []
    for i in range(n_profiles):
        mask = region_ids == i
        pts = all_points[mask]
        profiles.append({
            "id": i,
            "barycenter": pts.mean(axis=0),
            "n_points": len(pts),
        })
    return profiles


def classify_open_profiles(profiles: List[Dict]
                           ) -> Tuple[List[int], List[int], Dict[int, str]]:
    """
    Classify profiles into anatomical labels.

    Rules:
    1. Trachea = lowest Z. Always present.
    2. Nostrils = pair with similar Z, separated in X. Always together.
    3. Mouth = whatever remains.

    Returns: (source_ids, target_ids, labels)
    """
    print(f"\n--- Classifying {len(profiles)} open profiles ---")
    for p in profiles:
        print(f"  Profile {p['id']}: barycenter={p['barycenter']}")

    if len(profiles) < 2:
        raise ValueError(f"Need at least 2 open profiles, got {len(profiles)}")

    labels = {}

    # Trachea = lowest Z
    z_values = [p["barycenter"][2] for p in profiles]
    trachea = profiles[int(np.argmin(z_values))]
    labels[trachea["id"]] = "Trachea"
    print(f"  Trachea: profile {trachea['id']} (Z={trachea['barycenter'][2]:.1f})")

    remaining = [p for p in profiles if p["id"] != trachea["id"]]

    # Look for nostril pair: similar Z, different X
    if len(remaining) >= 2:
        best_pair = None
        best_z_diff = float("inf")
        for i in range(len(remaining)):
            for j in range(i + 1, len(remaining)):
                z_diff = abs(remaining[i]["barycenter"][2] - remaining[j]["barycenter"][2])
                x_diff = abs(remaining[i]["barycenter"][0] - remaining[j]["barycenter"][0])
                if z_diff < best_z_diff and x_diff > 5.0:
                    best_z_diff = z_diff
                    best_pair = (i, j)

        if best_pair is not None and best_z_diff < 20.0:
            li, ri = best_pair
            if remaining[li]["barycenter"][0] > remaining[ri]["barycenter"][0]:
                li, ri = ri, li
            labels[remaining[li]["id"]] = "LeftNose"
            labels[remaining[ri]["id"]] = "RightNose"
            print(f"  LeftNose: profile {remaining[li]['id']}")
            print(f"  RightNose: profile {remaining[ri]['id']}")
            remaining = [p for p in remaining
                         if p["id"] not in (remaining[li]["id"], remaining[ri]["id"])]

    # Remainder = Mouth
    for p in remaining:
        labels[p["id"]] = "Mouth"
        print(f"  Mouth: profile {p['id']}")

    source_ids = [trachea["id"]]
    target_ids = [p["id"] for p in profiles if p["id"] != trachea["id"]]

    return source_ids, target_ids, labels


# ---------------------------------------------------------------------------
# VMTK pipeline
# ---------------------------------------------------------------------------

def extract_centerlines(surface: vtk.vtkPolyData,
                        profiles: List[Dict],
                        source_ids: List[int],
                        target_ids: List[int],
                        resampling_step: float = 2.0) -> vtk.vtkPolyData:
    """
    Single VMTK centerline call: 1 source (trachea) → N targets.
    Returns polydata with N centerline cells sharing a common trunk.
    """
    print(f"\n--- Extracting centerlines ---")

    # Cap surface (VMTK needs closed surface internally)
    capper = vmtkscripts.vmtkSurfaceCapper()
    capper.Surface = surface
    capper.Interactive = 0
    capper.Method = "simple"
    capper.Execute()
    capped = capper.Surface

    # Build flat coordinate lists
    source_points = []
    for sid in source_ids:
        bc = next(p["barycenter"] for p in profiles if p["id"] == sid)
        source_points.extend(bc.tolist())

    target_points = []
    for tid in target_ids:
        bc = next(p["barycenter"] for p in profiles if p["id"] == tid)
        target_points.extend(bc.tolist())

    print(f"  Source: {source_points}")
    print(f"  Targets ({len(target_ids)}): {target_points}")

    cl = vmtkscripts.vmtkCenterlines()
    cl.Surface = capped
    cl.SeedSelectorName = "pointlist"
    cl.SourcePoints = source_points
    cl.TargetPoints = target_points
    cl.AppendEndPoints = 1
    cl.Resampling = 1
    cl.ResamplingStepLength = resampling_step
    cl.Execute()

    centerlines = cl.Centerlines
    for i in range(centerlines.GetNumberOfCells()):
        print(f"  Line {i}: {centerlines.GetCell(i).GetNumberOfPoints()} points")
    return centerlines


def smooth_centerlines(centerlines: vtk.vtkPolyData) -> vtk.vtkPolyData:
    print("\n--- Smoothing centerlines ---")
    geom = vmtkscripts.vmtkCenterlineGeometry()
    geom.Centerlines = centerlines
    geom.Smoothing = 1
    geom.NumberOfSmoothingIterations = 30
    geom.SmoothingFactor = 0.2
    geom.OutputSmoothed = 1
    geom.Execute()
    return geom.Centerlines


def extract_branches(centerlines: vtk.vtkPolyData) -> vtk.vtkPolyData:
    print("\n--- Extracting branches ---")
    ext = vmtkscripts.vmtkBranchExtractor()
    ext.Centerlines = centerlines
    ext.RadiusArrayName = "MaximumInscribedSphereRadius"
    ext.Execute()
    branched = ext.Centerlines

    gids = vtk_to_numpy(branched.GetCellData().GetArray("GroupIds"))
    blanking = vtk_to_numpy(branched.GetCellData().GetArray("Blanking"))
    print(f"  GroupIds: {np.unique(gids)}")
    print(f"  Branches (Blanking=0): {np.sum(blanking == 0)}, "
          f"Bifurcations (Blanking=1): {np.sum(blanking == 1)}")
    return branched


def clip_surface(surface: vtk.vtkPolyData,
                 branched_centerlines: vtk.vtkPolyData) -> vtk.vtkPolyData:
    print("\n--- Clipping surface ---")
    clip = vmtkscripts.vmtkBranchClipper()
    clip.Surface = surface
    clip.Centerlines = branched_centerlines
    clip.GroupIdsArrayName = "GroupIds"
    clip.Execute()
    clipped = clip.Surface
    print(f"  Clipped: {clipped.GetNumberOfCells()} cells")
    return clipped


# ---------------------------------------------------------------------------
# Centerline cell extraction (raw — no reassembly)
# ---------------------------------------------------------------------------

def extract_centerline_cell(centerlines: vtk.vtkPolyData, cell_idx: int) -> np.ndarray:
    """Extract points from a single centerline cell. Returns Nx3 array."""
    cell = centerlines.GetCell(cell_idx)
    pts = []
    for i in range(cell.GetNumberOfPoints()):
        pid = cell.GetPointId(i)
        pts.append(centerlines.GetPoint(pid))
    return np.array(pts)


def trim_centerline(points: np.ndarray, trim_n: int = 5) -> np.ndarray:
    """Trim points from both ends of a centerline to stay inside the lumen."""
    if len(points) <= 2 * trim_n:
        return points
    return points[trim_n:-trim_n]


def cap_surface(surface: vtk.vtkPolyData) -> vtk.vtkPolyData:
    """Cap open profiles with flat triangulated faces."""
    capper = vmtkscripts.vmtkSurfaceCapper()
    capper.Surface = surface
    capper.Interactive = 0
    capper.Method = "simple"
    capper.Execute()
    capped = capper.Surface
    n_added = capped.GetNumberOfCells() - surface.GetNumberOfCells()
    print(f"  Capped surface: +{n_added} cap cells")
    return capped


def split_centerline_at_bifurcations(full_path: np.ndarray,
                                      branched_centerlines: vtk.vtkPolyData,
                                      cell_idx: int) -> np.ndarray:
    """
    Extract the segment-only portion of a full-path centerline by finding
    where the centerline diverges from the shared trunk (bifurcation point).

    The segment = points from the last bifurcation to the endpoint (nostril/mouth).
    Uses the branched centerlines' Blanking array to find bifurcation boundaries.
    """
    cell_data = branched_centerlines.GetCellData()
    group_ids = vtk_to_numpy(cell_data.GetArray("GroupIds"))
    blanking = vtk_to_numpy(cell_data.GetArray("Blanking"))
    centerline_ids = vtk_to_numpy(cell_data.GetArray("CenterlineIds"))

    # Find cells belonging to this centerline
    my_cells = []
    for ci in range(branched_centerlines.GetNumberOfCells()):
        if centerline_ids[ci] == cell_idx:
            my_cells.append(ci)

    if not my_cells:
        return full_path

    # Find the last bifurcation region (Blanking=1) along this centerline's tracts
    # The segment starts AFTER the last bifurcation
    last_bif_end_pt = None
    for ci in reversed(my_cells):
        if blanking[ci] == 1:
            # Get the last point of this bifurcation cell
            cell = branched_centerlines.GetCell(ci)
            last_pid = cell.GetPointId(cell.GetNumberOfPoints() - 1)
            last_bif_end_pt = np.array(branched_centerlines.GetPoint(last_pid))
            break

    if last_bif_end_pt is None:
        # No bifurcation found — return the whole path
        return full_path

    # Find the closest point on the full path to this bifurcation endpoint
    dists = np.linalg.norm(full_path - last_bif_end_pt, axis=1)
    split_idx = int(np.argmin(dists))

    # The segment is from the bifurcation point to the end (the nostril/mouth)
    segment = full_path[split_idx:]

    return segment


# ---------------------------------------------------------------------------
# Surface extraction by GroupIds
# ---------------------------------------------------------------------------

def collect_group_info(branched_centerlines: vtk.vtkPolyData) -> Dict:
    """Collect centerline_ids for each GroupId from branched centerlines."""
    cell_data = branched_centerlines.GetCellData()
    group_ids = vtk_to_numpy(cell_data.GetArray("GroupIds"))
    blanking = vtk_to_numpy(cell_data.GetArray("Blanking"))
    centerline_ids = vtk_to_numpy(cell_data.GetArray("CenterlineIds"))

    info = {}
    for i in range(len(group_ids)):
        gid = int(group_ids[i])
        if gid not in info:
            info[gid] = {"blanking": int(blanking[i]), "centerline_ids": set()}
        info[gid]["centerline_ids"].add(int(centerline_ids[i]))
    return info


def assign_surface_partitions(group_info: Dict,
                              profiles: List[Dict],
                              profile_labels: Dict[int, str],
                              branched_centerlines: vtk.vtkPolyData
                              ) -> Dict[str, List[int]]:
    """
    Assign GroupIds to anatomical partitions for surface splitting.

    - Terminal GroupIds → matched to profiles by endpoint proximity
    - GroupIds shared by ALL centerlines → DescendingAirway
    - Remaining → assigned to nearest inlet by centerline_id overlap
    """
    print(f"\n--- Assigning surface partitions ---")

    n_centerlines = max(max(info["centerline_ids"]) for info in group_info.values()) + 1

    # Find terminal GroupIds (Blanking=0) closest to each profile barycenter
    branch_gids = [gid for gid, info in group_info.items() if info["blanking"] == 0]

    # Extract endpoint coordinates for each branch GroupId
    cell_data = branched_centerlines.GetCellData()
    all_gids = vtk_to_numpy(cell_data.GetArray("GroupIds"))

    gid_endpoints = {}
    for gid in branch_gids:
        # Find cells belonging to this GroupId, get first and last points
        for ci in range(branched_centerlines.GetNumberOfCells()):
            if int(all_gids[ci]) == gid:
                cell = branched_centerlines.GetCell(ci)
                n = cell.GetNumberOfPoints()
                start = np.array(branched_centerlines.GetPoint(cell.GetPointId(0)))
                end = np.array(branched_centerlines.GetPoint(cell.GetPointId(n - 1)))
                if gid not in gid_endpoints:
                    gid_endpoints[gid] = {"start": start, "end": end}
                else:
                    gid_endpoints[gid]["end"] = end
                break

    # Map each profile to its closest terminal GroupId
    inlet_labels = {}  # gid -> label
    used_gids = set()
    for profile in profiles:
        pid = profile["id"]
        if pid not in profile_labels:
            continue
        label = profile_labels[pid]
        if label == "Trachea":
            label = "DescendingAirway"

        bc = profile["barycenter"]
        best_gid, best_dist = None, float("inf")
        for gid in branch_gids:
            if gid in used_gids:
                continue
            ep = gid_endpoints.get(gid)
            if ep is None:
                continue
            for pt in [ep["start"], ep["end"]]:
                d = np.linalg.norm(pt - bc)
                if d < best_dist:
                    best_dist = d
                    best_gid = gid

        if best_gid is not None:
            inlet_labels[best_gid] = label
            used_gids.add(best_gid)
            print(f"  {label}: GroupId {best_gid} (dist={best_dist:.1f})")

    # Build partitions using simple rule:
    # - Exclusive to 1 centerline → that inlet
    # - Shared by 2+ centerlines → DescendingAirway
    partitions = {}
    inlet_gids = set()
    for gid, label in inlet_labels.items():
        if label == "DescendingAirway":
            continue
        if label not in partitions:
            partitions[label] = []
        partitions[label].append(gid)
        inlet_gids.add(gid)

    # Everything else → DescendingAirway
    # (shared by all centerlines, shared by 2, or bifurcation regions)
    desc_gids = []
    for gid, info in group_info.items():
        if gid in inlet_gids:
            continue
        desc_gids.append(gid)
    partitions["DescendingAirway"] = sorted(desc_gids)

    for label, gids in partitions.items():
        print(f"  {label}: GroupIds {sorted(gids)}")

    return partitions


def extract_surface_by_group_ids(clipped_surface: vtk.vtkPolyData,
                                 group_ids: List[int],
                                 remove_fragments: bool = True) -> vtk.vtkPolyData:
    """Extract surface cells by GroupIds (checks PointData then CellData).
    If remove_fragments=True, keeps only the largest connected region."""
    gid_array = clipped_surface.GetPointData().GetArray("GroupIds")
    if gid_array is None:
        gid_array = clipped_surface.GetCellData().GetArray("GroupIds")
        if gid_array is None:
            print("  WARNING: No GroupIds on surface")
            return vtk.vtkPolyData()
        # CellData path
        surface_gids = vtk_to_numpy(gid_array)
        target = set(group_ids)
        ids = vtk.vtkIdTypeArray()
        ids.SetNumberOfComponents(1)
        for i in range(clipped_surface.GetNumberOfCells()):
            if int(surface_gids[i]) in target:
                ids.InsertNextValue(i)
    else:
        # PointData path — majority vote per cell
        point_gids = vtk_to_numpy(gid_array)
        target = set(group_ids)
        ids = vtk.vtkIdTypeArray()
        ids.SetNumberOfComponents(1)
        for i in range(clipped_surface.GetNumberOfCells()):
            cell = clipped_surface.GetCell(i)
            n = cell.GetNumberOfPoints()
            match = sum(1 for j in range(n) if int(point_gids[cell.GetPointId(j)]) in target)
            if match > n / 2:
                ids.InsertNextValue(i)

    sel_node = vtk.vtkSelectionNode()
    sel_node.SetFieldType(vtk.vtkSelectionNode.CELL)
    sel_node.SetContentType(vtk.vtkSelectionNode.INDICES)
    sel_node.SetSelectionList(ids)

    sel = vtk.vtkSelection()
    sel.AddNode(sel_node)

    ext = vtk.vtkExtractSelection()
    ext.SetInputData(0, clipped_surface)
    ext.SetInputData(1, sel)
    ext.Update()

    geom = vtk.vtkGeometryFilter()
    geom.SetInputData(ext.GetOutput())
    geom.Update()

    result = geom.GetOutput()

    # Remove stray disconnected fragments (vmtkbranchclipper artifact)
    if remove_fragments and result.GetNumberOfCells() > 0:
        conn = vtk.vtkPolyDataConnectivityFilter()
        conn.SetInputData(result)
        conn.SetExtractionModeToLargestRegion()
        conn.Update()
        cleaned = conn.GetOutput()
        removed = result.GetNumberOfCells() - cleaned.GetNumberOfCells()
        if removed > 0:
            print(f"  Extracted {result.GetNumberOfCells()} cells for GroupIds {group_ids}, "
                  f"removed {removed} stray fragments → {cleaned.GetNumberOfCells()}")
        else:
            print(f"  Extracted {result.GetNumberOfCells()} cells for GroupIds {group_ids}")
        return cleaned

    print(f"  Extracted {result.GetNumberOfCells()} cells for GroupIds {group_ids}")
    return result


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

class AirwayBranchDetector:

    def __init__(self, stl_path: str, resampling_step: float = 2.0):
        self.stl_path = Path(stl_path)
        self.resampling_step = resampling_step
        if not self.stl_path.exists():
            raise FileNotFoundError(f"STL not found: {self.stl_path}")

    def run(self, output_dir: str, subject_id: str = ""):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        prefix = f"{subject_id}_" if subject_id else ""

        print("=" * 70)
        print(f"AIRWAY BRANCH DETECTION: {self.stl_path.name}")
        print("=" * 70)

        # Step 1: Load and validate
        surface = read_stl_as_vtk(str(self.stl_path))
        fe = vtk.vtkFeatureEdges()
        fe.SetInputData(surface)
        fe.BoundaryEdgesOn()
        fe.FeatureEdgesOff()
        fe.ManifoldEdgesOff()
        fe.NonManifoldEdgesOff()
        fe.Update()
        n_boundary = fe.GetOutput().GetNumberOfCells()
        print(f"\n  Boundary edges: {n_boundary}")
        if n_boundary == 0:
            raise RuntimeError("Mesh has no open profiles.")

        # Step 2: Detect and classify open profiles
        profiles = detect_open_profiles(surface)
        source_ids, target_ids, profile_labels = classify_open_profiles(profiles)

        # Build target label list (order matches VMTK centerline cell order)
        target_labels = []
        for tid in target_ids:
            label = profile_labels[tid]
            target_labels.append(label)
        print(f"  Target order: {list(zip(target_ids, target_labels))}")

        # Step 3: Extract centerlines (single VMTK call)
        centerlines = extract_centerlines(
            surface, profiles, source_ids, target_ids, self.resampling_step)

        # Step 4: Smooth
        centerlines = smooth_centerlines(centerlines)

        # Step 5: Branch extraction + surface clipping
        branched = extract_branches(centerlines)
        clipped = clip_surface(surface, branched)

        # Step 6: Assign surface partitions by GroupIds
        group_info = collect_group_info(branched)
        surface_partitions = assign_surface_partitions(
            group_info, profiles, profile_labels, branched)

        # Step 7: Save everything
        print(f"\n{'='*70}")
        print(f"SAVING to {output_dir}")
        print(f"{'='*70}")

        n_cells = centerlines.GetNumberOfCells()
        saved = []

        # Extract and trim all full-path centerlines
        full_paths = {}  # label -> trimmed full path (trachea → branch)
        for i, label in enumerate(target_labels):
            if i >= n_cells:
                print(f"  WARNING: No centerline cell for {label}")
                continue
            pts = extract_centerline_cell(centerlines, i)
            full_paths[label] = trim_centerline(pts, trim_n=5)

        # Save combined trimmed centerlines
        all_parts = [centerline_to_vtp(pts) for pts in full_paths.values()]
        if all_parts:
            appender = vtk.vtkAppendPolyData()
            for part in all_parts:
                appender.AddInputData(part)
            appender.Update()
            write_vtp(appender.GetOutput(), str(output_dir / "all_centerlines.vtp"))

        # Save per-target: full-path centerline, segment centerline, surface mesh
        for label, full_pts in full_paths.items():
            part_dir = output_dir / label
            part_dir.mkdir(parents=True, exist_ok=True)

            # Full-path centerline (trachea → this branch)
            full_name = f"Trachea_to_{label}"
            write_vtp(centerline_to_vtp(full_pts),
                      str(part_dir / f"{prefix}{full_name}_centerline.vtp"))

            # Surface mesh
            surf_vtk = None
            if label in surface_partitions:
                surf_vtk = extract_surface_by_group_ids(clipped, surface_partitions[label])
                if surf_vtk.GetNumberOfCells() > 0:
                    write_stl_from_vtk(surf_vtk, str(part_dir / f"{prefix}{label}_mesh.stl"))
                    saved.append(label)

            # Segment centerline (split at last bifurcation point)
            i = target_labels.index(label)
            seg_pts = split_centerline_at_bifurcations(full_pts, branched, i)
            if len(seg_pts) > 1:
                write_vtp(centerline_to_vtp(seg_pts),
                          str(part_dir / f"{prefix}{label}_segment_centerline.vtp"))
                print(f"  {label} segment: {len(seg_pts)} pts "
                      f"(from {len(full_pts)} full path)")

        # Save DescendingAirway surface (no fragment removal — keep all shared regions)
        if "DescendingAirway" in surface_partitions:
            part_dir = output_dir / "DescendingAirway"
            part_dir.mkdir(parents=True, exist_ok=True)
            surf = extract_surface_by_group_ids(
                clipped, surface_partitions["DescendingAirway"],
                remove_fragments=False)
            if surf.GetNumberOfCells() > 0:
                write_stl_from_vtk(surf,
                                   str(part_dir / f"{prefix}DescendingAirway_mesh.stl"))
                saved.append("DescendingAirway")

        print(f"\nDone. Saved {len(saved)} partitions: {saved}")
        return saved


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Detect airway branches from a full airway STL using VMTK")
    parser.add_argument("stl_path", help="Path to full airway STL file")
    parser.add_argument("--output-dir", default="./branch_output")
    parser.add_argument("--subject-id", default="")
    parser.add_argument("--resampling-step", type=float, default=2.0)
    args = parser.parse_args()

    detector = AirwayBranchDetector(args.stl_path, args.resampling_step)
    detector.run(args.output_dir, args.subject_id)
