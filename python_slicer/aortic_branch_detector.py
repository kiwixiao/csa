#!/usr/bin/env python3
"""
Aortic Branch Detection using VMTK

Takes a single aortic STL (with open profiles at inlet + N outlets)
and automatically:
1. Classifies open profiles: inlet = largest area, outlets = rest
2. Labels outlets generically: Branch_1, Branch_2, ...
3. Extracts centerlines and splits surface into labeled branches

Uses the airway branch detector's VMTK pipeline but with aortic
profile classification. Maps Inlet→Trachea internally so the
parent's save logic (which expects "Trachea") works correctly.
Output folders are renamed after saving.

Designed to run inside the 'vmtk' conda environment.
"""

import os
import shutil
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple

from branch_detector import (
    detect_open_profiles,
    classify_open_profiles,
    read_stl_as_vtk,
    extract_centerlines,
    smooth_centerlines,
    extract_branches,
    clip_surface,
    collect_group_info,
    assign_surface_partitions,
    extract_centerline_cell,
    trim_centerline,
    centerline_to_vtp,
    split_centerline_at_bifurcations,
    extract_surface_by_group_ids,
    write_stl_from_vtk,
    write_vtp,
    AirwayBranchDetector,
)

import vtk
from vtk.util.numpy_support import vtk_to_numpy


def classify_open_profiles_aortic(
    profiles: List[Dict],
) -> Tuple[List[int], List[int], Dict[int, str]]:
    """
    Classify profiles for aortic geometry.

    Rules:
    1. Inlet = profile with largest boundary loop (most points = largest opening)
    2. Outlets = everything else, labeled Branch_1, Branch_2, ... by size

    Internally maps Inlet → "Trachea" so the parent's save logic
    (which expects "Trachea" as source) works. Output folders are
    renamed after the parent finishes saving.

    Returns: (source_ids, target_ids, labels)
    """
    print(f"\n--- Classifying {len(profiles)} open profiles (aortic) ---")
    for p in profiles:
        print(f"  Profile {p['id']}: barycenter={p['barycenter']}, "
              f"n_points={p['n_points']}")

    if len(profiles) < 2:
        raise ValueError(f"Need at least 2 open profiles, got {len(profiles)}")

    labels = {}

    # Inlet = largest boundary loop (most points)
    inlet = max(profiles, key=lambda p: p["n_points"])
    # Map to "Trachea" so parent's save logic works
    labels[inlet["id"]] = "Trachea"
    print(f"  Inlet: profile {inlet['id']} ({inlet['n_points']} pts) "
          f"→ mapped to 'Trachea' for VMTK compat")

    # Outlets = rest, sorted by size (descending)
    outlets = [p for p in profiles if p["id"] != inlet["id"]]
    outlets.sort(key=lambda p: p["n_points"], reverse=True)

    for i, p in enumerate(outlets):
        label = f"Branch_{i + 1}"
        labels[p["id"]] = label
        print(f"  {label}: profile {p['id']} ({p['n_points']} pts)")

    source_ids = [inlet["id"]]
    target_ids = [p["id"] for p in outlets]

    return source_ids, target_ids, labels


class AorticBranchDetector:
    """Aortic-specific branch detector.

    Uses the same VMTK pipeline as AirwayBranchDetector but with
    aortic profile classification. Maps Inlet→Trachea internally,
    then renames output folders after saving.
    """

    def __init__(self, stl_path: str, resampling_step: float = 2.0):
        self.stl_path = Path(stl_path)
        self.resampling_step = resampling_step
        if not self.stl_path.exists():
            raise FileNotFoundError(f"STL not found: {self.stl_path}")

    def run(self, output_dir: str, subject_id: str = ""):
        output_dir = Path(output_dir)
        prefix = f"{subject_id}_" if subject_id else ""

        # Monkey-patch classify_open_profiles for this run
        import branch_detector
        original_classify = branch_detector.classify_open_profiles
        branch_detector.classify_open_profiles = classify_open_profiles_aortic

        try:
            # Use parent's run() — it handles all VMTK + save logic correctly
            detector = AirwayBranchDetector(str(self.stl_path), self.resampling_step)
            result = detector.run(str(output_dir), subject_id)
        finally:
            # Restore original
            branch_detector.classify_open_profiles = original_classify

        # Rename output folders: DescendingAirway → Trunk
        desc_dir = output_dir / "DescendingAirway"
        trunk_dir = output_dir / "Trunk"
        if desc_dir.exists() and not trunk_dir.exists():
            shutil.move(str(desc_dir), str(trunk_dir))
            # Rename files inside
            for f in trunk_dir.iterdir():
                new_name = f.name.replace("DescendingAirway", "Trunk")
                if new_name != f.name:
                    f.rename(trunk_dir / new_name)
            print(f"  Renamed DescendingAirway → Trunk")

        # Rename centerline files: Trachea_to_ → Inlet_to_
        for branch_dir in output_dir.iterdir():
            if not branch_dir.is_dir():
                continue
            for f in branch_dir.iterdir():
                new_name = f.name.replace("Trachea_to_", "Inlet_to_")
                if new_name != f.name:
                    f.rename(branch_dir / new_name)

        print(f"\nAortic branch detection complete.")
        print(f"  Branches: {[d.name for d in sorted(output_dir.iterdir()) if d.is_dir()]}")
        return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Aortic branch detection using VMTK")
    parser.add_argument("stl_path", help="Input aortic STL file")
    parser.add_argument("--output-dir", default="./branches")
    parser.add_argument("--subject-id", default="")
    parser.add_argument("--resampling-step", type=float, default=2.0)
    args = parser.parse_args()

    detector = AorticBranchDetector(args.stl_path, args.resampling_step)
    detector.run(args.output_dir, args.subject_id)
