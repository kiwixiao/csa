#!/usr/bin/env python3
"""
Multi-frame bifurcation CSA pipeline.

Processes a full breathing cycle (21+ deformed frames) through:
1. Branch detection + septum refinement (frame 0 only)
2. Slicing all frames with per-region labeling
3. Mutual plane detection
4. Combined CSV + visualization

Usage:
    python python_slicer/pipeline_bifurcation.py ENT001/
    python python_slicer/pipeline_bifurcation.py ENT001/ --skip-interp
"""

import sys
import argparse
import glob
import logging
import time
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent))
from slicer.io_utils import read_stl, read_vtk_centerline
from slicer.septum_refine import refine_nasal_partition
from slice_bifurcation import (
    run_single_frame, discover_branches, setup_logging,
    generate_csa_plot, generate_summary, save_centerline_vtk,
)
from interpolate_frames import generate_motion_frames


log = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Multi-frame bifurcation CSA pipeline",
    )
    parser.add_argument("subject_dir", help="Subject directory (e.g., ENT001/)")
    parser.add_argument("--subject-id", default="")
    parser.add_argument("--skip-interp", action="store_true",
                        help="Skip interpolation if motion/ already populated")
    parser.add_argument("--start", type=float, default=0)
    parser.add_argument("--stop", type=float, default=2000)
    parser.add_argument("--step", type=float, default=100)
    parser.add_argument("--max-frames", type=int, default=0,
                        help="Max frames to process (0=all, e.g. 5 for testing)")
    parser.add_argument("--frames", type=str, default="",
                        help="Specific frames: '0,5,10' or '0-5' or '0,3,5-10'")
    args = parser.parse_args()

    subject_dir = Path(args.subject_dir).resolve()
    subject_id = args.subject_id or subject_dir.name
    output_dir = subject_dir / "csa_bifurcation"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging(str(output_dir), subject_id)
    pipeline_start = time.time()

    logger.info("=" * 60)
    logger.info("MULTI-FRAME BIFURCATION CSA PIPELINE")
    logger.info(f"  Subject: {subject_dir}")
    logger.info(f"  Subject ID: {subject_id}")
    logger.info("=" * 60)

    # ── Step 0: Branch detection (frame 0) ───────────────────
    logger.info(f"\n--- Step 0: Branch detection ---")
    branches_dir = subject_dir / "branches"

    if not branches_dir.exists():
        logger.info("  Running branch detection on frame 0...")
        from step0_branch_detect import main as branch_detect_main
        old_argv = sys.argv
        sys.argv = ["step0_branch_detect.py",
                     str(subject_dir / "surface" / "frame0.stl"),
                     "--output-dir", str(branches_dir),
                     "--subject-id", subject_id]
        branch_detect_main()
        sys.argv = old_argv
    else:
        logger.info(f"  Branches already exist: {branches_dir}")

    branches = discover_branches(branches_dir, subject_id)
    has_left = "LeftNose" in branches and branches["LeftNose"]["centerline"]
    has_right = "RightNose" in branches and branches["RightNose"]["centerline"]
    has_mouth = "Mouth" in branches and branches["Mouth"].get("mesh")

    if not has_left or not has_right:
        logger.error("Need both LeftNose and RightNose. Cannot proceed.")
        sys.exit(1)

    # ── Step 1: Septum refinement (frame 0) ────────────────
    logger.info(f"\n--- Step 1: Septum refinement ---")
    frame0_stl = subject_dir / "surface" / "frame0.stl"
    full_mesh = read_stl(str(frame0_stl))
    left_cl_f0 = read_vtk_centerline(str(branches["LeftNose"]["centerline"]))
    right_cl_f0 = read_vtk_centerline(str(branches["RightNose"]["centerline"]))

    mouth_mesh = None
    if has_mouth:
        mouth_mesh = read_stl(str(branches["Mouth"]["mesh"]))

    refine_result = refine_nasal_partition(
        left_cl_f0, right_cl_f0, full_mesh, mouth_mesh=mouth_mesh,
    )

    # Save refined outputs + midline CL to branches/ for interpolation
    refine_result["left_mesh"].export(str(output_dir / f"{subject_id}_LeftNose_refined.stl"))
    refine_result["right_mesh"].export(str(output_dir / f"{subject_id}_RightNose_refined.stl"))
    refine_result["desc_mesh"].export(str(output_dir / f"{subject_id}_DescendingAirway_refined.stl"))
    save_centerline_vtk(refine_result["merged_centerline"],
                        output_dir / f"{subject_id}_merged_centerline.vtk")
    # Also save midline CL to branches/ so interpolation can find it
    save_centerline_vtk(refine_result["merged_centerline"],
                        branches_dir / f"{subject_id}_Trachea_to_NoseMidline_centerline.vtk")
    logger.info(f"  Saved refined meshes + merged centerline (also in branches/)")

    # ── Step 2: Generate deformed frames ────────────────────
    motion_dir = subject_dir / "motion"
    stl_dir = motion_dir / "stl"

    n_existing = len(list(stl_dir.glob("*.stl"))) if stl_dir.exists() else 0

    if args.skip_interp and n_existing > 0:
        logger.info(f"\n--- Step 2: Skipping interpolation ({n_existing} STLs exist) ---")
    else:
        logger.info(f"\n--- Step 2: Generating deformed frames ---")
        generate_motion_frames(
            subject_dir=subject_dir,
            start=args.start, stop=args.stop, step=args.step,
        )
        n_existing = len(list(stl_dir.glob("*.stl")))

    logger.info(f"  Motion frames: {n_existing}")

    # ── Step 3: Slice all frames ────────────────────────────
    logger.info(f"\n--- Step 3: Slicing {n_existing} frames ---")

    stl_files = sorted(stl_dir.glob("*.stl"))

    # Frame selection: --frames "0,5,10" or "0-5" or --max-frames 5
    if args.frames:
        selected = set()
        for part in args.frames.split(","):
            part = part.strip()
            if "-" in part:
                a, b = part.split("-", 1)
                selected.update(range(int(a), int(b) + 1))
            else:
                selected.add(int(part))
        stl_files = [f for i, f in enumerate(stl_files) if i in selected]
        logger.info(f"  Selected frames: {sorted(selected)} ({len(stl_files)} files)")
    elif args.max_frames > 0:
        stl_files = stl_files[:args.max_frames]
        logger.info(f"  Limited to {args.max_frames} frames")
    cl_left_dir = motion_dir / "centerlines" / "leftnose"
    cl_right_dir = motion_dir / "centerlines" / "rightnose"
    cl_midline_dir = motion_dir / "centerlines" / "midline"

    all_frame_rows = []
    frame_plane_indices = defaultdict(lambda: defaultdict(set))

    for fi, stl_file in enumerate(stl_files):
        frame_name = stl_file.stem
        frame_output = output_dir / "frames" / frame_name
        logger.info(f"\n  Frame {fi+1}/{len(stl_files)}: {frame_name}")

        # Nose CLs: use frame 0 (nose doesn't move)
        left_cl = left_cl_f0
        right_cl = right_cl_f0

        # Check for morphed nose CLs (if morph_nose was True)
        left_cl_file = cl_left_dir / f"{frame_name}.vtk"
        right_cl_file = cl_right_dir / f"{frame_name}.vtk"
        if left_cl_file.exists():
            left_cl = read_vtk_centerline(str(left_cl_file))
        if right_cl_file.exists():
            right_cl = read_vtk_centerline(str(right_cl_file))

        # Midline CL: always use deformed version (descending airway moves)
        midline_cl_file = cl_midline_dir / f"{frame_name}.vtk"
        midline_cl = read_vtk_centerline(str(midline_cl_file)) if midline_cl_file.exists() else None

        # Slice this frame
        rows = run_single_frame(
            stl_path=str(stl_file),
            output_dir=str(frame_output),
            name=frame_name,
            refine_result=refine_result,
            left_cl=left_cl,
            right_cl=right_cl,
            midline_cl=midline_cl,
            log=logger,
        )

        # Add frame metadata
        for row in rows:
            row["frame_name"] = frame_name
            row["frame_index"] = fi
            # Track plane indices per region per frame
            frame_plane_indices[row["region"]][fi].add(row["plane_index"])

        all_frame_rows.extend(rows)

    logger.info(f"\n  Total rows across all frames: {len(all_frame_rows)}")

    # ── Step 4: Mutual plane detection ──────────────────────
    logger.info(f"\n--- Step 4: Mutual plane detection ---")

    n_frames = len(stl_files)
    mutual_indices = {}

    for region in ["DescendingAirway", "LeftNose", "RightNose"]:
        region_frames = frame_plane_indices[region]
        if not region_frames:
            mutual_indices[region] = set()
            continue

        # Find plane indices present in ALL frames
        all_indices = None
        for fi, indices in region_frames.items():
            if all_indices is None:
                all_indices = indices.copy()
            else:
                all_indices &= indices

        mutual_indices[region] = all_indices or set()
        logger.info(f"  {region}: {len(mutual_indices[region])} mutual planes "
                     f"(from {len(region_frames)} frames)")

        # Save mutual indices
        mut_file = output_dir / f"{region}_mutual_indices.txt"
        with open(mut_file, "w") as f:
            for idx in sorted(mutual_indices[region]):
                f.write(f"{idx}\n")

    # ── Step 5: Filter to mutual planes + combine ───────────
    logger.info(f"\n--- Step 5: Filter to mutual planes ---")

    df = pd.DataFrame(all_frame_rows)
    if len(df) == 0:
        logger.error("No data produced!")
        sys.exit(1)

    # Filter: keep only mutual plane indices
    mask = df.apply(
        lambda row: row["plane_index"] in mutual_indices.get(row["region"], set()),
        axis=1
    )
    df_mutual = df[mask].copy()
    logger.info(f"  Before filter: {len(df)} rows")
    logger.info(f"  After filter:  {len(df_mutual)} rows (mutual planes only)")

    for region in df_mutual["region"].unique():
        n = len(df_mutual[df_mutual["region"] == region])
        logger.info(f"    {region}: {n} rows")

    # Save CSVs
    csv_all = output_dir / f"{subject_id}-Data-all.csv"
    df.to_csv(csv_all, index=False)
    logger.info(f"  Saved: {csv_all.name} (all planes)")

    csv_mutual = output_dir / f"{subject_id}-Data.csv"
    df_mutual.to_csv(csv_mutual, index=False)
    logger.info(f"  Saved: {csv_mutual.name} (mutual planes only)")

    # ── Step 6: Visualize ───────────────────────────────────
    logger.info(f"\n--- Step 6: Visualization ---")

    plot_path = output_dir / f"{subject_id}-CSA_plot.png"
    generate_csa_plot(csv_mutual, plot_path, subject_id)
    logger.info(f"  Plot: {plot_path.name}")

    summary_path = output_dir / f"{subject_id}-summary.txt"
    generate_summary(csv_mutual, summary_path, subject_id, refine_result)
    logger.info(f"  Summary: {summary_path.name}")

    # ── Step 7: Post-processing (per-region analysis + whole-airway visuals) ──
    logger.info(f"\n--- Step 7: Post-processing ---")
    from postprocess_bifurcation import run_postprocessing
    run_postprocessing(output_dir, subject_id)

    total_time = time.time() - pipeline_start
    logger.info(f"\n{'='*60}")
    logger.info(f"DONE in {total_time:.1f}s")
    logger.info(f"  Frames: {n_existing}")
    logger.info(f"  Mutual planes: {sum(len(v) for v in mutual_indices.values())}")
    logger.info(f"  Output: {output_dir}")
    logger.info(f"{'='*60}")


if __name__ == "__main__":
    main()
