# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Cross-Sectional Area (CSA) analysis pipeline for airway and vascular geometry from CFD simulation data. Three pipelines: legacy single-branch airway, bifurcation-aware airway, and aortic arch with CoA clinical metrics. Pip-installable as `csa-slicer` package, designed to run inside the `mirtk` conda environment alongside the MIRTK morphing pipeline.

## Environment

The primary environment is the `mirtk` conda env (installed via [mirtk-install](https://github.com/kiwixiao/mirtk-install)). It includes all CSA dependencies plus the MIRTK registration toolkit.

```bash
conda activate mirtk
```

A standalone `csa` env can be created from `python_slicer/environment.yml` for development, but production runs use the `mirtk` env.

Branch detection requires a separate `vmtk` conda env (VMTK toolkit). The pipelines call it via `conda run -n vmtk python branch_detector.py`.

## Pip Package

```bash
# Install (editable, for development)
cd python_slicer && pip install -e .

# Install from GitHub (production)
pip install "csa-slicer @ git+https://github.com/kiwixiao/csa.git#subdirectory=python_slicer"
```

## CLI Commands

```bash
# Legacy single-branch airway pipeline (21 time points)
csa-legacy OSAMRI037 LeftNoseDecending

# Bifurcation airway pipeline
csa-bifurcation /path/to/subject/ --subject-id ENT001
csa-bifurcation /path/to/subject/ --subject-id ENT001 --frames 0        # single frame
csa-bifurcation /path/to/subject/ --subject-id ENT001 --skip-interp     # skip morphing

# Aortic pipeline (single time point)
csa-aortic /path/to/subject/ --subject-id COA001
csa-aortic /path/to/subject/ --subject-id COA001 --skip-branch-detect   # reuse branches

# Core slicer only
csa-slice OSAMRI037 LeftNoseDecending --workers 8

# Convert a named multi-solid STL (solid leftnose/rightnose/mouth/outlet/wall)
# into an open-profile mesh + anatomical seeds (standalone / debugging).
csa-uncap /path/to/named_airway.stl --out /path/to/subject/surface/_derived
```

## Running Tests

```bash
cd python_slicer
pytest tests/ -v
```

Tests must run from `python_slicer/` directory (imports use `slicer.*` paths). Tests cover geometry utilities, diameter calculations, unit conversion, and validation.

## Architecture

### Three Pipelines

| Pipeline | Entry Point | Input | Use Case |
|---|---|---|---|
| `csa-legacy` | `pipeline_full.py` | `{Partition}/FFD/stl/` + `vtk/` (21 time points) | Single-branch airway breathing cycle |
| `csa-bifurcation` | `pipeline_bifurcation.py` | `subject/surface/` + auto-detected `branches/` | Airway with nasal L/R split, septum detection |
| `csa-aortic` | `pipeline_aortic.py` | `subject/surface/` + auto-detected `branches/` | Aortic arch with CoA metrics |

All three converge on the same core measurement logic in `slicer/`. They differ in input handling (FFD time series vs branch-detected regions) and post-processing.

### Input Formats

A named multi-solid STL (ASCII, `solid leftnose/rightnose/mouth/outlet/wall … endsolid`
per region, as OpenFOAM `constant/triSurface` writes) dropped into `subject/surface/`
is auto-detected by `csa-bifurcation`. The adapter (`slicer/anatomical_input.py`)
splits it into an open-profile lumen mesh + `anatomical_seeds.json`, which seeds
branch detection directly instead of relying on geometric endpoint search. Names
are matched case-insensitively (with common aliases) and the adapter works with
whichever of leftnose/rightnose/mouth/outlet are present, as long as `wall` and
`outlet` both exist. Use `csa-uncap` standalone to run this conversion for
debugging/inspection without invoking the full pipeline.

### Core Library (`python_slicer/slicer/`)

| Module | Role |
|---|---|
| `plane_slicer.py` | `AirwaySlicer` - orchestrates slicing mesh along centerline |
| `mesh_intersection.py` | Plane-mesh intersection, cross-section triangulation |
| `geometry.py` | 3D rotations, plane normals, arc length, unit conversion |
| `measurements.py` | `DiameterCalculator` - area, perimeter, hydraulic/equivalent diameter, ellipse fitting |
| `validation.py` | Quality checks (centroid proximity, plane enclosure) |
| `io_utils.py` | STL/VTK file I/O, results serialization |
| `boundary_caps.py` | Boundary cap generation for closing open mesh sections |
| `septum_refine.py` | Nasal L/R septum refinement for bifurcation pipeline |
| `loop_filter.py` | Merges nested polygon loops ("doughnut" handling) |

### Key Specialized Modules

- `branch_detector.py` - VMTK-based automatic branch detection (requires vmtk env)
- `aortic_branch_detector.py` - Aortic arch branch identification (BCA, LCCA, LSA)
- `slice_bifurcation.py` - Shared region classification logic used by bifurcation and aortic pipelines
- `interpolate_frames.py` - Generates motion frames using MIRTK `transform-points` (requires mirtk env)
- `detect_incomplete_planes_optimized.py` - Two-stage plane validation (edge length + surface proximity)

### Upstream Connection: mirtk-install

The morphing pipeline (`mirtk-pipeline`) produces deformed meshes. After morphing, `mirtk-prepare-slicer` packages the output into the directory structure expected by CSA pipelines:

```
subject/
  registration/   # FFD transforms
  surface/        # Time-zero mesh (wall.stl or seg_0.stl)
  motion/         # Deformed STLs + centerlines (from mirtk-interpolate)
  branches/       # Auto-created by branch detection (step 0)
```

### Input Data (Legacy Pipeline)

```
{Partition}/FFD/stl/out_000000.stl ... out_002000.stl   (21 deformed meshes)
{Partition}/FFD/vtk/out_000000.vtk ... out_002000.vtk   (21 centerlines)
```

Files numbered in increments of 100 (000000 to 002000). STL and VTK counts must match.

## Key Measurements

Output CSV columns: `arc_length_mm`, `area_mm2`, `perimeter_mm`, `hydraulic_diameter_mm`, `equivalent_diameter_mm`, `major_axis_mm`, `minor_axis_mm`, `diameter_ratio`, `centroid_x/y/z`, `plane_index`, `is_valid`, `file_index`.

## Known Issues

- `MeshPlaneSlicer` containment check bug: skips containment for single loops, returns all loops instead of filtering. Known, not yet fixed.
