# Python Slicer - Airway Cross-Sectional Analysis

Python implementation of the MATLAB airway slicer pipeline for cross-sectional analysis of airway geometries.

## Overview

This package replaces the MATLAB slicer pipeline (`SlicerMasterCode.m`, `trachSlice_OSA_AB_3.m`, and `sublobeMulti_DG_AB4.m`) with a pure Python implementation. It performs cross-sectional slicing of 3D airway meshes (STL files) along centerlines (VTK files) and computes geometric measurements including area, diameter, and diameter ratios.

## Features

### Core Functionality
- ✅ **STL/VTK File I/O**: Read STL meshes and VTK centerlines
- ✅ **Plane-Mesh Intersection**: Slice 3D meshes with cutting planes
- ✅ **Cross-Section Triangulation**: Constrained Delaunay triangulation
- ✅ **Geometric Measurements**: Area, perimeter, hydraulic diameter, equivalent diameter
- ✅ **Diameter Analysis**: Major/minor axis fitting, diameter ratios
- ✅ **Nested Polygon Handling**: Proper detection and merging of holes ("doughnuts")

### Improvements Over MATLAB
- 🆕 **Convex Hull Validation**: Reject planes outside mesh bounds
- 🆕 **Adjacent Plane Check**: Prevent overlapping cutting planes
- 🆕 **CSV Output**: All measurements in CSV format (no MATLAB required)
- 🆕 **Structured Results**: Clear data classes vs loose variables
- 🆕 **Better Error Handling**: Skip bad planes gracefully with detailed logging
- 🆕 **Comprehensive Diameter Metrics**: Multiple diameter definitions
- 🆕 **Summary Reports**: Automatic generation of summary statistics

## Installation

### 1. Create Virtual Environment (Recommended)

```bash
cd /storage1/Qiwei/CSATemplate/python_slicer
python3 -m venv venv
source venv/bin/activate
```

### 2. Install Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### 3. Verify Installation

```bash
python -c "import trimesh, pyvista, shapely, scipy, pandas; print('All dependencies installed successfully!')"
```

## Usage

### Basic Usage

```bash
cd /storage1/Qiwei/CSATemplate/python_slicer
python main.py <subject_name> <partition>
```

### Examples

Process left nose descending partition:
```bash
python main.py DYMOSA801 LeftNoseDecending
```

Process right nose partition:
```bash
python main.py DYMOSA801 RightNose
```

Specify custom output directory:
```bash
python main.py DYMOSA801 LeftNoseDecending --output /path/to/output
```

### Expected Directory Structure

```
CSATemplate/
├── python_slicer/           # This directory
│   ├── main.py
│   ├── requirements.txt
│   ├── README.md
│   └── slicer/
│       ├── __init__.py
│       ├── io_utils.py
│       ├── geometry.py
│       ├── validation.py
│       ├── mesh_intersection.py
│       ├── plane_slicer.py
│       └── measurements.py
├── LeftNoseDecending/
│   └── FFD/
│       ├── stl/             # Input STL meshes
│       └── vtk/             # Input VTK centerlines
└── RightNose/
    └── FFD/
        ├── stl/
        └── vtk/
```

## Output Files

For each processed file pair, the pipeline generates:

1. **Individual Plane STLs**: `<basename>-Planes-<NNN>.stl`
   - Each cross-sectional plane as a separate STL file
   - Numbered sequentially (001, 002, ...)

2. **Combined Planes STL**: `<basename>-Planes-All.stl`
   - All cross-sections combined into single STL file
   - Useful for visualization

3. **Measurements CSV**: `<basename>-Data.csv`
   - Columns: arc_length_mm, area_mm2, perimeter_mm, hydraulic_diameter_mm, equivalent_diameter_mm, major_axis_mm, minor_axis_mm, diameter_ratio, centroid_x/y/z, is_valid
   - One row per valid cross-section

4. **Summary Files** (generated at end):
   - `<subject>_<partition>_results.pkl`: All results in Python pickle format
   - `<subject>_<partition>_all_measurements.csv`: Combined measurements from all files
   - `<subject>_<partition>_summary.csv`: Summary statistics per file

## Measurements

### Area Metrics
- **Cross-Sectional Area** (mm²): Area enclosed by boundary
- **Perimeter** (mm): Length of boundary

### Diameter Metrics
- **Hydraulic Diameter** (mm): `Dh = 4 × Area / Perimeter`
  - Commonly used in fluid dynamics
- **Equivalent Diameter** (mm): `D = sqrt(4 × Area / π)`
  - Diameter of circle with same area
- **Major Axis** (mm): Longest axis of fitted ellipse
- **Minor Axis** (mm): Shortest axis of fitted ellipse
- **Diameter Ratio**: `Major / Minor`
  - Quantifies elongation (1.0 = circular, >1.0 = elliptical)

### Position Metrics
- **Arc Length** (mm): Cumulative distance along centerline
- **Centroid** (x, y, z): Geometric center of cross-section

## Quality Control

The pipeline includes several quality checks:

1. **Plane Intersection Check**: Ensures adjacent planes don't intersect
2. **Centroid Proximity**: Validates centroid is close to centerline point (<10mm)
3. **Mesh Enclosure**: Checks if cross-section is enclosed by mesh (convex hull)
4. **Triangulation Quality**: Validates mesh quality after triangulation
5. **Area Validation**: Checks for reasonable area values

Planes failing quality checks are skipped and logged with reasons.

## Module Overview

### Core Modules

1. **io_utils.py**: File I/O operations
   - Read/write STL files
   - Read VTK centerlines
   - Write CSV measurements
   - Pickle serialization

2. **geometry.py**: 3D geometry operations
   - Plane normal calculations
   - 3D rotations (Rodrigues formula)
   - Arc length computation
   - Unit conversions

3. **validation.py**: Quality control
   - Convex hull checks
   - Centroid proximity validation
   - Plane intersection detection
   - Mesh quality checks

4. **mesh_intersection.py**: Core slicing logic
   - Plane-mesh intersection
   - Connected component detection
   - Nested polygon handling ("doughnuts")
   - Constrained Delaunay triangulation

5. **plane_slicer.py**: Main slicing algorithm
   - Orchestrates slicing pipeline
   - Manages quality checks
   - Exports results

6. **measurements.py**: Geometric measurements
   - Diameter calculations
   - Ellipse fitting
   - Summary statistics

## Comparison with MATLAB

| Feature | MATLAB | Python |
|---------|--------|--------|
| STL/VTK I/O | Custom readers | trimesh, pyvista |
| Mesh slicing | Manual implementation | trimesh.section() |
| Polygon ops | inpolygon | shapely |
| Triangulation | delaunayTriangulation | scipy.spatial.Delaunay |
| Rotations | AxelRot (Rodrigues) | scipy.spatial.transform.Rotation |
| Output | .mat + CSV | CSV + pickle |
| Validation | Basic proximity | Convex hull + multiple checks |
| Performance | ~2-5 min/subject | ~1-3 min/subject (estimated) |

## Troubleshooting

### Import Errors

If you get import errors, ensure all dependencies are installed:
```bash
pip install -r requirements.txt
```

### VTK File Issues

If VTK files fail to load, check the file format:
```bash
head -20 path/to/file.vtk
```
Should start with `# vtk DataFile Version ...`

### Mesh Quality Warnings

If you see many "no intersection" warnings:
- Check that STL and VTK files are properly aligned
- Verify centerline passes through the mesh
- Inspect mesh with visualization tool

### Empty Results

If no valid planes are produced:
- Check quality thresholds in validation.py
- Try disabling quality checks: modify `quality_checks=False` in main.py
- Visualize inputs with ParaView or similar

## Testing

### Run Tests
```bash
cd python_slicer
pytest tests/
```

### Compare with MATLAB
```bash
python tests/test_against_matlab.py --matlab-dir /path/to/matlab/outputs --python-dir /path/to/python/outputs
```

## Advanced Usage

### Import as Library

```python
from slicer import AirwaySlicer

# Create slicer
slicer = AirwaySlicer("mesh.stl", "centerline.vtk")

# Perform slicing
results = slicer.slice_along_centerline(quality_checks=True)

# Get measurements
profile = results.compute_diameter_profile()
df = profile.to_dataframe()

# Export
slicer.export_results(results, "output_base_name", "output_dir")
```

### Custom Quality Thresholds

Edit `validation.py` to modify thresholds:
- `validate_centroid_proximity`: threshold=10.0 (mm)
- `check_plane_enclosed_by_mesh`: tolerance=1.0 (mm)

### Disable Specific Checks

In `plane_slicer.py`, comment out specific quality checks as needed.

## Contributing

To add new features or fix bugs:

1. Create a new branch
2. Make changes
3. Add tests in `tests/`
4. Update this README
5. Submit pull request

## License

[Add license information]

## Authors

Implementation by Claude Code
Based on original MATLAB code from [original authors]

## Citation

If you use this code, please cite:
```
[Add citation]
```

## Contact

For questions or issues, contact [contact information]
