# Test Data Summary

## Available Test Data

### LeftNoseDecending
- **Location**: `../LeftNoseDecending/FFD/`
- **STL files**: 21 files in `stl/` (out_000000.stl through out_002000.stl)
- **VTK files**: 21 files in `vtk/` (out_000000.vtk through out_002000.vtk)
- **Subject**: OSAMRI037
- **Status**: ✅ Ready to process

### RightNose
- **Location**: `../RightNose/FFD/`
- **STL files**: 21 files in `stl/` (out_000000.stl through out_002000.stl)
- **VTK files**: 21 files in `vtk/` (out_000000.vtk through out_002000.vtk)
- **Subject**: OSAMRI037
- **Status**: ✅ Ready to process

## File Details

### STL Meshes
- **Format**: Binary or ASCII STL
- **Size**: ~3.1 MB per file
- **Content**: 3D airway mesh geometry at different time points/deformations
- **Vertices**: ~50,000-100,000 per mesh
- **Faces**: ~100,000-200,000 per mesh

### VTK Centerlines
- **Format**: ASCII VTK POLYDATA
- **Size**: ~33 KB per file
- **Content**: Centerline points through airway
- **Points**: ~148 points per centerline
- **Coordinates**: 3D (x, y, z) in mm

## Testing Steps

### 1. Quick I/O Test (Verify files can be loaded)

```bash
cd /storage1/Qiwei/CSATemplate/python_slicer
source venv/bin/activate
python test_io.py
```

Expected output:
```
✓ Successfully loaded STL
  - Vertices: XXXXX
  - Faces: XXXXX
✓ Successfully loaded VTK centerline
  - Points: 148
```

### 2. Single File Test (Test slicing on one file pair)

```bash
# Create a test script for single file
python -c "
from slicer import AirwaySlicer
slicer = AirwaySlicer(
    '../LeftNoseDecending/FFD/stl/out_000000.stl',
    '../LeftNoseDecending/FFD/vtk/out_000000.vtk'
)
results = slicer.slice_along_centerline(quality_checks=True)
print(f'Valid planes: {len(results.valid_sections)}')
print(f'Skipped planes: {len(results.skipped_planes)}')
"
```

### 3. Full Pipeline Test (Process all files)

```bash
# LeftNoseDecending
python main.py OSAMRI037 LeftNoseDecending

# RightNose
python main.py OSAMRI037 RightNose
```

Expected outputs:
- `LeftNoseDecendingSlicedSTLs/` directory with:
  - 21 sets of cross-section STL files
  - 21 measurement CSV files
  - 21 combined planes STL files
- `OSAMRI037_LeftNoseDecending_*.csv` summary files
- `OSAMRI037_LeftNoseDecending_results.pkl` combined results

## Validation Against MATLAB

### Expected Outputs from MATLAB
The MATLAB pipeline should have produced:
1. Individual plane STL files: `out_NNNNNN-Planes-XXX.stl`
2. Combined planes STL: `out_NNNNNN-Planes-All.stl`
3. Data CSV files: `out_NNNNNN-Data.csv`
4. MAT file: `OSAMRI037_LeftNoseDecending.mat`

### Comparison Steps

1. **Compare Areas**:
```bash
# Extract area column from MATLAB CSV
cut -d',' -f2 matlab_output/out_000000-Data.csv > matlab_areas.txt

# Extract area column from Python CSV
cut -d',' -f3 python_output/out_000000-Data.csv > python_areas.txt

# Compare
paste matlab_areas.txt python_areas.txt | awk '{print $1, $2, ($1-$2)/$1*100"%"}'
```

2. **Compare Number of Planes**:
```bash
# MATLAB
ls matlab_output/out_000000-Planes-*.stl | wc -l

# Python
ls python_output/out_000000-Planes-*.stl | wc -l
```

3. **Visual Comparison** (using ParaView or similar):
```bash
paraview python_output/out_000000-Planes-All.stl &
paraview matlab_output/out_000000-Planes-All.stl &
```

## Expected Performance

### Processing Time (estimated)
- **Single file pair**: ~5-15 seconds
- **Full subject (21 pairs)**: ~2-5 minutes
- **MATLAB baseline**: ~3-7 minutes

### Expected Success Rate
- **Valid planes per file**: 100-140 out of 148 centerline points
- **Common skip reasons**:
  - `no_intersection`: Plane doesn't cut mesh (ends of centerline)
  - `centroid_too_far`: Cross-section centroid far from centerline
  - `outside_mesh`: Centroid outside mesh bounds

## Troubleshooting

### Issue: "No such file or directory"
Check paths relative to python_slicer directory:
```bash
ls -la ../LeftNoseDecending/FFD/stl/
ls -la ../LeftNoseDecending/FFD/vtk/
```

### Issue: "Module not found"
Activate virtual environment:
```bash
source venv/bin/activate
```

### Issue: Many planes skipped
This is normal. Check the skip reasons in the output. Typical valid plane counts:
- Start of centerline: May skip first 5-10 planes (outside mesh)
- End of centerline: May skip last 5-10 planes (outside mesh)
- Valid middle section: Should have most planes valid

### Issue: Areas don't match MATLAB exactly
Expected tolerance:
- **Area differences**: < 1% (numerical precision)
- **Arc length differences**: < 0.1 mm (accumulated differences)
- **Diameter differences**: < 0.5 mm (ellipse fitting variations)

## Known Differences from MATLAB

1. **Triangulation Algorithm**: Python uses scipy.spatial.Delaunay, MATLAB uses delaunayTriangulation. Minor differences in triangle distribution.

2. **Polygon Simplification**: Different retry strategies when triangulation fails.

3. **Ellipse Fitting**: Python uses least-squares optimization, may give slightly different major/minor axes.

4. **Quality Checks**: Python has additional checks (convex hull) that may skip more planes.

## Next Steps After Testing

1. ✅ Verify I/O works (test_io.py)
2. ✅ Test single file pair
3. ✅ Process full subject
4. ✅ Compare with MATLAB outputs
5. ✅ Document any discrepancies
6. 🔧 Adjust quality thresholds if needed
7. 🔧 Fine-tune triangulation if needed
8. ✅ Validate measurements are clinically reasonable
