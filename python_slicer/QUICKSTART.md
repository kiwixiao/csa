# Quick Start Guide

## 1. Installation (First Time Only)

```bash
cd /storage1/Qiwei/CSATemplate/python_slicer
./setup.sh
```

This will:
- Create a Python virtual environment
- Install all required dependencies
- Verify installation

## 2. Activate Environment

Every time you want to use the slicer, activate the virtual environment:

```bash
cd /storage1/Qiwei/CSATemplate/python_slicer
source venv/bin/activate
```

## 3. Run the Slicer

### Basic Usage

```bash
python main.py <subject_name> <partition>
```

### Examples

Process left nose:
```bash
python main.py DYMOSA801 LeftNoseDecending
```

Process right nose:
```bash
python main.py DYMOSA801 RightNose
```

## 4. Check Results

Results will be in: `<partition>SlicedSTLs/`

Key output files:
- `*-Planes-*.stl` - Individual cross-section meshes
- `*-Planes-All.stl` - All cross-sections combined
- `*-Data.csv` - Measurements (area, diameter, etc.)
- `*_summary.csv` - Summary statistics

## 5. View Measurements

```bash
# View a specific measurement file
column -t -s, LeftNoseDecendingSlicedSTLs/*-Data.csv | head -20

# Or open in Excel/LibreOffice
libreoffice LeftNoseDecendingSlicedSTLs/*-Data.csv
```

## Troubleshooting

### Issue: "No module named 'slicer'"
**Solution**: Make sure you activated the virtual environment:
```bash
source venv/bin/activate
```

### Issue: "VTK directory not found"
**Solution**: Make sure you run from the `python_slicer` directory:
```bash
cd /storage1/Qiwei/CSATemplate/python_slicer
python main.py ...
```

### Issue: "No VTK or STL files found"
**Solution**: Check that your input files exist:
```bash
ls -la ../LeftNoseDecending/FFD/vtk/
ls -la ../LeftNoseDecending/FFD/stl/
```

### Issue: Many "SKIP" messages
**Solution**: This is normal for some planes. Check the summary at the end to see how many valid planes were found.

## Comparing with MATLAB

To compare outputs:

1. Run MATLAB version and save outputs
2. Run Python version
3. Compare CSV files:

```bash
# Compare areas
diff <(cut -d, -f2 matlab_output.csv) <(cut -d, -f3 python_output.csv)
```

## Running Tests

```bash
cd /storage1/Qiwei/CSATemplate/python_slicer
source venv/bin/activate
pytest tests/ -v
```

## Getting Help

```bash
python main.py --help
```

## Workflow Summary

```
1. Place STL files in: ../<partition>/FFD/stl/
2. Place VTK files in: ../<partition>/FFD/vtk/
3. cd python_slicer
4. source venv/bin/activate
5. python main.py <subject> <partition>
6. Check results in: <partition>SlicedSTLs/
```

## Next Steps

Once you have test data:
1. Run the slicer on test data
2. Compare outputs with MATLAB version
3. Verify measurements match
4. Adjust quality thresholds if needed (in validation.py)
5. Report any issues or discrepancies
