# Python Slicer for CSA Template Project

Python replacement for MATLAB airway slicer pipeline.

## 🚀 Quick Start (3 Steps)

### 1. Create conda environment (one-time setup)

```bash
cd /storage1/Qiwei/CSATemplate/python_slicer
./setup_conda.sh
```

### 2. Activate environment

```bash
conda activate csa
```

### 3. Run from project root

```bash
cd /storage1/Qiwei/CSATemplate
./run_slicer.sh OSAMRI037 LeftNoseDecending
```

**That's it!** Results will be in `LeftNoseDecendingSlicedSTLs/`

---

## 📋 What It Does

Processes airway meshes to compute cross-sectional measurements:

**Inputs:**
- STL files: `LeftNoseDecending/FFD/stl/*.stl` (21 deformed meshes)
- VTK files: `LeftNoseDecending/FFD/vtk/*.vtk` (21 centerlines)

**Outputs (per file):**
- Individual cross-section STL files
- Combined planes STL
- CSV with measurements: area, diameters, arc length, etc.

**Summary outputs:**
- Combined measurements CSV
- Summary statistics CSV
- Results pickle file

---

## 🗂️ Project Structure

```
CSATemplate/
├── run_slicer.sh              ← Run from here (project root)
├── test_conda_setup.sh        ← Test conda environment
├── CONDA_QUICKSTART.md        ← Quick reference
├── README_PYTHON_SLICER.md    ← This file
│
├── python_slicer/
│   ├── environment.yml        ← Conda environment spec
│   ├── setup_conda.sh         ← Create conda env
│   ├── main.py                ← Entry point
│   ├── slicer/                ← Python modules
│   │   ├── io_utils.py
│   │   ├── geometry.py
│   │   ├── validation.py
│   │   ├── mesh_intersection.py
│   │   ├── plane_slicer.py
│   │   └── measurements.py
│   └── README.md              ← Full documentation
│
├── LeftNoseDecending/FFD/
│   ├── stl/                   ← Input: 21 STL files
│   └── vtk/                   ← Input: 21 VTK files
│
└── RightNose/FFD/
    ├── stl/                   ← Input: 21 STL files
    └── vtk/                   ← Input: 21 VTK files
```

---

## 📊 Usage Examples

### Process single partition
```bash
conda activate csa
cd /storage1/Qiwei/CSATemplate
./run_slicer.sh OSAMRI037 LeftNoseDecending
```

### Process both partitions
```bash
conda activate csa
cd /storage1/Qiwei/CSATemplate
./run_slicer.sh OSAMRI037 LeftNoseDecending
./run_slicer.sh OSAMRI037 RightNose
```

### Run manually with custom options
```bash
conda activate csa
cd /storage1/Qiwei/CSATemplate
python python_slicer/main.py OSAMRI037 LeftNoseDecending --output ./custom_output_dir
```

---

## 🔧 Conda Environment Management

### Create environment
```bash
cd python_slicer
./setup_conda.sh
```

### Update environment
```bash
conda activate csa
conda env update -f python_slicer/environment.yml
```

### Remove environment
```bash
conda env remove -n csa
```

### Test environment
```bash
./test_conda_setup.sh
```

---

## ✅ Verification

After running, check outputs:

```bash
# Count generated files
ls LeftNoseDecendingSlicedSTLs/*.stl | wc -l

# View measurements
head -20 LeftNoseDecendingSlicedSTLs/out_000000-Data.csv

# View summary
cat OSAMRI037_LeftNoseDecending_summary.csv
```

---

## 📚 Documentation

- **CONDA_QUICKSTART.md** - Quick reference for conda
- **python_slicer/README.md** - Full technical documentation
- **python_slicer/QUICKSTART.md** - Troubleshooting guide
- **python_slicer/IMPLEMENTATION_COMPLETE.md** - Implementation details
- **python_slicer/INPUT_FILES_EXPLAINED.md** - Input file hierarchy

---

## 🆚 Comparison with MATLAB

| Feature | MATLAB | Python |
|---------|--------|--------|
| License | Requires MATLAB | Free (conda) |
| Setup | MATLAB + toolboxes | `conda activate csa` |
| Run | `matlab -r "..."` | `./run_slicer.sh ...` |
| Output | .mat + CSV | CSV + pickle |
| Speed | ~3-7 min | ~2-5 min (estimated) |
| Extensible | Difficult | Easy (Python modules) |

---

## 🐛 Troubleshooting

### "conda: command not found"
Install Anaconda or Miniconda first.

### "csa environment not found"
```bash
cd python_slicer
./setup_conda.sh
```

### "Module not found" errors
```bash
conda activate csa
conda env update -f python_slicer/environment.yml
```

### Many planes skipped
This is normal (10-20%). Check skip reasons in output.

### Results differ from MATLAB
Expected: < 1% difference in areas. See validation guide in python_slicer/README.md

---

## 💡 Tips

1. **Always run from project root** (`/storage1/Qiwei/CSATemplate`)
2. **Activate conda first**: `conda activate csa`
3. **Check test data exists**: `ls LeftNoseDecending/FFD/stl/*.stl`
4. **Use convenience script**: `./run_slicer.sh` handles paths automatically
5. **View progress**: Output shows plane-by-plane progress

---

## 📞 Support

For issues or questions, check:
1. `CONDA_QUICKSTART.md` - Common conda issues
2. `python_slicer/README.md` - Full documentation
3. `./test_conda_setup.sh` - Test your setup

---

## 🎯 Workflow Summary

```
┌─────────────────────────────────────┐
│  One-time: Create conda env        │
│  cd python_slicer                   │
│  ./setup_conda.sh                   │
└─────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────┐
│  Every time: Activate env           │
│  conda activate csa                 │
└─────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────┐
│  Run slicer from project root       │
│  cd /storage1/Qiwei/CSATemplate     │
│  ./run_slicer.sh SUBJECT PARTITION  │
└─────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────┐
│  Check results                      │
│  ls *SlicedSTLs/                    │
└─────────────────────────────────────┘
```

---

**Ready to use!** Run `./test_conda_setup.sh` to verify your setup, then use `./run_slicer.sh` to process data.
