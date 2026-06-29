# Quick Start with Conda Environment

## One-Time Setup

### Step 1: Create conda environment "csa"

```bash
cd /storage1/Qiwei/CSATemplate
cd python_slicer
./setup_conda.sh
```

This creates a conda environment with all dependencies installed.

---

## Every Time You Use the Slicer

### Step 2: Activate the environment

```bash
conda activate csa
```

### Step 3: Run from project root

```bash
cd /storage1/Qiwei/CSATemplate
./run_slicer.sh OSAMRI037 LeftNoseDecending
```

Or manually:
```bash
cd /storage1/Qiwei/CSATemplate
conda activate csa
python python_slicer/main.py OSAMRI037 LeftNoseDecending
```

---

## Process Both Partitions

```bash
conda activate csa
cd /storage1/Qiwei/CSATemplate

# Process left
./run_slicer.sh OSAMRI037 LeftNoseDecending

# Process right
./run_slicer.sh OSAMRI037 RightNose
```

---

## Check Results

Results will be in the current directory:

```bash
ls -la LeftNoseDecendingSlicedSTLs/
ls -la OSAMRI037_LeftNoseDecending_*.csv
```

---

## Deactivate Environment

When done:
```bash
conda deactivate
```

---

## Troubleshooting

### If conda environment already exists:
```bash
conda env remove -n csa
cd python_slicer
./setup_conda.sh
```

### If "module not found" errors:
```bash
conda activate csa
conda env update -f python_slicer/environment.yml
```

### Verify installation:
```bash
conda activate csa
python -c "import trimesh, pyvista, shapely, scipy, pandas; print('All modules OK')"
```

---

## Full Workflow

```bash
# One-time setup
cd /storage1/Qiwei/CSATemplate/python_slicer
./setup_conda.sh

# Every time you use it
cd /storage1/Qiwei/CSATemplate
conda activate csa
./run_slicer.sh OSAMRI037 LeftNoseDecending

# When done
conda deactivate
```
