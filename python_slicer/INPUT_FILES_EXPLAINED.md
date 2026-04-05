# Input Files Explanation

## ✅ YES - Python Implementation Reads the CORRECT Input Files!

Your Python slicer is **correctly configured** to read the same input files as the MATLAB slicer.

---

## The Complete Pipeline

### Stage 1: Original Source Files (You Provide These)

Located in: `LeftNoseDecending/FFD/` or `RightNose/FFD/`

```
OSAMRI037_LeftNoseDecending_Scaled.stl          (850 KB)  - Base mesh
OSAMRI037_LeftNoseDecending_Scaled_CL_smooth.vtp (24 KB)  - Base centerline
ffd_1.dof.gz, ffd_2.dof.gz, ..., ffd_6.dof.gz            - Deformation fields
motiontable.csv                                 (11 MB)   - Motion data
staticFrom4D__t10_n7_0.nii.gz                  (189 KB)  - 4D imaging
```

**Purpose**: These represent the airway at a reference state plus time-varying deformations

---

### Stage 2: Preprocessing (Steps 2-3 Generate These)

#### Step 2 (geoFromDof.sh):
- Applies FFD deformations to base mesh and centerline
- Creates 21 time points (0ms, 100ms, 200ms, ..., 2000ms)

**Generates**:
```
LeftNoseDecending/FFD/stl/
├── out_000000.stl  (3.1 MB) - Mesh at t=0ms
├── out_000100.stl  (3.1 MB) - Mesh at t=100ms
├── out_000200.stl  (3.1 MB) - Mesh at t=200ms
...
└── out_002000.stl  (3.1 MB) - Mesh at t=2000ms

LeftNoseDecending/FFD/vtp/
├── out_000000.vtp - Centerline at t=0ms
├── out_000100.vtp - Centerline at t=100ms
...
└── out_002000.vtp - Centerline at t=2000ms
```

#### Step 3 (vtpTovtk.sh):
- Converts VTP → VTK format (MATLAB compatibility)

**Generates**:
```
LeftNoseDecending/FFD/vtk/
├── out_000000.vtk  (33 KB) - Centerline at t=0ms (VTK format)
├── out_000100.vtk  (33 KB) - Centerline at t=100ms (VTK format)
...
└── out_002000.vtk  (33 KB) - Centerline at t=2000ms (VTK format)
```

---

### Stage 3: Slicer Input (What MATLAB & Python Read)

## 🎯 THE ACTUAL INPUT FILES FOR SLICER:

Both MATLAB and Python read from:

```
Input directories:
├── LeftNoseDecending/FFD/stl/*.stl    (21 deformed meshes)
└── LeftNoseDecending/FFD/vtk/*.vtk    (21 corresponding centerlines)
```

**MATLAB code** (SlicerMasterCode.m, lines 11-12):
```matlab
vtkfiles = ['../',path,'/FFD/vtk/*.vtk'];
stlfiles = ['../',path,'/FFD/stl/*.stl'];
```

**Python code** (main.py, lines 28-29):
```python
vtk_dir = Path(f"../{partition}/FFD/vtk/")
stl_dir = Path(f"../{partition}/FFD/stl/")
```

### ✅ PERFECT MATCH! Both read the same files.

---

## What Each File Pair Represents

| File Pair | Time Point | Description |
|-----------|-----------|-------------|
| out_000000.stl + out_000000.vtk | t=0ms | Reference state |
| out_000100.stl + out_000100.vtk | t=100ms | Airway deformed at 100ms |
| out_000200.stl + out_000200.vtk | t=200ms | Airway deformed at 200ms |
| ... | ... | ... |
| out_002000.stl + out_002000.vtk | t=2000ms | End of respiratory cycle |

Each pair represents the airway geometry at a specific time during breathing/swallowing motion.

---

## Your Current Data Status

### ✅ You Have Complete Data Ready!

```bash
LeftNoseDecending/FFD/
├── stl/  (21 files) ✓
├── vtk/  (21 files) ✓
└── vtp/  (21 files) ✓ (not needed by slicer, but available)

RightNose/FFD/
├── stl/  (21 files) ✓
├── vtk/  (21 files) ✓
└── vtp/  (21 files) ✓
```

**This is EXACTLY what the MATLAB slicer uses!**

Your Python implementation will read these same files and produce equivalent outputs.

---

## Summary

### Question: "Which files are being used by the MATLAB slicer?"

**Answer**: The MATLAB slicer reads from:
- `stl/` folder: 21 deformed mesh files (out_NNNNNN.stl)
- `vtk/` folder: 21 centerline files (out_NNNNNN.vtk)

These are **NOT** the original source files - they're **generated** by the preprocessing pipeline (Steps 2-3).

### Question: "Does Python use the same inputs?"

**Answer**: **YES!** The Python implementation is configured to read from the exact same directories:
- `../{partition}/FFD/stl/*.stl`
- `../{partition}/FFD/vtk/*.vtk`

### Ready to Run?

**YES!** Your Python slicer is correctly configured and ready to process the same inputs as MATLAB:

```bash
cd /storage1/Qiwei/CSATemplate/python_slicer
source venv/bin/activate
python main.py OSAMRI037 LeftNoseDecending
```

This will process all 21 STL/VTK pairs, just like MATLAB does. 🚀
