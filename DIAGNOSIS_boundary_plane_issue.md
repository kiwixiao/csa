# Diagnosis: Boundary Plane Intersection Failure

## Problem Statement
Frame 14, plane 0 and plane 147 fail to create valid cross-sections despite having EXACT boundary cap normals.

## Investigation Results

### 1. Normal Alignment (✓ PERFECT)
```
Start cap normal:  [-0.047810, 0.277850, 0.959434]
Plane 0 normal:    [-0.047810, 0.277850, 0.959434]
Angle difference:  0.000000° (EXACT match)

End cap normal:    [0.001318, -0.500344, 0.865826]
Plane 147 normal:  [0.001318, -0.500344, 0.865826]
Angle difference:  0.000000° (EXACT match)
```

**Conclusion**: Normals are PERFECT. Not the issue.

### 2. Position Analysis (✓ TECHNICALLY CORRECT, BUT...)
```
Centerline[0] → Start Cap Plane:
  Signed distance: -0.073768 mm
  Status: INSIDE mesh by 0.074 mm

Centerline[-1] → End Cap Plane:
  Signed distance: 0.004323 mm
  Status: ON the end cap plane (essentially 0)
```

**Conclusion**: Positions are technically correct (inside/on boundary), but EXTREMELY CLOSE to boundaries.

### 3. Slicing Results (✗ FAILURE)
```
Exported planes: plane_001.stl, plane_002.stl, ..., plane_146.stl
Missing: plane_000.stl (skipped - no intersection)
Missing: plane_147.stl (likely skipped or poor quality)
```

**Conclusion**: Plane 0 completely fails to intersect. Plane 147 may also have issues.

## Root Cause

**The planes are TOO CLOSE to boundary openings for reliable mesh intersection.**

- Centerline[0] is only **0.074 mm inside** the mesh start boundary
- Centerline[147] is essentially **ON** the end boundary (0.004 mm)

When you place a slicing plane this close to a boundary opening:
- There's barely any mesh "thickness" to intersect
- Trimesh's slicing algorithm creates degenerate/empty cross-sections
- Numerical precision issues near boundaries cause failures

### Analogy
Imagine trying to slice an orange at its very end - you'd just be grazing the peel with almost nothing to cut through. That's what's happening here.

## Why Boundary Cap Normals Alone Don't Solve It

Having the correct normal (parallel to boundary opening) ensures the plane won't cut through half-open geometry IF there's sufficient mesh to intersect. But if the plane position is too close to the boundary, there's simply not enough mesh material to create a valid cross-section.

## Recommended Solutions

### Option 1: Offset First/Last Planes Inward (RECOMMENDED)
Move the first and last plane positions 1-2mm inward from the boundary while keeping the boundary cap normals:

```python
# For first plane
offset_distance = 1.0  # mm
plane_0_position = centerline[0] + offset_distance * start_cap_normal

# For last plane
plane_N_position = centerline[-1] - offset_distance * end_cap_normal
```

This ensures:
- Planes still use boundary cap normals (won't cut through half-open geometry)
- Planes are far enough inside to create valid cross-sections
- We capture cross-sections very close to boundaries (within 1mm)

### Option 2: Skip Boundary Planes
Simply skip plane 0 and plane 147, starting from plane 1:
- Pros: Simple, no modification needed
- Cons: Lose cross-sectional data near boundaries

### Option 3: Use Multiple Centerline Points
Instead of using just centerline[0], average the first few points:
```python
plane_0_position = np.mean(centerline[0:3], axis=0)
```
- Pros: Naturally moves position inward
- Cons: Position may no longer align with specific centerline index

## Recommended Approach

Implement **Option 1** with 1.0-1.5mm offset:
1. Keep EXACT boundary cap normals for planes 0 and 147
2. Offset positions 1-1.5mm inward along the normal direction
3. This ensures valid intersections while maintaining correct orientation
4. The offset is small enough to still capture near-boundary geometry

## Expected Outcome

After implementing the offset:
- ✓ Plane 0 will create a valid cross-section ~1mm inside the start boundary
- ✓ Plane 147 will create a valid cross-section ~1mm inside the end boundary
- ✓ Both planes will have correct normals (parallel to boundaries)
- ✓ No more cutting through half-open geometry
- ✓ CSA measurements will be available for near-boundary regions

## Files for Reference

- `diagnose_centerline_position.py` - Diagnostic script
- `frame14_boundary_guided_test/` - Test output showing plane 0 missing
- `python_slicer/slicer/geometry.py:139-173` - Boundary normal blending code
