#!/usr/bin/env python3
"""
Quick test to verify Python implementation can read the test data
"""

import sys
from pathlib import Path

# Add slicer to path
sys.path.insert(0, str(Path(__file__).parent))

from slicer.io_utils import read_stl, read_vtk_centerline

def test_load_files():
    """Test loading actual data files"""

    print("="*80)
    print("Testing File I/O with Real Data")
    print("="*80)

    # Test LeftNoseDecending
    print("\n1. Testing LeftNoseDecending files...")

    stl_path = "../LeftNoseDecending/FFD/stl/out_000000.stl"
    vtk_path = "../LeftNoseDecending/FFD/vtk/out_000000.vtk"

    if Path(stl_path).exists():
        try:
            mesh = read_stl(stl_path)
            print(f"✓ Successfully loaded STL")
            print(f"  - Vertices: {len(mesh.vertices)}")
            print(f"  - Faces: {len(mesh.faces)}")
            print(f"  - Bounds: {mesh.bounds}")
            print(f"  - Is watertight: {mesh.is_watertight}")
        except Exception as e:
            print(f"✗ Failed to load STL: {e}")
            return False
    else:
        print(f"✗ STL file not found: {stl_path}")
        return False

    if Path(vtk_path).exists():
        try:
            centerline = read_vtk_centerline(vtk_path)
            print(f"✓ Successfully loaded VTK centerline")
            print(f"  - Points: {len(centerline)}")
            print(f"  - First point: {centerline[0]}")
            print(f"  - Last point: {centerline[-1]}")
            print(f"  - Total length: ~{sum([((centerline[i+1] - centerline[i])**2).sum()**0.5 for i in range(len(centerline)-1)]):.2f} mm")
        except Exception as e:
            print(f"✗ Failed to load VTK: {e}")
            return False
    else:
        print(f"✗ VTK file not found: {vtk_path}")
        return False

    # Test RightNose
    print("\n2. Testing RightNose files...")

    stl_path = "../RightNose/FFD/stl/out_000000.stl"
    vtk_path = "../RightNose/FFD/vtk/out_000000.vtk"

    if Path(stl_path).exists():
        try:
            mesh = read_stl(stl_path)
            print(f"✓ Successfully loaded STL")
            print(f"  - Vertices: {len(mesh.vertices)}")
            print(f"  - Faces: {len(mesh.faces)}")
        except Exception as e:
            print(f"✗ Failed to load STL: {e}")
            return False
    else:
        print(f"⚠ RightNose STL file not found (may not have test data)")

    if Path(vtk_path).exists():
        try:
            centerline = read_vtk_centerline(vtk_path)
            print(f"✓ Successfully loaded VTK centerline")
            print(f"  - Points: {len(centerline)}")
        except Exception as e:
            print(f"✗ Failed to load VTK: {e}")
            return False
    else:
        print(f"⚠ RightNose VTK file not found (may not have test data)")

    print("\n" + "="*80)
    print("✓ All I/O tests passed!")
    print("="*80)

    return True

if __name__ == "__main__":
    success = test_load_files()
    sys.exit(0 if success else 1)
