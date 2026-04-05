"""
Basic tests for Python Slicer modules
"""

import pytest
import numpy as np
from pathlib import Path

# Test imports
def test_imports():
    """Test that all modules can be imported"""
    from slicer import io_utils
    from slicer import geometry
    from slicer import validation
    from slicer import mesh_intersection
    from slicer import plane_slicer
    from slicer import measurements
    assert True


def test_geometry_arc_length():
    """Test arc length calculation"""
    from slicer.geometry import compute_arc_length

    # Simple test: straight line
    points = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [2, 0, 0],
        [3, 0, 0]
    ])

    arc_lengths = compute_arc_length(points)

    assert len(arc_lengths) == 4
    assert arc_lengths[0] == 0.0
    assert arc_lengths[1] == pytest.approx(1.0)
    assert arc_lengths[2] == pytest.approx(2.0)
    assert arc_lengths[3] == pytest.approx(3.0)


def test_geometry_plane_normal():
    """Test plane normal calculation"""
    from slicer.geometry import compute_plane_normal

    # Simple centerline along X-axis
    centerline = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [2, 0, 0]
    ])

    # Normal at middle point should be along X-axis
    normal = compute_plane_normal(centerline, 1, smooth=False)

    assert np.linalg.norm(normal) == pytest.approx(1.0)
    assert normal[0] > 0.99  # Mostly along X


def test_diameter_calculator():
    """Test diameter calculations"""
    from slicer.measurements import DiameterCalculator

    # Test hydraulic diameter
    area = 100.0  # mm^2
    perimeter = 40.0  # mm
    hyd_diam = DiameterCalculator.compute_hydraulic_diameter(area, perimeter)
    assert hyd_diam == pytest.approx(10.0)

    # Test equivalent diameter (circle with area 100 mm^2)
    equiv_diam = DiameterCalculator.compute_equivalent_diameter(area)
    expected = np.sqrt(4 * 100 / np.pi)
    assert equiv_diam == pytest.approx(expected)


def test_diameter_ratio():
    """Test diameter ratio calculation"""
    from slicer.measurements import DiameterCalculator

    ratio = DiameterCalculator.compute_diameter_ratio(10.0, 5.0)
    assert ratio == pytest.approx(2.0)

    # Circle case
    ratio = DiameterCalculator.compute_diameter_ratio(5.0, 5.0)
    assert ratio == pytest.approx(1.0)


def test_unit_conversion():
    """Test unit conversion from meters to millimeters"""
    from slicer.geometry import convert_units_if_needed

    # Data in meters (max arc_length < 1)
    arc_length_m = np.array([0.0, 0.1, 0.2, 0.3])
    area_m2 = np.array([0.0001, 0.0002, 0.0003, 0.0004])

    arc_length_mm, area_mm2 = convert_units_if_needed(arc_length_m, area_m2)

    assert arc_length_mm[1] == pytest.approx(100.0)  # 0.1 m = 100 mm
    assert area_mm2[1] == pytest.approx(200.0)  # 0.0002 m^2 = 200 mm^2


def test_validation_centroid_proximity():
    """Test centroid proximity validation"""
    from slicer.validation import validate_centroid_proximity

    centroid = np.array([1.0, 2.0, 3.0])
    centerline_point = np.array([1.5, 2.5, 3.5])

    # Should be valid (distance < 10mm)
    is_valid = validate_centroid_proximity(centroid, centerline_point, threshold=10.0)
    assert is_valid

    # Should be invalid
    far_point = np.array([100.0, 100.0, 100.0])
    is_valid = validate_centroid_proximity(centroid, far_point, threshold=10.0)
    assert not is_valid


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
