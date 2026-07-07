"""
Tests for Mouth region CSA row-building (Task 7).
"""

import numpy as np
from slice_bifurcation import build_region_row


def test_build_region_row_mouth():
    row = build_region_row(
        plane_idx=3, region="Mouth", arc_len=12.0, branch_arc=12.0,
        area=4.0, perimeter=8.0,
        centroid=np.array([1.0, 2.0, 3.0]), centerline_source="mouth")
    assert row["region"] == "Mouth"
    assert row["area_mm2"] == 4.0
    assert row["hydraulic_diameter_mm"] == 2.0   # 4*area/perimeter = 16/8
    assert row["centroid_z"] == 3.0
    assert row["is_valid"] is True
