import pytest
from pipeline_bifurcation import select_regions

def _b(*names):  # build a fake discover_branches dict
    d = {}
    for n in names:
        d[n] = {"centerline": f"/x/{n}.vtk", "mesh": f"/x/{n}.stl", "segment_centerline": None, "dir": "/x"}
    return d

def test_select_regions_nose_and_mouth():
    r = select_regions(_b("LeftNose", "RightNose", "Mouth", "DescendingAirway"))
    assert r["has_left"] and r["has_right"] and r["has_mouth"]
    assert set(r["inlets"]) == {"LeftNose", "RightNose", "Mouth"}

def test_select_regions_mouth_only():
    r = select_regions(_b("Mouth", "DescendingAirway"))
    assert r["inlets"] == ["Mouth"]
    assert not r["has_left"] and not r["has_right"] and r["has_mouth"]

def test_select_regions_requires_inlet_and_outlet():
    with pytest.raises(ValueError):
        select_regions(_b("DescendingAirway"))         # no inlet
    with pytest.raises(ValueError):
        select_regions(_b("LeftNose", "RightNose"))    # no outlet/descending

def test_select_regions_mouth_without_centerline():
    d = _b("LeftNose", "Mouth", "DescendingAirway")
    d["Mouth"]["centerline"] = None
    r = select_regions(d)
    assert not r["has_mouth"]
