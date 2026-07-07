import json
from pathlib import Path

import numpy as np
import trimesh
from slicer.anatomical_input import read_named_solids, normalize_cap_name, CAP_LABELS, classify_solids, Cap, is_named_multisolid_stl, adapt_named_stl

MULTISOLID = """solid wall
facet normal 0 0 0
 outer loop
  vertex 0 0 0
  vertex 1 0 0
  vertex 0 1 0
 endloop
endfacet
endsolid wall
solid outlet
facet normal 0 0 0
 outer loop
  vertex 0 0 5
  vertex 1 0 5
  vertex 0 1 5
 endloop
endfacet
endsolid outlet
"""

def test_read_named_solids_splits_blocks(tmp_path):
    p = tmp_path / "m.stl"
    p.write_text(MULTISOLID)
    solids = read_named_solids(str(p))
    assert set(solids.keys()) == {"wall", "outlet"}
    assert isinstance(solids["wall"], trimesh.Trimesh)
    assert len(solids["wall"].faces) == 1
    assert np.allclose(solids["outlet"].vertices[:, 2].mean(), 5.0)


def test_normalize_cap_name_case_and_aliases():
    assert normalize_cap_name("LeftNose") == "leftnose"
    assert normalize_cap_name("  RIGHT_NOSE ") == "rightnose"
    assert normalize_cap_name("outlet_1") == "outlet"
    assert normalize_cap_name("Trachea") == "outlet"
    assert normalize_cap_name("Mouth") == "mouth"
    assert normalize_cap_name("wall") == "wall"
    assert {"leftnose", "rightnose", "mouth", "outlet"} == CAP_LABELS


def test_classify_solids_identifies_wall_and_caps(tmp_path):
    # wall = many faces; outlet + leftnose = small caps; junk = unknown name
    stl = ["solid wall"]
    for i in range(6):  # 6 facets -> wall is the largest solid
        stl += ["facet normal 0 0 0", " outer loop",
                f"  vertex {i} 0 0", f"  vertex {i+1} 0 0", f"  vertex {i} 1 0",
                " endloop", "endfacet"]
    stl.append("endsolid wall")
    for name, z in [("outlet", 0.0), ("leftnose", 9.0), ("junkpatch", 3.0)]:
        stl += [f"solid {name}", "facet normal 0 0 1", " outer loop",
                f"  vertex 0 0 {z}", f"  vertex 1 0 {z}", f"  vertex 0 1 {z}",
                " endloop", "endfacet", f"endsolid {name}"]
    p = tmp_path / "m.stl"
    p.write_text("\n".join(stl) + "\n")

    wall, caps, unknown = classify_solids(read_named_solids(str(p)))
    assert len(wall.faces) == 6
    labels = sorted(c.label for c in caps)
    assert labels == ["leftnose", "outlet"]
    assert unknown == ["junkpatch"]
    outlet = next(c for c in caps if c.label == "outlet")
    assert outlet.centroid.shape == (3,)
    assert abs(np.linalg.norm(outlet.normal) - 1.0) < 1e-6


def test_classify_solids_maxface_fallback_and_area_weighting():
    # No solid named "wall" -> classify_solids must fall back to the solid
    # with the most faces. A box has 12 faces, easily beating a 2-face cap.
    lumen = trimesh.creation.box(extents=[10, 10, 10])

    # "outlet" cap built from two triangles of very different area, both in
    # the z=0 plane, so area-weighting of the centroid is actually exercised
    # (a single-facet cap would make the area-weighted average tautological).
    #   triangle 1: (0,0,0),(1,0,0),(0,1,0)   -> area 0.5, center x = 1/3
    #   triangle 2: (11,0,0),(15,0,0),(11,4,0) -> area 8.0 (16x larger), center x = 37/3
    # simple (unweighted) average of the two centers: x = (1/3 + 37/3) / 2 = 19/3 ~= 6.333
    # area-weighted average (weights 0.5 vs 8):      x = (0.5*(1/3) + 8*(37/3)) / 8.5 ~= 11.627
    vertices = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [11.0, 0.0, 0.0],
        [15.0, 0.0, 0.0],
        [11.0, 4.0, 0.0],
    ])
    faces = [[0, 1, 2], [3, 4, 5]]
    outlet = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    solids = {"lumen": lumen, "outlet": outlet}
    wall, caps, unknown = classify_solids(solids)

    assert wall is lumen  # max-faces fallback chosen since nothing is named "wall"
    assert unknown == []

    cap = next(c for c in caps if c.label == "outlet")
    assert cap.label == "outlet"

    simple_avg_x = 19 / 3  # ~= 6.333
    area_weighted_x = (0.5 * (1 / 3) + 8 * (37 / 3)) / 8.5  # ~= 11.627
    assert simple_avg_x < area_weighted_x
    # threshold strictly between the simple average and the area-weighted value
    assert cap.centroid[0] > 9.0

    assert abs(np.linalg.norm(cap.normal) - 1) < 1e-6


def _write_airway(tmp_path, names_z):
    """Write an ASCII multi-solid STL with a wall and named caps."""
    stl = ["solid wall"]
    for i in range(6):
        stl += ["facet normal 0 0 0", " outer loop",
                f"  vertex {i} 0 0", f"  vertex {i+1} 0 0", f"  vertex {i} 1 0",
                " endloop", "endfacet"]
    stl.append("endsolid wall")
    for name, z in names_z:
        stl += [f"solid {name}", "facet normal 0 0 1", " outer loop",
                f"  vertex 0 0 {z}", f"  vertex 1 0 {z}", f"  vertex 0 1 {z}",
                " endloop", "endfacet", f"endsolid {name}"]
    p = tmp_path / "airway.stl"
    p.write_text("\n".join(stl) + "\n")
    return p


def test_is_named_multisolid_true_and_binary_false(tmp_path):
    p = _write_airway(tmp_path, [("outlet", 0.0), ("mouth", 9.0)])
    assert is_named_multisolid_stl(str(p)) is True
    b = tmp_path / "bin.stl"
    b.write_bytes(b"solid binaryheader" + bytes(200))  # binary-ish, no endsolid lines
    assert is_named_multisolid_stl(str(b)) is False


def test_adapt_named_stl_writes_openprofile_and_seeds(tmp_path):
    p = _write_airway(tmp_path, [("outlet", 0.0), ("leftnose", 9.0), ("rightnose", 9.0)])
    out = tmp_path / "surface" / "_derived"
    op, seeds_path, unknown = adapt_named_stl(str(p), str(out))
    assert Path(op).exists() and Path(seeds_path).exists()
    data = json.loads(Path(seeds_path).read_text())
    labels = sorted(s["label"] for s in data["seeds"])
    assert labels == ["leftnose", "outlet", "rightnose"]
    assert unknown == []


def test_adapt_named_stl_requires_outlet_and_inlet(tmp_path):
    import pytest
    p = _write_airway(tmp_path, [("leftnose", 9.0), ("rightnose", 9.0)])  # no outlet
    with pytest.raises(ValueError):
        adapt_named_stl(str(p), str(tmp_path / "d"))
