import numpy as np
import trimesh
from slicer.anatomical_input import read_named_solids, normalize_cap_name, CAP_LABELS, classify_solids, Cap

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
