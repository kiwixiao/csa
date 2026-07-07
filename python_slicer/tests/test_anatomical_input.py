import numpy as np
import trimesh
from slicer.anatomical_input import read_named_solids, normalize_cap_name, CAP_LABELS

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
