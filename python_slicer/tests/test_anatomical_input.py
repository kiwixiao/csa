import numpy as np
import trimesh
from slicer.anatomical_input import read_named_solids

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
