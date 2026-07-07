import logging
import numpy as np
import trimesh
from slice_bifurcation import run_single_branch_frame

def test_single_branch_frame_mouth_only(tmp_path):
    a = trimesh.creation.cylinder(radius=3.0, height=40.0, sections=48)
    b = trimesh.creation.cylinder(radius=3.0, height=40.0, sections=48)
    b.apply_translation([30.0, 0.0, 0.0])                     # residual nose tube
    mesh = trimesh.util.concatenate([a, b])
    stl = tmp_path / "wall.stl"; mesh.export(str(stl))
    cl = np.array([[0.0, 0.0, z] for z in np.linspace(-15, 15, 11)])

    rows = run_single_branch_frame(str(stl), str(tmp_path / "out"), "F0", cl, "Mouth",
                                   logging.getLogger("t"))
    assert rows and all(r["region"] == "Mouth" for r in rows)
    # residual tube B is excluded, so area ~ pi*3^2 (~28), never doubled
    assert all(r["area_mm2"] < 40.0 for r in rows)
