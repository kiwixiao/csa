import logging
import numpy as np
import trimesh
from slice_bifurcation import slice_along_centerline

def _two_tubes():
    a = trimesh.creation.cylinder(radius=3.0, height=40.0, sections=48)   # axis z, at origin
    b = trimesh.creation.cylinder(radius=3.0, height=40.0, sections=48)
    b.apply_translation([30.0, 0.0, 0.0])                                 # residual "nose" tube
    return trimesh.util.concatenate([a, b])

def test_keep_only_containing_rejects_disjoint_loop():
    mesh = _two_tubes()
    centerline = np.array([[0.0, 0.0, z] for z in np.linspace(-15, 15, 11)])  # through tube A
    log = logging.getLogger("t"); log.addHandler(logging.NullHandler())

    kept_all = slice_along_centerline(mesh, centerline, None, log, keep_only_containing=False)
    kept_one = slice_along_centerline(mesh, centerline, None, log, keep_only_containing=True)

    assert max(len(loops) for _, _, loops in kept_all) >= 2     # both tubes cut
    assert kept_one and all(len(loops) == 1 for _, _, loops in kept_one)  # only tube A kept
