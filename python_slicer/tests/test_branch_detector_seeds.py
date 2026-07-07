"""
Tests for the seeded anatomy-aware classifier in branch_detector.py.

assign_labels_from_seeds is a pure function (no VMTK dependency) so it must
import cleanly even outside the vmtk conda env.
"""

import numpy as np
from branch_detector import assign_labels_from_seeds


def test_assign_labels_from_seeds_maps_nearest():
    profiles = [
        {"id": 0, "barycenter": np.array([0.0, 0.0, 0.0])},   # near outlet
        {"id": 1, "barycenter": np.array([-5.0, 0.0, 40.0])},  # near leftnose
        {"id": 2, "barycenter": np.array([5.0, 0.0, 40.0])},   # near rightnose
    ]
    seeds = [
        {"label": "outlet",    "centroid": [0.0, 0.0, 0.0],  "normal": [0, 0, 1]},
        {"label": "leftnose",  "centroid": [-5.0, 0.0, 41.0], "normal": [0, 0, 1]},
        {"label": "rightnose", "centroid": [5.0, 0.0, 41.0],  "normal": [0, 0, 1]},
    ]
    source_ids, target_ids, labels = assign_labels_from_seeds(profiles, seeds)
    assert source_ids == [0]
    assert sorted(target_ids) == [1, 2]
    assert labels == {0: "Trachea", 1: "LeftNose", 2: "RightNose"}


def test_assign_labels_from_seeds_requires_one_source():
    import pytest
    profiles = [{"id": 0, "barycenter": np.array([0, 0, 0.0])}]
    seeds = [{"label": "leftnose", "centroid": [0, 0, 0.0], "normal": [0, 0, 1]}]
    with pytest.raises(ValueError):
        assign_labels_from_seeds(profiles, seeds)
