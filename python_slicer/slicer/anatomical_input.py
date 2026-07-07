"""Adapter for named multi-solid STL input.

Reads an ASCII STL with multiple `solid <name> … endsolid` blocks (as OpenFOAM
constant/triSurface writes), where each opening is a named cap and the lumen is
the `wall` solid. Produces an open-profile mesh + anatomical seed points that
feed the existing branch-detection / slicing pipeline.
"""
import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import trimesh

log = logging.getLogger(__name__)


def read_named_solids(path: str) -> "dict[str, trimesh.Trimesh]":
    """Parse an ASCII multi-solid STL into {solid_name: trimesh.Trimesh}."""
    text = Path(path).read_text(errors="replace")
    solids = {}
    name = None
    verts = []
    for line in text.splitlines():
        s = line.strip()
        low = s.lower()
        if low.startswith("solid"):
            name = s[5:].strip() or "unnamed"
            verts = []
        elif low.startswith("vertex"):
            parts = s.split()
            verts.append([float(parts[1]), float(parts[2]), float(parts[3])])
        elif low.startswith("endsolid"):
            if name is not None and len(verts) >= 3:
                v = np.asarray(verts, dtype=float)
                faces = np.arange(len(v)).reshape(-1, 3)
                solids[name] = trimesh.Trimesh(vertices=v, faces=faces, process=True)
            name = None
            verts = []
    return solids
