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


CAP_LABELS = {"leftnose", "rightnose", "mouth", "outlet"}

_ALIASES = {
    "trachea": "outlet",
    "out": "outlet",
    "outlet1": "outlet",
    "descending": "outlet",
    "leftnostril": "leftnose",
    "lnose": "leftnose",
    "left": "leftnose",
    "nl": "leftnose",
    "rightnostril": "rightnose",
    "rnose": "rightnose",
    "right": "rightnose",
    "nr": "rightnose",
}


def normalize_cap_name(name: str) -> str:
    """Lowercase, strip non-alphanumerics, resolve common synonyms."""
    key = "".join(ch for ch in name.strip().lower() if ch.isalnum())
    return _ALIASES.get(key, key)


@dataclass
class Cap:
    label: str          # normalized: leftnose/rightnose/mouth/outlet
    raw_name: str       # original solid name
    centroid: np.ndarray  # (3,) area-weighted face-center
    normal: np.ndarray    # (3,) unit, area-weighted


def _area_weighted_centroid(mesh: trimesh.Trimesh) -> np.ndarray:
    """Compute area-weighted centroid of mesh faces."""
    w = mesh.area_faces
    if w.sum() <= 0:
        return mesh.vertices.mean(axis=0)
    return np.average(mesh.triangles_center, axis=0, weights=w)


def _area_weighted_normal(mesh: trimesh.Trimesh) -> np.ndarray:
    """Compute unit area-weighted normal of mesh."""
    w = mesh.area_faces
    n = np.average(mesh.face_normals, axis=0, weights=w) if w.sum() > 0 \
        else mesh.face_normals.mean(axis=0)
    norm = np.linalg.norm(n)
    return n / norm if norm > 0 else np.array([0.0, 0.0, 1.0])


def classify_solids(solids: "dict[str, trimesh.Trimesh]") -> tuple[trimesh.Trimesh, list[Cap], list[str]]:
    """Split named solids into (wall_mesh, [Cap], [unknown_solid_names])."""
    wall_key = next((k for k in solids if normalize_cap_name(k) == "wall"), None)
    if wall_key is None:
        wall_key = max(solids, key=lambda k: len(solids[k].faces))
    wall = solids[wall_key]

    caps, unknown = [], []
    for k, m in solids.items():
        if k == wall_key:
            continue
        label = normalize_cap_name(k)
        if label not in CAP_LABELS:
            unknown.append(k)
            continue
        caps.append(Cap(label=label, raw_name=k,
                        centroid=_area_weighted_centroid(m),
                        normal=_area_weighted_normal(m)))
    return wall, caps, unknown


def is_named_multisolid_stl(path: str) -> bool:
    """True iff the STL is ASCII with more than one `endsolid` line."""
    try:
        text = Path(path).read_text(errors="replace")
    except (OSError, UnicodeError):
        return False
    n = sum(1 for ln in text.splitlines() if ln.strip().lower().startswith("endsolid"))
    return n > 1


def adapt_named_stl(named_stl: str, out_dir: str) -> tuple[Path, Path, list[str]]:
    """named multi-solid STL -> (openprofile.stl, anatomical_seeds.json, unknown)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    wall, caps, unknown = classify_solids(read_named_solids(named_stl))
    if wall is None or len(wall.faces) == 0:
        raise ValueError(f"No 'wall' solid (lumen surface) found in {named_stl}")
    labels = {c.label for c in caps}
    if "outlet" not in labels:
        raise ValueError("No 'outlet' cap found — cannot set centerline source.")
    if not (labels - {"outlet"}):
        raise ValueError("No inlet caps (leftnose/rightnose/mouth) found.")
    if unknown:
        log.warning("Unrecognized cap solids skipped: %s", unknown)

    open_profile = out_dir / "openprofile.stl"
    wall.export(str(open_profile))

    seeds = [{"label": c.label, "centroid": c.centroid.tolist(),
              "normal": c.normal.tolist()} for c in caps]
    seeds_path = out_dir / "anatomical_seeds.json"
    seeds_path.write_text(json.dumps({"seeds": seeds, "unknown_solids": unknown}, indent=2))
    log.info("Adapter: wrote %s (%d caps: %s)",
             open_profile.name, len(caps), sorted(labels))
    return open_profile, seeds_path, unknown
