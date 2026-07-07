# Anatomical Named-STL Input Adapter — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a named multi-solid STL input format that extracts the lumen `wall` as an open-profile mesh plus anatomical cap seeds (`leftnose/rightnose/mouth/outlet`), so the existing bifurcation slicing pipeline runs on it unchanged, with a real Mouth region and adaptive handling of whichever caps are present.

**Architecture:** A new pure-Python adapter module (`slicer/anatomical_input.py`) parses the named solids, writes an open-profile STL + `anatomical_seeds.json`. `branch_detector.py` gains a seeded mode that uses those seeds instead of geometric guessing (and skips the aortic `Branch_N` rename). `pipeline_bifurcation.py` auto-detects the named STL, runs the adapter, relaxes the "both nostrils" guard, and slices the Mouth region.

**Tech Stack:** Python 3.10, trimesh, numpy, json (all already in the `mirtk`/`csa` env); VMTK only in the `vmtk` env for `branch_detector.py`. pytest for tests.

## Global Constraints

- The adapter and pipeline run in the `mirtk` (or `csa`) conda env; `branch_detector.py` runs in the `vmtk` env. **`branch_detector.py` must not import from `slicer.*`** (that package is not installed in the `vmtk` env) — the seeds JSON is the only cross-env artifact.
- No new dependencies. Use trimesh, numpy, json only.
- Tests run from `python_slicer/` (imports use `slicer.*` and top-level module paths). Run with `pytest tests/ -v`.
- Anatomical cap names are matched **case-insensitively** after alphanumeric normalization; the pipeline must adapt to whichever caps exist (never require a fixed set).
- Preserve existing behavior: closed-loop validity (`is_valid`), output CSV columns, per-region file layout. Do not touch the slicing/measurement math.
- **Commit messages must NOT contain a `Co-Authored-By` / Claude attribution line** (user preference).
- Follow existing code style (module-level functions, `print` in `branch_detector.py`, `logging` in pipeline/slice modules).

**Label mapping (canonical, used everywhere):**

| normalized cap label | pipeline label | role |
|---|---|---|
| `outlet` | `Trachea` (→ relabeled `DescendingAirway` by `assign_surface_partitions`) | centerline source |
| `leftnose` | `LeftNose` | target |
| `rightnose` | `RightNose` | target |
| `mouth` | `Mouth` | target |

---

### Task 1: Multi-solid ASCII STL reader

**Files:**
- Create: `python_slicer/slicer/anatomical_input.py`
- Test: `python_slicer/tests/test_anatomical_input.py`

**Interfaces:**
- Produces: `read_named_solids(path: str) -> dict[str, trimesh.Trimesh]` — one mesh per `solid <name> … endsolid` block, keyed by the raw solid name.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_anatomical_input.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd python_slicer && pytest tests/test_anatomical_input.py::test_read_named_solids_splits_blocks -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'slicer.anatomical_input'`

- [ ] **Step 3: Write minimal implementation**

```python
# slicer/anatomical_input.py
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


def read_named_solids(path):
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd python_slicer && pytest tests/test_anatomical_input.py::test_read_named_solids_splits_blocks -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python_slicer/slicer/anatomical_input.py python_slicer/tests/test_anatomical_input.py
git commit -m "Add multi-solid ASCII STL reader for anatomical input"
```

---

### Task 2: Cap-name normalization + alias table

**Files:**
- Modify: `python_slicer/slicer/anatomical_input.py`
- Test: `python_slicer/tests/test_anatomical_input.py`

**Interfaces:**
- Produces: `normalize_cap_name(name: str) -> str` (lowercased, alphanumeric-only, alias-resolved); module constant `CAP_LABELS = {"leftnose","rightnose","mouth","outlet"}`.

- [ ] **Step 1: Write the failing test**

```python
from slicer.anatomical_input import normalize_cap_name, CAP_LABELS

def test_normalize_cap_name_case_and_aliases():
    assert normalize_cap_name("LeftNose") == "leftnose"
    assert normalize_cap_name("  RIGHT_NOSE ") == "rightnose"
    assert normalize_cap_name("outlet_1") == "outlet"
    assert normalize_cap_name("Trachea") == "outlet"
    assert normalize_cap_name("Mouth") == "mouth"
    assert normalize_cap_name("wall") == "wall"
    assert {"leftnose", "rightnose", "mouth", "outlet"} == CAP_LABELS
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd python_slicer && pytest tests/test_anatomical_input.py::test_normalize_cap_name_case_and_aliases -v`
Expected: FAIL with `ImportError: cannot import name 'normalize_cap_name'`

- [ ] **Step 3: Write minimal implementation**

Add to `slicer/anatomical_input.py`:

```python
CAP_LABELS = {"leftnose", "rightnose", "mouth", "outlet"}

_ALIASES = {
    "trachea": "outlet", "out": "outlet", "outlet1": "outlet", "descending": "outlet",
    "leftnostril": "leftnose", "lnose": "leftnose", "left": "leftnose", "nl": "leftnose",
    "rightnostril": "rightnose", "rnose": "rightnose", "right": "rightnose", "nr": "rightnose",
}


def normalize_cap_name(name):
    """Lowercase, strip non-alphanumerics, resolve common synonyms."""
    key = "".join(ch for ch in name.strip().lower() if ch.isalnum())
    return _ALIASES.get(key, key)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd python_slicer && pytest tests/test_anatomical_input.py::test_normalize_cap_name_case_and_aliases -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python_slicer/slicer/anatomical_input.py python_slicer/tests/test_anatomical_input.py
git commit -m "Add case-insensitive cap-name normalization with alias table"
```

---

### Task 3: Cap classification + seed extraction

**Files:**
- Modify: `python_slicer/slicer/anatomical_input.py`
- Test: `python_slicer/tests/test_anatomical_input.py`

**Interfaces:**
- Consumes: `read_named_solids`, `normalize_cap_name`, `CAP_LABELS`.
- Produces:
  - `@dataclass Cap: label:str, raw_name:str, centroid:np.ndarray(3,), normal:np.ndarray(3,)`
  - `classify_solids(solids: dict) -> tuple[trimesh.Trimesh, list[Cap], list[str]]` → `(wall_mesh, caps, unknown_solid_names)`.

- [ ] **Step 1: Write the failing test**

```python
from slicer.anatomical_input import classify_solids, Cap, read_named_solids

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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd python_slicer && pytest tests/test_anatomical_input.py::test_classify_solids_identifies_wall_and_caps -v`
Expected: FAIL with `ImportError: cannot import name 'classify_solids'`

- [ ] **Step 3: Write minimal implementation**

Add to `slicer/anatomical_input.py`:

```python
@dataclass
class Cap:
    label: str          # normalized: leftnose/rightnose/mouth/outlet
    raw_name: str       # original solid name
    centroid: np.ndarray  # (3,) area-weighted face-center
    normal: np.ndarray    # (3,) unit, area-weighted


def _area_weighted_centroid(mesh):
    w = mesh.area_faces
    if w.sum() <= 0:
        return mesh.vertices.mean(axis=0)
    return np.average(mesh.triangles_center, axis=0, weights=w)


def _area_weighted_normal(mesh):
    w = mesh.area_faces
    n = np.average(mesh.face_normals, axis=0, weights=w) if w.sum() > 0 \
        else mesh.face_normals.mean(axis=0)
    norm = np.linalg.norm(n)
    return n / norm if norm > 0 else np.array([0.0, 0.0, 1.0])


def classify_solids(solids):
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd python_slicer && pytest tests/test_anatomical_input.py::test_classify_solids_identifies_wall_and_caps -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python_slicer/slicer/anatomical_input.py python_slicer/tests/test_anatomical_input.py
git commit -m "Classify named solids into wall + anatomical cap seeds"
```

---

### Task 4: Adapter orchestrator + multi-solid detector

**Files:**
- Modify: `python_slicer/slicer/anatomical_input.py`
- Test: `python_slicer/tests/test_anatomical_input.py`

**Interfaces:**
- Consumes: `read_named_solids`, `classify_solids`.
- Produces:
  - `is_named_multisolid_stl(path: str) -> bool` — True iff the file has >1 `endsolid` line (ASCII multi-solid; binary STL → False).
  - `adapt_named_stl(named_stl: str, out_dir: str) -> tuple[Path, Path, list[str]]` → `(openprofile_stl_path, seeds_json_path, unknown_names)`. Writes `openprofile.stl` (the wall) and `anatomical_seeds.json` (`{"seeds":[{"label","centroid":[x,y,z],"normal":[x,y,z]}], "unknown_solids":[...]}`). Raises `ValueError` if no wall, no `outlet`, or no inlet.

- [ ] **Step 1: Write the failing test**

```python
import json
from pathlib import Path
from slicer.anatomical_input import is_named_multisolid_stl, adapt_named_stl

def _write_airway(tmp_path, names_z):
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd python_slicer && pytest tests/test_anatomical_input.py -k "multisolid or adapt" -v`
Expected: FAIL with `ImportError: cannot import name 'is_named_multisolid_stl'`

- [ ] **Step 3: Write minimal implementation**

Add to `slicer/anatomical_input.py`:

```python
def is_named_multisolid_stl(path):
    """True iff the STL is ASCII with more than one `endsolid` line."""
    try:
        text = Path(path).read_text(errors="replace")
    except (OSError, UnicodeError):
        return False
    n = sum(1 for ln in text.splitlines() if ln.strip().lower().startswith("endsolid"))
    return n > 1


def adapt_named_stl(named_stl, out_dir):
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd python_slicer && pytest tests/test_anatomical_input.py -v`
Expected: PASS (all anatomical_input tests)

- [ ] **Step 5: Commit**

```bash
git add python_slicer/slicer/anatomical_input.py python_slicer/tests/test_anatomical_input.py
git commit -m "Add adapter orchestrator + named-multisolid detector"
```

---

### Task 5: Seeded mode in branch_detector + step0 passthrough

**Files:**
- Modify: `python_slicer/branch_detector.py` (add `assign_labels_from_seeds`; thread `seeds` through `AirwayBranchDetector.__init__`/`run`; gate the `Branch_N` rename; add `--seeds` CLI)
- Modify: `python_slicer/step0_branch_detect.py` (add `--seeds`, pass through to `branch_detector.py`)
- Test: `python_slicer/tests/test_branch_detector_seeds.py`

**Interfaces:**
- Produces (importable, pure, no VMTK): `assign_labels_from_seeds(profiles: list[dict], seeds: list[dict]) -> tuple[list[int], list[int], dict[int,str]]` — same return shape as `classify_open_profiles`. Each profile (`{"id", "barycenter": np.ndarray}`) is labeled by its nearest seed centroid; `outlet` → source (mapped to `"Trachea"`), others → targets (`LeftNose/RightNose/Mouth`). Raises `ValueError` unless exactly one source.
- Consumed by: `pipeline_bifurcation.py` (Task 6) via `step0_branch_detect.py --seeds`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_branch_detector_seeds.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd python_slicer && pytest tests/test_branch_detector_seeds.py -v`
Expected: FAIL with `ImportError: cannot import name 'assign_labels_from_seeds'`

(Note: `branch_detector.py` imports `vtk`/`vmtk` at module top. If pytest runs outside the `vmtk` env this import fails. Guard those imports — see Step 3 — so the pure function is importable anywhere.)

- [ ] **Step 3a: Make VMTK imports lazy so the module imports without vmtk**

In `branch_detector.py`, near the top where `import vtk` / `from vmtk import vmtkscripts` live, wrap them:

```python
try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    from vmtk import vmtkscripts
    _HAVE_VMTK = True
except ImportError:  # allows importing pure helpers (e.g. assign_labels_from_seeds) without vmtk
    _HAVE_VMTK = False
```

(If the existing imports are spread across lines, consolidate them into this guarded block.)

- [ ] **Step 3b: Add the pure seeded classifier**

Add near `classify_open_profiles` in `branch_detector.py`:

```python
import numpy as np  # already imported at top; keep single import

CAP_TO_PIPELINE_LABEL = {
    "outlet": "Trachea",       # relabeled DescendingAirway downstream
    "leftnose": "LeftNose",
    "rightnose": "RightNose",
    "mouth": "Mouth",
}


def assign_labels_from_seeds(profiles, seeds):
    """Label detected open profiles from anatomical seed centroids.

    profiles: [{"id", "barycenter": (3,) array}], seeds: [{"label","centroid",...}].
    Each profile takes the label of its nearest seed. Returns
    (source_ids, target_ids, labels) — same shape as classify_open_profiles.
    """
    seed_pts = np.asarray([s["centroid"] for s in seeds], dtype=float)
    seed_lbl = [s["label"] for s in seeds]

    labels, source_ids, target_ids = {}, [], []
    for p in profiles:
        d = np.linalg.norm(seed_pts - np.asarray(p["barycenter"], dtype=float), axis=1)
        cap = seed_lbl[int(np.argmin(d))]
        labels[p["id"]] = CAP_TO_PIPELINE_LABEL.get(cap, cap)
        (source_ids if cap == "outlet" else target_ids).append(p["id"])

    if len(source_ids) != 1:
        raise ValueError(
            f"Expected exactly one 'outlet' profile, got {len(source_ids)} "
            f"(labels={list(labels.values())})")
    return source_ids, target_ids, labels
```

- [ ] **Step 3c: Thread seeds into the detector and gate the rename**

In `AirwayBranchDetector.__init__` add a `seeds` parameter:

```python
    def __init__(self, stl_path: str, resampling_step: float = 2.0, seeds=None):
        self.stl_path = Path(stl_path)
        self.resampling_step = resampling_step
        self.seeds = seeds            # None -> geometric guessing; list -> seeded
        if not self.stl_path.exists():
            raise FileNotFoundError(f"STL not found: {self.stl_path}")
```

In `run()`, replace the classify call (currently line ~607):

```python
        # Step 2: Detect and classify open profiles
        profiles = detect_open_profiles(surface)
        if self.seeds is not None:
            source_ids, target_ids, profile_labels = assign_labels_from_seeds(
                profiles, self.seeds)
        else:
            source_ids, target_ids, profile_labels = classify_open_profiles(profiles)
```

Wrap the reorder/rename block (currently lines ~671-697, from the `# Reorder branches by arc-length…` comment through the three `rename_map.get(...)` reassignments) in a guard so it only runs for the geometric/aortic path. Concretely: insert `if self.seeds is None:` immediately above the existing line 671 comment, indent the entire existing block (lines 671-697) one level to the right (do not change any of its logic), and add the trailing comment. The result is:

```python
        # >>> existing lines 671-697 kept verbatim, indented under this guard <<<
        if self.seeds is None:
            # Reorder branches by arc-length position along the longest centerline
            longest_label = max(full_paths, key=lambda k: len(full_paths[k]))
            longest_pts = full_paths[longest_label]
            from scipy.spatial import cKDTree as _cKDTree
            main_tree = _cKDTree(longest_pts)
            branch_order = []
            main_diffs = np.linalg.norm(np.diff(longest_pts, axis=0), axis=1)
            main_arc = np.concatenate([[0], np.cumsum(main_diffs)])
            for label, pts in full_paths.items():
                _, idx = main_tree.query(pts[-1])
                branch_order.append((label, float(main_arc[idx])))
            branch_order.sort(key=lambda x: x[1])
            rename_map = {}
            for i, (old_label, arc) in enumerate(branch_order):
                new_label = f"Branch_{i+1}"
                rename_map[old_label] = new_label
                print(f"  {new_label} (was {old_label}): arc={arc:.1f}mm")
            full_paths = {rename_map.get(k, k): v for k, v in full_paths.items()}
            surface_partitions = {rename_map.get(k, k): v for k, v in surface_partitions.items()}
            target_labels = [rename_map.get(t, t) for t in target_labels]
        # When seeds are given, anatomical labels (LeftNose/RightNose/Mouth/
        # DescendingAirway) pass straight through — no Branch_N renaming.
```

- [ ] **Step 3d: Add `--seeds` to the branch_detector CLI**

Replace the CLI tail (lines ~755-767) so it loads seeds and passes them in:

```python
if __name__ == "__main__":
    import argparse, json
    parser = argparse.ArgumentParser(
        description="Detect airway branches from a full airway STL using VMTK")
    parser.add_argument("stl_path", help="Path to full airway STL file")
    parser.add_argument("--output-dir", default="./branch_output")
    parser.add_argument("--subject-id", default="")
    parser.add_argument("--resampling-step", type=float, default=2.0)
    parser.add_argument("--seeds", default="",
                        help="anatomical_seeds.json (enables seeded, anatomy-aware mode)")
    args = parser.parse_args()

    seeds = None
    if args.seeds:
        seeds = json.loads(Path(args.seeds).read_text())["seeds"]

    detector = AirwayBranchDetector(args.stl_path, args.resampling_step, seeds=seeds)
    detector.run(args.output_dir, args.subject_id)
```

- [ ] **Step 3e: Pass `--seeds` through step0**

In `step0_branch_detect.py`: add the argument (after `--verbose`, line ~69):

```python
    parser.add_argument("--seeds", default="",
                        help="Path to anatomical_seeds.json for seeded mode")
```

and extend the command (after the `--resampling-step` block, line ~98):

```python
    if args.seeds:
        cmd.extend(["--seeds", str(Path(args.seeds).resolve())])
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd python_slicer && pytest tests/test_branch_detector_seeds.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python_slicer/branch_detector.py python_slicer/step0_branch_detect.py python_slicer/tests/test_branch_detector_seeds.py
git commit -m "Add seeded anatomy-aware mode to branch detector; bypass Branch_N rename"
```

---

### Task 6: Pipeline auto-detects named STL + adaptive branch guard

**Files:**
- Modify: `python_slicer/pipeline_bifurcation.py` (`find_surface_stl` / Step 0 / guard / thread mouth centerline)
- Test: `python_slicer/tests/test_pipeline_adaptive.py`

**Interfaces:**
- Consumes: `slicer.anatomical_input.is_named_multisolid_stl`, `adapt_named_stl`; `discover_branches` (existing).
- Produces: `resolve_surface_and_seeds(subject_dir) -> tuple[Path, Path|None]` → `(open_profile_stl, seeds_json_or_None)`; `select_regions(branches) -> dict` with keys `has_left/has_right/has_mouth/inlets` and validation that raises `SystemExit`-free `ValueError` when no outlet or no inlet.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_pipeline_adaptive.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd python_slicer && pytest tests/test_pipeline_adaptive.py -v`
Expected: FAIL with `ImportError: cannot import name 'select_regions'`

- [ ] **Step 3a: Add `select_regions` and `resolve_surface_and_seeds`**

Add to `pipeline_bifurcation.py` (after `find_surface_stl`, line ~90):

```python
from slicer.anatomical_input import is_named_multisolid_stl, adapt_named_stl


def resolve_surface_and_seeds(subject_dir):
    """Return (open_profile_stl, seeds_json | None).

    If the surface STL is a named multi-solid file, run the adapter to derive an
    open-profile mesh + anatomical seeds under surface/_derived/. Otherwise pass
    the STL through unchanged (legacy geometric detection, no seeds).
    """
    surface_stl = find_surface_stl(subject_dir)
    if is_named_multisolid_stl(str(surface_stl)):
        derived = Path(subject_dir) / "surface" / "_derived"
        open_profile, seeds_json, unknown = adapt_named_stl(str(surface_stl), str(derived))
        if unknown:
            log.warning(f"  Adapter skipped unrecognized cap solids: {unknown}")
        log.info(f"  Named multi-solid STL detected -> {open_profile.name} + seeds")
        return open_profile, seeds_json
    return surface_stl, None


def select_regions(branches):
    """Validate + report which anatomical regions are present."""
    if "DescendingAirway" not in branches:
        raise ValueError("No DescendingAirway (outlet) branch found.")
    inlets = [r for r in ("LeftNose", "RightNose", "Mouth")
              if r in branches and branches[r].get("centerline")]
    if not inlets:
        raise ValueError("No inlet branch (LeftNose/RightNose/Mouth) found.")
    return {
        "inlets": inlets,
        "has_left": "LeftNose" in inlets,
        "has_right": "RightNose" in inlets,
        "has_mouth": ("Mouth" in branches and branches["Mouth"].get("mesh") is not None),
    }
```

- [ ] **Step 3b: Wire Step 0 to use resolver + seeds**

Replace the Step-0 block (lines ~124-153) so it feeds the resolved surface + seeds and uses the adaptive guard:

```python
    # ── Step 0: Branch detection (frame 0) ───────────────────
    logger.info(f"\n--- Step 0: Branch detection ---")
    branches_dir = subject_dir / "branches"
    surface_stl, seeds_json = resolve_surface_and_seeds(subject_dir)

    if not branches_dir.exists():
        logger.info("  Running branch detection on frame 0 (vmtk env)...")
        import subprocess
        script_dir = Path(__file__).parent
        cmd = [
            "conda", "run", "-n", "vmtk",
            "python", str(script_dir / "step0_branch_detect.py"),
            str(surface_stl),
            "--output-dir", str(branches_dir),
            "--subject-id", subject_id,
        ]
        if seeds_json is not None:
            cmd += ["--seeds", str(seeds_json)]
        result = subprocess.run(cmd, capture_output=False)
        if result.returncode != 0:
            logger.error("Branch detection failed!")
            sys.exit(1)
    else:
        logger.info(f"  Branches already exist: {branches_dir}")

    branches = discover_branches(branches_dir, subject_id)
    try:
        regions = select_regions(branches)
    except ValueError as e:
        logger.error(str(e))
        sys.exit(1)
    has_left, has_right, has_mouth = regions["has_left"], regions["has_right"], regions["has_mouth"]
    logger.info(f"  Regions present: {regions['inlets']} (+ DescendingAirway)")
    if not (has_left and has_right):
        logger.warning("  Fewer than two nostrils present; septum split will be "
                       "skipped for missing sides (see mouth-only handling).")
```

Also update `find_surface_stl` to ignore the derived dir so it never re-picks `_derived/openprofile.stl` (line ~85):

```python
    stl_files = [p for p in surface_dir.glob("*.stl") if p.parent.name != "_derived"]
```

(Glob is non-recursive, so this is a no-op today but guards against a future move; keep it for clarity.)

- [ ] **Step 3c: Thread the mouth centerline into refine_result**

After `refine_result = refine_nasal_partition(...)` (line ~166-168), attach the mouth centerline + mesh path so `run_single_frame` (Task 7) can slice it:

```python
    if has_mouth:
        refine_result["mouth_cl"] = read_vtk_centerline(str(branches["Mouth"]["centerline"]))
        refine_result["mouth_mesh_path"] = str(branches["Mouth"]["mesh"])
    else:
        refine_result["mouth_cl"] = None
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd python_slicer && pytest tests/test_pipeline_adaptive.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python_slicer/pipeline_bifurcation.py python_slicer/tests/test_pipeline_adaptive.py
git commit -m "Auto-detect named STL, adaptive region selection, seed passthrough"
```

---

### Task 7: Mouth region slicing in run_single_frame

**Files:**
- Modify: `python_slicer/slice_bifurcation.py` (`run_single_frame` — add mouth block; add `"Mouth"` to output-dir creation and the Step-4 region loop in `pipeline_bifurcation.py`)
- Test: `python_slicer/tests/test_mouth_region.py`

**Interfaces:**
- Consumes: `slice_along_centerline` (existing), `extract_submesh` (existing), `refine_result["mouth_cl"]`, mouth submesh (face label `0`).
- Produces: mouth rows appended to `run_single_frame`'s returned `all_rows` with `region="Mouth"`, `centerline_source="mouth"`, same column set as other regions.

- [ ] **Step 1: Write the failing test** (unit-test the row-builder helper in isolation)

Add a small pure helper so mouth logic is testable without a full mesh/VMTK run:

```python
# tests/test_mouth_region.py
import numpy as np
from slice_bifurcation import build_region_row

def test_build_region_row_mouth():
    row = build_region_row(
        plane_idx=3, region="Mouth", arc_len=12.0, branch_arc=12.0,
        area=4.0, perimeter=8.0,
        centroid=np.array([1.0, 2.0, 3.0]), centerline_source="mouth")
    assert row["region"] == "Mouth"
    assert row["area_mm2"] == 4.0
    assert row["hydraulic_diameter_mm"] == 2.0   # 4*area/perimeter = 16/8
    assert row["centroid_z"] == 3.0
    assert row["is_valid"] is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd python_slicer && pytest tests/test_mouth_region.py -v`
Expected: FAIL with `ImportError: cannot import name 'build_region_row'`

- [ ] **Step 3a: Extract the shared row-builder**

Add to `slice_bifurcation.py` (near the top of the slicing helpers, after `slice_along_midline`):

```python
def build_region_row(plane_idx, region, arc_len, branch_arc, area, perimeter,
                     centroid, centerline_source):
    """Assemble one CSA output row (shared by all regions)."""
    hyd = DiameterCalculator.compute_hydraulic_diameter(area, perimeter)
    equiv = DiameterCalculator.compute_equivalent_diameter(area)
    return {
        "plane_index": plane_idx, "region": region,
        "arc_length_mm": arc_len, "arc_length_branch_mm": branch_arc,
        "area_mm2": area, "perimeter_mm": perimeter,
        "hydraulic_diameter_mm": hyd, "equivalent_diameter_mm": equiv,
        "centroid_x": centroid[0], "centroid_y": centroid[1], "centroid_z": centroid[2],
        "centerline_source": centerline_source, "is_valid": True,
    }
```

- [ ] **Step 3b: Add the mouth slice block to `run_single_frame`**

In `run_single_frame`, add `"Mouth_Planes"` to the output-dir list (line ~467):

```python
    for d in ["DescendingAirway_Planes", "LeftNose_Planes", "RightNose_Planes", "Mouth_Planes"]:
        (output_dir / d).mkdir(parents=True, exist_ok=True)
```

After the right-nose block (before `log.info(f"  Frame {name}: ...")`, line ~573), add:

```python
    # ── Slice mouth (single passage, containment-filtered) ──
    mouth_count = 0
    mouth_cl = refine_result.get("mouth_cl")
    if mouth_cl is not None and full_face_labels is not None:
        mouth_mesh = extract_submesh(full_mesh, full_face_labels == 0)
        if len(mouth_mesh.faces) > 0:
            mouth_results = slice_along_centerline(
                mouth_mesh, mouth_cl, mouth_mesh, log, label_prefix="Mouth ")
            mouth_positions = np.array([pos for _, pos, _ in mouth_results])
            mouth_arc = compute_arc_length(mouth_positions) if len(mouth_positions) else np.array([])
            for k, (plane_idx, pos, loops) in enumerate(mouth_results):
                merged = merge_loops(loops)
                if merged is None:
                    continue
                arc_len = float(mouth_arc[k]) if k < len(mouth_arc) else 0.0
                stl_out = output_dir / "Mouth_Planes" / f"{name}-Mouth-{mouth_count:03d}.stl"
                write_cross_section_stl(str(stl_out), merged.vertices, merged.faces)
                mouth_count += 1
                all_rows.append(build_region_row(
                    plane_idx, "Mouth", arc_len, arc_len,
                    merged.area, merged.perimeter, merged.centroid, "mouth"))
```

Update the frame summary log line to include mouth:

```python
    log.info(f"  Frame {name}: Desc={desc_count}, Left={left_count}, "
             f"Right={right_count}, Mouth={mouth_count}")
```

- [ ] **Step 3c: Include Mouth in the pipeline's mutual-plane + All_Planes loops**

In `pipeline_bifurcation.py` line ~276 change the region list:

```python
    for region in ["DescendingAirway", "LeftNose", "RightNose", "Mouth"]:
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd python_slicer && pytest tests/test_mouth_region.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python_slicer/slice_bifurcation.py python_slicer/pipeline_bifurcation.py python_slicer/tests/test_mouth_region.py
git commit -m "Implement Mouth region CSA (single-passage slicing)"
```

---

### Task 8: `csa-uncap` standalone command + docs

**Files:**
- Create: `python_slicer/uncap_cli.py`
- Modify: `python_slicer/setup.py` (or `pyproject.toml` — whichever defines `console_scripts`)
- Modify: `python_slicer/tests/test_anatomical_input.py` (CLI smoke test)
- Modify: `CLAUDE.md` (document the new input format + command)

**Interfaces:**
- Consumes: `adapt_named_stl`.
- Produces: console entry `csa-uncap <named.stl> [--out DIR]` printing the derived open-profile + seeds paths.

- [ ] **Step 1: Write the failing test**

```python
def test_uncap_cli_smoke(tmp_path, monkeypatch, capsys):
    import uncap_cli
    p = _write_airway(tmp_path, [("outlet", 0.0), ("mouth", 9.0)])  # helper from Task 4 tests
    out = tmp_path / "out"
    monkeypatch.setattr("sys.argv", ["csa-uncap", str(p), "--out", str(out)])
    uncap_cli.main()
    captured = capsys.readouterr().out
    assert "openprofile.stl" in captured
    assert (out / "anatomical_seeds.json").exists()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd python_slicer && pytest tests/test_anatomical_input.py::test_uncap_cli_smoke -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'uncap_cli'`

- [ ] **Step 3a: Implement the CLI**

```python
# uncap_cli.py
"""csa-uncap: convert a named multi-solid STL into open-profile + seeds."""
import argparse
from pathlib import Path
from slicer.anatomical_input import adapt_named_stl


def main():
    ap = argparse.ArgumentParser(description="Extract open-profile mesh + anatomical "
                                             "seeds from a named multi-solid STL")
    ap.add_argument("named_stl", help="ASCII multi-solid STL (solid leftnose/... /wall)")
    ap.add_argument("--out", default="./_derived", help="Output dir (default ./_derived)")
    args = ap.parse_args()

    op, seeds, unknown = adapt_named_stl(args.named_stl, args.out)
    print(f"open-profile: {op}")
    print(f"seeds:        {seeds}")
    if unknown:
        print(f"WARNING: skipped unrecognized cap solids: {unknown}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3b: Register the console script**

In `setup.py`, add to the existing `entry_points["console_scripts"]` list an entry:

```python
        "csa-uncap = uncap_cli:main",
```

(Match the existing style of the other `csa-*` entries; if the project uses `pyproject.toml [project.scripts]`, add `csa-uncap = "uncap_cli:main"` there instead. Re-run `pip install -e .` after editing.)

- [ ] **Step 3c: Document in CLAUDE.md**

Under "## CLI Commands", add:

```markdown
# Convert a named multi-solid STL (solid leftnose/rightnose/mouth/outlet/wall)
# into an open-profile mesh + anatomical seeds (standalone / debugging).
csa-uncap /path/to/named_airway.stl --out /path/to/subject/surface/_derived
```

Under "### Three Pipelines" (or a new "Input Formats" note), record: a named
multi-solid STL dropped into `subject/surface/` is auto-detected by
`csa-bifurcation`, which runs the adapter and feeds anatomical seeds to branch
detection (case-insensitive names; adapts to whichever of leftnose/rightnose/
mouth/outlet are present).

- [ ] **Step 4: Run test to verify it passes**

Run: `cd python_slicer && pip install -e . && pytest tests/test_anatomical_input.py::test_uncap_cli_smoke -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python_slicer/uncap_cli.py python_slicer/setup.py python_slicer/tests/test_anatomical_input.py CLAUDE.md
git commit -m "Add csa-uncap command and document named-STL input format"
```

---

### Task 9: Full regression + end-to-end verification

**Files:**
- Test: run existing suite; manual end-to-end note.

- [ ] **Step 1: Run the whole test suite**

Run: `cd python_slicer && pytest tests/ -v`
Expected: PASS (new tests + existing `test_basic.py`; `test_branch_detector.py` may already be stale/skipped — do not regress it further).

- [ ] **Step 2: End-to-end smoke (requires vmtk env + a real named STL)**

On a subject whose `surface/` holds a named multi-solid STL:

```bash
conda activate mirtk
csa-bifurcation /path/to/subject/ --subject-id TEST --frames 0
```

Expected: `surface/_derived/openprofile.stl` + `anatomical_seeds.json` created; `branches/` contains anatomical folders (`LeftNose/RightNose/Mouth/DescendingAirway`, whichever present) — NOT `Branch_1/2/3`; per-region CSV includes a `Mouth` region when a mouth cap exists; `is_valid` column present and open slices excluded.

- [ ] **Step 3: Commit any fixes**

```bash
git add -A && git commit -m "Fixes from end-to-end verification of anatomical-STL input"
```

---

## Notes / scope boundaries

- **Mouth-only (no nostrils) and single-nostril configs:** the adapter, seeded detection, and `select_regions` all accept these, and `branch_detector` will produce a `Mouth`+`DescendingAirway` (or single-nostril) branch set. However, `refine_nasal_partition` (the septum split) assumes **both** nostrils. For a mouth-only airway there is no septum and no L/R split — the Mouth branch slices fine, but the descending-airway slicing currently derives its centerline from the nasal midline. Making the descending region fall back to the outlet→mouth centerline when <2 nostrils are present is a **follow-up** (a small branch in `run_single_frame` Step "Slice descending"), flagged here so it is not mistaken for covered. This plan fully delivers the both-nostril (+ optional mouth) case end-to-end.
- `branch_detector.py` legacy geometric path (open mesh, no seeds) keeps its current `Branch_N` behavior — unchanged and out of scope; the seeded path is the supported route for the new input.
