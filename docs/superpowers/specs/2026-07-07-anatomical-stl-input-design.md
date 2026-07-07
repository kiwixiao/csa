# Anatomical Named-STL Input Adapter — Design Spec

**Date:** 2026-07-07
**Status:** Draft for review
**Author:** brainstorming session (qiwei.xiao@cchmc.org + Claude)

## 1. Goal

Add support for a **new input format**: a watertight/capped airway solid exported as a
**named multi-solid STL**, where each opening is a separately named cap
(`leftnose`, `rightnose`, `mouth`, `outlet`) and the lumen is the `wall` solid.

The capped solid exists only to produce an **anatomy-aware open-profile mesh + seeds**
that feed the **existing bifurcation slicing pipeline unchanged**. No slicing, centerline,
septum-split, or CSA math changes. Two small pieces ride along because this input requires
them: a real **Mouth** region, and **adaptive** handling of whichever caps are present.

This is a front-end format adapter, not a new pipeline.

## 2. Input format

- **ASCII multi-solid STL** with `solid <name> … endsolid <name>` blocks (the format OpenFOAM
  `constant/triSurface` already writes, e.g. `ENT001_clean.stl` → `solid inlet`, `solid outlet_1`,
  `solid wall`). Binary STL cannot carry per-region names, so named input is ASCII by definition.
- Exactly one `wall` solid (the lumen surface, no caps) plus one solid per opening.
- Opening names are drawn from `{leftnose, rightnose, mouth, outlet}`, matched **case-insensitively**
  and after normalization (lowercase, strip spaces/underscores/hyphens). A small alias table maps
  common synonyms (e.g. `trachea`→`outlet`; `left_nose`/`lnose`→`leftnose`).
- **Not every opening is present.** Supported configurations (always exactly one `outlet`):
  - `outlet` + `mouth` (mouth-only, no nose)
  - `outlet` + `leftnose` + `rightnose` (nose-only, no mouth)
  - `outlet` + `leftnose` + `rightnose` + `mouth` (all)
  - `outlet` + a single nostril + optionally `mouth` (asymmetric; supported, warned)

## 3. Non-goals (out of scope)

- Detecting caps geometrically on a single-solid capped STL (coplanar-face clustering). Not needed:
  named solids already separate wall from caps. A future single-solid path can be added later.
- Any change to centerline generation (still VMTK), septum refinement, slicing, loop filtering,
  or CSA measurement.
- Multi-frame / FFD / motion wiring. The adapter produces a static open-profile + seeds; whatever
  the downstream pipeline does after that is unchanged.

## 4. Architecture

```
named multi-solid STL
        │
        ▼
[ NEW: AnatomicalCapAdapter ]              ← the only genuinely new module
   • parse named solids  → {wall, caps}
   • wall solid          → open-profile mesh   (caps simply not included)
   • each cap            → seed {label, centroid, normal}
        │
        ├─► open-profile STL  (surface for existing pipeline)
        └─► anatomical_seeds.json  (sidecar: [{label, centroid, normal}])
        │
        ▼
[ branch_detector.py — NEW seeded mode ]
   • detect_open_profiles(wall) finds the holes (unchanged)
   • IF seeds present: match each hole barycenter → nearest seed centroid → anatomical label
     (SKIP classify_open_profiles guessing; SKIP the Branch_N rename)
   • extract_centerlines: source = outlet seed, targets = inlet seeds  (VMTK "backward", unchanged)
   • partition surface, save branches/  with anatomical folder names
        │
        ▼
[ pipeline_bifurcation.py — adaptive region set ]
   • relax "both nostrils required" guard → "≥1 inlet"
   • region set = whichever caps exist  {LeftNose?, RightNose?, Mouth?, DescendingAirway}
        │
        ▼
[ septum split → per-region slice → CSA ]   ← UNCHANGED
   • closed-loop validity, centroid-proximity, plane-enclosure checks all preserved
```

### Label mapping (input cap name → existing pipeline name)

| Input cap  | Role             | Existing pipeline label |
|------------|------------------|-------------------------|
| `outlet`   | centerline source (caudal, trachea end) | `DescendingAirway` (shared trunk) |
| `leftnose` | centerline target | `LeftNose`  |
| `rightnose`| centerline target | `RightNose` |
| `mouth`    | centerline target | `Mouth`     |

## 5. Components

### 5.1 `AnatomicalCapAdapter` (new module, e.g. `slicer/anatomical_input.py`)

- `read_named_solids(path) -> dict[str, trimesh.Trimesh]`
  Custom parser splitting on `solid <name> … endsolid` (trimesh collapses multi-solid, so parse
  blocks ourselves). Returns one mesh per named solid.
- `classify_solids(solids) -> (wall_mesh, caps: list[Cap])`
  `wall` = the solid named `wall` (fallback: the solid with the most faces). Each remaining solid is
  a `Cap` with normalized+aliased anatomical `label`, `centroid` (area-weighted face centroid),
  and `normal` (area-weighted / PCA plane normal).
- `write_open_profile(wall, out_stl)` — writes the wall as the open-profile surface mesh.
- `write_seeds(caps, out_json)` — writes `[{label, centroid[x,y,z], normal[x,y,z]}]`.

### 5.2 `branch_detector.py` seeded mode (modify existing)

- New optional `--seeds <json>` argument threaded through `step0_branch_detect.py`.
- When seeds are present: after `detect_open_profiles(wall)`, match each detected hole barycenter to
  the nearest seed centroid (1:1, nearest-neighbor) and adopt that seed's anatomical label.
  Skip `classify_open_profiles()` entirely.
- `extract_centerlines`: source = the hole matched to `outlet`; targets = holes matched to inlets.
- **Do not run the `Branch_1/2/3` rename block** (the regression from commit `8bcafdf`) on the
  seeded/airway path — anatomical labels pass straight through to `branches/`. The rename stays only
  on the aortic path (guarded, e.g. by detector type or a `--positional-labels` flag).

### 5.3 `pipeline_bifurcation.py` adaptive regions (modify existing)

- Replace the `has_left and has_right` hard requirement (`:147–153`) with: require exactly one
  `DescendingAirway` (outlet) + at least one inlet. Log the detected configuration.
- Region loops (`:276` and the refine/slice steps) iterate over **present** regions only.
- Include `Mouth` in refine + slice (see 5.4). Missing regions are skipped cleanly, not errored.

### 5.4 Mouth region CSA (implement existing stub)

- `slice_bifurcation.py` currently stubs mouth slicing (~`:916-918`) and refine excludes it.
- Implement mouth as a first-class branch: it is a single passage (like the descending airway,
  not a bifurcation), so slice its sub-mesh with its own centerline using the existing
  `slice_along_centerline()` primitive (containment-filtered), producing the same per-region CSV /
  STL / plot outputs as LeftNose/RightNose/DescendingAirway.
- Mouth-only configs (no nose) must run end-to-end with no septum step.

## 6. Preserved behavior (explicitly unchanged)

- Closed-loop validity: `validation.py` closed-loop check (~`:284`), `is_valid` flag
  (`measurements.py:29`), stats/summary aggregate only `m.is_valid` planes (`measurements.py:301–333`).
  A slice that is not a closed polygon is flagged and excluded from CSA.
- Centroid-proximity and plane-enclosure validation.
- Loop filtering / doughnut merging, MeshPlaneSlicer containment behavior (including the known
  single-loop containment quirk — tracked separately, not touched here).
- Output CSV columns and per-region file layout.

## 7. Error handling

- No `wall` solid → hard error (need the lumen surface).
- No `outlet` (or alias) → hard error (no centerline source).
- Zero inlets → hard error.
- Unknown / unmatched cap name → **warn loudly and record in a manifest; never silently drop**
  (consistent with the "investigate and report, don't silently modify data" principle). The run
  proceeds with the recognized caps; the manifest lists what was skipped and why.
- Seed↔hole match ambiguity (a seed with no nearby hole, or two seeds contesting one hole) → error
  with the offending labels and distances, rather than a silent mislabel.

## 8. Testing

- **Unit — reader:** parse a synthetic ASCII multi-solid STL; case-insensitive + alias matching;
  each supported configuration (mouth-only, nose-only, all, asymmetric).
- **Unit — seeds:** cap centroid/normal correctness on a synthetic capped tube.
- **Unit — matching:** hole-barycenter ↔ seed-centroid nearest-neighbor assignment, including the
  ambiguity error path.
- **Integration — reader on real data:** run the reader on `ENT001_clean.stl` (real named multi-solid;
  names `inlet/outlet_1/wall`) to validate parsing + wall extraction; add an anatomically-named
  fixture (renamed patches or synthetic) to exercise the full label mapping.
- **Regression:** existing `slicer/` geometry, measurement, validation tests still pass unchanged.
- **End-to-end:** a mouth-only case and an all-openings case each produce per-region CSA with the
  closed-loop validity flags intact.

## 9. Resolved decisions

- **Integration:** pre-step feeding the existing pipeline (not a standalone pipeline, not a merged
  detector) — chosen for smallest, cleanest change reusing the trusted slicer.
- **Seeds transport:** sidecar `anatomical_seeds.json` read by `branch_detector` via `--seeds`
  (invisible at the CLI; keeps the "single STL in surface/" convention intact — JSON is not globbed
  as a surface).
- **Adapter invocation:** auto-detect inside `csa-bifurcation` — the pipeline inspects the surface
  STL and, if it is a named multi-solid file (ASCII with >1 `solid` block), runs the adapter
  automatically before Step 0. This delivers the "just a new format input, drop it in and run"
  experience. To avoid the two-STL glob ambiguity in `find_surface_stl`, the adapter writes the
  derived open-profile mesh + `anatomical_seeds.json` into `surface/_derived/`, and Step 0 is
  pointed at that explicit path (not the glob); the original named STL stays as the single file in
  `surface/`. An explicit `csa-uncap <stl>` subcommand is also exposed for running the adapter
  standalone (debugging / inspecting the extracted wall and seeds), but is not required in the
  normal flow.
