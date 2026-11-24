# cfd_postproc: CSA + CFD post-processing for breathing-airway CFD

This folder holds the post-processing that turns a breathing-airway CFD run (STAR-CCM+, moving mesh) plus its
image registration into a report:
- cross-sectional area (CSA) along the airway
- pressure and flow at every cross-section
- work of breathing (WOB) by region
- a report PDF/DOCX and a report video

It lives inside the `csa` repo because the CSA step (`python_slicer`, `csa-bifurcation`) is the first half of the chain.

Everything here is meant to run without an agent:
1. Run one setup command on a new machine.
2. Write one config file per subject.
3. Run one driver command.

The checks that need a person are listed in section 4, with what to look at.

```
4-D CT --MIRTK registration--> ffd_*.dof.gz (breathing motion; image-to-image, the surface is NOT an input)
                                      |
t0 airway surface (capped, named) ----+-- csa-bifurcation  (conda env mirtk; step 0 calls env vmtk)
                                      v
            motion/stl/out_<ms>.stl (surface every 25 ms) + CSA planes + <S>-Data.csv   = CSA from the REGISTRATION
                                      |
STAR-CCM+ export (EnSight .case) -----+-- run_subject.sh  (ParaView pvbatch + python)
                                      v
            pressure / flow / energy flux on every CSA plane -> plots -> WOB -> report PDF/.docx -> report video
```

The report's CSA comes from the registration surfaces. The CFD supplies pressure, flow and velocity on those same planes.

---

## 1. Setup on a new machine

Requirements: Linux x86_64, conda >= 23 on PATH, internet, about 10 GB.

```bash
git clone git@github.com:kiwixiao/csa.git && cd csa          # this repo (branch with cfd_postproc)
bash cfd_postproc/setup.sh                                   # idempotent; see below
#   ParaView 6.0.1 already installed elsewhere:  setup.sh --pvbatch /path/to/bin/pvbatch
#   conda installed but not on PATH:             source <conda>/etc/profile.d/conda.sh  first
source cfd_postproc/activate.sh                              # every new shell: conda env cfdpost + PVBATCH
bash cfd_postproc/check_env.sh                               # must print "all OK"
```

`setup.sh` installs the following (anything already present is left as it is):

| piece | how | versions the reference results were made with |
|---|---|---|
| conda env **`mirtk`** | `mirtk-install/install.sh` (cloned next to this repo), then `pip install -e python_slicer` from this repo | MIRTK 2024.03.20, Python 3.11, vtk 9.4.2, trimesh 4.11.5 (`env/mirtk_env.yml`) |
| conda env **`vmtk`** | `env/vmtk_env.yml` | vmtk 1.5.0, Python 3.10, vtk 9.2.6. `csa-bifurcation` calls it **by this name** (`conda run -n vmtk`). |
| conda env **`cfdpost`** | `env/cfdpost_env.yml` | Python 3.11, numpy 1.26.4, scipy 1.13.1, matplotlib 3.9.2, pandas 2.2.3, pillow 10.4.0, ffmpeg, node |
| **ParaView 6.0.1** | `--pvbatch PATH`, or found in `/opt`, or binary download to `~/opt`; path written to `activate.sh` as `PVBATCH` | `ParaView-6.0.1-MPI-Linux-Python3.12-x86_64`, its own Python 3.12 (numpy 1.26.4, scipy 1.15.0) |
| npm `docx` | `npm install` here (`package.json`) | 9.15.4 |

Why this isn't one environment:
- MIRTK and vmtk need different, incompatible VTK/Python builds.
- ParaView brings its own Python.

The pinned versions reproduce the reference numbers. Different versions (especially VTK/ParaView) can change results slightly.

How each script runs:
- **`pvb.sh <script.py>`** runs it with ParaView's Python, for anything that reads the CFD, cuts planes or renders. It clears PYTHONPATH/PYTHONHOME and uses `$PVBATCH`.
- **`python3`** from env `cfdpost` runs plots, tables and composition.
- **`conda run -n mirtk ...`** runs the CSA slicer and MIRTK.

---

## 2. Per-subject inputs

| input | example (ENT005) | made by |
|---|---|---|
| CFD export, EnSight gold, one step per 25 ms | `CFD/ENT005_Case/ENT005.case` | STAR-CCM+. Check the time list for restarts and duplicates (section 6). |
| flow + pressure profile tables | `CFD/profiles_from_sim/ENT005FlowProfile_fromSim.csv` | `python3 tools/extract_sim_tables.py <run.sim> <outdir>` reads them straight from the `.sim` (no STAR licence) |
| registration folder | `ENT5_noalign_be0.001/` (`ffd_1..N.dof.gz`, `seg_0.stl`, `img_0.nii.gz`) | MIRTK pipeline (`mirtk-pipeline`) |
| registration control-point CSV | `ENT5_noalign_be0.001.csv` (points x per-ms xyz) | MIRTK pipeline; the same points the CFD motion tables use |
| t0 capped surface | `ENT005/surface/seg_0_boundary.stl` (solids Leftnose, Rightnose, outlet, wall) | section 3 |

### Subject folder for the CSA slicer
```bash
conda run -n mirtk mirtk-prepare-slicer --reg-dir <reg_dir> --output-dir <ROOT> --subject <S>
rm <ROOT>/<S>/surface/seg_0.stl                     # replaced by the capped, named surface:
cp seg_0_boundary.stl <ROOT>/<S>/surface/
```

---

## 3. The t0 capped surface (once per subject; look at the review output)

The input is the airway surface at t0 in mm, with the face/mask removed: nose + nasal cavity + pharynx + trachea + bronchi, open at the face opening and at the bronchi.
- ENT005's came from the STAR initial surface with the CFD Face/Mask boundaries removed.
- The CFD wall export (Wall + Nose blocks of the first step) also works; it needs `--fill-pinholes-mm 3`.

```bash
pvb.sh geometry/nostril_sweep.py airway_noface_mm.vtp                  # -> *_nostril_plane.json; prints rings per depth
pvb.sh geometry/build_capped_surface.py airway_noface_mm.vtp seg_0_boundary_auto.stl \
       --nostril-plane airway_noface_mm_nostril_plane.json --nostril-depth 4.5 --carina-z <z of carina, mm>
       # trachea cut 10 mm above the carina, DecimatePro 90 %; [--fill-pinholes-mm 3] for CFD-wall exports
pvb.sh geometry/find_septal_holes.py seg_0_boundary_auto.stl --front-of-y 40 --json septal_holes.json
pvb.sh geometry/patch_septal_holes.py seg_0_boundary_auto.stl seg_0_boundary.stl --front-of-y 38
python3 geometry/drop_slivers.py seg_0_boundary.stl seg_0_boundary.stl   # if qa/facecheck.py ever fails
```

What to look at:
- `nostril_sweep` must find a depth with exactly two closed nostril rings. Use that depth or slightly deeper (ENT005: found 4.0, used 4.5).
- The builder must report exactly 3 open rings: 2 nostrils + outlet.
- `patch_septal_holes` must end with `0 meeting edges`, meaning the septum separates left from right. Septal holes are closed in this CSA surface only; the CFD is not changed. Report holes as a possible clinical finding.
- ENT005's carina z was 1670 mm. `--front-of-y` depends on orientation: in ENT005's frame, larger y is anterior.

---

## 4. Run a subject: `run_subject.sh`

Copy `subjects/ENT005.env` to `subjects/<S>.env` and edit it. The settings are:
- **Inputs:** `CASE`, `CSA` (= `SUBJ_DIR/csa_bifurcation/<S>-Data.csv`), `PROF`, `REG_CSV`, `BREATH_MS` (breath length, for `csa-bifurcation --stop`).
- **Fitted values:** `RHO` (density fitted by the profile check), `VMAX` (velocity colour max, from the `velscan` stage).
- **Regions:** `REGIONS`, `PATHS_ANAT`, `MAIN_PATH`, `WOB_REGIONS` (region boundaries in mm from the nostril; your choice).
- **Output:** `POST`, plus `CSA_PLOT` / `PE_IN` / `PLANES_DIR` (the cleaned tables under `$POST/csa_clean`).

```bash
cfd_postproc/run_subject.sh cfd_postproc/subjects/<S>.env all          # everything, in order
cfd_postproc/run_subject.sh cfd_postproc/subjects/<S>.env plots regions pdf   # or any stages by name
```

| stage | what | ENT005 time | output | **check** |
|---|---|---|---|---|
| `csa` | `csa-bifurcation`: VMTK centerlines → septum split → surface deformed to every 25 ms → planes | ~20 min | `SUBJ_DIR/csa_bifurcation/<S>-Data.csv`, `frames/`, `motion/stl/` | expected frame count (ENT005: 84) |
| `facecheck` | face labels == faces of every deformed frame | ~1 min | `POST/facecheck.json` | **must be OK**; otherwise the slicer silently froze CSA at t0. Fix: `geometry/drop_slivers.py`, re-run `csa`. |
| `motioncheck` | registration control points vs CFD wall at every step | ~5 min | `POST/motion_check/<S>_motion_check.png` | **PASS** (median ≤ 0.1 mm, worst ≤ 1 mm). FAIL means the CFD mesh did not follow the motion and the CFD results are invalid. |
| `velscan` | velocity percentiles per step | ~2 min | `VELSCAN` | set `VMAX` ≈ max of the 99.9th percentile |
| `extract` | total pressure / flow / energy flux on every CSA plane, per region | ~2.5 h | `POST/csa_pe/` | |
| `clean` | nostril end: reject off-centerline planes, add fill planes from the cap to the first good plane | ~40 min | `POST/csa_clean/` | the printed fill check (re-cut vs slicer, median ~0 %) |
| `velgrid` | velocity sampling grid (manifest only) | ~1 min | `POST/velocity_breath_final/manifest.json` | |
| `plots` | anatomy / pressure / WOB per nostril path | ~1 min | `POST/csa_anatomy/`, `PE_IN/plots*/` | |
| `regions` | WOB per region and per 25 mm segment | <1 min | `POST/wob_regions/` | |
| `pdf` | report, one page per instant | ~4 min | `POST/<S>_CFD_report.pdf`, `.docx` | |
| `vellight` | velocity volume render, every CFD step | ~1.6 h | `POST/velocity_breath_light/` | |
| `reportvideo` | the report page animated over the breath | ~5 min | `POST/<S>_report_breath.mp4` | |

Also available: `qa/run_profile_check.sh <case> <flow.csv> <pres.csv> <outdir> <S> [workers]`. It compares the CFD boundary flow and pressure with the tables and fits `RHO`. Run it once per CFD run.

---

## 5. Edited geometry

**CFD of an edited geometry** (ENT005 "Aryt", see `subjects/ENT005_Aryt.env`):
- It is a separate CFD case processed on the baseline's planes.
- Its edited span is taken from `recut_csa.py` + `apply_csa_recut.py`.
- Run `motioncheck` on it first.

**Registration edit test** (CSA of an edited t0 surface, all from the registration):
1. Build the edited t0 surface the same way as section 3.
2. Run `csa` + `facecheck` on a second subject folder with the same registration.
3. Compare on the baseline's planes:
```bash
compare/run_samplane.sh <baseline_subject_dir> <edited_subject_dir> <out_dir> [workers]   # pvbatch, ~2 min
python3 compare/plot_samplane.py --csv <out_dir>/samplane_base_vs_edit.csv --prof <flow profile> --out fig.png \
        [--edit-zone 143:170] [--labels "baseline,edited"] [--title "..."]
```
- Compare on the same planes, not by region label. The slicer's septum split can move by cm between two surfaces.
- On the baseline surface, the re-cut reproduces the slicer table: ratio 1.000 on 99 % of planes (ENT005).

---

## 6. Pitfalls found on ENT005 (check these on every subject)

| problem | symptom | check / fix |
|---|---|---|
| STAR macro `TimeStep` ≠ simulation dt | CFD wall drifts away from the registration; WOB wrong | `motioncheck`; macro `TimeStep = 1.0` for 1 ms tables |
| restarts / duplicates in the export | repeated or out-of-order times in the `.case` | dedupe on (BronchiOut Q, cap area), keep the later copy |
| slicer face count ≠ labels | CSA identical in every frame, no warning | `facecheck`; `geometry/drop_slivers.py` |
| `csa-bifurcation` without `--stop` | wrong frame count | `BREATH_MS` in the .env (`csa` stage passes it) |
| septal holes | no LeftNose planes / left and right mixed | section 3 (find + patch) |
| VMTK septum split is sensitive | a small surface change moves the left-nose/descending boundary | compare geometries on the same planes (section 5) |
| nostril planes cross the cap | 8–24 mm² "pocket" CSA at s < 8 mm | `clean` stage |
| plane-STL numbering | wrong plane normals | rank over `<S>-Data-all.csv` (all planes), as the scripts do |
| GPU volume render segfaults (VisRTX), rarely | one render step dies | `pdf` / `vellight` retry a crashed step once; a second failure stops the run |
| pvbatch in a conda/anaconda shell | `import paraview` fails | always go through `pvb.sh` |
| flow sign / pressure time axis | inhale and exhale swapped, pressure offset | `qa/run_profile_check.sh` |

---

## 7. Reference run: ENT005 reproduced from this folder (2026-11-23)

`/storage4/ENT005/repro_test/` holds a full `all` run from this folder on the original inputs (~5.7 h on 1 GPU + 16 cores),
compared file by file with the existing ENT005 results by `repro_test/compare_all.sh` -> `repro_test/REPRO_REPORT.md`
(`qa/compare_outputs.py` does the comparing; use it the same way for any re-run).

| part | result |
|---|---|
| CSA pipeline (`csa`) | identical: tables (after row sort), face labels, all 84 motion STLs byte-identical |
| `facecheck` / `motioncheck` / `velscan` | OK / same FAIL numbers as before (native ENT005 mesh does not breathe) / byte-identical |
| `extract` | differs only on planes affected by the 2026-11-02 plane-STL rank fix (the reference predates it); every other plane identical. Whole-path WOB 211.900 mJ and peak dP 875.8 Pa unchanged; nose-region WOB 29.20 -> 29.57 mJ |
| `clean` -> `plots` -> `regions`, fed the reference extraction | every table and figure byte-identical, with the current python **and** with the `cfdpost` env from `setup.sh` |
| renders | repeat to max 9/255 run-to-run; the 2026-11-02 frames differ by a faint shade inside the mask cone (cosmetic) |

---

## 8. Layout

```
cfd_postproc/
  README.md            this runbook
  setup.sh, check_env.sh, activate.sh (written by setup.sh), pvb.sh, package.json
  run_subject.sh       driver (stages above)
  subjects/*.env       per-subject settings (ENT005, ENT005_Aryt)
  *.py                 pipeline scripts called by run_subject.sh
  geometry/            t0 surface: nostril_sweep, build_capped_surface, find/patch_septal_holes, drop_slivers
  qa/                  facecheck, reg_traj + wall_track + motion_check_fig, profile check, velocity scan
  compare/             same-plane CSA comparison of two geometries
  tools/               extract_sim_tables (profiles from a STAR .sim)
  env/                 environment specs (cfdpost, vmtk) + exported mirtk env record
```

Claude Code users get the project skill `.claude/skills/cfd-postproc/SKILL.md`, which points to this runbook.

On the original machine, `/storage4/ENT_postproc` is a link to this folder, so old paths still work.
