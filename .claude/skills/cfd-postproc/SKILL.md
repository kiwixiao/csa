---
name: cfd-postproc
description: Run, reproduce or extend the airway CSA + CFD post-processing chain (csa-bifurcation CSA from the registration, then cfd_postproc/run_subject.sh for pressure, flow, WOB, report PDF and video). Use when setting up a new machine, processing a new subject, re-running ENT005-style results, checking a CFD run's mesh motion, or comparing an edited geometry.
---

# CFD post-processing (CSA + CFD)

The runbook is `cfd_postproc/README.md`. Read it first and follow it; it is written so a person can run the whole
chain without an agent, and the agent should do the same steps the same way.

Order of work:
1. New machine: `bash cfd_postproc/setup.sh`, `source cfd_postproc/activate.sh`, `bash cfd_postproc/check_env.sh` (must print `all OK`).
2. New subject: inputs (README section 2), t0 capped surface (section 3, show the user the review output), subject folder
   via `mirtk-prepare-slicer`, then copy `cfd_postproc/subjects/ENT005.env` to `<S>.env` and edit it.
3. Run: `cfd_postproc/run_subject.sh cfd_postproc/subjects/<S>.env all` (or named stages). It takes hours: run it in the
   background and report per stage.

Rules:
- Never skip `facecheck` or `motioncheck`. A facecheck failure means the CSA is frozen at t0; a motioncheck FAIL means
  the CFD mesh did not follow the registration and the CFD numbers are not valid. Stop and tell the user in both cases.
- Ask the user for the choices that are theirs: `WOB_REGIONS` boundaries, nostril depth / carina z if the sweep is
  ambiguous, whether septal holes are real.
- Compare two geometries on the same planes (`cfd_postproc/compare/`), never by the slicer's region labels.
- Pitfalls and their checks: README section 6. Reference outputs for ENT005: README section 7.
- To check a re-run against earlier results use `cfd_postproc/qa/compare_outputs.py` (file by file -> markdown), and
  explain every difference before calling it reproduced.
- Run pvbatch scripts only through `cfd_postproc/pvb.sh`; the CSA slicer only in conda env `mirtk`.
