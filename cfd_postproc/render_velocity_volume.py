#!/usr/bin/env pvbatch
# -*- coding: utf-8 -*-
"""
render_velocity_volume.py

General, OOM-safe ParaView pipeline that auto-detects a transient CFD case folder
(STAR-CCM+ EnSight Gold export, moving/deforming mesh) and volume-renders the
velocity magnitude flowing *inside* the duct/airway, one PNG per timestep.

Why this is OOM-safe for big 3D unstructured CFD data:
  * The classic ParaView "OOM / crash on big 3D volume" happens when you volume-render
    a large UNSTRUCTURED grid directly: it falls back to the projected-tetrahedra
    (unstructured) volume mapper, which allocates enormous intermediate buffers.
  * We instead (a) load ONLY the velocity arrays (not all ~36 exported fields),
    (b) ResampleToImage the unstructured grid onto a uniform vtkImageData whose
    resolution is capped, and (c) volume-render the *image* with the GPU smart
    volume mapper, which is stable and bounded in memory.
  * We process ONE timestep at a time and never load the whole time series.

Run:
  /opt/ParaView-6.0.1-MPI-Linux-Python3.12-x86_64/bin/pvbatch --force-offscreen-rendering \
      render_velocity_volume.py --case-dir /storage4/ENT003 [options]

Key options:
  --case-dir DIR     Folder to auto-scan for a *.case file (searched up to depth 2).
  --max-dim N        Max voxels along the longest axis of the resample grid (default 300).
  --max-voxels N     Hard cap on total resample voxels (default 40e6).
  --stride S         Render every S-th timestep (default 1 = all).
  --frames a:b       Restrict to timestep index range [a, b) (0-based).
  --field NAME       Scalar to render (default auto: VelocityMagnitude, else |Velocity|).
  --size WxH         Output image size (default 1280x960).
  --vmax V           Force top of color/opacity range (default: sampled global max).
  --preset NAME      Color preset (default "Turbo").
  --max-opacity O    Peak voxel opacity 0..1 (default 0.45; lower = see deeper inside).
  --bounds-samples N Timesteps sampled to estimate union bounds + vmax (default 12).
  --outdir DIR       Output dir (default <case-dir>/postproc/velocity_volume).
  --mp4              Also encode frames to velocity_volume.mp4 via ffmpeg if available.
  --dry-run          Detect + report only; render nothing.
"""

import argparse
import glob
import json
import os
import sys
import time

from paraview.simple import *  # noqa: F401,F403


# --------------------------------------------------------------------------- #
# Argument parsing (pvbatch passes script args after the script name)
# --------------------------------------------------------------------------- #
def parse_args(argv):
    p = argparse.ArgumentParser(description="OOM-safe velocity volume render for CFD EnSight cases")
    p.add_argument("--case-dir", required=True)
    p.add_argument("--case-file", default=None, help="Explicit .case file (skips auto-detect)")
    p.add_argument("--max-dim", type=int, default=300)
    p.add_argument("--max-voxels", type=float, default=40e6)
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--frames", default=None, help="index range a:b (0-based, b exclusive)")
    p.add_argument("--shard", default=None,
                   help="Round-robin worker shard 'w:N' for parallel multi-process rendering")
    p.add_argument("--field", default=None)
    p.add_argument("--size", default="1280x960")
    p.add_argument("--vmax", type=float, default=None)
    p.add_argument("--vmin", type=float, default=0.0)
    p.add_argument("--sampling-bounds", default=None,
                   help="Fixed resample box 'x0,x1,y0,y1,z0,z1' (keeps parallel workers consistent)")
    p.add_argument("--dims", default=None, help="Fixed resample grid 'nx,ny,nz'")
    p.add_argument("--preset", default="Turbo")
    p.add_argument("--max-opacity", type=float, default=0.45)
    p.add_argument("--fill-floor", type=float, default=0.14,
                   help="fill-lumen mode: opacity fraction at velocity 0 (near-wall). "
                        "Higher = flow fills the lumen out to the wall (default 0.14)")
    p.add_argument("--opacity-floor", type=float, default=0.01,
                   help="Fraction of vmax at which opacity starts rising above 0. Small = show "
                        "even slow flow so the lumen fills continuously (no disconnected gaps).")
    p.add_argument("--no-fill-lumen", dest="fill_lumen", action="store_false",
                   help="Disable the mask-based lumen fill (by default every in-mesh voxel is "
                        "shown, so velocity 0 renders as blue and the lumen fills completely).")
    p.set_defaults(fill_lumen=True)
    p.add_argument("--bounds-samples", type=int, default=12)
    p.add_argument("--outdir", default=None)
    p.add_argument("--azimuth", type=float, default=30.0)
    p.add_argument("--elevation", type=float, default=20.0)
    p.add_argument("--view", default=None,
                   help="Named camera view: iso|posX|negX|posZ|negZ|posY (overrides azimuth/elevation)")
    p.add_argument("--parallel-projection", action="store_true",
                   help="Orthographic projection (recommended for anatomical planes)")
    p.add_argument("--surface-opacity", type=float, default=0.0,
                   help="Overlay the duct wall surface at this opacity (0 = off)")
    p.add_argument("--surface-color", default="0.85,0.85,0.9")
    p.add_argument("--surface-source", default="stl:auto",
                   help="Wall surface: 'stl:<path>' (clean, static) | 'stl:auto' "
                        "(auto-find a *.stl near the case) | 'extract' (from the moving "
                        "mesh, but may show speckle from internal faces)")
    p.add_argument("--surface-scale", type=float, default=0.001,
                   help="Scale applied to an STL surface to match the mesh (STL in mm, mesh in m -> 0.001)")
    p.add_argument("--roll", type=float, default=0.0,
                   help="Roll the named view by this many degrees CCW about the view axis")
    p.add_argument("--outline", action="store_true", help="Draw the sampling bounding box outline")
    p.add_argument("--skip-frames", default="",
                   help="Comma-separated timestep indices to skip (corrupted moving-mesh export "
                        "steps). Get the list from scan_velocity_peak.py / velocity_peak.csv.")
    p.add_argument("--clamp-above", type=float, default=1.2,
                   help="Opacity safety net: voxels above this * vmax render transparent, so "
                        "stray garbage nodes never paint the volume. <=1 disables.")
    p.add_argument("--mp4", action="store_true")
    p.add_argument("--background", default=None,
                   help="R,G,B (0..1) render background; default dark. Use 1,1,1 for a white/PDF panel.")
    p.add_argument("--no-annotations", action="store_true",
                   help="Hide orientation axes + time annotation (clean panel for compositing)")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args(argv)


# Named anatomical camera views: (direction from center to camera, view-up).
# 'sagittal' = negX rolled 90 deg CCW (up = +Z): the chosen airway side view.
NAMED_VIEWS = {
    "iso":      ((1.0, 0.6, 1.0), (0, 1, 0)),
    "posX":     ((1.0, 0.0, 0.0), (0, 1, 0)),
    "negX":     ((-1.0, 0.0, 0.0), (0, 1, 0)),
    "posZ":     ((0.0, 0.0, 1.0), (0, 1, 0)),
    "negZ":     ((0.0, 0.0, -1.0), (0, 1, 0)),
    "posY":     ((0.0, 1.0, 0.0), (0, 0, -1)),
    "sagittal": ((-1.0, 0.0, 0.0), (0, 0, 1)),
}


def find_stl(case_dir):
    """Find a plausible wall STL near the case (prefer a full/mm one, not decimated)."""
    hits = []
    for depth in ("", "*/", "*/*/"):
        hits += glob.glob(os.path.join(case_dir, depth + "*.stl"))
    hits = sorted(set(os.path.realpath(h) for h in hits))
    if not hits:
        return None
    # prefer files hinting full-resolution millimetre geometry over decimated ones
    hits.sort(key=lambda p: ("decimated" in p.lower(), "mm" not in p.lower(), len(p)))
    return hits[0]


import builtins as _bi


def logline(msg, logfh=None):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    if logfh:
        logfh.write(line + "\n")
        logfh.flush()


def _restore_globals():
    """A ProgrammableFilter (used for 'part:' surface) executes `from ...numpy... import *`
    into this module's globals, shadowing builtins max/min. Our logger is named `logline`
    (numpy has no such name) so it is immune; this just restores max/min defensively."""
    g = globals()
    g["max"], g["min"] = _bi.max, _bi.min


# --------------------------------------------------------------------------- #
# Case detection
# --------------------------------------------------------------------------- #
def find_case_file(case_dir, explicit=None):
    if explicit:
        if not os.path.isfile(explicit):
            raise SystemExit("Given --case-file does not exist: %s" % explicit)
        return explicit
    hits = []
    for depth in ("", "*/", "*/*/"):
        hits += glob.glob(os.path.join(case_dir, depth + "*.case"))
    # de-dup, keep stable order, prefer shallowest
    seen, ordered = set(), []
    for h in hits:
        rp = os.path.realpath(h)
        if rp not in seen:
            seen.add(rp)
            ordered.append(h)
    if not ordered:
        raise SystemExit("No *.case file found under %s" % case_dir)
    ordered.sort(key=lambda x: x.count(os.sep))
    return ordered[0]


def pick_velocity_arrays(available):
    """Return (arrays_to_load, scalar_field_name, need_calc)."""
    load = []
    scalar = None
    have_vec = "Velocity" in available
    if have_vec:
        load.append("Velocity")
    if "VelocityMagnitude" in available:
        load.append("VelocityMagnitude")
        scalar = "VelocityMagnitude"
    # generic fallbacks for other cases
    if scalar is None:
        for cand in available:
            cl = cand.lower()
            if cl in ("velocitymagnitude", "vel_mag", "umag", "velmag"):
                load.append(cand)
                scalar = cand
                break
    if not have_vec and scalar is None:
        for cand in available:
            if cand.lower() in ("velocity", "u", "vel"):
                load.append(cand)
                have_vec = True
                break
    return load, scalar, (scalar is None and have_vec)


# --------------------------------------------------------------------------- #
# Grid sizing
# --------------------------------------------------------------------------- #
def compute_dims(bounds, max_dim, max_voxels):
    dx = max(bounds[1] - bounds[0], 1e-9)
    dy = max(bounds[3] - bounds[2], 1e-9)
    dz = max(bounds[5] - bounds[4], 1e-9)
    longest = max(dx, dy, dz)
    h = longest / float(max_dim)  # target ~cubic voxel edge
    nx = max(16, int(round(dx / h)))
    ny = max(16, int(round(dy / h)))
    nz = max(16, int(round(dz / h)))
    total = nx * ny * nz
    if total > max_voxels:
        s = (max_voxels / total) ** (1.0 / 3.0)
        nx = max(16, int(nx * s))
        ny = max(16, int(ny * s))
        nz = max(16, int(nz * s))
    return [nx, ny, nz]


def union_bounds(a, b):
    if a is None:
        return list(b)
    return [min(a[0], b[0]), max(a[1], b[1]),
            min(a[2], b[2]), max(a[3], b[3]),
            min(a[4], b[4]), max(a[5], b[5])]


def _median(xs):
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return 0.0
    m = n // 2
    return s[m] if n % 2 else 0.5 * (s[m - 1] + s[m])


def robust_vmax(step_max, logfh=None):
    """Pick a physical color ceiling from per-step maxima, rejecting spike outliers.

    step_max: list of (timestep_index, reported_max). Rejects any value more than
    12*MAD above the median (garbage spikes like 1e23), returns the max of the rest.
    """
    vals = [v for _, v in step_max if v == v]  # drop NaN
    if not vals:
        return 1.0
    med = _median(vals)
    mad = _median([abs(v - med) for v in vals])
    thresh = med + 12.0 * mad + 1e-12
    kept, rejected = [], []
    for i, v in step_max:
        (kept if v <= thresh else rejected).append((i, v))
    if rejected and logfh is not None:
        logline("Rejected %d spike step(s) from range calc: %s"
            % (len(rejected), ", ".join("idx%d=%.3g" % (i, v) for i, v in rejected)), logfh)
    vmax = max((v for _, v in kept), default=med)
    # final guard: if everything looked like a spike, fall back to the median
    if vmax > 1e6 and med < 1e6:
        vmax = med
    return vmax if vmax > 0 else 1.0


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    args = parse_args(sys.argv[1:])
    case_dir = os.path.abspath(args.case_dir)
    case_file = find_case_file(case_dir, args.case_file)

    outdir = args.outdir or os.path.join(case_dir, "postproc", "velocity_volume")
    frames_dir = os.path.join(outdir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    logfh = open(os.path.join(outdir, "render.log"), "a")

    logline("=" * 70, logfh)
    logline("Case dir : %s" % case_dir, logfh)
    logline("Case file: %s" % case_file, logfh)
    logline("Out dir  : %s" % outdir, logfh)

    # ---- open reader, detect fields/timesteps ----
    reader = OpenDataFile(case_file)
    available = list(reader.PointArrays.Available)
    load, scalar, need_calc = pick_velocity_arrays(available)
    if not load:
        raise SystemExit("Could not find any velocity array. Available: %s" % available)
    reader.PointArrays = load
    logline("Available point arrays (%d): %s" % (len(available), ", ".join(available)), logfh)
    logline("Loading only: %s" % ", ".join(load), logfh)

    tvals = list(reader.TimestepValues) if reader.TimestepValues else [0.0]
    n_t = len(tvals)
    logline("Timesteps: %d  (t=%.4g .. %.4g)" % (n_t, tvals[0], tvals[-1]), logfh)

    # ---- source of scalar field ----
    src = reader
    field = args.field or scalar
    if field is None and need_calc:
        calc = Calculator(Input=reader)
        calc.ResultArrayName = "VelocityMagnitude"
        calc.Function = "mag(Velocity)"
        calc.AttributeType = "Point Data"
        src = calc
        field = "VelocityMagnitude"
    if field is None:
        raise SystemExit("No scalar field resolved to color by.")
    logline("Rendering scalar field: %s" % field, logfh)

    # ---- select frame indices ----
    idx = list(range(0, n_t, max(1, args.stride)))
    if args.frames:
        a, b = args.frames.split(":")
        a = int(a) if a else 0
        b = int(b) if b else n_t
        idx = [i for i in idx if a <= i < b]
    # round-robin sharding for multi-worker parallel rendering: --shard "w:N" keeps
    # the frames at list positions where position % N == w (balances slow late frames
    # across workers). All shards MUST share fixed bounds/dims/vmax/camera for consistency.
    if args.shard:
        w, nsh = (int(x) for x in args.shard.split(":"))
        idx = [idx[k] for k in range(len(idx)) if k % nsh == w]
        logline("Shard %d/%d -> %d frames" % (w, nsh, len(idx)), logfh)
    logline("Frames to render: %d (stride=%d)" % (len(idx), args.stride), logfh)

    # ---- pass 1: union bounds + robust global vmax over a few sampled steps ----
    # NOTE: transient moving-mesh CFD exports can contain rare spike/garbage nodes
    # (e.g. a single node reading ~1e23). Using the raw max as the color range would
    # collapse the real flow (a few m/s) to invisible. We take each step's reported
    # max but reject outliers via median-absolute-deviation so vmax tracks the real
    # physical range. Override with --vmax if you want an exact ceiling.
    fixed_ub = [float(x) for x in args.sampling_bounds.split(",")] if args.sampling_bounds else None
    fixed_dims = [int(x) for x in args.dims.split(",")] if args.dims else None
    # We can skip the (expensive) sampling pass entirely only when bounds, dims AND
    # vmax are all pre-supplied -- that is how parallel frame-workers start fast and
    # stay perfectly consistent with each other.
    need_sampling = not (fixed_ub and fixed_dims and args.vmax is not None)

    ub = fixed_ub
    step_max = []
    if need_sampling:
        n_s = max(2, min(args.bounds_samples, len(idx)))
        sample_idx = [idx[int(round(k * (len(idx) - 1) / (n_s - 1)))] for k in range(n_s)]
        sample_idx = sorted(set(sample_idx))
        logline("Estimating union bounds + range from %d sampled steps..." % len(sample_idx), logfh)
        acc = None
        for i in sample_idx:
            t = tvals[i]
            src.UpdatePipeline(t)
            di = src.GetDataInformation()
            acc = union_bounds(acc, di.GetBounds())
            ai = di.GetPointDataInformation().GetArrayInformation(field)
            if ai is not None:
                rng = ai.GetComponentRange(-1 if ai.GetNumberOfComponents() > 1 else 0)
                step_max.append((i, rng[1]))
        if ub is None:
            pad = 0.01 * max(acc[1] - acc[0], acc[3] - acc[2], acc[5] - acc[4])
            ub = [acc[0] - pad, acc[1] + pad, acc[2] - pad,
                  acc[3] + pad, acc[4] - pad, acc[5] + pad]
    else:
        logline("Using fixed bounds/dims/vmax (no sampling pass).", logfh)

    if args.vmax is not None:
        vmax = args.vmax
    else:
        # NOTE: robust_vmax uses reported per-step maxima, which are UNRELIABLE here
        # (the EnSight reader can report inf for clean steps too). Prefer passing an
        # explicit --vmax computed from the numpy scan (velocity_peak.csv clean peak).
        vmax = robust_vmax(step_max, logfh)
    dims = fixed_dims if fixed_dims else compute_dims(ub, args.max_dim, args.max_voxels)
    logline("Union bounds: [%.4g,%.4g] [%.4g,%.4g] [%.4g,%.4g]" %
        (ub[0], ub[1], ub[2], ub[3], ub[4], ub[5]), logfh)
    logline("Resample dims: %s  (%.2fM voxels)  vmax=%.4g" %
        (dims, dims[0] * dims[1] * dims[2] / 1e6, vmax), logfh)

    manifest = {
        "case_file": case_file, "field": field, "n_timesteps": n_t,
        "time_values": tvals, "frames_index": idx, "union_bounds": ub,
        "resample_dims": dims, "vmin": args.vmin, "vmax": vmax,
        "loaded_arrays": load, "size": args.size, "preset": args.preset,
    }
    with open(os.path.join(outdir, "manifest.json"), "w") as mf:
        json.dump(manifest, mf, indent=2)

    if args.dry_run:
        logline("Dry run complete. No frames rendered.", logfh)
        return

    # ---- build resample + volume pipeline (fixed grid = stable animation) ----
    resample = ResampleToImage(Input=src)
    resample.UseInputBounds = 0
    resample.SamplingBounds = ub
    resample.SamplingDimensions = dims

    W, H = (int(x) for x in args.size.lower().split("x"))
    rv = GetActiveViewOrCreate("RenderView")
    rv.ViewSize = [W, H]
    rv.OrientationAxesVisibility = 0 if args.no_annotations else 1
    rv.Background = [float(x) for x in args.background.split(",")] if args.background else [0.05, 0.05, 0.08]
    try:
        rv.UseColorPaletteForBackground = 0
    except Exception:
        pass

    # Lumen fill: ResampleToImage tags every voxel with vtkValidPointMask (1 inside the
    # mesh, 0 outside). Use it so IN-mesh voxels are always shown (velocity 0 -> blue, the
    # whole lumen fills solid) while OUT-of-mesh space stays transparent. Outside voxels get
    # a large-negative sentinel that the opacity map drops; inside keeps its real velocity.
    SENT = -1.0e4
    cfield = field
    vol_src = resample
    if args.fill_lumen:
        calc = Calculator(Input=resample)
        calc.AttributeType = "Point Data"
        calc.ResultArrayName = "vmag_lumen"
        calc.Function = "vtkValidPointMask*%s + (vtkValidPointMask-1)*%g" % (field, -SENT)
        cfield = "vmag_lumen"
        vol_src = calc
        logline("Lumen fill ON: in-mesh voxels always shown (0 -> blue)", logfh)

    disp = Show(vol_src, rv)
    disp.SetRepresentationType("Volume")
    ColorBy(disp, ("POINTS", cfield))

    # color transfer function (colored over [vmin, vmax]; sentinel clamps to bottom color)
    ctf = GetColorTransferFunction(cfield)
    ctf.RescaleTransferFunction(args.vmin, vmax)
    try:
        ctf.ApplyPreset(args.preset, True)
    except Exception as e:
        logline("preset '%s' failed (%s), keeping default" % (args.preset, e), logfh)

    otf = GetOpacityTransferFunction(cfield)
    mo = args.max_opacity
    span = vmax - args.vmin
    if args.fill_lumen:
        # in-mesh velocity 0 is VISIBLE (blue); only the outside sentinel is transparent
        pts = [
            SENT, 0.0, 0.5, 0.0,
            args.vmin - 1e-6 * span, 0.0, 0.5, 0.0,
            args.vmin, mo * args.fill_floor, 0.5, 0.0,    # velocity 0 inside -> blue, shown (near-wall)
            args.vmin + 0.40 * span, mo * max(0.5, args.fill_floor + 0.25), 0.5, 0.0,
            args.vmin + 0.75 * span, mo * 0.78, 0.5, 0.0,
            vmax, mo, 0.5, 0.0,
        ]
        lo = SENT
    else:
        fl = max(1e-6, args.opacity_floor) * span
        pts = [
            args.vmin, 0.0, 0.5, 0.0,
            args.vmin + fl, mo * 0.16, 0.5, 0.0,
            args.vmin + 0.40 * span, mo * 0.45, 0.5, 0.0,
            args.vmin + 0.75 * span, mo * 0.75, 0.5, 0.0,
            vmax, mo, 0.5, 0.0,
        ]
        lo = args.vmin
    if args.clamp_above > 1.0:
        hi = args.clamp_above * vmax
        pts += [hi, 0.0, 0.5, 0.0]           # fade to transparent just above vmax
        otf.RescaleTransferFunction(lo, hi)
    else:
        otf.RescaleTransferFunction(lo, vmax)
    otf.Points = pts
    try:
        ctf.NanOpacity = 0.0
    except Exception:
        pass
    # opacity unit distance ~ voxel size keeps opacity meaningful regardless of dims
    vsize = min((ub[1] - ub[0]) / dims[0], (ub[3] - ub[2]) / dims[1], (ub[5] - ub[4]) / dims[2])
    try:
        disp.ScalarOpacityUnitDistance = vsize
    except Exception:
        pass

    # faint context: bounding outline (off by default; on with --outline)
    if args.outline:
        outline = Outline(Input=resample)
        od = Show(outline, rv)
        od.DiffuseColor = [0.4, 0.4, 0.5]
        od.Opacity = 0.3

    # optional semi-solid duct wall surface.
    #   'part:<n>'                -> EnSight boundary surface block <n> (e.g. part:1 = the
    #                                wall). MOVES with the mesh each timestep, clean (it is a
    #                                real boundary mesh, not a dirty volume-surface extraction).
    #   'stl:<path>' / 'stl:auto' -> clean watertight STL; STATIC reference (does not deform).
    #   'extract'                 -> surface of the moving mesh; deforms but speckles.
    surf = None
    surf_moving = False
    if args.surface_opacity > 0.0:
        src_spec = args.surface_source
        if src_spec.startswith("part:"):
            spec = src_spec[5:]
            surf = ProgrammableFilter(Input=reader)
            surf.OutputDataSetType = "vtkUnstructuredGrid"
            if spec == "all":
                # merge every boundary surface block (1..N-1); block 0 is the fluid volume
                surf.Script = (
                    "from vtkmodules.vtkFiltersCore import vtkAppendFilter\n"
                    "inp=self.GetInputDataObject(0,0)\n"
                    "app=vtkAppendFilter()\n"
                    "for i in range(1, inp.GetNumberOfBlocks()):\n"
                    "    b=inp.GetBlock(i)\n"
                    "    if b is not None: app.AddInputData(b)\n"
                    "app.Update()\n"
                    "self.GetOutputDataObject(0).ShallowCopy(app.GetOutput())")
                logline("Wall surface: ALL EnSight boundary parts merged (moving)", logfh)
            else:
                bidx = int(spec)
                surf.Script = ("inp=self.GetInputDataObject(0,0); "
                               "self.GetOutputDataObject(0).ShallowCopy(inp.GetBlock(%d))" % bidx)
                logline("Wall surface: EnSight boundary part %d (moving)" % bidx, logfh)
            surf_moving = True
        elif src_spec.startswith("stl:"):
            stl_path = src_spec[4:]
            if stl_path == "auto":
                stl_path = find_stl(case_dir)
            if not stl_path or not os.path.isfile(stl_path):
                logline("Surface: no STL found (%s); falling back to mesh extract." % stl_path, logfh)
                src_spec = "extract"
            else:
                raw = OpenDataFile(stl_path)
                surf = Transform(Input=raw)
                surf.Transform.Scale = [args.surface_scale] * 3
                logline("Wall surface: STL %s (scale %.4g, static)" % (stl_path, args.surface_scale), logfh)
        if src_spec == "extract":
            surf = ExtractSurface(Input=MergeBlocks(Input=reader))
            surf_moving = True
            logline("Wall surface: extracted from moving mesh (may show speckle)", logfh)

        sdisp = Show(surf, rv)
        _restore_globals()  # Show(ProgrammableFilter) may shadow our log/max/min
        sdisp.SetRepresentationType("Surface")
        try:
            ColorBy(sdisp, None)  # turn off scalar coloring (arrayed surfaces)
        except Exception:
            pass  # arrayless STL already renders solid via DiffuseColor
        sdisp.DiffuseColor = [float(x) for x in args.surface_color.split(",")]
        sdisp.AmbientColor = sdisp.DiffuseColor
        sdisp.Opacity = args.surface_opacity
        for prop, val in (("Ambient", 0.45), ("Diffuse", 0.6), ("Specular", 0.0)):
            try:
                setattr(sdisp, prop, val)
            except Exception:
                pass
        logline("Duct wall surface overlay enabled (opacity=%.2f)" % args.surface_opacity, logfh)

    # solution time annotation, bottom-left corner
    try:
        if args.no_annotations:
            raise RuntimeError("annotations disabled")
        atf = AnnotateTimeFilter(Input=reader)
        atf.Format = "t = {time:.3f} s"
        ad = Show(atf, rv)
        ad.FontSize = 14
        try:
            ad.WindowLocation = "Lower Left Corner"
        except Exception:
            ad.WindowLocation = "LowerLeftCorner"
    except Exception as e:
        logline("time annotation unavailable: %s" % e, logfh)

    # color legend (skipped for clean composite panels; we draw our own colorbar)
    disp.SetScalarBarVisibility(rv, not args.no_annotations)
    if not args.no_annotations:
        sb = GetScalarBar(ctf, rv)
        sb.Title = field
        sb.ComponentTitle = "m/s"

    # camera (fixed for whole sequence, since bounds are fixed)
    src.UpdatePipeline(tvals[idx[0]])
    resample.UpdatePipeline(tvals[idx[0]])
    if surf is not None and surf_moving:
        surf.UpdatePipeline(tvals[idx[0]])
        _restore_globals()
    center = [(ub[0] + ub[1]) / 2, (ub[2] + ub[3]) / 2, (ub[4] + ub[5]) / 2]
    rv.CenterOfRotation = center
    if args.parallel_projection:
        rv.CameraParallelProjection = 1
    if args.view and args.view in NAMED_VIEWS:
        d, up = NAMED_VIEWS[args.view]
        n = (d[0] ** 2 + d[1] ** 2 + d[2] ** 2) ** 0.5
        diag = ((ub[1] - ub[0]) ** 2 + (ub[3] - ub[2]) ** 2 + (ub[5] - ub[4]) ** 2) ** 0.5
        rv.CameraFocalPoint = center
        rv.CameraPosition = [center[j] + (d[j] / n) * diag * 2.0 for j in range(3)]
        rv.CameraViewUp = up
        ResetCamera(rv)
        if args.roll:
            GetActiveCamera().Roll(args.roll)
        logline("Camera: named view '%s'%s" % (args.view,
            (" roll %+g deg" % args.roll) if args.roll else ""), logfh)
    else:
        ResetCamera(rv)
        cam = GetActiveCamera()
        cam.Azimuth(args.azimuth)
        cam.Elevation(args.elevation)

    # ---- render loop, one timestep at a time ----
    skip_set = set(int(x) for x in args.skip_frames.split(",") if x.strip() != "")
    if skip_set:
        logline("Skip list (corrupted steps): %s" % sorted(skip_set), logfh)
    logline("Rendering up to %d frames at %dx%d ..." % (len(idx), W, H), logfh)
    t_start = time.time()
    ndigits = _bi.max(4, len(str(n_t)))  # _bi.max: immune to ProgrammableFilter numpy shadowing
    skipped = []
    rendered = 0
    for k, i in enumerate(idx):
        t = tvals[i]
        try:
            if i in skip_set:
                logline("  frame idx %d (t=%.3f) SKIPPED (in --skip-frames)" % (i, t), logfh)
                skipped.append(i)
                continue
            rv.ViewTime = t
            resample.UpdatePipeline(t)
            if surf is not None and surf_moving:
                surf.UpdatePipeline(t)
                _restore_globals()  # in case the part-surface ProgrammableFilter re-pollutes
            Render(rv)
            fn = os.path.join(frames_dir, "frame_%0*d.png" % (ndigits, i))
            SaveScreenshot(fn, rv, ImageResolution=[W, H])
            rendered += 1
            if k % 10 == 0 or k == len(idx) - 1:
                el = time.time() - t_start
                rate = (k + 1) / el if el > 0 else 0
                eta = (len(idx) - k - 1) / rate if rate > 0 else 0
                logline("  frame %d/%d (idx %d, t=%.3f)  %.1fs elapsed  ETA %.0fs"
                    % (k + 1, len(idx), i, t, el, eta), logfh)
        except Exception as e:
            logline("  !! frame idx %d (t=%.3f) FAILED: %s" % (i, t, e), logfh)

    if skipped:
        logline("Skipped %d corrupted frame(s): %s" % (len(skipped), skipped), logfh)
    logline("Done. %d frames rendered (%d skipped) in %s" % (rendered, len(skipped), frames_dir), logfh)

    # ---- optional mp4 ----
    if args.mp4:
        import subprocess
        import shutil
        ff = shutil.which("ffmpeg")
        if ff:
            mp4 = os.path.join(outdir, "velocity_volume.mp4")
            # glob input tolerates strided frame numbering / gaps
            cmd = [ff, "-y", "-framerate", "20", "-pattern_type", "glob",
                   "-i", os.path.join(frames_dir, "frame_*.png"),
                   "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", mp4]
            logline("Encoding mp4: %s" % " ".join(cmd), logfh)
            try:
                subprocess.run(cmd, check=True)
                logline("Wrote %s" % mp4, logfh)
            except Exception as e:
                logline("ffmpeg failed: %s" % e, logfh)
        else:
            logline("ffmpeg not found; skipping mp4.", logfh)

    logfh.close()


if __name__ == "__main__":
    main()
