#!/usr/bin/env python3
"""Compare a reproduction run with reference results, file by file (same relative paths).

  .csv  -> rows/columns, max |diff| per numeric column, max relative diff; IDENTICAL / NUMERIC (<= --rtol) / DIFFERENT
  .json -> equal after parsing          .npy -> max |diff|          .png -> fraction of pixels that differ
  other -> byte compare
Usage: compare_outputs.py --ref DIR --new DIR [--glob 'csa_pe/*.csv' ...] [--rtol 1e-6] [--md report.md]
Pairs of single files can be given as --pair REF NEW (repeatable).
"""
import argparse, filecmp, glob, json, os
import numpy as np, pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument("--ref"); ap.add_argument("--new"); ap.add_argument("--glob", nargs="*", default=[])
ap.add_argument("--pair", nargs=2, action="append", default=[]); ap.add_argument("--rtol", type=float, default=1e-6)
ap.add_argument("--md", default=None); ap.add_argument("--title", default="comparison")
a = ap.parse_args()


def cmp_csv(r, n):
    st, info = _cmp_csv(pd.read_csv(r), pd.read_csv(n))
    if st == "DIFFERENT":                                  # parallel writers: same rows in another order?
        A, B = pd.read_csv(r), pd.read_csv(n)
        if list(A.columns) == list(B.columns) and len(A) == len(B):
            key = lambda D: D.sort_values(list(D.columns), kind="mergesort").reset_index(drop=True)
            st2, info2 = _cmp_csv(key(A), key(B))
            if st2 != "DIFFERENT":
                return st2, info2 + " (row order differs)"
    return st, info


def _cmp_csv(A, B):
    if list(A.columns) != list(B.columns):
        return "DIFFERENT", "columns differ: %s vs %s" % (sorted(set(A.columns) ^ set(B.columns)), "")
    if len(A) != len(B):
        return "DIFFERENT", "rows %d vs %d" % (len(A), len(B))
    num = [c for c in A.columns if pd.api.types.is_numeric_dtype(A[c]) and pd.api.types.is_numeric_dtype(B[c])]
    oth = [c for c in A.columns if c not in num]
    bad_txt = [c for c in oth if not (A[c].astype(str).values == B[c].astype(str).values).all()]
    if not num:
        return ("IDENTICAL" if not bad_txt else "DIFFERENT"), "%d rows, text columns %s" % (len(A), "equal" if not bad_txt else "differ: %s" % bad_txt)
    x, y = A[num].to_numpy(float), B[num].to_numpy(float)
    both_nan = np.isnan(x) & np.isnan(y)
    d = np.where(both_nan, 0, np.abs(x - y)); scale = np.maximum(np.abs(x), np.abs(y)); scale[scale == 0] = 1
    rel = np.where(both_nan, 0, d / scale); worst = num[int(np.nanargmax(np.nanmax(rel, axis=0)))]
    st = "IDENTICAL" if np.nanmax(d) == 0 and not bad_txt else ("NUMERIC" if np.nanmax(rel) <= a.rtol and not bad_txt else "DIFFERENT")
    return st, "%d rows x %d cols; max |diff| %.3g, max rel %.3g (%s)%s" % (len(A), len(A.columns), np.nanmax(d), np.nanmax(rel), worst,
                                                                              "; text differs: %s" % bad_txt if bad_txt else "")


def cmp_png(r, n):
    from PIL import Image
    A, B = np.asarray(Image.open(r).convert("RGB")).astype(int), np.asarray(Image.open(n).convert("RGB")).astype(int)
    if A.shape != B.shape:
        return "DIFFERENT", "size %s vs %s" % (A.shape, B.shape)
    f = (np.abs(A - B).sum(2) > 30).mean()
    return ("IDENTICAL" if f == 0 else "NUMERIC" if f < 1e-3 else "DIFFERENT"), "%.4f%% pixels differ" % (100 * f)


def compare(r, n):
    if not os.path.exists(n):
        return "MISSING", "not in the new run"
    if filecmp.cmp(r, n, shallow=False):
        return "IDENTICAL", "byte-identical"
    ext = os.path.splitext(r)[1].lower()
    try:
        if ext == ".csv": return cmp_csv(r, n)
        if ext == ".json": return ("IDENTICAL" if json.load(open(r)) == json.load(open(n)) else "DIFFERENT"), "parsed json"
        if ext == ".npy":
            x, y = np.load(r), np.load(n)
            return ("IDENTICAL" if x.shape == y.shape and np.array_equal(x, y) else "DIFFERENT"), "shape %s vs %s" % (x.shape, y.shape)
        if ext == ".png": return cmp_png(r, n)
    except Exception as e:
        return "ERROR", repr(e)[:200]
    return "DIFFERENT", "bytes differ"


pairs = list(a.pair)
for g in a.glob:
    for r in sorted(glob.glob(os.path.join(a.ref, g), recursive=True)):
        if os.path.isfile(r):
            pairs.append((r, os.path.join(a.new, os.path.relpath(r, a.ref))))
rows = [(os.path.relpath(r, a.ref) if a.ref and r.startswith(a.ref) else r, *compare(r, n)) for r, n in pairs]
w = max([len(x[0]) for x in rows] + [4])
for f, st, info in rows:
    print("%-*s  %-9s  %s" % (w, f, st, info))
tally = {s: sum(1 for x in rows if x[1] == s) for s in ("IDENTICAL", "NUMERIC", "DIFFERENT", "MISSING", "ERROR")}
print("summary:", tally)
if a.md:
    with open(a.md, "a") as fh:
        fh.write("\n### %s\n\n| file | result | detail |\n|---|---|---|\n" % a.title)
        for f, st, info in rows:
            fh.write("| `%s` | %s | %s |\n" % (f, st, info.replace("|", "/")))
        fh.write("\n%s\n" % ", ".join("%s %d" % kv for kv in tally.items() if kv[1]))
