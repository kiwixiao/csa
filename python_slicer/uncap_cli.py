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
