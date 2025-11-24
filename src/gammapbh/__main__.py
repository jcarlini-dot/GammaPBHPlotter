from __future__ import annotations

import argparse
import sys
import csv
from typing import Iterable, Tuple, List

import numpy as np

from . import (
    __version__,
    discover_mass_folders,
    load_spectra_components,
)

def _flatten_mass_names(masses) -> List[str]:
    # Accept (nums, names) or a single list
    if isinstance(masses, tuple) and len(masses) == 2:
        a, b = masses
        if a and isinstance(a[0], (int, float)):
            # b should be names
            return [str(x) for x in b]
        else:
            # a should be names
            return [str(x) for x in a]
    if isinstance(masses, list):
        # If numeric, format; else keep as string
        if masses and isinstance(masses[0], (int, float)):
            return [f"{float(x):.1e}" for x in masses]
        return [str(x) for x in masses]
    return []

def _compute_total(comps: dict) -> np.ndarray | None:
    if "total" in comps:
        return np.asarray(comps["total"], dtype=float)
    pri = comps.get("total_primary")
    sec = comps.get("total_secondary")
    if pri is not None and sec is not None:
        return np.asarray(pri, dtype=float) + np.asarray(sec, dtype=float)
    return None

def _write_csv(path: str, E: np.ndarray, comps: dict):
    keys = ["E"] + sorted(comps.keys())
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for i in range(E.size):
            row = [E[i]] + [np.asarray(comps[k])[i] for k in keys[1:]]
            w.writerow(row)

def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gammapbh",
        description="GammaPBH: inspect/plot bundled PBH gamma-ray spectra.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    parser.add_argument("--list-masses", action="store_true",
                        help="List available mass folders and exit.")
    parser.add_argument("--mass", type=float,
                        help="Select a mass from the grid (e.g., 5.0e13).")
    parser.add_argument("--align", choices=["secondary", "primary", "union", "intersection"],
                        default="secondary", help="Energy grid alignment (default: secondary).")

    parser.add_argument("--plot", action="store_true",
                        help="Plot the total spectrum (requires matplotlib).")
    parser.add_argument("--save", type=str, default=None,
                        help="File path to save the plot (if --plot) or CSV (if --csv).")
    parser.add_argument("--csv", action="store_true",
                        help="Write full spectra (E + components) to CSV (use --save).")
    parser.add_argument("--head", type=int, default=10,
                        help="When not plotting or CSV, print first N rows (default: 10).")
    parser.add_argument("--no-show", action="store_true",
                        help="With --plot, do not open an interactive window.")

    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.list_masses:
        masses = discover_mass_folders()
        names = _flatten_mass_names(masses)
        if not names:
            print("No masses found.", file=sys.stderr)
            return 1
        print("Mass grid entries:")
        for s in names:
            print("  ", s)
        return 0

    if args.mass is None:
        # No subcommand chosen; show help and succeed
        parser.print_help()
        return 0

    # Load spectra
    E, comps = load_spectra_components(args.mass, align=args.align)
    E = np.asarray(E, dtype=float)

    if args.csv:
        if not args.save:
            print("Error: --csv requires --save <file.csv>", file=sys.stderr)
            return 2
        _write_csv(args.save, E, comps)
        print(f"Wrote CSV: {args.save}")
        return 0

    if args.plot:
        import matplotlib.pyplot as plt  # local import keeps CLI light if not plotting
        y = _compute_total(comps)
        if y is None:
            print("No 'total' spectrum found; cannot plot.", file=sys.stderr)
            return 2
        plt.figure()
        plt.loglog(E, y, label=f"total ({args.align})")
        plt.xlabel("Energy")
        plt.ylabel("dN/dE (arb.)")
        plt.title(f"PBH spectra, M = {args.mass:.3e}")
        plt.grid(True, which="both", ls=":")
        plt.legend()
        if args.save:
            plt.savefig(args.save, dpi=150, bbox_inches="tight")
            print(f"Saved plot: {args.save}")
        if not args.no_show:
            plt.show()
        return 0

    # Text preview (head)
    total = _compute_total(comps)
    keys = ["total"] if total is not None else sorted(comps.keys())
    n = min(args.head, E.size)
    print(f"E grid points: {E.size}; showing first {n} rows")
    print("Columns:", "E", *keys)
    for i in range(n):
        vals = [np.asarray(comps[k])[i] if k != "total" else total[i] for k in keys]
        print(E[i], *vals)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())