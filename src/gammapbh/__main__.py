# __main__.py
"""
Command-line entry point for the ``gammapbh`` package.

This module implements a thin, dependency-light CLI around the public library
functions that (a) discover the available PBH mass grid from the bundled data
and (b) load the precomputed spectral components for a selected PBH mass.
It is intentionally simple: the heavy lifting (data layout, interpolation,
component assembly) is done in the library; the CLI is a convenience layer for
listing masses, previewing arrays on the console, exporting to CSV, and making
a quick log–log plot.

Typical usage
-------------
From a shell, after installing the package::

    # Show help and available options
    python -m gammapbh -h

    # List the mass grid available in the bundled data
    python -m gammapbh --list-masses

    # Preview the first 10 rows (default) of the spectra for a mass
    python -m gammapbh --mass 3.0e15

    # Export all columns (E and components) to CSV
    python -m gammapbh --mass 3.0e15 --csv --save spectra_3e15.csv

    # Make a quick plot of the total spectrum (no interactive window)
    python -m gammapbh --mass 3.0e15 --plot --no-show --save spectra_3e15.png

Notes
-----
* The CLI does **not** download or generate new data. It only reads the
  pre-bundled tables and composes them.
* Energy grid alignment can be controlled via ``--align``; by default, the
  secondary grid is used.
* The CSV writer includes **all** component columns present for the chosen
  mass; column names mirror the underlying data keys (e.g. ``total_primary``,
  ``total_secondary``, etc.).

Exit status
-----------
0 on success; non-zero on error (e.g. missing data, invalid arguments).
"""

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
    """
    Normalize a heterogeneous "masses" descriptor into a list of string names.

    The discovery helper may return either:
      * a single list (of numbers or strings), or
      * a 2-tuple ``(nums, names)``.

    This utility converts those shapes into a flat ``List[str]`` of displayable
    folder names (e.g. ``["5.0e13", "1.0e14", ...]``).

    Parameters
    ----------
    masses :
        Either a list of masses (numbers or strings) or a tuple
        ``(numbers, names)``.

    Returns
    -------
    List[str]
        Mass identifiers formatted as strings, suitable for display.
    """
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
    """
    Compute (or retrieve) the total spectrum from component columns.

    The function prefers an explicit ``"total"`` key if present. Otherwise,
    it adds ``"total_primary"`` and ``"total_secondary"`` if both exist.
    If neither path is possible, returns ``None``.

    Parameters
    ----------
    comps : dict
        Mapping of component name -> 1D array (NumPy-coercible).

    Returns
    -------
    numpy.ndarray | None
        The total spectrum array (same shape as energy grid), or ``None`` if
        an unambiguous total cannot be determined.
    """
    if "total" in comps:
        return np.asarray(comps["total"], dtype=float)
    pri = comps.get("total_primary")
    sec = comps.get("total_secondary")
    if pri is not None and sec is not None:
        return np.asarray(pri, dtype=float) + np.asarray(sec, dtype=float)
    return None


def _write_csv(path: str, E: np.ndarray, comps: dict):
    """
    Write energy and all component arrays to a CSV file.

    The first row is a header of column names: ``E`` followed by the sorted
    component keys. Each subsequent row corresponds to a single energy bin.

    Parameters
    ----------
    path : str
        Output CSV file path.
    E : numpy.ndarray
        1D energy grid (same length as each component column).
    comps : dict
        Mapping of component name -> 1D array (NumPy-coercible).
    """
    keys = ["E"] + sorted(comps.keys())
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for i in range(E.size):
            row = [E[i]] + [np.asarray(comps[k])[i] for k in keys[1:]]
            w.writerow(row)


def main(argv: Iterable[str] | None = None) -> int:
    """
    CLI entry point.

    Parameters
    ----------
    argv : Iterable[str] | None
        Sequence of command-line arguments *excluding* the program name.
        If ``None`` (default), arguments are taken from ``sys.argv``.

    Returns
    -------
    int
        Process exit code (0 on success).

    Command-line options
    --------------------
    --list-masses
        List available mass folders and exit.
    --mass FLOAT
        Select a PBH mass from the grid (e.g., ``5.0e13``).
    --align {secondary,primary,union,intersection}
        Choose the energy grid alignment (default: ``secondary``).
    --plot
        Plot the total spectrum (requires Matplotlib).
    --csv
        Write full spectra (E + components) to CSV (requires ``--save``).
    --save PATH
        Path to write the plot (PNG, PDF, etc.) or CSV, depending on flags.
    --head N
        Number of preview rows to print when not plotting or exporting (default: 10).
    --no-show
        With ``--plot``, do not open an interactive window (useful in batch runs).

    Examples
    --------
    Preview data::

        python -m gammapbh --mass 3e15 --head 5

    Export to CSV::

        python -m gammapbh --mass 3e15 --csv --save out.csv

    Save a plot without showing a window::

        python -m gammapbh --mass 3e15 --plot --no-show --save out.png
    """
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
