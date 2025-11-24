# src/gammapbh/__main__.py
"""
Module entry point so `python -m gammapbh [args]` works.

This simply forwards argv to `gammapbh.cli.main`, so things like:
  - python -m gammapbh --help
  - python -m gammapbh monochromatic --help
  - python -m gammapbh  (interactive CLI, if your main supports it)
behave exactly the same as running the `gammapbh` console script.
"""

from __future__ import annotations
import sys


def _entry() -> int:
    # Import inside the function so import errors show as runtime errors
    # when the module is invoked (mirrors how console_scripts work).
    from .cli import main  # your argparse/click entry point

    # Forward only the user args (skip the interpreter/module bits)
    # Expect `main(argv: list[str] | None)` to handle parsing and return
    # an int exit code or None (treated as 0 here).
    rv = main(sys.argv[1:])
    return int(rv) if isinstance(rv, int) else 0


if __name__ == "__main__":
    raise SystemExit(_entry())
