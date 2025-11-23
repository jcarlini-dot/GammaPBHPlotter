"""Entry-point for `python -m gammapbh` with proper --help behavior."""
from __future__ import annotations
import sys

HELP_FLAGS = {"-h", "--help", "/h", "/?"}
HELP_TEXT = """\
GammaPBHPlotter — PBH Spectrum Tool (v1.1.3)

Usage:
  python -m gammapbh [--help]

Description:
  Runs an interactive TUI to compute/plot PBH gamma-ray spectra
  from precomputed databases (BlackHawk primary/secondary + FSR + IFA).
"""

try:
    from .cli import main as cli_main  # launches TUI
except Exception:
    from . import main as cli_main      # type: ignore[attr-defined]

def _entry() -> int:
    if any(flag in sys.argv[1:] for flag in HELP_FLAGS):
        print(HELP_TEXT); return 0
    rv = cli_main()
    return int(rv or 0)

if __name__ == "__main__":
    raise SystemExit(_entry())