"""Module entry-point for `python -m gammapbh`."""
from __future__ import annotations
import sys

HELP_FLAGS = {"-h", "--help", "/h", "/?"}

try:
    from .cli import main as cli_main  # zero-arg TUI entrypoint
except Exception:
    from . import main as cli_main  # type: ignore[attr-defined]

def _entry() -> int:
    if any(flag in sys.argv[1:] for flag in HELP_FLAGS):
        print("GammaPBHPlotter — PBH Spectrum Tool (v1.1.3)\n")
        print("Usage:\n  python -m gammapbh [--help]\n")
        print("Description:\n  Runs an interactive TUI to compute/plot PBH gamma-ray spectra")
        print("  from precomputed databases (BlackHawk primary/secondary + FSR + IFA).")
        return 0
    rv = cli_main()  # DO NOT pass argv; your main() is zero-arg
    return int(rv or 0)

if __name__ == "__main__":
    raise SystemExit(_entry())