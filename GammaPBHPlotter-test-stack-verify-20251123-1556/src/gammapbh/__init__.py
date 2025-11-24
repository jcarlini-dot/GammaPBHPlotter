"""
gammapbh — utilities to load BlackHawk tables, compose PBH gamma spectra,
and plot/average spectra for various PBH mass distributions.

Units
-----
Mass: grams (g)
Energy: MeV
Spectral density: dN/dE [MeV^-1 s^-1]
"""

# Resolve package version with a safe fallback for dev/editable installs
try:
    from importlib.metadata import version as _pkg_version  # Py3.8+
except Exception:  # pragma: no cover
    _pkg_version = None  # very old Python; fallback below

try:
    __version__ = _pkg_version("gammapbh") if _pkg_version else "0+local"
except Exception:  # PackageNotFoundError etc. during local dev
    __version__ = "0+local"

# TEMPORARY: re-export non-interactive helpers from cli.py.
# (Recommended: move these into a new core.py and import from .core instead.)
from .cli import (
    discover_mass_folders,
    load_spectra_components,
    mass_function,
    mass_function_exact,
    mass_function_lognormal,
)

__all__ = [
    "discover_mass_folders",
    "load_spectra_components",
    "mass_function",
    "mass_function_exact",
    "mass_function_lognormal",
    "__version__",
]
