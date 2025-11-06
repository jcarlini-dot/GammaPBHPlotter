# src/gammapbh/__init__.py
"""
gammapbh — utilities to load BlackHawk tables, compose PBH gamma spectra,
and plot/average spectra for various PBH mass distributions.

Units
-----
Mass: grams (g)
Energy: MeV
Spectral density: dN/dE [MeV^-1 s^-1]
"""

from importlib.metadata import version as _v

__all__ = [
    "load_spectra_components",
    "discover_mass_folders",
    "mass_function",
    "mass_function_exact",
    "mass_function_lognormal",
    # add other public functions/classes you want users to import
]
__version__ = _v("gammapbh")

# Re-export key APIs from submodules (so docs & users can do `from gammapbh import ...`)
from .cli import (  # or better: move non-CLI helpers to a library module and import from there
    discover_mass_folders,
    load_spectra_components,
    mass_function,
    mass_function_exact,
    mass_function_lognormal,
)

__all__ = ["cli"]
__version__ = "1.1.2"


