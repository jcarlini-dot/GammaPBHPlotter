"""
gammapbh — utilities to load BlackHawk tables, compose PBH gamma spectra,
and compute/plot spectra for various PBH mass distributions.

Units:
  Mass: grams (g)
  Energy: MeV
  Spectral density: dN/dE [MeV^-1 s^-1]
"""
from importlib.metadata import version as _v
__version__ = _v("gammapbh")

# Public programmatic API (import from .api so defaults are applied)
from .api import (
    discover_mass_folders,
    load_spectra_components,
    load_spectra_components_full,  # <-- add this
    mass_function,
    mass_function_exact,
    mass_function_lognormal,
)

__all__ = [
    "discover_mass_folders",
    "load_spectra_components",
    "load_spectra_components_full",  # <-- and list it here
    "mass_function",
    "mass_function_exact",
    "mass_function_lognormal",
    "__version__",
]
