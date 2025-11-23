# src/gammapbh/__init__.py
"""
gammapbh — utilities to load BlackHawk tables, compose PBH gamma spectra,
and compute/plot spectra for various PBH mass distributions.

Mass in grams (g); Energy in MeV; Spectral density dN/dE [MeV^-1 s^-1]
"""
from importlib.metadata import version as _v

__version__ = _v("gammapbh")

# Re-export the programmatic functions so users can do:
#   from gammapbh import load_spectra_components, mass_function_lognormal, ...
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
