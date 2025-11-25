# api.py
"""
Public, importable API for the ``gammapbh`` package.

This module provides a stable, documented interface for:
  * Discovering which PBH masses are available in the bundled dataset.
  * Loading spectral components for a selected mass and returning them in a
    normalized, consistent shape: ``(E, components_dict)``.
  * Accessing mass-distribution helper functions (Gaussian, exact lists,
    log-normal, etc.) as simple pass-throughs.

Design
------
The actual I/O and composition logic lives in the underlying CLI/library
utilities (imported from ``.cli``). This module focuses on a clean, tidy
surface for library users:

* ``discover_mass_folders(data_dir=None)`` returns the available mass grid.
* ``load_spectra_components(mass, data_dir=None, align='secondary')`` always
  returns a 1D energy grid and a dict of aligned component arrays.
* Pass-throughs ``mass_function*`` expose mass distribution helpers.

Units
-----
Unless otherwise documented by the dataset, the conventional units are:

* Mass: grams (g)
* Energy: MeV
* Spectral density: dN/dE [MeV^-1 s^-1]
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Any, Tuple, List

import numpy as np

from .cli import (
    discover_mass_folders as _discover_mass_folders,
    load_spectra_components as _load_spectra_components,
    mass_function as _mass_function,
    mass_function_exact as _mass_function_exact,
    mass_function_lognormal as _mass_function_lognormal,
)

_PKG_DIR = Path(__file__).parent
_DATA_DIR = _PKG_DIR / "blackhawk_data"


def _resolve_data_dir(data_dir: Optional[str | Path]) -> Path:
    """
    Resolve the base data directory.

    Parameters
    ----------
    data_dir : str | pathlib.Path | None
        If provided, this path is used. Otherwise, the packaged data directory
        is returned.

    Returns
    -------
    pathlib.Path
        Resolved directory containing mass subfolders.
    """
    return Path(data_dir) if data_dir else _DATA_DIR


def _normalize_lists(masses_obj) -> Tuple[List[float], List[str]]:
    """
    Normalize a discovery result to parallel ``(numbers, names)`` lists.

    The discovery helper may yield either:
      * a 2-tuple ``(nums, names)``, or
      * a single list of numbers or strings.

    This function converts those shapes into a pair of ``List[float]`` (best
    effort; non-parsable names become ``nan``) and ``List[str]`` names.

    Parameters
    ----------
    masses_obj :
        Discovery result as returned by the underlying helper.

    Returns
    -------
    (List[float], List[str])
        Parallel lists: numeric masses (where possible) and string labels.
    """
    # Handle (nums, names) or just one list
    if isinstance(masses_obj, tuple) and len(masses_obj) == 2:
        a, b = masses_obj
        if a and isinstance(a[0], (int, float)):
            nums, names = list(map(float, a)), list(map(str, b))
        else:
            nums, names = list(map(float, b)), list(map(str, a))
        return nums, names
    if isinstance(masses_obj, list) and masses_obj:
        if isinstance(masses_obj[0], (int, float)):
            nums = list(map(float, masses_obj))
            names = [f"{x:.1e}".replace("E", "e") for x in nums]
            return nums, names
        if isinstance(masses_obj[0], str):
            names = list(map(str, masses_obj))
            nums = []
            for s in names:
                try:
                    nums.append(float(s))
                except Exception:
                    nums.append(float("nan"))
            return nums, names
    return [], []


def _mass_to_folder_name(mass: float | str, masses_obj) -> str:
    """
    Convert a user request (float or string) into a canonical folder name.

    The function consults the discovered grid to preserve exact names where
    available (e.g. ``"3.0e15"``). If not found, it formats the numeric value.

    Parameters
    ----------
    mass : float | str
        Requested mass; may be a float or string (e.g. ``"3e15"``).
    masses_obj :
        Discovery result (tuple or list) used to look up canonical names.

    Returns
    -------
    str
        Canonical folder name for the requested mass.
    """
    nums, names = _normalize_lists(masses_obj)
    if isinstance(mass, str):
        if mass in names:
            return mass
        try:
            m = float(mass)
            return f"{m:.1e}".replace("E", "e")
        except Exception:
            return mass
    try:
        idx = nums.index(float(mass))
        return names[idx]
    except Exception:
        return f"{float(mass):.1e}".replace("E", "e")


# ---------- Public API ----------

def discover_mass_folders(data_dir: Optional[str | Path] = None):
    """
    Discover the available PBH mass grid from the dataset.

    Parameters
    ----------
    data_dir : str | pathlib.Path | None
        Optional override for the data root. If ``None`` (default), the
        built-in packaged data path is used.

    Returns
    -------
    object
        Whatever the underlying helper returns (a list or a ``(nums, names)``
        tuple). Use with ``_normalize_lists`` if you need parallel lists.
    """
    root = _resolve_data_dir(data_dir)
    try:
        return _discover_mass_folders(root)
    except TypeError:
        # Backward-compat: older helper may not accept an explicit root
        return _discover_mass_folders()


def load_spectra_components_full(mass: float | str, data_dir: Optional[str | Path] = None):
    """
    Load the *raw* spectral data structure for a given mass.

    This returns the unmodified output of the lower-level loader (dict/tuple),
    which may contain separate primary/secondary energy axes and component
    arrays on their native grids.

    Parameters
    ----------
    mass : float | str
        Target mass (e.g. ``3.0e15``). May be numeric or a folder-style string.
    data_dir : str | pathlib.Path | None
        Optional override for the data root.

    Returns
    -------
    object
        The underlying loader's native output (implementation detail).
    """
    root = _resolve_data_dir(data_dir)
    masses_obj = discover_mass_folders(root)
    folder_name = _mass_to_folder_name(mass, masses_obj)
    mass_dir = root / str(folder_name)
    return _load_spectra_components(mass_dir)


def _resample(y: np.ndarray, x_old: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    """
    1D linear interpolation helper with zero fill outside the original domain.

    Parameters
    ----------
    y : numpy.ndarray
        Values on the original grid ``x_old``.
    x_old : numpy.ndarray
        Original, strictly increasing grid.
    x_new : numpy.ndarray
        Target grid.

    Returns
    -------
    numpy.ndarray
        Values resampled onto ``x_new``.
    """
    return np.interp(x_new, x_old, y, left=0.0, right=0.0)


def _normalize_to_E_comps_from_dict(out: dict, align: str):
    """
    Normalize a dict-shaped loader output to a single energy grid + components.

    Expected keys are ``"energy_primary"`` and/or ``"energy_secondary"`` with
    component arrays suffixed by ``"_primary"`` or ``"_secondary"``. The
    selected alignment is applied to return a uniform energy grid and aligned
    component arrays.

    If those explicit keys are absent, the function tries common fallbacks
    (``"E"``, ``"energy"``, ``"energies"``, ``"E_MeV"``) and assumes all
    components are already on that grid.

    Parameters
    ----------
    out : dict
        Loader output.
    align : {'primary','secondary','union','intersection'}
        Energy alignment rule.

    Returns
    -------
    (numpy.ndarray, dict[str, numpy.ndarray])
        ``(E, components_dict)`` with any of ``total_primary``, ``total_secondary``,
        and ``total`` synthesized when possible.
    """
    # Expect explicit energy axes
    if ("energy_primary" not in out) or ("energy_secondary" not in out):
        # Try generic keys as fallback
        for k in ("E", "energy", "energies", "E_MeV"):
            if k in out:
                E = np.asarray(out[k], dtype=float)
                comps = {kk: np.asarray(v, dtype=float) for kk, v in out.items() if kk != k}
                # Uniform lengths check
                for kk, vv in comps.items():
                    if vv.shape != E.shape:
                        comps[kk] = _resample(vv, np.linspace(E.min(), E.max(), vv.size), E)
                return E, comps
        raise ValueError("Dictionary output lacks recognizable energy keys.")

    Ep = np.asarray(out["energy_primary"], dtype=float)
    Es = np.asarray(out["energy_secondary"], dtype=float)

    if align == "primary":
        E = Ep
    elif align == "secondary":
        E = Es
    elif align == "union":
        E = np.union1d(Ep, Es)
    elif align == "intersection":
        # true intersection (not union!) — sorted, unique
        E = np.intersect1d(Ep, Es)
    else:
        raise ValueError(f"Unknown align='{align}'. Use 'primary','secondary','union','intersection'.")

    comps: dict[str, np.ndarray] = {}
    for k, v in out.items():
        if k in ("energy_primary", "energy_secondary"):
            continue
        arr = np.asarray(v, dtype=float)
        # Decide source grid by suffix
        if k.endswith("_primary"):
            src = Ep
        elif k.endswith("_secondary"):
            src = Es
        else:
            # Heuristic: pick the matching length; else assume target already
            if arr.shape == Ep.shape:
                src = Ep
            elif arr.shape == Es.shape:
                src = Es
            else:
                src = E
        comps[k] = arr if src is E else _resample(arr, src, E)

    # Totals (if prim/sec components exist)
    prim_keys = [k for k in comps if k.endswith("_primary")]
    sec_keys  = [k for k in comps if k.endswith("_secondary")]
    if prim_keys:
        comps["total_primary"] = np.sum([comps[k] for k in prim_keys], axis=0)
    if sec_keys:
        comps["total_secondary"] = np.sum([comps[k] for k in sec_keys], axis=0)
    if prim_keys and sec_keys:
        comps["total"] = comps["total_primary"] + comps["total_secondary"]
    elif "total" not in comps and (prim_keys or sec_keys):
        comps["total"] = comps.get("total_primary", np.zeros_like(E)) + comps.get("total_secondary", np.zeros_like(E))

    return E, comps


def _normalize_to_E_comps(out, align: str):
    """
    Normalize any supported loader output to the canonical ``(E, comps)`` form.

    Parameters
    ----------
    out : object
        Either a tuple/list ``(E, components_dict)`` or a dict as described in
        ``_normalize_to_E_comps_from_dict``.
    align : {'primary','secondary','union','intersection'}
        Energy alignment rule. Only used for the dict case.

    Returns
    -------
    (numpy.ndarray, dict[str, numpy.ndarray])
        Canonical tuple ``(E, components_dict)``.
    """
    # tuple/list case
    if isinstance(out, (tuple, list)) and len(out) >= 2:
        E, comps = out[0], out[1]
        return np.asarray(E, dtype=float), {k: np.asarray(v, dtype=float) for k, v in comps.items()}
    # dict case
    if isinstance(out, dict):
        return _normalize_to_E_comps_from_dict(out, align)
    raise TypeError(f"Could not normalize loader output to (E, components). Got: {type(out)!r}")


def load_spectra_components(mass: float | str, data_dir: Optional[str | Path] = None, align: str = "secondary"):
    """
    Load spectra for ``mass`` and **always** return ``(E, components_dict)``.

    This is the user-friendly, shape-stable loader. It applies the requested
    ``align`` rule to reconcile primary/secondary grids when necessary.

    Parameters
    ----------
    mass : float | str
        Target mass (e.g. ``3.0e15``). May be numeric or a folder-style string.
    data_dir : str | pathlib.Path | None
        Optional override for the data root.
    align : {'secondary','primary','union','intersection'}, default 'secondary'
        Energy alignment rule applied when the underlying data expose separate
        grids.

    Returns
    -------
    (numpy.ndarray, dict[str, numpy.ndarray])
        Energy grid and component arrays aligned to that grid.
    """
    out = load_spectra_components_full(mass, data_dir)
    return _normalize_to_E_comps(out, align=align)


# Pass-throughs

def mass_function(*args: Any, **kwargs: Any):
    """
    Pass-through to the underlying mass function helper.

    Parameters
    ----------
    *args, **kwargs :
        Forwarded verbatim to the implementation in ``.cli``.

    Returns
    -------
    Any
        Whatever the underlying helper returns.
    """
    return _mass_function(*args, **kwargs)


def mass_function_exact(*args: Any, **kwargs: Any):
    """
    Pass-through to the exact (discrete) mass function helper.

    Parameters
    ----------
    *args, **kwargs :
        Forwarded verbatim to the implementation in ``.cli``.

    Returns
    -------
    Any
        Whatever the underlying helper returns.
    """
    return _mass_function_exact(*args, **kwargs)


def mass_function_lognormal(*args: Any, **kwargs: Any):
    """
    Pass-through to the log-normal mass function helper.

    Parameters
    ----------
    *args, **kwargs :
        Forwarded verbatim to the implementation in ``.cli``.

    Returns
    -------
    Any
        Whatever the underlying helper returns.
    """
    return _mass_function_lognormal(*args, **kwargs)
