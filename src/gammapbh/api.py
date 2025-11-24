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
    return Path(data_dir) if data_dir else _DATA_DIR

def _normalize_lists(masses_obj) -> Tuple[List[float], List[str]]:
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
    root = _resolve_data_dir(data_dir)
    try:
        return _discover_mass_folders(root)
    except TypeError:
        return _discover_mass_folders()

def load_spectra_components_full(mass: float | str, data_dir: Optional[str | Path] = None):
    """
    Return the *raw* return from the CLI loader (dict, tuple, etc.).
    """
    root = _resolve_data_dir(data_dir)
    masses_obj = discover_mass_folders(root)
    folder_name = _mass_to_folder_name(mass, masses_obj)
    mass_dir = root / str(folder_name)
    return _load_spectra_components(mass_dir)

def _resample(y: np.ndarray, x_old: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    return np.interp(x_new, x_old, y, left=0.0, right=0.0)

def _normalize_to_E_comps_from_dict(out: dict, align: str):
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
    Always return exactly (E, components_dict).
    align: 'secondary' (default), 'primary', 'union', or 'intersection'
    """
    out = load_spectra_components_full(mass, data_dir)
    return _normalize_to_E_comps(out, align=align)

# Pass-throughs
def mass_function(*args: Any, **kwargs: Any):
    return _mass_function(*args, **kwargs)

def mass_function_exact(*args: Any, **kwargs: Any):
    return _mass_function_exact(*args, **kwargs)

def mass_function_lognormal(*args: Any, **kwargs: Any):
    return _mass_function_lognormal(*args, **kwargs)