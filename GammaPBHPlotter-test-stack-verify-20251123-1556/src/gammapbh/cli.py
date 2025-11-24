# src/gammapbh/cli.py
"""
GammaPBHPlotter — interactive CLI to analyze and visualize Hawking-radiation
gamma-ray spectra of primordial black holes (PBHs).

This module provides:
  - Monochromatic spectra visualization for selected PBH masses.
  - Distributed spectra from physically motivated mass PDFs:
      - Gaussian collapse (Press–Schechter–like).
      - Non-Gaussian collapse (Biagetti et al. formulation).
      - Log-normal mass function.
  - A custom-equation mass PDF tool that lets users enter f(m) directly.
  - A viewer for previously saved runs (with spectrum overlays and
    per-selection mass histograms, including analytic/KDE overlays).

All user-facing plotting is log–log with stable zero-flooring in linear space
to avoid numerical warnings. Interpolations are performed in (logM, logE) space
with linear/cubic bivariate splines, and inflight-annihilation tails are
sanity-trimmed to prevent staircase artifacts in the rightmost bins.

Conventions
-----------
- Masses are in grams [g].
- Energies are in MeV.
- Spectra are per energy [MeV^-1 s^-1].
- “E² dN/dE” overlays are used for SED-style views.
"""

from __future__ import annotations

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.special import erf
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import trapezoid
from types import SimpleNamespace
from colorama import Fore, Style

try:
    # works when invoked as `python -m gammapbh.cli` or via installed entry-point
    from . import __version__  # type: ignore
except Exception:
    try:
        # works when invoked as a plain script `python path/to/cli.py`
        from gammapbh import __version__  # type: ignore
    except Exception:
        __version__ = "dev"

def pause(msg="Press Enter to continue…"):
    # Only pause if running interactively
    if sys.stdin.isatty():
        input(msg)

# … then in view_previous_spectra(), keep your `pause()` calls unchanged.
# Under pytest (non-tty), pause() becomes a no-op and won’t consume feeder in

# ---------------------------
# Matplotlib/NumPy basics
# ---------------------------
plt.rcParams.update({'font.size': 12})
# Suppress harmless warnings when we intentionally clamp underflows to ~0
np.seterr(divide='ignore', invalid='ignore')


# ---------------------------
# Paths (package-internal only)
# ---------------------------
def _resolve_data_dir() -> str:
    """
    Resolve the *package-internal* data directory that contains BlackHawk tables.

    Returns
    -------
    str
        Absolute path to the packaged `blackhawk_data` directory.

    Notes
    -----
    - We do not permit user-provided paths here; reproducibility requires the
      tables bundled with the installed package to be used.
    """
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(pkg_dir, "blackhawk_data")


def _resolve_results_root() -> str:
    """
    Resolve the *package-internal* results directory used for all outputs.

    Returns
    -------
    str
        Absolute path to the packaged `results` directory.

    Raises
    ------
    RuntimeError
        If the directory cannot be created or written to.
    """
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    dest = os.path.join(pkg_dir, "results")
    os.makedirs(dest, exist_ok=True)
    # Writability quick check: create and remove a tiny temp file
    try:
        test = os.path.join(dest, ".writetest.tmp")
        with open(test, "w") as fh:
            fh.write("ok")
        os.remove(test)
    except Exception as e:
        raise RuntimeError(f"Results directory is not writable: {dest}\n{e}")
    return dest


DATA_DIR     = _resolve_data_dir()
RESULTS_DIR  = _resolve_results_root()

MONO_RESULTS_DIR   = os.path.join(RESULTS_DIR, "monochromatic")
CUSTOM_RESULTS_DIR = os.path.join(RESULTS_DIR, "custom_equation")
GAUSS_RESULTS_DIR  = os.path.join(RESULTS_DIR, "gaussian")
NGAUSS_RESULTS_DIR = os.path.join(RESULTS_DIR, "non_gaussian")
LOGN_RESULTS_DIR   = os.path.join(RESULTS_DIR, "lognormal")

for _d in (MONO_RESULTS_DIR, CUSTOM_RESULTS_DIR, GAUSS_RESULTS_DIR, NGAUSS_RESULTS_DIR, LOGN_RESULTS_DIR):
    os.makedirs(_d, exist_ok=True)


# ---------------------------
# Labels
# ---------------------------
GAUSSIAN_METHOD     = "Gaussian collapse"
NON_GAUSSIAN_METHOD = "Non-Gaussian Collapse"
LOGNORMAL_METHOD    = "Log-Normal Distribution"


# ---------------------------
# Required files within each mass folder
# ---------------------------
REQUIRED_FILES = [
    "instantaneous_primary_spectra.txt",
    "instantaneous_secondary_spectra.txt",
    "inflight_annihilation_prim.txt",
    "inflight_annihilation_sec.txt",
    "final_state_radiation_prim.txt",
    "final_state_radiation_sec.txt",
]


# ---------------------------
# Back navigation support
# ---------------------------
class BackRequested(Exception):
    """Raised when the user enters 'b' or 'back' to return to the prior screen."""
    pass


# ---------------------------
# Discovery helpers
# ---------------------------
def discover_mass_folders(data_dir: str) -> tuple[list[float], list[str]]:
    """
    Discover valid mass folders within `data_dir` that contain all required files.

    Parameters
    ----------
    data_dir : str
        Absolute or relative path to the BlackHawk data directory.

    Returns
    -------
    (list[float], list[str])
        A pair (masses, names) sorted by mass, where `masses[i]` corresponds to
        directory name `names[i]`.

    Notes
    -----
    - Folders are expected to be named as a float mass in grams (e.g., "1.00e+16").
    - Only folders containing the full REQUIRED_FILES set are returned.
    """
    masses, names = [], []
    try:
        for name in os.listdir(data_dir):
            p = os.path.join(data_dir, name)
            if not os.path.isdir(p):
                continue
            try:
                m = float(name)
            except ValueError:
                continue
            if all(os.path.isfile(os.path.join(p, f)) for f in REQUIRED_FILES):
                masses.append(m); names.append(name)
    except FileNotFoundError:
        return [], []
    if not masses:
        return [], []
    order = np.argsort(masses)
    return [float(masses[i]) for i in order], [names[i] for i in order]


# ---------------------------
# CLI + parsing helpers
# ---------------------------
def info(msg: str) -> None:
    """Print an informational (cyan) line."""
    print(Fore.CYAN + "ℹ " + msg + Style.RESET_ALL)


def warn(msg: str) -> None:
    """Print a warning (yellow) line."""
    print(Fore.YELLOW + "⚠ " + msg + Style.RESET_ALL)


def err(msg: str) -> None:
    """Print an error (red) line."""
    print(Fore.RED + "✖ " + msg + Style.RESET_ALL)


def user_input(
    prompt: str,
    *,
    allow_back: bool = False,
    allow_exit: bool = True,
    allow_comment: bool = True,
) -> str:
    """
    Wrapper around :func:`input` that also understands navigation and, optionally,
    ignores comment-style lines.

    Parameters
    ----------
    prompt : str
        Text to display when requesting input from the user.

    allow_back : bool, optional
        If True, the special commands ``"b"`` or ``"back"`` (case-insensitive)
        will *not* be returned as normal input. Instead, a :class:`BackRequested`
        exception is raised, allowing the caller to return to a previous menu.

    allow_exit : bool, optional
        If True, the special commands ``"q"`` or ``"exit"`` (case-insensitive)
        cause the program to terminate immediately via :func:`sys.exit`. This is
        used for global “quit” shortcuts in the interactive CLI.

    allow_comment : bool, optional
        If True (default), any line whose first non-whitespace character is
        ``'#'`` is treated as a comment and ignored. This is specifically to make
        copy–pasted transcripts from the documentation or README work smoothly:
        users can paste sequences like::

            gammapbh
            # Example A — Monochromatic spectra
            # 1) Select "Monochromatic spectra"
            1
            # 2) Enter PBH masses (g)
            3.14e15, 1.4e14

        and all ``'# …'`` lines will be skipped until a non-comment line is read.

    Returns
    -------
    str
        The raw user input (with leading/trailing whitespace stripped) that is
        *not* a navigation/exit command and is not suppressed as a comment.

    Raises
    ------
    BackRequested
        If the user enters ``"b"`` or ``"back"`` and ``allow_back`` is True.

    SystemExit
        If the user enters ``"q"`` or ``"exit"`` and ``allow_exit`` is True.
        This is the normal way to terminate the top-level CLI loop.

    Notes
    -----
    - All comparisons for navigation/exit/comment detection are done on the
      raw line **before** it is returned to the caller.
    - The function loops until it encounters a line that is not treated as a
      comment or navigation/exit command, so callers always receive a
      meaningful payload string.
    """
    while True:
        txt = input(prompt)

        # Treat comment-only lines as no-ops so that README examples with '#'
        # can be copy–pasted directly into an interactive session.
        if allow_comment and txt.lstrip().startswith("#"):
            continue

        txt = txt.strip()
        low = txt.lower()

        if allow_exit and low in ("exit", "q"):
            print("Exiting software.")
            sys.exit(0)

        if allow_back and low in ("b", "back"):
            raise BackRequested()

        return txt



def list_saved_runs(base_dir: str) -> list[str]:
    """
    List child directories beneath `base_dir`.

    Parameters
    ----------
    base_dir : str
        Root directory containing saved runs.

    Returns
    -------
    list[str]
        Sorted child directory names (no files).
    """
    try:
        return sorted(d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)))
    except FileNotFoundError:
        return []


def snap_to_available(mval: float, available: list[float], tol: float = 1e-12) -> float | None:
    """
    If `mval` is essentially equal (in log space) to one of `available` masses,
    return that available mass. Otherwise return None.

    Parameters
    ----------
    mval : float
        Desired mass value [g].
    available : list[float]
        Pre-rendered masses available.
    tol : float
        Allowed absolute difference in ln-space for a snap.

    Returns
    -------
    float or None
    """
    if not available:
        return None
    log_m = np.log(mval)
    log_available = np.log(np.array(available))
    diffs = np.abs(log_available - log_m)
    idx = np.argmin(diffs)
    return available[idx] if diffs[idx] < tol else None


def parse_float_list_verbose(
    s: str,
    *,
    name: str = "value",
    bounds: tuple[float | None, float | None] | None = None,
    allow_empty: bool = False,
    positive_only: bool = False,
    strict_gt: bool = False,
    strict_lt: bool = False,
) -> list[float]:
    """
    Parse a comma-separated list of floats with verbose validation.

    Parameters
    ----------
    s : str
        Input string, e.g. "1e15, 2e15".
    name : str
        Friendly name used in warning messages.
    bounds : (float|None, float|None) or None
        Inclusive (lo, hi) bounds if provided.
    allow_empty : bool
        If False and parsing yields nothing, a warning is printed.
    positive_only : bool
        If True, keep only values > 0.
    strict_gt : bool
        If True and bounds[0] not None: enforce v > lo; else v ≥ lo.
    strict_lt : bool
        If True and bounds[1] not None: enforce v < hi; else v ≤ hi.

    Returns
    -------
    list[float]
        Validated, de-duplicated floats (first occurrence kept).
    """
    if (s is None or s.strip() == ""):
        if not allow_empty:
            warn(f"No {name}s provided.")
        return []
    vals, seen = [], set()
    lo, hi = (bounds or (None, None))
    for tok in s.split(","):
        t = tok.strip()
        if not t:
            continue
        try:
            v = float(t)
        except Exception:
            warn(f"Skipping token '{t}': {name} is not a valid number.")
            continue
        if positive_only and v <= 0:
            warn(f"Skipping {name} {v:g}: must be > 0.")
            continue
        if lo is not None:
            if (strict_gt and not (v > lo)) or (not strict_gt and not (v >= lo)):
                cmp = ">" if strict_gt else "≥"
                warn(f"Skipping {name} {v:g}: must be {cmp} {lo:g}.")
                continue
        if hi is not None:
            if (strict_lt and not (v < hi)) or (not strict_lt and not (v <= hi)):
                cmp = "<" if strict_lt else "≤"
                warn(f"Skipping {name} {v:g}: must be {cmp} {hi:g}.")
                continue
        if v in seen:
            warn(f"Duplicate {name} {v:g}: keeping first, skipping this one.")
            continue
        vals.append(v); seen.add(v)
    if not vals and not allow_empty:
        warn(f"No usable {name}s parsed.")
    return vals


# ---------------------------
# PDFs (collapse space)
# ---------------------------
def delta_l(mass_ratio: np.ndarray, kappa: float, delta_c: float, gamma: float) -> np.ndarray:
    """
    Convert mass ratio to the linear threshold δ_l used in collapse models.

    Parameters
    ----------
    mass_ratio : ndarray
        Dimensionless M/M_peak (or equivalent model-specific scaling).
    kappa : float
    delta_c : float
    gamma : float

    Returns
    -------
    ndarray
        δ_l(mass_ratio) with the analytic mapping and a safe clip for the sqrt argument.
    """
    y = (mass_ratio / kappa)**(1.0 / gamma)
    arg = 64 - 96 * (delta_c + y)
    arg = np.clip(arg, 0.0, None)
    return (8 - np.sqrt(arg)) / 6


def mass_function(delta_l_val: np.ndarray, sigma_x: float, delta_c: float, gamma: float) -> np.ndarray:
    """
    Gaussian-collapse proxy mass function in δ_l-space.

    Parameters
    ----------
    delta_l_val : ndarray
        δ_l grid.
    sigma_x : float
        Collapse dispersion parameter.
    delta_c : float
        Critical collapse threshold.
    gamma : float
        Shape parameter.

    Returns
    -------
    ndarray
        Unnormalized mass function (shape same as input).
    """
    term1 = 1.0 / (np.sqrt(2 * np.pi) * sigma_x)
    term2 = np.exp(-delta_l_val**2 / (2 * sigma_x**2))
    term3 = delta_l_val - (3/8) * delta_l_val**2 - delta_c
    term4 = gamma * np.abs(1 - (3/4) * delta_l_val)
    return term1 * term2 * term3 / term4


def mass_function_exact(
    delta_l_val: np.ndarray,
    sigma_X: float,
    sigma_Y: float,
    delta_c: float,
    gamma: float
) -> np.ndarray:
    """
    Non-Gaussian mass function (Biagetti et al., Eq. 20 shape—up to constants),
    mapped into δ_l with a Jacobian consistent with the collapse mapping used above.

    Parameters
    ----------
    delta_l_val : ndarray
        δ_l grid.
    sigma_X : float
        Dispersion along X-direction.
    sigma_Y : float
        Dispersion along Y-direction (often tied to sigma_X via ratio).
    delta_c : float
        Critical threshold.
    gamma : float
        Shape parameter for the mapping.

    Returns
    -------
    ndarray
        Unnormalized mass function (shape same as input).
    """
    A = sigma_X**2 + (sigma_Y * delta_l_val)**2
    exp_pref = np.exp(-1.0 / (2.0 * sigma_Y**2))
    term1 = 2.0 * sigma_Y * np.sqrt(A)
    inner_exp = np.exp(sigma_X**2 / (2.0 * sigma_Y**2 * (sigma_X**2 + 2.0 * (sigma_Y * delta_l_val)**2)))
    erf_arg = sigma_X * np.sqrt(2.0) / np.sqrt(A)  # stable
    term2 = np.sqrt(2.0 * np.pi) * sigma_X * inner_exp * erf(erf_arg)
    bracket = term1 + term2
    norm = exp_pref * sigma_X / (2.0 * np.pi * A**1.5)
    jacobian = ((delta_l_val - 0.375 * delta_l_val**2 - delta_c) /
                (gamma * np.abs(1.0 - 0.75 * delta_l_val)))
    return norm * bracket * jacobian


def mass_function_lognormal(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """
    Standard log-normal PDF in variable x.

    Parameters
    ----------
    x : ndarray
        Positive support (will be clipped below to avoid divide-by-zero).
    mu : float
        Mean in ln-space.
    sigma : float
        Std. dev. in ln-space (must be > 0).

    Returns
    -------
    ndarray
        Log-normal PDF values at x.
    """
    x_clipped = np.clip(x, 1e-16, None)
    return (1.0 / (x_clipped * sigma * np.sqrt(2 * np.pi))
            * np.exp(- (np.log(x_clipped) - mu)**2 / (2 * sigma**2)))


# ---------------------------
# Data loaders
# ---------------------------
def load_data(filepath: str, skip_header: int = 0) -> np.ndarray:
    """
    A thin wrapper around `numpy.genfromtxt` with explicit header skipping.

    Parameters
    ----------
    filepath : str
        Path to file.
    skip_header : int
        Number of header lines to skip.

    Returns
    -------
    ndarray
        Parsed numeric array.

    Raises
    ------
    FileNotFoundError
        If file does not exist.
    ValueError
        If `genfromtxt` fails due to column inconsistency.
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    return np.genfromtxt(filepath, skip_header=skip_header)


def load_xy_lenient(filepath: str, skip_header: int = 0, min_cols: int = 2) -> np.ndarray:
    """
    Robustly load at least two numeric columns from a whitespace/CSV-like text file,
    skipping blank lines, comment lines, and any lines with fewer than `min_cols` tokens.

    This specifically fixes files where the first data row contains a single integer
    (e.g., a length or counter), followed by proper 2-column numeric rows; vanilla
    `genfromtxt` would lock onto the one-column width and then error.

    Parameters
    ----------
    filepath : str
        Path to the file to read.
    skip_header : int, optional
        Number of initial lines to skip unconditionally.
    min_cols : int, optional
        Minimum number of numeric columns required to accept a line (default 2).

    Returns
    -------
    ndarray
        Array of shape (N, >=min_cols). Only the first `min_cols` columns are guaranteed.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If no usable numeric rows are found.

    Notes
    -----
    - Treats lines starting with '#' as comments.
    - Replaces commas with spaces to tolerate CSV-ish files.
    - Silently skips lines that fail float conversion or are too short.
    """
    rows = []
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        for i, raw in enumerate(fh):
            if i < skip_header:
                continue
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            line = line.replace(",", " ")
            parts = [p for p in line.split() if p]
            if len(parts) < min_cols:
                continue
            try:
                nums = [float(parts[j]) for j in range(min_cols)]
            except Exception:
                continue
            rows.append(nums)
    if not rows:
        raise ValueError(f"No usable numeric rows with ≥{min_cols} columns in {filepath}")
    return np.asarray(rows, dtype=float)


def load_spectra_components(directory: str) -> dict[str, np.ndarray]:
    """
    Load and align spectral components for a given mass-directory.

    Files expected in `directory`
    -----------------------------
    instantaneous_primary_spectra.txt
        Columns: E(GeV)  dN/dE (GeV^-1 s^-1)  [we later convert E to MeV and flux to MeV^-1 s^-1]
    instantaneous_secondary_spectra.txt
        Columns: E(MeV)  dN/dE (MeV^-1 s^-1)
    inflight_annihilation_prim.txt
        Typically two columns (E(MeV), rate) but may contain a spurious single-number line first.
    inflight_annihilation_sec.txt
        Same caveat as above.
    final_state_radiation_prim.txt
        Typically two columns, sometimes with one header line to skip.
    final_state_radiation_sec.txt
        Typically two columns, sometimes with one header line to skip.

    Returns
    -------
    dict[str, ndarray]
        Keys:
            energy_primary, energy_secondary,
            direct_gamma_primary, direct_gamma_secondary,
            IFA_primary, IFA_secondary,
            FSR_primary, FSR_secondary

    Notes
    -----
    - This function now uses `load_xy_lenient` for IFA/FSR files to survive files with
      leading single-value rows.
    """
    primary   = load_data(os.path.join(directory, "instantaneous_primary_spectra.txt"),   skip_header=2)[123:]
    secondary = load_data(os.path.join(directory, "instantaneous_secondary_spectra.txt"), skip_header=1)

    # lenient loads for files that sometimes start with a single-number row
    IFA_prim  = load_xy_lenient(os.path.join(directory, "inflight_annihilation_prim.txt"))
    IFA_sec   = load_xy_lenient(os.path.join(directory, "inflight_annihilation_sec.txt"))
    FSR_prim  = load_xy_lenient(os.path.join(directory, "final_state_radiation_prim.txt"), skip_header=1)
    FSR_sec   = load_xy_lenient(os.path.join(directory, "final_state_radiation_sec.txt"),  skip_header=1)

    E_prim = primary[:, 0] * 1e3  # convert GeV → MeV
    E_sec  = secondary[:, 0]      # already in MeV

    return {
        'energy_primary':         E_prim,
        'energy_secondary':       E_sec,
        'direct_gamma_primary':   primary[:, 1] / 1e3,  # GeV^-1 → MeV^-1
        'direct_gamma_secondary': secondary[:, 1],
        'IFA_primary':            np.interp(E_prim, IFA_prim[:, 0], IFA_prim[:, 1], left=0.0, right=0.0),
        'IFA_secondary':          np.interp(E_sec,  IFA_sec[:, 0],  IFA_sec[:, 1],  left=0.0, right=0.0),
        'FSR_primary':            np.interp(E_prim, FSR_prim[:, 0], FSR_prim[:, 1]),
        'FSR_secondary':          np.interp(E_sec,  FSR_sec[:, 0],  FSR_sec[:, 1]),
    }


# ---------------------------
# Monochromatic
# ---------------------------
def generate_monochromatic_for_mass(target_mass: float, data_dir: str, out_dir: str) -> str:
    """
    Generate (or more precisely, assemble) a monochromatic spectrum file for the
    nearest available pre-rendered mass to `target_mass`.

    Parameters
    ----------
    target_mass : float
        Desired PBH mass [g].
    data_dir : str
        Directory containing the BlackHawk mass folders.
    out_dir : str
        Directory to write the output TXT file.

    Returns
    -------
    str
        Path to the saved monochromatic spectrum file. Columns:
            E_gamma(MeV), TotalSpectrum(MeV^-1 s^-1)

    Notes
    -----
    - We re-compute the *total* as Direct + Secondary + IFA + FSR aligned onto the
      primary energy grid for consistency with plotting routines.
    - The output file name encodes the requested mass, not the snapped mass, to
      reflect the user's intention; the data inside reflects the snapped folder.
    """
    masses, names = discover_mass_folders(data_dir)
    if not masses:
        raise RuntimeError("No valid mass folders found to generate monochromatic spectrum.")
    snap = snap_to_available(target_mass, masses)
    if snap is None:
        # choose nearest in log-space
        log_t = np.log(target_mass)
        idx = int(np.argmin(np.abs(np.log(masses) - log_t)))
        snap = masses[idx]
    idx_snap = np.where(np.isclose(masses, snap, rtol=0, atol=0))[0][0]
    sub = os.path.join(data_dir, names[idx_snap])
    S = load_spectra_components(sub)

    # align everything on the primary grid
    E = S['energy_primary']
    total = (
        S['direct_gamma_primary']
        + np.interp(E, S['energy_secondary'], S['direct_gamma_secondary'], left=0, right=0)
        + S['IFA_primary'] + np.interp(E, S['energy_secondary'], S['IFA_secondary'], left=0, right=0)
        + S['FSR_primary'] + np.interp(E, S['energy_secondary'], S['FSR_secondary'], left=0, right=0)
    )

    out_name = os.path.join(out_dir, f"{target_mass:.2e}_mono_generated.txt")
    np.savetxt(out_name, np.column_stack((E, total)),
               header="E_gamma(MeV)   TotalSpectrum (MeV^-1 s^-1)", fmt="%.10e")
    return out_name


def monochromatic_spectra() -> None:
    """
    Interactive tool to plot one or more monochromatic spectra.

    Flow
    ----
    1) Discover available pre-rendered masses.
    2) Ask user to enter a comma-separated list of target masses.
    3) Build logM–logE splines across the full mass grid to allow interpolation
       for off-grid masses as needed.
    4) For each requested mass, plot component curves and total; then offer to
       save selected spectra into the monochromatic results folder.
    """
    masses, names = discover_mass_folders(DATA_DIR)
    if not masses:
        warn(f"No valid mass folders found under: {DATA_DIR}")
        return
    MIN_MASS, MAX_MASS = min(masses), max(masses)

    try:
        masses_str = user_input(
            f"Enter PBH masses (g) to simulate (comma-separated; allowed range [{MIN_MASS:.2e}, {MAX_MASS:.2e}]): ",
            allow_back=True
        )
    except BackRequested:
        return

    mass_list = []
    if masses_str.strip():
        for tok in masses_str.split(','):
            t = tok.strip()
            if not t:
                continue
            try:
                mval = float(t)
            except Exception:
                warn(f"Skipping mass token '{t}': not a number.")
                continue
            if not (MIN_MASS <= mval <= MAX_MASS):
                warn(f"Skipping mass {mval:.3e} g: outside allowed range [{MIN_MASS:.2e}, {MAX_MASS:.2e}].")
                continue
            mass_list.append(mval)
    if not mass_list:
        warn("No valid masses provided. Returning to menu.")
        return

    info("Pre-loading pre-rendered components …")
    first_S = load_spectra_components(os.path.join(DATA_DIR, names[0]))
    E_ref   = first_S['energy_primary']
    N_E     = len(E_ref)
    N_M     = len(masses)

    direct_mat     = np.zeros((N_M, N_E))
    secondary_mat  = np.zeros((N_M, N_E))
    inflight_mat   = np.zeros((N_M, N_E))
    finalstate_mat = np.zeros((N_M, N_E))
    Emax_ifa       = np.zeros(N_M)

    for i, m in enumerate(masses):
        sub = os.path.join(DATA_DIR, names[i])
        S = load_spectra_components(sub)
        direct_mat[i]     = S['direct_gamma_primary']
        secondary_mat[i]  = np.interp(E_ref, S['energy_secondary'], S['direct_gamma_secondary'], left=0, right=0)
        inflight_mat[i]   = S['IFA_primary'] + np.interp(E_ref, S['energy_secondary'], S['IFA_secondary'], left=0, right=0)
        finalstate_mat[i] = S['FSR_primary'] + np.interp(E_ref, S['energy_secondary'], S['FSR_secondary'], left=0, right=0)
        p = load_xy_lenient(os.path.join(sub, "inflight_annihilation_prim.txt"))
        s = load_xy_lenient(os.path.join(sub, "inflight_annihilation_sec.txt"))
        Emax_ifa[i] = max(p[:,0].max() if p.size else 0, s[:,0].max() if s.size else 0)

    logM_all = np.log(masses)
    logE     = np.log(E_ref)
    tiny     = 1e-300

    ld = np.log(np.where(direct_mat>tiny,     direct_mat,     tiny))
    ls = np.log(np.where(secondary_mat>tiny,  secondary_mat,  tiny))
    li = np.log(np.where(inflight_mat>tiny,   inflight_mat,   tiny))
    lf = np.log(np.where(finalstate_mat>tiny, finalstate_mat, tiny))

    spline_direct     = RectBivariateSpline(logM_all, logE, ld, kx=1, ky=3, s=0)
    spline_secondary  = RectBivariateSpline(logM_all, logE, ls, kx=1, ky=3, s=0)
    spline_inflight   = RectBivariateSpline(logM_all, logE, li, kx=1, ky=3, s=0)
    spline_finalstate = RectBivariateSpline(logM_all, logE, lf, kx=1, ky=3, s=0)
    info("Built splines (linear in logM, cubic in logE).")

    all_data = []
    for mval in mass_list:
        snapped = snap_to_available(mval, masses)
        if snapped is not None:
            i = np.where(np.isclose(masses, snapped, rtol=0, atol=0))[0][0]
            kind = 'pre-rendered'
            d = direct_mat[i].copy()
            s = secondary_mat[i].copy()
            it= inflight_mat[i].copy()
            f = finalstate_mat[i].copy()
        else:
            kind = 'interpolated'
            idx_up = int(np.searchsorted(masses, mval, side='left'))
            idx_low = max(0, idx_up-1)
            idx_up  = min(idx_up, N_M-1)
            Ecut = min(Emax_ifa[idx_low], Emax_ifa[idx_up])
            logm = np.log(mval)
            d  = np.exp(spline_direct(logm, logE, grid=False))
            s  = np.exp(spline_secondary(logm, logE, grid=False))
            it = np.exp(spline_inflight(logm, logE, grid=False))
            f  = np.exp(spline_finalstate(logm, logE, grid=False))
            # Guard tails in inflight
            for k in range(len(it)-1, 0, -1):
                if np.isclose(it[k], it[k-1], rtol=1e-8):
                    it[k] = 0.0
                else:
                    break
            log10i = np.log10(np.where(it>0, it, tiny))
            for j in range(1, len(log10i)):
                if log10i[j] - log10i[j-1] < -50:
                    it[j:] = 0.0
                    break
            it[E_ref >= Ecut] = 0.0

        tot = d + s + it + f
        tol = 1e-299
        for arr in (d, s, it, f, tot):
            arr[arr < tol] = 0.0

        # Plot components and total
        plt.figure(figsize=(10,7))
        if np.any(d>0):  plt.plot(E_ref[d>0],  d[d>0],  label="Direct Hawking", lw=2)
        if np.any(s>0):  plt.plot(E_ref[s>0],  s[s>0],  label="Secondary",     lw=2, linestyle='--')
        if np.any(it>0): plt.plot(E_ref[it>0], it[it>0], label="Inflight",      lw=2)
        if np.any(f>0):  plt.plot(E_ref[f>0],  f[f>0],  label="Final State",   lw=2)
        if np.any(tot>0):plt.plot(E_ref[tot>0],tot[tot>0],'k.', label="Total Spectrum")
        plt.xlabel(r'$E_\gamma$ (MeV)')
        plt.ylabel(r'$dN_\gamma/dE_\gamma$ (MeV$^{-1}$ s$^{-1}$)')
        plt.xscale('log'); plt.yscale('log')
        peak_total = tot.max() if tot.size else 1e-20
        plt.ylim(peak_total/1e3, peak_total*1e1)
        plt.xlim(0.5, 5000.0)
        plt.grid(True, which='both', linestyle='--')
        plt.legend()
        plt.title(f'Components for {mval:.2e} g ({kind})')
        plt.tight_layout()
        plt.show()
        plt.close()

        all_data.append({
            'mass': mval, 'kind': kind, 'E': E_ref.copy(),
            'direct': d.copy(), 'secondary': s.copy(),
            'inflight': it.copy(), 'finalstate': f.copy(),
            'total': tot.copy()
        })

    # Overlaid E² dN/dE plot across all requested masses
    if all_data:
        fig = plt.figure(figsize=(10,7))
        summed = np.zeros_like(all_data[0]['E'])
        peaks = []
        for entry in all_data:
            Ecur = entry['E']; tot = entry['total']; valid = tot>0
            if np.any(valid):
                plt.plot(Ecur[valid], Ecur[valid]**2 * tot[valid], lw=2,
                         label=f"{entry['mass']:.2e} g ({entry['kind']})")
                summed += tot
                peaks.append((Ecur[valid]**2 * tot[valid]).max())
        vs = summed > 0
        plt.plot(all_data[0]['E'][vs], all_data[0]['E'][vs]**2 * summed[vs],
                 'k:', lw=3, label="Summed")
        ymax_o = max(peaks) * 1e1
        ymin_o = ymax_o / 1e3
        plt.xlabel(r'$E_\gamma$ (MeV)')
        plt.ylabel(r'$E^2 dN_\gamma/dE_\gamma$ (MeV s$^{-1}$)')
        plt.xscale('log'); plt.yscale('log')
        plt.xlim(0.5, 5000.0); plt.ylim(ymin_o, ymax_o)
        plt.grid(True, which='both', linestyle='--')
        plt.legend()
        plt.title('Total Hawking Radiation Spectra (E²·dN/dE)')
        plt.tight_layout()
        plt.show()
        plt.close(fig)

        sv = user_input("Save any spectra? (y/n): ", allow_back=False, allow_exit=True).strip().lower()
        if sv in ['y', 'yes']:
            print("Select spectra by index to save (single file each):")
            for idx, e in enumerate(all_data, start=1):
                print(f" {idx}: {e['mass']:.2e} g ({e['kind']})")
            choice = user_input("Enter comma-separated indices (e.g. 1,3,5) or '0' to save ALL: ",
                                allow_back=False, allow_exit=True).strip().lower()
            if choice == '0':
                picks = list(range(1, len(all_data)+1))
            else:
                try:
                    picks = [int(x) for x in choice.split(',')]
                except ValueError:
                    err("Invalid indices; skipping save.")
                    picks = []
            for i in picks:
                if 1 <= i <= len(all_data):
                    e = all_data[i - 1]
                    mass_label = f"{e['mass']:.2e}"
                    filename = os.path.join(MONO_RESULTS_DIR, f"{mass_label}_spectrum.txt")
                    data_cols = np.column_stack((
                        e['E'],
                        e['direct'], e['secondary'], e['inflight'], e['finalstate'], e['total']
                    ))
                    header = "E_gamma(MeV)    Direct    Secondary    Inflight    FinalState    Total (MeV^-1 s^-1)"
                    np.savetxt(filename, data_cols, header=header, fmt="%e")
                    print(f"Saved → {filename}")


# ---------------------------
# Right-edge spike trimming helper (kept for reference)
# ---------------------------
def _trim_right_spike(
    x_line: np.ndarray,
    y_line: np.ndarray,
    up_thresh: float = 1.35,
    down_thresh: float = 0.35,
    max_trim_frac: float = 0.10
) -> int:
    """
    Heuristic to trim a suspicious final spike/drop on the right edge of a curve.

    Parameters
    ----------
    x_line : ndarray
        X grid (unused in logic; provided for potential future use).
    y_line : ndarray
        Y values to inspect.
    up_thresh : float
        If y[-1]/y[-2] > up_thresh, treat as spike.
    down_thresh : float
        If y[-1]/y[-2] < down_thresh, treat as plunge.
    max_trim_frac : float
        Do not trim more than this fraction of the array length.

    Returns
    -------
    int
        New usable length index (exclusive). Caller may slice up to this index.
    """
    y = np.asarray(y_line, dtype=float)
    n = y.size
    if n < 3:
        return n - 1
    y_nm1, y_nm2 = y[-1], y[-2]
    if not (np.isfinite(y_nm1) and np.isfinite(y_nm2)) or y_nm2 == 0:
        return n - 1
    ratio = y_nm1 / max(y_nm2, 1e-300)
    if (ratio <= up_thresh) and (ratio >= down_thresh):
        return n - 1
    max_trim = max(3, int(max_trim_frac * n))
    j = n - 1
    trimmed = 0
    if ratio > up_thresh:
        while (j > 1 and trimmed < max_trim and np.isfinite(y[j]) and np.isfinite(y[j-1]) and
               (y[j] / max(y[j-1], 1e-300) > up_thresh)):
            j -= 1; trimmed += 1
        return max(j, 2)
    while (j > 1 and trimmed < max_trim and np.isfinite(y[j]) and np.isfinite(y[j-1]) and
           (y[j] / max(y[j-1], 1e-300) < down_thresh)):
        j -= 1; trimmed += 1
    return max(j, 2)


# ---------------------------
# Distributed (Gaussian collapse / Non-Gaussian / Lognormal)
# ---------------------------
def distributed_spectrum(distribution_method: str) -> None:
    """
    Generate distributed spectra using one of the supported PBH mass distributions.

    Parameters
    ----------
    distribution_method : str
        One of:
            - GAUSSIAN_METHOD ("Gaussian collapse")
            - NON_GAUSSIAN_METHOD ("Non-Gaussian Collapse")
            - LOGNORMAL_METHOD ("Log-Normal Distribution")

    Interactive Flow
    ----------------
    1) Prompt for one or more peak masses (must lie within available pre-rendered grid).
    2) Prompt for target sample size N.
    3) Prompt for distribution-specific width parameter(s) (σ / σ_X / σ in ln-space).
    4) Sample masses, accumulate average spectra via log–log splines (with IFA tail guards).
    5) Plot dN/dE and E² dN/dE overlays across all chosen parameter sets.
    6) For each set, plot its mass histogram with a counts-scaled analytic PDF overlay.
    7) Offer to save results into a unique directory under the method-specific results root.

    Notes
    -----
    - For Non-Gaussian, we enforce 0.04 ≤ σ_X ≤ 0.16 and set σ_Y/σ_X = 0.75 (typical choice).
    - For Log-Normal, we interpret the user's σ as the ln-space standard deviation and choose
      μ such that the mode equals the requested peak (μ_eff = ln(peak) + σ²).
    - All interpolation occurs in (logM, logE) space; inflight annihilation tails are trimmed.
    """
    is_g  = (distribution_method == GAUSSIAN_METHOD)
    is_ng = (distribution_method == NON_GAUSSIAN_METHOD)
    is_ln = (distribution_method == LOGNORMAL_METHOD)

    masses, names = discover_mass_folders(DATA_DIR)
    if not masses:
        warn(f"No valid mass folders found under: {DATA_DIR}")
        return
    MIN_MASS, MAX_MASS = min(masses), max(masses)

    try:
        pstr = user_input(
            f"Enter peak PBH masses (g) (comma-separated; each must be within [{MIN_MASS:.2e}, {MAX_MASS:.2e}]): ",
            allow_back=True, allow_exit=True
        )
    except BackRequested:
        return

    peaks = parse_float_list_verbose(pstr, name="peak mass (g)", bounds=(MIN_MASS, MAX_MASS), allow_empty=False)
    if not peaks:
        warn("No valid peaks; returning.")
        return

    try:
        nstr = user_input("Enter target N (integer, e.g. 1000): ",
                          allow_back=True, allow_exit=True)
    except BackRequested:
        return

    try:
        N_target = int(nstr)
        if N_target <= 0:
            err("N must be > 0. Returning.")
            return
    except Exception:
        err("Invalid N (not an integer). Returning.")
        return

    # collapse parameters (shared constants used in the literature fitting)
    kappa, gamma_p, delta_c = 3.3, 0.36, 0.59

    # read parameter lists
    param_sets = []
    if is_g:
        try:
            sstr = user_input("Enter σ list for Gaussian collapse (comma-separated; each must be within [0.03, 0.255]): ",
                              allow_back=True, allow_exit=True).strip()
        except BackRequested:
            return
        sigmas = parse_float_list_verbose(sstr, name="σ", bounds=(0.03, 0.255), allow_empty=False)
        if not sigmas:
            warn("No valid σ for Gaussian; returning.")
            return
        for sx in sigmas:
            param_sets.append({"sigma_x": sx})

    elif is_ng:
        try:
            sx_str = user_input("Enter σ_X list for Non-Gaussian collapse (comma-separated; σ must be within [0.04, 0.16]): ",
                                allow_back=True, allow_exit=True).strip()
        except BackRequested:
            return
        sigmas_X = parse_float_list_verbose(sx_str, name="σ_X", bounds=(0.04, 0.16), allow_empty=False)
        if not sigmas_X:
            warn("No valid σ for Non-Gaussian; returning.")
            return
        for sX in sigmas_X:
            param_sets.append({"sigma_X": sX, "ratio": 0.75})

    else:  # is_ln
        try:
            sig_str = user_input("Enter σ list (log-space std) for Log-Normal (comma-separated; each > 0): ",
                                 allow_back=True, allow_exit=True).strip()
        except BackRequested:
            return
        sigmas_ln = parse_float_list_verbose(sig_str, name="σ", bounds=(1e-12, None), allow_empty=False, strict_gt=True)
        if not sigmas_ln:
            warn("No valid σ for Log-Normal; returning.")
            return
        for sln in sigmas_ln:
            param_sets.append({"sigma_ln": sln})

    # pre-load all component matrices on a shared energy grid
    first = load_spectra_components(os.path.join(DATA_DIR, names[0]))
    E_grid = first['energy_primary']
    logE = np.log(E_grid)
    N_M = len(masses)

    direct_mat     = np.zeros((N_M, len(E_grid)))
    secondary_mat  = np.zeros_like(direct_mat)
    inflight_mat   = np.zeros_like(direct_mat)
    final_mat      = np.zeros_like(direct_mat)
    Emax_ifa       = np.zeros(N_M)

    for i, m in enumerate(masses):
        sub = os.path.join(DATA_DIR, names[i])
        S = load_spectra_components(sub)
        direct_mat[i]    = S['direct_gamma_primary']
        secondary_mat[i] = np.interp(E_grid, S['energy_secondary'], S['direct_gamma_secondary'], left=0, right=0)
        inflight_mat[i]  = S['IFA_primary'] + np.interp(E_grid, S['energy_secondary'], S['IFA_secondary'], left=0, right=0)
        final_mat[i]     = S['FSR_primary'] + np.interp(E_grid, S['energy_secondary'], S['FSR_secondary'], left=0, right=0)

        p = load_xy_lenient(os.path.join(sub, "inflight_annihilation_prim.txt"))
        s = load_xy_lenient(os.path.join(sub, "inflight_annihilation_sec.txt"))
        Emax_ifa[i] = max(p[:,0].max() if p.size else 0, s[:,0].max() if s.size else 0)

    logM_all = np.log(masses)
    floor = 1e-300

    ld = np.log(np.where(direct_mat    > floor, direct_mat,     floor))
    ls = np.log(np.where(secondary_mat > floor, secondary_mat,  floor))
    li = np.log(np.where(inflight_mat  > floor, inflight_mat,   floor))
    lf = np.log(np.where(final_mat     > floor, final_mat,      floor))

    sp_d = RectBivariateSpline(logM_all, logE, ld, kx=1, ky=3, s=0)
    sp_s = RectBivariateSpline(logM_all, logE, ls, kx=1, ky=3, s=0)
    sp_i = RectBivariateSpline(logM_all, logE, li, kx=1, ky=3, s=0)
    sp_f = RectBivariateSpline(logM_all, logE, lf, kx=1, ky=3, s=0)

    results: list[dict] = []

    for params in param_sets:

        if is_g:
            sigma_x = params["sigma_x"]
            x = np.linspace(0.001, 1.30909, 2000)
            mf = mass_function(delta_l(x, 3.3, 0.59, 0.36), sigma_x, 0.59, 0.36)
            label_param = f"σ={sigma_x:.3g}"
            mf = np.where(np.isfinite(mf) & (mf > 0), mf, 0.0)
            if mf.sum() <= 0:
                warn(f"Underlying PDF vanished for σ={sigma_x:g}; skipping.")
                continue
            probabilities = mf / mf.sum()
            r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]

        elif is_ng:
            sigma_X = params["sigma_X"]; ratio = params["ratio"]; sigma_Y = ratio * sigma_X
            x = np.linspace(0.001, 1.30909, 2000)
            mf = mass_function_exact(delta_l(x, 3.3, 0.59, 0.36), sigma_X, sigma_Y, 0.59, 0.36)
            label_param = f"σX={sigma_X:.3g}"
            mf = np.where(np.isfinite(mf) & (mf > 0), mf, 0.0)
            if mf.sum() <= 0:
                warn(f"Underlying PDF vanished for σ_X={sigma_X:g}; skipping.")
                continue
            probabilities = mf / mf.sum()
            r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]

        else:  # is_ln
            sigma_ln = params["sigma_ln"]
            label_param = f"σ={sigma_ln:.3g}"

        for peak in peaks:
            sum_d = np.zeros_like(E_grid); sum_s = np.zeros_like(E_grid)
            sum_i = np.zeros_like(E_grid); sum_f = np.zeros_like(E_grid)
            md    = []

            bar = tqdm(total=N_target, desc=f"Sampling  peak {peak:.2e}  [{label_param}]", unit="BH")

            if is_ln:
                mu_eff = np.log(peak) + sigma_ln**2
                try:
                    masses_drawn = np.random.lognormal(mean=mu_eff, sigma=sigma_ln, size=N_target)
                except Exception as e:
                    err(f"Sampling error (lognormal, peak {peak:.3e}, σ={sigma_ln:g}): {e}. Skipping.")
                    bar.close()
                    continue
                for mraw in masses_drawn:
                    md.append(float(mraw))
                    if mraw < MIN_MASS or mraw > MAX_MASS:
                        d_vals = s_vals = i_vals = f_vals = np.zeros_like(E_grid)
                    else:
                        try:
                            snap = snap_to_available(mraw, masses)
                            mval = snap if snap else mraw
                            idx_up = int(np.searchsorted(masses, mval, side='left'))
                            idx_low = max(0, idx_up-1)
                            idx_up  = min(idx_up, N_M-1)
                            Ecut = min(Emax_ifa[idx_low], Emax_ifa[idx_up])
                            logm = np.log(mval)
                            d_vals = np.exp(sp_d(logm, logE, grid=False))
                            s_vals = np.exp(sp_s(logm, logE, grid=False))
                            i_vals = np.exp(sp_i(logm, logE, grid=False))
                            f_vals = np.exp(sp_f(logm, logE, grid=False))
                        except Exception as e:
                            warn(f"Interpolation error at mass {mraw:.3e} g: {e}. Skipping draw.")
                            d_vals = s_vals = i_vals = f_vals = np.zeros_like(E_grid)
                        # guard inflight tails
                        for j in range(len(i_vals)-1,0,-1):
                            if np.isclose(i_vals[j], i_vals[j-1], rtol=1e-8): i_vals[j] = 0.0
                            else: break
                        log10i = np.log10(np.where(i_vals>0, i_vals, floor))
                        for j in range(1,len(log10i)):
                            if log10i[j] - log10i[j-1] < -50:
                                i_vals[j:] = 0.0; break
                        i_vals[E_grid >= Ecut] = 0.0
                    sum_d += d_vals; sum_s += s_vals; sum_i += i_vals; sum_f += f_vals
                    bar.update(1)

            else:
                scale = peak / r_mode
                for _ in range(N_target):
                    r = np.random.choice(x, p=probabilities)
                    mraw = r * scale
                    md.append(mraw)
                    if mraw < MIN_MASS or mraw > MAX_MASS:
                        d_vals = s_vals = i_vals = f_vals = np.zeros_like(E_grid)
                    else:
                        try:
                            snap = snap_to_available(mraw, masses)
                            mval = snap if snap else mraw
                            idx_up = int(np.searchsorted(masses, mval, side='left'))
                            idx_low = max(0, idx_up-1)
                            idx_up  = min(idx_up, N_M-1)
                            Ecut = min(Emax_ifa[idx_low], Emax_ifa[idx_up])
                            logm = np.log(mval)
                            d_vals = np.exp(sp_d(logm, logE, grid=False))
                            s_vals = np.exp(sp_s(logm, logE, grid=False))
                            i_vals = np.exp(sp_i(logm, logE, grid=False))
                            f_vals = np.exp(sp_f(logm, logE, grid=False))
                        except Exception as e:
                            warn(f"Interpolation error at mass {mraw:.3e} g: {e}. Skipping draw.")
                            d_vals = s_vals = i_vals = f_vals = np.zeros_like(E_grid)
                        # guard inflight tails
                        for j in range(len(i_vals)-1,0,-1):
                            if np.isclose(i_vals[j], i_vals[j-1], rtol=1e-8): i_vals[j] = 0.0
                            else: break
                        log10i = np.log10(np.where(i_vals>0, i_vals, floor))
                        for j in range(1,len(log10i)):
                            if log10i[j] - log10i[j-1] < -50:
                                i_vals[j:] = 0.0; break
                        i_vals[E_grid >= Ecut] = 0.0
                    sum_d += d_vals; sum_s += s_vals; sum_i += i_vals; sum_f += f_vals
                    bar.update(1)

            bar.close()

            avg_d = sum_d / N_target; avg_s = sum_s / N_target
            avg_i = sum_i / N_target; avg_f = sum_f / N_target
            avg_tot = avg_d + avg_s + avg_i + avg_f
            tol = 1e-299
            for arr in (avg_d, avg_s, avg_i, avg_f, avg_tot):
                arr[arr < tol] = 0.0

            results.append({
                "method": ("gaussian" if is_g else "non_gaussian" if is_ng else "lognormal"),
                "peak": peak,
                "params": params.copy(),
                "E": E_grid.copy(),
                "spectrum": avg_tot.copy(),
                "mdist": md[:],
                "label_param": label_param,
                "nsamp": N_target
            })

    if not results:
        return

    # dN/dE overlays
    fig1 = plt.figure(figsize=(10,7))
    peaks_dn = []
    for r in results:
        E = r["E"]; sp = r["spectrum"]; m = sp > 0
        plt.plot(E[m], sp[m], lw=2,
                 label=f"{distribution_method} {r['peak']:.1e}_{r['label_param'].replace('σ=','').replace('σX=','')}")
        peaks_dn.append(sp.max())
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel(r'$E_\gamma$ (MeV)'); plt.ylabel(r'$dN_\gamma/dE_\gamma$')
    if peaks_dn: plt.ylim(min(peaks_dn)/1e3, max(peaks_dn)*10)
    plt.xlim(0.5, 5e3); plt.grid(True, which='both', linestyle='--'); plt.legend()
    plt.title("Comparison: dN/dE"); plt.tight_layout(); plt.show(); plt.close(fig1)

    # E^2 dN/dE overlays
    fig2 = plt.figure(figsize=(10,7))
    peaks_e2 = []
    for r in results:
        E = r["E"]; sp = r["spectrum"]; m = sp > 0
        plt.plot(E[m], E[m]**2 * sp[m], lw=2,
                 label=f"{distribution_method} {r['peak']:.1e}_{r['label_param'].replace('σ=','').replace('σX=','')}")
        peaks_e2.append((E[m]**2 * sp[m]).max() if np.any(m) else 0.0)
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel(r'$E_\gamma$ (MeV)'); plt.ylabel(r'$E^2\,dN_\gamma/dE_\gamma$')
    if peaks_e2: plt.ylim(min(peaks_e2)/1e3, max(peaks_e2)*10)
    plt.xlim(0.5, 5e3); plt.grid(True, which='both', linestyle='--'); plt.legend()
    plt.title("Comparison: $E^2$ dN/dE"); plt.tight_layout(); plt.show(); plt.close(fig2)

    # Histograms + theoretical mass-PDF overlays (in counts space)
    def _hist_common_bins(md: np.ndarray):
        """Compute histogram bins via Freedman–Diaconis, clamped to [1,50]."""
        md = np.asarray(md, dtype=float)
        md = md[np.isfinite(md)]
        if md.size < 2 or (md.size > 0 and md.min() == md.max()):
            center = md[0] if md.size else 0.0
            eps = abs(center)*1e-9 if center != 0 else 1e-9
            return 1, (center - eps, center + eps), None
        q25, q75 = np.percentile(md, [25, 75])
        iqr = q75 - q25
        if iqr > 0:
            bw = 2 * iqr * md.size ** (-1/3)
            k  = int(np.clip(np.ceil((md.max() - md.min()) / bw), 1, 50))
        else:
            k  = int(np.clip(np.sqrt(md.size), 1, 50))
        return k, None, md

    for r in results:
        method = r["method"]
        figH = plt.figure(figsize=(10,6))

        md = np.asarray(r["mdist"], dtype=float)
        md = md[np.isfinite(md)]

        k, fixed_range, md_safe = _hist_common_bins(md)
        if fixed_range is not None:
            _, bins, _ = plt.hist(md, bins=1, range=fixed_range, alpha=0.7, edgecolor='k',
                                  label=f'{distribution_method} samples ({r["label_param"]})')
        else:
            _, bins, _ = plt.hist(md_safe, bins=k, alpha=0.7, edgecolor='k',
                                  label=f'{distribution_method} samples ({r["label_param"]})')

        bin_widths = (bins[1:] - bins[:-1])
        ref_width  = float(np.median(bin_widths)) if bin_widths.size else 1.0

        if method == "gaussian":
            sigma_x = r["params"]["sigma_x"]
            x = np.linspace(0.001, 1.30909, 2000)
            mf = mass_function(delta_l(x, 3.3, 0.59, 0.36), sigma_x, 0.59, 0.36)
            mf = np.where(np.isfinite(mf) & (mf > 0), mf, 0.0)
            if mf.sum() > 0:
                probabilities = mf / mf.sum()
                r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]
                scale = r["peak"] / r_mode
                dx = x[1] - x[0]; dm = dx * scale
                pdf_mass = probabilities / dm
                m_line = x * scale
                mask = (m_line >= bins[0]) & (m_line <= bins[-1]) & np.isfinite(pdf_mass) & (pdf_mass > 0)
                if np.any(mask):
                    y_line = pdf_mass[mask] * ref_width * len(r["mdist"])
                    plt.plot(m_line[mask], y_line, 'r--', lw=2, zorder=3, label='Underlying PDF (counts)')

        elif method == "non_gaussian":
            sigma_X = r["params"]["sigma_X"]; ratio = 0.75; sigma_Y = ratio * sigma_X
            x = np.linspace(0.001, 1.30909, 2000)
            mf = mass_function_exact(delta_l(x, 3.3, 0.59, 0.36), sigma_X, sigma_Y, 0.59, 0.36)
            mf = np.where(np.isfinite(mf) & (mf > 0), mf, 0.0)
            if mf.sum() > 0:
                probabilities = mf / mf.sum()
                r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]
                scale = r["peak"] / r_mode
                dx = x[1] - x[0]; dm = dx * scale
                pdf_mass = probabilities / dm
                m_line = x * scale
                mask = (m_line >= bins[0]) & (m_line <= bins[-1]) & np.isfinite(pdf_mass) & (pdf_mass > 0)
                if np.any(mask):
                    y_line = pdf_mass[mask] * ref_width * len(r["mdist"])
                    plt.plot(m_line[mask], y_line, 'r--', lw=2, zorder=3, label='Underlying PDF (counts)')

        else:  # lognormal
            sigma_ln = r["params"]["sigma_ln"]; mu_eff = np.log(r["peak"]) + sigma_ln**2
            mlo_tail = np.exp(mu_eff - 6.0*sigma_ln); mhi_tail = np.exp(mu_eff + 6.0*sigma_ln)
            m_plot = np.logspace(np.log10(min(bins[0], mlo_tail)), np.log10(max(bins[-1], mhi_tail)), 2000)
            pdf = (1.0/(m_plot*sigma_ln*np.sqrt(2*np.pi))) * np.exp( - (np.log(m_plot)-mu_eff)**2 / (2*sigma_ln**2) )
            y_plot = pdf * ref_width * len(r["mdist"])
            plt.plot(m_plot, y_plot, 'r--', lw=2, zorder=3, label='Underlying PDF (counts)')
            plt.legend(title=f"σ={sigma_ln:.3f}")

        plt.xlabel('Simulated PBH Mass (g)')
        plt.ylabel('Count')
        plt.title(f'Mass Distribution & PDF for Peak {r["peak"]:.2e} g')
        plt.grid(True, which='both', linestyle='--')
        plt.legend()
        plt.tight_layout()
        plt.show(); plt.close(figH)

    # === Save distributed results ===
    try:
        tosave = user_input("Save distributed results? (y/n): ",
                            allow_back=True, allow_exit=True).strip().lower()
    except BackRequested:
        tosave = 'n'
    if tosave in ('y', 'yes'):
        for r in results:
            method = r["method"]
            if method == "gaussian":
                base = GAUSS_RESULTS_DIR
                tag  = f"peak_{r['peak']:.2e}_{r['label_param'].replace('=','')}_N{r['nsamp']}"
            elif method == "non_gaussian":
                base = NGAUSS_RESULTS_DIR
                tag  = f"peak_{r['peak']:.2e}_{r['label_param'].replace('=','')}_N{r['nsamp']}"
            else:
                base = LOGN_RESULTS_DIR
                tag  = f"peak_{r['peak']:.2e}_{r['label_param'].replace('=','')}_N{r['nsamp']}"
            outdir = os.path.join(base, tag)
            k = 1; unique = outdir
            while os.path.exists(unique):
                unique = f"{outdir}_{k}"; k += 1
            os.makedirs(unique, exist_ok=True)
            np.savetxt(os.path.join(unique, "distributed_spectrum.txt"),
                       np.column_stack((r["E"], r["spectrum"])),
                       header="E_gamma(MeV)   TotalSpectrum", fmt="%.10e")
            np.savetxt(os.path.join(unique, "mass_distribution.txt"),
                       np.asarray(r["mdist"], dtype=float),
                       header="Sampled masses (g)", fmt="%.12e")
            print(f"Saved → {unique}")


# ---------------------------
# Helpers for Custom Equation: safe eval + variable prompting
# ---------------------------
def _build_safe_numpy_namespace() -> SimpleNamespace:
    """
    Build a restricted numpy-like namespace exposing only safe math functions.

    Returns
    -------
    SimpleNamespace
        Object exposing e.g. log, exp, sqrt, sin/cos/tan, etc.
    """
    safe_np = SimpleNamespace(
        log=np.log, log10=np.log10, log1p=np.log1p, exp=np.exp, sqrt=np.sqrt, power=np.power,
        sin=np.sin, cos=np.cos, tan=np.tan, arctan=np.arctan,
        abs=np.abs, minimum=np.minimum, maximum=np.maximum, clip=np.clip, erf=erf,
        pi=np.pi, e=np.e
    )
    return safe_np


SAFE_FUNCS = {
    "log","log10","log1p","exp","sqrt","pow","sin","cos","tan","arctan",
    "abs","minimum","maximum","clip","erf","pi","e","m","np","numpy"
}


def _detect_custom_variables(expr: str) -> list[str]:
    """
    Detect identifiers in a user expression that are not known safe names.

    Parameters
    ----------
    expr : str
        RHS expression in variable `m` (grams).

    Returns
    -------
    list[str]
        Sorted names of variables that require values from the user.

    Notes
    -----
    - Greek letters like 'μ','α','β' are supported as identifiers.
    - Strings are stripped to avoid false positives.
    """
    expr_wo_strings = re.sub(r"(\".*?\"|'.*?')", "", expr)
    tokens = set(re.findall(r"\b[^\W\d]\w*\b", expr_wo_strings, flags=re.UNICODE))
    unknown = sorted([t for t in tokens if t not in SAFE_FUNCS])
    return unknown


def _prompt_variable_values(var_names: list[str]) -> dict[str, float]:
    """
    Prompt the user for each variable value. Accepts numeric expressions
    using pi, e, and np.*.

    Parameters
    ----------
    var_names : list[str]
        Variables needing values.

    Returns
    -------
    dict[str, float]
        Mapping from name → float value.

    Raises
    ------
    BackRequested
        If the user backs out.
    SystemExit
        If the user exits.
    """
    vals: dict[str, float] = {}
    safe_np = _build_safe_numpy_namespace()
    num_ctx = {"__builtins__": None, "pi": np.pi, "e": np.e, "np": safe_np, "numpy": safe_np}
    for name in var_names:
        while True:
            try:
                s = user_input(f"Enter value for variable '{name}': ",
                               allow_back=True, allow_exit=True).strip()
                val = eval(s, {"__builtins__": None}, num_ctx)
                val = float(val)
                vals[name] = val
                break
            except BackRequested:
                raise
            except SystemExit:
                raise
            except Exception:
                err("Could not parse value. Use a number or an expression like '1e16' or '2*np.pi'. Try again.")
    return vals


# ---------------------------
# Custom Mass PDF from user-entered EQUATION
# ---------------------------
def custom_equation_pdf_tool() -> None:
    """
    Build a PBH mass PDF from a user-entered equation f(m, params...), normalize it per gram,
    sample N masses, accumulate ONLY the TOTAL spectrum, then show:

      (1) total dN/dE,
      (2) total E^2 dN/dE,
      (3) mass histogram (counts) with analytic PDF scaled to counts (log bins).

    Saved outputs (if requested)
    ----------------------------
    - equation.txt           : The expression and any variable values (commented).
    - samples_sorted.txt     : Sorted sampled masses (g).
    - distributed_spectrum.txt : Columns: E_gamma(MeV), TotalSpectrum

    Notes
    -----
    - The expression must be a *right-hand-side* function of `m` (no "f(m)=" prefix needed).
    - Allowed functions: subset from numpy (log/exp/sqrt/sin/cos/tan/arctan/abs/clip/min/max/erf).
    - Variables unknown to the safe namespace will be auto-detected and prompted for.
    """
    # Discover data domain
    masses, names = discover_mass_folders(DATA_DIR)
    if masses:
        M_MIN, M_MAX = min(masses), max(masses)
    else:
        M_MIN, M_MAX = 5e13, 1e19

    N_BINS = 50

    def log_edges(a, b, k):
        return np.logspace(np.log10(a), np.log10(b), k + 1)

    def safe_eval_on_grid(expr, m_grid, user_vars):
        safe_np = _build_safe_numpy_namespace()
        safe = {
            "m": m_grid,
            "log": np.log, "log10": np.log10, "log1p": np.log1p,
            "exp": np.exp, "sqrt": np.sqrt, "pow": np.power,
            "sin": np.sin, "cos": np.cos, "tan": np.tan, "arctan": np.arctan,
            "abs": np.abs, "minimum": np.minimum, "maximum": np.maximum, "clip": np.clip,
            "erf": erf, "pi": np.pi, "e": np.e,
            "np": safe_np, "numpy": safe_np
        }
        safe.update(user_vars)
        try:
            y = eval(expr, {"__builtins__": None}, safe)
        except BackRequested:
            raise
        except SystemExit:
            raise
        except Exception as e:
            raise ValueError(f"Could not evaluate expression: {e}")
        y = np.asarray(y, dtype=float)
        if y.size == 1:
            y = np.full_like(m_grid, float(y))
        if y.shape != m_grid.shape:
            raise ValueError("Expression did not return an array of the same shape as m.")
        return y

    def cdf_from_pdf(m, pdf):
        cdf = np.empty_like(pdf)
        cdf[0] = 0.0
        dm = np.diff(m)
        cdf[1:] = np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * dm)
        total = cdf[-1]
        if not np.isfinite(total) or total <= 0:
            raise ValueError("PDF integrates to non-positive value.")
        cdf /= total
        return cdf

    # ---- read the equation & prompt variables ----
    print("\n=== Custom Equation Mass PDF ===")
    print("Domain: m in [{:.2e}, {:.2e}] g".format(M_MIN, M_MAX))
    print("Enter a Python expression for your PDF f(m) using 'm' in grams and any constants/variables you define.")
    print("Examples:")
    print("f(m) = (m/mp)**(-(alpha0 + beta*log(m/mp))) / m")
    print("f(m) = exp(-m/5e17) / m")
    try:
        expr = user_input("f(m) = ", allow_back=True, allow_exit=True).strip()
    except BackRequested:
        return

    # If someone pastes "f(m) = ..." or "fm = ...", strip the prefix anyway.
    expr = re.sub(r'^\s*(?:f\s*\(\s*m\s*\)|fm)\s*=\s*', '', expr, flags=re.IGNORECASE)

    # Detect custom variables (excluding allowed function names, m, pi, e, np, numpy)
    vars_needed = _detect_custom_variables(expr)
    user_vars = {}
    if vars_needed:
        info(f"Variables detected: {', '.join(vars_needed)}")
        try:
            user_vars = _prompt_variable_values(vars_needed)
        except BackRequested:
            return

    # ---- build normalized PDF on a fine m-grid ----
    m_grid = np.logspace(np.log10(M_MIN), np.log10(M_MAX), 20000)
    try:
        f = safe_eval_on_grid(expr, m_grid, user_vars)
    except BackRequested:
        return
    except ValueError as e:
        err(str(e))
        return

    f = np.clip(f, 0.0, None)
    area = trapezoid(f, m_grid)
    if not np.isfinite(area) or area <= 0.0:
        err("Your f(m) is nonpositive or non-integrable over the domain.")
        return
    pdf = f / area  # per gram
    cdf = cdf_from_pdf(m_grid, pdf)

    # ---- ask for N ----
    try:
        n_default = 1000
        n_str = user_input(f"Enter target N (integer, e.g. 1000):  ",
                           allow_back=True, allow_exit=True).strip()
        if n_str == "":
            N = n_default
        else:
            N = int(n_str)
            if N <= 0:
                err("N must be > 0.")
                return
    except BackRequested:
        return
    except Exception:
        err("Invalid N (must be a positive integer).")
        return

    # ---- pre-load spectral grids & splines ----
    if not masses:
        err("No valid mass folders found under the data directory.")
        return

    first = load_spectra_components(os.path.join(DATA_DIR, names[0]))
    E_grid = first['energy_primary']
    logE   = np.log(E_grid)
    N_M    = len(masses)

    direct_mat     = np.zeros((N_M, len(E_grid)))
    secondary_mat  = np.zeros_like(direct_mat)
    inflight_mat   = np.zeros_like(direct_mat)
    final_mat      = np.zeros_like(direct_mat)
    Emax_ifa       = np.zeros(N_M)

    for i, m in enumerate(masses):
        sub = os.path.join(DATA_DIR, names[i])
        S = load_spectra_components(sub)
        direct_mat[i]    = S['direct_gamma_primary']
        secondary_mat[i] = np.interp(E_grid, S['energy_secondary'], S['direct_gamma_secondary'], left=0, right=0)
        inflight_mat[i]  = S['IFA_primary'] + np.interp(E_grid, S['energy_secondary'], S['IFA_secondary'], left=0, right=0)
        final_mat[i]     = S['FSR_primary'] + np.interp(E_grid, S['energy_secondary'], S['FSR_secondary'], left=0, right=0)

        p = load_xy_lenient(os.path.join(sub, "inflight_annihilation_prim.txt"))
        s = load_xy_lenient(os.path.join(sub, "inflight_annihilation_sec.txt"))
        Emax_ifa[i] = max(p[:,0].max() if p.size else 0, s[:,0].max() if s.size else 0)

    logM_all = np.log(masses)
    floor = 1e-300
    ld = np.log(np.where(direct_mat    > floor, direct_mat,     floor))
    ls = np.log(np.where(secondary_mat > floor, secondary_mat,  floor))
    li = np.log(np.where(inflight_mat  > floor, inflight_mat,   floor))
    lf = np.log(np.where(final_mat     > floor, final_mat,      floor))

    sp_d = RectBivariateSpline(logM_all, logE, ld, kx=1, ky=3, s=0)
    sp_s = RectBivariateSpline(logM_all, logE, ls, kx=1, ky=3, s=0)
    sp_i = RectBivariateSpline(logM_all, logE, li, kx=1, ky=3, s=0)
    sp_f = RectBivariateSpline(logM_all, logE, lf, kx=1, ky=3, s=0)

    # ---- sample masses via inverse CDF and accumulate ONLY total spectrum ----
    rng = np.random.default_rng()
    u = rng.random(N)
    samples = np.interp(u, cdf, m_grid)
    samples.sort()

    sum_tot = np.zeros_like(E_grid)

    bar = tqdm(total=N, desc="Sampling custom PDF", unit="BH")
    for mraw in samples:
        if mraw < masses[0] or mraw > masses[-1]:
            bar.update(1)
            continue
        try:
            snap = snap_to_available(mraw, masses)
            mval = snap if snap else mraw
            idx_up  = int(np.searchsorted(masses, mval, side='left'))
            idx_low = max(0, idx_up-1)
            idx_up  = min(idx_up, len(masses)-1)
            Ecut    = min(Emax_ifa[idx_low], Emax_ifa[idx_up])
            logm    = np.log(mval)
            d_vals  = np.exp(sp_d(logm, logE, grid=False))
            s_vals  = np.exp(sp_s(logm, logE, grid=False))
            i_vals  = np.exp(sp_i(logm, logE, grid=False))
            f_vals  = np.exp(sp_f(logm, logE, grid=False))
        except Exception:
            bar.update(1)
            continue

        # trim inflight tails (stability)
        for j in range(len(i_vals)-1, 0, -1):
            if np.isclose(i_vals[j], i_vals[j-1], rtol=1e-8):
                i_vals[j] = 0.0
            else:
                break
        log10i = np.log10(np.where(i_vals > 0, i_vals, floor))
        for j in range(1, len(log10i)):
            if log10i[j] - log10i[j-1] < -50:
                i_vals[j:] = 0.0
                break
        i_vals[E_grid >= Ecut] = 0.0

        sum_tot += (d_vals + s_vals + i_vals + f_vals)
        bar.update(1)
    bar.close()

    avg_tot = sum_tot / max(N, 1)
    avg_tot[avg_tot < 1e-299] = 0.0

    # ---- FIGURE A: total dN/dE ----
    msk = avg_tot > 0
    plt.figure(figsize=(10, 7))
    plt.plot(E_grid[msk], avg_tot[msk], lw=2, label="Total spectrum")
    plt.xscale('log'); plt.yscale('log')
    plt.xlim(0.5, 5e3)
    if np.any(msk):
        peak = avg_tot[msk].max()
        plt.ylim(peak/1e3, peak*10)
    plt.xlabel(r'$E_\gamma$ (MeV)')
    plt.ylabel(r'$dN_\gamma/dE_\gamma$ (MeV$^{-1}$ s$^{-1}$)')
    plt.grid(True, which='both', linestyle='--')
    plt.legend()
    plt.title("Custom Equation — Total $dN/dE$")
    plt.tight_layout()
    plt.show()

    # ---- FIGURE B: total E^2 dN/dE ----
    plt.figure(figsize=(10, 7))
    if np.any(msk):
        plt.plot(E_grid[msk], (E_grid[msk]**2) * avg_tot[msk], lw=2, label="Total")
        peak_e2 = ((E_grid[msk]**2) * avg_tot[msk]).max()
        plt.ylim(peak_e2/1e3, peak_e2*10)
    plt.xscale('log'); plt.yscale('log')
    plt.xlim(0.5, 5e3)
    plt.xlabel(r'$E_\gamma$ (MeV)')
    plt.ylabel(r'$E^2\,dN_\gamma/dE_\gamma$ (MeV s$^{-1}$)')
    plt.grid(True, which='both', linestyle='--')
    plt.legend()
    plt.title("Custom Equation — Total $E^2 dN/dE$")
    plt.tight_layout()
    plt.show()

    # ---- FIGURE C: Mass histogram (counts) + SMOOTH analytic PDF scaled to counts ----
    edges = log_edges(masses[0], masses[-1], N_BINS)
    plt.figure(figsize=(10, 6))
    # Blue = sampled counts per (log) bin
    plt.hist(samples, bins=edges, density=False, alpha=0.6, edgecolor='k',
             label=f"Sampled counts per bin (N={N})")

    # Orange = smooth line proportional to expected counts/bin for log bins:
    # expected counts in a narrow log bin: N * pdf(m) * m * d(ln m)
    dln = (np.log(masses[-1]) - np.log(masses[0])) / N_BINS
    counts_line = N * pdf * m_grid * dln
    plt.plot(m_grid, counts_line, lw=2.5, label="Analytic PDF (scaled to counts)")
    plt.xscale("log")
    plt.xlabel("Mass m (g)")
    plt.ylabel("Count per bin")
    plt.title("Custom Equation — Mass Histogram (counts) + Smooth PDF overlay")
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ---- Save exactly 3 files for custom: equation, mass distribution, distributed spectrum ----
    try:
        sv = user_input("\nSave this custom spectrum? (y/n): ",
                        allow_back=True, allow_exit=True).strip().lower()
    except BackRequested:
        return
    if sv in ('y', 'yes'):
        median_mass = float(np.median(samples)) if samples.size else 0.0
        folder = f"{median_mass:.2e}_custom_eq"
        outdir = os.path.join(CUSTOM_RESULTS_DIR, folder)
        base = outdir; k = 1
        while os.path.exists(outdir):
            outdir = f"{base}_{k}"; k += 1
        os.makedirs(outdir, exist_ok=True)

        with open(os.path.join(outdir, "equation.txt"), "w", encoding="utf-8") as fh:
            if user_vars:
                fh.write("# Variables:\n")
                for kname, kval in user_vars.items():
                    fh.write(f"# {kname} = {kval:.10e}\n")
            fh.write(expr + "\n")
        np.savetxt(os.path.join(outdir, "samples_sorted.txt"), samples,
                   header="Simulated masses (g), sorted ascending", fmt="%.12e")
        np.savetxt(os.path.join(outdir, "distributed_spectrum.txt"),
                   np.column_stack((E_grid, avg_tot)),
                   header="E_gamma(MeV)   TotalSpectrum", fmt="%.10e")
        print(f"Saved → {outdir}")


# ---------------------------
# View previous spectra (with queue)
# ---------------------------
def view_previous_spectra() -> None:
    """
    View previously saved spectra with a queue:
      - Selecting items adds them to the queue only.
      - Press '0' to plot ALL queued items: spectra first (dN/dE, E^2 dN/dE), then histograms.
      - Queue auto-clears after plotting.
    """
    # --- allowed mono input range (use discovered data domain) ---
    masses_all, names_all = discover_mass_folders(DATA_DIR)
    if masses_all:
        M_MIN_MONO, M_MAX_MONO = min(masses_all), max(masses_all)
    else:
        M_MIN_MONO, M_MAX_MONO = 5e13, 1e19

    cat_map = {
        '1': ("Monochromatic Distribution", MONO_RESULTS_DIR,  None,                          "mono"),
        '2': (GAUSSIAN_METHOD,             GAUSS_RESULTS_DIR,  "distributed_spectrum.txt",    "gaussian"),
        '3': (NON_GAUSSIAN_METHOD,         NGAUSS_RESULTS_DIR, "distributed_spectrum.txt",    "non_gaussian"),
        '4': (LOGNORMAL_METHOD,            LOGN_RESULTS_DIR,   "distributed_spectrum.txt",    "lognormal"),
        '5': ("Custom equation (user-defined mass PDF)", CUSTOM_RESULTS_DIR, "distributed_spectrum.txt", "custom"),
    }

    # ---------- helpers: text cleanup & equation parsing ----------
    def _strip_invisibles(s: str) -> str:
        for ch in (
            "\ufeff","\u200b","\u200c","\u200d","\u2060","\u200e","\u200f",
            "\u202a","\u202b","\u202c","\u202d","\u202e","\u202f",
            "\u00a0","\r"
        ):
            s = s.replace(ch, "")
        return s

    def _normalize_expr_line(s: str) -> str:
        s = _strip_invisibles(s.strip())
        s = re.sub(r'^\s*(?:f\s*\(\s*m\s*\)|fm)\s*=\s*','',s,flags=re.IGNORECASE)
        s = (s.replace('^','**')
               .replace('×','*')
               .replace('·','*')
               .replace('÷','/')
               .replace('−','-')
               .replace('—','-')
               .replace('–','-')
               .replace('“','"').replace('”','"')
               .replace('’',"'").replace('‘',"'"))
        out, in_sin, in_dbl = [], False, False
        for ch in s:
            if ch == "'" and not in_dbl:
                in_sin = not in_sin
            elif ch == '"' and not in_sin:
                in_dbl = not in_dbl
            if ch == '#' and not in_sin and not in_dbl:
                break
            out.append(ch)
        return ''.join(out).strip()

    def _read_equation_file(run_dir: str) -> tuple[str, dict[str, float]]:
        eq_path = os.path.join(run_dir, "equation.txt")
        try:
            try:
                lines = open(eq_path,"r",encoding="utf-8-sig").readlines()
            except UnicodeDecodeError:
                lines = open(eq_path,"r",encoding="latin-1").readlines()
        except Exception as e:
            raise RuntimeError(f"Cannot read equation.txt: {e}")

        user_vars: dict[str, float] = {}
        expr = None
        for raw in lines:
            s = _strip_invisibles(raw).strip()
            if not s:
                continue
            if s.startswith("#"):
                if "=" in s[1:]:
                    try:
                        k,v = s[1:].split("=",1)
                        user_vars[k.strip()] = float(_strip_invisibles(v).strip())
                    except Exception:
                        pass
                continue
            norm = _normalize_expr_line(s)
            if norm:
                expr = norm
        if not expr:
            for raw in reversed(lines):
                s = raw.strip()
                if s and not s.lstrip().startswith("#"):
                    s = _normalize_expr_line(s)
                    if s:
                        expr = s
                        break
        if not expr:
            raise RuntimeError("No custom equation found in equation.txt.")
        return expr, user_vars

    def _safe_eval_on_grid(expr: str, m_grid: np.ndarray, user_vars: dict[str, float]) -> np.ndarray:
        safe_np = _build_safe_numpy_namespace()
        safe = {
            "m": m_grid,
            "log": np.log, "log10": np.log10, "log1p": np.log1p,
            "exp": np.exp, "sqrt": np.sqrt, "pow": np.power,
            "sin": np.sin, "cos": np.cos, "tan": np.tan, "arctan": np.arctan,
            "abs": np.abs, "minimum": np.minimum, "maximum": np.maximum, "clip": np.clip,
            "erf": erf, "pi": np.pi, "e": np.e,
            "np": _build_safe_numpy_namespace(), "numpy": _build_safe_numpy_namespace()
        }
        safe.update(user_vars)
        y = eval(expr, {"__builtins__": None}, safe)
        y = np.asarray(y, dtype=float)
        if y.size == 1:
            y = np.full_like(m_grid, float(y))
        if y.shape != m_grid.shape:
            raise ValueError("Expression did not return an array of the same shape as m.")
        return y

    # ---------- parse run_name for peak and sigma ----------
    def _extract_peak_sigma(run_name: str, kind: str) -> tuple[float | None, float | None, str | None]:
        peak_val = None
        sigma_val = None
        sigma_str = None
        m_peak = re.search(r"peak_([0-9.+\-eE]+)", run_name)
        if m_peak:
            try:
                peak_val = float(m_peak.group(1))
            except Exception:
                peak_val = None
        if kind == "non_gaussian":
            m_sigx = re.search(r"σX([0-9.]+)", run_name)
            if m_sigx:
                try:
                    sigma_val = float(m_sigx.group(1))
                except Exception:
                    sigma_val = None
                sigma_str = f"σX={sigma_val:.3g}" if sigma_val is not None else "σX=?"
        else:
            m_sig = re.search(r"σ([0-9.]+)", run_name)
            if m_sig:
                try:
                    sigma_val = float(m_sig.group(1))
                except Exception:
                    sigma_val = None
                sigma_str = f"σ={sigma_val:.3g}" if sigma_val is not None else "σ=?"
        return peak_val, sigma_val, sigma_str

    # ---------- plotting helpers (updated sizes) ----------
    def _plot_dn(results: list[tuple[str, tuple[np.ndarray, np.ndarray]]]) -> None:
        if not results:
            return
        fig = plt.figure(figsize=(10,7))
        peaks = []
        for lab, (E, S) in results:
            msk = S > 0
            plt.plot(E[msk], S[msk], lw=2, label=lab)
            if np.any(msk):
                peaks.append(S[msk].max())
        plt.xscale('log'); plt.yscale('log')
        plt.xlabel(r'$E_\gamma$ (MeV)')
        plt.ylabel(r'$dN_\gamma/dE_\gamma$')
        if peaks:
            plt.ylim(min(peaks)/1e3, max(peaks)*10)
        plt.xlim(0.5, 5e3)
        plt.grid(True, which='both', linestyle='--')
        plt.legend()
        plt.title("Comparison: dN/dE")
        plt.tight_layout()
        plt.show()
        plt.close(fig)

    def _plot_e2(results: list[tuple[str, tuple[np.ndarray, np.ndarray]]]) -> None:
        if not results:
            return
        fig = plt.figure(figsize=(10,7))
        peaks = []
        for lab, (E, S) in results:
            msk = S > 0
            plt.plot(E[msk], (E[msk]**2)*S[msk], lw=2, label=lab)
            if np.any(msk):
                peaks.append(((E[msk]**2)*S[msk]).max())
        plt.xscale('log'); plt.yscale('log')
        plt.xlabel(r'$E_\gamma$ (MeV)')
        plt.ylabel(r'$E^2\,dN_\gamma/dE_\gamma$')
        if peaks:
            plt.ylim(min(peaks)/1e3, max(peaks)*10)
        plt.xlim(0.5, 5e3)
        plt.grid(True, which='both', linestyle='--')
        plt.legend()
        plt.title("Comparison: $E^2$ dN/dE")
        plt.tight_layout()
        plt.show()
        plt.close(fig)

    def _hist_gaussian(samples, peak, sigma_x, title_prefix):
        """Histogram + counts-scaled Gaussian-collapse PDF overlay (bigger figure)."""
        md = np.asarray(samples, dtype=float)
        md = md[np.isfinite(md)]
        if md.size == 0:
            return
        plt.figure(figsize=(10,6))  # <-- bigger now
        if md.size < 2 or md.min() == md.max():
            center = md[0]
            eps = abs(center)*1e-9 if center != 0 else 1e-9
            _, bins, _ = plt.hist(md, bins=1, range=(center-eps, center+eps),
                                  alpha=0.7, edgecolor='k', label='samples')
        else:
            q25, q75 = np.percentile(md, [25, 75])
            iqr = q75 - q25
            if iqr > 0:
                bw = 2 * iqr * md.size ** (-1/3)
                k  = int(np.clip(np.ceil((md.max() - md.min()) / bw), 1, 50))
            else:
                k  = int(np.clip(np.sqrt(md.size), 1, 50))
            _, bins, _ = plt.hist(md, bins=k, alpha=0.7, edgecolor='k', label='samples')

        x  = np.linspace(0.001, 1.30909, 2000)
        mf = mass_function(delta_l(x,3.3,0.59,0.36), sigma_x, 0.59, 0.36)
        mf = np.where(np.isfinite(mf) & (mf>0), mf, 0.0)
        if mf.sum() > 0:
            probabilities = mf / mf.sum()
            r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]
            scale = peak / r_mode if peak is not None else 1.0
            dx = x[1] - x[0]
            dm = dx * scale
            pdf_mass = probabilities / dm
            m_line = x * scale
            bin_widths = (bins[1:] - bins[:-1])
            ref_width  = float(np.median(bin_widths)) if bin_widths.size else 1.0
            mask = ((m_line >= bins[0]) & (m_line <= bins[-1]) &
                    np.isfinite(pdf_mass) & (pdf_mass > 0))
            if np.any(mask):
                y_line = pdf_mass[mask] * ref_width * len(md)
                plt.plot(m_line[mask], y_line, 'r--', lw=2, zorder=3,
                         label='Underlying PDF (counts)')

        plt.xlabel('Simulated PBH Mass (g)')
        plt.ylabel('Count')
        plt.title(f'{title_prefix} — Mass Distribution & PDF overlay')
        plt.grid(True, which='both', linestyle='--')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def _hist_nongaussian(samples, peak, sigma_X, title_prefix):
        """Histogram + counts-scaled Non-Gaussian PDF overlay (bigger figure)."""
        md = np.asarray(samples, dtype=float)
        md = md[np.isfinite(md)]
        if md.size == 0:
            return
        plt.figure(figsize=(10,6))  # <-- bigger now
        if md.size < 2 or md.min() == md.max():
            center = md[0]
            eps = abs(center)*1e-9 if center != 0 else 1e-9
            _, bins, _ = plt.hist(md, bins=1, range=(center-eps, center+eps),
                                  alpha=0.7, edgecolor='k', label='samples')
        else:
            q25, q75 = np.percentile(md, [25, 75])
            iqr = q75 - q25
            if iqr > 0:
                bw = 2 * iqr * md.size ** (-1/3)
                k  = int(np.clip(np.ceil((md.max() - md.min()) / bw), 1, 50))
            else:
                k  = int(np.clip(np.sqrt(md.size), 1, 50))
            _, bins, _ = plt.hist(md, bins=k, alpha=0.7, edgecolor='k', label='samples')

        x = np.linspace(0.001, 1.30909, 2000)
        sigma_Y = 0.75 * (sigma_X if sigma_X is not None else 0.0)
        mf = mass_function_exact(delta_l(x,3.3,0.59,0.36),
                                 sigma_X if sigma_X is not None else 0.0,
                                 sigma_Y,
                                 0.59, 0.36)
        mf = np.where(np.isfinite(mf) & (mf>0), mf, 0.0)
        if mf.sum() > 0:
            probabilities = mf / mf.sum()
            r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]
            scale = peak / r_mode if peak is not None else 1.0
            dx = x[1] - x[0]
            dm = dx * scale
            pdf_mass = probabilities / dm
            m_line = x * scale
            bin_widths = (bins[1:] - bins[:-1])
            ref_width  = float(np.median(bin_widths)) if bin_widths.size else 1.0
            mask = ((m_line >= bins[0]) & (m_line <= bins[-1]) &
                    np.isfinite(pdf_mass) & (pdf_mass > 0))
            if np.any(mask):
                y_line = pdf_mass[mask] * ref_width * len(md)
                plt.plot(m_line[mask], y_line, 'r--', lw=2, zorder=3,
                         label='Underlying PDF (counts)')

        plt.xlabel('Simulated PBH Mass (g)')
        plt.ylabel('Count')
        plt.title(f'{title_prefix} — Mass Distribution & PDF overlay')
        plt.grid(True, which='both', linestyle='--')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def _hist_lognormal(samples, peak, sigma_ln, title_prefix):
        """Histogram + counts-scaled Log-normal PDF overlay (bigger figure)."""
        md = np.asarray(samples, dtype=float)
        md = md[np.isfinite(md)]
        if md.size == 0:
            return
        plt.figure(figsize=(10,6))  # <-- bigger now
        if md.size < 2 or md.min() == md.max():
            center = md[0]
            eps = abs(center)*1e-9 if center != 0 else 1e-9
            _, bins, _ = plt.hist(md, bins=1, range=(center-eps, center+eps),
                                  alpha=0.7, edgecolor='k', label='samples')
        else:
            q25, q75 = np.percentile(md, [25, 75])
            iqr = q75 - q25
            if iqr > 0:
                bw = 2 * iqr * md.size ** (-1/3)
                k  = int(np.clip(np.ceil((md.max() - md.min()) / bw), 1, 50))
            else:
                k  = int(np.clip(np.sqrt(md.size), 1, 50))
            _, bins, _ = plt.hist(md, bins=k, alpha=0.7, edgecolor='k', label='samples')

        bin_widths = (bins[1:] - bins[:-1])
        ref_width  = float(np.median(bin_widths)) if bin_widths.size else 1.0

        mu_eff = None
        if peak is not None and sigma_ln is not None:
            mu_eff = np.log(peak) + sigma_ln**2

        if mu_eff is not None:
            mlo_tail = np.exp(mu_eff - 6.0*sigma_ln)
            mhi_tail = np.exp(mu_eff + 6.0*sigma_ln)
            m_plot = np.logspace(np.log10(min(bins[0], mlo_tail)),
                                 np.log10(max(bins[-1], mhi_tail)),
                                 2000)
            pdf = (1.0/(m_plot*sigma_ln*np.sqrt(2*np.pi))) * np.exp(
                - (np.log(m_plot)-mu_eff)**2 / (2*sigma_ln**2)
            )
            y_plot = pdf * ref_width * len(md)
            plt.plot(m_plot, y_plot, 'r--', lw=2, zorder=3,
                     label='Underlying PDF (counts)')

        plt.xlabel('Simulated PBH Mass (g)')
        plt.ylabel('Count')
        plt.title(f'{title_prefix} — Mass Distribution & PDF overlay')
        plt.grid(True, which='both', linestyle='--')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def _hist_custom(run_dir, title_prefix):
        """Histogram for a custom-equation run (size already large at 10×6 in generator)."""
        spath = os.path.join(run_dir, "samples_sorted.txt")
        try:
            samples = np.loadtxt(spath)
        except Exception:
            return
        samples = np.asarray(samples, dtype=float)
        samples = samples[np.isfinite(samples)]
        if samples.size == 0:
            return

        try:
            expr, user_vars = _read_equation_file(run_dir)
        except Exception:
            expr = None
            user_vars = {}

        masses_all2, _ = discover_mass_folders(DATA_DIR)
        if masses_all2:
            M_MIN, M_MAX = min(masses_all2), max(masses_all2)
        else:
            M_MIN, M_MAX = 5e13, 1e19

        pdf = None
        if expr is not None:
            m_grid = np.logspace(np.log10(M_MIN), np.log10(M_MAX), 20000)
            try:
                f = _safe_eval_on_grid(expr, m_grid, user_vars)
            except Exception:
                f = None
            if f is not None:
                f = np.clip(f, 0.0, None)
                area = trapezoid(f, m_grid)
                if np.isfinite(area) and area > 0:
                    pdf = f / area

        N_BINS = 50
        edges = np.logspace(np.log10(M_MIN), np.log10(M_MAX), N_BINS + 1)
        plt.figure(figsize=(10,6))
        plt.hist(samples, bins=edges, density=False, alpha=0.6, edgecolor='k',
                 label=f"Sampled counts per bin (N={len(samples)})")

        if pdf is not None:
            dln = (np.log(M_MAX) - np.log(M_MIN)) / N_BINS
            counts_line = len(samples) * pdf * m_grid * dln
            plt.plot(m_grid, counts_line, lw=2.5, label="Analytic PDF (scaled to counts)")

        plt.xscale("log")
        plt.xlabel("Mass m (g)")
        plt.ylabel("Count per bin")
        plt.title(f"{title_prefix} — Mass Histogram (counts) + Smooth PDF")
        plt.grid(True, which='both', linestyle='--', alpha=0.5)
        plt.legend()
        plt.tight_layout()
        plt.show()

    # ---------- UI ----------
    def _print_menu():
        print("\nView Previous — choose:")
        print(" 1: Monochromatic Distribution")
        print(f" 2: {GAUSSIAN_METHOD}")
        print(f" 3: {NON_GAUSSIAN_METHOD}")
        print(f" 4: {LOGNORMAL_METHOD}")
        print(f" 5: Custom equation (user-defined mass PDF)")
        print(" 0: Plot all Queued | b: Back | q: Quit")

    # queue elements:
    queue: list[dict] = []

    # Preload matrices/splines once for fast monochromatic interpolation in this view
    mono_ready = False
    mono_E = None
    sp_d = sp_s = sp_i = sp_f = None
    logE = None
    Emax_ifa = None
    floor = 1e-300
    N_M = 0

    def _ensure_mono_interpolator():
        nonlocal mono_ready, mono_E, sp_d, sp_s, sp_i, sp_f, logE, Emax_ifa, N_M
        if mono_ready:
            return
        if not masses_all:
            raise RuntimeError("No valid mass folders found under data directory.")
        first = load_spectra_components(os.path.join(DATA_DIR, names_all[0]))
        mono_E = first['energy_primary']
        logE = np.log(mono_E)
        N_M = len(masses_all)

        direct_mat     = np.zeros((N_M, len(mono_E)))
        secondary_mat  = np.zeros_like(direct_mat)
        inflight_mat   = np.zeros_like(direct_mat)
        final_mat      = np.zeros_like(direct_mat)
        Emax_ifa       = np.zeros(N_M)

        for i, m in enumerate(masses_all):
            sub = os.path.join(DATA_DIR, names_all[i])
            S = load_spectra_components(sub)
            direct_mat[i]    = S['direct_gamma_primary']
            secondary_mat[i] = np.interp(mono_E, S['energy_secondary'], S['direct_gamma_secondary'], left=0, right=0)
            inflight_mat[i]  = S['IFA_primary'] + np.interp(mono_E, S['energy_secondary'], S['IFA_secondary'], left=0, right=0)
            final_mat[i]     = S['FSR_primary'] + np.interp(mono_E, S['energy_secondary'], S['FSR_secondary'], left=0, right=0)
            p = load_xy_lenient(os.path.join(sub, "inflight_annihilation_prim.txt"))
            s = load_xy_lenient(os.path.join(sub, "inflight_annihilation_sec.txt"))
            Emax_ifa[i] = max(p[:,0].max() if p.size else 0, s[:,0].max() if s.size else 0)

        logM_all = np.log(masses_all)
        ld = np.log(np.where(direct_mat    > floor, direct_mat,     floor))
        ls = np.log(np.where(secondary_mat > floor, secondary_mat,  floor))
        li = np.log(np.where(inflight_mat  > floor, inflight_mat,   floor))
        lf = np.log(np.where(final_mat     > floor, final_mat,      floor))

        sp_d = RectBivariateSpline(logM_all, logE, ld, kx=1, ky=3, s=0)
        sp_s = RectBivariateSpline(logM_all, logE, ls, kx=1, ky=3, s=0)
        sp_i = RectBivariateSpline(logM_all, logE, li, kx=1, ky=3, s=0)
        sp_f = RectBivariateSpline(logM_all, logE, lf, kx=1, ky=3, s=0)
        mono_ready = True

    while True:
        _print_menu()
        try:
            choice = user_input("Choice: ", allow_back=True, allow_exit=True).strip().lower()
        except BackRequested:
            return
        # Handle feeders that don't raise BackRequested / SystemExit
        if choice in ("b", "back"):
                    return
        if choice in ("q", "exit"):
            return  # don't sys.exit() here; returning keeps tests happy

        # ----- plot all queued -----
        if choice == '0':
            if not queue:
                warn("Queue is empty.")
                continue
            plot_pack = [(item['label'], (item['E'], item['S'])) for item in queue]
            _plot_dn(plot_pack)
            _plot_e2(plot_pack)

            for item in queue:
                h = item.get('hist')
                if not h:  # monochromatic has no histogram
                    continue
                kind = h['kind']
                run_dir = h['run_dir']
                label_for_hist = item['label']
                peak_val = h.get('peak')
                sigma_val = h.get('sigma')  # for non_gaussian this is sigma_X
                try:
                    if kind in ("gaussian", "non_gaussian", "lognormal"):
                        md_path = os.path.join(run_dir, "mass_distribution.txt")
                        md = np.loadtxt(md_path)
                        if kind == "gaussian":
                            _hist_gaussian(md,
                                           peak_val if peak_val is not None else np.median(md),
                                           sigma_val if sigma_val is not None else 0.05,
                                           label_for_hist)
                        elif kind == "non_gaussian":
                            _hist_nongaussian(md,
                                              peak_val if peak_val is not None else np.median(md),
                                              sigma_val if sigma_val is not None else 0.08,
                                              label_for_hist)
                        else:  # lognormal
                            _hist_lognormal(md,
                                            peak_val if peak_val is not None else np.median(md),
                                            sigma_val if sigma_val is not None else 1.0,
                                            label_for_hist)
                    elif kind == "custom":
                        _hist_custom(run_dir, label_for_hist)
                except FileNotFoundError as e:
                    warn(f"Histogram inputs missing in {run_dir}: {e}")
                except Exception as e:
                    warn(f"Histogram reconstruction failed for {run_dir}: {e}")

            queue.clear()
            info("Queue cleared.")
            continue

        if choice not in cat_map:
            warn("Invalid choice; try again.")
            continue

        label_group, root, spec_file, kind = cat_map[choice]

        # ----- NEW Monochromatic branch (no listing, no rounding; interpolate & queue) -----
        if kind == "mono":
            print(f"Enter PBH masses (g) to QUEUE for monochromatic plots (range [{M_MIN_MONO:.2e}, {M_MAX_MONO:.2e}]).")
            try:
                mstr = user_input("Masses (comma-separated): ", allow_back=True, allow_exit=True).strip()
            except BackRequested:
                continue
            if not mstr:
                continue

            # build interpolators once
            try:
                _ensure_mono_interpolator()
            except Exception as e:
                err(f"Cannot prepare interpolator: {e}")
                continue

            req_masses = parse_float_list_verbose(
                mstr, name="mass (g)", bounds=(M_MIN_MONO, M_MAX_MONO), allow_empty=False
            )
            if not req_masses:
                continue

            for mval in req_masses:
                try:
                    # no snapping: always interpolate in (logM, logE)
                    idx_up  = int(np.searchsorted(masses_all, mval, side='left'))
                    idx_low = max(0, idx_up-1)
                    idx_up  = min(idx_up, len(masses_all)-1)
                    Ecut    = min(Emax_ifa[idx_low], Emax_ifa[idx_up])
                    logm    = np.log(mval)
                    d_vals  = np.exp(sp_d(logm, logE, grid=False))
                    s_vals  = np.exp(sp_s(logm, logE, grid=False))
                    i_vals  = np.exp(sp_i(logm, logE, grid=False))
                    f_vals  = np.exp(sp_f(logm, logE, grid=False))
                    # guard inflight tails
                    for j in range(len(i_vals)-1,0,-1):
                        if np.isclose(i_vals[j], i_vals[j-1], rtol=1e-8): i_vals[j] = 0.0
                        else: break
                    log10i = np.log10(np.where(i_vals>0, i_vals, floor))
                    for j in range(1,len(log10i)):
                        if log10i[j] - log10i[j-1] < -50:
                            i_vals[j:] = 0.0; break
                    i_vals[mono_E >= Ecut] = 0.0
                    T = d_vals + s_vals + i_vals + f_vals
                    T[T < 1e-299] = 0.0

                    queue.append({
                        'label': f"Monochromatic {mval:.2e} g (interp)",
                        'E': mono_E.copy(),
                        'S': T.copy(),
                        'hist': None
                    })
                    info(f"Queued Monochromatic {mval:.2e} g (interpolated)")
                except Exception as e:
                    warn(f"Failed to queue {mval:.2e} g: {e}")
            continue

        # ----- Distributed branches (unchanged, except bigger hist funcs will be used later) -----
        try:
            subdirs = [
                d for d in sorted(os.listdir(root))
                if os.path.isdir(os.path.join(root, d))
            ]
        except FileNotFoundError:
            subdirs = []

        # Pretty listing entries
        pretty_entries = []
        for d in subdirs:
            peak_val, sigma_val, sigma_str = _extract_peak_sigma(d, kind)
            if kind == "gaussian":
                pretty = (f"{GAUSSIAN_METHOD} peak {peak_val:.2e} g ({sigma_str})"
                          if (peak_val is not None and sigma_str is not None) else f"{GAUSSIAN_METHOD} {d}")
            elif kind == "non_gaussian":
                pretty = (f"{NON_GAUSSIAN_METHOD} peak {peak_val:.2e} g ({sigma_str})"
                          if (peak_val is not None and sigma_str is not None) else f"{NON_GAUSSIAN_METHOD} {d}")
            elif kind == "lognormal":
                pretty = (f"{LOGNORMAL_METHOD} peak {peak_val:.2e} g ({sigma_str})"
                          if (peak_val is not None and sigma_str is not None) else f"{LOGNORMAL_METHOD} {d}")
            elif kind == "custom":
                pretty = d
            else:
                pretty = d
            pretty_entries.append((d, pretty, peak_val, sigma_val))

        print(f"Available in {label_group}:")
        for i, (_, pretty, _, _) in enumerate(pretty_entries, start=1):
            print(f" {i}: {pretty}")

        sel = user_input(
            "Enter indices to QUEUE (comma-separated): ",
            allow_back=True, allow_exit=True
        ).strip()
        if not sel:
            continue

        try:
            idxs = [int(x) for x in sel.split(",") if x.strip()]
        except Exception:
            warn("Invalid indices input.")
            continue

        for i_sel in idxs:
            if not (1 <= i_sel <= len(pretty_entries)):
                warn(f"Index {i_sel} out of range.")
                continue

            run_name, pretty_label, peak_val, sigma_val = pretty_entries[i_sel - 1]
            run_dir  = os.path.join(root, run_name)
            spec_path = os.path.join(run_dir, spec_file) if spec_file else None

            try:
                if spec_path and os.path.isfile(spec_path):
                    arr = np.loadtxt(spec_path)
                    if arr.ndim >= 2 and arr.shape[1] >= 2:
                        E = arr[:,0]
                        S = arr[:,1]
                        if kind == "gaussian":
                            plot_label = (f"{GAUSSIAN_METHOD} peak {peak_val:.2e} g (σ={sigma_val:.3g})"
                                          if peak_val is not None and sigma_val is not None else
                                          f"{GAUSSIAN_METHOD} {run_name}")
                        elif kind == "non_gaussian":
                            plot_label = (f"{NON_GAUSSIAN_METHOD} peak {peak_val:.2e} g (σX={sigma_val:.3g})"
                                          if peak_val is not None and sigma_val is not None else
                                          f"{NON_GAUSSIAN_METHOD} {run_name}")
                        elif kind == "lognormal":
                            plot_label = (f"{LOGNORMAL_METHOD} peak {peak_val:.2e} g (σ={sigma_val:.3g})"
                                          if peak_val is not None and sigma_val is not None else
                                          f"{LOGNORMAL_METHOD} {run_name}")
                        elif kind == "custom":
                            plot_label = run_name
                        else:
                            plot_label = run_name

                        queue.append({
                            'label': plot_label,
                            'E': E,
                            'S': S,
                            'hist': {
                                'kind': kind,
                                'run_dir': run_dir,
                                'peak': peak_val,
                                'sigma': sigma_val
                            }
                        })
                        info(f"Queued {plot_label}")
                    else:
                        warn(f"{spec_file} malformed in {run_name}.")
                else:
                    warn(f"No '{spec_file}' found in {run_name}; skipping queue.")
            except Exception as e:
                warn(f"Failed to queue {run_name}: {e}")

    # ---------- UI ----------
    def _print_menu():
        print("\nView Previous — choose:")
        print(" 1: Monochromatic Distribution")
        print(f" 2: {GAUSSIAN_METHOD}")
        print(f" 3: {NON_GAUSSIAN_METHOD}")
        print(f" 4: {LOGNORMAL_METHOD}")
        print(f" 5: Custom equation (user-defined mass PDF)")
        print(" 0: Plot all Queued | b: Back | q: Quit")

    # queue elements:
    # {
    #    'label': str (pretty label for overlay / hist title),
    #    'E': ndarray,
    #    'S': ndarray,
    #    'hist': {'kind': kind, 'run_dir': run_dir, 'peak': float?, 'sigma': float?} or None
    # }
    queue: list[dict] = []

    while True:
        _print_menu()
        try:
            choice = user_input("Choice: ", allow_back=True, allow_exit=True).strip().lower()
        except BackRequested:
            return

        # ----- plot all queued -----
        if choice == '0':
            if not queue:
                warn("Queue is empty.")
                continue

            # spectra first
            plot_pack = [(item['label'], (item['E'], item['S'])) for item in queue]
            _plot_dn(plot_pack)
            _plot_e2(plot_pack)

            # histograms second
            for item in queue:
                h = item.get('hist')
                if not h:
                    continue
                kind = h['kind']
                run_dir = h['run_dir']
                label_for_hist = item['label']
                peak_val = h.get('peak')
                sigma_val = h.get('sigma')  # for non_gaussian this is sigma_X

                try:
                    if kind in ("gaussian", "non_gaussian", "lognormal"):
                        md_path = os.path.join(run_dir, "mass_distribution.txt")
                        md = np.loadtxt(md_path)

                        if kind == "gaussian":
                            _hist_gaussian(md,
                                           peak_val if peak_val is not None else np.median(md),
                                           sigma_val if sigma_val is not None else 0.05,
                                           label_for_hist)

                        elif kind == "non_gaussian":
                            _hist_nongaussian(md,
                                              peak_val if peak_val is not None else np.median(md),
                                              sigma_val if sigma_val is not None else 0.08,
                                              label_for_hist)

                        else:  # lognormal
                            _hist_lognormal(md,
                                            peak_val if peak_val is not None else np.median(md),
                                            sigma_val if sigma_val is not None else 1.0,
                                            label_for_hist)

                    elif kind == "custom":
                        _hist_custom(run_dir, label_for_hist)

                except FileNotFoundError as e:
                    warn(f"Histogram inputs missing in {run_dir}: {e}")
                except Exception as e:
                    warn(f"Histogram reconstruction failed for {run_dir}: {e}")

            # clear queue
            queue.clear()
            info("Queue cleared.")
            continue

        if choice not in cat_map:
            warn("Invalid choice; try again.")
            continue

        label_group, root, spec_file, kind = cat_map[choice]

        # ----- Monochromatic branch -----
        if kind == "mono":
            runs = []
            try:
                for fn in sorted(os.listdir(root)):
                    if fn.lower().endswith(".txt"):
                        runs.append(("file", fn))
            except FileNotFoundError:
                pass

            print(f"Available in {label_group}:")
            for i,(_,fn) in enumerate(runs, start=1):
                print(f" {i}: {fn}")

            print("\nYou can also request a target mass (in grams) to generate the nearest pre-rendered mono spectrum.")
            sel = user_input(
                "Enter indices to QUEUE (comma-separated) OR a mass (e.g. 1e15), or Enter to cancel: ",
                allow_back=True, allow_exit=True
            ).strip()
            if not sel:
                continue

            # numeric mass path?
            try:
                mass_try = float(sel)
            except Exception:
                mass_try = None

            if mass_try is not None:
                if not (M_MIN_MONO <= mass_try <= M_MAX_MONO):
                    warn(f"Mass outside allowed view window [{M_MIN_MONO:.2e}, {M_MAX_MONO:.2e}].")
                    continue
                try:
                    fname = generate_monochromatic_for_mass(mass_try, DATA_DIR, MONO_RESULTS_DIR)
                    arr = np.loadtxt(fname)
                    E = arr[:,0]
                    T = arr[:,1]
                    queue.append({
                        'label': f"Monochromatic {mass_try:.2e} g",
                        'E': E,
                        'S': T,
                        'hist': None
                    })
                    info(f"Queued Monochromatic {mass_try:.2e} g → {os.path.basename(fname)}")
                except Exception as e:
                    err(f"Could not generate/queue mono spectrum: {e}")
                continue

            # index list path
            try:
                idxs = [int(x) for x in sel.split(",") if x.strip()]
            except Exception:
                warn("Invalid indices input.")
                continue

            for i in idxs:
                if 1 <= i <= len(runs):
                    _, fn = runs[i-1]
                    path = os.path.join(root, fn)
                    try:
                        arr = np.loadtxt(path)
                        E = arr[:,0]
                        T = arr[:,1] if arr.ndim > 1 and arr.shape[1] >= 2 else np.zeros_like(E)
                        queue.append({
                            'label': f"Monochromatic {fn}",
                            'E'   : E,
                            'S'   : T,
                            'hist': None
                        })
                        info(f"Queued {fn}")
                    except Exception as e:
                        warn(f"Failed to read {fn}: {e}")
            continue

        # ----- Distributed branches -----
        try:
            subdirs = [
                d for d in sorted(os.listdir(root))
                if os.path.isdir(os.path.join(root, d))
            ]
        except FileNotFoundError:
            subdirs = []

        # Build pretty listing entries
        pretty_entries = []
        for d in subdirs:
            peak_val, sigma_val, sigma_str = _extract_peak_sigma(d, kind)

            if kind == "gaussian":
                if peak_val is not None and sigma_str is not None:
                    pretty = f"{GAUSSIAN_METHOD} peak {peak_val:.2e} g ({sigma_str})"
                else:
                    pretty = f"{GAUSSIAN_METHOD} {d}"

            elif kind == "non_gaussian":
                if peak_val is not None and sigma_str is not None:
                    pretty = f"{NON_GAUSSIAN_METHOD} peak {peak_val:.2e} g ({sigma_str})"
                else:
                    pretty = f"{NON_GAUSSIAN_METHOD} {d}"

            elif kind == "lognormal":
                if peak_val is not None and sigma_str is not None:
                    pretty = f"{LOGNORMAL_METHOD} peak {peak_val:.2e} g ({sigma_str})"
                else:
                    pretty = f"{LOGNORMAL_METHOD} {d}"

            elif kind == "custom":
                # Just show folder name (no equation)
                pretty = d

            else:
                pretty = d

            pretty_entries.append((d, pretty, peak_val, sigma_val))

        print(f"Available in {label_group}:")
        for i, (_, pretty, _, _) in enumerate(pretty_entries, start=1):
            print(f" {i}: {pretty}")

        sel = user_input(
            "Enter indices to QUEUE (comma-separated): ",
            allow_back=True, allow_exit=True
        ).strip()
        if not sel:
            continue

        try:
            idxs = [int(x) for x in sel.split(",") if x.strip()]
        except Exception:
            warn("Invalid indices input.")
            continue

        for i_sel in idxs:
            if not (1 <= i_sel <= len(pretty_entries)):
                warn(f"Index {i_sel} out of range.")
                continue

            run_name, pretty_label, peak_val, sigma_val = pretty_entries[i_sel - 1]
            run_dir  = os.path.join(root, run_name)
            spec_path = os.path.join(run_dir, spec_file) if spec_file else None

            try:
                if spec_path and os.path.isfile(spec_path):
                    arr = np.loadtxt(spec_path)
                    if arr.ndim >= 2 and arr.shape[1] >= 2:
                        E = arr[:,0]
                        S = arr[:,1]

                        # Build final label for plots:
                        if kind == "gaussian":
                            if peak_val is not None and sigma_val is not None:
                                plot_label = f"{GAUSSIAN_METHOD} peak {peak_val:.2e} g (σ={sigma_val:.3g})"
                            else:
                                plot_label = f"{GAUSSIAN_METHOD} {run_name}"

                        elif kind == "non_gaussian":
                            if peak_val is not None and sigma_val is not None:
                                plot_label = f"{NON_GAUSSIAN_METHOD} peak {peak_val:.2e} g (σX={sigma_val:.3g})"
                            else:
                                plot_label = f"{NON_GAUSSIAN_METHOD} {run_name}"

                        elif kind == "lognormal":
                            if peak_val is not None and sigma_val is not None:
                                plot_label = f"{LOGNORMAL_METHOD} peak {peak_val:.2e} g (σ={sigma_val:.3g})"
                            else:
                                plot_label = f"{LOGNORMAL_METHOD} {run_name}"

                        elif kind == "custom":
                            plot_label = run_name

                        else:
                            plot_label = run_name

                        queue.append({
                            'label': plot_label,
                            'E': E,
                            'S': S,
                            'hist': {
                                'kind': kind,
                                'run_dir': run_dir,
                                'peak': peak_val,
                                'sigma': sigma_val
                            }
                        })
                        info(f"Queued {plot_label}")
                    else:
                        warn(f"{spec_file} malformed in {run_name}.")
                else:
                    warn(f"No '{spec_file}' found in {run_name}; skipping queue.")
            except Exception as e:
                warn(f"Failed to queue {run_name}: {e}")


# ---------------------------
# UI
# ---------------------------
from colorama import Fore, Style, init as colorama_init
colorama_init(autoreset=True)

def show_start_screen() -> None:
    """
    Print the program banner and helpful usage hints.
    """
    width = 56  # inner width of the box
    top = "╔" + "═" * width + "╗"
    bot = "╚" + "═" * width + "╝"
    title = "GammaPBHPlotter: PBH Spectrum Tool"
    ver   = f"Version {__version__}"

    print("\n" + Fore.CYAN + Style.BRIGHT + top)
    print(        Fore.CYAN + Style.BRIGHT + f"║{title.center(width)}║")
    print(        Fore.CYAN + Style.BRIGHT + f"║{ver.center(width)}║")
    print(        Fore.CYAN + Style.BRIGHT + bot + Style.RESET_ALL)
    print()
    print("Analyze and visualize Hawking radiation spectra of primordial black holes.\n")
    print(Fore.YELLOW + "📄 Associated Publication:" + Style.RESET_ALL)
    print("   John Carlini & Ilias Cholis — Particle Astrophysics Research\n")
    print("At any prompt: 'b' = back, 'q' = quit.")

def main() -> None:
    """
    Entry point for the interactive CLI loop.

    Menu
    ----
    1: Monochromatic spectra
    2: Distributed spectra (Gaussian collapse)
    3: Distributed spectra (Non-Gaussian Collapse)
    4: Distributed spectra (Log-Normal Distribution)
    5: Distributed spectra (Custom mass PDF)
    6: View previous spectra
    0: Exit
    """
    show_start_screen()
    while True:
        print("\nSelect:")
        print("1: Monochromatic spectra")
        print(f"2: Distributed spectra ({GAUSSIAN_METHOD})")
        print(f"3: Distributed spectra ({NON_GAUSSIAN_METHOD})")
        print(f"4: Distributed spectra ({LOGNORMAL_METHOD})")
        print("5: Distributed spectra (Custom mass PDF)")
        print("6: View previous spectra")
        print("0: Exit")
        choice = user_input("Choice: ", allow_back=False, allow_exit=True).strip().lower()
        if choice == '1':
            monochromatic_spectra()
        elif choice == '2':
            distributed_spectrum(GAUSSIAN_METHOD)
        elif choice == '3':
            distributed_spectrum(NON_GAUSSIAN_METHOD)
        elif choice == '4':
            distributed_spectrum(LOGNORMAL_METHOD)
        elif choice == '5':
            custom_equation_pdf_tool()
        elif choice == '6':
            view_previous_spectra()
        elif choice in ['0','exit','q']:
            print("Goodbye.")
            break
        else:
            print("Invalid; try again.")


if __name__ == '__main__':
    try:
        main()
    except BackRequested:
        # If a BackRequested was thrown at the top-level, just exit cleanly.
        pass
    except Exception:
        import traceback
        traceback.print_exc()
        try:
            input("\nAn error occurred. Press Enter to exit…")
        except Exception:
            pass
