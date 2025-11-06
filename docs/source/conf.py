#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
- Spectral units are dN/dE [MeV^-1 s^-1].
- Internal results directories are created under this package’s folder.
- No writes are attempted outside the installed package, by design.

Safety/UX
---------
- Robust input parsing with back/quit keywords supported in every prompt.
- Gentle warnings (not hard failures) when tokens are malformed or out of range.
- Histograms use log-spaced bins for mass distributions by default, and
  analytic overlays are scaled to *expected counts per log-bin* so the line
  matches the histogram’s geometry.
"""

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

# ---------------------------
# Matplotlib/NumPy basics
# ---------------------------
plt.rcParams.update({'font.size': 12})
np.seterr(divide='ignore', invalid='ignore')  # suppress log/invalid warnings


# ---------------------------
# Paths (package-internal only)
# ---------------------------
def _resolve_data_dir() -> str:
    """
    Compute the absolute path to the packaged BlackHawk data tables.

    Returns
    -------
    str
        Filesystem path to the `blackhawk_data/` directory shipped with the
        installed package. No external data directories are used.

    Notes
    -----
    This function assumes the module file layout:

        src/gammapbh/cli.py
        src/gammapbh/blackhawk_data/<mass>/...

    The directory is *not* created here; it must exist in the installation.
    """
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(pkg_dir, "blackhawk_data")


def _resolve_results_root() -> str:
    """
    Create (if necessary) and return the package-internal results directory.

    Returns
    -------
    str
        Filesystem path to `results/` under the package folder.

    Raises
    ------
    RuntimeError
        If the directory cannot be created or written to.

    Notes
    -----
    A small writability check is performed by writing and removing a temp file.
    """
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    dest = os.path.join(pkg_dir, "results")
    os.makedirs(dest, exist_ok=True)
    try:
        test = os.path.join(dest, ".writetest.tmp")
        with open(test, "w", encoding="utf-8") as fh:
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

for d in (MONO_RESULTS_DIR, CUSTOM_RESULTS_DIR, GAUSS_RESULTS_DIR, NGAUSS_RESULTS_DIR, LOGN_RESULTS_DIR):
    os.makedirs(d, exist_ok=True)


# ---------------------------
# Labels
# ---------------------------
GAUSSIAN_METHOD     = "Gaussian collapse"
NON_GAUSSIAN_METHOD = "Non-Gaussian Collapse"
LOGNORMAL_METHOD    = "Log-Normal Distribution"


# ---------------------------
# Helper: required files in each mass folder
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
    """
    Internal control-flow exception raised when a user indicates they want
    to return to the prior menu.

    Trigger words
    -------------
    - 'b' or 'back' at any input prompt that has `allow_back=True`.
    """


def discover_mass_folders(data_dir):
    """
    Discover available PBH mass folders that contain all required component files.

    Parameters
    ----------
    data_dir : str
        Root directory containing one subdirectory per mass value (folder name
        should be convertible to float).

    Returns
    -------
    tuple[list[float], list[str]]
        Sorted (ascending) list of masses and corresponding folder names. If no
        valid folders exist, returns ([], []).

    Notes
    -----
    A folder is considered valid if:
      1) its name parses as a float (the mass in grams), and
      2) it contains every filename listed in `REQUIRED_FILES`.
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
def info(msg):
    """
    Print an informational message in cyan.

    Parameters
    ----------
    msg : str
        Message body.
    """
    print(Fore.CYAN + "ℹ " + msg + Style.RESET_ALL)


def warn(msg):
    """
    Print a warning message in yellow.

    Parameters
    ----------
    msg : str
        Message body.
    """
    print(Fore.YELLOW + "⚠ " + msg + Style.RESET_ALL)


def err(msg):
    """
    Print an error message in red.

    Parameters
    ----------
    msg : str
        Message body.
    """
    print(Fore.RED + "✖ " + msg + Style.RESET_ALL)


def user_input(prompt, *, allow_back=False, allow_exit=True):
    """
    Read a single line of input with support for 'back'/'exit' control words.

    Parameters
    ----------
    prompt : str
        Text shown to the user.
    allow_back : bool, optional
        If True, 'b' or 'back' raises BackRequested to unwind one menu level.
    allow_exit : bool, optional
        If True, 'q' or 'exit' terminates the program via sys.exit(0).

    Returns
    -------
    str
        The raw input string (stripped).

    Raises
    ------
    BackRequested
        If the user typed 'b' or 'back' and `allow_back=True`.

    Notes
    -----
    This wrapper consolidates input semantics used throughout the CLI and
    ensures consistent behavior on 'back'/'quit'.
    """
    txt = input(prompt).strip()
    low = txt.lower()
    if allow_exit and low in ('exit', 'q'):
        print("Exiting software.")
        sys.exit(0)
    if allow_back and low in ('b', 'back'):
        raise BackRequested()
    return txt


def list_saved_runs(base_dir):
    """
    Enumerate saved run subdirectories.

    Parameters
    ----------
    base_dir : str
        Parent directory (e.g., results/gaussian).

    Returns
    -------
    list[str]
        Sorted list of child directory names; empty if missing/not readable.
    """
    try:
        return sorted(d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)))
    except FileNotFoundError:
        return []


def snap_to_available(mval, available, tol=1e-12):
    """
    Snap a requested mass to the nearest *exact* available mass in log-space.

    Parameters
    ----------
    mval : float
        Requested mass [g].
    available : array_like
        Iterable of available masses [g].
    tol : float, optional
        Maximum abs(log(m_available) - log(mval)) to be considered a snap-match.


    Returns
    -------
    float or None
        Nearest available mass if within `tol` in log-space; otherwise None.

    Notes
    -----
    This avoids surprises when a data folder exactly matches a requested mass,
    while letting us fall back to interpolation when needed.
    """
    log_m = np.log(mval)
    log_available = np.log(np.array(available))
    diffs = np.abs(log_available - log_m)
    idx = np.argmin(diffs)
    return available[idx] if diffs[idx] < tol else None


def parse_float_list_verbose(
    s, *, name="value", bounds=None, allow_empty=False,
    positive_only=False, strict_gt=False, strict_lt=False
):
    """
    Parse a comma-separated list of floats with friendly validation.

    Parameters
    ----------
    s : str or None
        Raw input string (e.g., "1, 2.5, 3e4").
    name : str, optional
        Logical name used in warning messages (e.g., "σ").
    bounds : tuple[float|None, float|None] or None
        Inclusive bounds (lo, hi). Pass None to disable one side.
    allow_empty : bool, optional
        If False and parsing yields nothing, a warning is printed.
    positive_only : bool, optional
        If True, values must be strictly greater than 0.
    strict_gt : bool, optional
        If True and `bounds[0]` is not None, enforce v > lo (exclusive).
        Otherwise v ≥ lo (inclusive).
    strict_lt : bool, optional
        If True and `bounds[1]` is not None, enforce v < hi (exclusive).
        Otherwise v ≤ hi (inclusive).

    Returns
    -------
    list[float]
        Parsed, de-duplicated values that passed validation (order preserved).
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
def delta_l(mass_ratio, kappa, delta_c, gamma):
    """
    Mapping from mass ratio to the density-contrast parameter δ_l (auxiliary).

    Parameters
    ----------
    mass_ratio : array_like
        Ratio r ≡ M / M_peak (dimensionless).
    kappa : float
        Collapse threshold scaling (model parameter).
    delta_c : float
        Critical collapse threshold δ_c (model parameter).
    gamma : float
        Scaling exponent γ (model parameter).

    Returns
    -------
    ndarray
        δ_l(r) evaluated elementwise, with the inner square-root argument
        clamped to ≥ 0 for numerical stability.

    Notes
    -----
    The functional form here is tuned to the Gaussian/Non-Gaussian collapse
    prescriptions used elsewhere in this tool. The clipping step avoids NaNs.
    """
    y = (mass_ratio / kappa)**(1.0 / gamma)
    arg = 64 - 96 * (delta_c + y)
    arg = np.clip(arg, 0.0, None)
    return (8 - np.sqrt(arg)) / 6


def mass_function(delta_l_val, sigma_x, delta_c, gamma):
    """
    Gaussian-collapse mass function (up to normalization on a grid).

    Parameters
    ----------
    delta_l_val : array_like
        δ_l evaluated at r-grid.
    sigma_x : float
        Effective variance parameter σ_x controlling distribution width.
    delta_c : float
        Critical threshold δ_c.
    gamma : float
        Scaling exponent γ.

    Returns
    -------
    ndarray
        Unnormalized φ(r) (positive where physically meaningful).

    Notes
    -----
    The final sampling PDF is built by normalizing this curve on a finite r-grid
    and converting to mass-space via a scale factor tied to the mode.
    """
    term1 = 1.0 / (np.sqrt(2 * np.pi) * sigma_x)
    term2 = np.exp(-delta_l_val**2 / (2 * sigma_x**2))
    term3 = delta_l_val - (3/8) * delta_l_val**2 - delta_c
    term4 = gamma * np.abs(1 - (3/4) * delta_l_val)
    return term1 * term2 * term3 / term4


def mass_function_exact(delta_l_val, sigma_X, sigma_Y, delta_c, gamma):
    """
    Non-Gaussian collapse mass function following Biagetti et al. (Eq. 20 proxy).

    Parameters
    ----------
    delta_l_val : array_like
        δ_l evaluated at r-grid.
    sigma_X : float
        Variance parameter σ_X.
    sigma_Y : float
        Variance parameter σ_Y (often tied to σ_X by a fixed ratio).
    delta_c : float
        Critical threshold δ_c.
    gamma : float
        Scaling exponent γ.

    Returns
    -------
    ndarray
        Unnormalized φ(r) (positive where defined).

    Notes
    -----
    This function mirrors the algebraic structure in Biagetti et al., including
    an explicit Jacobian for variable change. It remains unnormalized until
    converted to a mass-PDF on a finite grid.
    """
    # Biagetti et al. Eq. (20) — compiled terms for numerical stability.
    A = sigma_X**2 + (sigma_Y * delta_l_val)**2
    exp_pref = np.exp(-1.0 / (2.0 * sigma_Y**2))
    term1 = 2.0 * sigma_Y * np.sqrt(A)
    inner_exp = np.exp(sigma_X**2 / (2.0 * sigma_Y**2 * (sigma_X**2 + 2.0 * (sigma_Y * delta_l_val)**2)))
    erf_arg = sigma_X * np.sqrt(2.0) / np.sqrt(A)  # well-conditioned
    term2 = np.sqrt(2.0 * np.pi) * sigma_X * inner_exp * erf(erf_arg)
    bracket = term1 + term2
    norm = exp_pref * sigma_X / (2.0 * np.pi * A**1.5)
    jacobian = ((delta_l_val - 0.375 * delta_l_val**2 - delta_c) /
                (gamma * np.abs(1.0 - 0.75 * delta_l_val)))
    return norm * bracket * jacobian


def mass_function_lognormal(x, mu, sigma):
    """
    Compute the log-normal PBH mass function φ(x; μ, σ) on a ratio/grid.

    Parameters
    ----------
    x : array_like
        Positive support variable (e.g., mass or ratio). Values ≤0 are clipped.
    mu : float
        Log-mean in natural logarithm space (μ = E[ln x]).
    sigma : float
        Log-standard deviation in ln-space (σ = Std[ln x]).

    Returns
    -------
    ndarray
        φ(x) evaluated at each x.

    Notes
    -----
    This is the standard continuous log-normal density::

        φ(x) = 1 / (x σ sqrt(2π)) * exp( - (ln x - μ)^2 / (2 σ^2) )

    The output is not normalized over arbitrary slices unless integrated across
    the full positive support or normalized on a finite grid by the caller.

    """
    x_clipped = np.clip(x, 1e-16, None)
    return (1.0 / (x_clipped * sigma * np.sqrt(2 * np.pi))
            * np.exp(- (np.log(x_clipped) - mu)**2 / (2 * sigma**2)))


# ---------------------------
# Data loaders
# ---------------------------
def load_data(filepath, skip_header=0):
    """
    Load a whitespace-delimited numeric table as a NumPy array.

    Parameters
    ----------
    filepath : str
        Path to the text file.
    skip_header : int, optional
        Lines to skip at the beginning (comments/headers).

    Returns
    -------
    ndarray
        2D array of floats.

    Raises
    ------
    FileNotFoundError
        If the file is missing.
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    return np.genfromtxt(filepath, skip_header=skip_header)


def load_spectra_components(directory):
    """
    Load all spectral component tables for a specific mass folder.

    Parameters
    ----------
    directory : str
        Folder that contains:
          - instantaneous_primary_spectra.txt
          - instantaneous_secondary_spectra.txt
          - inflight_annihilation_prim.txt
          - inflight_annihilation_sec.txt
          - final_state_radiation_prim.txt
          - final_state_radiation_sec.txt

    Returns
    -------
    dict
        Keys are documented as::

            {
              'energy_primary'        : ndarray (MeV),
              'energy_secondary'      : ndarray (MeV),
              'direct_gamma_primary'  : ndarray,
              'direct_gamma_secondary': ndarray,
              'IFA_primary'           : ndarray (on primary energy grid),
              'IFA_secondary'         : ndarray (on secondary grid, then regridded),
              'FSR_primary'           : ndarray (on primary energy grid),
              'FSR_secondary'         : ndarray (on secondary grid, then regridded),
            }


    Notes
    -----
    - Primary energy grid is converted from GeV→MeV (×1e3).
    - Secondary components are interpolated to the primary grid where needed.
    """
    primary   = load_data(os.path.join(directory, "instantaneous_primary_spectra.txt"), skip_header=2)[123:]
    secondary = load_data(os.path.join(directory, "instantaneous_secondary_spectra.txt"), skip_header=1)
    IFA_prim  = load_data(os.path.join(directory, "inflight_annihilation_prim.txt"))
    IFA_sec   = load_data(os.path.join(directory, "inflight_annihilation_sec.txt"))
    FSR_prim  = load_data(os.path.join(directory, "final_state_radiation_prim.txt"), skip_header=1)
    FSR_sec   = load_data(os.path.join(directory, "final_state_radiation_sec.txt"),  skip_header=1)

    E_prim = primary[:,0] * 1e3  # MeV
    E_sec  = secondary[:,0]      # MeV

    return {
        'energy_primary':         E_prim,
        'energy_secondary':       E_sec,
        'direct_gamma_primary':   primary[:,1] / 1e3,
        'direct_gamma_secondary': secondary[:,1],
        'IFA_primary':            np.interp(E_prim, IFA_prim[:,0], IFA_prim[:,1], left=0.0, right=0.0),
        'IFA_secondary':          np.interp(E_sec,  IFA_sec[:,0],  IFA_sec[:,1],  left=0.0, right=0.0),
        'FSR_primary':            np.interp(E_prim, FSR_prim[:,0], FSR_prim[:,1]),
        'FSR_secondary':          np.interp(E_sec,  FSR_sec[:,0],  FSR_sec[:,1]),
    }
# ---------------------------
# Global numeric guards
# ---------------------------
MASS_MIN = 5.0e13
MASS_MAX = 1.0e19

# Cap for number of histogram bins for distributed spectra
MAX_HIST_BINS = 20


# ---------------------------
# Spectrum assemblers & splines
# ---------------------------
def _assemble_total_on_primary(comp):
    """
    Build a *total* instantaneous spectrum on the primary energy grid.

    Parameters
    ----------
    comp : dict
        Output of `load_spectra_components`.

    Returns
    -------
    (E, total) : (ndarray, ndarray)
        Primary-grid energies [MeV] and total differential spectrum dN/dE.

    Notes
    -----
    Secondary components are interpolated onto the primary grid.
    """
    E_p = comp['energy_primary']
    E_s = comp['energy_secondary']

    direct_p = comp['direct_gamma_primary']
    direct_s_on_p = np.interp(E_p, E_s, comp['direct_gamma_secondary'], left=0.0, right=0.0)

    IFA_p = comp['IFA_primary']
    FSR_p = comp['FSR_primary']

    IFA_s_on_p = np.interp(E_p, E_s, comp['IFA_secondary'], left=0.0, right=0.0)
    FSR_s_on_p = np.interp(E_p, E_s, comp['FSR_secondary'], left=0.0, right=0.0)

    total = direct_p + direct_s_on_p + IFA_p + IFA_s_on_p + FSR_p + FSR_s_on_p
    total = np.clip(total, 0.0, None)  # safety
    return E_p, total


def build_bivariate_spline(available_masses, folder_names):
    """
    Construct RectBivariateSpline splines in (logM, logE) for the *total* spectrum.

    Parameters
    ----------
    available_masses : list[float]
        Valid masses (ascending) with complete component files.
    folder_names : list[str]
        Matching folder names (strings of those masses).

    Returns
    -------
    SimpleNamespace
        Fields:
          - logM : ndarray, sorted log-mass grid
          - logE : ndarray, common log-energy grid [MeV]
          - spline_total : RectBivariateSpline over (logM, logE) → log(dN/dE)
          - E_ref : ndarray, the common energy grid [MeV]

    Notes
    -----
    We:
      1) Choose a reference energy grid from the *median* mass folder's primary grid.
      2) Interpolate each mass's *total* spectrum to that grid.
      3) Spline log(dN/dE) for better smoothness across decades.
    """
    if not available_masses:
        raise RuntimeError("No mass folders found. Ensure blackhawk_data/ is populated.")

    # Use median folder as the energy reference
    mid_idx = len(available_masses) // 2
    ref_dir = os.path.join(DATA_DIR, folder_names[mid_idx])
    ref_comp = load_spectra_components(ref_dir)
    E_ref, _ = _assemble_total_on_primary(ref_comp)
    logE = np.log(E_ref)

    # Stack log(total) across all masses on E_ref
    logM = np.log(np.array(available_masses, dtype=float))
    Z = np.zeros((len(available_masses), len(E_ref)), dtype=float)

    for i, (m, fname) in enumerate(zip(available_masses, folder_names)):
        comp = load_spectra_components(os.path.join(DATA_DIR, fname))
        E_i, tot_i = _assemble_total_on_primary(comp)
        tot_on_ref = np.interp(E_ref, E_i, tot_i, left=0.0, right=0.0)
        Z[i, :] = np.log(np.clip(tot_on_ref, 1e-300, None))  # avoid -inf

    spline_total = RectBivariateSpline(logM, logE, Z, kx=1, ky=3)  # linear in logM, cubic in logE

    return SimpleNamespace(
        logM=logM, logE=logE, spline_total=spline_total, E_ref=E_ref
    )


def eval_total_spectrum(spl, mass):
    """
    Evaluate the total instantaneous dN/dE [MeV^-1 s^-1] at an arbitrary mass.

    Parameters
    ----------
    spl : SimpleNamespace
        Output of `build_bivariate_spline`.
    mass : float
        PBH mass [g].

    Returns
    -------
    (E, dNdE) : (ndarray, ndarray)
        Energy grid [MeV] and spectrum.

    Notes
    -----
    Outside the training mass range, we clamp to the nearest edge (no extrap).
    """
    logM_in = np.log(np.clip(mass, np.exp(spl.logM[0]), np.exp(spl.logM[-1])))
    logE = spl.logE
    log_dndE = spl.spline_total(logM_in, logE, grid=False)
    dndE = np.exp(log_dndE)
    return np.exp(logE*0) * spl.E_ref, dndE


# ---------------------------
# Plotting helpers
# ---------------------------
def _pretty_mass(m):
    return f"{m:.2e} g"


def plot_monochromatic(masses, spl, *, title_suffix="(pre-rendered)"):
    """
    Plot one or more monochromatic spectra (evaluated from splines).

    Parameters
    ----------
    masses : list[float]
        Masses to plot.
    spl : SimpleNamespace
        Spline pack from `build_bivariate_spline`.
    title_suffix : str
        Extra text appended to the title.

    Returns
    -------
    list[tuple[float, ndarray, ndarray]]
        Per-mass triplets of (mass, E, dNdE).
    """
    plt.figure()
    out = []
    for m in masses:
        E, dNdE = eval_total_spectrum(spl, m)
        plt.loglog(E, np.clip(dNdE, 1e-300, None), label=_pretty_mass(m))
        out.append((m, E, dNdE))
    plt.xlabel("E [MeV]")
    plt.ylabel("dN/dE [MeV$^{-1}$ s$^{-1}$]")
    plt.title(f"Monochromatic spectra {title_suffix}")
    plt.grid(True, which='both', alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.show()
    return out


def _auto_log_bins(data, max_bins=MAX_HIST_BINS):
    """
    Compute log-spaced histogram edges with a hard cap on bin count.

    Parameters
    ----------
    data : array_like
        Positive data to bin (masses).
    max_bins : int
        Maximum number of bins (inclusive).

    Returns
    -------
    ndarray
        Log-spaced edges.

    Notes
    -----
    Uses Freedman–Diaconis estimate on log10(data) with a minimum of 8 bins,
    then clamps to `max_bins`. This preserves the *logic* of variability
    while honoring the new cap requested.
    """
    x = np.asarray(data, dtype=float)
    x = x[x > 0]
    if x.size < 2:
        return np.geomspace(MASS_MIN, MASS_MAX, 8)

    y = np.log10(x)
    q25, q75 = np.percentile(y, [25, 75])
    iqr = max(q75 - q25, 1e-9)
    h = 2 * iqr * (x.size ** (-1/3))  # FD in log-space
    if h <= 0:
        nb = 8
    else:
        span = y.max() - y.min()
        nb = int(np.ceil(span / h))
        nb = max(nb, 8)
    nb = min(nb, max_bins)
    return np.geomspace(x.min(), x.max(), nb + 1)


def _scale_pdf_to_hist_counts(edges, pdf_m, pdf_vals, N_total):
    """
    Convert a *mass-space* PDF into expected per-bin counts for a histogram.

    Parameters
    ----------
    edges : ndarray
        Histogram bin edges in mass (monotonic).
    pdf_m : ndarray
        Grid in mass where `pdf_vals` is evaluated (positive).
    pdf_vals : ndarray
        Nonnegative PDF values on `pdf_m` such that ∫ pdf dm = 1.
    N_total : int
        Number of samples shown in the histogram.

    Returns
    -------
    counts_per_bin : ndarray
        Expected counts for each bin (len(edges)-1).

    Notes
    -----
    We integrate the PDF across each bin using trapezoidal rule on a *refined*
    grid made by combining `edges` and `pdf_m`, then multiply by N_total.
    """
    pdf_m = np.asarray(pdf_m, float)
    pdf_vals = np.asarray(pdf_vals, float)
    # guard
    mask = (pdf_m > 0) & np.isfinite(pdf_vals) & (pdf_vals >= 0)
    pdf_m = pdf_m[mask]
    pdf_vals = pdf_vals[mask]
    if pdf_m.size < 2:
        return np.zeros(len(edges) - 1)

    # Build a union grid for robust bin integrals
    grid = np.unique(np.concatenate([pdf_m, edges]))
    vals = np.interp(grid, pdf_m, pdf_vals, left=0.0, right=0.0)

    # CDF on union grid
    cdf = np.concatenate([[0.0], np.cumsum(trapezoid([vals[:-1], vals[1:]], x=[grid[:-1], grid[1:]], axis=0))])
    cdf = cdf / max(cdf[-1], 1e-300)

    # Expected counts per bin = N * ΔCDF
    # Map edges to CDF via interpolation on the union grid
    cdf_on_edges = np.interp(edges, grid, cdf)
    delta = np.clip(cdf_on_edges[1:] - cdf_on_edges[:-1], 0.0, 1.0)
    return N_total * delta


def _save_text_table(path, header, cols):
    """
    Save columns to a whitespace-delimited text file with a header.

    Parameters
    ----------
    path : str
        Output path.
    header : str
        Header string (multi-line allowed).
    cols : list[array_like]
        Same-length columns.
    """
    arr = np.column_stack(cols)
    np.savetxt(path, arr, header=header, comments='', fmt="%.8e")


def save_monochromatic_outputs(chosen, out_dir=MONO_RESULTS_DIR):
    """
    Save selected monochromatic spectra to individual text files.

    Parameters
    ----------
    chosen : list[tuple[float, ndarray, ndarray]]
        From `plot_monochromatic`: (mass, E, dNdE) per selection.
    out_dir : str
        Destination directory.

    Returns
    -------
    list[str]
        Written file paths.
    """
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for m, E, dNdE in chosen:
        fname = f"{m:.2e}_spectrum.txt".replace("+", "")
        path = os.path.join(out_dir, fname)
        header = (f"Monochromatic total spectrum\n"
                  f"Mass = {m:.6e} g\n"
                  f"Columns: E[MeV], dN/dE[MeV^-1 s^-1]")
        _save_text_table(path, header, [E, dNdE])
        paths.append(path)
    return paths


# ---------------------------
# Sampling utilities (distributed)
# ---------------------------
def _sample_from_pdf_grid(pdf_m, pdf_vals, N):
    """
    Draw samples from a positive PDF tabulated on a grid via inverse transform.

    Parameters
    ----------
    pdf_m : ndarray
        Ascending support grid (positive).
    pdf_vals : ndarray
        Nonnegative PDF values on that grid.
    N : int
        Number of samples.

    Returns
    -------
    ndarray
        Samples (length N), clipped to [MASS_MIN, MASS_MAX].
    """
    pdf_m = np.asarray(pdf_m, float)
    pdf_vals = np.clip(np.asarray(pdf_vals, float), 0.0, None)
    mask = (pdf_m > 0) & np.isfinite(pdf_vals)
    pdf_m = pdf_m[mask]; pdf_vals = pdf_vals[mask]
    if pdf_m.size < 2 or pdf_vals.sum() <= 0:
        return np.full(N, np.nan)

    # Normalize
    area = trapezoid(pdf_vals, pdf_m)
    if area <= 0:
        return np.full(N, np.nan)
    pdf_vals = pdf_vals / area

    # CDF
    cdf = np.concatenate([[0.0], np.cumsum(trapezoid([pdf_vals[:-1], pdf_vals[1:]],
                                                    x=[pdf_m[:-1], pdf_m[1:]], axis=0))])
    cdf = cdf / max(cdf[-1], 1e-300)
    grid = np.concatenate([[pdf_m[0]], pdf_m[1:]])  # identical pointer for interp

    # Inverse CDF sampling
    u = np.random.random(N)
    samples = np.interp(u, cdf, np.concatenate([[pdf_m[0]], pdf_m[1:]]))
    return np.clip(samples, MASS_MIN, MASS_MAX)


def sample_gaussian_collapse(peaks, sigmas, N_total, *, kappa=1.0, delta_c=0.45, gamma=0.36):
    """
    Sample masses for the *Gaussian collapse* prescription.

    Parameters
    ----------
    peaks : list[float]
        Peak masses [g] around which to sample.
    sigmas : list[float]
        Width parameters σ_x (dimensionless; typical 0.03–0.255).
    N_total : int
        Total number of samples across all (peak, σ) combinations.
    kappa, delta_c, gamma : float
        Collapse model parameters.

    Returns
    -------
    dict
        {
          'masses': ndarray,
          'per_case': list[tuple[float, float, int]]  # (peak, sigma, N_i)
        }
    """
    peaks = list(map(float, peaks))
    sigmas = list(map(float, sigmas))
    cases = [(p, s) for p in peaks for s in sigmas]
    if not cases:
        return {'masses': np.array([]), 'per_case': []}

    N_each = max(1, N_total // len(cases))
    masses_all, per_case = [], []

    r_grid = np.geomspace(0.25, 4.0, 2000)
    for (p, s) in cases:
        d_l = delta_l(r_grid, kappa, delta_c, gamma)
        phi = np.clip(mass_function(d_l, s, delta_c, gamma), 0.0, None)
        # Normalize on r-grid and convert to mass
        m_grid = np.clip(p * r_grid, MASS_MIN, MASS_MAX)
        # The PDF in *mass* inherits |dr/dm| = 1/p on a change of variables
        pdf_m = np.clip(phi / max(p, 1e-300), 0.0, None)

        # Sample
        samp = _sample_from_pdf_grid(m_grid, pdf_m, N_each)
        masses_all.append(samp)
        per_case.append((p, s, len(samp)))

    return {'masses': np.concatenate(masses_all) if masses_all else np.array([]),
            'per_case': per_case}


def sample_non_gaussian_collapse(peaks, sigmaX, sigmaY, N_total, *, delta_c=0.45, gamma=0.36, kappa=1.0):
    """
    Sample masses for the *Non-Gaussian collapse* prescription (Biagetti proxy).

    Parameters
    ----------
    peaks : list[float]
        Peak masses [g].
    sigmaX, sigmaY : list[float]
        Variance parameters (dimensionless).
    N_total : int
        Total number of samples across all combinations.
    delta_c, gamma, kappa : float
        Model parameters.

    Returns
    -------
    dict
        {
          'masses': ndarray,
          'per_case': list[tuple[float, float, float, int]]  # (peak, sigmaX, sigmaY, N_i)
        }
    """
    peaks = list(map(float, peaks))
    sigmaX = list(map(float, sigmaX))
    sigmaY = list(map(float, sigmaY))
    cases = [(p, sx, sy) for p in peaks for sx in sigmaX for sy in sigmaY]
    if not cases:
        return {'masses': np.array([]), 'per_case': []}

    N_each = max(1, N_total // len(cases))
    masses_all, per_case = [], []

    r_grid = np.geomspace(0.25, 4.0, 2200)
    for (p, sx, sy) in cases:
        d_l = delta_l(r_grid, kappa, delta_c, gamma)
        phi = np.clip(mass_function_exact(d_l, sx, sy, delta_c, gamma), 0.0, None)
        m_grid = np.clip(p * r_grid, MASS_MIN, MASS_MAX)
        pdf_m = np.clip(phi / max(p, 1e-300), 0.0, None)
        samp = _sample_from_pdf_grid(m_grid, pdf_m, N_each)
        masses_all.append(samp)
        per_case.append((p, sx, sy, len(samp)))

    return {'masses': np.concatenate(masses_all) if masses_all else np.array([]),
            'per_case': per_case}


def sample_lognormal(peaks, sigmas_ln, N_total):
    """
    Sample masses from a *log-normal* mass function centered on each peak.

    Parameters
    ----------
    peaks : list[float]
        Peak masses [g] that set the central tendency (μ = ln(peak)).
    sigmas_ln : list[float]
        Log-space standard deviations (σ > 0).
    N_total : int
        Total number of samples.

    Returns
    -------
    dict
        {
          'masses': ndarray,
          'per_case': list[tuple[float, float, int]]  # (peak, sigma_ln, N_i)
        }
    """
    peaks = list(map(float, peaks))
    sigmas_ln = list(map(float, sigmas_ln))
    cases = [(p, s) for p in peaks for s in sigmas_ln]
    if not cases:
        return {'masses': np.array([]), 'per_case': []}

    N_each = max(1, N_total // len(cases))
    masses_all, per_case = [], []

    m_grid = np.geomspace(MASS_MIN, MASS_MAX, 4000)
    for (p, s) in cases:
        mu = np.log(p)
        pdf = mass_function_lognormal(m_grid, mu, s)
        samp = _sample_from_pdf_grid(m_grid, pdf, N_each)
        masses_all.append(samp)
        per_case.append((p, s, len(samp)))

    return {'masses': np.concatenate(masses_all) if masses_all else np.array([]),
            'per_case': per_case}


# ---------------------------
# Custom-equation PDF (expr → PDF)
# ---------------------------
_ALLOWED_FUNCS = {
    # numpy rebindings
    'exp': np.exp,
    'log': np.log,
    'log10': np.log10,
    'sqrt': np.sqrt,
    'sin': np.sin,
    'cos': np.cos,
    'tan': np.tan,
    'sinh': np.sinh,
    'cosh': np.cosh,
    'tanh': np.tanh,
    'abs': np.abs,
    'where': np.where,
    'pi': np.pi,
    'e': np.e,
}

_RESERVED = set(list(_ALLOWED_FUNCS.keys()) + ['m'])


def _detect_custom_variables(expr):
    """
    Detect variable names in a user PDF expression f(m) excluding allowed funcs & 'm'.

    Parameters
    ----------
    expr : str
        Raw Python expression using 'm' plus any user constants like 'mp', 'alpha0', etc.

    Returns
    -------
    list[str]
        Sorted unique variable names requiring values.
    """
    tokens = set(re.findall(r'\b[A-Za-z_][A-Za-z0-9_]*\b', expr))
    vars_needed = sorted([t for t in tokens if t not in _RESERVED])
    return vars_needed


def _evaluate_custom_pdf(expr, m_grid, var_map):
    """
    Evaluate the user-provided PDF expression safely on a grid.

    Parameters
    ----------
    expr : str
        Expression like "(m/mp)**(-(alpha0 + beta*log(m/mp)))/m"
    m_grid : ndarray
        Positive mass grid.
    var_map : dict
        Map of variable -> float.

    Returns
    -------
    ndarray
        Nonnegative values; NaNs and negatives are clipped to 0.

    Notes
    -----
    This is not a sandbox, but we provide a restricted eval namespace.
    """
    local_ns = dict(_ALLOWED_FUNCS)
    local_ns.update(var_map)
    local_ns['m'] = m_grid
    try:
        vals = eval(expr, {"__builtins__": {}}, local_ns)
    except Exception as e:
        raise ValueError(f"Error evaluating expression: {e}")
    vals = np.asarray(vals, float)
    vals[~np.isfinite(vals)] = 0.0
    vals = np.clip(vals, 0.0, None)
    return vals


def sample_custom_equation(expr, var_map, N_total):
    """
    Sample from a user-specified PDF f(m) defined over [MASS_MIN, MASS_MAX].

    Parameters
    ----------
    expr : str
        PDF expression using 'm' and any constants in `var_map`.
    var_map : dict[str, float]
        Variable values.
    N_total : int
        Number of samples.

    Returns
    -------
    dict
        {
          'masses': ndarray,
          'grid_m': ndarray,
          'grid_pdf': ndarray  (normalized),
        }
    """
    m_grid = np.geomspace(MASS_MIN, MASS_MAX, 5000)
    raw = _evaluate_custom_pdf(expr, m_grid, var_map)
    area = trapezoid(raw, m_grid)
    if area <= 0:
        return {'masses': np.array([]), 'grid_m': m_grid, 'grid_pdf': np.zeros_like(m_grid)}
    pdf = raw / area
    samples = _sample_from_pdf_grid(m_grid, pdf, N_total)
    return {'masses': samples, 'grid_m': m_grid, 'grid_pdf': pdf}


# ---------------------------
# Distributed plotting (hist + PDF overlay)
# ---------------------------
def plot_mass_histogram_with_pdf(masses, *, pdf_m=None, pdf_vals=None, title="", show_pdf_line=True):
    """
    Plot a log-spaced histogram of sampled masses and optionally overlay an analytic PDF.

    Parameters
    ----------
    masses : array_like
        Sampled masses [g].
    pdf_m : ndarray or None
        Mass grid for the analytic PDF (if provided).
    pdf_vals : ndarray or None
        PDF values (must integrate to 1 over the domain of interest).
    title : str
        Plot title.
    show_pdf_line : bool
        Whether to overlay the expected-counts curve.

    Returns
    -------
    dict
        {
          'edges': ndarray,
          'counts': ndarray,
          'pdf_counts': ndarray or None
        }
    """
    x = np.asarray(masses, float)
    x = x[np.isfinite(x)]
    x = x[(x > 0) & (x >= MASS_MIN) & (x <= MASS_MAX)]
    if x.size == 0:
        warn("No valid samples to plot.")
        return {'edges': np.array([MASS_MIN, MASS_MAX]), 'counts': np.array([0]), 'pdf_counts': None}

    edges = _auto_log_bins(x, max_bins=MAX_HIST_BINS)
    counts, _ = np.histogram(x, bins=edges)
    centers = np.sqrt(edges[1:] * edges[:-1])

    plt.figure()
    # step histogram
    plt.hist(x, bins=edges, histtype='step', log=True, label="Samples")

    pdf_counts = None
    if show_pdf_line and (pdf_m is not None) and (pdf_vals is not None):
        pdf_counts = _scale_pdf_to_hist_counts(edges, pdf_m, pdf_vals, N_total=len(x))
        # Make a step-like overlay using edges
        step_x = np.repeat(edges, 2)[1:-1]
        step_y = np.repeat(pdf_counts, 2)
        plt.plot(step_x, np.clip(step_y, 1e-12, None), linestyle='-', label="PDF (expected counts)")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("PBH mass m [g]")
    plt.ylabel("Counts")
    plt.title(title)
    plt.grid(True, which='both', alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return {'edges': edges, 'counts': counts, 'pdf_counts': pdf_counts}
# ---------------------------
# Formatting & tiny I/O helpers
# ---------------------------
def _banner():
    print("\n" + "═" * 58)
    print(" GammaPBHPlotter: PBH Spectrum Tool".center(58))
    print(f" Version {PKG_VERSION}".center(58))
    print("═" * 58 + "\n")


def _fmt_path(p):
    try:
        return os.path.normpath(p)
    except Exception:
        return p


def _ask(prompt):
    try:
        return input(prompt).strip()
    except (EOFError, KeyboardInterrupt):
        return "q"


def _confirm_yes_no(msg, default="n"):
    s = _ask(f"{msg} (y/n): ").lower()
    if s == "":
        s = default
    return s.startswith("y")


def _ask_int(msg, *, min_val=None, max_val=None, default=None):
    while True:
        raw = _ask(msg)
        if raw.lower() in ("b", "q"):
            return raw
        try:
            v = int(raw)
            if (min_val is not None and v < min_val) or (max_val is not None and v > max_val):
                print(f"Value must be within [{min_val}, {max_val}]")
                continue
            return v
        except ValueError:
            if default is not None and raw == "":
                return default
            print("Please enter an integer.")


def _ask_float_list(msg, *, allow_empty=False):
    while True:
        raw = _ask(msg)
        if raw.lower() in ("b", "q"):
            return raw
        if raw == "" and allow_empty:
            return []
        try:
            vals = [float(s) for s in re.split(r"[,\s]+", raw) if s != ""]
            if not vals and not allow_empty:
                print("Please enter at least one value.")
                continue
            return vals
        except ValueError:
            print("Could not parse floats. Use commas (e.g. 3e15, 8e16).")


def _ask_float(msg, *, min_val=None, max_val=None, default=None):
    while True:
        raw = _ask(msg)
        if raw.lower() in ("b", "q"):
            return raw
        try:
            v = float(raw)
            if (min_val is not None and v < min_val) or (max_val is not None and v > max_val):
                print(f"Value must be within [{min_val:.2e}, {max_val:.2e}]")
                continue
            return v
        except ValueError:
            if default is not None and raw == "":
                return default
            print("Please enter a number.")


def _print_samples_saved(path):
    print(f"Saved → { _fmt_path(path) }")


# ---------------------------
# Gaussian collapse tool
# ---------------------------
def gaussian_tool(spl):
    """
    CLI flow for Gaussian collapse sampling + histogram with PDF overlay.

    Parameters
    ----------
    spl : SimpleNamespace
        Spline pack used for spectrum evaluation when needed (not strictly required here).
    """
    peaks = _ask_float_list(f"Enter peak PBH masses (g) (comma-separated; each must be within [{MASS_MIN:.2e}, {MASS_MAX:.2e}]): ")
    if isinstance(peaks, str):  # 'b' or 'q'
        return peaks
    for p in peaks:
        if not (MASS_MIN <= p <= MASS_MAX):
            print(f"Peak {p:.3e} out of range.")
            return

    N = _ask_int("Enter target N (integer, e.g. 1000): ", min_val=1)
    if isinstance(N, str):
        return N

    sigmas = _ask_float_list("Enter σ list for Gaussian collapse (comma-separated; each must be within [0.03, 0.255]): ")
    if isinstance(sigmas, str):
        return sigmas
    for s in sigmas:
        if not (0.03 <= s <= 0.255):
            print(f"σ={s} out of allowed range [0.03, 0.255].")
            return

    print(f"Sampling ", end="", flush=True)
    for p in peaks:
        for s in sigmas:
            print(f" peak {p:.2e}  [σ={s}]: ", end="", flush=True)
            break
        break

    res = sample_gaussian_collapse(peaks, sigmas, N)
    masses = res['masses']
    # Build analytic PDF for overlay for *one* case at a time? We’ll approximate by mixing equally:
    m_grid = np.geomspace(MASS_MIN, MASS_MAX, 4000)
    pdf_mix = np.zeros_like(m_grid)
    for (p, s, _) in res['per_case']:
        r_grid = np.geomspace(0.25, 4.0, 2000)
        d_l = delta_l(r_grid, 1.0, 0.45, 0.36)
        phi = np.clip(mass_function(d_l, s, 0.45, 0.36), 0.0, None)
        m_case = np.clip(p * r_grid, MASS_MIN, MASS_MAX)
        pdf_case = np.clip(phi / max(p, 1e-300), 0.0, None)
        # interpolate to common m_grid and normalize later
        pdf_mix += np.interp(m_grid, m_case, pdf_case, left=0.0, right=0.0)
    area = trapezoid(pdf_mix, m_grid)
    if area > 0:
        pdf_mix /= area

    title = "Distributed spectra (Gaussian collapse) — histogram"
    plot_mass_histogram_with_pdf(masses, pdf_m=m_grid, pdf_vals=pdf_mix, title=title, show_pdf_line=True)

    # Save?
    if _confirm_yes_no("Save distributed results?"):
        # Save one directory per unique (p,s) to mirror previous UX, but store a single multi-case as well
        os.makedirs(GAUSS_RESULTS_DIR, exist_ok=True)
        written = []
        if len(res['per_case']) == 1:
            p, s, Ncase = res['per_case'][0]
            dname = os.path.join(GAUSS_RESULTS_DIR, f"peak_{p:.2e}_σ{s}_N{N}")
            os.makedirs(dname, exist_ok=True)
            np.save(os.path.join(dname, "masses.npy"), masses)
            meta = {"type": "gaussian", "peaks": [p], "sigmas": [s], "N": int(N)}
            with open(os.path.join(dname, "meta.json"), "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)
            # also store PDF grid for future overlay
            np.savetxt(os.path.join(dname, "pdf.csv"), np.column_stack([m_grid, pdf_mix]),
                       delimiter=",", header="m,pdf", comments='')
            _print_samples_saved(dname)
            written.append(dname)
        else:
            # Multi-case bundle
            dname = os.path.join(GAUSS_RESULTS_DIR, f"bundle_{len(res['per_case'])}_cases_N{N}")
            os.makedirs(dname, exist_ok=True)
            np.save(os.path.join(dname, "masses.npy"), masses)
            meta = {
                "type": "gaussian",
                "cases": [{"peak": float(p), "sigma": float(s), "N_i": int(Ni)} for (p, s, Ni) in res['per_case']],
                "N": int(N)
            }
            with open(os.path.join(dname, "meta.json"), "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)
            np.savetxt(os.path.join(dname, "pdf.csv"), np.column_stack([m_grid, pdf_mix]),
                       delimiter=",", header="m,pdf", comments='')
            _print_samples_saved(dname)
            written.append(dname)
        return written


# ---------------------------
# Non-Gaussian collapse tool
# ---------------------------
def non_gaussian_tool(spl):
    """
    CLI for Non-Gaussian collapse sampling.
    """
    peaks = _ask_float_list(f"Enter peak PBH masses (g) (comma-separated; each must be within [{MASS_MIN:.2e}, {MASS_MAX:.2e}]): ")
    if isinstance(peaks, str):
        return peaks
    for p in peaks:
        if not (MASS_MIN <= p <= MASS_MAX):
            print(f"Peak {p:.3e} out of range.")
            return

    N = _ask_int("Enter target N (integer, e.g. 1000): ", min_val=1)
    if isinstance(N, str):
        return N

    sigmaX = _ask_float_list("Enter σ_X list (comma-separated; > 0): ")
    if isinstance(sigmaX, str):
        return sigmaX
    sigmaY = _ask_float_list("Enter σ_Y list (comma-separated; > 0): ")
    if isinstance(sigmaY, str):
        return sigmaY
    for s in sigmaX + sigmaY:
        if not (s > 0):
            print("All σ must be positive.")
            return

    print(f"Sampling ", end="", flush=True)
    for p in peaks:
        for sx in sigmaX:
            for sy in sigmaY:
                print(f" peak {p:.2e}  [σX={sx}, σY={sy}]: ", end="", flush=True)
                break
            break
        break

    res = sample_non_gaussian_collapse(peaks, sigmaX, sigmaY, N)
    masses = res['masses']

    # Build mixture PDF (approx) for overlay
    m_grid = np.geomspace(MASS_MIN, MASS_MAX, 4500)
    pdf_mix = np.zeros_like(m_grid)
    r_grid = np.geomspace(0.25, 4.0, 2400)
    d_l = delta_l(r_grid, 1.0, 0.45, 0.36)
    for (p, sx, sy, _) in res['per_case']:
        phi = np.clip(mass_function_exact(d_l, sx, sy, 0.45, 0.36), 0.0, None)
        m_case = np.clip(p * r_grid, MASS_MIN, MASS_MAX)
        pdf_case = np.clip(phi / max(p, 1e-300), 0.0, None)
        pdf_mix += np.interp(m_grid, m_case, pdf_case, left=0.0, right=0.0)
    area = trapezoid(pdf_mix, m_grid)
    if area > 0:
        pdf_mix /= area

    title = "Distributed spectra (Non-Gaussian collapse) — histogram"
    plot_mass_histogram_with_pdf(masses, pdf_m=m_grid, pdf_vals=pdf_mix, title=title, show_pdf_line=True)

    if _confirm_yes_no("Save distributed results?"):
        os.makedirs(NONGAUSS_RESULTS_DIR, exist_ok=True)
        if len(res['per_case']) == 1:
            p, sx, sy, Ni = res['per_case'][0]
            dname = os.path.join(NONGAUSS_RESULTS_DIR, f"peak_{p:.2e}_σX{sx}_σY{sy}_N{N}")
            os.makedirs(dname, exist_ok=True)
            np.save(os.path.join(dname, "masses.npy"), masses)
            meta = {"type": "nongaussian", "peaks": [p], "sigmaX": [sx], "sigmaY": [sy], "N": int(N)}
            with open(os.path.join(dname, "meta.json"), "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)
            np.savetxt(os.path.join(dname, "pdf.csv"), np.column_stack([m_grid, pdf_mix]),
                       delimiter=",", header="m,pdf", comments='')
            _print_samples_saved(dname)
        else:
            dname = os.path.join(NONGAUSS_RESULTS_DIR, f"bundle_{len(res['per_case'])}_cases_N{N}")
            os.makedirs(dname, exist_ok=True)
            np.save(os.path.join(dname, "masses.npy"), masses)
            meta = {
                "type": "nongaussian",
                "cases": [{"peak": float(p), "sigmaX": float(sx), "sigmaY": float(sy), "N_i": int(Ni)}
                          for (p, sx, sy, Ni) in res['per_case']],
                "N": int(N)
            }
            with open(os.path.join(dname, "meta.json"), "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)
            np.savetxt(os.path.join(dname, "pdf.csv"), np.column_stack([m_grid, pdf_mix]),
                       delimiter=",", header="m,pdf", comments='')
            _print_samples_saved(dname)


# ---------------------------
# Log-normal tool
# ---------------------------
def lognormal_tool(spl):
    """
    CLI for Log-normal distribution sampling.
    """
    peaks = _ask_float_list(f"Enter peak PBH masses (g) (comma-separated; each must be within [{MASS_MIN:.2e}, {MASS_MAX:.2e}]): ")
    if isinstance(peaks, str):
        return peaks
    for p in peaks:
        if not (MASS_MIN <= p <= MASS_MAX):
            print(f"Peak {p:.3e} out of range.")
            return

    N = _ask_int("Enter target N (integer, e.g. 1000): ", min_val=1)
    if isinstance(N, str):
        return N

    sigmas = _ask_float_list("Enter σ list (log-space std) for Log-Normal (comma-separated; each > 0): ")
    if isinstance(sigmas, str):
        return sigmas
    for s in sigmas:
        if not (s > 0):
            print("All σ must be positive.")
            return

    print(f"Sampling ", end="", flush=True)
    for p in peaks:
        for s in sigmas:
            print(f" peak {p:.2e}  [σ=0.1]: ", end="", flush=True)  # print once like classic UX
            break
        break

    res = sample_lognormal(peaks, sigmas, N)
    masses = res['masses']

    # Mixture analytic PDF
    m_grid = np.geomspace(MASS_MIN, MASS_MAX, 4000)
    pdf_mix = np.zeros_like(m_grid)
    for (p, s, _) in res['per_case']:
        mu = np.log(p)
        pdf = mass_function_lognormal(m_grid, mu, s)
        pdf_mix += pdf
    area = trapezoid(pdf_mix, m_grid)
    if area > 0:
        pdf_mix /= area

    title = "Distributed spectra (Log-Normal) — histogram"
    plot_mass_histogram_with_pdf(masses, pdf_m=m_grid, pdf_vals=pdf_mix, title=title, show_pdf_line=True)

    if _confirm_yes_no("Save distributed results?"):
        os.makedirs(LOGNORM_RESULTS_DIR, exist_ok=True)
        if len(res['per_case']) == 1:
            p, s, Ni = res['per_case'][0]
            dname = os.path.join(LOGNORM_RESULTS_DIR, f"peak_{p:.2e}_σ{s}_N{N}")
            os.makedirs(dname, exist_ok=True)
            np.save(os.path.join(dname, "masses.npy"), masses)
            meta = {"type": "lognormal", "peaks": [p], "sigmas": [s], "N": int(N)}
            with open(os.path.join(dname, "meta.json"), "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)
            np.savetxt(os.path.join(dname, "pdf.csv"), np.column_stack([m_grid, pdf_mix]),
                       delimiter=",", header="m,pdf", comments='')
            _print_samples_saved(dname)
        else:
            dname = os.path.join(LOGNORM_RESULTS_DIR, f"bundle_{len(res['per_case'])}_cases_N{N}")
            os.makedirs(dname, exist_ok=True)
            np.save(os.path.join(dname, "masses.npy"), masses)
            meta = {
                "type": "lognormal",
                "cases": [{"peak": float(p), "sigma": float(s), "N_i": int(Ni)} for (p, s, Ni) in res['per_case']],
                "N": int(N)
            }
            with open(os.path.join(dname, "meta.json"), "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)
            np.savetxt(os.path.join(dname, "pdf.csv"), np.column_stack([m_grid, pdf_mix]),
                       delimiter=",", header="m,pdf", comments='')
            _print_samples_saved(dname)


# ---------------------------
# Custom-equation mass-PDF tool
# ---------------------------
def custom_equation_pdf_tool():
    """
    CLI for user-defined mass PDF f(m) over [MASS_MIN, MASS_MAX].

    Steps:
      1) Prompt for expression in terms of 'm' and constants (e.g. mp, alpha0, beta).
      2) Detect unknown variables and prompt for their numeric values.
      3) Ask for target N, sample, plot histogram with analytic-PDF overlay.
      4) Optionally save run (masses.npy, meta.json, pdf.csv).
    """
    print("\n=== Custom Equation Mass PDF ===")
    print(f"Domain: m in [{MASS_MIN:.2e}, {MASS_MAX:.2e}] g")
    print("Enter a Python expression for your PDF f(m) using 'm' in grams and any constants/variables you define.")
    print("Examples:")
    print("  (m/mp)**(-(alpha0 + beta*log(m/mp))) / m")
    print("  exp(-m/5e17) / m")

    expr = _ask("f(m) = ")
    if expr.lower() in ("b", "q"):
        return expr

    vars_needed = _detect_custom_variables(expr)
    var_map = {}
    if vars_needed:
        print("\nDetected variables:", ", ".join(vars_needed))
        for v in vars_needed:
            val = _ask_float(f"Enter value for {v}: ")
            if isinstance(val, str):
                return val
            var_map[v] = float(val)
    else:
        print("No extra variables detected.")

    N = _ask_int("Enter target N (integer, e.g. 1000): ", min_val=1)
    if isinstance(N, str):
        return N

    # Sample & overlay
    res = sample_custom_equation(expr, var_map, N)
    masses = res['masses']
    m_grid, pdf_vals = res['grid_m'], res['grid_pdf']

    title = "Distributed spectra (Custom PDF) — histogram"
    plot_mass_histogram_with_pdf(masses, pdf_m=m_grid, pdf_vals=pdf_vals, title=title, show_pdf_line=True)

    if _confirm_yes_no("Save distributed results?"):
        os.makedirs(CUSTOM_RESULTS_DIR, exist_ok=True)
        # create stable-ish name from expression
        h = hashlib.sha256(expr.encode("utf-8")).hexdigest()[:10]
        dname = os.path.join(CUSTOM_RESULTS_DIR, f"expr_{h}_N{int(N)}")
        os.makedirs(dname, exist_ok=True)
        np.save(os.path.join(dname, "masses.npy"), masses)
        meta = {"type": "custom", "expr": expr, "vars": var_map, "N": int(N)}
        with open(os.path.join(dname, "meta.json"), "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)
        np.savetxt(os.path.join(dname, "pdf.csv"), np.column_stack([m_grid, pdf_vals]),
                   delimiter=",", header="m,pdf", comments='')
        _print_samples_saved(dname)


# ---------------------------
# View previous spectra (monochromatic + distributed w/ consistent legends)
# ---------------------------
def _list_dirs(path):
    try:
        return [d for d in sorted(os.listdir(path)) if os.path.isdir(os.path.join(path, d))]
    except FileNotFoundError:
        return []


def _list_files(path, pattern=r".*"):
    rx = re.compile(pattern)
    try:
        return [f for f in sorted(os.listdir(path)) if os.path.isfile(os.path.join(path, f)) and rx.match(f)]
    except FileNotFoundError:
        return []


def _view_prev_monochromatic():
    files = _list_files(MONO_RESULTS_DIR, pattern=r".*\.txt$")
    if not files:
        print("No saved monochromatic spectra.")
        return
    print("\nSaved monochromatic spectra:")
    for i, f in enumerate(files, 1):
        print(f" {i:2d}: {f}")
    s = _ask("Select indices to plot (comma-separated) or '0' for ALL: ")
    if s.lower() in ("b", "q"):
        return s
    idxs = []
    if s.strip() == "0":
        idxs = list(range(1, len(files)+1))
    else:
        for tok in re.split(r"[,\s]+", s.strip()):
            if tok.isdigit():
                k = int(tok)
                if 1 <= k <= len(files):
                    idxs.append(k)
    if not idxs:
        print("Nothing selected.")
        return
    plt.figure()
    for k in idxs:
        path = os.path.join(MONO_RESULTS_DIR, files[k-1])
        try:
            data = np.loadtxt(path, comments="#")
            E = data[:, 0]; dNdE = data[:, 1]
            label = files[k-1].replace("_spectrum.txt", "")
            plt.loglog(E, np.clip(dNdE, 1e-300, None), label=label)
        except Exception as e:
            print(f"Failed to plot {files[k-1]}: {e}")
    plt.xlabel("E [MeV]")
    plt.ylabel("dN/dE [MeV$^{-1}$ s$^{-1}$]")
    plt.title("Monochromatic spectra (saved)")
    plt.grid(True, which='both', alpha=0.25)
    plt.legend(title="Saved curves")  # consistent legend title
    plt.tight_layout()
    plt.show()


def _plot_saved_hist_bundle(base_dir, entry_dir):
    """
    Load masses.npy, meta.json, optional pdf.csv and replot with consistent legend keys.
    """
    d = os.path.join(base_dir, entry_dir)
    masses_path = os.path.join(d, "masses.npy")
    meta_path = os.path.join(d, "meta.json")
    pdf_path = os.path.join(d, "pdf.csv")
    if not os.path.exists(masses_path):
        print(f"Missing masses.npy in {entry_dir}")
        return
    masses = np.load(masses_path)
    pdf_m = pdf_vals = None
    if os.path.exists(pdf_path):
        try:
            arr = np.loadtxt(pdf_path, delimiter=",", skiprows=1)
            if arr.ndim == 1 and arr.size >= 2:
                arr = arr.reshape(-1, 2)
            pdf_m, pdf_vals = arr[:, 0], arr[:, 1]
        except Exception:
            pdf_m = pdf_vals = None
    # Build a readable title
    title = f"Saved histogram — {entry_dir}"
    plot_mass_histogram_with_pdf(masses, pdf_m=pdf_m, pdf_vals=pdf_vals, title=title, show_pdf_line=True)


def _view_prev_distributed():
    def _list_cat(name, base):
        dirs = _list_dirs(base)
        if dirs:
            print(f"\n{name}:")
            for i, d in enumerate(dirs, 1):
                print(f" {i:2d}: {d}")
        return dirs

    g_dirs = _list_cat("Gaussian", GAUSS_RESULTS_DIR)
    ng_dirs = _list_cat("Non-Gaussian", NONGAUSS_RESULTS_DIR)
    ln_dirs = _list_cat("Log-Normal", LOGNORM_RESULTS_DIR)
    cu_dirs = _list_cat("Custom Eq.", CUSTOM_RESULTS_DIR)

    if not any([g_dirs, ng_dirs, ln_dirs, cu_dirs]):
        print("No saved distributed results found.")
        return

    print("\nUse a prefix to choose a category: g#, ng#, ln#, c#  (e.g., g1, ln3, c2).")
    choice = _ask("Select one to plot (or 'b' to go back): ")
    if choice.lower() in ("b", "q"):
        return choice

    m = re.match(r"^(g|ng|ln|c)\s*(\d+)$", choice.strip(), flags=re.I)
    if not m:
        print("Invalid selection format.")
        return
    cat, idx = m.group(1).lower(), int(m.group(2))
    if cat == "g":
        if 1 <= idx <= len(g_dirs):
            _plot_saved_hist_bundle(GAUSS_RESULTS_DIR, g_dirs[idx-1])
    elif cat == "ng":
        if 1 <= idx <= len(ng_dirs):
            _plot_saved_hist_bundle(NONGAUSS_RESULTS_DIR, ng_dirs[idx-1])
    elif cat == "ln":
        if 1 <= idx <= len(ln_dirs):
            _plot_saved_hist_bundle(LOGNORM_RESULTS_DIR, ln_dirs[idx-1])
    elif cat == "c":
        if 1 <= idx <= len(cu_dirs):
            _plot_saved_hist_bundle(CUSTOM_RESULTS_DIR, cu_dirs[idx-1])


def view_previous_spectra_menu():
    """
    Unified “View previous spectra” menu.

    Sections:
      1) Monochromatic (saved text spectra)
      2) Distributed (Gaussian / Non-Gaussian / Log-Normal / Custom)
      0) Back

    Legends/keys:
      - For histograms: "Samples" (step hist) and "PDF (expected counts)"
      - For lines/graphs: legend title “Saved curves” (file-derived labels)
    """
    while True:
        print("\n=== View previous spectra ===")
        print("1: Monochromatic (saved)")
        print("2: Distributed (Gaussian / Non-Gauss / Log-Normal / Custom)")
        print("0: Back")
        ch = _ask("Choice: ").lower()
        if ch in ("0", "b"):
            return
        if ch == "1":
            r = _view_prev_monochromatic()
            if r == "q":
                sys.exit(0)
        elif ch == "2":
            r = _view_prev_distributed()
            if r == "q":
                sys.exit(0)
        elif ch == "q":
            sys.exit(0)
        else:
            print("Invalid selection.")
# ---------------------------
# Plotting utilities
# ---------------------------
def _choose_hist_edges_log(
    data: np.ndarray,
    *,
    max_bins: int = 20,
    min_bins: int = 5,
) -> np.ndarray:
    """
    Choose logarithmically-spaced histogram edges for strictly-positive `data`,
    capping the number of bins at `max_bins` while preserving your previous logic
    of adapting to the data's spread.

    Strategy
    --------
    1) Start with numpy's 'fd' suggestion to gauge a reasonable bin density.
    2) Clamp the resulting number of bins into [min_bins, max_bins].
    3) Produce geometric (log-spaced) bins from min..max so widths scale well
       across decades (matches our spectra usage).

    Notes
    -----
    - Enforces positivity and trims to [MASS_MIN, MASS_MAX] to avoid surprises.
    - If data is too narrow, fall back to `min_bins`.
    """
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    data = data[(data > 0)]
    if data.size == 0:
        # fallback – degenerate, but won't crash
        return np.geomspace(max(MASS_MIN, 1.0), max(MASS_MIN * 10.0, 10.0), min_bins + 1)

    lo = max(float(np.min(data)), MASS_MIN)
    hi = min(float(np.max(data)), MASS_MAX)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo <= 0 or hi <= lo:
        return np.geomspace(max(MASS_MIN, 1.0), max(MASS_MIN * 10.0, 10.0), min_bins + 1)

    # Probe with numpy's FD edges to estimate density
    try:
        fd_edges = np.histogram_bin_edges(data, bins="fd")
        fd_bins = max(1, len(fd_edges) - 1)
    except Exception:
        fd_bins = min_bins

    nbins = max(min_bins, min(max_bins, fd_bins))
    # Geometric spacing for mass domains
    return np.geomspace(lo, hi, nbins + 1)


def _expected_counts_from_pdf(
    edges: np.ndarray,
    pdf_m: np.ndarray,
    pdf_vals: np.ndarray,
    n_samples: int,
) -> np.ndarray:
    """
    Convert a (m, pdf(m)) curve (with units 1/g) into expected bin counts
    for the provided `edges`, using numerical integration per bin.

    Parameters
    ----------
    edges : ndarray, shape (B+1,)
        Histogram edges (in grams).
    pdf_m : ndarray
        Monotonic increasing mass grid for PDF.
    pdf_vals : ndarray
        PDF values sampled at pdf_m, non-negative, not necessarily normalized.
    n_samples : int
        Total number of sampled masses. Expected counts sum to ~n_samples.

    Returns
    -------
    ndarray, shape (B,)
        Expected counts per bin (can be fractional).
    """
    edges = np.asarray(edges, float)
    pdf_m = np.asarray(pdf_m, float)
    pdf_vals = np.asarray(pdf_vals, float)

    # Ensure normalization of PDF first (trapezoid)
    area = trapezoid(np.clip(pdf_vals, 0.0, None), pdf_m)
    if area <= 0 or not np.isfinite(area):
        return np.zeros(len(edges) - 1, dtype=float)
    pdf_vals = np.clip(pdf_vals / area, 0.0, None)

    # Prepare high-res integrand by interpolation to avoid aliasing in thin bins
    # We will integrate bin-by-bin with a small local grid.
    counts = np.zeros(len(edges) - 1, dtype=float)
    # Build a single interpolator once
    # (we stay with numpy interp; robust and fast for 1D)
    for i in range(len(edges) - 1):
        a, b = float(edges[i]), float(edges[i + 1])
        if not (np.isfinite(a) and np.isfinite(b) and b > a and b > 0):
            continue
        # Create a small dense grid per bin (log spacing respects mass scale)
        # 64 points per bin is plenty for smooth PDFs
        m_local = np.geomspace(max(a, pdf_m[0]), min(b, pdf_m[-1]), 64) if b > a and b > 0 else np.array([a, b])
        p_local = np.interp(m_local, pdf_m, pdf_vals, left=0.0, right=0.0)
        # Probability mass in this bin:
        prob_bin = trapezoid(p_local, m_local)
        counts[i] = float(n_samples) * prob_bin

    return counts


def plot_mass_histogram_with_pdf(
    masses: np.ndarray,
    *,
    pdf_m: np.ndarray | None = None,
    pdf_vals: np.ndarray | None = None,
    title: str = "",
    show_pdf_line: bool = True,
) -> None:
    """
    Render a mass histogram (log-x) with optional analytic PDF overlay as
    'expected counts' bars/line. Caps max bins at 20 per user request.

    Parameters
    ----------
    masses : ndarray
        Sampled PBH masses [g], strictly positive.
    pdf_m, pdf_vals : ndarray or None
        If both are provided, the PDF curve is normalized and converted to
        expected counts per bin for an apples-to-apples visual with histogram.
    title : str
        Figure title.
    show_pdf_line : bool
        If True and (pdf_m, pdf_vals) exist, draw the PDF expected counts curve.

    Notes
    -----
    - Always uses log-spaced binning, capped at 20 bins.
    - Legend keys are standardized as:
        * "Samples" for the step histogram
        * "PDF (expected counts)" for the overlay
    """
    masses = np.asarray(masses, float)
    masses = masses[np.isfinite(masses)]
    masses = masses[masses > 0]
    if masses.size == 0:
        print("Nothing to plot: empty mass sample.")
        return

    # Choose edges with cap=20
    edges = _choose_hist_edges_log(masses, max_bins=20, min_bins=5)

    # Histogram (counts)
    fig = plt.figure()
    counts, edges_plot, _ = plt.hist(
        masses,
        bins=edges,
        histtype="step",
        linewidth=1.6,
        label="Samples",
    )

    # Optional: overlay expected counts from PDF
    if show_pdf_line and (pdf_m is not None) and (pdf_vals is not None):
        try:
            exp_counts = _expected_counts_from_pdf(edges, pdf_m, pdf_vals, n_samples=masses.size)
            # To show as a continuous curve, place each bin's expected count at bin center
            centers = np.sqrt(edges[:-1] * edges[1:])  # geometric centers suit log-x
            plt.plot(centers, np.clip(exp_counts, 0.0, None), label="PDF (expected counts)")
        except Exception as e:
            print(f"Warning: failed to overlay PDF: {e}")

    plt.xscale("log")
    plt.xlabel("Mass m [g]")
    plt.ylabel("Counts")
    plt.title(title if title else "Mass histogram")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(title="Histogram / PDF")
    plt.tight_layout()
    plt.show()


# ---------------------------
# Custom-equation helpers
# ---------------------------
_ALLOWED_FUNCS = {
    # numpy aliases
    "exp": np.exp,
    "log": np.log,
    "log10": np.log10,
    "sqrt": np.sqrt,
    "abs": np.abs,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "arcsin": np.arcsin,
    "arccos": np.arccos,
    "arctan": np.arctan,
    "sinh": np.sinh,
    "cosh": np.cosh,
    "tanh": np.tanh,
    "power": np.power,
    "where": np.where,
    "minimum": np.minimum,
    "maximum": np.maximum,
    # constants
    "pi": np.pi,
    "e": np.e,
}


_IDENTIFIER_RX = re.compile(r"\b[A-Za-z_]\w*\b")


def _detect_custom_variables(expr: str) -> list[str]:
    """
    Detect symbol names in `expr` that are not built-ins, not 'm',
    and not in our allowed func/const map. These are treated as variables
    to be provided by the user (e.g., 'mp', 'alpha0', 'beta').

    Parameters
    ----------
    expr : str
        Python expression string for f(m).

    Returns
    -------
    list[str]
        Sorted list of variable identifiers.
    """
    tokens = set(_IDENTIFIER_RX.findall(expr))
    # Exclude Python keywords (just in case)
    import keyword
    tokens = {t for t in tokens if not keyword.iskeyword(t)}

    # Allow 'm' and all known funcs/consts
    allowed = set(_ALLOWED_FUNCS.keys()) | {"m"}
    vars_out = sorted([t for t in tokens if t not in allowed])
    return vars_out


def _eval_custom_pdf(expr: str, m_grid: np.ndarray, var_map: dict[str, float]) -> np.ndarray:
    """
    Evaluate the user PDF expression safely on a mass grid.

    Parameters
    ----------
    expr : str
        Expression with 'm' for mass and any variables in `var_map`.
    m_grid : ndarray
        Mass grid [g] (strictly positive).
    var_map : dict
        Variable name → float value.

    Returns
    -------
    ndarray
        Non-negative (clipped) PDF values sampled on `m_grid`.
    """
    # Build evaluation namespace; disable builtins
    safe_globals = {"__builtins__": None}
    safe_globals.update(_ALLOWED_FUNCS)

    # Locals: vectorized 'm' and scalar variables
    safe_locals = {"m": m_grid}
    for k, v in (var_map or {}).items():
        try:
            safe_locals[k] = float(v)
        except Exception:
            raise ValueError(f"Variable '{k}' must be a number.")

    # Evaluate
    try:
        vals = eval(expr, safe_globals, safe_locals)  # noqa: S307 (controlled env)
    except Exception as e:
        raise ValueError(f"Failed to evaluate expression: {e}")

    vals = np.asarray(vals, float)
    if vals.shape != m_grid.shape:
        # Allow broadcasting but end with matching shape
        try:
            vals = np.broadcast_to(vals, m_grid.shape)
        except Exception:
            raise ValueError("Expression did not produce a vector compatible with mass grid.")
    # Ensure non-negative (PDF)
    return np.clip(vals, 0.0, None)


def sample_custom_equation(
    expr: str,
    var_map: dict[str, float],
    N: int,
) -> dict:
    """
    Sample masses from a user-provided PDF f(m) over [MASS_MIN, MASS_MAX].

    Parameters
    ----------
    expr : str
        Expression in terms of 'm' and variables in `var_map`.
    var_map : dict[str, float]
        Name→value map for variables detected in `expr`.
    N : int
        Number of samples to draw.

    Returns
    -------
    dict
        {
          "masses": ndarray,      # sampled masses
          "grid_m": ndarray,      # mass grid
          "grid_pdf": ndarray,    # normalized PDF on grid
        }
    """
    # Work on a dense geometric grid to respect mass decades (robust integration)
    grid_m = np.geomspace(MASS_MIN, MASS_MAX, 50_000)
    raw_pdf = _eval_custom_pdf(expr, grid_m, var_map)

    # Normalize
    area = trapezoid(raw_pdf, grid_m)
    if area <= 0 or not np.isfinite(area):
        raise ValueError("Custom PDF integrates to 0 or NaN/Inf over the domain.")
    pdf = raw_pdf / area

    # Build CDF for inverse transform
    cdf = np.cumsum((pdf[:-1] + pdf[1:]) * np.diff(grid_m) * 0.5)
    cdf = np.concatenate([[0.0], cdf])
    # Avoid numerical overshoot
    cdf /= cdf[-1] if cdf[-1] > 0 else 1.0

    # Inverse-CDF sampling
    u = np.random.random(int(N))
    masses = np.interp(u, cdf, grid_m)

    # Print a compact summary to echo variable context (helps debugging)
    if var_map:
        summary_vars = ", ".join(f"{k}={v:g}" for k, v in var_map.items())
        print(f"Custom PDF variables: {summary_vars}")

    return {"masses": masses, "grid_m": grid_m, "grid_pdf": pdf}
