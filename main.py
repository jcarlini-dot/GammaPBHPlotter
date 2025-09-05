#!/usr/bin/env python3
import importlib
import subprocess
import sys
import os

# ---------------------------
# Auto-install
# ---------------------------
def _ensure_packages(modules):
    for import_name, pkg_name in modules.items():
        try:
            importlib.import_module(import_name)
        except ImportError:
            print(f"⏳ Installing missing package: {pkg_name} …")
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg_name, "--quiet"])

_REQUIRED = {
    "colorama":   "colorama",
    "numpy":      "numpy",
    "matplotlib": "matplotlib",
    "tqdm":       "tqdm",
    "scipy":      "scipy",
}
_ensure_packages(_REQUIRED)
del _ensure_packages, _REQUIRED

# ---------------------------
# Imports
# ---------------------------
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.special import erf
from scipy.interpolate import RectBivariateSpline
from colorama import Fore, Style

plt.rcParams.update({'font.size': 12})
np.seterr(divide='ignore', invalid='ignore')  # suppress log warnings

# ---------------------------
# Paths
# ---------------------------
BASE_DIR         = os.path.dirname(os.path.abspath(__file__))
DATA_DIR         = os.path.join(BASE_DIR, "blackhawk_data")
RESULTS_DIR      = os.path.join(BASE_DIR, "results")
MONO_RESULTS_DIR = os.path.join(RESULTS_DIR, "monochromatic")
os.makedirs(MONO_RESULTS_DIR, exist_ok=True)

# ---------------------------
# Labels
# ---------------------------
GAUSSIAN_METHOD     = "Gaussian collapse"  # collapse model (σ≤0.255 enforced)
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
    """Raised when the user enters 'b' or 'back' to return to the prior screen."""
    pass

def discover_mass_folders(data_dir):
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
# CLI + parsing helpers (no defaults; verbose skipping)
# ---------------------------
def info(msg):  print(Fore.CYAN   + "ℹ " + msg + Style.RESET_ALL)
def warn(msg):  print(Fore.YELLOW + "⚠ " + msg + Style.RESET_ALL)
def err(msg):   print(Fore.RED    + "✖ " + msg + Style.RESET_ALL)

def user_input(prompt, *, allow_back=False, allow_exit=True):
    """
    Wrapper for input() that optionally allows 'b'/'back' to raise BackRequested.
    'exit' or 'q' will terminate the program if allow_exit is True.
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
    try:
        return sorted(d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)))
    except FileNotFoundError:
        return []

def snap_to_available(mval, available, tol=1e-12):
    log_m = np.log(mval)
    log_available = np.log(np.array(available))
    diffs = np.abs(log_available - log_m)
    idx = np.argmin(diffs)
    return available[idx] if diffs[idx] < tol else None

def parse_float_list_verbose(s, *, name="value", bounds=None, allow_empty=False,
                             positive_only=False, strict_gt=False, strict_lt=False):
    if (s is None or s.strip() == ""):
        if not allow_empty:
            warn(f"No {name}s provided.")
        return []
    vals, seen = [], set()
    lo, hi = (bounds or (None, None))
    for tok in s.split(","):
        t = tok.strip()
        if not t: continue
        try:
            v = float(t)
        except Exception:
            warn(f"Skipping token '{t}': {name} is not a valid number."); continue
        if positive_only and v <= 0:
            warn(f"Skipping {name} {v:g}: must be > 0."); continue
        if lo is not None:
            if (strict_gt and not (v > lo)) or (not strict_gt and not (v >= lo)):
                cmp = ">" if strict_gt else "≥"
                warn(f"Skipping {name} {v:g}: must be {cmp} {lo:g}."); continue
        if hi is not None:
            if (strict_lt and not (v < hi)) or (not strict_lt and not (v <= hi)):
                cmp = "<" if strict_lt else "≤"
                warn(f"Skipping {name} {v:g}: must be {cmp} {hi:g}."); continue
        if v in seen:
            warn(f"Duplicate {name} {v:g}: keeping first, skipping this one."); continue
        vals.append(v); seen.add(v)
    if not vals and not allow_empty:
        warn(f"No usable {name}s parsed.")
    return vals

# ---------------------------
# PDFs (collapse space)
# ---------------------------
def delta_l(mass_ratio, kappa, delta_c, gamma):
    y = (mass_ratio / kappa)**(1.0 / gamma)
    arg = 64 - 96 * (delta_c + y)
    arg = np.clip(arg, 0.0, None)
    return (8 - np.sqrt(arg)) / 6

def mass_function(delta_l_val, sigma_x, delta_c, gamma):
    term1 = 1.0 / (np.sqrt(2 * np.pi) * sigma_x)
    term2 = np.exp(-delta_l_val**2 / (2 * sigma_x**2))
    term3 = delta_l_val - (3/8) * delta_l_val**2 - delta_c
    term4 = gamma * np.abs(1 - (3/4) * delta_l_val)
    return term1 * term2 * term3 / term4

def mass_function_exact(delta_l_val, sigma_X, sigma_Y, delta_c, gamma):
    # Biagetti et al. Eq. (20)
    A = sigma_X**2 + (sigma_Y * delta_l_val)**2
    exp_pref = np.exp(-1.0 / (2.0 * sigma_Y**2))
    term1 = 2.0 * sigma_Y * np.sqrt(A)
    inner_exp = np.exp(sigma_X**2 / (2.0 * sigma_Y**2 * (sigma_X**2 + 2.0 * (sigma_Y * delta_l_val)**2)))
    erf_arg = sigma_X * np.sqrt(2.0 * sigma_Y) / np.sqrt(A)
    term2 = np.sqrt(2.0 * np.pi) * sigma_X * inner_exp * erf(erf_arg)
    bracket = term1 + term2
    norm = exp_pref * sigma_X / (2.0 * np.pi * A**1.5)
    jacobian = ((delta_l_val - 0.375 * delta_l_val**2 - delta_c) /
                (gamma * np.abs(1.0 - 0.75 * delta_l_val)))
    return norm * bracket * jacobian

def mass_function_lognormal(x, mu, sigma):
    x_clipped = np.clip(x, 1e-16, None)
    return (1.0 / (x_clipped * sigma * np.sqrt(2 * np.pi))
            * np.exp(- (np.log(x_clipped) - mu)**2 / (2 * sigma**2)))

# ---------------------------
# Data loaders
# ---------------------------
def load_data(filepath, skip_header=0):
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    return np.genfromtxt(filepath, skip_header=skip_header)

def load_spectra_components(directory):
    primary   = load_data(os.path.join(directory, "instantaneous_primary_spectra.txt"), skip_header=2)[123:]
    secondary = load_data(os.path.join(directory, "instantaneous_secondary_spectra.txt"), skip_header=1)
    IFA_prim  = load_data(os.path.join(directory, "inflight_annihilation_prim.txt"))
    IFA_sec   = load_data(os.path.join(directory, "inflight_annihilation_sec.txt"))
    FSR_prim  = load_data(os.path.join(directory, "final_state_radiation_prim.txt"), skip_header=1)
    FSR_sec   = load_data(os.path.join(directory, "final_state_radiation_sec.txt"), skip_header=1)

    E_prim = primary[:,0] * 1e3  # MeV
    E_sec  = secondary[:,0]      # MeV

    return {
        'energy_primary':         E_prim,
        'energy_secondary':       E_sec,
        'direct_gamma_primary':   primary[:,1] / 1e3,
        'direct_gamma_secondary': secondary[:,1],
        'IFA_primary':            np.interp(E_prim, IFA_prim[:,0], IFA_prim[:,1], left=0.0, right=0.0),
        'IFA_secondary':          np.interp(E_sec,   IFA_sec[:,0],  IFA_sec[:,1],  left=0.0, right=0.0),
        'FSR_primary':            np.interp(E_prim, FSR_prim[:,0], FSR_prim[:,1]),
        'FSR_secondary':          np.interp(E_sec,   FSR_sec[:,0],  FSR_sec[:,1])
    }

# ---------------------------
# Monochromatic
# ---------------------------
def monochromatic_spectra():
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
            if not t: continue
            try:
                mval = float(t)
            except Exception:
                warn(f"Skipping mass token '{t}': not a number."); continue
            if not (MIN_MASS <= mval <= MAX_MASS):
                warn(f"Skipping mass {mval:.3e} g: outside allowed range [{MIN_MASS:.2e}, {MAX_MASS:.2e}]."); continue
            mass_list.append(mval)
    if not mass_list:
        warn("No valid masses provided. Returning to menu."); return

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
        p = np.genfromtxt(os.path.join(sub, "inflight_annihilation_prim.txt"))
        s = np.genfromtxt(os.path.join(sub, "inflight_annihilation_sec.txt"))
        Emax_ifa[i] = max(p[:,0].max() if p.size else 0, s[:,0].max() if s.size else 0)

    logM_all = np.log(masses)
    logE     = np.log(E_ref)
    tiny     = 1e-300

    ld = np.log(np.where(direct_mat>tiny, direct_mat, tiny))
    ls = np.log(np.where(secondary_mat>tiny, secondary_mat, tiny))
    li = np.log(np.where(inflight_mat>tiny, inflight_mat, tiny))
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
            for k in range(len(it)-1,0,-1):
                if np.isclose(it[k], it[k-1], rtol=1e-8): it[k] = 0.0
                else: break
            log10i = np.log10(np.where(it>0, it, tiny))
            for j in range(1, len(log10i)):
                if log10i[j] - log10i[j-1] < -50:
                    it[j:] = 0.0; break
            it[E_ref >= Ecut] = 0.0

        tot = d + s + it + f
        tol = 1e-299
        for arr in (d, s, it, f, tot): arr[arr < tol] = 0.0

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
        plt.show(block=True)
        plt.close()

        all_data.append({'mass': mval,'kind': kind,'E': E_ref.copy(),
                         'direct': d.copy(),'secondary': s.copy(),
                         'inflight': it.copy(),'finalstate': f.copy(),
                         'total': tot.copy()})

    if all_data:
        fig = plt.figure(figsize=(10,7))
        summed = np.zeros_like(all_data[0]['E']); peaks = []
        for entry in all_data:
            Ecur = entry['E']; tot = entry['total']; valid = tot>0
            if np.any(valid):
                plt.plot(Ecur[valid], Ecur[valid]**2 * tot[valid], lw=2,
                         label=f"{entry['mass']:.2e} g ({entry['kind']})")
                summed += tot
                peaks.append((Ecur[valid]**2 * tot[valid]).max())
        vs = summed>0
        plt.plot(all_data[0]['E'][vs], all_data[0]['E'][vs]**2 * summed[vs],
                 'k:', lw=3, label="Summed")
        ymax_o = max(peaks) * 1e1; ymin_o = ymax_o / 1e3
        plt.xlabel(r'$E_\gamma$ (MeV)'); plt.ylabel(r'$E^2 dN_\gamma/dE_\gamma$ (MeV s$^{-1}$)')
        plt.xscale('log'); plt.yscale('log'); plt.xlim(0.5, 5000.0); plt.ylim(ymin_o, ymax_o)
        plt.grid(True, which='both', linestyle='--'); plt.legend()
        plt.title('Total Hawking Radiation Spectra (E²·dN/dE)'); plt.tight_layout()
        plt.show(block=True); plt.close(fig)

        sv = user_input("Save any spectra? (y/n): ").strip().lower()
        if sv in ['y', 'yes']:
            print("Select spectra by index to save (single file each):")
            for idx, e in enumerate(all_data, start=1):
                print(f" {idx}: {e['mass']:.2e} g ({e['kind']})")
            choice = user_input("Enter comma-separated indices (e.g. 1,3,5) or '0' to save ALL: ").strip().lower()
            if choice == '0':
                picks = list(range(1, len(all_data)+1))
            else:
                try:
                    picks = [int(x) for x in choice.split(',')]
                except ValueError:
                    err("Invalid indices; skipping save."); picks = []
            for i in picks:
                if 1 <= i <= len(all_data):
                    e = all_data[i - 1]
                    mass_label = f"{e['mass']:.2e}"
                    filename = os.path.join(MONO_RESULTS_DIR, f"{mass_label}_spectrum.txt")
                    data_cols = np.column_stack((e['E'], e['direct'], e['secondary'], e['inflight'], e['finalstate'], e['total']))
                    header = "E_gamma(MeV)    Direct    Secondary    Inflight    FinalState    Total (MeV^-1 s^-1)"
                    np.savetxt(filename, data_cols, header=header, fmt="%e")
                    print(f"Saved → {filename}")

# ---------------------------
# NEW: Robust right-edge spike/cliff trimming
# ---------------------------
def _trim_right_spike(x_line, y_line, up_thresh=1.35, down_thresh=0.35, max_trim_frac=0.10):
    """
    Trim only a *terminal anomaly* at the far right of a curve:
      - Up-spike  : final point jumps way UP ( y[-1] / y[-2] > up_thresh )
      - Down-cliff: final point drops way DOWN ( y[-1] / y[-2] < down_thresh )
    Walk left while that anomaly persists, but trim at most max_trim_frac of points.
    Return the last index to keep (inclusive). If no anomaly, keep all points.
    """
    y = np.asarray(y_line, dtype=float)
    n = y.size
    if n < 3:
        return n - 1

    y_nm1, y_nm2 = y[-1], y[-2]
    if not (np.isfinite(y_nm1) and np.isfinite(y_nm2)) or y_nm2 == 0:
        return n - 1

    ratio = y_nm1 / max(y_nm2, 1e-300)
    # No anomaly → keep all
    if (ratio <= up_thresh) and (ratio >= down_thresh):
        return n - 1

    max_trim = max(3, int(max_trim_frac * n))
    j = n - 1
    trimmed = 0

    if ratio > up_thresh:
        # Upward spike at the very end
        while (
            j > 1 and trimmed < max_trim and
            np.isfinite(y[j]) and np.isfinite(y[j-1]) and
            (y[j] / max(y[j-1], 1e-300) > up_thresh)
        ):
            j -= 1
            trimmed += 1
        return max(j, 2)

    # Downward cliff at the very end
    while (
        j > 1 and trimmed < max_trim and
        np.isfinite(y[j]) and np.isfinite(y[j-1]) and
        (y[j] / max(y[j-1], 1e-300) < down_thresh)
    ):
        j -= 1
        trimmed += 1
    return max(j, 2)

# ---------------------------
# Distributed (Gaussian collapse / Non-Gaussian / Lognormal)
# ---------------------------
def distributed_spectrum(distribution_method):
    """
    LOGNORMAL_METHOD: unbounded mass sampling (mode-centered).
    NON_GAUSSIAN: enforce 0.04 ≤ σ_X ≤ 0.16 (with σ_Y/σ_X = 0.75).
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
            allow_back=True
        )
    except BackRequested:
        return

    peaks = parse_float_list_verbose(pstr, name="peak mass (g)", bounds=(MIN_MASS, MAX_MASS), allow_empty=False)
    if not peaks:
        warn("No valid peaks; returning."); return

    try:
        nstr = user_input("Enter target N (integer, e.g. 1000): ", allow_back=True)
    except BackRequested:        
        return

    try:
        N_target = int(nstr)
        if N_target <= 0: err("N must be > 0. Returning."); return
    except Exception:
        err("Invalid N (not an integer). Returning."); return

    kappa, gamma_p, delta_c = 3.3, 0.36, 0.59
    x_collapse = np.linspace(0.001, 1.30909, 2000)

    param_sets = []
    if is_g:
        try:
            sstr = user_input("Enter σ list for Gaussian collapse (comma-separated; each must be within [0.03, 0.255]): ",
                              allow_back=True).strip()
        except BackRequested:
            return
        sigmas = parse_float_list_verbose(sstr, name="σ", bounds=(0.03, 0.255), allow_empty=False)
        if not sigmas: warn("No valid σ for Gaussian; returning."); return
        for sx in sigmas: param_sets.append({"sigma_x": sx})
    elif is_ng:
        try:
            sx_str = user_input("Enter σ_X list for Non-Gaussian collapse (comma-separated; σ must be within [0.04, 0.16]): ",
                                allow_back=True).strip()
        except BackRequested:
            return
        sigmas_X = parse_float_list_verbose(sx_str, name="σ_X", bounds=(0.04, 0.16), allow_empty=False)
        if not sigmas_X: warn("No valid σ for Non-Gaussian; returning."); return
        for sX in sigmas_X: param_sets.append({"sigma_X": sX, "ratio": 0.75})
    else:
        try:
            sig_str = user_input("Enter σ list (log-space std) for Log-Normal (comma-separated; each > 0): ",
                                 allow_back=True).strip()
        except BackRequested:
            return
        sigmas_ln = parse_float_list_verbose(sig_str, name="σ", bounds=(1e-12, None), allow_empty=False)
        if not sigmas_ln: warn("No valid σ for Log-Normal; returning."); return
        for sln in sigmas_ln: param_sets.append({"sigma_ln": sln})

    first = load_spectra_components(os.path.join(DATA_DIR, names[0]))
    E_grid = first['energy_primary']; logE = np.log(E_grid)
    N_M = len(masses)
    direct_mat     = np.zeros((N_M,len(E_grid)))
    secondary_mat  = np.zeros_like(direct_mat)
    inflight_mat   = np.zeros_like(direct_mat)
    final_mat      = np.zeros_like(direct_mat)
    Emax_ifa       = np.zeros(N_M)

    for i,m in enumerate(masses):
        sub = os.path.join(DATA_DIR, names[i])
        S = load_spectra_components(sub)
        direct_mat[i]    = S['direct_gamma_primary']
        secondary_mat[i] = np.interp(E_grid, S['energy_secondary'], S['direct_gamma_secondary'])
        inflight_mat[i]  = S['IFA_primary'] + np.interp(E_grid, S['energy_secondary'], S['IFA_secondary'])
        final_mat[i]     = S['FSR_primary'] + np.interp(E_grid, S['energy_secondary'], S['FSR_secondary'])
        p = np.genfromtxt(os.path.join(sub, "inflight_annihilation_prim.txt"))
        s = np.genfromtxt(os.path.join(sub, "inflight_annihilation_sec.txt"))
        Emax_ifa[i] = max(p[:,0].max() if p.size else 0, s[:,0].max() if s.size else 0)

    logM_all = np.log(masses)
    floor = 1e-300
    ld = np.log(np.where(direct_mat>floor, direct_mat, floor))
    ls = np.log(np.where(secondary_mat>floor, secondary_mat, floor))
    li = np.log(np.where(inflight_mat>floor, inflight_mat, floor))
    lf = np.log(np.where(final_mat>floor, final_mat, floor))
    sp_d = RectBivariateSpline(logM_all, logE, ld, kx=1, ky=3, s=0)
    sp_s = RectBivariateSpline(logM_all, logE, ls, kx=1, ky=3, s=0)
    sp_i = RectBivariateSpline(logM_all, logE, li, kx=1, ky=3, s=0)
    sp_f = RectBivariateSpline(logM_all, logE, lf, kx=1, ky=3, s=0)

    results = []

    for params in param_sets:

        if is_g:
            sigma_x = params["sigma_x"]
            x = x_collapse
            mf = mass_function(delta_l(x, kappa, delta_c, gamma_p), sigma_x, delta_c, gamma_p)
            label_param = f"σ={sigma_x:.3g}"
            mf = np.where(np.isfinite(mf) & (mf > 0), mf, 0.0)
            if mf.sum() <= 0: warn(f"Underlying PDF vanished for σ={sigma_x:g}; skipping."); continue
            probabilities = mf / mf.sum()
            r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]

        elif is_ng:
            sigma_X = params["sigma_X"]; ratio = params["ratio"]; sigma_Y = ratio * sigma_X
            x = x_collapse
            mf = mass_function_exact(delta_l(x, kappa, delta_c, gamma_p), sigma_X, sigma_Y, delta_c, gamma_p)
            label_param = f"σX={sigma_X:.3g}"
            mf = np.where(np.isfinite(mf) & (mf > 0), mf, 0.0)
            if mf.sum() <= 0: warn(f"Underlying PDF vanished for σ_X={sigma_X:g}; skipping."); continue
            probabilities = mf / mf.sum()
            r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]

        else:
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
                    err(f"Sampling error (lognormal, peak {peak:.3e}, σ={sigma_ln:g}): {e}. Skipping."); bar.close(); continue
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
                        for j in range(len(i_vals)-1,0,-1):
                            if np.isclose(i_vals[j], i_vals[j-1], rtol=1e-8): i_vals[j] = 0.0
                            else: break
                        log10i = np.log10(np.where(i_vals>0, i_vals, floor))
                        for j in range(1,len(log10i)):
                            if log10i[j] - log10i[j-1] < -50: i_vals[j:] = 0.0; break
                        i_vals[E_grid>=Ecut] = 0.0
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
                        for j in range(len(i_vals)-1,0,-1):
                            if np.isclose(i_vals[j], i_vals[j-1], rtol=1e-8): i_vals[j] = 0.0
                            else: break
                        log10i = np.log10(np.where(i_vals>0, i_vals, floor))
                        for j in range(1,len(log10i)):
                            if log10i[j] - log10i[j-1] < -50: i_vals[j:] = 0.0; break
                        i_vals[E_grid>=Ecut] = 0.0
                    sum_d += d_vals; sum_s += s_vals; sum_i += i_vals; sum_f += f_vals
                    bar.update(1)

            bar.close()

            avg_d = sum_d / N_target; avg_s = sum_s / N_target
            avg_i = sum_i / N_target; avg_f = sum_f / N_target
            avg_tot = avg_d + avg_s + avg_i + avg_f
            tol = 1e-299
            for arr in (avg_d, avg_s, avg_i, avg_f, avg_tot): arr[arr < tol] = 0.0

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

    if results:
        # dN/dE overlays
        fig1 = plt.figure(figsize=(10,7))
        peaks_dn = []
        for r in results:
            E = r["E"]; sp = r["spectrum"]; m = sp > 0
            plt.plot(E[m], sp[m], lw=2,
                     label=f"{distribution_method} {r['peak']:.1e}_{r['label_param'].replace('σ=','').replace('σX=','')}")
            peaks_dn.append(sp.max())
        plt.xscale('log'); plt.yscale('log'); plt.xlabel(r'$E_\gamma$ (MeV)'); plt.ylabel(r'$dN_\gamma/dE_\gamma$')
        if peaks_dn: plt.ylim(min(peaks_dn)/1e3, max(peaks_dn)*10)
        plt.xlim(0.5, 5e3); plt.grid(True, which='both', linestyle='--'); plt.legend()
        plt.title("Comparison: dN/dE"); plt.tight_layout(); plt.show(block=True); plt.close(fig1)

        # E^2 dN/dE overlays
        fig2 = plt.figure(figsize=(10,7))
        peaks_e2 = []
        for r in results:
            E = r["E"]; sp = r["spectrum"]; m = sp > 0
            plt.plot(E[m], E[m]**2 * sp[m], lw=2,
                     label=f"{distribution_method} {r['peak']:.1e}_{r['label_param'].replace('σ=','').replace('σX=','')}")
            peaks_e2.append((E[m]**2 * sp[m]).max() if np.any(m) else 0.0)
        plt.xscale('log'); plt.yscale('log'); plt.xlabel(r'$E_\gamma$ (MeV)'); plt.ylabel(r'$E^2\,dN_\gamma/dE_\gamma$')
        if peaks_e2: plt.ylim(min(peaks_e2)/1e3, max(peaks_e2)*10)
        plt.xlim(0.5, 5e3); plt.grid(True, which='both', linestyle='--'); plt.legend()
        plt.title("Comparison: $E^2$ dN/dE"); plt.tight_layout(); plt.show(block=True); plt.close(fig2)

        # Histograms + PDF overlays with robust end-trim
        for r in results:
            method = r["method"]
            figH = plt.figure(figsize=(10,6))

            # --- safe histogram for tiny samples / zero-range data ---
            md = np.asarray(r["mdist"], dtype=float)
            md = md[np.isfinite(md)]

            if md.size < 2 or (md.size > 0 and md.min() == md.max()):
                center = md[0] if md.size else 0.0
                eps = abs(center)*1e-9 if center != 0 else 1e-9
                counts, bins, _ = plt.hist(
                    md, bins=1, range=(center - eps, center + eps),
                    alpha=0.7, edgecolor='k',
                    label=f'{distribution_method} samples ({r["label_param"]})'
                )
            else:
                q25, q75 = np.percentile(md, [25, 75])
                iqr = q75 - q25
                if iqr > 0:
                    bw = 2 * iqr * md.size ** (-1/3)  # Freedman–Diaconis
                    k = int(np.clip(np.ceil((md.max() - md.min()) / bw), 1, 50))
                else:
                    k = int(np.clip(np.sqrt(md.size), 1, 50))
                counts, bins, _ = plt.hist(
                    md, bins=k, alpha=0.7, edgecolor='k',
                    label=f'{distribution_method} samples ({r["label_param"]})'
                )

            bin_widths  = (bins[1:] - bins[:-1])
            ref_width   = float(np.median(bin_widths)) if bin_widths.size else 1.0

            if method == "gaussian":
                sigma_x = r["params"]["sigma_x"]
                x = x_collapse
                mf = mass_function(delta_l(x, kappa, delta_c, gamma_p), sigma_x, delta_c, gamma_p)
                mf = np.where(np.isfinite(mf) & (mf > 0), mf, 0.0)
                if mf.sum() > 0:
                    probabilities = mf / mf.sum()
                    r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]
                    scale = r["peak"] / r_mode
                    dx = x[1] - x[0]; dm = dx * scale
                    pdf_mass = probabilities / dm
                    m_line = x * scale
                    mask_hist = (m_line >= bins[0]) & (m_line <= bins[-1])
                    m_plot = m_line[mask_hist]
                    y_plot = (pdf_mass * ref_width * len(r["mdist"]))[mask_hist]
                    keep_idx = _trim_right_spike(m_plot, y_plot, up_thresh=1.35, down_thresh=0.35, max_trim_frac=0.10)
                    plt.plot(m_plot[:keep_idx+1], y_plot[:keep_idx+1], 'r--', lw=2, label='Underlying PDF (counts)')

            elif method == "non_gaussian":
                sigma_X = r["params"]["sigma_X"]; ratio = 0.75; sigma_Y = ratio * sigma_X
                x = x_collapse
                mf = mass_function_exact(delta_l(x, kappa, delta_c, gamma_p), sigma_X, sigma_Y, delta_c, gamma_p)
                mf = np.where(np.isfinite(mf) & (mf > 0), mf, 0.0)
                if mf.sum() > 0:
                    probabilities = mf / mf.sum()
                    r_mode = x[np.argmax(mf)] if np.any(mf) else x[len(x)//2]
                    scale = r["peak"] / r_mode
                    dx = x[1] - x[0]; dm = dx * scale
                    pdf_mass = probabilities / dm
                    m_line = x * scale
                    mask_hist = (m_line >= bins[0]) & (m_line <= bins[-1])
                    m_plot = m_line[mask_hist]
                    y_plot = (pdf_mass * ref_width * len(r["mdist"]))[mask_hist]
                    keep_idx = _trim_right_spike(m_plot, y_plot, up_thresh=1.35, down_thresh=0.35, max_trim_frac=0.10)
                    plt.plot(m_plot[:keep_idx+1], y_plot[:keep_idx+1], 'r--', lw=2, label='Underlying PDF (counts)')

            else:
                sigma_ln = r["params"]["sigma_ln"]; mu_eff = np.log(r["peak"]) + sigma_ln**2
                mlo_tail = np.exp(mu_eff - 6.0*sigma_ln); mhi_tail = np.exp(mu_eff + 6.0*sigma_ln)
                mlo = min(bins[0],  mlo_tail); mhi = max(bins[-1], mhi_tail)
                m_plot = np.logspace(np.log10(mlo), np.log10(mhi), 2000)
                pdf = (1.0/(m_plot*sigma_ln*np.sqrt(2*np.pi))) * np.exp( - (np.log(m_plot)-mu_eff)**2 / (2*sigma_ln**2) )
                y_plot = pdf * ref_width * len(r["mdist"])
                plt.plot(m_plot, y_plot, 'r--', lw=2, label='Underlying PDF (counts)')
                plt.legend(title=f"σ={sigma_ln:.3f}")

            plt.xlabel('Simulated PBH Mass (g)')
            plt.ylabel('Count')
            plt.title(f'Mass Distribution & PDF for Peak {r["peak"]:.2e} g')
            plt.grid(True, which='both', linestyle='--')
            plt.legend()
            plt.tight_layout()
            plt.show(block=True); plt.close(figH)

    if not results: return

    sv = user_input("\nSave distributed results? (y/n): ").lower()
    if sv in ['y','yes']:
        print("Select spectra by index to save:")
        for i, r in enumerate(results, 1):
            disp = f"{r['peak']:.2e} g  [{r['label_param']}]"
            print(f" {i}: {disp}")
        idxs = user_input("Enter comma-separated (i.e. 1,3,5,...) or '0' to save ALL listed spectra: ")
        if idxs.strip() == '0':
            picks = list(range(1, len(results)+1))
        else:
            try:
                picks = [int(x) for x in idxs.split(',')]
            except Exception:
                err("Invalid indices; skipping save."); picks = []
        for idx in picks:
            if 1 <= idx <= len(results):
                r = results[idx-1]
                method = r["method"]; pk = r["peak"]
                E = r["E"]; spec = r["spectrum"]; md = r["mdist"]
                if method in ("gaussian", "non_gaussian"):
                    sigma_str = f"{r['params'].get('sigma_x', r['params'].get('sigma_X')):.3g}"
                    subdir = f"{pk:.1e}_{sigma_str}"
                else:
                    subdir = f"{pk:.1e}_{r['params']['sigma_ln']:.3g}"
                base_sub = ("gaussian" if method=="gaussian" else "non_gaussian" if method=="non_gaussian" else "lognormal")
                od = os.path.join(RESULTS_DIR, base_sub, subdir)
                os.makedirs(od, exist_ok=True)
                np.savetxt(os.path.join(od, "distributed_spectrum.txt"),
                           np.column_stack((E, spec)),
                           header="E_gamma(MeV)   Spectrum", fmt="%e")
                sorted_md = sorted(md); idxs_arr  = np.arange(1, len(sorted_md)+1)
                head = f"Mean mass (g): {np.mean(md):.2e}"
                np.savetxt(os.path.join(od, "mass_distribution.txt"),
                           np.column_stack((idxs_arr, sorted_md)),
                           header=head, fmt="%d %e")
                print(f"Saved → {od}")

# ---------------------------
# Monochromatic helper
# ---------------------------
def generate_monochromatic_for_mass(mval, DATA_DIR, MONO_RESULTS_DIR, tol=1e-12):
    masses, names = discover_mass_folders(DATA_DIR)
    if not masses:
        raise RuntimeError("No valid mass folders in data directory.")
    snapped = snap_to_available(mval, masses, tol=tol)
    exact   = snapped is not None
    mkey    = float(snapped if exact else mval)
    kind    = "pre-rendered" if exact else "interpolated"
    subdir  = f"{mkey:.1e}_{kind}"
    outdir  = os.path.join(MONO_RESULTS_DIR, subdir)
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, "Spectrum_components_and_total.txt")

    if not os.path.isfile(outfile):
        first_name = names[0]
        first = load_spectra_components(os.path.join(DATA_DIR, first_name))
        E_ref = first['energy_primary']

        if exact:
            idx_exact = np.where(np.isclose(masses, snapped, rtol=0, atol=0))[0][0]
            S = load_spectra_components(os.path.join(DATA_DIR, names[idx_exact]))
            d = S['direct_gamma_primary']
            s = np.interp(E_ref, S['energy_secondary'], S['direct_gamma_secondary'], left=0, right=0)
            i = S['IFA_primary'] + np.interp(E_ref, S['energy_secondary'], S['IFA_secondary'], left=0, right=0)
            f = S['FSR_primary'] + np.interp(E_ref, S['energy_secondary'], S['FSR_secondary'], left=0, right=0)
        else:
            mats = {k: [] for k in ('direct','secondary','inflight','finalstate')}
            Emax = []
            for i, m in enumerate(masses):
                sub = os.path.join(DATA_DIR, names[i])
                S = load_spectra_components(sub)
                mats['direct'].append(S['direct_gamma_primary'])
                mats['secondary'].append(np.interp(E_ref, S['energy_secondary'], S['direct_gamma_secondary'], left=0, right=0))
                mats['inflight'].append(S['IFA_primary'] + np.interp(E_ref, S['energy_secondary'], S['IFA_secondary'], left=0, right=0))
                mats['finalstate'].append(S['FSR_primary'] + np.interp(E_ref, S['energy_secondary'], S['FSR_secondary'], left=0, right=0))
                p = np.genfromtxt(os.path.join(sub, "inflight_annihilation_prim.txt"))
                q = np.genfromtxt(os.path.join(sub, "inflight_annihilation_sec.txt"))
                Emax.append(max(p[:,0].max() if p.size else 0, q[:,0].max() if q.size else 0))
            tiny = 1e-300
            logM_all = np.log(masses)
            logE     = np.log(E_ref)
            ld = np.log(np.where(np.vstack(mats['direct'])>tiny, mats['direct'], tiny))
            ls = np.log(np.where(np.vstack(mats['secondary'])>tiny, mats['secondary'], tiny))
            li = np.log(np.where(np.vstack(mats['inflight'])>tiny, mats['inflight'], tiny))
            lf = np.log(np.where(np.vstack(mats['finalstate'])>tiny, mats['finalstate'], tiny))
            sp_d = RectBivariateSpline(logM_all, logE, ld, kx=1, ky=3, s=0)
            sp_s = RectBivariateSpline(logM_all, logE, ls, kx=1, ky=3, s=0)
            sp_i = RectBivariateSpline(logM_all, logE, li, kx=1, ky=3, s=0)
            sp_f = RectBivariateSpline(logM_all, logE, lf, kx=1, ky=3, s=0)
            logm = np.log(mval)
            d = np.exp(sp_d(logm, logE, grid=False))
            s = np.exp(sp_s(logm, logE, grid=False))
            i = np.exp(sp_i(logm, logE, grid=False))
            f = np.exp(sp_f(logm, logE, grid=False))

        tot = d + s + i + f
        for arr in (d, s, i, f, tot): arr[arr < 1e-299] = 0.0
        np.savetxt(outfile, np.column_stack((E_ref, d, s, i, f, tot)),
                   header="E_gamma(MeV) Direct Secondary Inflight FinalState Total", fmt="%e")
    return outfile

# ---------------------------
# View previous
# ---------------------------
def view_previous_spectra():
    cat_map = {
        '2': (GAUSSIAN_METHOD,     os.path.join(RESULTS_DIR,'gaussian')),
        '3': (NON_GAUSSIAN_METHOD, os.path.join(RESULTS_DIR,'non_gaussian')),
        '4': (LOGNORMAL_METHOD,    os.path.join(RESULTS_DIR,'lognormal'))
    }
    sels = []
    while True:
        print("\nView Previous — choose:")
        print(" 1: Monochromatic Distribution")
        print(" 2: Gaussian collapse")
        print(" 3: Non-Gaussian collapse")
        print(" 4: Log-Normal distribution")
        print(" 0: Plot all Queued | b: Back")
        c = user_input("Choice: ").strip().lower()
        if c in ('b','back'): return
        if c == '1':
            try:
                masses = user_input("Enter mass(es) to queue (comma-separated): ", allow_back=True)
            except BackRequested:
                # Return to the view-previous menu
                continue
            for mstr in masses.split(','):
                try:
                    mval = float(mstr)
                    path = generate_monochromatic_for_mass(mval, DATA_DIR, MONO_RESULTS_DIR)
                    data = np.loadtxt(path, skiprows=1)
                    E, tot = data[:,0], data[:,-1]
                    sels.append((f"Mono {mval:.2e}", E, tot))
                    print(f"Queued Mono {mval:.2e}")
                except Exception as e:
                    warn(f"Skipping token '{mstr.strip()}': {e}")
            continue
        if c in cat_map:
            lbl, base = cat_map[c]
            runs = list_saved_runs(base)
            if not runs:
                print("None saved yet."); continue
            print(f"Available in {lbl}:")
            for i,run in enumerate(runs, 1):
                print(f" {i}: {run}")
            picks = user_input("Enter indices to queue (comma-separated): ")
            for idx in picks.split(','):
                try:
                    run = runs[int(idx)-1]
                    fn  = os.path.join(base, run, "distributed_spectrum.txt")
                    data = np.loadtxt(fn, skiprows=1)
                    E, tot = data[:,0], data[:,1]
                    sels.append((f"{lbl} {run}", E, tot))
                    print(f"Queued {lbl} {run}")
                except Exception as e:
                    warn(f"Skipping selection '{idx.strip()}': {e}")
            continue
        if c == '0': break
        print("Invalid choice.")

    if not sels:
        print("None queued."); return

    peaks_dn  = [sp.max() for _,_,sp in sels]
    ymax1 = max(peaks_dn)*10; ymin1 = min(peaks_dn)/1e3
    peaks_e2  = [ (E**2 * sp).max() for _,E,sp in sels ]
    ymax2 = max(peaks_e2)*10; ymin2 = min(peaks_e2)/1e3

    figA = plt.figure(figsize=(10,7))
    for name, E, sp in sels:
        m = (sp>0); plt.plot(E[m], sp[m], lw=2, label=name)
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel(r'$E_\gamma$ (MeV)'); plt.ylabel(r'$dN_\gamma/dE_\gamma$')
    plt.ylim(ymin1, ymax1); plt.xlim(0.5, 5e3)
    plt.grid(True, which='both', linestyle='--'); plt.legend()
    plt.title("Comparison: dN/dE"); plt.tight_layout(); plt.show(block=True); plt.close(figA)

    figB = plt.figure(figsize=(10,7))
    for name, E, sp in sels:
        m = (sp>0); plt.plot(E[m], E[m]**2 * sp[m], lw=2, label=name)
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel(r'$E_\gamma$ (MeV)'); plt.ylabel(r'$E^2\,dN_\gamma/dE_\gamma$')
    plt.ylim(ymin2, ymax2); plt.xlim(0.5, 5e3)
    plt.grid(True, which='both', linestyle='--'); plt.legend()
    plt.title("Comparison: $E^2$ dN/dE"); plt.tight_layout(); plt.show(block=True); plt.close(figB)

# ---------------------------
# UI
# ---------------------------
def show_start_screen():
    print("\n" + Fore.CYAN + Style.BRIGHT + "╔════════════════════════════════════════════════════════╗")
    print(        Fore.CYAN + Style.BRIGHT + "║               GammaPBHPlotter: PBH Spectrum Tool       ║")
    print(        Fore.CYAN + Style.BRIGHT + "║                        Version 1.0.0                   ║")
    print(        Fore.CYAN + Style.BRIGHT + "╚════════════════════════════════════════════════════════╝")
    print()
    print("Analyze and visualize Hawking radiation spectra of primordial black holes.\n")
    print(Fore.YELLOW + "📄 Associated Publication:" + Style.RESET_ALL)
    print("   John Carlini & Ilias Cholis — Particle Astrophysics Research\n")
    print("Type 'b' or 'back' to return to main menu or 'q' to exit.")

def main():
    show_start_screen()
    while True:
        print("\nSelect:")
        print("1: Monochromatic spectra")
        print(f"2: Distributed spectra ({GAUSSIAN_METHOD})")
        print(f"3: Distributed spectra ({NON_GAUSSIAN_METHOD})")
        print(f"4: Distributed spectra ({LOGNORMAL_METHOD})")
        print("5: View previous spectra")
        print("0: Exit")
        choice = user_input("Choice: ").strip().lower()
        if choice=='1':
            monochromatic_spectra()
        elif choice=='2':
            distributed_spectrum(GAUSSIAN_METHOD)
        elif choice=='3':
            distributed_spectrum(NON_GAUSSIAN_METHOD)
        elif choice=='4':
            distributed_spectrum(LOGNORMAL_METHOD)
        elif choice=='5':
            view_previous_spectra()
        elif choice in ['0','exit','q']:
            print("Goodbye."); break
        else:
            print("Invalid; try again.")

if __name__=='__main__':
    try:
        main()
    except Exception:
        import traceback
        traceback.print_exc()
        input("\nAn error occurred. Press Enter to exit…")
