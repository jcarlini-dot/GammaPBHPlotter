import os
import math
import shutil
import types
import builtins
import pathlib as _pl

import numpy as np
import pytest

# Force non-interactive matplotlib backend for all tests
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

# --- tiny synthetic BlackHawk dataset helpers --------------------------------

REQUIRED_FILES = [
    "instantaneous_primary_spectra.txt",
    "instantaneous_secondary_spectra.txt",
    "inflight_annihilation_prim.txt",
    "inflight_annihilation_sec.txt",
    "final_state_radiation_prim.txt",
    "final_state_radiation_sec.txt",
]

def _write_primary(path):
    # 2 header lines, then >= 124 rows because cli slices [123:]
    # Energies in GeV; we'll cover 0.001–5 GeV (~1–5000 MeV after conversion)
    N = 300
    Egev = np.linspace(0.001, 5.0, N)
    # Simple, smooth spectrum ~ E^-2 (per GeV); CLI converts to MeV^-1
    y = 1.0 / np.maximum(Egev, 1e-9)**2
    with open(path, "w") as f:
        f.write("# E(GeV)  dN/dE (GeV^-1 s^-1)\n")
        f.write("# synthetic test data\n")
        for e, v in zip(Egev, y):
            f.write(f"{e:.6e} {v:.6e}\n")

def _write_secondary(path):
    # 1 header, energies in MeV over same span
    N = 300
    Emev = np.linspace(1.0, 5000.0, N)
    # Smooth shape ~ E^-1 * exp(-E/3000)
    y = (1.0 / np.maximum(Emev, 1.0)) * np.exp(-Emev / 3000.0)
    with open(path, "w") as f:
        f.write("# E(MeV)  dN/dE (MeV^-1 s^-1)\n")
        for e, v in zip(Emev, y):
            f.write(f"{e:.6e} {v:.6e}\n")

def _write_xy_lenient(path, with_singleton_header=False):
    # Optionally start with a lone number to exercise load_xy_lenient
    with open(path, "w") as f:
        if with_singleton_header:
            f.write("123\n")
        E = np.linspace(1.0, 5000.0, 120)
        y = 1e-5 * np.exp(-E / 4000.0)
        for e, v in zip(E, y):
            f.write(f"{e:.6e} {v:.6e}\n")

def _make_mass_dir(root, mass_label):
    d = os.path.join(root, mass_label)
    os.makedirs(d, exist_ok=True)
    _write_primary(os.path.join(d, "instantaneous_primary_spectra.txt"))
    _write_secondary(os.path.join(d, "instantaneous_secondary_spectra.txt"))
    _write_xy_lenient(os.path.join(d, "inflight_annihilation_prim.txt"), with_singleton_header=True)
    _write_xy_lenient(os.path.join(d, "inflight_annihilation_sec.txt"), with_singleton_header=True)
    _write_xy_lenient(os.path.join(d, "final_state_radiation_prim.txt"))
    _write_xy_lenient(os.path.join(d, "final_state_radiation_sec.txt"))
    return d

@pytest.fixture(scope="session")
def synthetic_blackhawk(tmp_path_factory):
    """
    Creates a minimal blackhawk_data tree with 3 mass folders.
    """
    root = tmp_path_factory.mktemp("bh_data")
    # Use masses that are far enough apart to exercise interpolation
    for mass in ("1.00e+16", "3.00e+16", "1.00e+17"):
        _make_mass_dir(str(root), mass)
    return str(root)

@pytest.fixture
def temp_results(tmp_path):
    out = tmp_path / "results"
    out.mkdir(parents=True, exist_ok=True)
    # Subdirs mirroring the CLI layout
    for sub in ("monochromatic", "gaussian", "non_gaussian", "lognormal", "custom_equation"):
        (out / sub).mkdir(exist_ok=True)
    return str(out)

@pytest.fixture(autouse=True)
def no_show(monkeypatch):
    # Prevent GUI/interactive windows during tests
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
    yield

@pytest.fixture
def cli_module(synthetic_blackhawk, temp_results, monkeypatch):
    """
    Import gammapbh.cli with DATA_DIR/RESULTS_DIR pointed at our temp fixtures.
    """
    import importlib
    import sys
    # Ensure project src is importable if running from repo root
    sys.path.insert(0, str(_pl.Path("src").resolve()))
    cli = importlib.import_module("gammapbh.cli")

    # Redirect globals to our temp fixtures
    monkeypatch.setattr(cli, "DATA_DIR", synthetic_blackhawk)
    monkeypatch.setattr(cli, "RESULTS_DIR", temp_results)
    monkeypatch.setattr(cli, "MONO_RESULTS_DIR", os.path.join(temp_results, "monochromatic"))
    monkeypatch.setattr(cli, "GAUSS_RESULTS_DIR", os.path.join(temp_results, "gaussian"))
    monkeypatch.setattr(cli, "NGAUSS_RESULTS_DIR", os.path.join(temp_results, "non_gaussian"))
    monkeypatch.setattr(cli, "LOGN_RESULTS_DIR", os.path.join(temp_results, "lognormal"))
    monkeypatch.setattr(cli, "CUSTOM_RESULTS_DIR", os.path.join(temp_results, "custom_equation"))
    return cli

class InputFeeder:
    """
    Feeds predetermined responses to cli.user_input.
    Any unexpected prompt will raise StopIteration for easier debugging.
    """
    def __init__(self, answers):
        self._it = iter(answers)

    def __call__(self, prompt, *, allow_back=False, allow_exit=True):
        try:
            return next(self._it)
        except StopIteration:
            # If the CLI asks for more than we expected, fail loudly
            raise AssertionError(f"No more canned answers for prompt: {prompt!r}")

@pytest.fixture
def feed_inputs(monkeypatch, cli_module):
    """
    Factory that installs a feeder on gammapbh.cli.user_input.
    Use: feeder = feed_inputs(["ans1", "ans2"]); then call your CLI function.
    """
    def _install(answers):
        feeder = InputFeeder(answers)
        monkeypatch.setattr(cli_module, "user_input", feeder)
        return feeder
    return _install
