# GammaPBHPlotter Version 1.1.4 (22 November 2025)
-----------------------------------
By John Carlini (jcarlini@oakland.edu) and Ilias Cholis (cholis@oakland.edu)

INTRODUCTION
-----------------------------------
The most recent version of this program can be obtained from: <br>
https://test.pypi.org/project/gammapbh <br>
https://zenodo.org/records/16944093 <br>
https://pypi.org/project/gammapbh/ <br>

This Python package is designed to simulate, display, and record the Hawking gamma-ray differential spectra per unit time (d^2 Nγ/(dEγ dt)) of primordial black holes (PBHs) in units of inverse megaelectron volts per second. The mass range of simulated PBHs is between 5×10^13 and 1×10^19 grams. It does this through a combination of interpolating direct Hawking radiation (DHR) spectral data from the existing software BlackHawk, as well as computations of the final state radiation (FSR) from electrons and positrons, and the energy produced by the annihilation of said positrons with electrons in the interstellar medium, referred to as inflight annihilation (IFA).

This software was designed for use by physicists and astronomers as both a comprehensive and user-friendly means of modelling different distributions of PBHs. These results can be compared to any excess gamma-rays detected from certain regions of space. Matches could be used as evidence not only for the presence of PBHs, but their number, density, and distribution. 

DISTRIBUTION METHODS
-----------------------------------
Gamma-ray Hawking Spectra can be generated in one of many distribution methods. In monochromatic distribution all black holes possess an identical mass and by extension an identical gamma-ray Hawking spectrum. If the mass of the simulated black hole happens to align with one of the 56 pre-rendered spectra generated via BlackHawk (for DHR) and our calculations of the FSR and IFA, then the resulting spectrum is presented as is. For any other mass within the appropriate range, a new spectrum is interpolated based on that existing data.

All other means of distribution are simulated by producing a number of randomly generated black hole masses according to a probability density function. The individual spectra of each black hole are then simulated in a similar manner to the monochromatic method before being added up to produce a final result. The number of simulated black holes is inputted by the user, but a sample size of at least 1000 is a recommended minimum for accuracy. It is also worth noting that in order to better understand the average contribution of singular black holes and account for different needed sample sizes, the distributed spectra this software produces have been divided in amplitude by their sample size. So if a user were to generate a sample size of 1000, it is important to remember that the results seen and saved are only an average per black hole and would be 3 orders of magnitude lower in amplitude than the total radiation that particular number of black holes of those masses would actually produce.

What differentiates these distribution methods is what specific probability density function they follow. The Gaussian and non-Gaussian collapse PDFs are based upon a model of PBH formation and early universe structure from the paper "The Formation Probability of Primordial Black Holes" by Matteo Biagetti et al. Due to the specific limitations of that model, it only remains accurate for a limited range of values for the standard deviation (σ). Only values within that range may be simulated by this software.

0.03 < σ < 0.255 for the case of the Gaussian collapse. <br>
0.04 < σ < 0.16 for the case of non-Gaussian collapse

For the lognormal distribution, it is a simpler and more malleable model which can accommodate values of standard deviations as long as σ > 0. That being said, values of σ=2 or lower are recommended for utility as any higher of a spread would most of the distribution lying outside of our mass range.

Custom distributions allow a user to enter any probability density function required in the form of a Python expression f(m). The variable "m" refers to the PBH mass in grams. Identifiers for Numpy "np." are not required and expressions like "sin", "log", "exp" etc can be entered in directly. The resulting function will always be normalized to 1 for the purposes of calculating probability. Constants may be entered into the equation, but the user will be prompted after entering to define the definition for each of them.

Running these simulations will take additional time proportional to the sample sized used. Since monochromatic distributions only have a sample size of 1, they are effectively instant.

RUNNING SIMULATIONS
-----------------------------------
No matter which of the distribution methods is needed, the act of simulating them is a similar process for users. Upon selecting their desired method, the program will send a text prompt to input the mass values most likely to appear in the distribution (referred to as your peak masses). Peak masses can be entered individually or as a series of comma separated numbers in scientific notation. I.e. 1e15, 2.5e16, 3.75e17, etc. Since masses can only be generated within the limits of 5×10^13 and 1×10^19 g, placing a peak mass too close to said limits while using a non-monochromatic distribution method can cause a significant number of the black holes to fall outside the range, have their spectra treated as 0, and cause the overall data to lose accuracy. This is less of an issue when the peak mass is near the higher end of masses (5e18 g) as masses above that value have such low values of Hawking radiation that counting them as 0 is an imperceptible difference in almost all use cases. Additionally, outside of a monochromatic distribution, the peak mass does not always coincide with the average mass. For this reason, the mean mass of simulated black holes is provided as well in graphs as well as in saved results.

If a peak mass is entered individually, two graphs will appear. One will show the number of gamma-ray photons emitted per unit energy and unit time (or dNγ/dEγ) in units of Inverse Megaelectron Volts and inverse seconds on the y-axis and Energy (E) in units of Megaelectron volts. The next graph is opened by closing the previous one and displays that same data in units of Megaelectron Volts per second. That is done by multiplying the y axis (dNγ/dEγ) of each data point by the x axis squared (E²). The first graph provides data in a form more useful for the simulation of Hawking radiation, while the second provides data in the form of the luminosity in Megaelectron Volts per second. If multiple masses are entered, the resulting spectra are presented in separate graphs. Each graph will appear once the previous one is closed in the order they were entered. Once through, all spectra as well as their cumulative sum will be presented in one final graph in units of MeV s^-1.

SAVING RESULTS
-----------------------------------
Once the final graph produced by any simulation has been closed, the program will give the user a y/n prompt of whether they would like to save their results or not. If "y" is entered, an indexed list of the entered peak masses will appear alongside a prompt asking the user which simulations to save. This task is done by entering a single number, a comma separated list of numbers, or simply pressing "0" to save all of them. Once finished, the results will automatically be saved as a .txt file in a destination folder named "results" and a specific subfolder depending on the method used to generate them. If "n" is selected, the user returns to the main menu.

Monochromatic distribution 	=	".../results/monochromatic" <br>
Gaussian distribution 		=	".../results/gaussian" <br>
Lognormal distribution 		=	".../results/lognormal" <br>
non-Gaussian collapse		=	".../results/non_gaussian" <br>
custom distribution		=	"...results/custom_equation" <br>

Within the appropriate subfolder, the results are saved as another subfolder named after the peak mass used to generate them with three significant figures. For example, a gaussian distribution with peak mass 3.1415e15 grams would be saved under ".../results/gaussian/3.14e+15". Be careful to back up your files when performing multiple simulations of identical or sufficiently close masses via the same method. You may overwrite your previous data. 

Spectra generated from monochromatic distribution provide the spectra for each individual component of the spectrum (Direct Primary, Direct Secondary, Inflight Annihilation, and Final State Radiation) all in their own columns of the same file named "spectrum components". Spectra from the other three methods instead produce two files. One called "distributed_spectrum" which includes a one column spectrum of the total hawking gamma-ray spectrum as seen in the graph. Additionally, there is also the "mass_distribution' file which lists all the masses generated by the simulation as well as the average mass.

VIEWING PREVIOUS SPECTRA
-----------------------------------
If it is desired to compare spectra from different PBH mass distributions, the "view previous spectra" feature on the main menu is provided. Once selected, the user is presented with a screen similar to the main menu. The user first selects the type of PBH mass distribution, i.e. monochromatic, Gaussian, non-Gaussian, or lognormal. Then the user selects a peak mass. For the monochromatic distribution, the user may input any mass within the allowed range. For any of the other three cases, the program provides an indexed list of saved spectral files. Once all the desired file(s) are selected, the user will see a message which writes "→ Queued: {Method} {Peak Mass}" or "→ Queued: Gaussian Distribution 3.14e+15" to use the earlier example. This can be done multiple times for multiple different PBH distribution types. Once everything the user wishes to graph is selected, the user needs only to press 0 from the "previous spectra menu" to view all of them in two graphs. One of them is in units of MeV^-1 s^-1, the other in MeV s^-1. 


REQUIREMENTS  
----------------------------------- 
All that is needed to run GammaPBHPlotter is the following:
- Python 3.9 or newer  
- Internet connection (for first-time installation of dependencies)  

These modules will be automatically installed by pip if not already present:
- colorama
- numpy  
- matplotlib 
- tqdm  
- scipy 

INSTALLATION STEPS  
-----------------------------------  

Option A — Recommended (via pip):  
```text
	pip install gammapbh  
```
You can then run the text user-interface directly from your terminal with:
```text
	gammapbh-tui  
```
or equivalently:
```text
	python -m gammapbh-tui  
```
Option B — Manual build (from source):
```text
	git clone https://github.com/jcarlini-dot/GammaPBHPlotter 
	cd GammaPBHPlotter  
	python -m pip install .  
```
To verify a successful installation: 
```text
	python -c "import gammapbh, importlib.metadata as md; print(gammapbh.__version__, md.version('gammapbh'))"
```

TUI EXAMPLE RUNS
-----------------------------------
To test if the package is successfully installed and operational, please copy and paste these blocks of inputs into your device's command prompt
Example A — Monochromatic spectra
```text
  #1) Start Package TUI 
	gammapbh
  #2) Pick to generate monochromatic spectra 
	1  
  #3) Enter masses (g) within the available grid
	3.14e15, 1.4e14
  #4) When prompted, choose to save
	y
  #5) Choose to save all masses
	0
  #Outputs:
    #results/monochromatic/3.14e+15_spectrum.txt
    #results/monochromatic/1.40e+14_spectrum.txt
    #Output File Columns: E_gamma(MeV)  Direct  Secondary  Inflight  FinalState  Total
```
Example B — Log-normal distributed spectrum
  ```text
  #1) Start Package  
	gammapbh
  #2) Pick to generate lognormal spectra
	4
  #3) Enter the Peak PBH mass (g) within available grid 
	3e16
  #4) Enter target "N"
	2000
  #5) Enter σ: 
	0.6
  #6) Save results when prompted
	y
  #7) Choose to save one spectrum
	1
  #Outputs:
    #results/lognormal/peak_3.00e+16_σ0.6_N2000/distributed_spectrum.txt   (E_gamma(MeV), TotalSpectrum)
        #Output File Columns: (E_gamma(MeV), TotalSpectrum)
    #results/lognormal/peak_3.00e+16_σ0.6_N2000/mass_distribution.txt 
        #Output File Columns: (N sampled masses in g)
```
## Quick Start (CLI & API)
If a user wishes to not use GammaPBHplotter via the interactive TUI, then other options include:
- a **non-interactive CLI** (great for scripts/automation), and  
- the **Python API** (for notebooks or custom pipelines).  

All examples assume **v1.1.4** is installed and on your `PATH`.

> **Units:** energy `E` is in **MeV**. Differential spectra use **MeV⁻¹ s⁻¹**; the SED-style view uses **MeV s⁻¹** via \(E^2\,\mathrm{d}N/\mathrm{d}E\).

### Command-line (non-interactive)

```bash
# Show help and available options
python -m gammapbh --help

# List the available PBH mass grid entries bundled with the package
python -m gammapbh --list-masses

# Write a CSV of energy + all components (aligned to the "secondary" grid)
python -m gammapbh --mass 5e13 --align secondary --csv --save out.csv

# Save a plot of the total spectrum (no GUI window)
python -m gammapbh --mass 5e13 --align secondary --plot --save out.png --no-show
```

**Notes**
- `--align` controls the common energy grid:
  - `secondary` *(default)*, `primary`, `union` (all unique points), or `intersection` (only overlapping points).
- `--csv` writes **E** plus each component column (e.g., `direct_gamma_*`, `IFA_*`, `FSR_*`, `total_*`, and `total` if present).
- `--plot` plots the **total** spectrum on the chosen grid; combine with `--save` and `--no-show` for headless runs.

### Python API

```python
import numpy as np
import matplotlib.pyplot as plt
import gammapbh as gp

# Discover the bundled mass grid and pick one mass
masses = gp.discover_mass_folders()
# masses can be a tuple (numbers, names) or a flat list; pick the first numeric value robustly:
if isinstance(masses, tuple) and len(masses) == 2 and masses[0]:
    mass = float(masses[0][0]) if isinstance(masses[0][0], (int, float)) else float(masses[1][0])
elif isinstance(masses, list) and masses:
    mass = float(masses[0]) if isinstance(masses[0], (int, float)) else float(masses[0])
else:
    raise RuntimeError("No masses found in the bundled grid.")

# Load spectra aligned to the 'secondary' grid (default). Returns energy array E and a dict of components.
E, comps = gp.load_spectra_components(mass, align="secondary")  # also: "primary", "union", "intersection"

# Access the total spectrum (if not directly provided, sum primary+secondary totals)
total = comps.get("total")
if total is None:
    tp, ts = comps.get("total_primary"), comps.get("total_secondary")
    if tp is not None and ts is not None:
        total = np.asarray(tp) + np.asarray(ts)
    else:
        raise RuntimeError("No total spectrum available in components.")

# Quick plot (log–log)
plt.figure()
plt.loglog(E, total, label=f"total, M={mass:.3e} g (secondary grid)")
plt.xlabel("E [MeV]")
plt.ylabel("dN/dE [MeV$^{-1}$ s$^{-1}$]")
plt.grid(True, which="both", ls=":")
plt.legend()
plt.tight_layout()
plt.show()

# (Optional) Write a simple CSV (E + all component columns)
import csv
with open("spectrum_out.csv", "w", newline="") as f:
    w = csv.writer(f)
    cols = ["E"] + sorted(comps.keys())
    w.writerow(cols)
    for i in range(len(E)):
        w.writerow([E[i]] + [np.asarray(comps[k])[i] for k in cols[1:]])
```

> **Tip:** For distributed (stochastic) runs done via the TUI, results depend on random sampling size/seed. The CLI/API examples above operate on the bundled per-mass spectra and are deterministic for a given mass and alignment.


## Tests

This project uses `pytest` for the test suite, with optional tools for coverage, style, and types. Below are quick-start instructions for running and extending tests on Linux/macOS and Windows PowerShell.

### 1) Set up a virtual environment

**Linux / macOS (bash/zsh)**  
    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install --upgrade pip

**Windows PowerShell**  
    py -3 -m venv .venv
    .\.venv\Scripts\Activate.ps1
    python -m pip install --upgrade pip

### 2) Install the package and test dependencies

If the project exposes extras (recommended):

**Either**  
    pip install -e ".[dev]"        # includes test/linters/types/etc.  
**Or**  
    pip install -e ".[test]"       # just test deps

If you prefer explicit files:

    pip install -e .
    pip install -r requirements-dev.txt

> Tip: If you use `matplotlib` in tests, the suite sets a non-interactive backend (Agg). You can also export `MPLBACKEND=Agg` (bash) or `$env:MPLBACKEND='Agg'` (PowerShell) before running tests.

### 3) Run the full test suite

    pytest -q

To see verbose output:

    pytest -vv

To stop on first failure:

    pytest -x

### 4) Run subsets (markers / paths)

Run only unit tests in a directory:

    pytest -q tests/unit

Run only tests marked as “fast” (see markers below):

    pytest -m fast -q

Exclude slow tests:

    pytest -m "not slow" -q

Run a single test file or a single test:

    pytest tests/test_spectra.py -q
    pytest tests/test_spectra.py::test_primary_component_shapes -q

### 5) Coverage

Terminal coverage with missing lines:

    pytest --cov=gammapbh --cov-report=term-missing

Create an HTML report at `htmlcov/index.html`:

    pytest --cov=gammapbh --cov-report=term-missing --cov-report=html

Open the HTML report in your browser.

### 6) Style and type checks (optional but recommended)

If configured:

    ruff check .
    black --check .
    mypy src/ tests/

You can auto-fix style issues:

    ruff check . --fix
    black .

### 7) CLI smoke tests

Ensure the CLI responds and help prints without error:

    python -m gammapbh.cli --help
    gammapbh --help

Run a specific CLI path:

    pytest tests/test_cli.py -q

> For CLI snapshot tests, prefer comparing normalized text or small JSON outputs to avoid flakiness.

### 8) Plot/image tests (optional)

For deterministic plot tests, use `matplotlib`’s Agg backend. If you use image comparison:

- Prefer lightweight assertions on figure metadata (axes count, labels present, line counts) rather than pixel-exact diffs.
- If using `pytest-mpl`, store baseline images under `tests/baseline/` and run:

      pytest --mpl -q

- Keep numeric seeds fixed (e.g., `rng = np.random.default_rng(0)`) for any stochastic examples.

### 9) Floating-point comparisons

When asserting numerical equality, use tolerances:

    import numpy as np
    assert np.allclose(got, expected, rtol=1e-10, atol=0.0)

Choose tolerances appropriate to the algorithm and units. Prefer relative (`rtol`) for scale-invariant checks and absolute (`atol`) when values can be near zero.

### 10) Markers and skips

Declare markers in `pytest.ini`:

    [pytest]
    addopts = -ra
    testpaths = tests
    markers =
        slow: long-running tests (use -m "not slow" to skip)
        net: tests that require network access
        cli: command-line interface tests
        plot: plotting/image tests

Use them in tests:

    import pytest

    @pytest.mark.slow
    def test_heavy_integration():
        ...

    @pytest.mark.skipif(not HAVE_DATA, reason="requires local data assets")
    def test_with_optional_data():
        ...

### 11) Test layout & naming

Recommended structure:

    tests/
      unit/
        test_loaders.py
        test_mass_functions.py
        test_interpolation.py
      cli/
        test_cli.py
      plot/
        test_plots.py
      data/              # tiny fixtures only; larger assets should be generated or downloaded in CI

Pytest discovers tests matching `test_*.py` or `*_test.py`.

### 12) Writing new tests (examples)

Unit test for a mass PDF normalization:

    import numpy as np
    from gammapbh import mass_function_lognormal

    def test_lognormal_normalizes_to_one():
        m = np.logspace(13, 20, 2000)  # grams
        pdf = mass_function_lognormal(m, mu=1e15, sigma=0.5)
        # integrate in log-space to reduce bias
        area = np.trapz(pdf, x=np.log(m))
        assert np.allclose(area, 1.0, rtol=1e-3)

Interpolation monotonicity:

    from gammapbh import load_spectra_components

    def test_interpolation_monotonic_energy_grid(tmp_path):
        comp = load_spectra_components(mass_g=3e15, data_root="blackhawk_data")
        e = comp.energy_mev
        assert (np.diff(e) > 0).all()

CLI smoke test:

    import subprocess, sys

    def test_cli_help_runs():
        code = subprocess.call([sys.executable, "-m", "gammapbh.cli", "--help"])
        assert code == 0

### 13) Doctests (optional)

You can validate code examples in the README or docstrings:

    pytest --doctest-glob="*.md" README.md
    pytest --doctest-modules src/gammapbh

Keep examples fast and deterministic.

### 14) Continuous Integration (GitHub Actions example)

Create `.github/workflows/tests.yml`:

    name: tests
    on:
      push:
      pull_request:
    jobs:
      test:
        runs-on: ubuntu-latest
        steps:
          - uses: actions/checkout@v4
          - uses: actions/setup-python@v5
            with:
              python-version: "3.11"
          - name: Install
            run: |
              python -m pip install --upgrade pip
              pip install -e ".[dev]" || pip install -e . && pip install -r requirements-dev.txt
          - name: Run tests
            env:
              MPLBACKEND: Agg
            run: |
              pytest -q --cov=gammapbh --cov-report=term-missing

### 15) Tips for reliable tests

- Avoid network and large file dependencies; use tiny fixtures or generate data on the fly.
- Seed all randomness.
- Keep slow tests behind a `@pytest.mark.slow` marker.
- Assert on invariants (shapes, monotonicity, normalization) in addition to exact numbers.
- Prefer black-box API tests at the module boundary, and targeted unit tests for numerics.

If anything fails or is unclear, run with `-vv -ra` for more detail and include the command, OS, Python version, and full traceback when filing an issue.

INCLUDED FILES  
-----------------------------------  

Top-Level Project Structure:  
	GammaPBHPlotter/  
	│  
	├── pyproject.toml           (Build configuration for pip and PyPI)  
	├── LICENSE                  (GNU GPL v3 license)  
	├── README.txt               (This documentation file)  
	├── CITATION.cff             (Citation metadata for Zenodo and GitHub)  
	├── CHANGELOG.md             (Version history and updates)  
	│  
	├── src/  
	│   └── gammapbh/  
	│       ├── __init__.py  
	│       ├── __main__.py          (Enables 'python -m gammapbh')  
	│       ├── cli.py               (Primary program logic and CLI interface)  
	│       ├── blackhawk_data/      (Tables from the BlackHawk software)  
	│       └── results/             (Default save directory, auto-created)  
	│  
	└── tests/         (Unit tests and validation scripts)  


FILE AND FOLDER DESCRIPTIONS  
-----------------------------------  
	src/gammapbh/cli.py		= Core logic for spectrum interpolation, PDF sampling, and visualization.  
	src/gammapbh/__main__.py	= Enables running the software with "python -m gammapbh".  
	src/gammapbh/blackhawk_data/	= Precomputed spectral tables from BlackHawk (Arbey & Auffinger 2019, 2021).  
	src/gammapbh/results/		= Auto-generated output directory for saved spectra; initially empty.  
	LICENSE				= Terms of redistribution under GNU GPL v3.  
	CITATION.cff			= Citation metadata for Zenodo and scholarly references.  
	README.txt			= Full user guide (this file).  
	CHANGELOG.md			= Summarized update history for each release.  


VERSION HISTORY
-----------------------------------  

- **v1.0.0** – 2025-08-25  
  - First public release as an executable.

- **v1.1.0** – 2025-10-28  
  - Released as a Python package on PyPI.  
  - Added the ability to enter custom PDF equations.

- **v1.1.1** – 2025-10-29  
  - Fixed issues related to faulty data upload.

- **v1.1.3** – 2025-11-04  
  - Fixed bug involving unenforced mass limits while viewing prior monochromatic spectra.  
  - Tweaked plot labeling for increased legibility.  
  - Added the ability to view previously generated histograms.

- **v1.1.4** – 2025-11-22  
  - Reformatted paper and README for compatibility.  
  - Comment inputs beginning with a pound sign (`#`) are ignored, allowing example runs to be directly copy-pasted into a terminal.
  - Added Non-interactive CLI (`python -m gammapbh`)
  - Added Public Python API** (`gammapbh.api`)
  - Added Test suite (pytest)

Acknowledgements
-----------------------------------  

The program has been tested on windows 11, Mac, and Linux devices.

If you use GammaPBHPlotter to write a paper, please cite:

linktocitation.placeholder

As well as the paper published for the BlackHawk software.

A. Arbey and J. Auffinger, Eur. Phys. J. C79 (2019) 693, arXiv:1905.04268 [gr-qc]
A. Arbey and J. Auffinger, Eur. Phys. J. C81 (2021) 910, arXiv:2108.02737 [gr-qc]

And if you use the gaussian or non-gaussian collapse for your paper, please cite Biagetti et al.

M. Biagetti, V. De Luca, G. Franciolini, A. Kehagias and A. Riotto, Phys. Lett. B 820 (2021) 136602, arXiv:2105.07810 [astro-ph.CO].

LICENSE
-----------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any 
    later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    See <http://www.gnu.org/licenses/>.  
