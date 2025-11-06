# GammaPBHPlotter Version 1.1.2 (29 October 2025)
-----------------------------------
By John Carlini (jcarlini@oakland.edu) and Ilias Cholis (cholis@oakland.edu)

INTRODUCTION
-----------------------------------
The most recent version of this program can be obtained from: 
https://test.pypi.org/project/gammapbh
https://zenodo.org/records/16944093
https://test.pypi.org/project/gammapbh

This Python-coded software is designed to simulate, display, and record the Hawking gamma-ray differential spectra per unit time (d^2 Nγ/(dEγ dt)) of black holes in units of inverse megaelectron volts per second. The mass range of simulated black holes is between 5×10^13 and 1×10^19 grams. It does this through a combination of interpolating direct Hawking radiation (DHR) spectral data from the existing software BlackHawk, as well as computations of the final state radiation (FSR) from electrons and positrons, and the energy produced by the annihilation of said positrons with electrons in the interstellar medium, referred to as inflight annihilation (IFA).

DISTRIBUTION METHODS
-----------------------------------
Gamma-ray Hawking Spectra can be generated in one of many distribution methods. In monochromatic distribution all black holes possess an identical mass and by extension an identical gamma-ray Hawking spectrum. If the mass of the simulated black hole happens to align with one of the 56 pre-rendered spectra generated via BlackHawk (for DHR) and our calculations of the FSR and IFA, then the resulting spectrum is presented as is. For any other mass within the appropriate range, a new spectrum is interpolated based on that existing data.

All other means of distribution are simulated by producing a number of randomly generated black hole masses according to a probability density function. The individual spectra of each black hole are then simulated in a similar manner to the monochromatic method before being added up to produce a final result. The number of simulated black holes is inputted by the user, but a sample size of at least 1000 is a recommended minimum for accuracy. It is also worth noting that in order to better understand the average contribution of singular black holes and account for different needed sample sizes, the distributed spectra this software produces have been divided in amplitude by their sample size. So if a user were to generate a sample size of 1000, it is important to remember that the results seen and saved are only an average per black hole and would be 3 orders of magnitude lower in amplitude than the total radiation that particular number of black holes of those masses would actually produce.

What differentiates these distribution methods is what specific probability density function they follow. The Gaussian and non-Gaussian collapse PDFs are based upon a model of PBH formation and early universe structure from the paper "The Formation Probability of Primordial Black Holes" by Matteo Biagetti et al. Due to the specific limitations of that model, it only remains accurate for a limited range of values for the standard deviation (σ). Only values within that range may be simulated by this software.

0.03 < σ < 0.255 for the case of the Gaussian collapse.
0.04 < σ < 0.16 for the case of non-Gaussian collapse

For the lognormal distribution, it is a simpler and more malleable model which can accommodate values of standard deviations as long as σ > 0. That being said, values of σ=2 or lower are recommended for utility as any higher of a spread would most of the distribution lying outside of our mass range.

Custom distributions allow a user to enter any probability density function required in the form of a Python expression f(m). The variable "m" refers to the PBH mass in grams. Identifiers for Numpy "np." are not required and expressions like "sin", "log", "exp" etc can be entered in directly. The resulting function will always be normalized to 1 for the purposes of calculating probability. Constants may be entered into the equation, but the user will be prompted after entering to define the definition for each of them.

Running these simulations will take additional time proportional to the sample sized used. Since monochromatic distributions only have a sample size of 1, they are effectively instant.

RUNNING SIMULATIONS
-----------------------------------
No matter which of the distribution methods is needed, the act of simulating them is a similar process for users. Upon selecting their desired method, the program will send a text prompt to input the mass values most likely to appear in the distribution (referred to as your peak masses). Peak masses can be entered individually or as a series of comma separated numbers in scientific notation. I.e. 1e15, 2.5e16, 3.75e17, etc. Since masses can only be generated within the limits of 5e13 and 5e18 g, placing a peak mass too close to said limits while using a non-monochromatic distribution method can cause a significant number of the black holes to fall outside the range, have their spectra treated as 0, and cause the overall data to lose accuracy. This is less of an issue when the peak mass is near the higher end of masses (5e18 g) as masses above that value have such low values of Hawking radiation that counting them as 0 is an imperceptible difference in almost all use cases. Additionally, outside of a monochromatic distribution, the peak mass does not always coincide with the average mass. For this reason, the mean mass of simulated black holes is provided as well in graphs as well as in saved results.

If a peak mass is entered individually, two graphs will appear. One will show the number of gamma-ray photons emitted per unit energy and unit time (or dNγ/dEγ) in units of Inverse Megaelectron Volts and inverse seconds on the y-axis and Energy (E) in units of Megaelectron volts. The next graph is opened by closing the previous one and displays that same data in units of Megaelectron Volts per second. That is done by multiplying the y axis (dNγ/dEγ) of each data point by the x axis squared (E²). The first graph provides data in a form more useful for the simulation of Hawking radiation, while the second provides data in the form of the luminosity in Megaelectron Volts per second. If multiple masses are entered, the resulting spectra are presented in separate graphs. Each graph will appear once the previous one is closed in the order they were entered. Once through, all spectra as well as their cumulative sum will be presented in one final graph in units of MeV s^-1.

SAVING RESULTS
-----------------------------------
Once the final graph produced by any simulation has been closed, the program will give the user a y/n prompt of whether they would like to save their results or not. If "y" is entered, an indexed list of the entered peak masses will appear alongside a prompt asking the user which simulations to save. This task is done by entering a single number, a comma separated list of numbers, or simply pressing "0" to save all of them. Once finished, the results will automatically be saved as a .txt file in a destination folder named "results" and a specific subfolder depending on the method used to generate them. If "n" is selected, the user returns to the main menu.

Monochromatic distribution 	=	".../results/monochromatic"
Gaussian distribution 		=	".../results/gaussian"
Lognormal distribution 		=	".../results/lognormal"
non-Gaussian collapse		=	".../results/non_gaussian"
custom distribution		=	"...results/custom_equation"

Within the appropriate subfolder, the results are saved as another subfolder named after the peak mass used to generate them with three significant figures. For example, a gaussian distribution with peak mass 3.1415e15 grams would be saved under ".../results/gaussian/3.14e+15". Be careful to back up your files when performing multiple simulations of identical or sufficiently close masses via the same method. You may overwrite your previous data. 

Spectra generated from monochromatic distribution provide the spectra for each individual component of the spectrum (Direct Primary, Direct Secondary, Inflight Annihilation, and Final State Radiation) all in their own columns of the same file named "spectrum components". Spectra from the other three methods instead produce two files. One called "distributed_spectrum" which includes a one column spectrum of the total hawking gamma-ray spectrum as seen in the graph. Additionally, there is also the "mass_distribution' file which lists all the masses generated by the simulation as well as the average mass.

VIEWING PREVIOUS SPECTRA
-----------------------------------
If it is desired to compare spectra from different PBH mass distributions, the "view previous spectra" feature on the main menu is provided. Once selected, the user is presented with a screen similar to the main menu. The user first selects the type of PBH mass distribution, i.e. monochromatic, Gaussian, non-Gaussian, or lognormal. Then the user selects a peak mass. For the monochromatic distribution, the user may input any mass within the allowed range. For any of the other three cases, the program provides an indexed list of saved spectral files. Once all the desired file(s) are selected, the user will see a message which writes "→ Queued: {Method} {Peak Mass}" or "→ Queued: Gaussian Distribution 3.14e+15" to use the earlier example. This can be done multiple times for multiple different PBH distribution types. Once everything the user wishes to graph is selected, the user needs only to press 0 from the "previous spectra menu" to view all of them in two graphs. One of them is in units of MeV^-1 s^-1, the other in MeV s^-1. 

The program has been tested on windows 11, Mac, and Linux devices.

If you use GammaPBHPlotter to write a paper, please cite:

linktocitation.placeholder

As well as the paper published for the BlackHawk software.

A. Arbey and J. Auffinger, Eur. Phys. J. C79 (2019) 693, arXiv:1905.04268 [gr-qc]
A. Arbey and J. Auffinger, Eur. Phys. J. C81 (2021) 910, arXiv:2108.02737 [gr-qc]

And if you use the gaussian or non-gaussian collapse for your paper, please cite Biagetti et al.

M. Biagetti, V. De Luca, G. Franciolini, A. Kehagias and A. Riotto, Phys. Lett. B 820 (2021) 136602, arXiv:2105.07810 [astro-ph.CO].

arXiv:2105.07810

INSTALLATION  
-----------------------------------  

GammaPBHPlotter v1.1.0 is now distributed as a Python package and can be installed directly from PyPI using pip.  

REQUIREMENTS  
-----------------------------------  
- Python 3.9 or newer  
- Internet connection (for first-time installation of dependencies)  

The following modules will be automatically installed by pip if not already present:  
	colorama  
	numpy  
	matplotlib  
	tqdm  
	scipy  

INSTALLATION STEPS  
-----------------------------------  

Option A — Recommended (via pip):  
	pip install gammapbh  

You can then run the program directly from your terminal with:  
	gammapbh  

or equivalently:  
	python -m gammapbh  

Option B — Manual build (from source):  
	git clone https://github.com/jcarlini-dot/GammaPBHPlotter 
	cd GammaPBHPlotter  
	python -m pip install .  

To verify a successful installation:  
	python -c "import gammapbh, importlib.metadata as md; print(gammapbh.__version__, md.version('gammapbh'))"  

Both commands should print 1.1.0.  


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
	└── tests/ (optional)            (Unit tests or validation scripts)  


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



HISTORY
-----------------------------------

v1.0.0 - 08/25/2025
 	- First public release as an executable.
v1.1.0 - 10/28/2025
 	- Released as a Python Package on PyPi
 	- Added the ability to enter custom PDF equations.
v1.1.1 - 10/29/2025
 	-Fixed issues regarding faulty data upload.
v1.1.2 - 11/04/2025
	-Fixed bug involving unenforced mass limitations while viewing prior monochromatic spectra.
	-Tweaked plot labelling for increased legibility.
	-Added the ability to view previously generated histograms.


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