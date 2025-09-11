## Modeling the emission components

The gamma-ray spectrum of a PBH within the relevant mass range consists
of four primary components; direct/primary Hawking radiation, secondary
radiation, final-state radiation, and in-flight annihilation.

Direct Hawking radiation accounts for all kinematically allowed
elementary particles formed at the event horizon, including gamma-ray
photons. Secondary radiation originates from the decay of unstable
particles and contributes significantly at lower energies. We rely on
`BlackHawk` to evaluate the gamma-ray primary and secondary spectral
components. `BlackHawk` uses `PYTHIA` for the modeling of the
hadronization and decay processes leading to the secondary spectra.
Final-state radiation originates from relativistic electrons and
positrons and has a differential spectrum given by Eq. [1](#eq:FSRRate):

<a id="eq:FSRRate"></a>

$$
\frac{dN_{\gamma}^{\mathrm{FSR}}}{dE_{\gamma}}
=\frac{\alpha}{2\pi}\int dE_{e^{+}}\,
\frac{dN_{e^{+}}}{dE_{e^{+}}}
\left(\frac{2}{E_{\gamma}}+\frac{E_{\gamma}}{E_{e^{+}}^{2}}-\frac{2}{E_{e^{+}}}\right)
\left[\ln\!\left(\frac{2E_{e^{+}}+(E_{e^{+}}-E_{\gamma})}{m_{e}^{2}}\right)-1\right].
$$

where $\alpha=137.037$ is the fine-structure constant, $E_{e^{+}}$
is the kinetic energy of a given positron ($e^{+}$), $E_{\gamma}$ is
the energy of the emitted photon, $m_{e}=0.511\ \mathrm{MeV}$ is the rest
mass of the electron, and $\frac{dN_{e^{+}}}{dE_{e^{+}}}$ is the
differential spectrum of emitted electrons/positrons. In addition to the
previously mentioned components, gamma rays can be produced through
pair-annihilation of positrons with interstellar-medium electrons. This
is known as in-flight annihilation and its differential spectrum is given
by Eq. [2](#eq:IARate):

<a id="eq:IARate"></a>

$$
\frac{dN_{\gamma}^{\mathrm{IA}}}{dE_{\gamma}}
=\frac{\pi\alpha^{2}n_{H}}{m_{e}}
\int_{m_{e}}^{\infty} dE_{e^{+}}\,
\frac{dN_{e^{+}}}{dE_{e^{+}}}
\int_{E_{\min}}^{E_{e^{+}}}
\frac{dE}{dE/dx}\,
\frac{P_{E_{e^{+}}\to E}}{E^{2}-m_{e}^{2}}
\left(
-2-\frac{(E+m_{e})\!\left[m_{e}^{2}(E+m_{e})+E_{\gamma}^{2}(E+3m_{e})-E_{\gamma}(E+m_{e})(E+3m_{e})\right]}
{E_{\gamma}^{2}(E-E_{\gamma}+m_{e})^{2}}
\right).
$$

We take $n_{H}=1\ \mathrm{cm}^{-3}$ as the density of interstellar
medium hydrogen (and by extension electrons). $E_{e^{+}}$ is again the
initial positron total energy, $E$ is the final positron total energy,
$dE/dx$ is the rate of positron energy lost per path via the
Bethe–Bloch formula, $E_{\gamma}$ is the resulting photon energy from
annihilation, and $P_{E_{e^{+}}\to E}$ is the probability of a
particular positron of a given initial and final energy to decay. This
probability matrix can be calculated as Eq. [3](#eq:Pmatrix):

<a id="eq:Pmatrix"></a>

$$
P_{E_{e^{+}}\to E}
=\exp\!\Biggl(
-\,n_{H}\int_{E}^{E_{e^{+}}}\sigma_{\mathrm{ann}}(E')\,\frac{dE'}{dx}\,dE'
\Biggr),
$$

where $\sigma_{\mathrm{ann}}$ is the annihilation cross section for
positrons of a given energy.

In Fig. 1, we give the individual gamma-ray spectral components as well
as their sum for a PBH of mass $3\times10^{15}$ grams.

![Monochromatic spectrum](figures/monochromatic.png)

**Figure 1.** The total gamma-ray spectrum of a $3\times10^{15}$ grams PBH as well as its components.

## PBH Mass Distribution

Users can calculate the gamma-ray spectra from four types of PBH mass
distributions. Those are, i) a monochromatic distribution with a mass to
be set in the range of $5\times10^{13}$ to $1\times10^{19}$ grams,
ii) a Gaussian distribution of PBH masses originating from a Gaussian
distribution of density perturbations, iii) a more realistic
non-Gaussian PBH mass distribution, and iv) a log-normal
distribution of PBH masses. In Fig. 2, we give the gamma-ray spectra
from monochromatic and Gaussian PBH mass distributions.

![Spectrum comparison](figures/spectrum_comparison.png)

**Figure 2.** The total gamma-ray spectrum per PBH, from a PBH of mass $3\times10^{15}$ grams (blue line) and from a Gaussian distribution of density perturbations leading to a distribution with mean mass $3\times10^{15}$ grams. $\sigma$ refers to the standard deviation of the initial density perturbations.
