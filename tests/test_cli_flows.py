import os
import numpy as np
import pytest

@pytest.mark.slow
def test_view_previous_queue_overlay(cli_module, feed_inputs):
    """
    Create a fake saved GAUSSIAN run, then drive 'view_previous_spectra' to:
      - choose Gaussian category
      - queue the first run
      - plot all (0)
      - back out (b)
    """
    # Make a fake saved run with distributed_spectrum.txt and mass_distribution.txt
    run_dir = os.path.join(cli_module.GAUSS_RESULTS_DIR, "peak_3.00e+16_σ0.05_N50")
    os.makedirs(run_dir, exist_ok=True)
    # Energy grid & spectrum
    E = np.logspace(np.log10(1.0), np.log10(5000.0), 200)
    S = 1e-9 * (E**-1.2) * np.exp(-E/4000.0)
    np.savetxt(os.path.join(run_dir, "distributed_spectrum.txt"),
               np.column_stack([E, S]),
               header="E_gamma(MeV)   TotalSpectrum", fmt="%.6e")
    # Mass samples
    md = np.random.lognormal(mean=np.log(3e16)+0.05**2, sigma=0.05, size=50)
    np.savetxt(os.path.join(run_dir, "mass_distribution.txt"), md, header="Sampled masses (g)", fmt="%.8e")

    # Drive the menu: pick "Gaussian" (=2), queue first item ("1"), plot all ("0"), then back ("b")
    feeder = feed_inputs(["2", "1", "0", "b"])
    cli_module.view_previous_spectra()
    assert True

@pytest.mark.slow
def test_custom_equation_tool(cli_module, feed_inputs):
    # f(m) = 1/m over domain, small N, do not save
    feeder = feed_inputs(["1/m", "10", "n"])
    cli_module.custom_equation_pdf_tool()
    assert True

def test_view_previous_quit_text(cli_module, feed_inputs):
    feeder = feed_inputs(["q"])
    # Should return cleanly (no exception) when user types "q"
    cli_module.view_previous_spectra()

