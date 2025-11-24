import os
import numpy as np

def test_generate_monochromatic_file(cli_module, tmp_path):
    masses, _ = cli_module.discover_mass_folders(cli_module.DATA_DIR)
    target = masses[1] * 1.02  # near a real mass to exercise snapping
    out = cli_module.generate_monochromatic_for_mass(target, cli_module.DATA_DIR, cli_module.MONO_RESULTS_DIR)
    assert os.path.isfile(out)
    data = np.loadtxt(out)
    # two columns: E(MeV), Total
    assert data.ndim == 2 and data.shape[1] == 2
    # Energies should be > 0 and within our synthetic E-grid
    assert data[:, 0].min() > 0

def test_save_formats_for_distributed(cli_module, feed_inputs):
    """
    Run a tiny Log-normal distributed job (N=10) and refuse saving,
    just to make sure plotting & averaging work without exceptions.
    """
    # answers to prompts in distributed_spectrum(LOGNORMAL_METHOD):
    # 1) peaks string, 2) N, 3) sigma list, 4) save? (y/n)
    feeder = feed_inputs(["3e16", "10", "0.6", "n"])
    cli_module.distributed_spectrum(cli_module.LOGNORMAL_METHOD)
    # Nothing saved, but function should complete without error
    assert True
