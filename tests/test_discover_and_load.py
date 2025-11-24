import numpy as np

def test_discover_has_entries(data_root):
    import gammapbh as gp
    masses = gp.discover_mass_folders()     # should work without args
    assert masses, "discover_mass_folders() returned empty."

def test_load_components_shapes(sample_mass):
    """
    gp.load_spectra_components() must return (E, comps_dict).
    E must be 1D and each component array must match len(E).
    """
    import gammapbh as gp
    E, comps = gp.load_spectra_components(sample_mass)  # default align
    E = np.asarray(E)
    assert E.ndim == 1 and E.size > 0

    assert isinstance(comps, dict) and len(comps) > 0
    for k, v in comps.items():
        a = np.asarray(v)
        assert a.shape == E.shape, f"Component {k} length != E length"

def test_full_loader_energy_keys(sample_mass):
    """
    The raw loader can return a dict; if so, it should carry at least one energy axis key.
    """
    import gammapbh as gp
    raw = gp.load_spectra_components_full(sample_mass)
    if isinstance(raw, dict):
        keys = set(raw.keys())
        assert any(k in keys for k in ("E","energy","energies","E_MeV","energy_primary","energy_secondary"))