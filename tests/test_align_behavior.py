import numpy as np
import gammapbh as gp

def _first_mass(masses):
    if isinstance(masses, tuple) and len(masses)==2 and masses[0]:
        return float(masses[0][0]) if isinstance(masses[0][0], (int,float)) else float(masses[1][0])
    return float(masses[0]) if isinstance(masses, list) else pytest.skip("no mass found")

def test_align_modes_monotonic_and_sizes():
    masses = gp.discover_mass_folders()
    m = _first_mass(masses)

    E_sec, comps_sec = gp.load_spectra_components(m, align="secondary")
    E_pri, comps_pri = gp.load_spectra_components(m, align="primary")
    E_uni, comps_uni = gp.load_spectra_components(m, align="union")
    E_int, comps_int = gp.load_spectra_components(m, align="intersection")

    for E in (E_sec, E_pri, E_uni, E_int):
        e = np.asarray(E)
        assert e.ndim == 1 and e.size > 0
        assert np.all(np.diff(e) > 0), "E must be strictly increasing"

    assert E_uni.size >= max(E_sec.size, E_pri.size)
    assert E_int.size <= min(E_sec.size, E_pri.size)