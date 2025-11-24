import numpy as np
import gammapbh as gp

def _first_mass(masses):
    if isinstance(masses, tuple) and len(masses)==2 and masses[0]:
        return float(masses[0][0]) if isinstance(masses[0][0], (int,float)) else float(masses[1][0])
    return float(masses[0]) if isinstance(masses, list) else None

def _check_totals(E, comps):
    import numpy as np
    E = np.asarray(E)
    kset = set(comps.keys())
    assert "total" in kset, "normalized API should expose 'total' spectra"
    # If the breaking into prim/sec is available, total should be their sum
    if {"total_primary","total_secondary"} <= kset:
        assert np.allclose(comps["total"], comps["total_primary"] + comps["total_secondary"], rtol=1e-8, atol=1e-12)

def test_totals_secondary_align():
    masses = gp.discover_mass_folders()
    m = _first_mass(masses)
    assert m is not None
    E, comps = gp.load_spectra_components(m, align="secondary")
    _check_totals(E, comps)

def test_totals_primary_align():
    masses = gp.discover_mass_folders()
    m = _first_mass(masses)
    assert m is not None
    E, comps = gp.load_spectra_components(m, align="primary")
    _check_totals(E, comps)