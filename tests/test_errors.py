import pytest, gammapbh as gp

def test_invalid_mass_raises():
    with pytest.raises(Exception):
        gp.load_spectra_components(1.23456789e42)  # not in grid