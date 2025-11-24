import pytest, pathlib
import gammapbh as gp

def _first_mass_value(masses):
    """
    Normalize discover_mass_folders() output (tuple or list) to one float mass.
    """
    if isinstance(masses, tuple) and len(masses) == 2:
        a, b = masses
        # prefer numeric list if present
        if a and isinstance(a[0], (int, float)):
            return float(a[0])
        if b and isinstance(b[0], (int, float)):
            return float(b[0])
        # else parse the first string
        for lst in (a, b):
            if lst and isinstance(lst[0], str):
                try:
                    return float(lst[0])
                except Exception:
                    pass
    if isinstance(masses, list) and masses:
        try:
            return float(masses[0])
        except Exception:
            pass
    pytest.skip("Could not normalize a candidate mass from discover_mass_folders().")

@pytest.fixture(scope="session")
def sample_mass():
    masses = gp.discover_mass_folders()
    return _first_mass_value(masses)

@pytest.fixture(scope="session")
def data_root():
    root = pathlib.Path(gp.__file__).parent / "blackhawk_data"
    if not root.exists():
        pytest.skip("Bundled blackhawk_data not found; skipping data-driven tests.")
    return root