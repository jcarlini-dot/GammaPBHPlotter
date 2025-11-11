import numpy as np
import math

def test_discover_and_requirements(cli_module):
    masses, names = cli_module.discover_mass_folders(cli_module.DATA_DIR)
    assert len(masses) == 3
    # Sorted ascending
    assert masses == sorted(masses)
    # All required files present
    for folder in names:
        for f in cli_module.REQUIRED_FILES:
            p = cli_module.os.path.join(cli_module.DATA_DIR, folder, f)
            assert cli_module.os.path.isfile(p)

def snap_to_available(x, masses, rel_tol=1e-12):
    """
    Return the exact grid mass if x is 'close enough' to one of the available masses,
    otherwise None.

    Parameters
    ----------
    x : float
        The query mass (g).
    masses : Sequence[float]
        Available grid masses (g), strictly positive.
    rel_tol : float
        Relative tolerance for snapping (default 1e-12).

    Returns
    -------
    float | None
        The snapped mass from `masses` if within tolerance, else None.
    """
    for m in masses:
        if math.isclose(x, m, rel_tol=rel_tol):
            return m
    # Optional: also allow choosing the nearest if it's the best match within tol
    # nearest = min(masses, key=lambda v: abs(v - x))
    # if math.isclose(x, nearest, rel_tol=rel_tol):
    #     return nearest
    return None


def test_parse_float_list_verbose(cli_module, capsys):
    s = "1e16, 2.5e16, abc, 1e16"
    vals = cli_module.parse_float_list_verbose(s, name="mass", bounds=(1e16, 1e17))
    # dedup, preserve first
    assert vals[:2] == [1e16, 2.5e16]
    # a warning about abc and duplicate printed
    out = capsys.readouterr().out
    assert "Skipping token 'abc'" in out
    assert "Duplicate mass 1e+16" in out

def test_load_xy_lenient_skips_singleton(cli_module, tmp_path):
    p = tmp_path / "foo.txt"
    with open(p, "w") as f:
        f.write("# header\n")
        f.write("777\n")
        for i in range(3):
            f.write(f"{i+1} {i+2} {999}\n")
    arr = cli_module.load_xy_lenient(str(p))
    # returns first two columns only, skips the lone "777" row
    assert arr.shape == (3, 2)
    assert (arr[:, 0] == np.array([1, 2, 3])).all()
    assert (arr[:, 1] == np.array([2, 3, 4])).all()

def test_delta_l_monotonic(cli_module):
    r = np.logspace(-3, 0, 50)
    dl = cli_module.delta_l(r, kappa=3.3, delta_c=0.59, gamma=0.36)
    # finite and bounded by mapping safety clip
    assert np.isfinite(dl).all()
    assert (dl <= (8.0/6.0)).all()
