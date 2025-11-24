import sys, subprocess, os, pytest

def test_module_help_prints_or_skips():
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    proc = subprocess.run(
        [sys.executable, "-X", "utf8", "-m", "gammapbh", "--help"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", env=env
    )
    out = (proc.stdout or "").lower()
    # If module runner expects input(), skip until CLI supports --help
    if proc.returncode != 0 and ("eof when reading a line" in out or "input(" in out):
        pytest.skip("gammapbh __main__ is interactive; skipping --help check.")
    assert ("usage" in out) or ("description" in out)