import sys, subprocess, os

def test_module_help_prints_and_exits():
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    proc = subprocess.run(
        [sys.executable, "-X", "utf8", "-m", "gammapbh", "--help"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", env=env
    )
    assert proc.returncode == 0
    out = (proc.stdout or "").lower()
    assert ("usage" in out) or ("description" in out)