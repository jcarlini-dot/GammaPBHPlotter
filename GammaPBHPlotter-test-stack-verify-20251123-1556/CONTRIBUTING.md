## Contributing to GammaPBHPlotter



Thanks for considering a contribution! This document explains how to:

- contribute code or docs,

- report bugs and request features,

- run tests locally.



> By contributing, you agree your contributions will be licensed under the repository’s license (see `LICENSE`, GNU GPL-3.0-or-later).



## Code of Conduct

This project follows the Contributor Covenant. See `CODE\_OF\_CONDUCT.md`.



## How to Get Help / Ask Questions

- Usage questions: open a \*\*GitHub Discussion\*\* (if enabled) or a \*\*question\*\*-labeled issue.

- Bug reports: open a \*\*Bug report\*\* issue (template provided).

- Security issues: email the maintainers (see `SECURITY.md` if present) rather than filing a public issue.



## Reporting Bugs

Please include:

- steps to reproduce,

- expected vs. actual behavior,

- your OS/Python versions,

- `pip show gammapbh`, and

- minimal data / command sequence (when possible).



## Requesting Features

Open a \*\*Feature request\*\* issue and explain:

- the problem you’re trying to solve,

- proposed changes, and

- any alternatives you considered.



## Development Setup

```bash

git clone https://github.com/jcarlini-dot/GammaPBHPlotter.git

cd GammaPBHPlotter

python -m venv .venv \&\& . .venv/Scripts/activate  # Windows PowerShell

pip install -e .

pip install pytest


