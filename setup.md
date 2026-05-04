# Setup Guide

This guide covers the shared environment and dependencies for both pipelines. Each pipeline also has a local setup guide:

- clustering/Setup.md
- regression/Setup.md

## Prerequisites

- uv (recommended) or conda
- R 4.0 or newer
- Git
- A C/C++ build toolchain for R packages

## System dependencies (new machine)

Install system tools before creating the Python environment.

Linux (Debian/Ubuntu):

```bash
sudo apt update
sudo apt install -y build-essential curl git r-base r-base-dev \
	libssl-dev libcurl4-openssl-dev libxml2-dev
```

Linux (Fedora/RHEL):

```bash
sudo dnf install -y gcc gcc-c++ make curl git R \
	openssl-devel libcurl-devel libxml2-devel
```

macOS:

```bash
xcode-select --install
brew install r git
```

Windows:

- Install R 4.0+ from https://cran.r-project.org/bin/windows/base/
- Install Rtools (needed to compile some R packages):
	https://cran.r-project.org/bin/windows/Rtools/
- Install Git from https://git-scm.com/download/win

If pyarrow fails to install, install system Arrow dependencies for your OS and retry.

## 1. Create and activate a uv environment (recommended)

Install uv if needed:

Linux/macOS:

```bash
curl -Ls https://astral.sh/uv/install.sh | sh
```

Windows PowerShell:

```powershell
irm https://astral.sh/uv/install.ps1 | iex
```

Create and activate the environment, then install the root requirements:

```bash
uv venv --python 3.13
source .venv/bin/activate
uv pip install -r requirements.txt
```

## 2. Alternative: conda environment

```bash
conda create -n synthetic_data python=3.13
conda activate synthetic_data
pip install --upgrade pip
pip install -r requirements.txt
```

## 3. Install R packages

The R scripts will install missing packages automatically on first run. If you prefer to preinstall them, run:

```bash
Rscript -e "install.packages(c('jsonlite','mvtnorm','synthpop','arrow','MultiDiscreteRNG'), repos='https://cloud.r-project.org')"
```

Notes:
- clustering uses jsonlite, mvtnorm, synthpop, arrow
- regression uses jsonlite, mvtnorm, MultiDiscreteRNG, synthpop, arrow

## 4. Make scripts executable (if needed)

```bash
chmod +x clustering/run_all.sh
chmod +x regression/run_all.sh
```

## 5. Smoke test (recommended)

Clustering quick test:
1. Edit clustering/config/config.json
2. Reduce the grid to one value each (p, k, separation, rho, distribution)
3. Set simulation.n to 1 and simulation.m to 2
4. Run:

```bash
cd clustering
./run_all.sh
```

Regression quick test:
1. Edit regression/config/config.json
2. Set simulation.M to 10 (and optionally reduce N and p lists)
3. Run:

```bash
cd regression
./run_all.sh
```

Restore the original values after testing.

## 6. Verify the Python environment (optional)

```bash
python - <<'PY'
import numpy, pandas, pyarrow, sklearn, scipy, statsmodels
print("OK", numpy.__version__, pandas.__version__)
PY
```

If this fails, re-run the environment setup and check the troubleshooting section.

## Troubleshooting

- uv not found: install uv (see step 1) and restart your shell.
- no active venv: run source .venv/bin/activate before ./run_all.sh.
- conda not found: ensure conda is installed and on PATH, then restart the terminal.
- Rscript not found: install R 4.0+ and verify Rscript is available on PATH.
- Python import errors: confirm the environment is active and requirements.txt is installed.
- arrow install issues: install system dependencies for Arrow and retry pip install.
