# Project Setup & Environment Guide

This project uses a specific **Conda** environment named `synthetic_data` running **Python 3.13** on WSL/Linux.

## 1\. Prerequisites: Installing Conda

If you do not have Conda installed (command `conda` not found), you need to install it first. We recommend **Miniconda** (a lightweight version) for WSL/Linux.

**Run these commands in your terminal:**

```bash
# 1. Download the installer
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh

# 2. Run the install script
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3

# 3. Initialize Conda (and restart your shell afterwards)
~/miniconda3/bin/conda init bash
source ~/.bashrc
```

## 2\. First-Time Installation

If you are setting this up for the first time (or on a new machine), follow these steps to create the environment and install libraries:

### Step A: Create the Environment

We use `conda-forge` to ensure the latest Python compatibility.

```bash
conda create -n synthetic_data python=3.13 -c conda-forge
```

### Step B: Activate

```bash
conda activate synthetic_data
```

### Step C: Install Dependencies

This project uses `pip` for package management within the Conda environment.

```bash
pip install -r requirements.txt
```

## 3\. Daily Usage (Quick Start)

Once the environment is set up, you only need to run this when you start working:

```bash
conda activate synthetic_data
```

To verify you are in the correct environment:

```bash
python --version
# Should output Python 3.13.x
```

## 4\. Maintenance: Adding New Libraries

**Crucial Step:**
If you install a **new** package during development (e.g., `pip install matplotlib`), you **must** update the requirements file immediately so others (or your future self) have the exact same setup.

Run this command after every new installation:

```bash
pip list --format=freeze > requirements.txt
```

## 5\. Deactivation

When you are finished working:

```bash
conda deactivate
```