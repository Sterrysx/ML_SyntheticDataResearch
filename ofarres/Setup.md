# Setup Guide - Synthetic Data Clustering Pipeline

**Complete installation and environment setup instructions**

---

## ⚠️ Prerequisites

Before starting, ensure you have:
1. **Conda** (Miniconda or Anaconda) installed
2. **R** (version 4.0 or higher)
3. **Python** 3.9 or higher (will be installed via conda)
4. **Git** (for cloning/updating the repository)

---

## 📋 Installation Steps

### Step 1: Create Conda Environment

**IMPORTANT**: You must create and activate the conda environment **BEFORE** running any installation scripts.

```bash
# Navigate to your repository root
cd /path/to/ML_SyntheticDataResearch

# Create a new conda environment with Python 3.13
conda create -n synthetic_data python=3.13

# Activate the environment
conda activate synthetic_data
```

**Verify activation:**
```bash
# You should see (synthetic_data) at the beginning of your terminal prompt
# Example: (synthetic_data) user@machine:~$

# Check Python version
python --version  # Should show Python 3.9.x
```

---

### Step 2: Navigate to Pipeline Directory

```bash
cd ofarres/
```

---

### Step 3: Run Installation Script

The `install_all.sh` script will:
1. Install all required R packages (jsonlite, mvtnorm, synthpop, dplyr)
2. Install all required Python packages (from `../requirements.txt`)
3. Verify the installation automatically

```bash
./installer_scripts/install_all.sh
```

**What happens during installation:**
- R packages are installed from CRAN
- Python packages are installed via pip
- Installation verification runs automatically
- Any errors will be reported immediately

**Expected output:**
```
==========================================
Installing Pipeline Dependencies
==========================================

✓ Conda environment active: synthetic_data

==========================================
Step 1: Installing R packages
==========================================
...
✅ All packages verified successfully!

==========================================
Step 2: Installing Python packages
==========================================
...
Successfully installed pandas-2.0.0 numpy-1.23.0 ...

==========================================
Step 3: Verifying installation
==========================================
...
✅ All dependencies verified!

==========================================
✅ Installation complete!
==========================================

You can now run the pipeline:
  ./run_all.sh
```

---

### Step 4: Verify Setup (Optional)

If you want to manually verify your setup later:

```bash
./installer_scripts/verify_setup.sh
```

This checks:
- R and Python installations
- All required R packages
- All required Python packages
- File permissions
- Configuration files

---

## 🚀 Running the Pipeline

After successful setup:

```bash
# Make sure you're in the ofarres/ directory
cd ofarres/

# Ensure conda environment is activated
conda activate synthetic_data

# Run the complete pipeline (~15-20 minutes)
./run_all.sh
```

---

## 🔧 Manual Installation (If Automated Script Fails)

### Install R Packages Manually

```bash
Rscript installer_scripts/install_requirements.R
```

Or in R console:
```r
install.packages(c("jsonlite", "mvtnorm", "synthpop", "dplyr"))
```

### Install Python Packages Manually

```bash
# Make sure conda environment is activated
conda activate synthetic_data

# Install from requirements file
pip install -r requirements.txt
```

---

## 🐛 Troubleshooting

### Issue: "conda: command not found"

**Solution**: Install Miniconda first
```bash
# Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run installer
bash Miniconda3-latest-Linux-x86_64.sh

# Restart terminal or source .bashrc
source ~/.bashrc
```

---

### Issue: "R is not installed"

**Solution**: Install R

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install r-base r-base-dev
```

**macOS (with Homebrew):**
```bash
brew install r
```

**Windows:**
- Download from https://cran.r-project.org/
- Install using the GUI installer

---

### Issue: "ERROR: No conda environment is active!"

**Solution**: Activate the environment before running install_all.sh
```bash
conda activate synthetic_data
./scripts/install_all.sh
```

---

### Issue: R package installation fails with "permission denied"

**Solution**: Don't use sudo with R package installation. Packages will install to your user library.

If prompted about creating a personal library, select **Yes**.

---

### Issue: Python package installation fails

**Solution**: Ensure pip is up to date
```bash
conda activate synthetic_data
pip install --upgrade pip
pip install -r ../requirements.txt
```

---

### Issue: "verify_setup.sh: Permission denied"

**Solution**: Make scripts executable
```bash
chmod +x installer_scripts/install_all.sh
chmod +x installer_scripts/verify_setup.sh
chmod +x run_all.sh
```

---

## 📦 What Gets Installed?

### R Packages
- **jsonlite** (1.8.0+) - JSON configuration parsing
- **mvtnorm** (1.1.0+) - Multivariate normal distribution generation
- **synthpop** (1.7.0+) - Synthetic data generation via CART
- **dplyr** (1.0.0+) - Data manipulation

### Python Packages (Key Dependencies)
- **pandas** (2.0.0+) - Data structures and analysis
- **numpy** (1.23.0+) - Numerical computing
- **scikit-learn** (1.2.0+) - Machine learning algorithms
- **matplotlib** (3.6.0+) - Plotting and visualization
- **seaborn** (0.12.0+) - Statistical visualization
- **scipy** (1.10.0+) - Scientific computing
- **jupyter** (1.0.0+) - Interactive notebooks
- **tqdm** (4.64.0+) - Progress bars

*See `../requirements.txt` for complete list*

---

## 🔄 Updating the Environment

If you need to update packages:

```bash
# Activate environment
conda activate synthetic_data

# Update Python packages
pip install --upgrade -r ../requirements.txt

# Update R packages
Rscript -e "update.packages(ask = FALSE)"
```

---

## 🗑️ Removing the Environment

If you need to start fresh:

```bash
# Deactivate if currently active
conda deactivate

# Remove environment
conda env remove -n synthetic_data

# Recreate from scratch
conda create -n synthetic_data python=3.9
conda activate synthetic_data
cd ofarres/
./scripts/install_all.sh
```

---

## ✅ Post-Installation Checklist

Before running the pipeline, verify:

- [ ] Conda environment `synthetic_data` is created and activated
- [ ] `installer_scripts/install_all.sh` completed without errors
- [ ] `installer_scripts/verify_setup.sh` shows all ✓ checks passed
- [ ] You are in the `ofarres/` directory
- [ ] `run_all.sh` has execute permissions (`chmod +x run_all.sh`)
- [ ] `config/config.json` exists in the `ofarres/config/` directory

---

## 📚 Next Steps

After successful setup:

1. **Review configuration**: Check `config/config.json` for simulation parameters
2. **Read pipeline documentation**: See [README.md](README.md) for details
3. **Run the pipeline**: Execute `./run_all.sh`
4. **Analyze results**: Outputs will be in `03_analysis/`

---

## 💡 Tips

- **Always activate the conda environment** before working with the pipeline
- **Don't use sudo** for R or Python package installations
- **Keep the environment clean** - only install required packages
- **Check verify_setup.sh** if anything behaves unexpectedly
- **Read error messages** - they usually indicate exactly what's missing

---

## 📞 Getting Help

If you encounter issues:

1. Run `./installer_scripts/verify_setup.sh` to diagnose the problem
2. Check the **Troubleshooting** section above
3. Review error messages carefully
4. Ensure conda environment is activated
5. Try manual installation steps

---

**Last Updated**: January 3, 2026  
**Conda Environment**: `synthetic_data`  
**Python Version**: 3.9+  
**R Version**: 4.0+
