#!/bin/bash

# ==============================================================================
# MASTER EXECUTION SCRIPT: Synthetic Data Regression Research Pipeline
# ==============================================================================
# Purpose: Execute the complete regression simulation pipeline (linear + logistic)
# Usage:   cd regression && ./run_all.sh
# ==============================================================================

set -e  # Exit on any error

# ==============================================================================
# ZOMBIE PROCESS CLEANUP TRAP
# Automatically hunts down and kills orphaned workers if you hit Ctrl+C
# ==============================================================================
cleanup() {
    pkill -f "sd_worker.R" 2>/dev/null || true
    pkill -f "od_linear_worker.R" 2>/dev/null || true
    pkill -f "od_logistic_worker.R" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

# Environment detection: prefer an active venv (uv), otherwise fall back to conda
if [ -n "${VIRTUAL_ENV:-}" ]; then
    echo "✓ Using active virtual environment: ${VIRTUAL_ENV}"
elif command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    if conda env list | grep -q "synthetic_data"; then
        conda activate synthetic_data
        echo "✓ Activated conda environment: synthetic_data"
    else
        echo "❌ ERROR: conda environment 'synthetic_data' not found"
        echo "Create it or use uv (from repo root):"
        echo "  uv venv --python 3.13"
        echo "  source .venv/bin/activate"
        echo "  uv pip install -r requirements.txt"
        exit 1
    fi
else
    echo "❌ ERROR: no active virtual environment and conda is not installed"
    echo "Activate a venv (uv) or install conda."
    exit 1
fi

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1"; }

record_timing_history() {
    mkdir -p results
    TIMING_CSV="results/pipeline_timing_history.csv"
    TIMING_JSON="results/latest_pipeline_timing.json"

    RUN_TIMESTAMP_BARCELONA="$(TZ=Europe/Madrid date +%Y-%m-%d-%H-%M-%S)" \
    CONFIG_HASH="$CONFIG_HASH" \
    CONFIG_CHANGED="$CONFIG_CHANGED" \
    SIM_M="$SIM_M" \
    NUM_N="$NUM_N" \
    NUM_P="$NUM_P" \
    NUM_RHO="$NUM_RHO" \
    NUM_SIGMA2="$NUM_SIGMA2" \
    NUM_P1="$NUM_P1" \
    OD_LINEAR_CONT_SCENARIOS="$OD_LINEAR_CONT_SCENARIOS" \
    OD_LINEAR_BIN_SCENARIOS="$OD_LINEAR_BIN_SCENARIOS" \
    OD_LOGISTIC_CONT_SCENARIOS="$OD_LOGISTIC_CONT_SCENARIOS" \
    OD_LOGISTIC_BIN_SCENARIOS="$OD_LOGISTIC_BIN_SCENARIOS" \
    OD_LINEAR_CONT_TARGET_ROWS="$OD_LINEAR_CONT_TARGET_ROWS" \
    OD_LINEAR_BIN_TARGET_ROWS="$OD_LINEAR_BIN_TARGET_ROWS" \
    OD_LOGISTIC_CONT_TARGET_ROWS="$OD_LOGISTIC_CONT_TARGET_ROWS" \
    OD_LOGISTIC_BIN_TARGET_ROWS="$OD_LOGISTIC_BIN_TARGET_ROWS" \
    SD_LINEAR_EXPECTED_FILES="$SD_LINEAR_EXPECTED_FILES" \
    SD_LOGISTIC_EXPECTED_FILES="$SD_LOGISTIC_EXPECTED_FILES" \
    SD_LINEAR_TARGET_ROWS="$SD_LINEAR_TARGET_ROWS" \
    SD_LOGISTIC_TARGET_ROWS="$SD_LOGISTIC_TARGET_ROWS" \
    LR_LINEAR_EXPECTED_MODELS="$LR_LINEAR_EXPECTED_MODELS" \
    LR_LOGISTIC_EXPECTED_MODELS="$LR_LOGISTIC_EXPECTED_MODELS" \
    STEP0_DURATION_S="$STEP0_DURATION_S" \
    OD_LINEAR_CONT_DURATION_S="$OD_LINEAR_CONT_DURATION_S" \
    OD_LINEAR_BIN_DURATION_S="$OD_LINEAR_BIN_DURATION_S" \
    OD_LINEAR_TOTAL_DURATION_S="$OD_LINEAR_TOTAL_DURATION_S" \
    OD_LOGISTIC_CONT_DURATION_S="$OD_LOGISTIC_CONT_DURATION_S" \
    OD_LOGISTIC_BIN_DURATION_S="$OD_LOGISTIC_BIN_DURATION_S" \
    OD_LOGISTIC_TOTAL_DURATION_S="$OD_LOGISTIC_TOTAL_DURATION_S" \
    SD_DURATION_S="$SD_DURATION_S" \
    SD_LINEAR_CONT_DURATION_S="$SD_LINEAR_CONT_DURATION_S" \
    SD_LINEAR_BIN_DURATION_S="$SD_LINEAR_BIN_DURATION_S" \
    SD_LOGISTIC_CONT_DURATION_S="$SD_LOGISTIC_CONT_DURATION_S" \
    SD_LOGISTIC_BIN_DURATION_S="$SD_LOGISTIC_BIN_DURATION_S" \
    LR_LINEAR_DURATION_S="$LR_LINEAR_DURATION_S" \
    LR_LOGISTIC_DURATION_S="$LR_LOGISTIC_DURATION_S" \
    PIPELINE_TOTAL_DURATION_S="$TOTAL_RUNTIME" \
    OD_LINEAR_SKIPPED="$OD_LINEAR_SKIPPED" \
    OD_LOGISTIC_SKIPPED="$OD_LOGISTIC_SKIPPED" \
    SD_SKIPPED="$SD_SKIPPED" \
    LR_LINEAR_SKIPPED="$LR_LINEAR_SKIPPED" \
    LR_LOGISTIC_SKIPPED="$LR_LOGISTIC_SKIPPED" \
    TIMING_CSV="$TIMING_CSV" \
    TIMING_JSON="$TIMING_JSON" \
    $PYTHON_CMD - <<'PY'
import csv
import json
import os

fields = [
    "run_timestamp_barcelona",
    "config_hash",
    "config_changed",
    "sim_m",
    "num_n",
    "num_p",
    "num_rho",
    "num_sigma2",
    "num_p1",
    "od_linear_cont_scenarios",
    "od_linear_bin_scenarios",
    "od_logistic_cont_scenarios",
    "od_logistic_bin_scenarios",
    "od_linear_cont_target_rows",
    "od_linear_bin_target_rows",
    "od_logistic_cont_target_rows",
    "od_logistic_bin_target_rows",
    "sd_linear_expected_files",
    "sd_logistic_expected_files",
    "sd_linear_target_rows",
    "sd_logistic_target_rows",
    "lr_linear_expected_models",
    "lr_logistic_expected_models",
    "step0_duration_s",
    "od_linear_cont_duration_s",
    "od_linear_bin_duration_s",
    "od_linear_total_duration_s",
    "od_logistic_cont_duration_s",
    "od_logistic_bin_duration_s",
    "od_logistic_total_duration_s",
    "sd_duration_s",
    "sd_linear_cont_duration_s",
    "sd_linear_bin_duration_s",
    "sd_logistic_cont_duration_s",
    "sd_logistic_bin_duration_s",
    "lr_linear_duration_s",
    "lr_logistic_duration_s",
    "pipeline_total_duration_s",
    "od_linear_skipped",
    "od_logistic_skipped",
    "sd_skipped",
    "lr_linear_skipped",
    "lr_logistic_skipped",
]

row = {
    "run_timestamp_barcelona": os.environ["RUN_TIMESTAMP_BARCELONA"],
    "config_hash": os.environ["CONFIG_HASH"],
    "config_changed": os.environ["CONFIG_CHANGED"],
    "sim_m": os.environ["SIM_M"],
    "num_n": os.environ["NUM_N"],
    "num_p": os.environ["NUM_P"],
    "num_rho": os.environ["NUM_RHO"],
    "num_sigma2": os.environ["NUM_SIGMA2"],
    "num_p1": os.environ["NUM_P1"],
    "od_linear_cont_scenarios": os.environ["OD_LINEAR_CONT_SCENARIOS"],
    "od_linear_bin_scenarios": os.environ["OD_LINEAR_BIN_SCENARIOS"],
    "od_logistic_cont_scenarios": os.environ["OD_LOGISTIC_CONT_SCENARIOS"],
    "od_logistic_bin_scenarios": os.environ["OD_LOGISTIC_BIN_SCENARIOS"],
    "od_linear_cont_target_rows": os.environ["OD_LINEAR_CONT_TARGET_ROWS"],
    "od_linear_bin_target_rows": os.environ["OD_LINEAR_BIN_TARGET_ROWS"],
    "od_logistic_cont_target_rows": os.environ["OD_LOGISTIC_CONT_TARGET_ROWS"],
    "od_logistic_bin_target_rows": os.environ["OD_LOGISTIC_BIN_TARGET_ROWS"],
    "sd_linear_expected_files": os.environ["SD_LINEAR_EXPECTED_FILES"],
    "sd_logistic_expected_files": os.environ["SD_LOGISTIC_EXPECTED_FILES"],
    "sd_linear_target_rows": os.environ["SD_LINEAR_TARGET_ROWS"],
    "sd_logistic_target_rows": os.environ["SD_LOGISTIC_TARGET_ROWS"],
    "lr_linear_expected_models": os.environ["LR_LINEAR_EXPECTED_MODELS"],
    "lr_logistic_expected_models": os.environ["LR_LOGISTIC_EXPECTED_MODELS"],
    "step0_duration_s": os.environ["STEP0_DURATION_S"],
    "od_linear_cont_duration_s": os.environ["OD_LINEAR_CONT_DURATION_S"],
    "od_linear_bin_duration_s": os.environ["OD_LINEAR_BIN_DURATION_S"],
    "od_linear_total_duration_s": os.environ["OD_LINEAR_TOTAL_DURATION_S"],
    "od_logistic_cont_duration_s": os.environ["OD_LOGISTIC_CONT_DURATION_S"],
    "od_logistic_bin_duration_s": os.environ["OD_LOGISTIC_BIN_DURATION_S"],
    "od_logistic_total_duration_s": os.environ["OD_LOGISTIC_TOTAL_DURATION_S"],
    "sd_duration_s": os.environ["SD_DURATION_S"],
    "sd_linear_cont_duration_s": os.environ["SD_LINEAR_CONT_DURATION_S"],
    "sd_linear_bin_duration_s": os.environ["SD_LINEAR_BIN_DURATION_S"],
    "sd_logistic_cont_duration_s": os.environ["SD_LOGISTIC_CONT_DURATION_S"],
    "sd_logistic_bin_duration_s": os.environ["SD_LOGISTIC_BIN_DURATION_S"],
    "lr_linear_duration_s": os.environ["LR_LINEAR_DURATION_S"],
    "lr_logistic_duration_s": os.environ["LR_LOGISTIC_DURATION_S"],
    "pipeline_total_duration_s": os.environ["PIPELINE_TOTAL_DURATION_S"],
    "od_linear_skipped": os.environ["OD_LINEAR_SKIPPED"],
    "od_logistic_skipped": os.environ["OD_LOGISTIC_SKIPPED"],
    "sd_skipped": os.environ["SD_SKIPPED"],
    "lr_linear_skipped": os.environ["LR_LINEAR_SKIPPED"],
    "lr_logistic_skipped": os.environ["LR_LOGISTIC_SKIPPED"],
}

csv_path = os.environ["TIMING_CSV"]
json_path = os.environ["TIMING_JSON"]
csv_exists = os.path.exists(csv_path)

with open(csv_path, "a", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    if not csv_exists:
        writer.writeheader()
    writer.writerow(row)

with open(json_path, "w", encoding="utf-8") as handle:
    json.dump(row, handle, indent=2)
PY
}

# Banner
echo "=============================================================================="
echo "        SYNTHETIC DATA REGRESSION RESEARCH PIPELINE (LINEAR + LOGISTIC)      "
echo "=============================================================================="
echo ""

PIPELINE_START_TIME=$(date +%s)
log_info "Starting full regression pipeline..."
echo ""

STEP0_DURATION_S=0
OD_LINEAR_CONT_DURATION_S=0
OD_LINEAR_BIN_DURATION_S=0
OD_LINEAR_TOTAL_DURATION_S=0
OD_LOGISTIC_CONT_DURATION_S=0
OD_LOGISTIC_BIN_DURATION_S=0
OD_LOGISTIC_TOTAL_DURATION_S=0
SD_DURATION_S=0
SD_LINEAR_CONT_DURATION_S=0
SD_LINEAR_BIN_DURATION_S=0
SD_LOGISTIC_CONT_DURATION_S=0
SD_LOGISTIC_BIN_DURATION_S=0
LR_LINEAR_DURATION_S=0
LR_LOGISTIC_DURATION_S=0

OD_LINEAR_SKIPPED=0
OD_LOGISTIC_SKIPPED=0
SD_SKIPPED=0
LR_LINEAR_SKIPPED=0
LR_LOGISTIC_SKIPPED=0

# ==============================================================================
# STEP 0: Validate Environment & SMART CACHE MANAGEMENT
# ==============================================================================
log_info "STEP 0: Validating environment and checking config consistency..."
STEP0_START_TIME=$(date +%s)

if [ ! -f "config/config.json" ]; then
    log_error "config/config.json not found!"
    exit 1
fi
log_success "config/config.json found"

# Compute hash of current config.json
CONFIG_HASH=$(md5sum config/config.json | awk '{print $1}')
HASH_FILE="data/.config_hash"
CONFIG_CHANGED=false

if [ -f "$HASH_FILE" ]; then
    PREV_HASH=$(cat "$HASH_FILE")
    if [ "$CONFIG_HASH" != "$PREV_HASH" ]; then
        CONFIG_CHANGED=true
        log_warning "config.json has CHANGED since last run. Nuking all previous data..."
        rm -f data/original/*.parquet 2>/dev/null || true
        rm -f data/synthetic/*.parquet 2>/dev/null || true
        rm -f data/*_packed.parquet 2>/dev/null || true
        rm -f results/aggregated_*.parquet 2>/dev/null || true
        rm -f results/summary_table.csv 2>/dev/null || true
        find results/figures -type f -delete 2>/dev/null || true
        find results/figures_final -type f -delete 2>/dev/null || true
        log_success "Environment cleaned (config changed)."
    else
        log_success "config.json unchanged — will resume from where we left off."
    fi
else
    # No hash file yet — this is either a true first run or a pre-existing run
    # that never had hashing. NEVER nuke existing data; just save the hash.
    log_info "No previous config hash found. Saving current config fingerprint."
    log_success "Existing data (if any) will be kept — resuming."
fi

# Save the current config hash for future runs
mkdir -p data
echo "$CONFIG_HASH" > "$HASH_FILE"

# Check Python
if command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
elif command -v python &> /dev/null; then
    PYTHON_CMD="python"
else
    log_error "Python is not installed or not in PATH"
    exit 1
fi

# Pre-compute expected counts from config
SIM_M_TYPE=$($PYTHON_CMD - <<'PY'
import json
c = json.load(open('config/config.json'))
print(type(c['simulation']['M']).__name__)
PY
)
if [ "$SIM_M_TYPE" != "int" ] && [ "$SIM_M_TYPE" != "float" ]; then
    log_error "config.simulation.M must be a single numeric value per run."
    exit 1
fi

EXPECTED_LINEAR_CONTINUOUS=$(python3 -c "
import json; c=json.load(open('config/config.json'))
if 'continuous' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']))
else:
    print(0)
")
EXPECTED_LINEAR_BINARY=$(python3 -c "
import json; c=json.load(open('config/config.json'))
if 'binary' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']) * len(c['parameters']['p1']))
else:
    print(0)
")
EXPECTED_OD_LINEAR=$((EXPECTED_LINEAR_CONTINUOUS + EXPECTED_LINEAR_BINARY))
EXPECTED_LOGISTIC_CONTINUOUS=$(python3 -c "
import json; c=json.load(open('config/config.json'))
if 'continuous' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']))
else:
    print(0)
")
EXPECTED_LOGISTIC_BINARY=$(python3 -c "
import json; c=json.load(open('config/config.json'))
if 'binary' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']) * len(c['parameters']['p1']))
else:
    print(0)
")
EXPECTED_OD_LOGISTIC=$((EXPECTED_LOGISTIC_CONTINUOUS + EXPECTED_LOGISTIC_BINARY))
CONT_METHODS=$(python3 -c "import json; c=json.load(open('config/config.json')); print(', '.join(m.upper() for m in c['synthesis']['continuous_methods']))")
BIN_METHODS=$(python3 -c "import json; c=json.load(open('config/config.json')); print(', '.join(m.upper() for m in c['synthesis']['binary_methods']))")
eval "$($PYTHON_CMD - <<'PY'
import json
c = json.load(open('config/config.json'))
N = c['simulation']['N']
p = c['simulation']['p']
rho = c['parameters']['rho']
sigma2 = c['parameters']['sigma_2']
p1 = c['parameters']['p1']
M = c['simulation']['M']

num_n = len(N)
num_p = len(p)
num_rho = len(rho)
num_sigma2 = len(sigma2)
num_p1 = len(p1)

od_linear_cont_scenarios = num_n * num_p * num_rho * num_sigma2 if 'continuous' in c['simulation']['var_type'] else 0
od_linear_bin_scenarios = num_n * num_p * num_rho * num_sigma2 * num_p1 if 'binary' in c['simulation']['var_type'] else 0
od_logistic_cont_scenarios = od_linear_cont_scenarios
od_logistic_bin_scenarios = od_linear_bin_scenarios

sum_n = sum(N)
od_linear_cont_target_rows = M * sum_n * num_p * num_rho * num_sigma2 if od_linear_cont_scenarios else 0
od_linear_bin_target_rows = M * sum_n * num_p * num_rho * num_sigma2 * num_p1 if od_linear_bin_scenarios else 0
od_logistic_cont_target_rows = od_linear_cont_target_rows
od_logistic_bin_target_rows = od_linear_bin_target_rows

sd_linear_expected_files = od_linear_cont_scenarios * 2 + od_linear_bin_scenarios * 4
sd_logistic_expected_files = od_logistic_cont_scenarios * 4 + od_logistic_bin_scenarios * 2

sd_linear_target_rows = M * sum_n * ((num_p * num_rho * num_sigma2 * 2) + (num_p * num_rho * num_sigma2 * num_p1 * 4))
sd_logistic_target_rows = M * sum_n * ((num_p * num_rho * num_sigma2 * 4) + (num_p * num_rho * num_sigma2 * num_p1 * 2))

lr_linear_expected_models = sd_linear_expected_files * M * 2
lr_logistic_expected_models = sd_logistic_expected_files * M * 2

for k, v in [
    ('SIM_M', M),
    ('NUM_N', num_n),
    ('NUM_P', num_p),
    ('NUM_RHO', num_rho),
    ('NUM_SIGMA2', num_sigma2),
    ('NUM_P1', num_p1),
    ('OD_LINEAR_CONT_SCENARIOS', od_linear_cont_scenarios),
    ('OD_LINEAR_BIN_SCENARIOS', od_linear_bin_scenarios),
    ('OD_LOGISTIC_CONT_SCENARIOS', od_logistic_cont_scenarios),
    ('OD_LOGISTIC_BIN_SCENARIOS', od_logistic_bin_scenarios),
    ('OD_LINEAR_CONT_TARGET_ROWS', od_linear_cont_target_rows),
    ('OD_LINEAR_BIN_TARGET_ROWS', od_linear_bin_target_rows),
    ('OD_LOGISTIC_CONT_TARGET_ROWS', od_logistic_cont_target_rows),
    ('OD_LOGISTIC_BIN_TARGET_ROWS', od_logistic_bin_target_rows),
    ('SD_LINEAR_EXPECTED_FILES', sd_linear_expected_files),
    ('SD_LOGISTIC_EXPECTED_FILES', sd_logistic_expected_files),
    ('SD_LINEAR_TARGET_ROWS', sd_linear_target_rows),
    ('SD_LOGISTIC_TARGET_ROWS', sd_logistic_target_rows),
    ('LR_LINEAR_EXPECTED_MODELS', lr_linear_expected_models),
    ('LR_LOGISTIC_EXPECTED_MODELS', lr_logistic_expected_models),
]:
    print(f'{k}={v}')
PY
)"
STEP0_END_TIME=$(date +%s)
STEP0_DURATION_S=$((STEP0_END_TIME - STEP0_START_TIME))

echo ""

# ==============================================================================
# STEP 1a: Generate Linear Original Data (R)
# ==============================================================================
log_info "STEP 1a: Generating Linear Original Data (OD_linear) ..."
echo ""

cd 01_generate_OD

EXISTING_OD_LINEAR=$(ls -1 ../data/original/OD_linear_*.parquet 2>/dev/null | wc -l)

if [ "$EXISTING_OD_LINEAR" -ge "$EXPECTED_OD_LINEAR" ]; then
    OD_LINEAR_SKIPPED=1
    log_success "All $EXPECTED_OD_LINEAR linear OD file(s) already exist — skipping Step 1a."
    OD_LINEAR_COUNT=$EXISTING_OD_LINEAR
else
    if [ "$EXISTING_OD_LINEAR" -gt 0 ]; then
        log_info "Resuming linear OD generation: $EXISTING_OD_LINEAR / $EXPECTED_OD_LINEAR already complete."
    fi
    START_TIME=$(date +%s)

    log_info "-> Phase 1a-i: Generating Continuous Data (linear)..."
    PHASE_START_TIME=$(date +%s)
    Rscript 01_generate_original_data_linear.R continuous
    PHASE_END_TIME=$(date +%s)
    OD_LINEAR_CONT_DURATION_S=$((PHASE_END_TIME - PHASE_START_TIME))

    if [ $? -ne 0 ]; then
        log_error "Phase 1a-i (Continuous linear) OD generation failed!"
        exit 1
    fi

    echo ""
    log_info "-> Phase 1a-ii: Generating Binary Data (linear)..."
    PHASE_START_TIME=$(date +%s)
    Rscript 01_generate_original_data_linear.R binary
    PHASE_END_TIME=$(date +%s)
    OD_LINEAR_BIN_DURATION_S=$((PHASE_END_TIME - PHASE_START_TIME))

    if [ $? -ne 0 ]; then
        log_error "Phase 1a-ii (Binary linear) OD generation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    OD_LINEAR_TOTAL_DURATION_S=$DURATION
    OD_LINEAR_COUNT=$(ls -1 ../data/original/OD_linear_*.parquet 2>/dev/null | wc -l)
    log_success "Linear OD generation complete in ${DURATION}s (${OD_LINEAR_COUNT} file(s))"
fi

echo ""

# ==============================================================================
# STEP 1b: Generate Logistic Original Data (R)
# ==============================================================================
log_info "STEP 1b: Generating Logistic Original Data (OD_logistic) ..."
echo ""

EXISTING_OD_LOGISTIC=$(ls -1 ../data/original/OD_logistic_*.parquet 2>/dev/null | wc -l)

if [ "$EXISTING_OD_LOGISTIC" -ge "$EXPECTED_OD_LOGISTIC" ]; then
    OD_LOGISTIC_SKIPPED=1
    log_success "All $EXPECTED_OD_LOGISTIC logistic OD file(s) already exist — skipping Step 1b."
    OD_LOGISTIC_COUNT=$EXISTING_OD_LOGISTIC
else
    if [ "$EXISTING_OD_LOGISTIC" -gt 0 ]; then
        log_info "Resuming logistic OD generation: $EXISTING_OD_LOGISTIC / $EXPECTED_OD_LOGISTIC already complete."
    fi
    START_TIME=$(date +%s)

    log_info "-> Phase 1b-i: Generating Continuous Data (logistic)..."
    PHASE_START_TIME=$(date +%s)
    Rscript 01_generate_original_data_logistic.R continuous
    PHASE_END_TIME=$(date +%s)
    OD_LOGISTIC_CONT_DURATION_S=$((PHASE_END_TIME - PHASE_START_TIME))

    if [ $? -ne 0 ]; then
        log_error "Phase 1b-i (Continuous logistic) OD generation failed!"
        exit 1
    fi

    echo ""
    log_info "-> Phase 1b-ii: Generating Binary Data (logistic)..."
    PHASE_START_TIME=$(date +%s)
    Rscript 01_generate_original_data_logistic.R binary
    PHASE_END_TIME=$(date +%s)
    OD_LOGISTIC_BIN_DURATION_S=$((PHASE_END_TIME - PHASE_START_TIME))

    if [ $? -ne 0 ]; then
        log_error "Phase 1b-ii (Binary logistic) OD generation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    OD_LOGISTIC_TOTAL_DURATION_S=$DURATION
    OD_LOGISTIC_COUNT=$(ls -1 ../data/original/OD_logistic_*.parquet 2>/dev/null | wc -l)
    log_success "Logistic OD generation complete in ${DURATION}s (${OD_LOGISTIC_COUNT} file(s))"
fi

OD_COUNT=$((OD_LINEAR_COUNT + OD_LOGISTIC_COUNT))

cd ..
echo ""

# ==============================================================================
# STEP 2: Generate Synthetic Data for ALL OD files (R)
# ==============================================================================
# Case A (linear+cont): 2 arms; Case C (linear+bin): 4 arms
# Case D (logistic+cont): 4 arms; Case B (logistic+bin): 2 arms
EXPECTED_SD_LINEAR=$((EXPECTED_LINEAR_CONTINUOUS * 2 + EXPECTED_LINEAR_BINARY * 4))
EXPECTED_SD_LOGISTIC=$((EXPECTED_LOGISTIC_CONTINUOUS * 4 + EXPECTED_LOGISTIC_BINARY * 2))
EXPECTED_SD=$((EXPECTED_SD_LINEAR + EXPECTED_SD_LOGISTIC))
EXISTING_SD=$(ls -1 data/synthetic/SD_*.parquet 2>/dev/null | wc -l)

log_info "STEP 2: Generating Synthetic Data (SD) [cont: ${CONT_METHODS}, bin: ${BIN_METHODS}] ..."
log_info "  OD total: ${OD_COUNT} (${OD_LINEAR_COUNT} linear + ${OD_LOGISTIC_COUNT} logistic) -> ${EXPECTED_SD} SD files expected"
echo ""

if [ "$EXISTING_SD" -ge "$EXPECTED_SD" ]; then
    SD_SKIPPED=1
    log_success "All $EXPECTED_SD SD file(s) already exist — skipping Step 2."
    SD_COUNT=$EXISTING_SD
else
    if [ "$EXISTING_SD" -gt 0 ]; then
        log_info "Resuming SD generation: $EXISTING_SD / $EXPECTED_SD already complete."
    fi

    cd 02_generate_SD
    START_TIME=$(date +%s)
    Rscript 02_generate_synthetic_data.R

    if [ $? -ne 0 ]; then
        log_error "SD generation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    SD_DURATION_S=$DURATION
    if [ -f "../results/latest_sd_timing_breakdown.json" ]; then
        eval "$($PYTHON_CMD - <<'PY'
import json
from pathlib import Path

path = Path("../results/latest_sd_timing_breakdown.json")
payload = json.loads(path.read_text())
scaled = payload.get("scaled_wall_s", {})
for shell_name, json_name in [
    ("SD_LINEAR_CONT_DURATION_S", "linear_cont"),
    ("SD_LINEAR_BIN_DURATION_S", "linear_bin"),
    ("SD_LOGISTIC_CONT_DURATION_S", "logistic_cont"),
    ("SD_LOGISTIC_BIN_DURATION_S", "logistic_bin"),
]:
    value = float(scaled.get(json_name, 0.0))
    print(f"{shell_name}={int(round(value))}")
PY
)"
    fi
    SD_COUNT=$(ls -1 ../data/synthetic/SD_*.parquet 2>/dev/null | wc -l)
    log_success "SD generation complete in ${DURATION}s (${SD_COUNT} file(s))"
    cd ..
fi
echo ""

# ==============================================================================
# STEP 3a: Parallel Linear (OLS) Regression (Python)
# ==============================================================================
log_info "STEP 3a: Running Parallel OLS Regressions (linear branch) ..."
echo ""

if [ -f "results/aggregated_model_metrics.parquet" ]; then
    LR_LINEAR_SKIPPED=1
    log_success "Linear regression results already exist — skipping Step 3a."
else
    cd 03_regression_analysis
    START_TIME=$(date +%s)
    $PYTHON_CMD 03_regression_linear.py

    if [ $? -ne 0 ]; then
        log_error "Linear regression evaluation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    LR_LINEAR_DURATION_S=$DURATION
    log_success "Linear regression evaluation complete in ${DURATION}s"
    cd ..
fi
echo ""

# ==============================================================================
# STEP 3b: Parallel Logistic Regression (Python)
# ==============================================================================
log_info "STEP 3b: Running Parallel Logistic Regressions (logistic branch) ..."
echo ""

if [ -f "results/aggregated_logistic_metrics.parquet" ]; then
    LR_LOGISTIC_SKIPPED=1
    log_success "Logistic regression results already exist — skipping Step 3b."
else
    cd 03_regression_analysis
    START_TIME=$(date +%s)
    $PYTHON_CMD 03_regression_logistic.py

    if [ $? -ne 0 ]; then
        log_error "Logistic regression evaluation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    LR_LOGISTIC_DURATION_S=$DURATION
    log_success "Logistic regression evaluation complete in ${DURATION}s"
    cd ..
fi
echo ""

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
echo "=============================================================================="
echo "                    PIPELINE EXECUTION COMPLETE                              "
echo "=============================================================================="
echo ""
log_success "All steps completed successfully!"
echo ""
echo "📁 Output Locations:"
echo "   ├─ Config:            config/config.json"
echo "   ├─ Original Data:     data/original/ (${OD_COUNT} files: ${OD_LINEAR_COUNT} linear + ${OD_LOGISTIC_COUNT} logistic)"
echo "   ├─ Synthetic Data:    data/synthetic/ (${SD_COUNT} files)"
echo "   ├─ Linear Results:    results/aggregated_model_metrics.parquet"
echo "   └─ Logistic Results:  results/aggregated_logistic_metrics.parquet"
echo ""

PIPELINE_END_TIME=$(date +%s)
TOTAL_RUNTIME=$((PIPELINE_END_TIME - PIPELINE_START_TIME))
record_timing_history
RUNTIME_MODEL_STATUS="not_run"
RUNTIME_REPORT_MODE="${RUN_ALL_RUNTIME_REPORT_MODE:-full}"
if [ "$RUNTIME_REPORT_MODE" != "skip" ] && [ -f "04_evaluation/fit_runtime_model.py" ]; then
    log_info "Refreshing runtime approximation model..."
    echo ""
    if $PYTHON_CMD 04_evaluation/fit_runtime_model.py --no-prompt; then
        RUNTIME_MODEL_STATUS="ok"
    else
        RUNTIME_MODEL_STATUS="failed"
        log_warning "Runtime approximation model refresh failed. See results/latest_runtime_model.txt"
    fi
fi
MINUTES=$(((TOTAL_RUNTIME % 3600) / 60))
SECONDS=$((TOTAL_RUNTIME % 60))

if [ $MINUTES -gt 0 ]; then
    log_success "Total pipeline runtime: ${MINUTES}m ${SECONDS}s (${TOTAL_RUNTIME} seconds)"
else
    log_success "Total pipeline runtime: ${SECONDS}s"
fi
echo ""
echo "⏱️ Timing breakdown saved to:"
echo "   ├─ results/pipeline_timing_history.csv"
echo "   └─ results/latest_pipeline_timing.json"
if [ "$RUNTIME_MODEL_STATUS" = "ok" ]; then
    echo "📈 Runtime approximation saved to:"
    echo "   ├─ results/latest_runtime_model.json"
    echo "   └─ results/latest_runtime_model.txt"
fi
echo "=============================================================================="
