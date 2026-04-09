#!/bin/bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

echo "=============================================================================="
echo "                         PIPELINE RUNTIME ESTIMATOR                           "
echo "=============================================================================="
echo ""

python3 04_evaluation/fit_runtime_model.py --estimator
