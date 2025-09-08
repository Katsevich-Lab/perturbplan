#!/bin/bash

#$ -N run_obtain_pilot_parameter
#$ -cwd
#$ -pe openmp 4
#$ -l m_mem_free=10G        # Request 16 GB memory *per core*
#$ -j y                     # Join stdout and stderr

set -euo pipefail

echo "Running install helper..."

PROJECT_ROOT="/home/mnt/weka/yihuihe/code/perturbplan/"
# SCRIPT_PATH="$(readlink -f "$0")"
# SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
# PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "Now in package root: $PROJECT_ROOT"

cd "$PROJECT_ROOT"

Rscript "data-raw/process_reference_expression_data.R"
