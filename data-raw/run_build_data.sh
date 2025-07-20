#!/bin/bash

#$ -N run_obtain_pilot_parameter
#$ -cwd
#$ -pe openmp 4
#$ -l m_mem_free=16G        # Request 16 GB memory *per core*
#$ -j y                     # Join stdout and stderr

echo "Running install helper..."
cd "$(dirname "$0")/.."
Rscript data_raw/build_reference_expression_data.R