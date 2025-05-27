# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

PerturbPlan is an R package for experimental design and power analysis for perturb-seq experiments (CRISPR-based single-cell perturbation experiments). It combines R and C++ code for efficient statistical computations.

## Common Development Commands

### Build and Check
```bash
# Build the package
R CMD build .

# Check the package (replace version as needed)
R CMD check perturbplan_0.0.1.tar.gz

# Install the package locally
R CMD INSTALL .
```

### Development Workflow
```r
# Load package during development
devtools::load_all()

# Run all tests
devtools::test()

# Run a specific test file
devtools::test(filter = "test-library_computation")

# Generate documentation from roxygen2 comments
devtools::document()

# Check package without building
devtools::check()
```

### Shiny App
```r
# Launch the interactive app (runs app-prototype.R)
perturbplan::launch_app()
```

## Architecture

### Core Components

1. **Power Analysis Pipeline** (`R/power.R`, `R/power_help.R`)
   - `compute_power_posthoc()`: Main function for post-hoc power analysis
   - Integrates with C++ implementations for performance

2. **Parameter Estimation** (`R/parameter_estimation.R`, `R/parameter_estimation_help.R`)
   - `library_estimation()`: Estimates parameters from existing data
   - `library_computation()`: Computes QC-aware library statistics

3. **C++ Performance Layer** (`src/`)
   - `BH_cutoff.cpp`: Benjamini-Hochberg multiple testing corrections
   - `compute_distribution_teststat_*.cpp`: Test statistic distributions
   - Uses Rcpp for seamless R/C++ integration

4. **Quality Control** (`R/QC_computation.R`)
   - Implements pairwise QC checks for perturb-seq data
   - Filters based on minimum non-zero cell counts

5. **Input Validation** (`R/check.R`)
   - Comprehensive validation functions for all major operations
   - Ensures data consistency across the pipeline

### Data Flow

1. User provides:
   - Cell counts per gRNA
   - Baseline expression statistics
   - Discovery pairs (perturbation-gene combinations to test)
   - Analysis parameters (control group, test side, QC thresholds)

2. The package:
   - Validates inputs
   - Computes QC-aware library sizes
   - Calculates test statistic distributions using C++
   - Estimates power for each perturbation-gene pair

3. Returns:
   - Individual power for each pair
   - Expected total discoveries

### Key Design Decisions

- **Rcpp Integration**: C++ code handles computationally intensive operations (distribution calculations, multiple testing corrections)
- **Modular Design**: Separate functions for parameter estimation, QC computation, and power analysis allow flexible workflows
- **Shiny Interface**: Provides non-programmatic access via `inst/shiny/app-prototype.R`

## Testing

The package uses testthat (edition 3) with helper functions in `tests/testthat/helper-*.R` for test data generation. Tests compare analytical computations against simulations to ensure accuracy.