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
# Launch the interactive app (runs app.R)
perturbplan::launch_app()
```

## Architecture

### Core Components

1. **Power Analysis Pipeline** (`R/power_plan.R`, `R/plan_help.R`)
   - `calculate_power_grid()`: Main function for heatmap power analysis  
   - `calculate_power_curves()`: Detailed power curves for selected conditions
   - `compute_power_grid_efficient()`: Efficient grid-based power analysis using C++ Monte Carlo
   - `compute_power_posthoc()`: Main function for post-hoc power analysis
   - Integrates with C++ implementations for performance

2. **Parameter Estimation** (`R/parameter_estimation.R`, `R/parameter_estimation_help.R`)
   - `library_estimation()`: Estimates parameters from existing data
   - `library_computation()`: Computes QC-aware library statistics

3. **C++ Performance Layer** (`src/`)
   - `BH_cutoff.cpp`: Benjamini-Hochberg multiple testing corrections
   - `compute_distribution_teststat_*.cpp`: Test statistic distributions
   - `power_curves.cpp`: Monte Carlo integration for power curves
   - Key C++ functions: `compute_monte_carlo_teststat_cpp()`, `compute_fc_curve_cpp()`, `compute_expression_curve_cpp()`
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
   - Perturbation-gene pairs to analyze (Random mode or Custom CSV with `grna_target`, `response_id` columns)
   - Analysis parameters (control group, test side, QC thresholds)

2. The package:
   - Validates inputs and CSV format (if Custom mode)
   - Handles gene multiplicity using weighted sampling for duplicate genes in pairs
   - Computes QC-aware library sizes
   - Calculates test statistic distributions using C++
   - Estimates power for each perturbation-gene pair

3. Returns:
   - Individual power for each pair
   - Expected total discoveries

### Key Design Decisions

- **Rcpp Integration**: C++ code handles computationally intensive operations (distribution calculations, multiple testing corrections)
- **Modular Design**: Separate functions for parameter estimation, QC computation, and power analysis allow flexible workflows
- **Shiny Interface**: Provides non-programmatic access via `inst/shiny/app.R`
- **C++ Optimization**: Monte Carlo loops implemented in C++ for significant performance improvements
- **Weighted Sampling**: Preserves gene multiplicity from perturbation-gene pairs using efficient weighted sampling instead of row duplication

## Performance

The package has been optimized for computational efficiency:

- **C++ Monte Carlo**: `.compute_power_plan_efficient()` replaces the older R-based `.compute_underspecified_power_efficient()` with C++ implementations
- **Batch Processing**: Monte Carlo samples processed in batch using `compute_monte_carlo_teststat_cpp()`
- **Efficient Curves**: Power curves computed using optimized C++ functions (`compute_fc_curve_cpp`, `compute_expression_curve_cpp`)
- **Memory-Efficient Sampling**: Uses weighted sampling for gene multiplicity instead of duplicating rows, reducing memory usage while preserving statistical correctness

## Shiny Application Features

### Perturbation-Gene Pairs Analysis

The Shiny app provides an intuitive interface for specifying perturbation-gene pairs:

- **Random Mode**: Randomly samples genes from the baseline expression dataset
- **Custom Mode**: Accepts CSV files with user-specified perturbation-gene pairs
  - Required format: CSV with `grna_target` and `response_id` columns
  - `response_id` must contain Ensembl gene IDs (e.g., ENSG00000141510)
  - Preserves gene multiplicity: genes appearing in multiple pairs get proportional weight in power calculations
  - Example file: `inst/extdata/sample_pairs.csv`

### UI Organization

The application features a streamlined two-tab structure:

1. **Overall Power**: Unified tab with sub-tabs for "Heatmap" and "Slice" views
2. **Drill-down Power**: Detailed power curve analysis for selected experimental conditions

Analysis choices are ordered for logical workflow:
1. **Perturbation-gene pairs to analyze**: Random/Custom dropdown
2. **Minimum TPM threshold**: Gene expression filtering
3. **Test side**: Left (knockdown), Right (overexpression)
4. **Control group**: Complement cells vs Non-targeting cells  
5. **FDR target level**: Multiple testing correction threshold

### Power Visualization Interface

#### Overall Power Tab
- **Heatmap Sub-tab**: Interactive power heatmap with click-to-select functionality
  - Drill-down controls for cells, reads per cell, or both (tiles)
  - Context-sensitive sidebar showing only relevant controls
- **Slice Sub-tab**: Line plots showing power curves for selected heatmap slices
  - Conditional display: shows instruction message when no slices selected
  - Interactive point selection with multiple selection support

#### Drill-down Power Tab
Provides detailed power curve analysis with:

- **Tabbed Interface**: Separate tabs for "Expression" and "Fold Change" plots
- **Display Options**: Control box with three visualization modes:
  - "All together": All experimental designs on a single plot with:
    - Color representing number of cells
    - Linetype and point shape representing reads per cell
    - Legends positioned on the right for optimal space usage
    - Clean legend labels (no redundant "reads/cell" text)
  - "Facet over cells": Horizontal panels separated by cell count, colored by reads per cell
  - "Facet over reads per cell": Horizontal panels separated by reads per cell, colored by cell count
- **Interactive Features**: 
  - ggside marginal histograms showing distribution of expression/fold change values
  - Points added to all line plots for better data visibility
  - Square aspect ratio panels for optimal viewing
  - Consistent 570px box heights across all tabs
- **Performance Optimization**: Default 10×10 grid (instead of 20×20) for faster computation

### File Validation

- Validates CSV format and required columns
- Provides clear error messages for format issues
- Shows loading status: "Loaded X pairs (Y unique genes)"
- Warns about genes filtered out due to low TPM

### Excel Download Organization

The results Excel file is organized with numbered sheets for logical reading:

1. **1_Parameters**: Analysis settings and input parameters
2. **2_Power_Grid**: Main heatmap results (cells × reads per cell power grid)
3. **3_Gene_List**: Input gene list (if custom pairs provided)
4. **4_Selected_Designs**: Information about drill-down selections
5. **5_Fold_Change_Power**: Detailed fold change power curves
6. **6_Expression_Power**: Detailed expression (TPM) power curves

Each sheet uses logical column ordering:
- **Design** column shows "cells × reads" format for easy identification
- **Cells** and **Reads_per_Cell** as separate numeric columns for analysis
- **Data columns** (Expression_TPM, Fold_Change, Power) follow design info
- **Clear naming**: Descriptive column headers without redundancy

## Testing

The package uses testthat (edition 3) with helper functions in `tests/testthat/helper-*.R` for test data generation. Tests compare analytical computations against simulations to ensure accuracy.

## Known Issues

Current R CMD check warnings that need attention:

- **Missing Imports**: Need to declare imports for `Matrix`, `sceptre`, `shiny` packages
- **Namespace Issues**: Missing imports for standard R functions (`setNames`, `read.csv`, `as`)
- **Hidden Files**: `.claude` directory should be added to `.Rbuildignore`

To fix namespace issues, add to NAMESPACE:
```r
importFrom("methods", "as")
importFrom("stats", "setNames")
importFrom("utils", "read.csv")
```

## Parameter Naming Convention

**IMPORTANT**: Use `tpm_threshold` instead of `tmp_threshold` everywhere in the package.

- All function parameters should use `tpm_threshold`
- All variable names should use `tpm_threshold`  
- All documentation should reference `tpm_threshold`
- UI inputs should use `"tpm_threshold"` as input ID

This ensures consistency across the entire codebase and avoids confusion between "TPM" (Transcripts Per Million) and "tmp" (temporary).

## Custom Baseline Expression Upload

The Shiny application supports uploading custom baseline expression data instead of using the default biological system data (K562). This allows users to perform power analysis with their own expression profiles.

### Using Custom Baseline Expression

1. **Navigate to "Pilot data choice" tab** in the sidebar
2. **Select "Custom"** for baseline expression
3. **Upload an RDS file** with the required structure (see below)
4. **Proceed with analysis** - all power calculations will use your custom data

### Required RDS File Structure

The RDS file must contain a list with exactly two elements:

```r
custom_baseline <- list(
  baseline_expression = data.frame(
    response_id = c("ENSG00000141510", "ENSG00000157764", ...),    # Ensembl gene IDs
    relative_expression = c(1.23e-05, 4.56e-06, ...),             # TPM/1e6 scale
    expression_size = c(0.45, 1.23, ...)                          # Dispersion parameters
  ),
  expression_dispersion_curve = function(v) {                     # Dispersion function
    pmax(0.01, 0.1 + 0.5 / sqrt(v))
  }
)

# Save as RDS file
saveRDS(custom_baseline, "my_custom_baseline.rds")
```

### Data Requirements

**baseline_expression data frame:**
- **response_id**: Character vector of gene IDs (preferably Ensembl format: ENSGXXXXXXXXXXX)
- **relative_expression**: Numeric vector of expression levels on TPM/1e6 scale (i.e., raw TPM divided by 1,000,000)
- **expression_size**: Numeric vector of positive dispersion parameters 
- **No missing values** in any column
- **Unique gene IDs** (duplicates will be removed, keeping first occurrence)

**expression_dispersion_curve function:**
- Must accept a numeric vector and return dispersion values of the same length
- Should model the mean-variance relationship in your expression data
- Example: `function(v) pmax(0.01, 0.1 + 0.5 / sqrt(v))`

### Creating Custom Baseline Files

#### Method 1: From Existing Data
```r
# Load the package and default data
library(perturbplan)
default_data <- extract_baseline_expression("K562")

# Subset or modify as needed
my_baseline <- default_data$baseline_expression[1:2000, ]  # Use first 2000 genes

# Create custom baseline structure
custom_baseline <- list(
  baseline_expression = my_baseline,
  expression_dispersion_curve = default_data$expression_dispersion_curve
)

# Save as RDS
saveRDS(custom_baseline, "my_custom_baseline.rds")
```

#### Method 2: From Scratch
Use the example script at `inst/extdata/create_custom_baseline_example.R` for guidance on creating baseline data from your own expression measurements.

### File Validation

The application automatically validates uploaded RDS files and provides detailed error messages for:
- Incorrect file structure or missing elements
- Invalid data types or value ranges
- Missing values or duplicate gene IDs
- File size limits (50MB maximum)
- R version compatibility issues

### Integration with Analysis Workflow

Custom baseline expression data integrates seamlessly with all analysis features:
- **Compatible with both Random and Custom gene list modes**
- **Works with all analysis parameters** (test side, control group, FDR levels)
- **Included in Excel downloads** with clear documentation of data source
- **Supports all visualization features** (heatmaps, power curves, drill-down analysis)

### Performance Considerations

- **File size**: Keep RDS files under 50MB for optimal performance
- **Gene count**: 1,000-10,000 genes typically provide good balance of comprehensiveness and speed
- **Memory usage**: Large datasets may require more RAM for analysis

### Example Files

Pre-built example files are available:
- `inst/extdata/example_custom_baseline.rds`: Subset of default K562 data (1,000 genes)
- `inst/extdata/create_custom_baseline_example.R`: Script for creating custom files

## Custom Library Parameters Upload

The Shiny application supports uploading custom library parameters instead of using the default biological system data (K562). This allows users to perform power analysis with their own UMI saturation and PCR bias parameters.

### Using Custom Library Parameters

1. **Navigate to "Pilot data choice" tab** in the sidebar
2. **Expand "Library size parameters" section**
3. **Select "Custom"** for library parameters
4. **Upload an RDS file** with the required structure (see below)
5. **Proceed with analysis** - all power calculations will use your custom parameters

### Required RDS File Structure

The RDS file must contain a list with exactly two elements:

```r
custom_library <- list(
  UMI_per_cell = 15000,    # Maximum UMI per cell parameter (positive numeric)
  variation = 0.25         # Variation parameter for PCR bias (positive numeric)
)

# Save as RDS file
saveRDS(custom_library, "my_custom_library.rds")
```

### Data Requirements

**Library parameters:**
- **UMI_per_cell**: Maximum UMI per cell parameter from saturation curve fitting (typically 1000-50000)
- **variation**: Variation parameter characterizing PCR amplification bias (typically 0.1-1.0)
- **Both parameters** must be positive single numeric values
- **No missing values** allowed

### Data Validation

The application automatically validates uploaded RDS files and provides detailed error messages for:
- Incorrect file structure or missing elements
- Invalid data types or value ranges
- Missing values or non-finite numbers
- File size limits (50MB maximum)
- R version compatibility issues

### Creating Custom Library Files

#### Method 1: From Default Parameters
```r
# Load the package and default data
library(perturbplan)
default_library <- extract_library_info("K562")

# Modify parameters as needed
custom_library <- list(
  UMI_per_cell = default_library$UMI_per_cell * 1.5,  # 50% higher capacity
  variation = default_library$variation * 0.8         # 20% lower bias
)

# Save as RDS
saveRDS(custom_library, "modified_k562_library.rds")
```

#### Method 2: From Your Own Measurements
```r
# Based on your saturation curve fitting results
custom_library <- list(
  UMI_per_cell = 18000,    # Your measured UMI capacity
  variation = 0.22         # Your measured PCR bias
)

# Validate before saving
validation_result <- validate_custom_library_rds(custom_library)
if (validation_result$valid) {
  saveRDS(custom_library, "my_measured_library.rds")
} else {
  cat("Validation errors:", paste(validation_result$errors, collapse = ", "))
}
```

### Integration with Analysis Workflow

Custom library parameters integrate seamlessly with all analysis features:
- **Compatible with both Random and Custom gene list modes**
- **Works with all analysis parameters** (test side, control group, FDR levels)
- **Included in Excel downloads** with clear documentation of parameter source
- **Supports all visualization features** (heatmaps, power curves, drill-down analysis)

### Performance Considerations

- **File size**: Keep RDS files small for optimal performance
- **Parameter ranges**: Extreme values may affect analysis speed or accuracy
- **Memory usage**: Large parameter differences from defaults may require more computation time

### Example Files

Pre-built example files are available:
- `inst/extdata/example_custom_library.rds`: Simple example with moderate parameters
- `inst/extdata/create_custom_library_example.R`: Script for creating custom files

### Summary Display

When custom library parameters are loaded, the application displays:
"Loaded custom library parameters  
UMI per cell: XX,XXX  
Variation: X.XXXe-XX"

## Development Notes

- **Function Migration**: `.compute_underspecified_power_efficient()` has been replaced with `.compute_power_plan_efficient()` for better performance
- **C++ Priority**: When possible, use C++ implementations over R loops for computationally intensive operations
- **Grid Analysis**: Use `compute_power_grid_efficient()` for systematic power analysis across experimental conditions
- **ggside Faceting**: When using ggplot2 faceting with ggside histograms, convert numeric faceting variables to factors explicitly to avoid "Can't combine factor and double" errors. Use `factor()` with proper levels and labels before `facet_grid()`.

## Git Workflow Requirements

**IMPORTANT**: When commit and push is requested, the **entire repository** should be committed and pushed, not just specific changes.

- **Complete Sync**: After commit and push, there should be **no difference** between the local directory and remote repository
- **Clean Working Tree**: `git status` should show a clean working tree after pushing
- **Full Commit**: Use `git add .` to stage all changes before committing, unless specifically instructed to commit only particular files
- **Repository Consistency**: The remote repository should always reflect the complete current state of the local development environment

This ensures repository consistency and prevents issues with uncommitted changes being left behind during development sessions.

## Shiny UI Tab Modification Guidelines

**IMPORTANT**: When making modifications to tabs in the Shiny UI, always ensure changes align with existing patterns:

- **Tab Header Colors**: Follow the established color scheme defined in CSS selectors (e.g., `#exp-header`, `#perturbation-header`, `#analysis-header`, `#effects-header`)
- **Collapsibility**: Maintain the collapsible functionality with proper JavaScript integration
- **CSS Consistency**: Update all relevant CSS selectors and JavaScript arrays when adding/removing/modifying tabs
- **Pattern Matching**: New tabs should follow the exact same structure as existing tabs:
  - Header styling with hover effects
  - Chevron icons with proper rotation
  - Content containers with consistent padding and background
  - Display states (`display: none` for collapsed, `display: block` for expanded)

When adding or modifying tabs, check:
1. CSS selectors in `ui_styles.R` include the new tab IDs
2. JavaScript arrays (`allSections`, `allChevrons`) are updated
3. Initial state setup includes the new tab
4. Color scheme and hover effects match existing tabs