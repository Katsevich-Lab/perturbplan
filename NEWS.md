# perturbplan 0.2.0

## Major Changes

### Parameter Naming Standardization
- **BREAKING CHANGE**: Renamed `reads_per_cell` and `raw_reads_per_cell` parameters to `sequenced_reads_per_cell` in some functions
  - This change clarifies that the parameter refers to raw sequencing reads (before mapping), not mapped reads
  - Affects: `compute_power_plan()`, `compute_power_posthoc()`, `cost_power_computation()`, `find_optimal_cost_design()`, and related functions
- Improved internal handling of mapping efficiency with automatic conversion and rounding

### New Vignettes and Documentation
- **New vignette**: "Prepare Data For Web App" (`preprocess-reference.Rmd`)
  - Comprehensive guide for preprocessing Cell Ranger outputs
  - Documentation for `reference_data_preprocessing_10x()` and `reference_data_processing()`
  - Examples for both Perturb-seq and TAP-seq workflows
  - Mathematical details of negative binomial expression model and read-UMI saturation curve
- Enhanced `prospective-power` vignette with updated parameter names
- Improved documentation for all pilot data preprocessing functions

### Enhanced Pilot Data Preprocessing
- Added `gene_list` parameter to `reference_data_processing()` for TAP-seq targeted gene panel support
- Created example TAP-seq pilot data with targeted gene list
- Improved handling of different prior values for variation parameter (`D2_rough`)
- Enhanced example raw data to reduce numerical singularity issues

### Bug Fixes
- Fixed Matrix package warnings for sparse matrix operations
- Resolved devtools::check() issues
- Fixed parameter passing in cost optimization functions
- Improved numerical stability in dispersion estimation

### Website and Documentation Improvements
- Redesigned pkgdown site structure and homepage
- Added hex logo
- Rebuilt all documentation with updated parameter names
- Improved function examples and cross-references
- Enhanced vignette figures and visualizations

## Testing
- All tests updated for new parameter naming convention
- Enhanced test coverage for parameter validation
- All R CMD check tests passing (0 errors, 0 warnings)

---

# perturbplan 0.1.0

## New Features

### Cost Optimization Functions
- `cost_power_computation()`: Comprehensive power analysis with cost minimization across multiple experimental parameters
- `find_optimal_cost_design()`: Binary search optimization to find experimental designs meeting power targets
- `obtain_fixed_variable_constraining_cost()`: Helper function for cost-constrained experimental design

### Enhanced Parameter Support
- Extended `minimizing_variable` parameter validation to support:
  - `"TPM_threshold"` - Expression threshold optimization
  - `"minimum_fold_change"` - Effect size optimization
  - `"cells_per_target"` - Cell count optimization (cost_power_computation only)
  - `"reads_per_cell"` - Sequencing depth optimization (cost_power_computation only)
  - `"cost"` - Total cost optimization
- Differentiated validation between `cost_power_computation` and `find_optimal_cost_design`

### Documentation Improvements
- Added comprehensive `@examples` sections to all major functions
- Enhanced function documentation with detailed parameter descriptions
- Added `@keywords internal` annotations for helper functions

## Bug Fixes

### Critical Fixes
- **Fixed minimum_fold_change < 0.6 bug**: Added safety check and fallback mechanism when cell range identification fails due to small effect sizes
- **Parameter naming standardization**: Changed all `tmp_threshold` variables to `TPM_threshold` throughout codebase to prevent typos

### Validation Improvements
- Enhanced input validation for cost optimization functions
- Added NSE-safe dplyr operations using `.data$` notation
- Improved error messages and parameter checking

## Performance Improvements
- Optimized switch statement logic for new minimizing variables
- Efficient handling of cost calculations across experimental parameter grids
- Enhanced C++ integration for computationally intensive operations

## Package Structure
- Updated `.Rbuildignore` and `.gitignore` for pkgdown support
- Comprehensive test coverage for all new cost optimization features
- Maintained backward compatibility with existing power analysis functions

---

**Breaking Changes**: None - all new features are additive and maintain full backward compatibility.

**Dependencies**: No new dependencies added for core functionality.
