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