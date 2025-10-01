# Test for input validation functions
# Validates input checking functions properly catch invalid inputs

library(testthat)

# Set up test data for validation tests
setup_validation_test_data <- function() {
  # Create valid baseline expression data
  baseline_expression_stats <- data.frame(
    response_id = paste0("ENSG", sprintf("%011d", 1:20)),
    relative_expression = exp(rnorm(20, mean = log(1e-5), sd = 0.5)),
    expression_size = runif(20, min = 0.5, max = 2.0)
  )

  # Create valid library parameters
  library_parameters <- list(
    UMI_per_cell = 12000,
    variation = 0.3
  )

  list(
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters
  )
}

test_that("input_check_compute_power_planvalidates TPM_threshold parameter", {
  test_data <- setup_validation_test_data()

  # Test invalid TPM_threshold type
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = "invalid",
      minimum_fold_change = 0.8,
      cells_per_target = 1000,
      reads_per_cell = 8000,
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10
    ),
    regexp = "must be numeric or the string 'varying'"
  )

  # Test negative tmp_threshold
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = -5,
      minimum_fold_change = 0.8,
      cells_per_target = 1000,
      reads_per_cell = 8000,
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10
    ),
    regexp = "`TPM_threshold` values must be non-negative!"
  )

  # Test valid inputs - should not throw error
  expect_silent(
    input_check_compute_power_plan(
      TPM_threshold = 10,
      minimum_fold_change = 0.8,
      cells_per_target = 1000,
      reads_per_cell = 8000,
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10
    )
  )
})

test_that("input_check_compute_power_planvalidates minimum_fold_change parameter", {
  test_data <- setup_validation_test_data()

  # Test invalid minimum_fold_change type
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 10,
      minimum_fold_change = "invalid",
      cells_per_target = 1000,
      reads_per_cell = 8000,
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10
    ),
    regexp = "must be numeric or the string 'varying'"
  )

  # Test negative minimum_fold_change
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 10,
      minimum_fold_change = -0.5,
      cells_per_target = 1000,
      reads_per_cell = 8000,
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10
    ),
    regexp = "values must be positive"
  )
})

test_that("input_check_compute_power_planvalidates baseline expression data", {
  test_data <- setup_validation_test_data()

  # Test missing baseline_expression_stats
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 10,
      minimum_fold_change = 0.8,
      cells_per_target = 1000,
      reads_per_cell = 8000,
      baseline_expression_stats = NULL,
      library_parameters = test_data$library_parameters,
      num_targets = 10
    ),
    regexp = "must be a data frame"
  )

  # Test invalid baseline_expression_stats structure
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 10,
      minimum_fold_change = 0.8,
      cells_per_target = 1000,
      reads_per_cell = 8000,
      baseline_expression_stats = data.frame(wrong_column = 1:5),
      library_parameters = test_data$library_parameters,
      num_targets = 10
    ),
    regexp = "must contain columns: response_id, relative_expression, expression_size"
  )
})

test_that("input_check_compute_power_planvalidates library parameters", {
  test_data <- setup_validation_test_data()

  # Test missing library_parameters
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 10,
      minimum_fold_change = 0.8,
      cells_per_target = 1000,
      reads_per_cell = 8000,
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = NULL,
      num_targets = 10
    ),
    regexp = "must be a list"
  )

  # Test invalid library_parameters structure
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 10,
      minimum_fold_change = 0.8,
      cells_per_target = 1000,
      reads_per_cell = 8000,
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = list(wrong_param = 100),
      num_targets = 10
    ),
    regexp = "must contain elements: UMI_per_cell, variation"
  )
})

test_that("input_check_cost_power_computation validates minimizing_variable parameter", {
  test_data <- setup_validation_test_data()

  # Test invalid minimizing_variable
  expect_error(
    input_check_cost_power_computation(
      minimizing_variable = "invalid_param",
      fixed_variable = list(minimum_fold_change = 0.8),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10,
      power_target = 0.6,
      cost_constraint = 10000
    ),
    regexp = "must be one of: TPM_threshold, minimum_fold_change"
  )

  # Test valid minimizing_variable - should not throw error
  expect_silent(
    input_check_cost_power_computation(
      minimizing_variable = "TPM_threshold",
      fixed_variable = list(minimum_fold_change = 0.8),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10,
      power_target = 0.6,
      cost_constraint = 10000
    )
  )
})

test_that("input_check_cost_power_computation validates power_target parameter", {
  test_data <- setup_validation_test_data()

  # Test invalid power_target range
  expect_error(
    input_check_cost_power_computation(
      minimizing_variable = "TPM_threshold",
      fixed_variable = list(minimum_fold_change = 0.8),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10,
      power_target = 1.5,  # > 1
      cost_constraint = 10000
    ),
    regexp = "must be a numeric value in \\(0,1\\)"
  )

  # Test negative power_target
  expect_error(
    input_check_cost_power_computation(
      minimizing_variable = "TPM_threshold",
      fixed_variable = list(minimum_fold_change = 0.8),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10,
      power_target = -0.1,
      cost_constraint = 10000
    ),
    regexp = "must be a numeric value in \\(0,1\\)"
  )
})

test_that("input_check_cost_power_computation validates cost_constraint parameter", {
  test_data <- setup_validation_test_data()

  # Test negative cost_constraint
  expect_error(
    input_check_cost_power_computation(
      minimizing_variable = "TPM_threshold",
      fixed_variable = list(minimum_fold_change = 0.8),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10,
      power_target = 0.6,
      cost_constraint = -1000
    ),
    regexp = "must be NULL or a positive numeric value"
  )

  # Test NULL cost_constraint - should not throw error
  expect_silent(
    input_check_cost_power_computation(
      minimizing_variable = "TPM_threshold",
      fixed_variable = list(minimum_fold_change = 0.8),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10,
      power_target = 0.6,
      cost_constraint = NULL
    )
  )
})

test_that("input_check_cost_power_computation validates fixed_variable parameter", {
  test_data <- setup_validation_test_data()

  # Test missing fixed_variable
  expect_error(
    input_check_cost_power_computation(
      minimizing_variable = "TPM_threshold",
      fixed_variable = NULL,
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10,
      power_target = 0.6,
      cost_constraint = 10000
    ),
    regexp = "must be a list"
  )

  # Test empty fixed_variable
  expect_error(
    input_check_cost_power_computation(
      minimizing_variable = "TPM_threshold",
      fixed_variable = list(),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10,
      power_target = 0.6,
      cost_constraint = 10000
    ),
    regexp = "must contain 'minimum_fold_change'"
  )

  # Test invalid fixed_variable parameter names
  expect_error(
    input_check_cost_power_computation(
      minimizing_variable = "TPM_threshold",
      fixed_variable = list(invalid_param = 0.8),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10,
      power_target = 0.6,
      cost_constraint = 10000
    ),
    regexp = "must contain 'minimum_fold_change'"
  )
})

test_that("input validation functions work with varying parameters", {
  test_data <- setup_validation_test_data()

  # Test with "varying" parameters - should not throw error
  expect_silent(
    input_check_compute_power_plan(
      TPM_threshold = "varying",
      minimum_fold_change = "varying",
      cells_per_target = "varying",
      reads_per_cell = "varying",
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10
    )
  )

  # Test with mixed numeric and varying parameters - should not throw error
  expect_silent(
    input_check_compute_power_plan(
      TPM_threshold = c(5, 10, 15),
      minimum_fold_change = "varying",
      cells_per_target = 1000,
      reads_per_cell = c(6000, 8000, 10000),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 10
    )
  )
})
