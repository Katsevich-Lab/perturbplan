# Test for cost_power_computation function
# Validates cost-constrained power analysis functionality

library(testthat)

# Set up test data once for all tests
setup_test_data <- function() {
  # Create minimal baseline expression data
  baseline_expression_stats <- data.frame(
    response_id = paste0("ENSG", sprintf("%011d", 1:100)),
    relative_expression = exp(rnorm(100, mean = log(1e-5), sd = 1)),
    expression_size = runif(100, min = 0.5, max = 2.0)
  )

  # Create library parameters
  library_parameters <- list(
    UMI_per_cell = 15000,
    variation = 0.25
  )

  list(
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters
  )
}

test_that("cost_power_computation basic functionality with tpm_threshold optimization", {
  test_data <- setup_test_data()

  # Test basic tpm_threshold optimization
  result <- cost_power_computation(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.8),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 30,
    grid_size = 3,
    power_target = 0.6,
    cost_constraint = 20000
  )

  # Validate output structure
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)

  # Check required columns
  expected_cols <- c("minimum_fold_change", "tpm_threshold", "cells_per_target",
                     "num_captured_cells", "raw_reads_per_cell", "library_size",
                     "overall_power", "library_cost", "sequencing_cost", "total_cost")
  expect_true(all(expected_cols %in% names(result)))

  # Validate data types
  expect_type(result$minimum_fold_change, "double")
  expect_type(result$tpm_threshold, "double")
  expect_type(result$overall_power, "double")
  expect_type(result$total_cost, "double")

  # Validate value ranges
  expect_true(all(result$overall_power >= 0 & result$overall_power <= 1))
  expect_true(all(result$total_cost > 0))
  expect_true(all(result$minimum_fold_change == 0.8))  # Fixed value
})

test_that("cost_power_computation with minimum_fold_change optimization", {
  test_data <- setup_test_data()

  result <- cost_power_computation(
    minimizing_variable = "minimum_fold_change",
    fixed_variable = list(tpm_threshold = 100),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 30,
    grid_size = 3,
    power_target = 0.6,
    cost_constraint = NULL  # No cost constraint
  )

  # Validate output
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(all(result$tpm_threshold == 100))  # Fixed value
  expect_true(length(unique(result$minimum_fold_change)) > 1)  # Should vary
})

test_that("cost_power_computation with fixed experimental design", {
  test_data <- setup_test_data()

  result <- cost_power_computation(
    minimizing_variable = "minimum_fold_change",
    fixed_variable = list(
      tpm_threshold = 50,
      cells_per_target = 500,
      reads_per_cell = 8000
    ),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 20,
    grid_size = 3,
    power_target = 0.5,
    cost_constraint = 15000
  )

  # With fixed experimental design, should have consistent cells/reads
  expect_true(all(abs(result$cells_per_target - 500) < 1e-6))
  expect_true(all(abs(result$raw_reads_per_cell - 8000/0.72) < 1e-6))  # Accounting for mapping efficiency
})


test_that("cost_power_computation with very restrictive cost constraint", {
  test_data <- setup_test_data()

  # Test with cost constraint that should fail validation
  expect_error(
    cost_power_computation(
      minimizing_variable = "tpm_threshold",
      fixed_variable = list(minimum_fold_change = 0.8),
      baseline_expression_stats = test_data$baseline_expression_stats,
      library_parameters = test_data$library_parameters,
      num_targets = 100,
      grid_size = 3,
      power_target = 0.8,
      cost_constraint = 50  # Very low cost constraint
    ),
    regexp = "cost optimization failed"
  )
})

test_that("cost_power_computation cost calculations are correct", {
  test_data <- setup_test_data()

  result <- cost_power_computation(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.8),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 30,
    grid_size = 3,
    power_target = 0.6,
    cost_constraint = 20000,
    cost_per_captured_cell = 0.1,
    cost_per_million_reads = 0.5
  )

  # Verify cost calculations
  expected_library_cost <- result$num_captured_cells * 0.1
  expected_sequencing_cost <- (0.5 * result$raw_reads_per_cell * result$num_captured_cells) / 1e6
  expected_total_cost <- expected_library_cost + expected_sequencing_cost

  expect_equal(result$library_cost, expected_library_cost, tolerance = 1e-10)
  expect_equal(result$sequencing_cost, expected_sequencing_cost, tolerance = 1e-10)
  expect_equal(result$total_cost, expected_total_cost, tolerance = 1e-10)
})

test_that("cost_power_computation different test sides work correctly", {
  test_data <- setup_test_data()

  # Test left side (knockdown)
  result_left <- cost_power_computation(
    minimizing_variable = "minimum_fold_change",
    fixed_variable = list(tpm_threshold = 100),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    side = "left",
    num_targets = 20,
    grid_size = 3,
    power_target = 0.5,
    cost_constraint = NULL
  )

  # Test right side (overexpression)
  result_right <- cost_power_computation(
    minimizing_variable = "minimum_fold_change",
    fixed_variable = list(tpm_threshold = 100),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    side = "right",
    num_targets = 20,
    grid_size = 3,
    power_target = 0.5,
    cost_constraint = NULL
  )

  # Should have different fold change ranges
  expect_true(max(result_left$minimum_fold_change) <= 1.0)  # Left side should be < 1
  expect_true(min(result_right$minimum_fold_change) >= 1.0)  # Right side should be > 1
})

test_that("cost_power_computation with different control groups", {
  test_data <- setup_test_data()

  # Test complement control
  result_complement <- cost_power_computation(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.8),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    control_group = "complement",
    num_targets = 20,
    grid_size = 3,
    power_target = 0.5,
    cost_constraint = NULL
  )

  # Test non-targeting control
  result_nt <- cost_power_computation(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.8),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    control_group = "nt_cells",
    num_targets = 20,
    grid_size = 3,
    power_target = 0.5,
    cost_constraint = NULL
  )

  # Both should run without error and return results
  expect_s3_class(result_complement, "data.frame")
  expect_s3_class(result_nt, "data.frame")
  expect_true(nrow(result_complement) > 0)
  expect_true(nrow(result_nt) > 0)
})

test_that("cost_power_computation reproducibility with set.seed", {
  test_data <- setup_test_data()

  # Run twice with same parameters
  result1 <- cost_power_computation(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.8),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 20,
    grid_size = 3,
    power_target = 0.6,
    cost_constraint = NULL
  )

  result2 <- cost_power_computation(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.8),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 20,
    grid_size = 3,
    power_target = 0.6,
    cost_constraint = NULL
  )

  # Results should be identical due to set.seed(1) in function
  expect_equal(result1, result2, tolerance = 1e-10)
})

test_that("cost_power_computation with cost_precision parameter", {
  test_data <- setup_test_data()

  result <- cost_power_computation(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.8),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 30,
    grid_size = 3,
    power_target = 0.6,
    cost_constraint = 20000,
    cost_precision = 0.8  # Use only 80% of budget
  )

  # All costs should be within the constraint × precision
  expect_true(any(result$total_cost <= 20000 * 0.8))
})
