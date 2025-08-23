# Test for find_optimal_cost_design function
# Validates cost-optimal experimental design identification

library(testthat)

# Set up test data once for all tests
setup_optimal_cost_test_data <- function() {
  # Create minimal baseline expression data
  baseline_expression_stats <- data.frame(
    response_id = paste0("ENSG", sprintf("%011d", 1:30)),
    relative_expression = exp(rnorm(30, mean = log(1e-5), sd = 0.5)),
    expression_size = runif(30, min = 0.5, max = 2.0)
  )

  # Create library parameters
  library_parameters <- list(
    UMI_per_cell = 12000,
    variation = 0.3
  )

  # Create mock cost-power data frame (simulating output from cost_power_computation)
  set.seed(12345)
  cost_power_df <- expand.grid(
    TPM_threshold = c(5, 10, 15, 20),
    cells_per_target = c(500, 1000, 1500),
    raw_reads_per_cell = c(5000, 8000, 12000)
  ) |>
    dplyr::mutate(
      # Add other required columns
      minimum_fold_change = 0.8,
      num_captured_cells = 100 * cells_per_target / 10, # MOI = 10, num_targets = 100
      library_size = 8000, # Fixed for simplicity
      library_cost = 0.086 * num_captured_cells,
      sequencing_cost = 0.374 * raw_reads_per_cell * num_captured_cells / 1e6,
      total_cost = library_cost + sequencing_cost,
      # Simulate realistic power values (higher for more cells/reads, lower for higher thresholds)
      overall_power = pmax(0.1, pmin(0.95,
        0.5 + 0.3 * log10(cells_per_target / 500) +
        0.2 * log10(raw_reads_per_cell / 5000) -
        0.15 * (TPM_threshold - 5) / 15 + rnorm(36, 0, 0.05)))
    )

  list(
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters,
    cost_power_df = cost_power_df
  )
}

test_that("find_optimal_cost_design basic functionality", {
  test_data <- setup_optimal_cost_test_data()

  result <- find_optimal_cost_design(
    cost_power_df = test_data$cost_power_df,
    minimizing_variable = "TPM_threshold",
    power_target = 0.8,
    power_precision = 0.05
  )

  # Validate output structure
  expect_type(result, "list")
  expect_length(result, 2)
  expect_named(result, c("optimal_cost_power_df", "optimal_cost_grid"))

  # Validate optimal_cost_power_df
  expect_s3_class(result$optimal_cost_power_df, "data.frame")
  required_cols_power <- c("TPM_threshold", "overall_power", "total_cost",
                          "cells_per_target", "raw_reads_per_cell", "minimum_cost")
  expect_true(all(required_cols_power %in% names(result$optimal_cost_power_df)))

  # Validate optimal_cost_grid
  expect_s3_class(result$optimal_cost_grid, "data.frame")
  required_cols_grid <- c("TPM_threshold", "minimum_cost", "raw_reads_per_cell", "total_cost")
  expect_true(all(required_cols_grid %in% names(result$optimal_cost_grid)))

  # Check that optimal_cost_grid has the flattened structure from unnest
  expect_true("raw_reads_per_cell" %in% names(result$optimal_cost_grid))
  expect_true("cells_per_target" %in% names(result$optimal_cost_grid))
})

test_that("find_optimal_cost_design with minimum_fold_change minimizing variable", {
  test_data <- setup_optimal_cost_test_data()

  # Create test data with minimum_fold_change as varying parameter
  cost_power_df_fc <- test_data$cost_power_df |>
    dplyr::select(-minimum_fold_change) |>
    dplyr::slice(rep(1:36, 3)) |>
    dplyr::mutate(
      minimum_fold_change = rep(c(0.6, 0.7, 0.8), nrow(test_data$cost_power_df)),
      # Adjust power based on fold change (lower fold change = easier detection = higher power)
      overall_power = pmax(0.1, pmin(0.95, overall_power + 0.1 * (0.8 - minimum_fold_change)))
    )

  result <- find_optimal_cost_design(
    cost_power_df = cost_power_df_fc,
    minimizing_variable = "minimum_fold_change",
    power_target = 0.7,
    power_precision = 0.05
  )

  # Should group by minimum_fold_change instead of TPM_threshold
  expect_true("minimum_fold_change" %in% names(result$optimal_cost_power_df))
  expect_true("minimum_fold_change" %in% names(result$optimal_cost_grid))

  # Check that we have multiple fold change levels
  expect_gt(length(unique(result$optimal_cost_power_df$minimum_fold_change)), 1)
})

test_that("find_optimal_cost_design cost calculations", {
  test_data <- setup_optimal_cost_test_data()

  result <- find_optimal_cost_design(
    cost_power_df = test_data$cost_power_df,
    minimizing_variable = "TPM_threshold",
    power_target = 0.6,  # Lower target to get more data
    power_precision = 0.1,
    cost_per_captured_cell = 0.1,
    cost_per_million_reads = 0.5
  )

  # Check that minimum_cost is reasonable
  expect_true(all(result$optimal_cost_power_df$minimum_cost > 0))

  # Cost precision should be 1% of minimum cost
  expected_precision <- result$optimal_cost_power_df$minimum_cost / 100
  expect_equal(result$optimal_cost_power_df$cost_precision, expected_precision,
               tolerance = 1e-10)

  # Check cost calculations in the unnested structure
  if (nrow(result$optimal_cost_grid) > 0) {
    grid_sample <- result$optimal_cost_grid

    # Verify cost calculations
    expected_lib_cost <- 0.1 * grid_sample$num_captured_cells
    expected_seq_cost <- 0.5 * grid_sample$raw_reads_per_cell * grid_sample$num_captured_cells / 1e6
    expected_total <- expected_lib_cost + expected_seq_cost

    expect_equal(grid_sample$library_cost, expected_lib_cost, tolerance = 1e-6)
    expect_equal(grid_sample$sequencing_cost, expected_seq_cost, tolerance = 1e-6)
    expect_equal(grid_sample$total_cost, expected_total, tolerance = 1e-6)
  }
})

test_that("find_optimal_cost_design MOI and num_targets effects", {
  test_data <- setup_optimal_cost_test_data()

  # Test with different MOI values
  result_moi5 <- find_optimal_cost_design(
    cost_power_df = test_data$cost_power_df,
    minimizing_variable = "TPM_threshold",
    power_target = 0.6,
    power_precision = 0.1,
    MOI = 5,  # Lower MOI
    num_targets = 100
  )

  result_moi20 <- find_optimal_cost_design(
    cost_power_df = test_data$cost_power_df,
    minimizing_variable = "TPM_threshold",
    power_target = 0.6,
    power_precision = 0.1,
    MOI = 20,  # Higher MOI
    num_targets = 100
  )

  # Lower MOI should result in more captured cells for same cells_per_target
  # But this mainly affects the cost_grid calculations
  expect_s3_class(result_moi5$optimal_cost_power_df, "data.frame")
  expect_s3_class(result_moi20$optimal_cost_power_df, "data.frame")
})

test_that("find_optimal_cost_design different cost_grid_size", {
  test_data <- setup_optimal_cost_test_data()

  result_small_grid <- find_optimal_cost_design(
    cost_power_df = test_data$cost_power_df,
    minimizing_variable = "TPM_threshold",
    power_target = 0.6,
    power_precision = 0.1,
    cost_grid_size = 10  # Small grid
  )

  result_large_grid <- find_optimal_cost_design(
    cost_power_df = test_data$cost_power_df,
    minimizing_variable = "TPM_threshold",
    power_target = 0.6,
    power_precision = 0.1,
    cost_grid_size = 50  # Larger grid
  )

  # Both should work but may have different numbers of grid points
  expect_s3_class(result_small_grid$optimal_cost_grid, "data.frame")
  expect_s3_class(result_large_grid$optimal_cost_grid, "data.frame")
})

test_that("find_optimal_cost_design parameter range validation", {
  test_data <- setup_optimal_cost_test_data()

  result <- find_optimal_cost_design(
    cost_power_df = test_data$cost_power_df,
    minimizing_variable = "TPM_threshold",
    power_target = 0.6,
    power_precision = 0.1
  )

  # Check parameter ranges make sense
  if (nrow(result$optimal_cost_grid) > 0) {
    ranges <- result$optimal_cost_grid |>
      dplyr::select(min_cells, max_cells, min_reads, max_reads)

    expect_true(all(ranges$min_cells <= ranges$max_cells))
    expect_true(all(ranges$min_reads <= ranges$max_reads))
    expect_true(all(ranges$min_cells > 0))
    expect_true(all(ranges$min_reads > 0))
  }
})


test_that("find_optimal_cost_design column renaming", {
  test_data <- setup_optimal_cost_test_data()

  result <- find_optimal_cost_design(
    cost_power_df = test_data$cost_power_df,
    minimizing_variable = "TPM_threshold",
    power_target = 0.6,
    power_precision = 0.1
  )

  # Check that raw_reads_per_cell column exists (no renaming)
  expect_true("raw_reads_per_cell" %in% names(result$optimal_cost_power_df))
  expect_false("reads_per_cell" %in% names(result$optimal_cost_power_df))

  # Same check for cost grid
  expect_true("raw_reads_per_cell" %in% names(result$optimal_cost_grid))
  expect_false("reads_per_cell" %in% names(result$optimal_cost_grid))
})
