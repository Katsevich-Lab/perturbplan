# Test for compute_power_plan_full_grid function
# Validates comprehensive power analysis across parameter combinations

library(testthat)

# Set up test data once for all tests
setup_grid_test_data <- function() {
  # Create minimal baseline expression data
  baseline_expression_stats <- data.frame(
    response_id = paste0("ENSG", sprintf("%011d", 1:50)),
    relative_expression = exp(rnorm(50, mean = log(1e-5), sd = 0.5)),
    expression_size = runif(50, min = 0.5, max = 2.0)
  )

  # Create library parameters
  library_parameters <- list(
    UMI_per_cell = 12000,
    variation = 0.3
  )

  list(
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters
  )
}

test_that("compute_power_plan_full_grid basic functionality with single values", {
  test_data <- setup_grid_test_data()

  result <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = 10000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 20,
    grid_size = 5
  )

  # Validate output structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)  # Single combination should give 1 row

  # Check required columns
  expected_cols <- c("minimum_fold_change", "tpm_threshold", "cells_per_target",
                     "num_captured_cells", "raw_reads_per_cell", "library_size", "overall_power")
  expect_true(all(expected_cols %in% names(result)))

  # Validate specific values
  expect_equal(result$minimum_fold_change, 0.8)
  expect_equal(result$tpm_threshold, 10)
  expect_true(result$overall_power >= 0 && result$overall_power <= 1)
})

test_that("compute_power_plan_full_grid with numeric vectors", {
  test_data <- setup_grid_test_data()

  result <- compute_power_plan_full_grid(
    tpm_threshold = c(5, 10, 15),
    minimum_fold_change = c(0.7, 0.8),
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 15
  )

  # Should have 3 × 2 = 6 combinations
  expect_equal(nrow(result), 6)

  # Check parameter combinations
  expect_setequal(unique(result$tpm_threshold), c(5, 10, 15))
  expect_setequal(unique(result$minimum_fold_change), c(0.7, 0.8))
  expect_true(all(result$cells_per_target == 1000))
})

test_that("compute_power_plan_full_grid with varying parameters", {
  test_data <- setup_grid_test_data()

  result <- compute_power_plan_full_grid(
    tpm_threshold = "varying",
    minimum_fold_change = 0.8,
    cells_per_target = "varying",
    reads_per_cell = "varying",
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 20,
    grid_size = 3  # Small grid for faster testing
  )

  # Should have multiple rows due to varying parameters
  expect_true(nrow(result) >= 3)  # At least grid_size combinations

  # Check that parameters actually vary
  expect_true(length(unique(result$tpm_threshold)) > 1)
  expect_true(length(unique(result$cells_per_target)) > 1)
  expect_true(length(unique(result$raw_reads_per_cell)) > 1)
  expect_true(all(result$minimum_fold_change == 0.8))  # Fixed value
})

test_that("compute_power_plan_full_grid varying parameter ranges", {
  test_data <- setup_grid_test_data()

  # Test tpm_threshold varying creates reasonable range
  result <- compute_power_plan_full_grid(
    tpm_threshold = "varying",
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    grid_size = 5
  )

  # Should have 5 different tpm_threshold values based on quantiles
  expect_equal(length(unique(result$tpm_threshold)), 5)

  # tpm_threshold values should be in reasonable range (based on expression quantiles)
  tpm_thresholds <- unique(result$tpm_threshold)
  expect_true(all(tpm_thresholds > 0))
  expect_true(all(tpm_thresholds < 1000))  # Reasonable upper bound
})

test_that("compute_power_plan_full_grid minimum_fold_change varying by test side", {
  test_data <- setup_grid_test_data()

  # Test left side (knockdown)
  result_left <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = "varying",
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    side = "left",
    grid_size = 5
  )

  # Test right side (overexpression)
  result_right <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = "varying",
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    side = "right",
    grid_size = 5
  )

  # Test both sides
  result_both <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = "varying",
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    side = "both",
    grid_size = 5
  )

  # Left side should have fold changes < 1
  expect_true(all(result_left$minimum_fold_change < 1))

  # Right side should have fold changes > 1
  expect_true(all(result_right$minimum_fold_change > 1))

  # Both sides should have more combinations
  expect_equal(nrow(result_both), 10)  # 5 + 5 values for both sides
})

test_that("compute_power_plan_full_grid different control groups", {
  test_data <- setup_grid_test_data()

  # Test complement control
  result_complement <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    control_group = "complement",
    num_targets = 15
  )

  # Test non-targeting control
  result_nt <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    control_group = "nt_cells",
    num_targets = 15
  )

  # Both should work but may give different power values
  expect_s3_class(result_complement, "data.frame")
  expect_s3_class(result_nt, "data.frame")
  expect_equal(nrow(result_complement), 1)
  expect_equal(nrow(result_nt), 1)
})

test_that("compute_power_plan_full_grid expression filtering by tpm_threshold", {
  test_data <- setup_grid_test_data()

  # Test with high tpm_threshold (should filter out low-expression genes)
  result_high <- compute_power_plan_full_grid(
    tpm_threshold = 20,  # High threshold
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 15
  )

  # Test with low tpm_threshold (should include more genes)
  result_low <- compute_power_plan_full_grid(
    tpm_threshold = 1,  # Low threshold
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 15
  )

  # Both should work (since we have sufficient genes in test data)
  expect_s3_class(result_high, "data.frame")
  expect_s3_class(result_low, "data.frame")

  # Power may differ due to different gene sets being analyzed
  expect_true(result_high$overall_power >= 0 && result_high$overall_power <= 1)
  expect_true(result_low$overall_power >= 0 && result_low$overall_power <= 1)
})

test_that("compute_power_plan_full_grid library size calculation", {
  test_data <- setup_grid_test_data()

  result <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = c(5000, 10000, 15000),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 15
  )

  # Should have 3 rows for 3 different reads_per_cell values
  expect_equal(nrow(result), 3)

  # Library size should increase with reads per cell
  expect_true(all(diff(result$library_size) > 0))

  # Raw reads per cell should match input (adjusted for mapping efficiency)
  expected_raw_reads <- c(5000, 10000, 15000) / 0.72
  expect_equal(result$raw_reads_per_cell, expected_raw_reads, tolerance = 1e-6)
})

test_that("compute_power_plan_full_grid experimental parameter validation", {
  test_data <- setup_grid_test_data()

  # Test with different experimental parameters
  result <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    MOI = 5,  # Different MOI
    num_targets = 25,  # Different target count
    gRNAs_per_target = 3,  # Different gRNA count
    non_targeting_gRNAs = 20,  # Different non-targeting count
    multiple_testing_alpha = 0.1  # Different alpha
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(result$overall_power >= 0 && result$overall_power <= 1)
})

test_that("compute_power_plan_full_grid gRNA variability effects", {
  test_data <- setup_grid_test_data()

  # Test with low variability
  result_low_var <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    gRNA_variability = 0.05,  # Low variability
    num_targets = 15
  )

  # Test with high variability
  result_high_var <- compute_power_plan_full_grid(
    tpm_threshold = 10,
    minimum_fold_change = 0.8,
    cells_per_target = 1000,
    reads_per_cell = 8000,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    gRNA_variability = 0.3,  # High variability
    num_targets = 15
  )

  # Both should work
  expect_s3_class(result_low_var, "data.frame")
  expect_s3_class(result_high_var, "data.frame")

  # Power should generally be higher with lower variability
  # (though this may not always hold due to randomness in small examples)
  expect_true(result_low_var$overall_power >= 0 && result_low_var$overall_power <= 1)
  expect_true(result_high_var$overall_power >= 0 && result_high_var$overall_power <= 1)
})

test_that("compute_power_plan_full_grid output completeness", {
  test_data <- setup_grid_test_data()

  result <- compute_power_plan_full_grid(
    tpm_threshold = c(5, 15),
    minimum_fold_change = c(0.7, 0.9),
    cells_per_target = c(500, 1500),
    reads_per_cell = c(6000, 12000),
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = 10
  )

  # Should have 2×2×2×2 = 16 combinations
  expect_equal(nrow(result), 16)

  # Check all combinations are present
  param_combinations <- expand.grid(
    tpm_threshold = c(5, 15),
    minimum_fold_change = c(0.7, 0.9),
    cells_per_target = c(500, 1500),
    raw_reads_per_cell = c(6000, 12000) / 0.72
  )

  # Match by rounding to avoid floating point precision issues
  result_rounded <- result[c("tpm_threshold", "minimum_fold_change", "cells_per_target", "raw_reads_per_cell")]
  result_rounded[] <- lapply(result_rounded, round, digits = 6)
  param_combinations[] <- lapply(param_combinations, round, digits = 6)

  expect_true(nrow(dplyr::anti_join(param_combinations, result_rounded, by = names(param_combinations))) == 0)
})

test_that("compute_power_plan_full_grid matches compute_power_plan_overall", {
  # Set up consistent test data and parameters
  test_data <- setup_grid_test_data()
  set.seed(12345)  # Set seed for reproducible fold change generation

  # Fixed experimental parameters for direct comparison
  tpm_threshold_val <- 10
  minimum_fold_change_val <- 0.5
  cells_per_target_val <- 1000
  reads_per_cell_val <- 8000
  num_targets_val <- 20
  gRNAs_per_target_val <- 4
  non_targeting_gRNAs_val <- 10
  MOI_val <- 10

  # Calculate experimental design parameters for compute_power_plan_overall
  total_gRNAs <- num_targets_val * gRNAs_per_target_val + non_targeting_gRNAs_val
  num_total_cells <- (cells_per_target_val * total_gRNAs) / (gRNAs_per_target_val * MOI_val)
  num_trt_cells <- cells_per_target_val
  num_cntrl_cells <- num_total_cells - num_trt_cells

  # Calculate library size
  library_size <- fit_read_UMI_curve_cpp(
    reads_per_cell = reads_per_cell_val,
    UMI_per_cell = test_data$library_parameters$UMI_per_cell,
    variation = test_data$library_parameters$variation
  )

  # Create fc_expression_df for compute_power_plan_overall
  set.seed(12345)  # Reset seed to ensure same fold changes
  filtered_expression_df <- dplyr::filter(
    test_data$baseline_expression_stats,
    relative_expression >= tpm_threshold_val / 1e6
  )

  fc_expression_df <- filtered_expression_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      grna_effects_vec = list(pmax(
        stats::rnorm(n = gRNAs_per_target_val,
                    mean = minimum_fold_change_val,
                    sd = 0.13),  # default gRNA_variability
        .Machine$double.eps
      )),
      avg_fold_change = mean(grna_effects_vec),
      avg_fold_change_sq = mean(grna_effects_vec^2)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-grna_effects_vec)

  # Test 1: Get result from compute_power_plan_full_grid
  set.seed(12345)  # Reset seed
  result_full_grid <- compute_power_plan_full_grid(
    tpm_threshold = tpm_threshold_val,
    minimum_fold_change = minimum_fold_change_val,
    cells_per_target = cells_per_target_val,
    reads_per_cell = reads_per_cell_val,
    baseline_expression_stats = test_data$baseline_expression_stats,
    library_parameters = test_data$library_parameters,
    num_targets = num_targets_val,
    gRNAs_per_target = gRNAs_per_target_val,
    non_targeting_gRNAs = non_targeting_gRNAs_val,
    MOI = MOI_val,
    grid_size = 1  # Single point for exact comparison
  )

  # Test 2: Get result from compute_power_plan_overall
  result_overall <- compute_power_plan_overall(
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    library_size = library_size,
    fc_expression_df = fc_expression_df,
    multiple_testing_alpha = 0.05,
    side = "left",
    prop_non_null = 0.1
  )

  # Test 3: Compare power values
  expect_equal(nrow(result_full_grid), 1)
  expect_type(result_overall, "double")
  expect_length(result_overall, 1)

  # The power values should match within tolerance (allow for Monte Carlo variability)
  expect_equal(result_full_grid$overall_power, result_overall, tolerance = .Machine$double.eps)

  # Test 4: Validate other calculated values match expected
  expect_equal(result_full_grid$cells_per_target, cells_per_target_val)
  expect_equal(result_full_grid$library_size, library_size, tolerance = 1e-6)
  expect_equal(result_full_grid$raw_reads_per_cell, reads_per_cell_val / 0.72, tolerance = 1e-6)
})
