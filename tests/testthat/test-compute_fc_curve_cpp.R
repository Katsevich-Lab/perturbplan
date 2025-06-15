# Test for compute_fc_curve_cpp function
# Validates fold change power curve computation against individual calculations

library(testthat)

test_that("compute_fc_curve_cpp matches individual power calculations", {
  
  # Use Gasperini data for realistic testing
  gas_data <- readRDS(system.file('extdata/baseline_expression', 'Gasperini_expression.rds', package = 'perturbplan'))
  baseline_df <- gas_data$baseline_expression
  
  # Sample representative genes across expression spectrum
  set.seed(42)
  n_samples <- 10  # Small for efficient testing
  sampled_indices <- sample(nrow(baseline_df), n_samples)
  
  fc_expression_df <- data.frame(
    relative_expression = baseline_df$relative_expression[sampled_indices],
    expression_size = baseline_df$expression_size[sampled_indices]
  )
  
  # Test with small grid for efficiency but cover key points
  fc_grid <- c(0.5, 1.0, 2.0)
  library_size <- 10000
  num_trt_cells <- 300
  num_cntrl_cells <- 700
  side <- "both"
  cutoff <- 0.05
  
  # Get curve result
  curve_result <- compute_fc_curve_cpp(
    fc_output_grid = fc_grid,
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    side = side,
    cutoff = cutoff
  )
  
  # Verify output structure
  expect_equal(nrow(curve_result), length(fc_grid))
  expect_equal(ncol(curve_result), 2)
  expect_true(all(c("fold_change", "power") %in% names(curve_result)))
  expect_equal(curve_result$fold_change, fc_grid)
  expect_true(all(curve_result$power >= 0 & curve_result$power <= 1))
  expect_true(all(is.finite(curve_result$power)))
  
  # Manually compute each point and compare
  tolerance <- 1e-12  # Very tight tolerance for exact consistency
  
  for (i in 1:length(fc_grid)) {
    fc_val <- fc_grid[i]
    
    # Manually compute what the curve function should do for this FC
    # This replicates the exact logic inside compute_fc_curve_cpp
    mc_expression_means <- library_size * fc_expression_df$relative_expression
    
    fc_means <- numeric(n_samples)
    fc_sds <- numeric(n_samples)
    
    # Compute test statistics for this FC across all Monte Carlo expression samples
    for (j in 1:n_samples) {
      individual_result <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = fc_val,
        expression_mean = mc_expression_means[j],
        expression_size = fc_expression_df$expression_size[j],
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = num_trt_cells  # Single gRNA case
      )
      fc_means[j] <- individual_result$mean
      fc_sds[j] <- individual_result$sd
    }
    
    # Compute rejection probabilities and average power
    fc_powers <- rejection_computation_cpp(fc_means, fc_sds, side, cutoff)
    manual_power <- mean(fc_powers)
    
    # Compare with curve result
    expect_equal(curve_result$power[i], manual_power, tolerance = tolerance,
                info = paste("Power mismatch at FC =", fc_val, "- curve:", 
                           round(curve_result$power[i], 6), "manual:", round(manual_power, 6)))
  }
  
  # Test with different experimental parameters to ensure consistency
  test_params <- list(
    list(lib_size = 5000, trt_cells = 200, ctrl_cells = 400, side = "right"),
    list(lib_size = 15000, trt_cells = 500, ctrl_cells = 1000, side = "left"),
    list(lib_size = 20000, trt_cells = 800, ctrl_cells = 800, side = "both")
  )
  
  for (params in test_params) {
    # Use smaller sample for efficiency in parameter testing
    small_fc_expression_df <- fc_expression_df[1:5, ]
    small_fc_grid <- c(0.8, 1.5)  # Just two points
    
    # Get curve result with these parameters
    param_curve_result <- compute_fc_curve_cpp(
      fc_output_grid = small_fc_grid,
      fc_expression_df = small_fc_expression_df,
      library_size = params$lib_size,
      num_trt_cells = params$trt_cells,
      num_cntrl_cells = params$ctrl_cells,
      side = params$side,
      cutoff = cutoff
    )
    
    # Manually verify first point for these parameters
    fc_val <- small_fc_grid[1]
    mc_expression_means <- params$lib_size * small_fc_expression_df$relative_expression
    
    fc_means <- numeric(nrow(small_fc_expression_df))
    fc_sds <- numeric(nrow(small_fc_expression_df))
    
    for (j in 1:nrow(small_fc_expression_df)) {
      individual_result <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = fc_val,
        expression_mean = mc_expression_means[j],
        expression_size = small_fc_expression_df$expression_size[j],
        num_trt_cells = params$trt_cells,
        num_cntrl_cells = params$ctrl_cells,
        num_cells = params$trt_cells
      )
      fc_means[j] <- individual_result$mean
      fc_sds[j] <- individual_result$sd
    }
    
    fc_powers <- rejection_computation_cpp(fc_means, fc_sds, params$side, cutoff)
    manual_power <- mean(fc_powers)
    
    expect_equal(param_curve_result$power[1], manual_power, tolerance = tolerance,
                info = paste("Parameter test failed for lib_size:", params$lib_size,
                           "cells:", params$trt_cells, params$ctrl_cells, "side:", params$side))
  }
  
  # Test edge case: single fold change value
  single_fc_result <- compute_fc_curve_cpp(
    fc_output_grid = c(1.5),
    fc_expression_df = fc_expression_df[1:3, ],  # Small subset
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    side = side,
    cutoff = cutoff
  )
  
  expect_equal(nrow(single_fc_result), 1)
  expect_equal(single_fc_result$fold_change[1], 1.5)
  expect_true(single_fc_result$power[1] >= 0 && single_fc_result$power[1] <= 1)
  expect_true(is.finite(single_fc_result$power[1]))
  
  # Test edge case: single expression sample
  single_expr_df <- data.frame(
    relative_expression = baseline_df$relative_expression[1],
    expression_size = baseline_df$expression_size[1]
  )
  
  single_expr_result <- compute_fc_curve_cpp(
    fc_output_grid = fc_grid,
    fc_expression_df = single_expr_df,
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    side = side,
    cutoff = cutoff
  )
  
  expect_equal(nrow(single_expr_result), length(fc_grid))
  expect_true(all(is.finite(single_expr_result$power)))
  expect_true(all(single_expr_result$power >= 0 & single_expr_result$power <= 1))
  
  # Print summary for manual inspection
  cat("\n--- FC Curve Test Summary ---")
  cat("\nFold changes tested:", paste(fc_grid, collapse = ", "))
  cat("\nPower values:", paste(round(curve_result$power, 4), collapse = ", "))
  cat("\nExpression samples:", n_samples)
  cat("\nAll consistency checks passed with tolerance:", tolerance)
})