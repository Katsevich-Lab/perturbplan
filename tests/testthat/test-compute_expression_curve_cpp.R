# Test for compute_expression_curve_cpp function
# Validates expression power curve computation against individual calculations

library(testthat)

test_that("compute_expression_curve_cpp matches individual calculations with mock dispersion", {
  
  # Create a simple, predictable dispersion function for testing
  simple_dispersion <- function(relative_expr) {
    # Simple linear relationship: size = 0.1 + 10 * relative_expr
    return(0.1 + 10 * relative_expr)
  }
  
  # Test data with fixed fold changes
  fc_expression_df <- data.frame(
    fold_change = c(0.5, 1.0, 2.0, 0.8, 1.5)
  )
  
  # Simple expression grid
  expr_grid <- c(1e-5, 5e-5, 1e-4)
  library_size <- 10000
  num_trt_cells <- 400
  num_cntrl_cells <- 600
  side <- "both"
  cutoff <- 0.05
  
  # Get curve result
  curve_result <- compute_expression_curve_cpp(
    expr_output_grid = expr_grid,
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    expression_dispersion_curve = simple_dispersion,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    side = side,
    cutoff = cutoff
  )
  
  # Verify output structure
  expect_equal(nrow(curve_result), length(expr_grid))
  expect_equal(ncol(curve_result), 2)
  expect_true(all(c("relative_expression", "power") %in% names(curve_result)))
  expect_equal(curve_result$relative_expression, expr_grid)
  expect_true(all(curve_result$power >= 0 & curve_result$power <= 1))
  expect_true(all(is.finite(curve_result$power)))
  
  # Manually compute each point and compare
  tolerance <- 1e-12  # Very tight tolerance for exact consistency
  
  for (i in 1:length(expr_grid)) {
    expr_val <- expr_grid[i]
    expr_mean <- library_size * expr_val
    expr_size <- simple_dispersion(expr_val)  # We control this!
    
    # Compute test statistics manually for this expression across all fold changes
    expr_means <- numeric(nrow(fc_expression_df))
    expr_sds <- numeric(nrow(fc_expression_df))
    
    for (j in 1:nrow(fc_expression_df)) {
      individual_result <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = fc_expression_df$fold_change[j],
        expression_mean = expr_mean,
        expression_size = expr_size,
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = num_trt_cells
      )
      expr_means[j] <- individual_result$mean
      expr_sds[j] <- individual_result$sd
    }
    
    # Compute rejection probabilities and average power
    expr_powers <- rejection_computation_cpp(expr_means, expr_sds, side, cutoff)
    manual_power <- mean(expr_powers)
    
    # Compare with curve result
    expect_equal(curve_result$power[i], manual_power, tolerance = tolerance,
                info = paste("Power mismatch at expression =", expr_val, "- curve:", 
                           round(curve_result$power[i], 6), "manual:", round(manual_power, 6)))
  }
})

test_that("compute_expression_curve_cpp matches individual calculations with real Gasperini dispersion", {
  
  # Load real dispersion function from Gasperini data
  gas_data <- readRDS(system.file('extdata/baseline_expression', 'Gasperini_expression.rds', package = 'perturbplan'))
  real_dispersion <- gas_data$expression_dispersion_curve
  
  # Test data with fixed fold changes
  fc_expression_df <- data.frame(
    fold_change = c(0.5, 1.0, 2.0, 0.8, 1.5)
  )
  
  # Expression grid covering realistic range
  expr_grid <- c(1e-5, 5e-5, 1e-4)
  library_size <- 10000
  num_trt_cells <- 400
  num_cntrl_cells <- 600
  side <- "both"
  cutoff <- 0.05
  
  # Get curve result
  curve_result <- compute_expression_curve_cpp(
    expr_output_grid = expr_grid,
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    expression_dispersion_curve = real_dispersion,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    side = side,
    cutoff = cutoff
  )
  
  # Verify output structure
  expect_equal(nrow(curve_result), length(expr_grid))
  expect_equal(ncol(curve_result), 2)
  expect_true(all(c("relative_expression", "power") %in% names(curve_result)))
  expect_equal(curve_result$relative_expression, expr_grid)
  expect_true(all(curve_result$power >= 0 & curve_result$power <= 1))
  expect_true(all(is.finite(curve_result$power)))
  
  # Manually compute each point and compare
  tolerance <- 1e-12  # Very tight tolerance for exact consistency
  
  for (i in 1:length(expr_grid)) {
    expr_val <- expr_grid[i]
    expr_mean <- library_size * expr_val
    expr_size <- real_dispersion(expr_val)  # Use real Gasperini dispersion!
    
    # Compute test statistics manually for this expression across all fold changes
    expr_means <- numeric(nrow(fc_expression_df))
    expr_sds <- numeric(nrow(fc_expression_df))
    
    for (j in 1:nrow(fc_expression_df)) {
      individual_result <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = fc_expression_df$fold_change[j],
        expression_mean = expr_mean,
        expression_size = expr_size,
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = num_trt_cells
      )
      expr_means[j] <- individual_result$mean
      expr_sds[j] <- individual_result$sd
    }
    
    # Compute rejection probabilities and average power
    expr_powers <- rejection_computation_cpp(expr_means, expr_sds, side, cutoff)
    manual_power <- mean(expr_powers)
    
    # Compare with curve result
    expect_equal(curve_result$power[i], manual_power, tolerance = tolerance,
                info = paste("Power mismatch at expression =", expr_val, "- curve:", 
                           round(curve_result$power[i], 6), "manual:", round(manual_power, 6)))
  }
  
  # Test biological expectation: power should generally increase with expression
  # (Allow some noise due to dispersion curve complexity)
  if (length(expr_grid) > 1) {
    power_trend <- diff(curve_result$power)
    positive_trend_fraction <- sum(power_trend > 0) / length(power_trend)
    expect_true(positive_trend_fraction >= 0.5,
               info = paste("Expected mostly increasing power trend, got", 
                           round(positive_trend_fraction, 2), "positive fraction"))
  }
  
  # Test different experimental parameters
  test_params <- list(
    list(lib_size = 5000, trt_cells = 200, ctrl_cells = 400, side = "right"),
    list(lib_size = 15000, trt_cells = 500, ctrl_cells = 1000, side = "left")
  )
  
  for (params in test_params) {
    # Use smaller sample for efficiency
    small_expr_grid <- expr_grid[1:2]  # Just two points
    small_fc_df <- data.frame(fold_change = fc_expression_df$fold_change[1:3])  # Just three fold changes
    
    param_curve_result <- compute_expression_curve_cpp(
      expr_output_grid = small_expr_grid,
      fc_expression_df = small_fc_df,
      library_size = params$lib_size,
      expression_dispersion_curve = real_dispersion,
      num_trt_cells = params$trt_cells,
      num_cntrl_cells = params$ctrl_cells,
      side = params$side,
      cutoff = cutoff
    )
    
    # Verify first point manually for these parameters
    expr_val <- small_expr_grid[1]
    expr_mean <- params$lib_size * expr_val
    expr_size <- real_dispersion(expr_val)
    
    expr_means <- numeric(nrow(small_fc_df))
    expr_sds <- numeric(nrow(small_fc_df))
    
    for (j in 1:nrow(small_fc_df)) {
      individual_result <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = small_fc_df$fold_change[j],
        expression_mean = expr_mean,
        expression_size = expr_size,
        num_trt_cells = params$trt_cells,
        num_cntrl_cells = params$ctrl_cells,
        num_cells = params$trt_cells
      )
      expr_means[j] <- individual_result$mean
      expr_sds[j] <- individual_result$sd
    }
    
    expr_powers <- rejection_computation_cpp(expr_means, expr_sds, params$side, cutoff)
    manual_power <- mean(expr_powers)
    
    expect_equal(param_curve_result$power[1], manual_power, tolerance = tolerance,
                info = paste("Parameter test failed for lib_size:", params$lib_size,
                           "cells:", params$trt_cells, params$ctrl_cells, "side:", params$side))
  }
  
  # Print summary for manual inspection
  cat("\n--- Expression Curve Test Summary ---")
  cat("\nExpression levels tested:", paste(expr_grid, collapse = ", "))
  cat("\nPower values:", paste(round(curve_result$power, 4), collapse = ", "))
  cat("\nFold change samples:", nrow(fc_expression_df))
  cat("\nUsing real Gasperini dispersion curve")
  cat("\nAll consistency checks passed with tolerance:", tolerance)
})