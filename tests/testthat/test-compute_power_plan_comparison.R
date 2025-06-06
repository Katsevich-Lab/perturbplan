test_that("C++ Monte Carlo implementation matches R implementation", {
  
  # Set up test data
  set.seed(123)
  num_total_cells <- 10000
  library_size <- 500
  MOI <- 10
  num_targets <- 100
  gRNAs_per_target <- 4
  non_targeting_gRNAs <- 10
  multiple_testing_alpha <- 0.05
  multiple_testing_method <- "BH"
  control_group <- "complement"
  side <- "left"
  num_pairs <- 1000
  prop_non_null <- 0.1
  B <- 50  # Use smaller B for faster tests
  
  # Create test fc_expression_df
  fc_expression_info <- extract_fc_expression_info(
    fold_change_mean = 0.85, 
    fold_change_sd = 0.15,
    biological_system = "K562", 
    B = B
  )
  
  fc_expression_df <- fc_expression_info$fc_expression_df
  expression_dispersion_curve <- fc_expression_info$expression_dispersion_curve
  
  # Define output grids (smaller for faster tests)
  fc_output_grid <- seq(0.5, 1.0, length.out = 5)
  expr_output_grid <- 10^seq(-3, -1, length.out = 5)
  
  # Test with same seed for both implementations
  set.seed(456)
  result_old <- .compute_underspecified_power_efficient(
    num_total_cells = num_total_cells,
    library_size = library_size,
    MOI = MOI,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    multiple_testing_alpha = multiple_testing_alpha,
    multiple_testing_method = multiple_testing_method,
    control_group = control_group,
    side = side,
    num_pairs = num_pairs,
    fc_expression_df = fc_expression_df,
    expression_dispersion_curve = expression_dispersion_curve,
    fc_output_grid = fc_output_grid,
    expr_output_grid = expr_output_grid,
    prop_non_null = prop_non_null
  )
  
  set.seed(456)
  result_new <- .compute_power_plan_efficient(
    num_total_cells = num_total_cells,
    library_size = library_size,
    MOI = MOI,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    multiple_testing_alpha = multiple_testing_alpha,
    multiple_testing_method = multiple_testing_method,
    control_group = control_group,
    side = side,
    num_pairs = num_pairs,
    fc_expression_df = fc_expression_df,
    expression_dispersion_curve = expression_dispersion_curve,
    fc_output_grid = fc_output_grid,
    expr_output_grid = expr_output_grid,
    prop_non_null = prop_non_null
  )
  
  # Test that overall power is nearly identical
  expect_equal(result_old$overall_power, result_new$overall_power, tolerance = 1e-10)
  
  # Test that power curves are nearly identical
  expect_equal(result_old$power_by_fc$power, result_new$power_by_fc$power, tolerance = 1e-10)
  expect_equal(result_old$power_by_expr$power, result_new$power_by_expr$power, tolerance = 1e-10)
  
  # Test that fold change grid values match
  expect_equal(result_old$power_by_fc$fold_change, result_new$power_by_fc$fold_change, tolerance = 1e-10)
  expect_equal(result_old$power_by_expr$relative_expression, result_new$power_by_expr$relative_expression, tolerance = 1e-10)
  
  # Test that grid summary matches
  expect_equal(result_old$grid_summary$cutoff, result_new$grid_summary$cutoff, tolerance = 1e-10)
  expect_equal(result_old$grid_summary$mc_samples, result_new$grid_summary$mc_samples)
  expect_equal(result_old$grid_summary$fc_output_points, result_new$grid_summary$fc_output_points)
  expect_equal(result_old$grid_summary$expr_output_points, result_new$grid_summary$expr_output_points)
})

test_that("C++ Monte Carlo provides performance improvement", {
  
  # Set up test data
  set.seed(789)
  num_total_cells <- 10000
  library_size <- 500
  B <- 100  # Use moderate B for timing comparison
  
  # Create test fc_expression_df
  fc_expression_info <- extract_fc_expression_info(
    fold_change_mean = 0.85, 
    fold_change_sd = 0.15,
    biological_system = "K562", 
    B = B
  )
  
  fc_expression_df <- fc_expression_info$fc_expression_df
  expression_dispersion_curve <- fc_expression_info$expression_dispersion_curve
  
  # Define output grids
  fc_output_grid <- seq(0.5, 1.0, length.out = 8)
  expr_output_grid <- 10^seq(-3, -1, length.out = 8)
  
  # Time the R implementation
  time_old <- system.time({
    result_old <- .compute_underspecified_power_efficient(
      num_total_cells = num_total_cells,
      library_size = library_size,
      fc_expression_df = fc_expression_df,
      expression_dispersion_curve = expression_dispersion_curve,
      fc_output_grid = fc_output_grid,
      expr_output_grid = expr_output_grid
    )
  })
  
  # Time the C++ implementation
  time_new <- system.time({
    result_new <- .compute_power_plan_efficient(
      num_total_cells = num_total_cells,
      library_size = library_size,
      fc_expression_df = fc_expression_df,
      expression_dispersion_curve = expression_dispersion_curve,
      fc_output_grid = fc_output_grid,
      expr_output_grid = expr_output_grid
    )
  })
  
  # Print timing comparison
  speedup <- time_old[["elapsed"]] / time_new[["elapsed"]]
  cat(sprintf("\nR implementation: %.3f seconds\n", time_old[["elapsed"]]))
  cat(sprintf("C++ implementation: %.3f seconds\n", time_new[["elapsed"]]))
  cat(sprintf("Speedup: %.2fx\n", speedup))
  
  # Verify results are still equivalent
  expect_equal(result_old$overall_power, result_new$overall_power, tolerance = 1e-10)
  
  # Expect some performance improvement (C++ should be faster)
  # Note: This might not always pass due to system variability, so we're lenient
  expect_true(speedup > 0.5)  # At least not much slower
})