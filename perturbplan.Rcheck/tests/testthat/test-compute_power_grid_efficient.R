library(testthat)
library(dplyr)

test_that("compute_test_stat_clean and compute_distribution_teststat_fixed_es_cpp give identical results for single gRNA case", {
  
  # Test parameters for single gRNA case
  num_trt_cells <- 100
  num_cntrl_cells <- 150  
  expression_mean <- 5.0
  expression_size <- 2.5
  fold_change_mean <- 0.8
  
  # Call R function
  result_r <- compute_test_stat_clean(
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    expression_mean = expression_mean,
    expression_size = expression_size,
    fold_change_mean = fold_change_mean
  )
  
  # Call C++ function  
  result_cpp <- compute_distribution_teststat_fixed_es_cpp(
    fold_change = fold_change_mean,
    expression_mean = expression_mean,
    expression_size = expression_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    num_cells = num_trt_cells  # For single gRNA case, all cells have same fold change
  )
  
  # Test that results are nearly identical (allowing for small numerical differences)
  expect_equal(result_r[["mean"]], result_cpp[["mean"]], tolerance = 1e-12)
  expect_equal(result_r[["sd"]], result_cpp[["sd"]], tolerance = 1e-12)
})

test_that("C++ function handles vector inputs correctly", {
  
  # Test parameters for multi-gRNA case
  num_trt_cells <- 300  # Total across all gRNAs
  num_cntrl_cells <- 150  
  expression_mean <- 5.0
  expression_size <- 2.5
  
  # Multi-gRNA case
  fold_change <- c(0.8, 0.85, 0.9)  # 3 gRNAs
  num_cells <- c(100, 100, 100)     # Equal cells per gRNA
  
  result_multi <- compute_distribution_teststat_fixed_es_cpp(
    fold_change = fold_change,
    expression_mean = expression_mean,
    expression_size = expression_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    num_cells = num_cells
  )
  
  # Test that result has expected structure
  expect_true(is.list(result_multi))
  expect_true("mean" %in% names(result_multi))
  expect_true("sd" %in% names(result_multi))
  expect_true(is.numeric(result_multi$mean))
  expect_true(is.numeric(result_multi$sd))
  expect_true(length(result_multi$mean) == 1)
  expect_true(length(result_multi$sd) == 1)
  
  # Results should be reasonable
  expect_true(is.finite(result_multi$mean))
  expect_true(is.finite(result_multi$sd))
  expect_true(result_multi$sd > 0)
})

test_that("C++ and R functions give different results for multi-gRNA case", {
  
  # Test that the more general C++ function gives different results 
  # when there are multiple gRNAs with different fold changes
  
  num_trt_cells <- 200
  num_cntrl_cells <- 150  
  expression_mean <- 5.0
  expression_size <- 2.5
  
  # Single gRNA equivalent using weighted average
  fold_change_avg <- 0.85
  result_r <- compute_test_stat_clean(
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    expression_mean = expression_mean,
    expression_size = expression_size,
    fold_change_mean = fold_change_avg
  )
  
  # Multi-gRNA case with heterogeneous fold changes  
  fold_change <- c(0.7, 1.0)  # Very different fold changes
  num_cells <- c(100, 100)    # Equal cells per gRNA
  
  result_cpp <- compute_distribution_teststat_fixed_es_cpp(
    fold_change = fold_change,
    expression_mean = expression_mean,
    expression_size = expression_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    num_cells = num_cells
  )
  
  # Results should be different due to variance across gRNAs
  # The mean should be very close (both functions use same weighted average)
  expect_equal(result_r[["mean"]], result_cpp[["mean"]], tolerance = 1e-12)
  # But the standard deviation should be different (C++ accounts for gRNA heterogeneity)
  expect_false(isTRUE(all.equal(result_r[["sd"]], result_cpp[["sd"]], tolerance = 1e-6)))
})

test_that("compute_power_grid_efficient produces identical results to compute_power_grid_separated", {
  
  # Set up simple test case
  cells_reads_df <- data.frame(
    num_total_cells = c(10000, 20000),
    reads_per_cell = c(500, 1000)
  )
  
  # Use small parameters for fast testing
  test_params <- list(
    num_targets = 10,
    gRNAs_per_target = 2,
    non_targeting_gRNAs = 5,
    num_pairs = 100,
    fdr_target = 0.1,
    fc_mean = 0.8,
    fc_sd = 0.1,
    prop_non_null = 0.2,
    B = 50,  # Small for speed
    fc_curve_points = 5,
    expr_curve_points = 5
  )
  
  # Run both implementations with same random seed
  set.seed(42)
  result_original <- compute_power_grid_separated(
    cells_reads_df = cells_reads_df,
    num_targets = test_params$num_targets,
    gRNAs_per_target = test_params$gRNAs_per_target,
    non_targeting_gRNAs = test_params$non_targeting_gRNAs,
    num_pairs = test_params$num_pairs,
    fdr_target = test_params$fdr_target,
    fc_mean = test_params$fc_mean,
    fc_sd = test_params$fc_sd,
    prop_non_null = test_params$prop_non_null,
    B = test_params$B,
    fc_curve_points = test_params$fc_curve_points,
    expr_curve_points = test_params$expr_curve_points
  )
  
  set.seed(42)
  result_efficient <- compute_power_grid_efficient(
    cells_reads_df = cells_reads_df,
    num_targets = test_params$num_targets,
    gRNAs_per_target = test_params$gRNAs_per_target,
    non_targeting_gRNAs = test_params$non_targeting_gRNAs,
    num_pairs = test_params$num_pairs,
    fdr_target = test_params$fdr_target,
    fc_mean = test_params$fc_mean,
    fc_sd = test_params$fc_sd,
    prop_non_null = test_params$prop_non_null,
    B = test_params$B,
    fc_curve_points = test_params$fc_curve_points,
    expr_curve_points = test_params$expr_curve_points
  )
  
  # Test that both results have the same structure
  expect_identical(names(result_original), names(result_efficient))
  expect_identical(nrow(result_original), nrow(result_efficient))
  expect_identical(ncol(result_original), ncol(result_efficient))
  
  # Test that overall power values are very close (allowing for small numerical differences)
  expect_equal(result_original$overall_power, result_efficient$overall_power, tolerance = 1e-10)
  
  # Test that power curves have same structure
  expect_identical(names(result_original$power_by_fc[[1]]), names(result_efficient$power_by_fc[[1]]))
  expect_identical(names(result_original$power_by_expr[[1]]), names(result_efficient$power_by_expr[[1]]))
  
  # Test that fold change grids are identical
  expect_equal(result_original$power_by_fc[[1]]$fold_change, 
               result_efficient$power_by_fc[[1]]$fold_change, tolerance = 1e-14)
  expect_equal(result_original$power_by_fc[[2]]$fold_change, 
               result_efficient$power_by_fc[[2]]$fold_change, tolerance = 1e-14)
  
  # Test that expression grids are identical  
  expect_equal(result_original$power_by_expr[[1]]$relative_expression, 
               result_efficient$power_by_expr[[1]]$relative_expression, tolerance = 1e-14)
  expect_equal(result_original$power_by_expr[[2]]$relative_expression, 
               result_efficient$power_by_expr[[2]]$relative_expression, tolerance = 1e-14)
  
  # Test that power values in curves are very close
  expect_equal(result_original$power_by_fc[[1]]$power, 
               result_efficient$power_by_fc[[1]]$power, tolerance = 1e-10)
  expect_equal(result_original$power_by_fc[[2]]$power, 
               result_efficient$power_by_fc[[2]]$power, tolerance = 1e-10)
  
  expect_equal(result_original$power_by_expr[[1]]$power, 
               result_efficient$power_by_expr[[1]]$power, tolerance = 1e-10)
  expect_equal(result_original$power_by_expr[[2]]$power, 
               result_efficient$power_by_expr[[2]]$power, tolerance = 1e-10)
})

test_that("functions are accessible in package namespace", {
  # Test that the core functions are accessible
  expect_true(exists("compute_test_stat_clean"))
  expect_true(exists("compute_distribution_teststat_fixed_es_cpp"))
})