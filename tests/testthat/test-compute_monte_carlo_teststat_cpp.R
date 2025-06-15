# Test for compute_monte_carlo_teststat_cpp function
# Validates Monte Carlo test statistic computation against individual calls

library(testthat)

test_that("compute_monte_carlo_teststat_cpp is consistent with individual function calls", {
  
  # Create test data with diverse scenarios
  fc_expression_df <- data.frame(
    fold_change = c(0.5, 0.8, 1.0, 1.2, 2.0, 3.0),
    relative_expression = c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4),
    expression_size = c(0.1, 0.5, 1.0, 2.0, 5.0, 10.0)
  )
  
  # Test parameters covering realistic experimental ranges
  library_size <- 10000
  num_trt_cells <- 500
  num_cntrl_cells <- 1000
  
  # Call batch function
  batch_result <- compute_monte_carlo_teststat_cpp(
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells
  )
  
  # Verify output structure
  expect_equal(length(batch_result$means), nrow(fc_expression_df))
  expect_equal(length(batch_result$sds), nrow(fc_expression_df))
  expect_true(all(is.finite(batch_result$means)))
  expect_true(all(is.finite(batch_result$sds)))
  expect_true(all(batch_result$sds > 0))  # Standard deviations should be positive
  
  # Compare each scenario with individual function calls
  tolerance <- 1e-12  # Very tight tolerance for exact consistency
  
  for (i in 1:nrow(fc_expression_df)) {
    # Compute expression mean for this scenario
    expression_mean <- fc_expression_df$relative_expression[i] * library_size
    
    # Call individual function
    individual_result <- compute_distribution_teststat_fixed_es_cpp(
      fold_change = fc_expression_df$fold_change[i],
      expression_mean = expression_mean,
      expression_size = fc_expression_df$expression_size[i],
      num_trt_cells = num_trt_cells,
      num_cntrl_cells = num_cntrl_cells,
      num_cells = num_trt_cells  # Single gRNA case
    )
    
    # Compare batch result with individual result
    expect_equal(batch_result$means[i], individual_result$mean, tolerance = tolerance)
    expect_equal(batch_result$sds[i], individual_result$sd, tolerance = tolerance)
  }
  
  # Test different experimental parameter combinations
  test_params <- list(
    list(lib_size = 5000, trt_cells = 300, ctrl_cells = 700),
    list(lib_size = 15000, trt_cells = 800, ctrl_cells = 1200),
    list(lib_size = 20000, trt_cells = 1000, ctrl_cells = 2000)
  )
  
  for (params in test_params) {
    batch_result_param <- compute_monte_carlo_teststat_cpp(
      fc_expression_df = fc_expression_df,
      library_size = params$lib_size,
      num_trt_cells = params$trt_cells,
      num_cntrl_cells = params$ctrl_cells
    )
    
    # Verify consistency for each scenario with these parameters
    for (i in 1:min(3, nrow(fc_expression_df))) {  # Test first 3 scenarios for efficiency
      expression_mean <- fc_expression_df$relative_expression[i] * params$lib_size
      
      individual_result_param <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = fc_expression_df$fold_change[i],
        expression_mean = expression_mean,
        expression_size = fc_expression_df$expression_size[i],
        num_trt_cells = params$trt_cells,
        num_cntrl_cells = params$ctrl_cells,
        num_cells = params$trt_cells
      )
      
      expect_equal(batch_result_param$means[i], individual_result_param$mean, 
                  tolerance = tolerance)
      expect_equal(batch_result_param$sds[i], individual_result_param$sd, 
                  tolerance = tolerance)
    }
  }
})