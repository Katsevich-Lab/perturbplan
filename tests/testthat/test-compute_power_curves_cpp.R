library(testthat)
library(dplyr)

# Helper function to create test data
create_test_data <- function(n_samples = 50) {
  set.seed(123)  # For reproducible tests
  
  fc_expression_df <- data.frame(
    fold_change = rnorm(n_samples, mean = 0.8, sd = 0.1),
    relative_expression = runif(n_samples, min = 0.0005, max = 0.005),
    expression_size = runif(n_samples, min = 8, max = 15)
  )
  
  # Simple expression dispersion curve for testing
  expression_dispersion_curve <- function(x) 10 + 5/x
  
  list(
    fc_expression_df = fc_expression_df,
    expression_dispersion_curve = expression_dispersion_curve,
    fc_output_grid = c(0.6, 0.8, 1.0, 1.2, 1.4),
    expr_output_grid = c(0.001, 0.002, 0.003, 0.004, 0.005),
    library_size = 10000,
    num_trt_cells = 100,
    num_cntrl_cells = 200,
    side = "left",
    cutoff = 0.05
  )
}

test_that("compute_fc_curve_cpp gives identical results to R implementation", {
  
  test_data <- create_test_data()
  
  # Call R function
  result_r <- compute_fc_curve(
    fc_output_grid = test_data$fc_output_grid,
    fc_expression_df = test_data$fc_expression_df,
    library_size = test_data$library_size,
    num_trt_cells = test_data$num_trt_cells,
    num_cntrl_cells = test_data$num_cntrl_cells,
    side = test_data$side,
    cutoff = test_data$cutoff
  )
  
  # Call C++ function
  result_cpp <- compute_fc_curve_cpp(
    fc_output_grid = test_data$fc_output_grid,
    fc_expression_df = test_data$fc_expression_df,
    library_size = test_data$library_size,
    num_trt_cells = test_data$num_trt_cells,
    num_cntrl_cells = test_data$num_cntrl_cells,
    side = test_data$side,
    cutoff = test_data$cutoff
  )
  
  # Test that results are nearly identical
  expect_equal(result_r$fold_change, result_cpp$fold_change, tolerance = 1e-12)
  expect_equal(result_r$power, result_cpp$power, tolerance = 1e-12)
  
  # Test data frame structure
  expect_s3_class(result_cpp, "data.frame")
  expect_named(result_cpp, c("fold_change", "power"))
  expect_equal(nrow(result_cpp), length(test_data$fc_output_grid))
})

test_that("compute_expression_curve_cpp gives identical results to R implementation", {
  
  test_data <- create_test_data()
  
  # Call R function
  result_r <- compute_expression_curve(
    expr_output_grid = test_data$expr_output_grid,
    fc_expression_df = test_data$fc_expression_df,
    library_size = test_data$library_size,
    expression_dispersion_curve = test_data$expression_dispersion_curve,
    num_trt_cells = test_data$num_trt_cells,
    num_cntrl_cells = test_data$num_cntrl_cells,
    side = test_data$side,
    cutoff = test_data$cutoff
  )
  
  # Call C++ function
  result_cpp <- compute_expression_curve_cpp(
    expr_output_grid = test_data$expr_output_grid,
    fc_expression_df = test_data$fc_expression_df,
    library_size = test_data$library_size,
    expression_dispersion_curve = test_data$expression_dispersion_curve,
    num_trt_cells = test_data$num_trt_cells,
    num_cntrl_cells = test_data$num_cntrl_cells,
    side = test_data$side,
    cutoff = test_data$cutoff
  )
  
  # Test that results are nearly identical
  expect_equal(result_r$relative_expression, result_cpp$relative_expression, tolerance = 1e-12)
  expect_equal(result_r$power, result_cpp$power, tolerance = 1e-12)
  
  # Test data frame structure
  expect_s3_class(result_cpp, "data.frame")
  expect_named(result_cpp, c("relative_expression", "power"))
  expect_equal(nrow(result_cpp), length(test_data$expr_output_grid))
})

test_that("C++ functions handle different test sides correctly", {
  
  test_data <- create_test_data()
  sides <- c("left", "right", "both")
  
  for (side in sides) {
    # Test FC curve
    result_r_fc <- compute_fc_curve(
      fc_output_grid = test_data$fc_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = side,
      cutoff = test_data$cutoff
    )
    
    result_cpp_fc <- compute_fc_curve_cpp(
      fc_output_grid = test_data$fc_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = side,
      cutoff = test_data$cutoff
    )
    
    expect_equal(result_r_fc$power, result_cpp_fc$power, tolerance = 1e-12)
    
    # Test expression curve
    result_r_expr <- compute_expression_curve(
      expr_output_grid = test_data$expr_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      expression_dispersion_curve = test_data$expression_dispersion_curve,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = side,
      cutoff = test_data$cutoff
    )
    
    result_cpp_expr <- compute_expression_curve_cpp(
      expr_output_grid = test_data$expr_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      expression_dispersion_curve = test_data$expression_dispersion_curve,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = side,
      cutoff = test_data$cutoff
    )
    
    expect_equal(result_r_expr$power, result_cpp_expr$power, tolerance = 1e-12)
  }
})

test_that("C++ functions handle different cutoff values correctly", {
  
  test_data <- create_test_data()
  cutoffs <- c(0.01, 0.05, 0.1, 0.2)
  
  for (cutoff in cutoffs) {
    # Test FC curve
    result_r_fc <- compute_fc_curve(
      fc_output_grid = test_data$fc_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = test_data$side,
      cutoff = cutoff
    )
    
    result_cpp_fc <- compute_fc_curve_cpp(
      fc_output_grid = test_data$fc_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = test_data$side,
      cutoff = cutoff
    )
    
    expect_equal(result_r_fc$power, result_cpp_fc$power, tolerance = 1e-12)
    
    # Test expression curve  
    result_r_expr <- compute_expression_curve(
      expr_output_grid = test_data$expr_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      expression_dispersion_curve = test_data$expression_dispersion_curve,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = test_data$side,
      cutoff = cutoff
    )
    
    result_cpp_expr <- compute_expression_curve_cpp(
      expr_output_grid = test_data$expr_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      expression_dispersion_curve = test_data$expression_dispersion_curve,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = test_data$side,
      cutoff = cutoff
    )
    
    expect_equal(result_r_expr$power, result_cpp_expr$power, tolerance = 1e-12)
  }
})

test_that("C++ functions handle edge cases correctly", {
  
  # Test with minimal data
  minimal_data <- data.frame(
    fold_change = 0.8,
    relative_expression = 0.001,
    expression_size = 10
  )
  
  # Test single grid point
  result_r_fc <- compute_fc_curve(
    fc_output_grid = 0.9,
    fc_expression_df = minimal_data,
    library_size = 5000,
    num_trt_cells = 50,
    num_cntrl_cells = 100,
    side = "left",
    cutoff = 0.05
  )
  
  result_cpp_fc <- compute_fc_curve_cpp(
    fc_output_grid = 0.9,
    fc_expression_df = minimal_data,
    library_size = 5000,
    num_trt_cells = 50,
    num_cntrl_cells = 100,
    side = "left",
    cutoff = 0.05
  )
  
  expect_equal(result_r_fc$power, result_cpp_fc$power, tolerance = 1e-12)
  expect_equal(nrow(result_cpp_fc), 1)
  
  # Test expression curve with single point
  expr_curve_func <- function(x) 12
  
  result_r_expr <- compute_expression_curve(
    expr_output_grid = 0.002,
    fc_expression_df = minimal_data,
    library_size = 5000,
    expression_dispersion_curve = expr_curve_func,
    num_trt_cells = 50,
    num_cntrl_cells = 100,
    side = "left",
    cutoff = 0.05
  )
  
  result_cpp_expr <- compute_expression_curve_cpp(
    expr_output_grid = 0.002,
    fc_expression_df = minimal_data,
    library_size = 5000,
    expression_dispersion_curve = expr_curve_func,
    num_trt_cells = 50,
    num_cntrl_cells = 100,
    side = "left",
    cutoff = 0.05
  )
  
  expect_equal(result_r_expr$power, result_cpp_expr$power, tolerance = 1e-12)
  expect_equal(nrow(result_cpp_expr), 1)
})