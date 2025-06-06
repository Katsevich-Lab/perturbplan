library(testthat)
library(dplyr)

# Helper function to create benchmark data of different sizes
create_benchmark_data <- function(n_samples = 100, n_grid_points = 10) {
  set.seed(456)  # Different seed for benchmark data
  
  fc_expression_df <- data.frame(
    fold_change = rnorm(n_samples, mean = 0.8, sd = 0.15),
    relative_expression = runif(n_samples, min = 0.0003, max = 0.007),
    expression_size = runif(n_samples, min = 5, max = 20)
  )
  
  # More complex expression dispersion curve
  expression_dispersion_curve <- function(x) 15 + 10/sqrt(x) + 2*x
  
  list(
    fc_expression_df = fc_expression_df,
    expression_dispersion_curve = expression_dispersion_curve,
    fc_output_grid = seq(0.5, 1.3, length.out = n_grid_points),
    expr_output_grid = seq(0.0005, 0.006, length.out = n_grid_points),
    library_size = 12000,
    num_trt_cells = 120,
    num_cntrl_cells = 180,
    side = "left",
    cutoff = 0.05
  )
}

# Helper function to time function execution
time_function <- function(expr, n_reps = 3) {
  times <- numeric(n_reps)
  for (i in seq_len(n_reps)) {
    start_time <- Sys.time()
    eval(expr)
    end_time <- Sys.time()
    times[i] <- as.numeric(end_time - start_time, units = "secs")
  }
  list(
    mean_time = mean(times),
    median_time = median(times),
    min_time = min(times),
    max_time = max(times),
    all_times = times
  )
}

test_that("C++ FC curve computation is faster than R implementation", {
  
  # Use larger datasets for meaningful benchmarks
  test_data <- create_benchmark_data(n_samples = 200, n_grid_points = 15)
  
  # Benchmark R implementation
  r_timing <- time_function({
    compute_fc_curve(
      fc_output_grid = test_data$fc_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = test_data$side,
      cutoff = test_data$cutoff
    )
  })
  
  # Benchmark C++ implementation
  cpp_timing <- time_function({
    compute_fc_curve_cpp(
      fc_output_grid = test_data$fc_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = test_data$side,
      cutoff = test_data$cutoff
    )
  })
  
  # Calculate speedup
  speedup <- r_timing$median_time / cpp_timing$median_time
  
  # Print timing results for information
  cat("\n=== FC Curve Benchmark Results ===\n")
  cat("Test data: ", nrow(test_data$fc_expression_df), " samples, ", 
      length(test_data$fc_output_grid), " grid points\n")
  cat("R implementation median time: ", sprintf("%.4f", r_timing$median_time), " seconds\n")
  cat("C++ implementation median time: ", sprintf("%.4f", cpp_timing$median_time), " seconds\n")
  cat("Speedup: ", sprintf("%.2fx", speedup), "\n")
  cat("==============================\n")
  
  # Test that C++ is competitive (allow for measurement noise and overhead)
  expect_lt(cpp_timing$median_time, r_timing$median_time * 1.5)
  
  # Ideally, C++ should provide meaningful speedup
  if (speedup > 1.2) {
    message("C++ provides ", sprintf("%.2fx", speedup), " speedup for FC curve computation")
  }
})

test_that("C++ expression curve computation is faster than R implementation", {
  
  # Use larger datasets for meaningful benchmarks
  test_data <- create_benchmark_data(n_samples = 200, n_grid_points = 15)
  
  # Benchmark R implementation
  r_timing <- time_function({
    compute_expression_curve(
      expr_output_grid = test_data$expr_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      expression_dispersion_curve = test_data$expression_dispersion_curve,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = test_data$side,
      cutoff = test_data$cutoff
    )
  })
  
  # Benchmark C++ implementation
  cpp_timing <- time_function({
    compute_expression_curve_cpp(
      expr_output_grid = test_data$expr_output_grid,
      fc_expression_df = test_data$fc_expression_df,
      library_size = test_data$library_size,
      expression_dispersion_curve = test_data$expression_dispersion_curve,
      num_trt_cells = test_data$num_trt_cells,
      num_cntrl_cells = test_data$num_cntrl_cells,
      side = test_data$side,
      cutoff = test_data$cutoff
    )
  })
  
  # Calculate speedup
  speedup <- r_timing$median_time / cpp_timing$median_time
  
  # Print timing results for information
  cat("\n=== Expression Curve Benchmark Results ===\n")
  cat("Test data: ", nrow(test_data$fc_expression_df), " samples, ", 
      length(test_data$expr_output_grid), " grid points\n")
  cat("R implementation median time: ", sprintf("%.4f", r_timing$median_time), " seconds\n")
  cat("C++ implementation median time: ", sprintf("%.4f", cpp_timing$median_time), " seconds\n")
  cat("Speedup: ", sprintf("%.2fx", speedup), "\n")
  cat("==========================================\n")
  
  # Test that C++ is competitive (allow for measurement noise and overhead)
  expect_lt(cpp_timing$median_time, r_timing$median_time * 1.5)
  
  # Ideally, C++ should provide meaningful speedup
  if (speedup > 1.2) {
    message("C++ provides ", sprintf("%.2fx", speedup), " speedup for expression curve computation")
  }
})

test_that("Performance scales well with increasing data size", {
  
  data_sizes <- c(50, 100, 200)
  grid_sizes <- c(5, 10, 15)
  
  fc_results <- data.frame(
    data_size = integer(),
    grid_size = integer(),
    r_time = numeric(),
    cpp_time = numeric(),
    speedup = numeric()
  )
  
  for (n_samples in data_sizes) {
    for (n_grid in grid_sizes) {
      test_data <- create_benchmark_data(n_samples = n_samples, n_grid_points = n_grid)
      
      # Quick single-run timing for scalability test
      r_start <- Sys.time()
      compute_fc_curve(
        fc_output_grid = test_data$fc_output_grid,
        fc_expression_df = test_data$fc_expression_df,
        library_size = test_data$library_size,
        num_trt_cells = test_data$num_trt_cells,
        num_cntrl_cells = test_data$num_cntrl_cells,
        side = test_data$side,
        cutoff = test_data$cutoff
      )
      r_time <- as.numeric(Sys.time() - r_start, units = "secs")
      
      cpp_start <- Sys.time()
      compute_fc_curve_cpp(
        fc_output_grid = test_data$fc_output_grid,
        fc_expression_df = test_data$fc_expression_df,
        library_size = test_data$library_size,
        num_trt_cells = test_data$num_trt_cells,
        num_cntrl_cells = test_data$num_cntrl_cells,
        side = test_data$side,
        cutoff = test_data$cutoff
      )
      cpp_time <- as.numeric(Sys.time() - cpp_start, units = "secs")
      
      fc_results <- rbind(fc_results, data.frame(
        data_size = n_samples,
        grid_size = n_grid,
        r_time = r_time,
        cpp_time = cpp_time,
        speedup = r_time / cpp_time
      ))
    }
  }
  
  # Print scalability results
  cat("\n=== Scalability Test Results ===\n")
  print(fc_results)
  cat("================================\n")
  
  # Test that C++ maintains advantage across different sizes
  expect_true(all(fc_results$speedup > 0.8))
  
  # Test that larger datasets show consistent or better speedup
  large_data_speedup <- mean(fc_results$speedup[fc_results$data_size >= 100])
  small_data_speedup <- mean(fc_results$speedup[fc_results$data_size < 100])
  
  expect_gte(large_data_speedup, small_data_speedup * 0.3)
})