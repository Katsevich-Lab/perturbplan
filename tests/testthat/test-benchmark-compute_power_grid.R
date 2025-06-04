library(testthat)
library(dplyr)

test_that("benchmark compute_power_grid_efficient vs compute_power_grid_separated", {
  
  # Skip benchmark on CI or if not explicitly requested
  skip_if_not(identical(Sys.getenv("RUN_BENCHMARKS"), "true"), "Benchmarks not requested")

  # Set up reasonable test case for demonstration
  cells_reads_df <- data.frame(
    num_total_cells = c(15000, 25000),
    reads_per_cell = c(1000, 1500)
  )

  # Use modest parameters for reasonable test execution time
  benchmark_params <- list(
    num_targets = 20,
    gRNAs_per_target = 3,
    non_targeting_gRNAs = 5,
    num_pairs = 200,
    fdr_target = 0.05,
    fc_mean = 0.8,
    fc_sd = 0.1,
    prop_non_null = 0.15,
    B = 100,  # Reasonable for testing
    fc_curve_points = 8,
    expr_curve_points = 8
  )

  cat("\n=== Performance Benchmark ===\n")
  cat("Test parameters:\n")
  cat("- Cells/reads combinations:", nrow(cells_reads_df), "\n")
  cat("- Monte Carlo samples (B):", benchmark_params$B, "\n")
  cat("- Targets:", benchmark_params$num_targets, "\n")
  cat("- gRNAs per target:", benchmark_params$gRNAs_per_target, "\n")
  cat("- Curve points (FC x Expr):", benchmark_params$fc_curve_points, "x", benchmark_params$expr_curve_points, "\n")

  # Calculate expected test statistic computations
  total_computations_per_condition <- benchmark_params$B * (1 + benchmark_params$fc_curve_points + benchmark_params$expr_curve_points)
  total_computations <- nrow(cells_reads_df) * total_computations_per_condition
  cat("- Expected test statistic calls:", total_computations, "\n\n")

  # Benchmark original R implementation
  cat("Benchmarking R implementation (compute_power_grid_separated)...\n")
  time_r_start <- Sys.time()

  set.seed(123)
  result_r <- compute_power_grid_separated(
    cells_reads_df = cells_reads_df,
    num_targets = benchmark_params$num_targets,
    gRNAs_per_target = benchmark_params$gRNAs_per_target,
    non_targeting_gRNAs = benchmark_params$non_targeting_gRNAs,
    num_pairs = benchmark_params$num_pairs,
    fdr_target = benchmark_params$fdr_target,
    fc_mean = benchmark_params$fc_mean,
    fc_sd = benchmark_params$fc_sd,
    prop_non_null = benchmark_params$prop_non_null,
    B = benchmark_params$B,
    fc_curve_points = benchmark_params$fc_curve_points,
    expr_curve_points = benchmark_params$expr_curve_points
  )

  time_r_end <- Sys.time()
  time_r <- as.numeric(difftime(time_r_end, time_r_start, units = "secs"))

  # Benchmark C++ implementation
  cat("Benchmarking C++ implementation (compute_power_grid_efficient)...\n")
  time_cpp_start <- Sys.time()

  set.seed(123)
  result_cpp <- compute_power_grid_efficient(
    cells_reads_df = cells_reads_df,
    num_targets = benchmark_params$num_targets,
    gRNAs_per_target = benchmark_params$gRNAs_per_target,
    non_targeting_gRNAs = benchmark_params$non_targeting_gRNAs,
    num_pairs = benchmark_params$num_pairs,
    fdr_target = benchmark_params$fdr_target,
    fc_mean = benchmark_params$fc_mean,
    fc_sd = benchmark_params$fc_sd,
    prop_non_null = benchmark_params$prop_non_null,
    B = benchmark_params$B,
    fc_curve_points = benchmark_params$fc_curve_points,
    expr_curve_points = benchmark_params$expr_curve_points
  )

  time_cpp_end <- Sys.time()
  time_cpp <- as.numeric(difftime(time_cpp_end, time_cpp_start, units = "secs"))

  # Calculate performance metrics
  speedup <- time_r / time_cpp
  percent_improvement <- ((time_r - time_cpp) / time_r) * 100
  computations_per_sec_r <- total_computations / time_r
  computations_per_sec_cpp <- total_computations / time_cpp

  # Display results
  cat("\n=== Benchmark Results ===\n")
  cat(sprintf("R implementation:   %.3f seconds (%.0f computations/sec)\n", time_r, computations_per_sec_r))
  cat(sprintf("C++ implementation: %.3f seconds (%.0f computations/sec)\n", time_cpp, computations_per_sec_cpp))
  cat(sprintf("Speedup:            %.2fx faster\n", speedup))
  cat(sprintf("Performance gain:   %.1f%% improvement\n", percent_improvement))
  cat(sprintf("Time saved:         %.3f seconds\n", time_r - time_cpp))

  # Verify results are still identical
  cat("\n=== Correctness Verification ===\n")
  power_diff <- abs(result_r$overall_power - result_cpp$overall_power)
  max_power_diff <- max(abs(unlist(result_r$power_by_fc[[1]]$power) - unlist(result_cpp$power_by_fc[[1]]$power)))

  cat(sprintf("Overall power difference:     %.2e\n", max(power_diff)))
  cat(sprintf("Max power curve difference:   %.2e\n", max_power_diff))

  # Test assertions
  expect_true(time_cpp < time_r, info = "C++ implementation should be faster")
  expect_true(speedup > 1.0, info = "Should show measurable speedup")
  expect_equal(result_r$overall_power, result_cpp$overall_power, tolerance = 1e-10,
               info = "Results should remain identical")
})

test_that("micro-benchmark individual test statistic functions", {

  # Skip benchmark on CI or if not explicitly requested
  skip_if_not(identical(Sys.getenv("RUN_BENCHMARKS"), "true"), "Benchmarks not requested")

  # Test parameters
  num_trt_cells <- 150
  num_cntrl_cells <- 200
  expression_mean <- 5.0
  expression_size <- 2.5
  fold_change_mean <- 0.8

  # Number of iterations for micro-benchmark
  n_iterations <- 10000

  cat("\n=== Micro-benchmark: Test Statistic Functions ===\n")
  cat("Iterations:", n_iterations, "\n")
  cat("Parameters: num_trt_cells =", num_trt_cells, ", num_cntrl_cells =", num_cntrl_cells, "\n")
  cat("            expression_mean =", expression_mean, ", fold_change =", fold_change_mean, "\n\n")

  # Benchmark R function
  cat("Benchmarking R function (compute_test_stat_clean)...\n")
  time_r_start <- Sys.time()

  for (i in 1:n_iterations) {
    result_r <- compute_test_stat_clean(
      num_trt_cells = num_trt_cells,
      num_cntrl_cells = num_cntrl_cells,
      expression_mean = expression_mean,
      expression_size = expression_size,
      fold_change_mean = fold_change_mean
    )
  }

  time_r_end <- Sys.time()
  time_r_micro <- as.numeric(difftime(time_r_end, time_r_start, units = "secs"))

  # Benchmark C++ function
  cat("Benchmarking C++ function (compute_distribution_teststat_fixed_es_cpp)...\n")
  time_cpp_start <- Sys.time()

  for (i in 1:n_iterations) {
    result_cpp <- compute_distribution_teststat_fixed_es_cpp(
      fold_change = fold_change_mean,
      expression_mean = expression_mean,
      expression_size = expression_size,
      num_trt_cells = num_trt_cells,
      num_cntrl_cells = num_cntrl_cells,
      num_cells = num_trt_cells
    )
  }

  time_cpp_end <- Sys.time()
  time_cpp_micro <- as.numeric(difftime(time_cpp_end, time_cpp_start, units = "secs"))

  # Calculate micro-benchmark metrics
  speedup_micro <- time_r_micro / time_cpp_micro
  time_per_call_r <- (time_r_micro / n_iterations) * 1000  # milliseconds
  time_per_call_cpp <- (time_cpp_micro / n_iterations) * 1000  # milliseconds

  cat("\n=== Micro-benchmark Results ===\n")
  cat(sprintf("R function:   %.3f seconds total (%.4f ms/call)\n", time_r_micro, time_per_call_r))
  cat(sprintf("C++ function: %.3f seconds total (%.4f ms/call)\n", time_cpp_micro, time_per_call_cpp))
  cat(sprintf("Speedup:      %.2fx faster per function call\n", speedup_micro))

  # Verify results are identical
  expect_equal(result_r[["mean"]], result_cpp[["mean"]], tolerance = 1e-12)
  expect_equal(result_r[["sd"]], result_cpp[["sd"]], tolerance = 1e-12)

  expect_true(speedup_micro > 1.0, info = "C++ should be faster than R for individual calls")

  cat("=== Micro-benchmark Complete ===\n\n")
})
