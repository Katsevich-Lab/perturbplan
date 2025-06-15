# Test for compute_distribution_teststat_fixed_es_cpp function
# Compares C++ implementation against R score_test implementation

library(testthat)

test_that("compute_distribution_teststat_fixed_es_cpp gives similar results to score_test", {

  # Load Gasperini baseline expression data
  gas_data <- readRDS(system.file('extdata/baseline_expression', 'Gasperini_expression.rds', package = 'perturbplan'))
  baseline_df <- gas_data$baseline_expression

  # Transform relative expression to mean expression with library size 10000
  library_size <- 10000
  baseline_df$mean_expression <- baseline_df$relative_expression * library_size

  # Sample 10 representative gene pairs across expression range
  set.seed(123)
  sampled_indices <- sample(nrow(baseline_df), 10)
  test_genes <- baseline_df[sampled_indices, ]

  # Test parameters
  fold_change_values <- c(0.5, 1, 2)
  num_trt_cells_values <- 500
  num_cntrl_cells_values <- 1000
  n_simulations <- 5000
  tolerance <- 0.1 # 10% tolerance

  # Test across different parameter combinations
  for (i in 1:nrow(test_genes)) {
    gene <- test_genes[i, ]
    expression_mean <- gene$mean_expression
    expression_size <- gene$expression_size

    for (fold_change in fold_change_values) {
      for (num_trt_cells in num_trt_cells_values) {
        for (num_cntrl_cells in num_cntrl_cells_values) {

          # Test C++ implementation (returns list with mean and sd)
          cpp_result <- compute_distribution_teststat_fixed_es_cpp(
            fold_change = fold_change,
            expression_mean = expression_mean,
            expression_size = expression_size,
            num_trt_cells = num_trt_cells,
            num_cntrl_cells = num_cntrl_cells,
            num_cells = num_trt_cells
          )

          # Test R implementation using score_test
          r_results <- numeric(n_simulations)

          for (sim in 1:n_simulations) {
            # Generate treatment group data (with fold change effect)
            trt_mean <- expression_mean * fold_change
            trt_counts <- rnbinom(num_trt_cells, size = expression_size, mu = trt_mean)

            # Generate control group data
            ctrl_counts <- rnbinom(num_cntrl_cells, size = expression_size, mu = expression_mean)

            # Combine data
            Y <- c(trt_counts, ctrl_counts)
            X <- c(rep(1, num_trt_cells), rep(0, num_cntrl_cells))

            # Compute test statistic using R implementation
            r_results[sim] <- score_test(X, Y, expression_size)
          }

          # Extract C++ results
          cpp_mean <- cpp_result$mean
          cpp_sd <- cpp_result$sd

          # Compute R results
          r_mean <- mean(r_results, na.rm = TRUE)
          r_sd <- sd(r_results, na.rm = TRUE)

          # Compare means
          mean_diff <- abs(cpp_mean - r_mean)
          expect_true(mean_diff < tolerance,
                     info = paste("Mean difference too large for gene", i,
                                "fold_change", fold_change, "cells", num_trt_cells, num_cntrl_cells,
                                "- C++:", round(cpp_mean, 4), "R:", round(r_mean, 4),
                                "diff:", round(mean_diff, 4)))

          # Compare standard deviations
          sd_diff <- abs(cpp_sd - r_sd)
          expect_true(sd_diff < tolerance,
                     info = paste("SD difference too large for gene", i,
                                "fold_change", fold_change, "cells", num_trt_cells, num_cntrl_cells,
                                "- C++:", round(cpp_sd, 4), "R:", round(r_sd, 4),
                                "diff:", round(sd_diff, 4)))
        }
      }
    }
  }
})
