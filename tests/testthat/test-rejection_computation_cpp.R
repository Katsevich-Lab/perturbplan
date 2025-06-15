# Test for rejection_computation_cpp function
# Validates vectorized rejection probability computation against R implementations

library(testthat)

test_that("rejection_computation_cpp gives correct results compared to R implementation", {
  
  # Load Gasperini baseline expression data
  gas_data <- readRDS(system.file('extdata/baseline_expression', 'Gasperini_expression.rds', package = 'perturbplan'))
  baseline_df <- gas_data$baseline_expression
  
  # Sample representative genes across expression range
  set.seed(42)
  n_test_genes <- 20
  sampled_indices <- sample(nrow(baseline_df), n_test_genes)
  test_genes <- baseline_df[sampled_indices, ]
  
  # Simulate test statistics using realistic parameters
  library_size <- 10000
  test_genes$mean_expression <- test_genes$relative_expression * library_size
  
  # Generate test statistic means and standard deviations
  # Simulate various experimental scenarios with different effect sizes
  fold_changes <- c(0.5, 0.8, 1.0, 1.2, 2.0)  # Range of biological effects
  num_cells_pairs <- list(c(300, 500), c(500, 1000), c(1000, 1500))  # (treatment, control)
  
  test_scenarios <- expand.grid(
    gene_idx = 1:min(5, nrow(test_genes)),  # Test subset of genes
    fold_change = fold_changes,
    cell_pair_idx = 1:length(num_cells_pairs),
    stringsAsFactors = FALSE
  )
  
  # Generate test statistics for each scenario
  test_stats <- data.frame(
    mean_test_stat = numeric(nrow(test_scenarios)),
    sd_test_stat = numeric(nrow(test_scenarios))
  )
  
  for (i in 1:nrow(test_scenarios)) {
    gene <- test_genes[test_scenarios$gene_idx[i], ]
    fold_change <- test_scenarios$fold_change[i]
    cells <- num_cells_pairs[[test_scenarios$cell_pair_idx[i]]]
    num_trt_cells <- cells[1]
    num_ctrl_cells <- cells[2]
    
    # Use the C++ function we already tested to get test statistic distribution
    cpp_result <- compute_distribution_teststat_fixed_es_cpp(
      fold_change = fold_change,
      expression_mean = gene$mean_expression,
      expression_size = gene$expression_size,
      num_trt_cells = num_trt_cells,
      num_cntrl_cells = num_ctrl_cells,
      num_cells = num_trt_cells
    )
    
    test_stats$mean_test_stat[i] <- cpp_result$mean
    test_stats$sd_test_stat[i] <- cpp_result$sd
  }
  
  # Test parameters
  cutoff_values <- c(0.01, 0.05, 0.1, 0.2)
  test_sides <- c("left", "right", "both")
  tolerance <- 1e-10  # Very tight tolerance for numerical precision
  
  # Test each combination of cutoff and side
  for (cutoff in cutoff_values) {
    for (side in test_sides) {
      
      # Test C++ implementation
      cpp_probs <- rejection_computation_cpp(
        mean_list = test_stats$mean_test_stat,
        sd_list = test_stats$sd_test_stat,
        side = side,
        cutoff = cutoff
      )
      
      # Compute R equivalent manually
      r_probs <- numeric(length(test_stats$mean_test_stat))
      
      if (side == "left") {
        # P(Z < threshold) where Z ~ N(mean, sd)
        threshold <- qnorm(cutoff)
        for (j in 1:length(r_probs)) {
          r_probs[j] <- pnorm(threshold, mean = test_stats$mean_test_stat[j], 
                             sd = test_stats$sd_test_stat[j])
        }
        
      } else if (side == "right") {
        # P(Z > threshold) where Z ~ N(mean, sd)
        threshold <- qnorm(1 - cutoff)
        for (j in 1:length(r_probs)) {
          r_probs[j] <- pnorm(threshold, mean = test_stats$mean_test_stat[j], 
                             sd = test_stats$sd_test_stat[j], lower.tail = FALSE)
        }
        
      } else if (side == "both") {
        # P(|Z| > threshold) = P(Z > threshold) + P(Z < -threshold)
        threshold_upper <- qnorm(1 - cutoff/2)
        threshold_lower <- qnorm(cutoff/2)
        for (j in 1:length(r_probs)) {
          r_probs[j] <- pnorm(threshold_upper, mean = test_stats$mean_test_stat[j], 
                             sd = test_stats$sd_test_stat[j], lower.tail = FALSE) +
                       pnorm(threshold_lower, mean = test_stats$mean_test_stat[j], 
                             sd = test_stats$sd_test_stat[j], lower.tail = TRUE)
        }
      }
      
      # Compare C++ and R results
      max_diff <- max(abs(cpp_probs - r_probs))
      expect_true(max_diff < tolerance,
                 info = paste("Maximum difference too large for side =", side,
                            "cutoff =", cutoff, "- max diff:", formatC(max_diff, format = "e")))
      
      # Test specific properties
      expect_true(all(cpp_probs >= 0 & cpp_probs <= 1),
                 info = paste("Probabilities not in [0,1] for side =", side, "cutoff =", cutoff))
      
      expect_equal(length(cpp_probs), length(test_stats$mean_test_stat),
                  info = paste("Output length mismatch for side =", side, "cutoff =", cutoff))
    }
  }
  
  # Test edge cases
  
  # Test with extreme values
  extreme_means <- c(-10, -1, 0, 1, 10)
  extreme_sds <- c(0.1, 1, 5)
  
  for (mean_val in extreme_means) {
    for (sd_val in extreme_sds) {
      cpp_result <- rejection_computation_cpp(
        mean_list = mean_val,
        sd_list = sd_val,
        side = "both",
        cutoff = 0.05
      )
      
      r_result <- pnorm(qnorm(0.975), mean = mean_val, sd = sd_val, lower.tail = FALSE) +
                 pnorm(qnorm(0.025), mean = mean_val, sd = sd_val, lower.tail = TRUE)
      
      expect_equal(cpp_result[1], r_result, tolerance = tolerance,
                  info = paste("Extreme value test failed for mean =", mean_val, "sd =", sd_val))
    }
  }
  
  # Test error handling
  expect_error(rejection_computation_cpp(c(1, 2), c(1), "both", 0.05),
              "mean_list and sd_list must have identical length")
  
  expect_error(rejection_computation_cpp(c(1), c(1), "invalid_side", 0.05),
              "side must be 'left', 'right', 'both', or 'two.sided'")
})

test_that("rejection_computation_cpp handles vectorized inputs correctly", {
  
  # Test vectorization with varying parameters
  n_tests <- 100
  set.seed(123)
  
  test_means <- rnorm(n_tests, mean = 0, sd = 2)
  test_sds <- runif(n_tests, min = 0.5, max = 3)
  
  # Test all sides
  for (side in c("left", "right", "both")) {
    cpp_results <- rejection_computation_cpp(
      mean_list = test_means,
      sd_list = test_sds,
      side = side,
      cutoff = 0.05
    )
    
    expect_equal(length(cpp_results), n_tests)
    expect_true(all(cpp_results >= 0 & cpp_results <= 1))
    expect_true(all(is.finite(cpp_results)))
  }
})

test_that("rejection_computation_cpp matches two.sided and both aliases", {
  
  # Test that "both" and "two.sided" give identical results
  test_means <- c(-2, -1, 0, 1, 2)
  test_sds <- c(0.5, 1, 1.5, 2, 2.5)
  cutoff <- 0.05
  
  results_both <- rejection_computation_cpp(test_means, test_sds, "both", cutoff)
  results_two_sided <- rejection_computation_cpp(test_means, test_sds, "two.sided", cutoff)
  
  expect_equal(results_both, results_two_sided, tolerance = 1e-15)
})