library(testthat)

test_that("compute_monte_carlo_teststat_new_cpp validates inputs", {
  # Create valid test data
  test_df <- data.frame(
    relative_expression = c(0.001, 0.002),
    expression_size = c(1.2, 0.8),
    avg_fold_change = c(0.8, 1.5),
    avg_fold_change_sq = c(0.7, 2.5)
  )
  
  # Test with missing columns
  incomplete_df <- data.frame(
    relative_expression = c(0.001, 0.002),
    expression_size = c(1.2, 0.8)
    # Missing avg_fold_change and avg_fold_change_sq
  )
  
  expect_error(
    compute_monte_carlo_teststat_new_cpp(
      fc_expression_df = incomplete_df,
      library_size = 5000,
      num_trt_cells = 100,
      num_cntrl_cells = 200
    ),
    "Index out of bounds"
  )
  
  # Test with negative library size
  expect_error(
    compute_monte_carlo_teststat_new_cpp(
      fc_expression_df = test_df,
      library_size = -1000,
      num_trt_cells = 100,
      num_cntrl_cells = 200
    ),
    "Expression parameters must be positive"
  )
  
  # Test with zero cell counts
  expect_error(
    compute_monte_carlo_teststat_new_cpp(
      fc_expression_df = test_df,
      library_size = 5000,
      num_trt_cells = 0,
      num_cntrl_cells = 200
    ),
    "Cell counts must be positive"
  )
})

test_that("compute_monte_carlo_teststat_new_cpp produces correct output structure", {
  # Create test data
  n_samples <- 5
  test_df <- data.frame(
    relative_expression = c(0.001, 0.002, 0.0005, 0.003, 0.0015),
    expression_size = c(1.2, 0.8, 2.0, 1.5, 0.9),
    avg_fold_change = c(0.8, 1.5, 0.6, 1.2, 0.9),
    avg_fold_change_sq = c(0.7, 2.5, 0.4, 1.5, 0.85)
  )
  
  result <- compute_monte_carlo_teststat_new_cpp(
    fc_expression_df = test_df,
    library_size = 5000,
    num_trt_cells = 100,
    num_cntrl_cells = 200
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("means", "sds") %in% names(result)))
  
  # Check dimensions
  expect_equal(length(result$means), n_samples)
  expect_equal(length(result$sds), n_samples)
  
  # Check all values are finite
  expect_true(all(is.finite(result$means)))
  expect_true(all(is.finite(result$sds)))
  
  # Check all SDs are positive
  expect_true(all(result$sds > 0))
})

test_that("compute_monte_carlo_teststat_new_cpp matches individual function calls", {
  # Create test data
  test_df <- data.frame(
    relative_expression = c(0.001, 0.002, 0.0005),
    expression_size = c(1.2, 0.8, 2.0),
    avg_fold_change = c(0.8, 1.5, 0.6),
    avg_fold_change_sq = c(0.7, 2.5, 0.4)
  )
  
  library_size <- 5000
  num_trt_cells <- 100
  num_cntrl_cells <- 200
  
  # Test batch function
  batch_result <- compute_monte_carlo_teststat_new_cpp(
    fc_expression_df = test_df,
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells
  )
  
  # Test individual function calls
  individual_means <- numeric(nrow(test_df))
  individual_sds <- numeric(nrow(test_df))
  
  for (i in 1:nrow(test_df)) {
    individual_result <- compute_distribution_teststat_random_es_cpp(
      num_trt_cell = num_trt_cells,
      num_cntrl_cell = num_cntrl_cells,
      expression_mean = library_size * test_df$relative_expression[i],
      expression_size = test_df$expression_size[i],
      avg_fold_change = test_df$avg_fold_change[i],
      avg_fold_change_sq = test_df$avg_fold_change_sq[i]
    )
    individual_means[i] <- individual_result$mean
    individual_sds[i] <- individual_result$sd
  }
  
  # Compare results
  expect_equal(batch_result$means, individual_means, tolerance = 1e-10)
  expect_equal(batch_result$sds, individual_sds, tolerance = 1e-10)
})

test_that("compute_monte_carlo_teststat_new_cpp handles edge cases", {
  # Test with single sample
  single_df <- data.frame(
    relative_expression = 0.001,
    expression_size = 1.2,
    avg_fold_change = 0.8,
    avg_fold_change_sq = 0.7
  )
  
  result1 <- compute_monte_carlo_teststat_new_cpp(
    fc_expression_df = single_df,
    library_size = 5000,
    num_trt_cells = 100,
    num_cntrl_cells = 200
  )
  
  expect_equal(length(result1$means), 1)
  expect_equal(length(result1$sds), 1)
  expect_true(is.finite(result1$means[1]))
  expect_true(is.finite(result1$sds[1]))
  expect_true(result1$sds[1] > 0)
  
  # Test with very small library size
  result2 <- compute_monte_carlo_teststat_new_cpp(
    fc_expression_df = single_df,
    library_size = 1,
    num_trt_cells = 100,
    num_cntrl_cells = 200
  )
  
  expect_true(is.finite(result2$means[1]))
  expect_true(is.finite(result2$sds[1]))
  expect_true(result2$sds[1] > 0)
  
  # Test with fold change = 1 (no effect)
  no_effect_df <- data.frame(
    relative_expression = 0.001,
    expression_size = 1.2,
    avg_fold_change = 1.0,
    avg_fold_change_sq = 1.0
  )
  
  result3 <- compute_monte_carlo_teststat_new_cpp(
    fc_expression_df = no_effect_df,
    library_size = 5000,
    num_trt_cells = 100,
    num_cntrl_cells = 200
  )
  
  expect_equal(result3$means[1], 0, tolerance = 1e-10)  # Should be zero when no effect
  expect_true(result3$sds[1] > 0)
})

test_that("compute_monte_carlo_teststat_new_cpp handles different parameter ranges", {
  # Test with various parameter combinations
  test_cases <- list(
    list(
      df = data.frame(
        relative_expression = c(0.001, 0.01),
        expression_size = c(0.5, 5.0),
        avg_fold_change = c(0.5, 2.0),
        avg_fold_change_sq = c(0.3, 4.5)
      ),
      library_size = 1000,
      num_trt_cells = 50,
      num_cntrl_cells = 100
    ),
    list(
      df = data.frame(
        relative_expression = c(0.0001, 0.1),
        expression_size = c(0.1, 10.0),
        avg_fold_change = c(0.1, 5.0),
        avg_fold_change_sq = c(0.011, 25.5)  # Fixed: 0.011 > 0.1^2 = 0.01, 25.5 > 5.0^2 = 25.0
      ),
      library_size = 10000,
      num_trt_cells = 200,
      num_cntrl_cells = 500
    )
  )
  
  for (i in seq_along(test_cases)) {
    case <- test_cases[[i]]
    
    result <- compute_monte_carlo_teststat_new_cpp(
      fc_expression_df = case$df,
      library_size = case$library_size,
      num_trt_cells = case$num_trt_cells,
      num_cntrl_cells = case$num_cntrl_cells
    )
    
    # Check basic validity
    expect_true(all(is.finite(result$means)), info = paste("Case", i, "- finite means"))
    expect_true(all(is.finite(result$sds)), info = paste("Case", i, "- finite sds"))
    expect_true(all(result$sds > 0), info = paste("Case", i, "- positive sds"))
    expect_equal(length(result$means), nrow(case$df), info = paste("Case", i, "- length"))
  }
})

test_that("compute_monte_carlo_teststat_new_cpp second moment constraint", {
  # Test with valid second moment constraint
  valid_df <- data.frame(
    relative_expression = c(0.001, 0.002),
    expression_size = c(1.2, 0.8),
    avg_fold_change = c(0.8, 1.5),
    avg_fold_change_sq = c(0.8^2, 1.5^2)  # Equal to first moment squared
  )
  
  result1 <- compute_monte_carlo_teststat_new_cpp(
    fc_expression_df = valid_df,
    library_size = 5000,
    num_trt_cells = 100,
    num_cntrl_cells = 200
  )
  
  expect_true(all(is.finite(result1$means)))
  expect_true(all(is.finite(result1$sds)))
  expect_true(all(result1$sds > 0))
  
  # Test with invalid second moment constraint (should fail)
  invalid_df <- data.frame(
    relative_expression = c(0.001, 0.002),
    expression_size = c(1.2, 0.8),
    avg_fold_change = c(0.8, 1.5),
    avg_fold_change_sq = c(0.5, 1.0)  # Less than first moment squared
  )
  
  expect_error(
    compute_monte_carlo_teststat_new_cpp(
      fc_expression_df = invalid_df,
      library_size = 5000,
      num_trt_cells = 100,
      num_cntrl_cells = 200
    ),
    "avg_fold_change_sq must be >= avg_fold_change\\^2"
  )
})