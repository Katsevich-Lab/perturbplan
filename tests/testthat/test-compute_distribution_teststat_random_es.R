library(testthat)

test_that("compute_distribution_teststat_random_es_cpp validates inputs", {
  # Test negative cell counts
  expect_error(
    compute_distribution_teststat_random_es_cpp(
      num_trt_cell = -1,
      num_cntrl_cell = 100,
      expression_mean = 0.001,
      expression_size = 1.2,
      avg_fold_change = 0.8,
      avg_fold_change_sq = 0.7
    ),
    "Cell counts must be positive"
  )
  
  expect_error(
    compute_distribution_teststat_random_es_cpp(
      num_trt_cell = 50,
      num_cntrl_cell = 0,
      expression_mean = 0.001,
      expression_size = 1.2,
      avg_fold_change = 0.8,
      avg_fold_change_sq = 0.7
    ),
    "Cell counts must be positive"
  )
  
  # Test negative expression parameters
  expect_error(
    compute_distribution_teststat_random_es_cpp(
      num_trt_cell = 50,
      num_cntrl_cell = 100,
      expression_mean = -0.001,
      expression_size = 1.2,
      avg_fold_change = 0.8,
      avg_fold_change_sq = 0.7
    ),
    "Expression parameters must be positive"
  )
  
  expect_error(
    compute_distribution_teststat_random_es_cpp(
      num_trt_cell = 50,
      num_cntrl_cell = 100,
      expression_mean = 0.001,
      expression_size = -1.2,
      avg_fold_change = 0.8,
      avg_fold_change_sq = 0.7
    ),
    "Expression parameters must be positive"
  )
  
  # Test invalid fold change parameters
  expect_error(
    compute_distribution_teststat_random_es_cpp(
      num_trt_cell = 50,
      num_cntrl_cell = 100,
      expression_mean = 0.001,
      expression_size = 1.2,
      avg_fold_change = -0.8,
      avg_fold_change_sq = 0.7
    ),
    "Fold change parameters must be positive"
  )
  
  # Test second moment constraint violation
  expect_error(
    compute_distribution_teststat_random_es_cpp(
      num_trt_cell = 50,
      num_cntrl_cell = 100,
      expression_mean = 0.001,
      expression_size = 1.2,
      avg_fold_change = 0.8,
      avg_fold_change_sq = 0.5  # This violates avg_fold_change_sq >= avg_fold_change^2
    ),
    "avg_fold_change_sq must be >= avg_fold_change\\^2"
  )
})

test_that("compute_distribution_teststat_random_es_cpp matches R implementation", {
  # Test parameters
  test_cases <- list(
    list(
      num_trt_cell = 50,
      num_cntrl_cell = 100,
      expression_mean = 0.001,
      expression_size = 1.2,
      avg_fold_change = 0.8,
      avg_fold_change_sq = 0.7
    ),
    list(
      num_trt_cell = 200,
      num_cntrl_cell = 200,
      expression_mean = 0.01,
      expression_size = 2.0,
      avg_fold_change = 1.5,
      avg_fold_change_sq = 2.5
    ),
    list(
      num_trt_cell = 20,
      num_cntrl_cell = 80,
      expression_mean = 0.0001,
      expression_size = 0.5,
      avg_fold_change = 2.0,
      avg_fold_change_sq = 4.2
    )
  )
  
  for (i in seq_along(test_cases)) {
    params <- test_cases[[i]]
    
    # Test C++ function
    cpp_result <- compute_distribution_teststat_random_es_cpp(
      num_trt_cell = params$num_trt_cell,
      num_cntrl_cell = params$num_cntrl_cell,
      expression_mean = params$expression_mean,
      expression_size = params$expression_size,
      avg_fold_change = params$avg_fold_change,
      avg_fold_change_sq = params$avg_fold_change_sq
    )
    
    # Test R function
    r_result <- compute_distribution_teststat_random_es(
      num_trt_cell = params$num_trt_cell,
      num_cntrl_cell = params$num_cntrl_cell,
      expression_mean = params$expression_mean,
      expression_size = params$expression_size,
      avg_fold_change = params$avg_fold_change,
      avg_fold_change_sq = params$avg_fold_change_sq
    )
    
    # Compare results with tolerance for floating point precision
    expect_equal(cpp_result$mean, unname(r_result[[1]]["mean"]), tolerance = 1e-10,
                 info = paste("Test case", i, "- mean mismatch"))
    expect_equal(cpp_result$sd, unname(r_result[[1]]["sd"]), tolerance = 1e-10,
                 info = paste("Test case", i, "- sd mismatch"))
  }
})

test_that("compute_distribution_teststat_random_es_cpp returns valid results", {
  # Test with typical parameters
  result <- compute_distribution_teststat_random_es_cpp(
    num_trt_cell = 50,
    num_cntrl_cell = 100,
    expression_mean = 0.001,
    expression_size = 1.2,
    avg_fold_change = 0.8,
    avg_fold_change_sq = 0.7
  )
  
  # Check result structure
  expect_true(is.list(result))
  expect_true(all(c("mean", "sd") %in% names(result)))
  
  # Check result values are finite and valid
  expect_true(is.finite(result$mean))
  expect_true(is.finite(result$sd))
  expect_true(result$sd > 0)
})

test_that("compute_distribution_teststat_random_es_cpp handles edge cases", {
  # Test with very small expression mean
  result1 <- compute_distribution_teststat_random_es_cpp(
    num_trt_cell = 50,
    num_cntrl_cell = 100,
    expression_mean = 1e-10,
    expression_size = 1.2,
    avg_fold_change = 0.8,
    avg_fold_change_sq = 0.7
  )
  expect_true(is.finite(result1$mean))
  expect_true(is.finite(result1$sd))
  expect_true(result1$sd > 0)
  
  # Test with very large expression size (low dispersion)
  result2 <- compute_distribution_teststat_random_es_cpp(
    num_trt_cell = 50,
    num_cntrl_cell = 100,
    expression_mean = 0.001,
    expression_size = 1000,
    avg_fold_change = 0.8,
    avg_fold_change_sq = 0.7
  )
  expect_true(is.finite(result2$mean))
  expect_true(is.finite(result2$sd))
  expect_true(result2$sd > 0)
  
  # Test with fold change = 1 (no effect)
  result3 <- compute_distribution_teststat_random_es_cpp(
    num_trt_cell = 50,
    num_cntrl_cell = 100,
    expression_mean = 0.001,
    expression_size = 1.2,
    avg_fold_change = 1.0,
    avg_fold_change_sq = 1.0
  )
  expect_true(is.finite(result3$mean))
  expect_true(is.finite(result3$sd))
  expect_true(result3$sd > 0)
  expect_equal(result3$mean, 0, tolerance = 1e-10)  # Should be zero when no effect
})

test_that("compute_distribution_teststat_random_es_cpp is consistent with second moment", {
  # Test that function respects the constraint that avg_fold_change_sq >= avg_fold_change^2
  
  # Case 1: avg_fold_change_sq = avg_fold_change^2 (no variance)
  result1 <- compute_distribution_teststat_random_es_cpp(
    num_trt_cell = 50,
    num_cntrl_cell = 100,
    expression_mean = 0.001,
    expression_size = 1.2,
    avg_fold_change = 0.8,
    avg_fold_change_sq = 0.8^2
  )
  
  # Case 2: avg_fold_change_sq > avg_fold_change^2 (some variance)
  result2 <- compute_distribution_teststat_random_es_cpp(
    num_trt_cell = 50,
    num_cntrl_cell = 100,
    expression_mean = 0.001,
    expression_size = 1.2,
    avg_fold_change = 0.8,
    avg_fold_change_sq = 0.7
  )
  
  # Both should be valid
  expect_true(is.finite(result1$mean))
  expect_true(is.finite(result1$sd))
  expect_true(result1$sd > 0)
  
  expect_true(is.finite(result2$mean))
  expect_true(is.finite(result2$sd))
  expect_true(result2$sd > 0)
  
  # The mean should be the same (doesn't depend on variance)
  expect_equal(result1$mean, result2$mean, tolerance = 1e-10)
})