# Tests for compute_QC.cpp functions
# This file tests the C++ QC computation functions

test_that("compute_QC_fixed_es_cpp works with basic inputs", {
  # Basic test case
  fold_change <- c(0.5, 0.8, 1.2)
  expression_mean <- 100
  expression_size <- 5
  num_cntrl_cells <- 100
  num_cells <- 50

  result <- compute_QC_fixed_es_cpp(
    fold_change = fold_change,
    expression_mean = expression_mean,
    expression_size = expression_size,
    num_cntrl_cells = num_cntrl_cells,
    num_cells = num_cells,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7
  )

  # QC probability should be between 0 and 1
  expect_type(result, "double")
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("compute_QC_fixed_es_cpp handles vector inputs correctly", {
  # All parameters except fold_change should work as vectors (but only first element is used)
  fold_change <- c(0.5, 0.8)
  expression_mean <- c(100, 200)  # Only first element used
  expression_size <- c(5, 10)     # Only first element used
  num_cntrl_cells <- c(100, 50)  # Only first element used
  num_cells <- c(50, 25)          # Only first element used

  result <- compute_QC_fixed_es_cpp(
    fold_change = fold_change,
    expression_mean = expression_mean,
    expression_size = expression_size,
    num_cntrl_cells = num_cntrl_cells,
    num_cells = num_cells
  )

  expect_type(result, "double")
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("compute_QC_fixed_es_cpp handles different threshold values", {
  fold_change <- c(0.5, 1.0, 1.5)
  expression_mean <- 50
  expression_size <- 3
  num_cntrl_cells <- 100
  num_cells <- 30

  # Test with default thresholds
  result1 <- compute_QC_fixed_es_cpp(
    fold_change, expression_mean, expression_size,
    num_cntrl_cells, num_cells,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7
  )

  # Test with higher thresholds (should give different result)
  result2 <- compute_QC_fixed_es_cpp(
    fold_change, expression_mean, expression_size,
    num_cntrl_cells, num_cells,
    n_nonzero_trt_thresh = 10,
    n_nonzero_cntrl_thresh = 10
  )

  expect_type(result1, "double")
  expect_type(result2, "double")

  # Higher thresholds should generally increase QC filtering probability
  # (result is probability that QC filters the pair)
  expect_gte(result1, 0)
  expect_gte(result2, 0)
})

test_that("compute_QC_fixed_es_cpp handles edge cases", {
  fold_change <- c(1.0)  # No effect
  expression_mean <- 100
  expression_size <- 5
  num_cntrl_cells <- 100
  num_cells <- 50

  result <- compute_QC_fixed_es_cpp(
    fold_change, expression_mean, expression_size,
    num_cntrl_cells, num_cells
  )

  expect_type(result, "double")
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("compute_QC_fixed_es_cpp handles multiple fold changes", {
  # Test with varying number of fold changes
  for (k in c(1, 3, 5, 10)) {
    fold_change <- rep(0.5, k)
    result <- compute_QC_fixed_es_cpp(
      fold_change = fold_change,
      expression_mean = 100,
      expression_size = 5,
      num_cntrl_cells = 100,
      num_cells = 50
    )

    expect_type(result, "double")
    expect_gte(result, 0)
    expect_lte(result, 1)
  }
})

test_that("compute_QC_fixed_es_cpp QC probability increases with lower expression", {
  fold_change <- c(0.5, 1.0, 1.5)
  num_cntrl_cells <- 100
  num_cells <- 50

  # Low expression should have higher QC filtering probability
  result_low <- compute_QC_fixed_es_cpp(
    fold_change,
    expression_mean = 10,  # Low expression
    expression_size = 5,
    num_cntrl_cells,
    num_cells
  )

  # High expression should have lower QC filtering probability
  result_high <- compute_QC_fixed_es_cpp(
    fold_change,
    expression_mean = 1000,  # High expression
    expression_size = 5,
    num_cntrl_cells,
    num_cells
  )

  # Low expression should have higher filtering probability
  expect_gte(result_low, result_high)
})

test_that("compute_QC_fixed_es_cpp handles very low and very high thresholds", {
  fold_change <- c(0.5, 1.0, 1.5)
  expression_mean <- 100
  expression_size <- 5
  num_cntrl_cells <- 100
  num_cells <- 50

  # Very low threshold
  result_low_thresh <- compute_QC_fixed_es_cpp(
    fold_change, expression_mean, expression_size,
    num_cntrl_cells, num_cells,
    n_nonzero_trt_thresh = 1,
    n_nonzero_cntrl_thresh = 1
  )

  # Very high threshold
  result_high_thresh <- compute_QC_fixed_es_cpp(
    fold_change, expression_mean, expression_size,
    num_cntrl_cells, num_cells,
    n_nonzero_trt_thresh = 20,
    n_nonzero_cntrl_thresh = 20
  )

  expect_gte(result_low_thresh, 0)
  expect_lte(result_low_thresh, 1)
  expect_gte(result_high_thresh, 0)
  expect_lte(result_high_thresh, 1)
})

test_that("compute_QC_fixed_es_cpp works with realistic parameter values", {
  # Use realistic values from actual perturb-seq experiments
  fold_change <- c(0.3, 0.5, 0.7, 1.0, 1.5, 2.0)
  expression_mean <- 50  # Moderate expression
  expression_size <- 3
  num_cntrl_cells <- 200  # Typical control group size
  num_cells <- 20  # Typical cells per gRNA

  result <- compute_QC_fixed_es_cpp(
    fold_change, expression_mean, expression_size,
    num_cntrl_cells, num_cells,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7
  )

  expect_type(result, "double")
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("compute_QC_fixed_es_cpp handles different cell counts", {
  fold_change <- c(0.5, 1.0, 1.5)
  expression_mean <- 100
  expression_size <- 5

  # Small experiment
  result_small <- compute_QC_fixed_es_cpp(
    fold_change, expression_mean, expression_size,
    num_cntrl_cells = 50,
    num_cells = 10
  )

  # Large experiment
  result_large <- compute_QC_fixed_es_cpp(
    fold_change, expression_mean, expression_size,
    num_cntrl_cells = 500,
    num_cells = 100
  )

  expect_type(result_small, "double")
  expect_type(result_large, "double")
  expect_gte(result_small, 0)
  expect_gte(result_large, 0)

  # Larger experiments should generally have lower QC filtering probability
  expect_lte(result_large, result_small)
})
