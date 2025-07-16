# Tests for obtain_expression_dispersion_curve function
library(testthat)

# Helper function to create test data
create_test_baseline_expression <- function(n_genes = 100, seed = 123) {
  set.seed(seed)
  data.frame(
    response_id = paste0("ENSG", sprintf("%011d", 1:n_genes)),
    relative_expression = sort(runif(n_genes, 1e-6, 1e-3)),
    expression_size = runif(n_genes, 0.1, 2.0),
    stringsAsFactors = FALSE
  )
}

test_that("obtain_expression_dispersion_curve works with valid input", {
  # Create test data
  baseline_expr <- create_test_baseline_expression(50)
  
  # Test function execution
  expect_silent({
    dispersion_curve <- obtain_expression_dispersion_curve(baseline_expr)
  })
  
  # Check that the function returns a function
  dispersion_curve <- obtain_expression_dispersion_curve(baseline_expr)
  expect_true(is.function(dispersion_curve))
  
  # Test that the returned function works
  test_expr_values <- c(1e-5, 5e-5, 1e-4, 5e-4)
  expect_silent({
    predicted_dispersion <- dispersion_curve(test_expr_values)
  })
  
  # Check output is numeric and has correct length
  predicted_dispersion <- dispersion_curve(test_expr_values)
  expect_true(is.numeric(predicted_dispersion))
  expect_equal(length(predicted_dispersion), length(test_expr_values))
  expect_true(all(predicted_dispersion > 0))
})

test_that("obtain_expression_dispersion_curve handles monotonicity", {
  # Create test data with known relationship
  baseline_expr <- data.frame(
    response_id = paste0("gene", 1:10),
    relative_expression = seq(1e-5, 1e-3, length.out = 10),
    expression_size = seq(0.5, 2.0, length.out = 10),
    stringsAsFactors = FALSE
  )
  
  dispersion_curve <- obtain_expression_dispersion_curve(baseline_expr)
  
  # Test monotonicity: higher expression should not lead to lower dispersion
  test_expr <- sort(runif(20, 1e-5, 1e-3))
  predicted_disp <- dispersion_curve(test_expr)
  
  # With isotonic regression, we expect non-decreasing dispersion
  expect_true(all(diff(predicted_disp) >= -1e-10))  # Allow for numerical precision
})

test_that("obtain_expression_dispersion_curve validates input data frame", {
  # Test with non-data.frame input
  expect_error(
    obtain_expression_dispersion_curve(c(1, 2, 3)),
    "baseline_expression must be a data frame"
  )
  
  expect_error(
    obtain_expression_dispersion_curve(matrix(1:6, nrow = 2)),
    "baseline_expression must be a data frame"
  )
})

test_that("obtain_expression_dispersion_curve validates required columns", {
  # Test missing relative_expression column
  baseline_expr_missing_rel <- data.frame(
    response_id = c("gene1", "gene2"),
    expression_size = c(0.5, 1.0)
  )
  expect_error(
    obtain_expression_dispersion_curve(baseline_expr_missing_rel),
    "Missing required columns: relative_expression"
  )
  
  # Test missing expression_size column
  baseline_expr_missing_size <- data.frame(
    response_id = c("gene1", "gene2"),
    relative_expression = c(1e-5, 2e-5)
  )
  expect_error(
    obtain_expression_dispersion_curve(baseline_expr_missing_size),
    "Missing required columns: expression_size"
  )
  
  # Test missing both columns
  baseline_expr_missing_both <- data.frame(
    response_id = c("gene1", "gene2")
  )
  expect_error(
    obtain_expression_dispersion_curve(baseline_expr_missing_both),
    "Missing required columns: relative_expression, expression_size"
  )
})

test_that("obtain_expression_dispersion_curve validates minimum data requirement", {
  # Test with single row (insufficient for curve fitting)
  baseline_expr_single <- data.frame(
    response_id = "gene1",
    relative_expression = 1e-5,
    expression_size = 0.5
  )
  expect_error(
    obtain_expression_dispersion_curve(baseline_expr_single),
    "baseline_expression must contain at least 2 rows for curve fitting"
  )
  
  # Test with empty data frame
  baseline_expr_empty <- data.frame(
    response_id = character(0),
    relative_expression = numeric(0),
    expression_size = numeric(0)
  )
  expect_error(
    obtain_expression_dispersion_curve(baseline_expr_empty),
    "baseline_expression must contain at least 2 rows for curve fitting"
  )
})

test_that("obtain_expression_dispersion_curve handles missing values", {
  # Test with NA in relative_expression (leaving only 1 valid row)
  baseline_expr_na_rel <- data.frame(
    response_id = c("gene1", "gene2", "gene3"),
    relative_expression = c(1e-5, NA, NA),
    expression_size = c(0.5, NA, 1.5)
  )
  expect_error(
    obtain_expression_dispersion_curve(baseline_expr_na_rel),
    "Not enough valid rows remaining after removing NAs for curve fitting"
  )
  
  # Test with NA in expression_size (leaving only 1 valid row)
  baseline_expr_na_size <- data.frame(
    response_id = c("gene1", "gene2", "gene3"),
    relative_expression = c(1e-5, NA, 3e-5),
    expression_size = c(0.5, NA, NA)
  )
  expect_error(
    obtain_expression_dispersion_curve(baseline_expr_na_size),
    "Not enough valid rows remaining after removing NAs for curve fitting"
  )
})

test_that("obtain_expression_dispersion_curve validates positive values", {
  # Test with zero/negative relative_expression
  baseline_expr_neg_rel <- data.frame(
    response_id = c("gene1", "gene2", "gene3"),
    relative_expression = c(-1e-5, 0, 3e-5),
    expression_size = c(0.5, 1.0, 1.5)
  )
  expect_error(
    obtain_expression_dispersion_curve(baseline_expr_neg_rel),
    "All relative_expression and expression_size values must be positive"
  )
  
  # Test with zero/negative expression_size
  baseline_expr_neg_size <- data.frame(
    response_id = c("gene1", "gene2", "gene3"),
    relative_expression = c(1e-5, 2e-5, 3e-5),
    expression_size = c(-0.5, 0, 1.5)
  )
  expect_error(
    obtain_expression_dispersion_curve(baseline_expr_neg_size),
    "All relative_expression and expression_size values must be positive"
  )
})

test_that("obtain_expression_dispersion_curve extrapolation behavior", {
  baseline_expr <- create_test_baseline_expression(20)
  dispersion_curve <- obtain_expression_dispersion_curve(baseline_expr)
  
  # Test extrapolation below range
  min_expr <- min(baseline_expr$relative_expression)
  below_range <- min_expr * 0.1
  pred_below <- dispersion_curve(below_range)
  expect_true(is.finite(pred_below))
  expect_true(pred_below > 0)
  
  # Test extrapolation above range
  max_expr <- max(baseline_expr$relative_expression)
  above_range <- max_expr * 10
  pred_above <- dispersion_curve(above_range)
  expect_true(is.finite(pred_above))
  expect_true(pred_above > 0)
  
  # With rule = 2, extrapolation should use constant values
  # Values below range should get minimum dispersion value
  # Values above range should get maximum dispersion value
})

test_that("obtain_expression_dispersion_curve integration with obtain_expression_information output", {
  # This test ensures compatibility with the actual output format
  # from obtain_expression_information
  
  # Simulate realistic output from obtain_expression_information
  simulated_output <- data.frame(
    response_id = paste0("ENSG", sprintf("%011d", sample(100000:999999, 50))),
    relative_expression = sort(10^runif(50, -6, -3)),  # Realistic TPM scale
    expression_size = rgamma(50, shape = 2, rate = 2),  # Realistic theta values
    stringsAsFactors = FALSE
  )
  
  # Should work without errors
  expect_silent({
    dispersion_curve <- obtain_expression_dispersion_curve(simulated_output)
  })
  
  # Test prediction on realistic expression values
  dispersion_curve <- obtain_expression_dispersion_curve(simulated_output)
  test_expressions <- 10^runif(10, -6, -3)
  
  expect_silent({
    predictions <- dispersion_curve(test_expressions)
  })
  
  predictions <- dispersion_curve(test_expressions)
  expect_true(all(is.finite(predictions)))
  expect_true(all(predictions > 0))
  expect_equal(length(predictions), length(test_expressions))
})