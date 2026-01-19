# Tests for library_size_curves.cpp functions
# This file tests the C++ saturation curve and library size range functions

# ============================================================================
# Tests for fit_read_UMI_curve_cpp
# ============================================================================

test_that("fit_read_UMI_curve_cpp works with basic inputs", {
  reads <- c(1000, 5000, 10000, 20000)
  UMI_per_cell <- 10000
  variation <- 0.3

  result <- fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation)

  # Should return numeric vector of same length
  expect_type(result, "double")
  expect_length(result, length(reads))

  # UMI counts should be positive and less than UMI_per_cell
  expect_true(all(result > 0))
  expect_true(all(result < UMI_per_cell))

  # UMI should increase with reads
  expect_true(all(diff(result) > 0))
})

test_that("fit_read_UMI_curve_cpp handles single value correctly", {
  reads <- 5000
  UMI_per_cell <- 15000
  variation <- 0.5

  result <- fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result > 0)
  expect_true(result < UMI_per_cell)
})

test_that("fit_read_UMI_curve_cpp validates UMI_per_cell", {
  reads <- c(1000, 5000)
  variation <- 0.3

  # Zero UMI_per_cell
  expect_error(
    fit_read_UMI_curve_cpp(reads, UMI_per_cell = 0, variation),
    "UMI_per_cell must be positive"
  )

  # Negative UMI_per_cell
  expect_error(
    fit_read_UMI_curve_cpp(reads, UMI_per_cell = -1000, variation),
    "UMI_per_cell must be positive"
  )
})

test_that("fit_read_UMI_curve_cpp validates variation", {
  reads <- c(1000, 5000)
  UMI_per_cell <- 10000

  # Zero variation
  expect_error(
    fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation = 0),
    "variation must be positive"
  )

  # Negative variation
  expect_error(
    fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation = -0.1),
    "variation must be positive"
  )
})

test_that("fit_read_UMI_curve_cpp validates reads_per_cell", {
  UMI_per_cell <- 10000
  variation <- 0.3

  # Negative reads
  expect_error(
    fit_read_UMI_curve_cpp(c(1000, -500), UMI_per_cell, variation),
    "reads_per_cell values must be non-negative"
  )
})

test_that("fit_read_UMI_curve_cpp handles zero reads", {
  reads <- c(0, 1000, 5000)
  UMI_per_cell <- 10000
  variation <- 0.3

  result <- fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation)

  # Zero reads should give zero UMI
  expect_equal(result[1], 0)
  expect_true(result[2] > 0)
  expect_true(result[3] > result[2])
})

test_that("fit_read_UMI_curve_cpp handles very high reads", {
  reads <- c(100000, 500000, 1000000)
  UMI_per_cell <- 10000
  variation <- 0.3

  result <- fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation)

  # At very high reads, should approach but not exceed UMI_per_cell
  expect_true(all(result < UMI_per_cell))
  expect_true(all(result > 0.9 * UMI_per_cell))  # Should be close to saturation
})

test_that("fit_read_UMI_curve_cpp handles various UMI_per_cell values", {
  reads <- c(5000, 10000, 20000)
  variation <- 0.3

  # Low UMI_per_cell
  result_low <- fit_read_UMI_curve_cpp(reads, UMI_per_cell = 5000, variation)
  expect_true(all(result_low < 5000))

  # Medium UMI_per_cell
  result_med <- fit_read_UMI_curve_cpp(reads, UMI_per_cell = 15000, variation)
  expect_true(all(result_med < 15000))

  # High UMI_per_cell
  result_high <- fit_read_UMI_curve_cpp(reads, UMI_per_cell = 50000, variation)
  expect_true(all(result_high < 50000))

  # Higher UMI_per_cell should give higher UMI counts
  expect_true(all(result_high > result_med))
  expect_true(all(result_med > result_low))
})

test_that("fit_read_UMI_curve_cpp handles various variation values", {
  reads <- c(5000, 10000, 20000)
  UMI_per_cell <- 10000

  # Low variation
  result_low_var <- fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation = 0.1)

  # Medium variation
  result_med_var <- fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation = 0.5)

  # High variation
  result_high_var <- fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation = 2.0)

  # All should be valid
  expect_true(all(result_low_var > 0))
  expect_true(all(result_med_var > 0))
  expect_true(all(result_high_var > 0))

  # All should be below saturation
  expect_true(all(result_low_var < UMI_per_cell))
  expect_true(all(result_med_var < UMI_per_cell))
  expect_true(all(result_high_var < UMI_per_cell))
})

test_that("fit_read_UMI_curve_cpp produces realistic saturation curves", {
  # Realistic K562 parameters from pilot data
  reads <- seq(1000, 50000, by = 5000)
  UMI_per_cell <- 60414
  variation <- 0.467

  result <- fit_read_UMI_curve_cpp(reads, UMI_per_cell, variation)

  # Check saturation behavior
  expect_true(all(result > 0))
  expect_true(all(result < UMI_per_cell))
  expect_true(all(diff(result) > 0))  # Monotonically increasing

  # At low reads, should have low UMI
  expect_true(result[1] < 0.2 * UMI_per_cell)

  # At high reads (50000), should have meaningful saturation
  # With high UMI_per_cell and moderate variation, 50k reads gives moderate saturation
  expect_true(result[length(result)] > 0.1 * UMI_per_cell)

  # Test with more reads to ensure higher saturation
  high_reads <- 200000
  high_result <- fit_read_UMI_curve_cpp(high_reads, UMI_per_cell, variation)
  expect_true(high_result > 0.5 * UMI_per_cell)
})

# ============================================================================
# Tests for identify_library_size_range_cpp
# ============================================================================

test_that("identify_library_size_range_cpp works with normal parameters", {
  result <- identify_library_size_range_cpp(
    experimental_platform = "K562",
    UMI_per_cell = 10000,
    variation = 0.3
  )

  expect_type(result, "list")
  expect_named(result, c("min_reads_per_cell", "max_reads_per_cell"))

  min_reads <- result$min_reads_per_cell
  max_reads <- result$max_reads_per_cell

  # Should return positive integers
  expect_true(min_reads > 0)
  expect_true(max_reads > 0)

  # Min should be less than max
  expect_true(min_reads < max_reads)

  # Should be reasonable values
  expect_true(min_reads >= 100)  # Minimum search starts at 100
  expect_true(max_reads <= 10 * 10000)  # Max upper bound is 10 * UMI_per_cell
})

test_that("identify_library_size_range_cpp validates UMI_per_cell", {
  # Zero UMI_per_cell
  expect_error(
    identify_library_size_range_cpp("K562", UMI_per_cell = 0, variation = 0.3),
    "UMI_per_cell must be positive"
  )

  # Negative UMI_per_cell
  expect_error(
    identify_library_size_range_cpp("K562", UMI_per_cell = -1000, variation = 0.3),
    "UMI_per_cell must be positive"
  )
})

test_that("identify_library_size_range_cpp validates variation", {
  # Negative variation
  expect_error(
    identify_library_size_range_cpp("K562", UMI_per_cell = 10000, variation = -0.1),
    "variation must be non-negative"
  )
})

test_that("identify_library_size_range_cpp handles various UMI_per_cell values", {
  variation <- 0.3

  # Low UMI_per_cell
  result_low <- identify_library_size_range_cpp("K562", 5000, variation)
  expect_true(result_low$min_reads_per_cell < result_low$max_reads_per_cell)

  # Medium UMI_per_cell
  result_med <- identify_library_size_range_cpp("K562", 15000, variation)
  expect_true(result_med$min_reads_per_cell < result_med$max_reads_per_cell)

  # High UMI_per_cell
  result_high <- identify_library_size_range_cpp("K562", 50000, variation)
  expect_true(result_high$min_reads_per_cell < result_high$max_reads_per_cell)

  # Higher UMI_per_cell should generally require more reads
  expect_true(result_high$max_reads_per_cell > result_low$max_reads_per_cell)
})

test_that("identify_library_size_range_cpp handles corner case with very high variation", {
  # Very high variation might make 98% saturation unreachable
  result <- identify_library_size_range_cpp(
    experimental_platform = "K562",
    UMI_per_cell = 10000,
    variation = 10.0  # Very high variation
  )

  # Should still return valid range
  expect_true(result$min_reads_per_cell > 0)
  expect_true(result$max_reads_per_cell > 0)
  expect_true(result$min_reads_per_cell < result$max_reads_per_cell)
})

# ============================================================================
# Tests for identify_reads_range_cpp (streamlined version)
# ============================================================================

test_that("identify_reads_range_cpp works with normal parameters", {
  result <- identify_reads_range_cpp(UMI_per_cell = 10000, variation = 0.3)

  expect_type(result, "list")
  expect_named(result, c("min_reads_per_cell", "max_reads_per_cell"))

  min_reads <- result$min_reads_per_cell
  max_reads <- result$max_reads_per_cell

  # Should return positive integers
  expect_true(min_reads > 0)
  expect_true(max_reads > 0)

  # Min should be less than max
  expect_true(min_reads < max_reads)
})

test_that("identify_reads_range_cpp validates inputs", {
  # Zero UMI_per_cell
  expect_error(
    identify_reads_range_cpp(UMI_per_cell = 0, variation = 0.3),
    "UMI_per_cell must be positive"
  )

  # Negative variation
  expect_error(
    identify_reads_range_cpp(UMI_per_cell = 10000, variation = -0.1),
    "variation must be non-negative"
  )
})

test_that("identify_reads_range_cpp matches identify_library_size_range_cpp", {
  UMI_per_cell <- 15000
  variation <- 0.4

  result1 <- identify_library_size_range_cpp("K562", UMI_per_cell, variation)
  result2 <- identify_reads_range_cpp(UMI_per_cell, variation)

  # Both should give same results (streamlined version just removes unused parameter)
  expect_equal(result1$min_reads_per_cell, result2$min_reads_per_cell)
  expect_equal(result1$max_reads_per_cell, result2$max_reads_per_cell)
})

# ============================================================================
# Tests for generate_reads_grid_cpp
# ============================================================================

test_that("generate_reads_grid_cpp works with default grid size", {
  result <- generate_reads_grid_cpp(
    experimental_platform = "K562",
    UMI_per_cell = 10000,
    variation = 0.3
  )

  # Should return numeric vector of length 10 (default grid_size)
  expect_type(result, "double")
  expect_length(result, 10)

  # Should be monotonically increasing
  expect_true(all(diff(result) >= 0))

  # All values should be positive
  expect_true(all(result > 0))
})

test_that("generate_reads_grid_cpp works with custom grid sizes", {
  UMI_per_cell <- 10000
  variation <- 0.3

  # Grid size 5
  result_5 <- generate_reads_grid_cpp("K562", UMI_per_cell, variation, grid_size = 5)
  expect_length(result_5, 5)

  # Grid size 20
  result_20 <- generate_reads_grid_cpp("K562", UMI_per_cell, variation, grid_size = 20)
  expect_length(result_20, 20)

  # Both should be monotonically increasing
  expect_true(all(diff(result_5) >= 0))
  expect_true(all(diff(result_20) >= 0))
})

test_that("generate_reads_grid_cpp grid endpoints match range", {
  UMI_per_cell <- 15000
  variation <- 0.4
  grid_size <- 10

  # Get range
  range_result <- identify_library_size_range_cpp("K562", UMI_per_cell, variation)

  # Get grid
  grid_result <- generate_reads_grid_cpp("K562", UMI_per_cell, variation, grid_size)

  # First and last grid points should match min and max from range
  expect_equal(grid_result[1], range_result$min_reads_per_cell)
  expect_equal(grid_result[grid_size], range_result$max_reads_per_cell)
})

# ============================================================================
# Tests for generate_reads_grid_streamlined_cpp
# ============================================================================

test_that("generate_reads_grid_streamlined_cpp works with default grid size", {
  result <- generate_reads_grid_streamlined_cpp(
    UMI_per_cell = 10000,
    variation = 0.3
  )

  # Should return numeric vector of length 10 (default grid_size)
  expect_type(result, "double")
  expect_length(result, 10)

  # Should be monotonically increasing
  expect_true(all(diff(result) >= 0))

  # All values should be positive
  expect_true(all(result > 0))
})

test_that("generate_reads_grid_streamlined_cpp matches generate_reads_grid_cpp", {
  UMI_per_cell <- 12000
  variation <- 0.35
  grid_size <- 15

  result1 <- generate_reads_grid_cpp("K562", UMI_per_cell, variation, grid_size)
  result2 <- generate_reads_grid_streamlined_cpp(UMI_per_cell, variation, grid_size)

  # Both should give same results
  expect_equal(result1, result2)
})

test_that("generate_reads_grid_streamlined_cpp works with various grid sizes", {
  UMI_per_cell <- 10000
  variation <- 0.3

  for (size in c(5, 10, 15, 20)) {
    result <- generate_reads_grid_streamlined_cpp(UMI_per_cell, variation, grid_size = size)
    expect_length(result, size)
    expect_true(all(diff(result) >= 0))
  }
})

# ============================================================================
# Integration tests
# ============================================================================

test_that("Saturation curve and range functions integrate correctly", {
  UMI_per_cell <- 20000
  variation <- 0.5

  # Get range
  range <- identify_reads_range_cpp(UMI_per_cell, variation)

  # Get UMI at min and max reads
  min_UMI <- fit_read_UMI_curve_cpp(range$min_reads_per_cell, UMI_per_cell, variation)
  max_UMI <- fit_read_UMI_curve_cpp(range$max_reads_per_cell, UMI_per_cell, variation)

  # Min should be around 10% saturation (allow some tolerance due to rounding)
  expect_true(min_UMI >= 0.08 * UMI_per_cell)
  expect_true(min_UMI <= 0.12 * UMI_per_cell)

  # Max should be around 98% saturation (or lower if not achievable)
  expect_true(max_UMI >= 0.90 * UMI_per_cell)
  expect_true(max_UMI <= UMI_per_cell)
})

test_that("Grid generation produces evenly spaced reads", {
  UMI_per_cell <- 15000
  variation <- 0.4
  grid_size <- 10

  grid <- generate_reads_grid_streamlined_cpp(UMI_per_cell, variation, grid_size)

  # Check spacing is approximately uniform
  diffs <- diff(grid)
  avg_diff <- mean(diffs)

  # All differences should be close to average (within 10% tolerance)
  expect_true(all(abs(diffs - avg_diff) / avg_diff < 0.1))
})

test_that("Real pilot data parameters work correctly", {
  # Test with actual K562 Gasperini parameters
  UMI_per_cell <- 60414
  variation <- 0.467

  # Get range
  range <- identify_reads_range_cpp(UMI_per_cell, variation)
  expect_true(range$min_reads_per_cell > 0)
  expect_true(range$max_reads_per_cell > range$min_reads_per_cell)

  # Generate grid
  grid <- generate_reads_grid_streamlined_cpp(UMI_per_cell, variation, 10)
  expect_length(grid, 10)
  expect_true(all(diff(grid) > 0))

  # Get UMI values
  UMI_values <- fit_read_UMI_curve_cpp(grid, UMI_per_cell, variation)
  expect_true(all(UMI_values > 0))
  expect_true(all(UMI_values < UMI_per_cell))
})
