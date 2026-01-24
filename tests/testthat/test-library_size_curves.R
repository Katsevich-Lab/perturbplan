# Tests for library_size_curves.cpp functions
# This file tests the C++ saturation curve and library size range functions

# Helper function to create ZTNB wrapper for testing
# Uses realistic parameter scaling based on K562 pilot data
create_test_wrapper <- function(UMI_per_cell, variation, reads_norm = NULL, n_cells = NULL) {
  # Use realistic defaults if not specified
  if (is.null(reads_norm)) {
    # Scale reads_norm proportionally to UMI_per_cell
    # K562 has reads_norm=23835, UMI=60414, ratio ~0.395
    reads_norm <- UMI_per_cell * 0.4
  }
  if (is.null(n_cells)) {
    # Use a reasonable default cell count
    n_cells <- 1000
  }

  list(
    method_used = "ZTNB",
    L = UMI_per_cell * n_cells,
    size = 1.0 / variation,
    mu = 0.4,  # Use realistic mu value (K562 uses 0.395)
    reads_norm = reads_norm,
    n_cells = n_cells,
    UMI_per_cell_at_saturation = UMI_per_cell
  )
}

# ============================================================================
# Tests for fit_read_UMI_curve_cpp (new interface with rSAC_fn_wrapper)
# ============================================================================

test_that("fit_read_UMI_curve_cpp works with basic inputs", {
  reads <- c(1000, 5000, 10000, 20000)
  UMI_per_cell <- 10000
  variation <- 0.3

  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  result <- fit_read_UMI_curve_cpp(reads, wrapper)

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

  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  result <- fit_read_UMI_curve_cpp(reads, wrapper)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result > 0)
  expect_true(result < UMI_per_cell)
})

test_that("fit_read_UMI_curve_cpp handles zero reads", {
  reads <- c(0, 1000, 5000)
  UMI_per_cell <- 10000
  variation <- 0.3

  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  result <- fit_read_UMI_curve_cpp(reads, wrapper)

  # Zero reads should give zero or near-zero UMI
  expect_true(result[1] < 0.01 * UMI_per_cell)
  expect_true(result[2] > 0)
  expect_true(result[3] > result[2])
})

test_that("fit_read_UMI_curve_cpp handles very high reads", {
  reads <- c(100000, 500000, 1000000)
  UMI_per_cell <- 10000
  variation <- 0.3

  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  result <- fit_read_UMI_curve_cpp(reads, wrapper)

  # At very high reads, should approach but not exceed UMI_per_cell
  expect_true(all(result < UMI_per_cell))
  expect_true(all(result > 0.9 * UMI_per_cell))  # Should be close to saturation
})

test_that("fit_read_UMI_curve_cpp handles various UMI_per_cell values", {
  reads <- c(5000, 10000, 20000)
  variation <- 0.3

  # Low UMI_per_cell
  wrapper_low <- create_test_wrapper(5000, variation)
  result_low <- fit_read_UMI_curve_cpp(reads, wrapper_low)
  expect_true(all(result_low < 5000))

  # Medium UMI_per_cell
  wrapper_med <- create_test_wrapper(15000, variation)
  result_med <- fit_read_UMI_curve_cpp(reads, wrapper_med)
  expect_true(all(result_med < 15000))

  # High UMI_per_cell
  wrapper_high <- create_test_wrapper(50000, variation)
  result_high <- fit_read_UMI_curve_cpp(reads, wrapper_high)
  expect_true(all(result_high < 50000))

  # Higher UMI_per_cell should give higher UMI counts
  expect_true(all(result_high > result_med))
  expect_true(all(result_med > result_low))
})

test_that("fit_read_UMI_curve_cpp handles various variation values", {
  reads <- c(5000, 10000, 20000)
  UMI_per_cell <- 10000

  # Low variation
  wrapper_low_var <- create_test_wrapper(UMI_per_cell, 0.1)
  result_low_var <- fit_read_UMI_curve_cpp(reads, wrapper_low_var)

  # Medium variation
  wrapper_med_var <- create_test_wrapper(UMI_per_cell, 0.5)
  result_med_var <- fit_read_UMI_curve_cpp(reads, wrapper_med_var)

  # High variation
  wrapper_high_var <- create_test_wrapper(UMI_per_cell, 2.0)
  result_high_var <- fit_read_UMI_curve_cpp(reads, wrapper_high_var)

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

  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  result <- fit_read_UMI_curve_cpp(reads, wrapper)

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
  high_result <- fit_read_UMI_curve_cpp(high_reads, wrapper)
  expect_true(high_result > 0.5 * UMI_per_cell)
})

# ============================================================================
# Tests for fit_read_UMI_curve_cpp as R-exported function
# ============================================================================

test_that("fit_read_UMI_curve_cpp works as exported R function", {
  reads <- c(1000, 5000, 10000, 20000)
  UMI_per_cell <- 10000
  variation <- 0.3

  # Use wrapper interface
  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  result <- fit_read_UMI_curve_cpp(reads, wrapper)

  expect_type(result, "double")
  expect_length(result, length(reads))
  expect_true(all(result > 0))
  expect_true(all(result < UMI_per_cell))
  expect_true(all(diff(result) > 0))
})

# ============================================================================
# Tests for identify_library_size_range_cpp (updated rSAC_fn_wrapper interface)
# ============================================================================

test_that("identify_library_size_range_cpp works with rSAC_fn_wrapper", {
  UMI_per_cell <- 10000
  variation <- 0.3
  wrapper <- create_test_wrapper(UMI_per_cell, variation)

  result <- identify_library_size_range_cpp("10x", wrapper)

  # Should return list with required elements
  expect_type(result, "list")
  expect_true(all(c("min_reads_per_cell", "max_reads_per_cell") %in% names(result)))

  # Min should be less than max
  expect_true(result$min_reads_per_cell < result$max_reads_per_cell)

  # Both should be positive integers
  expect_true(result$min_reads_per_cell > 0)
  expect_true(result$max_reads_per_cell > 0)

  # Min should correspond to ~10% saturation
  min_reads_vec <- c(result$min_reads_per_cell)
  min_UMI <- fit_read_UMI_curve_cpp(min_reads_vec, wrapper)
  expect_true(min_UMI > 0.05 * UMI_per_cell)  # At least 5% saturation
  expect_true(min_UMI < 0.15 * UMI_per_cell)  # At most 15% saturation

  # Max should correspond to ~80% saturation
  max_reads_vec <- c(result$max_reads_per_cell)
  max_UMI <- fit_read_UMI_curve_cpp(max_reads_vec, wrapper)
  expect_true(max_UMI > 0.7 * UMI_per_cell)  # At least 70% saturation
  expect_true(max_UMI < 0.9 * UMI_per_cell)  # At most 90% saturation
})

test_that("identify_library_size_range_cpp handles various UMI_per_cell values", {
  variation <- 0.3

  # Low UMI_per_cell
  wrapper_low <- create_test_wrapper(5000, variation)
  result_low <- identify_library_size_range_cpp("10x", wrapper_low)
  expect_true(result_low$max_reads_per_cell < 100000)

  # High UMI_per_cell
  wrapper_high <- create_test_wrapper(50000, variation)
  result_high <- identify_library_size_range_cpp("10x", wrapper_high)
  expect_true(result_high$max_reads_per_cell > result_low$max_reads_per_cell)
})

# ============================================================================
# Tests for identify_reads_range_cpp (streamlined version)
# ============================================================================

test_that("identify_reads_range_cpp works with rSAC_fn_wrapper", {
  UMI_per_cell <- 10000
  variation <- 0.3
  wrapper <- create_test_wrapper(UMI_per_cell, variation)

  result <- identify_reads_range_cpp(wrapper)

  # Should return list with required elements
  expect_type(result, "list")
  expect_true(all(c("min_reads_per_cell", "max_reads_per_cell") %in% names(result)))

  # Min should be less than max
  expect_true(result$min_reads_per_cell < result$max_reads_per_cell)

  # Both should be positive
  expect_true(result$min_reads_per_cell > 0)
  expect_true(result$max_reads_per_cell > 0)
})

test_that("identify_reads_range_cpp matches identify_library_size_range_cpp", {
  UMI_per_cell <- 15000
  variation <- 0.4
  wrapper <- create_test_wrapper(UMI_per_cell, variation)

  result1 <- identify_library_size_range_cpp("10x", wrapper)
  result2 <- identify_reads_range_cpp(wrapper)

  # Results should be identical (streamlined version just removes platform param)
  expect_equal(result1$min_reads_per_cell, result2$min_reads_per_cell)
  expect_equal(result1$max_reads_per_cell, result2$max_reads_per_cell)
})

# ============================================================================
# Tests for generate_reads_grid_cpp
# ============================================================================

test_that("generate_reads_grid_cpp generates correct grid size", {
  UMI_per_cell <- 10000
  variation <- 0.3
  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  grid_size <- 5

  result <- generate_reads_grid_cpp("10x", wrapper, grid_size)

  expect_type(result, "double")
  expect_length(result, grid_size)
})

test_that("generate_reads_grid_cpp creates evenly spaced grid", {
  UMI_per_cell <- 10000
  variation <- 0.3
  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  grid_size <- 10

  result <- generate_reads_grid_cpp("10x", wrapper, grid_size)

  # Should be monotonically increasing
  expect_true(all(diff(result) > 0))

  # First element should be near min_reads_per_cell
  range_result <- identify_library_size_range_cpp("10x", wrapper)
  expect_true(abs(result[1] - range_result$min_reads_per_cell) < 100)

  # Last element should be near max_reads_per_cell
  expect_true(abs(result[grid_size] - range_result$max_reads_per_cell) < 100)
})

# ============================================================================
# Tests for generate_reads_grid_streamlined_cpp
# ============================================================================

test_that("generate_reads_grid_streamlined_cpp generates correct grid", {
  UMI_per_cell <- 10000
  variation <- 0.3
  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  grid_size <- 8

  result <- generate_reads_grid_streamlined_cpp(wrapper, grid_size)

  expect_type(result, "double")
  expect_length(result, grid_size)
  expect_true(all(diff(result) > 0))
})

test_that("generate_reads_grid_streamlined_cpp matches generate_reads_grid_cpp", {
  UMI_per_cell <- 15000
  variation <- 0.4
  wrapper <- create_test_wrapper(UMI_per_cell, variation)
  grid_size <- 10

  result1 <- generate_reads_grid_cpp("10x", wrapper, grid_size)
  result2 <- generate_reads_grid_streamlined_cpp(wrapper, grid_size)

  # Results should be identical
  expect_equal(result1, result2)
})
