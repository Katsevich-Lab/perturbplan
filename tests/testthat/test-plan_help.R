library(testthat)

# Test fit_read_UMI_curve function
test_that("fit_read_UMI_curve works correctly", {
  
  # Test basic functionality with valid inputs
  result <- fit_read_UMI_curve(reads_per_cell = 1000, UMI_per_cell = 500, variation = 0.1)
  
  # Should return a numeric value
  expect_type(result, "double")
  expect_length(result, 1)
  
  # Result should be positive and less than UMI_per_cell
  expect_gt(result, 0)
  expect_lt(result, 500)
  
  # Test with different parameter values
  result2 <- fit_read_UMI_curve(reads_per_cell = 2000, UMI_per_cell = 1000, variation = 0.05)
  expect_type(result2, "double")
  expect_gt(result2, 0)
  expect_lt(result2, 1000)
  
  # Test edge cases
  # Very low reads per cell
  result_low <- fit_read_UMI_curve(reads_per_cell = 10, UMI_per_cell = 500, variation = 0.1)
  expect_gt(result_low, 0)
  expect_lt(result_low, 50)  # Should be much smaller
  
  # Very high reads per cell (should approach UMI_per_cell)
  result_high <- fit_read_UMI_curve(reads_per_cell = 10000, UMI_per_cell = 500, variation = 0.1)
  expect_gt(result_high, 400)  # Should be close to UMI_per_cell
  expect_lt(result_high, 500)
  
  # Test monotonicity: more reads should give more UMIs (up to saturation)
  reads_seq <- c(100, 200, 500, 1000)
  results <- sapply(reads_seq, function(r) fit_read_UMI_curve(r, 1000, 0.1))
  
  # Results should be increasing
  for (i in 2:length(results)) {
    expect_gt(results[i], results[i-1])
  }
  
  # Test with zero variation
  result_no_var <- fit_read_UMI_curve(reads_per_cell = 1000, UMI_per_cell = 500, variation = 0)
  expect_type(result_no_var, "double")
  expect_gt(result_no_var, 0)
  
  # Test vectorized inputs
  reads_vec <- c(500, 1000, 1500)
  results_vec <- fit_read_UMI_curve(reads_per_cell = reads_vec, UMI_per_cell = 1000, variation = 0.1)
  expect_length(results_vec, 3)
  expect_true(all(results_vec > 0))
  expect_true(all(results_vec < 1000))
})

# Test the mathematical properties of fit_read_UMI_curve
test_that("fit_read_UMI_curve has correct mathematical properties", {
  
  # Test that the function approaches UMI_per_cell as reads_per_cell increases
  UMI_per_cell <- 1000
  variation <- 0.1
  
  # Very high reads should give close to UMI_per_cell
  result_very_high <- fit_read_UMI_curve(reads_per_cell = 100000, UMI_per_cell = UMI_per_cell, variation = variation)
  expect_lt(abs(result_very_high - UMI_per_cell), UMI_per_cell * 0.05)  # Within 5% of UMI_per_cell
  
  # Test saturation curve behavior
  reads_range <- seq(100, 10000, by = 500)
  results <- sapply(reads_range, function(r) fit_read_UMI_curve(r, UMI_per_cell, variation))
  
  # Should be monotonically increasing
  expect_true(all(diff(results) > 0))
  
  # Should be concave (decreasing slope)
  differences <- diff(results)
  expect_true(all(diff(differences) < 0))
  
  # Test boundary conditions
  # When reads_per_cell = 0, result should be 0
  result_zero <- fit_read_UMI_curve(reads_per_cell = 0, UMI_per_cell = UMI_per_cell, variation = variation)
  expect_equal(result_zero, 0, tolerance = 1e-10)
})

# Test extract_library_info function (conditional on data availability)
test_that("extract_library_info function structure", {
  
  # Skip test if library data file doesn't exist
  library_path <- system.file("extdata/library_info", "Gasperini_library.rds", package = "perturbplan")
  
  if (file.exists(library_path)) {
    # Test with default biological system (K562)
    result <- extract_library_info()
    
    # Should return a list with specific elements
    expect_type(result, "list")
    expect_named(result, c("UMI_per_cell", "variation"))
    
    # Both elements should be numeric
    expect_type(result$UMI_per_cell, "double")
    expect_type(result$variation, "double")
    
    # Values should be positive
    expect_gt(result$UMI_per_cell, 0)
    expect_gt(result$variation, 0)
    
    # UMI_per_cell should be a reasonable value (typically hundreds to thousands)
    expect_gt(result$UMI_per_cell, 100)
    expect_lt(result$UMI_per_cell, 10000)
    
    # Variation should be a small positive value (typically < 1)
    expect_lt(result$variation, 1)
    
    # Test with explicit K562 system
    result_k562 <- extract_library_info(biological_system = "K562")
    expect_identical(result, result_k562)
    
    # Test that the function is deterministic (same inputs give same outputs)
    result2 <- extract_library_info()
    expect_identical(result, result2)
    
  } else {
    skip("Library data file not available - extract_library_info tests skipped")
  }
})

# Test integration with mock data
test_that("fit_read_UMI_curve works with realistic parameters", {
  
  # Use realistic parameter values
  UMI_per_cell <- 2000  # Typical value
  variation <- 0.15     # Typical variation
  
  # Test with realistic read depths
  reads_values <- c(500, 1000, 2000, 5000, 10000)
  
  results <- sapply(reads_values, function(reads) {
    fit_read_UMI_curve(
      reads_per_cell = reads,
      UMI_per_cell = UMI_per_cell,
      variation = variation
    )
  })
  
  # All results should be valid
  expect_true(all(is.finite(results)))
  expect_true(all(results > 0))
  expect_true(all(results < UMI_per_cell))
  
  # Results should be increasing with reads
  for (i in 2:length(results)) {
    expect_gt(results[i], results[i-1])
  }
  
  # At typical sequencing depth (2000 reads), should capture reasonable fraction
  result_typical <- fit_read_UMI_curve(reads_per_cell = 2000, UMI_per_cell = UMI_per_cell, variation = variation)
  expect_gt(result_typical, UMI_per_cell * 0.5)  # Should capture at least 50%
  expect_lt(result_typical, UMI_per_cell * 0.9)  # But not more than 90%
})

# Test function argument handling
test_that("fit_read_UMI_curve handles arguments correctly", {
  
  # Test with missing arguments (should work due to R's argument matching)
  expect_no_error(fit_read_UMI_curve(reads_per_cell = 1000, UMI_per_cell = 500, variation = 0.1))
  
  # Test with named arguments in different order
  result1 <- fit_read_UMI_curve(reads_per_cell = 1000, UMI_per_cell = 500, variation = 0.1)
  result2 <- fit_read_UMI_curve(variation = 0.1, UMI_per_cell = 500, reads_per_cell = 1000)
  expect_equal(result1, result2)
  
  # Test with different variation values
  results_var <- sapply(c(0, 0.05, 0.1, 0.2), function(v) {
    fit_read_UMI_curve(reads_per_cell = 1000, UMI_per_cell = 500, variation = v)
  })
  
  # Higher variation should generally give lower UMI counts (for same read depth)
  expect_gt(results_var[1], results_var[4])  # variation=0 > variation=0.2
})