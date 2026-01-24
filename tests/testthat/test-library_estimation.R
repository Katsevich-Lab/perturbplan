# Tests for library_estimation function
# This file tests the preseqR-based library parameter estimation

test_that("library_estimation returns correct structure", {
  # Create minimal synthetic QC data
  qc_data <- data.frame(
    num_reads = rep(1:10, times = c(100, 80, 60, 40, 30, 20, 15, 10, 8, 5)),
    UMI_id = 1:368,
    cell_id = rep(paste0("cell_", 1:50), length.out = 368),
    response_id = rep(paste0("gene_", 1:20), length.out = 368)
  )

  result <- library_estimation(qc_data)

  # Check structure
  expect_type(result, "list")
  expect_true("method_used" %in% names(result))
  expect_true("reads_norm" %in% names(result))
  expect_true("n_cells" %in% names(result))
  expect_true("UMI_per_cell_at_saturation" %in% names(result))

  # Check method is either ZTNB or RFA
  expect_true(result$method_used %in% c("ZTNB", "RFA"))

  # Check numeric values are positive
  expect_true(result$reads_norm > 0)
  expect_true(result$n_cells > 0)
  expect_true(result$UMI_per_cell_at_saturation > 0)
})

test_that("library_estimation ZTNB method includes correct parameters", {
  # Create data that should trigger ZTNB method (shape > 1)
  # Use more uniform read distribution
  qc_data <- data.frame(
    num_reads = rep(1:5, times = c(50, 60, 70, 65, 55)),
    UMI_id = 1:300,
    cell_id = rep(paste0("cell_", 1:30), length.out = 300),
    response_id = rep(paste0("gene_", 1:15), length.out = 300)
  )

  result <- library_estimation(qc_data)

  if (result$method_used == "ZTNB") {
    # Check ZTNB-specific parameters
    expect_true("L" %in% names(result))
    expect_true("size" %in% names(result))
    expect_true("mu" %in% names(result))

    expect_true(result$L > 0)
    expect_true(result$size > 0)
    expect_true(result$mu > 0)
  }
})

test_that("library_estimation RFA method includes correct parameters", {
  # Create data with highly skewed distribution (shape <= 1)
  # This should trigger RFA method
  qc_data <- data.frame(
    num_reads = rep(1:20, times = rev(seq(200, 10, length.out = 20))),
    UMI_id = 1:2100,
    cell_id = rep(paste0("cell_", 1:100), length.out = 2100),
    response_id = rep(paste0("gene_", 1:50), length.out = 2100)
  )

  result <- library_estimation(qc_data)

  if (result$method_used == "RFA") {
    # Check RFA-specific parameters
    expect_true("valid_estimator" %in% names(result))

    if (result$valid_estimator) {
      expect_true("coefs_real" %in% names(result))
      expect_true("coefs_imag" %in% names(result))
      expect_true("poles_real" %in% names(result))
      expect_true("poles_imag" %in% names(result))

      expect_type(result$coefs_real, "double")
      expect_type(result$coefs_imag, "double")
      expect_type(result$poles_real, "double")
      expect_type(result$poles_imag, "double")
    } else {
      expect_true("constant_value" %in% names(result))
      expect_true(result$constant_value > 0)
    }
  }
})

test_that("library_estimation output works with fit_read_UMI_curve_cpp", {
  # Create synthetic QC data
  qc_data <- data.frame(
    num_reads = rep(1:10, times = c(100, 80, 60, 40, 30, 20, 15, 10, 8, 5)),
    UMI_id = 1:368,
    cell_id = rep(paste0("cell_", 1:50), length.out = 368),
    response_id = rep(paste0("gene_", 1:20), length.out = 368)
  )

  # Get library parameters
  lib_params <- library_estimation(qc_data)

  # Test that we can predict library sizes
  reads <- c(1000, 5000, 10000, 20000)
  result <- fit_read_UMI_curve_cpp(reads, lib_params)

  # Check basic properties
  expect_type(result, "double")
  expect_length(result, length(reads))
  expect_true(all(result > 0))
  expect_true(all(diff(result) > 0))  # Monotonically increasing

  # UMI should not exceed saturation
  expect_true(all(result < lib_params$UMI_per_cell_at_saturation))
})

test_that("library_estimation handles different mt parameter", {
  qc_data <- data.frame(
    num_reads = rep(1:10, times = c(100, 80, 60, 40, 30, 20, 15, 10, 8, 5)),
    UMI_id = 1:368,
    cell_id = rep(paste0("cell_", 1:50), length.out = 368),
    response_id = rep(paste0("gene_", 1:20), length.out = 368)
  )

  # Test with different mt values
  result_10 <- library_estimation(qc_data, mt = 10)
  result_20 <- library_estimation(qc_data, mt = 20)
  result_30 <- library_estimation(qc_data, mt = 30)

  # All should return valid structures
  expect_type(result_10, "list")
  expect_type(result_20, "list")
  expect_type(result_30, "list")

  # UMI_per_cell_at_saturation should be similar but not necessarily identical
  expect_true(result_10$UMI_per_cell_at_saturation > 0)
  expect_true(result_20$UMI_per_cell_at_saturation > 0)
  expect_true(result_30$UMI_per_cell_at_saturation > 0)
})

test_that("library_estimation with real-like data produces reasonable values", {
  # Simulate more realistic perturb-seq data
  set.seed(42)
  n_umis <- 5000
  n_cells <- 100

  # Create read counts following a realistic distribution
  # Most UMIs have 1-3 reads, some have higher counts (PCR duplicates)
  read_counts <- sample(1:10, n_umis, replace = TRUE,
                       prob = c(0.5, 0.25, 0.15, 0.05, 0.03, 0.01, 0.005, 0.003, 0.001, 0.001))

  qc_data <- data.frame(
    num_reads = read_counts,
    UMI_id = 1:n_umis,
    cell_id = rep(paste0("cell_", 1:n_cells), length.out = n_umis),
    response_id = sample(paste0("gene_", 1:500), n_umis, replace = TRUE)
  )

  result <- library_estimation(qc_data)

  # Check reasonable ranges for perturb-seq data
  expect_true(result$UMI_per_cell_at_saturation >= n_umis / n_cells)  # At least observed UMIs per cell
  expect_true(result$UMI_per_cell_at_saturation <= 10 * n_umis / n_cells)  # Not unreasonably high
  expect_true(result$reads_norm > 0)
  expect_equal(result$n_cells, n_cells)
})
