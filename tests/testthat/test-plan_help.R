library(testthat)

# Test fit_read_UMI_curve_cpp function with rSAC_fn_wrapper
test_that("fit_read_UMI_curve_cpp works correctly with rSAC_fn_wrapper", {

  # Get library parameters from package data (rSAC_fn_wrapper format)
  pilot_data <- get_pilot_data_from_package("K562")
  rSAC_params <- pilot_data$library_parameters
  UMI_at_saturation <- rSAC_params$UMI_per_cell_at_saturation

  # Test basic functionality with valid inputs
  result <- fit_read_UMI_curve_cpp(reads_per_cell = 1000, rSAC_fn_wrapper = rSAC_params)

  # Should return a numeric value
  expect_type(result, "double")
  expect_length(result, 1)

  # Result should be positive and less than UMI_per_cell_at_saturation
  expect_gt(result, 0)
  expect_lt(result, UMI_at_saturation)

  # Test with different read depths
  result2 <- fit_read_UMI_curve_cpp(reads_per_cell = 5000, rSAC_fn_wrapper = rSAC_params)
  expect_type(result2, "double")
  expect_gt(result2, result)  # More reads should give more UMIs
  expect_lt(result2, UMI_at_saturation)

  # Test edge cases
  # Very low reads per cell
  result_low <- fit_read_UMI_curve_cpp(reads_per_cell = 10, rSAC_fn_wrapper = rSAC_params)
  expect_gt(result_low, 0)
  expect_lt(result_low, result)  # Should be much smaller

  # Very high reads per cell (should approach UMI_per_cell_at_saturation)
  result_high <- fit_read_UMI_curve_cpp(reads_per_cell = 500000, rSAC_fn_wrapper = rSAC_params)
  expect_gt(result_high, UMI_at_saturation * 0.8)  # Should be close to saturation
  expect_lt(result_high, UMI_at_saturation * 1.05)  # Within 5% of saturation

  # Test monotonicity: more reads should give more UMIs (up to saturation)
  reads_seq <- c(100, 500, 1000, 5000, 10000)
  results <- sapply(reads_seq, function(r) {
    fit_read_UMI_curve_cpp(reads_per_cell = r, rSAC_fn_wrapper = rSAC_params)
  })

  # Results should be increasing
  for (i in 2:length(results)) {
    expect_gt(results[i], results[i-1])
  }

  # Test vectorized inputs
  reads_vec <- c(500, 1000, 1500, 2000)
  results_vec <- fit_read_UMI_curve_cpp(reads_per_cell = reads_vec, rSAC_fn_wrapper = rSAC_params)
  expect_length(results_vec, length(reads_vec))
  expect_true(all(results_vec > 0))
  expect_true(all(results_vec < UMI_at_saturation))
})

# Test the mathematical properties of fit_read_UMI_curve_cpp
test_that("fit_read_UMI_curve_cpp has correct mathematical properties", {

  # Get library parameters from package data
  pilot_data <- get_pilot_data_from_package("K562")
  rSAC_params <- pilot_data$library_parameters
  UMI_at_saturation <- rSAC_params$UMI_per_cell_at_saturation

  # Very high reads should approach UMI_per_cell_at_saturation
  result_very_high <- fit_read_UMI_curve_cpp(reads_per_cell = 1000000, rSAC_fn_wrapper = rSAC_params)
  expect_lt(abs(result_very_high - UMI_at_saturation) / UMI_at_saturation, 0.1)  # Within 10%

  # Test saturation curve behavior
  reads_range <- seq(1000, 100000, by = 5000)
  results <- sapply(reads_range, function(r) {
    fit_read_UMI_curve_cpp(reads_per_cell = r, rSAC_fn_wrapper = rSAC_params)
  })

  # Should be monotonically increasing
  expect_true(all(diff(results) >= 0))  # Allow for saturation plateau

  # Test boundary conditions
  # When reads_per_cell = 0, result should be 0
  result_zero <- fit_read_UMI_curve_cpp(reads_per_cell = 0, rSAC_fn_wrapper = rSAC_params)
  expect_equal(result_zero, 0, tolerance = 1e-10)
})


# Test integration with realistic data
test_that("fit_read_UMI_curve_cpp works with realistic read depths", {

  # Get library parameters from package data
  pilot_data <- get_pilot_data_from_package("K562")
  rSAC_params <- pilot_data$library_parameters
  UMI_at_saturation <- rSAC_params$UMI_per_cell_at_saturation

  # Test with realistic read depths
  reads_values <- c(5000, 10000, 20000, 50000, 100000)

  results <- sapply(reads_values, function(reads) {
    fit_read_UMI_curve_cpp(reads_per_cell = reads, rSAC_fn_wrapper = rSAC_params)
  })

  # All results should be valid
  expect_true(all(is.finite(results)))
  expect_true(all(results > 0))
  expect_true(all(results <= UMI_at_saturation * 1.05))  # Allow small overshoot

  # Results should be increasing with reads (or plateau at saturation)
  for (i in 2:length(results)) {
    expect_gte(results[i], results[i-1])
  }

  # At typical sequencing depth (20000 reads), should capture reasonable fraction
  result_typical <- fit_read_UMI_curve_cpp(reads_per_cell = 20000, rSAC_fn_wrapper = rSAC_params)
  expect_gt(result_typical, UMI_at_saturation * 0.1)  # Should capture at least 10%
  expect_lt(result_typical, UMI_at_saturation * 1.1)  # Allow small overshoot
})

# ============================================================================
# Tests for get_pilot_data_from_package()
# ============================================================================

test_that("get_pilot_data_from_package loads K562 data correctly", {
  k562_data <- get_pilot_data_from_package("K562")

  # Test that it returns a list
  expect_type(k562_data, "list")

  # Test that required elements exist
  expect_true("baseline_expression_stats" %in% names(k562_data) || "baseline_expression" %in% names(k562_data))
  expect_true("library_parameters" %in% names(k562_data))

  # Test library_parameters structure (rSAC_fn_wrapper format)
  expect_type(k562_data$library_parameters, "list")
  expect_true("method_used" %in% names(k562_data$library_parameters))
  expect_true("UMI_per_cell_at_saturation" %in% names(k562_data$library_parameters))
  expect_true("reads_norm" %in% names(k562_data$library_parameters))
  expect_true("n_cells" %in% names(k562_data$library_parameters))

  # Test library parameter values
  expect_gt(k562_data$library_parameters$UMI_per_cell_at_saturation, 0)
  expect_gt(k562_data$library_parameters$reads_norm, 0)
  expect_gt(k562_data$library_parameters$n_cells, 0)
})

test_that("get_pilot_data_from_package works for different biological systems", {
  # Test various systems
  systems <- c("K562", "A549", "THP-1", "iPSC")

  for (sys in systems) {
    data <- get_pilot_data_from_package(sys)
    expect_type(data, "list")
    expect_true("library_parameters" %in% names(data))
  }
})

test_that("get_pilot_data_from_package errors on invalid system", {
  expect_error(
    get_pilot_data_from_package("InvalidSystem"),
    "Unsupported biological system"
  )
})

# ============================================================================
# Tests for validation functions
# ============================================================================

test_that("validate_custom_baseline accepts valid data", {
  # Create valid baseline data
  valid_data <- data.frame(
    response_id = c("ENSG00000141510", "ENSG00000157764", "ENSG00000186092"),
    relative_expression = c(1.5e-05, 2.3e-05, 5.7e-06),
    expression_size = c(0.5, 1.2, 0.8)
  )

  result <- validate_custom_baseline(valid_data)

  expect_true(result$valid)
  expect_equal(nrow(result$data), 3)
  expect_length(result$errors, 0)
})

test_that("validate_custom_baseline detects missing columns", {
  # Missing expression_size column
  invalid_data <- data.frame(
    response_id = c("ENSG00000141510", "ENSG00000157764"),
    relative_expression = c(1.5e-05, 2.3e-05)
  )

  result <- validate_custom_baseline(invalid_data)

  expect_false(result$valid)
  expect_true(any(grepl("Missing required columns", result$errors)))
})

test_that("validate_custom_baseline detects negative values", {
  # Negative relative_expression
  invalid_data <- data.frame(
    response_id = c("ENSG00000141510", "ENSG00000157764"),
    relative_expression = c(-1.5e-05, 2.3e-05),
    expression_size = c(0.5, 1.2)
  )

  result <- validate_custom_baseline(invalid_data)

  expect_false(result$valid)
  expect_true(any(grepl("non-negative", result$errors)))
})

test_that("validate_custom_baseline detects non-positive expression_size", {
  # Zero expression_size
  invalid_data <- data.frame(
    response_id = c("ENSG00000141510", "ENSG00000157764"),
    relative_expression = c(1.5e-05, 2.3e-05),
    expression_size = c(0, 1.2)
  )

  result <- validate_custom_baseline(invalid_data)

  expect_false(result$valid)
  expect_true(any(grepl("must be positive", result$errors)))
})

test_that("validate_custom_baseline handles duplicates", {
  # Duplicate gene IDs
  dup_data <- data.frame(
    response_id = c("ENSG00000141510", "ENSG00000141510", "ENSG00000157764"),
    relative_expression = c(1.5e-05, 2.3e-05, 5.7e-06),
    expression_size = c(0.5, 1.2, 0.8)
  )

  result <- validate_custom_baseline(dup_data)

  # Should be valid but with warnings, and duplicates removed
  expect_true(result$valid)
  expect_equal(nrow(result$data), 2)  # Duplicate removed
  expect_true(any(grepl("duplicate", result$warnings)))
})

test_that("validate_custom_baseline warns about non-Ensembl IDs", {
  # Non-Ensembl gene IDs
  non_ensembl_data <- data.frame(
    response_id = c("TP53", "BRCA1", "ENSG00000141510"),
    relative_expression = c(1.5e-05, 2.3e-05, 5.7e-06),
    expression_size = c(0.5, 1.2, 0.8)
  )

  result <- validate_custom_baseline(non_ensembl_data)

  expect_true(result$valid)
  expect_true(any(grepl("Ensembl format", result$warnings)))
})

test_that("validate_custom_baseline detects NA values", {
  # NA in relative_expression
  na_data <- data.frame(
    response_id = c("ENSG00000141510", "ENSG00000157764"),
    relative_expression = c(NA, 2.3e-05),
    expression_size = c(0.5, 1.2)
  )

  result <- validate_custom_baseline(na_data)

  expect_false(result$valid)
  expect_true(any(grepl("missing values", result$errors)))
})

test_that("validate_custom_baseline_rds accepts valid data frame", {
  # Create valid data frame (new format)
  valid_df <- data.frame(
    response_id = c("ENSG00000141510", "ENSG00000157764", "ENSG00000186092"),
    relative_expression = c(1.5e-05, 2.3e-05, 5.7e-06),
    expression_size = c(0.5, 1.2, 0.8)
  )

  result <- validate_custom_baseline_rds(valid_df)

  expect_true(result$valid)
  expect_length(result$errors, 0)
})

test_that("validate_custom_baseline_rds accepts valid list structure", {
  # Create valid list structure
  valid_list <- list(
    baseline_expression_stats = data.frame(
      response_id = c("ENSG00000141510", "ENSG00000157764"),
      relative_expression = c(1.5e-05, 2.3e-05),
      expression_size = c(0.5, 1.2)
    )
  )

  result <- validate_custom_baseline_rds(valid_list)

  expect_true(result$valid)
  expect_length(result$errors, 0)
})

test_that("validate_custom_library_rds accepts valid library parameters", {
  # Get valid library parameters from package data (rSAC_fn_wrapper format)
  pilot_data <- get_pilot_data_from_package("K562")
  valid_library <- pilot_data$library_parameters

  result <- validate_custom_library_rds(valid_library)

  expect_true(result$valid)
  expect_length(result$errors, 0)
  expect_equal(result$data$method_used, valid_library$method_used)
  expect_equal(result$data$UMI_per_cell_at_saturation, valid_library$UMI_per_cell_at_saturation)
})

test_that("validate_custom_library_rds detects missing elements", {
  # Missing required fields for rSAC_fn_wrapper format
  invalid_library <- list(
    method_used = "ZTNB",
    UMI_per_cell_at_saturation = 15000
    # Missing: reads_norm, n_cells
  )

  result <- validate_custom_library_rds(invalid_library)

  expect_false(result$valid)
  expect_true(any(grepl("Missing required elements", result$errors)))
})

test_that("validate_custom_library_rds detects invalid UMI_per_cell_at_saturation", {
  # Negative UMI_per_cell_at_saturation
  invalid_library <- list(
    method_used = "ZTNB",
    UMI_per_cell_at_saturation = -1000,
    reads_norm = 50000,
    n_cells = 10000
  )

  result <- validate_custom_library_rds(invalid_library)

  expect_false(result$valid)
  expect_true(any(grepl("must be positive", result$errors)))
})

test_that("validate_custom_library_rds warns about unusual values", {
  # Very low UMI_per_cell_at_saturation
  low_umi_library <- list(
    method_used = "ZTNB",
    UMI_per_cell_at_saturation = 500,
    reads_norm = 50000,
    n_cells = 10000
  )

  result <- validate_custom_library_rds(low_umi_library)

  expect_true(result$valid)
  expect_true(any(grepl("unusually low", result$warnings)))
})

test_that("validate_combined_pilot_data accepts valid combined data", {
  # Get valid library parameters from package data
  pilot_data <- get_pilot_data_from_package("K562")

  # Create valid combined data
  valid_combined <- list(
    baseline_expression_stats = data.frame(
      response_id = c("ENSG00000141510", "ENSG00000157764"),
      relative_expression = c(1.5e-05, 2.3e-05),
      expression_size = c(0.5, 1.2)
    ),
    library_parameters = pilot_data$library_parameters
  )

  result <- validate_combined_pilot_data(valid_combined)

  expect_true(result$valid)
  expect_length(result$errors, 0)
})

test_that("validate_combined_pilot_data detects missing components", {
  # Missing library_parameters
  invalid_combined <- list(
    baseline_expression_stats = data.frame(
      response_id = c("ENSG00000141510"),
      relative_expression = c(1.5e-05),
      expression_size = c(0.5)
    )
  )

  result <- validate_combined_pilot_data(invalid_combined)

  expect_false(result$valid)
  expect_true(any(grepl("library_parameters", result$errors)))
})

test_that("validate_combined_pilot_data validates nested components", {
  # Invalid baseline data within combined structure
  invalid_combined <- list(
    baseline_expression_stats = data.frame(
      response_id = c("ENSG00000141510"),
      relative_expression = c(-1.5e-05),  # Negative value
      expression_size = c(0.5)
    ),
    library_parameters = list(
      UMI_per_cell = 15000,
      variation = 0.25
    )
  )

  result <- validate_combined_pilot_data(invalid_combined)

  expect_false(result$valid)
  expect_true(any(grepl("non-negative", result$errors)))
})

# ============================================================================
# Tests for extract_fc_expression_info()
# ============================================================================

test_that("extract_fc_expression_info works with default parameters", {
  set.seed(123)

  result <- extract_fc_expression_info(
    minimum_fold_change = 0.8,
    gRNA_variability = 0.1,
    biological_system = "K562",
    B = 50,  # Small B for faster testing
    TPM_threshold = 10
  )

  # Test that it returns a list
  expect_type(result, "list")

  # Test that required elements exist
  expect_true("fc_expression_df" %in% names(result))
  expect_true("pilot_data" %in% names(result))

  # Test fc_expression_df is a data frame
  expect_s3_class(result$fc_expression_df, "data.frame")

  # Test that it has the expected number of rows
  expect_equal(nrow(result$fc_expression_df), 50)

  # Test that required columns exist
  expect_true("avg_fold_change" %in% colnames(result$fc_expression_df))
  expect_true("relative_expression" %in% colnames(result$fc_expression_df))
  expect_true("expression_size" %in% colnames(result$fc_expression_df))

  # Test that avg_fold_change values are around the specified minimum
  expect_true(all(result$fc_expression_df$avg_fold_change > 0))

  # Test that relative_expression is positive
  expect_true(all(result$fc_expression_df$relative_expression > 0))
})

test_that("extract_fc_expression_info works with gene_list", {
  set.seed(123)

  # Get some genes from K562 data
  k562_data <- get_pilot_data_from_package("K562")
  baseline_stats <- if ("baseline_expression_stats" %in% names(k562_data)) {
    k562_data$baseline_expression_stats
  } else {
    k562_data$baseline_expression$baseline_expression
  }

  # Take first 10 genes with reasonable expression
  gene_subset <- head(baseline_stats$response_id, 10)

  result <- extract_fc_expression_info(
    minimum_fold_change = 0.8,
    gRNA_variability = 0.1,
    biological_system = "K562",
    B = 20,  # Request 20 samples
    gene_list = gene_subset,
    TPM_threshold = 0  # Low threshold to ensure genes pass
  )

  # Test that result contains fc_expression_df
  expect_true("fc_expression_df" %in% names(result))
  expect_s3_class(result$fc_expression_df, "data.frame")

  # Test that result contains genes from gene_list
  expect_true(all(result$fc_expression_df$response_id %in% gene_subset))
})

test_that("extract_fc_expression_info respects TPM_threshold", {
  set.seed(123)

  # High threshold should filter out low-expression genes
  result_high <- extract_fc_expression_info(
    minimum_fold_change = 0.8,
    gRNA_variability = 0.1,
    biological_system = "K562",
    B = 50,
    TPM_threshold = 100  # High threshold
  )

  # All genes should have expression above threshold
  expect_true(all(result_high$fc_expression_df$relative_expression * 1e6 >= 100))
})

# ============================================================================
# Tests for identify_library_size_range() - internal function
# ============================================================================

test_that("identify_library_size_range returns valid range", {
  # Get library parameters from package data (rSAC_fn_wrapper format)
  pilot_data <- get_pilot_data_from_package("K562")
  library_params <- pilot_data$library_parameters

  # Access internal function with :::
  result <- perturbplan:::identify_library_size_range("10x Chromium v3", library_params)

  # Test that it returns a list with range values
  expect_type(result, "list")

  # Check for actual returned names from C++ function
  expect_true("min_reads_per_cell" %in% names(result))
  expect_true("max_reads_per_cell" %in% names(result))

  # Test that values are positive and min < max
  expect_gt(result$min_reads_per_cell, 0)
  expect_gt(result$max_reads_per_cell, 0)
  expect_lt(result$min_reads_per_cell, result$max_reads_per_cell)
})