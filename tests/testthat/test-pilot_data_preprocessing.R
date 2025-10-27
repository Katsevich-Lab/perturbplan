# Test for pilot_data_preprocessing.R functions
library(testthat)

# Get path to test data
extdata_path <- system.file("extdata", package = "perturbplan")

test_that("reference_data_preprocessing_10x aggregates data correctly", {
  # Run preprocessing with h5_rough = TRUE
  raw_data <- reference_data_preprocessing_10x(
    path_to_top_level_output = extdata_path,
    path_to_run_level_output = "cellranger_tiny",
    h5_rough = TRUE,
    skip_mapping_efficiency = FALSE
  )

  # Test that it returns a list with required elements
  expect_type(raw_data, "list")
  expect_true(all(c("response_matrix", "read_umi_table", "mapping_efficiency") %in% names(raw_data)))

  # Test response_matrix
  expect_s4_class(raw_data$response_matrix, "CsparseMatrix")
  expect_gt(nrow(raw_data$response_matrix), 0)
  expect_gt(ncol(raw_data$response_matrix), 0)

  # Test read_umi_table
  expect_s3_class(raw_data$read_umi_table, "data.frame")
  required_cols <- c("num_reads", "UMI_id", "cell_id", "response_id", "srr_idx")
  expect_true(all(required_cols %in% colnames(raw_data$read_umi_table)))

  # Test mapping_efficiency
  expect_type(raw_data$mapping_efficiency, "double")
  expect_length(raw_data$mapping_efficiency, 1)
  expect_gte(raw_data$mapping_efficiency, 0)
  expect_lte(raw_data$mapping_efficiency, 1)
})

test_that("reference_data_preprocessing_10x with skip_mapping_efficiency works", {
  # When h5_rough = FALSE and skip_mapping_efficiency = TRUE, mapping_efficiency should be NULL
  raw_data <- reference_data_preprocessing_10x(
    path_to_top_level_output = extdata_path,
    path_to_run_level_output = "cellranger_tiny",
    h5_rough = FALSE,
    skip_mapping_efficiency = TRUE
  )

  # Test that mapping_efficiency is NULL
  expect_null(raw_data$mapping_efficiency)

  # Other components should still exist
  expect_s4_class(raw_data$response_matrix, "CsparseMatrix")
  expect_s3_class(raw_data$read_umi_table, "data.frame")
})

test_that("reference_data_preprocessing_10x with h5_rough = FALSE works", {
  # Run preprocessing with h5_rough = FALSE
  raw_data <- reference_data_preprocessing_10x(
    path_to_top_level_output = extdata_path,
    path_to_run_level_output = "cellranger_tiny",
    h5_rough = FALSE,
    skip_mapping_efficiency = FALSE
  )

  # Test that it returns valid data
  expect_type(raw_data, "list")
  expect_s4_class(raw_data$response_matrix, "CsparseMatrix")
  expect_s3_class(raw_data$read_umi_table, "data.frame")

  # Mapping efficiency should be median of all runs
  expect_type(raw_data$mapping_efficiency, "double")
  expect_length(raw_data$mapping_efficiency, 1)
})

test_that("reference_data_preprocessing_10x handles warnings for missing directories", {
  # Test with non-existent directory in path_to_run_level_output
  expect_warning(
    reference_data_preprocessing_10x(
      path_to_top_level_output = extdata_path,
      path_to_run_level_output = c("cellranger_tiny", "nonexistent_dir"),
      h5_rough = TRUE
    ),
    "not found in the run directories"
  )
})

test_that("reference_data_preprocessing_10x errors with no valid directories", {
  # Test with only non-existent directories
  expect_error(
    reference_data_preprocessing_10x(
      path_to_top_level_output = extdata_path,
      path_to_run_level_output = "nonexistent_dir",
      h5_rough = TRUE
    ),
    "No valid run directories found"
  )
})

test_that("reference_data_processing extracts parameters correctly", {
  # Set seed for reproducibility
  set.seed(123)

  # First get raw data
  raw_data <- reference_data_preprocessing_10x(
    path_to_top_level_output = extdata_path,
    path_to_run_level_output = "cellranger_tiny",
    h5_rough = TRUE,
    skip_mapping_efficiency = FALSE
  )

  # Process the data
  pilot_data <- reference_data_processing(
    response_matrix = raw_data$response_matrix,
    read_umi_table = raw_data$read_umi_table,
    mapping_efficiency = raw_data$mapping_efficiency,
    gene_list = NULL,
    TPM_thres = 0.1,
    downsample_ratio = 0.6,
    D2_rough = 0.4,
    h5_only = FALSE,
    n_threads = 1
  )

  # Test that it returns a list with required elements
  expect_type(pilot_data, "list")
  expect_true(all(c("baseline_expression_stats", "library_parameters", "mapping_efficiency") %in% names(pilot_data)))

  # Test baseline_expression_stats
  expect_s3_class(pilot_data$baseline_expression_stats, "data.frame")
  baseline_cols <- c("response_id", "relative_expression", "expression_size")
  expect_true(all(baseline_cols %in% colnames(pilot_data$baseline_expression_stats)))

  # Test relative_expression properties
  expect_true(all(pilot_data$baseline_expression_stats$relative_expression > 0))
  expect_equal(sum(pilot_data$baseline_expression_stats$relative_expression), 1, tolerance = 1e-6)

  # Test expression_size is positive
  expect_true(all(pilot_data$baseline_expression_stats$expression_size > 0))

  # Test library_parameters
  expect_type(pilot_data$library_parameters, "list")
  expect_true(all(c("UMI_per_cell", "variation") %in% names(pilot_data$library_parameters)))

  # Test UMI_per_cell is positive
  expect_gt(pilot_data$library_parameters$UMI_per_cell, 0)

  # Test variation is between 0 and 1
  expect_gte(pilot_data$library_parameters$variation, 0)
  expect_lte(pilot_data$library_parameters$variation, 1)

  # Test mapping_efficiency
  expect_type(pilot_data$mapping_efficiency, "double")
  expect_gte(pilot_data$mapping_efficiency, 0)
  expect_lte(pilot_data$mapping_efficiency, 1)
})

test_that("reference_data_processing with h5_only = TRUE skips baseline expression", {
  # Set seed for reproducibility
  set.seed(123)

  # Get raw data
  raw_data <- reference_data_preprocessing_10x(
    path_to_top_level_output = extdata_path,
    path_to_run_level_output = "cellranger_tiny",
    h5_rough = TRUE
  )

  # Process with h5_only = TRUE
  pilot_data <- reference_data_processing(
    response_matrix = raw_data$response_matrix,
    read_umi_table = raw_data$read_umi_table,
    mapping_efficiency = raw_data$mapping_efficiency,
    h5_only = TRUE,
    n_threads = 1
  )

  # baseline_expression_stats should be NULL
  expect_null(pilot_data$baseline_expression_stats)

  # library_parameters and mapping_efficiency should still exist
  expect_type(pilot_data$library_parameters, "list")
  expect_true(all(c("UMI_per_cell", "variation") %in% names(pilot_data$library_parameters)))
})

test_that("reference_data_processing with gene_list filters correctly", {
  # Set seed for reproducibility
  set.seed(123)

  # Get raw data
  raw_data <- reference_data_preprocessing_10x(
    path_to_top_level_output = extdata_path,
    path_to_run_level_output = "cellranger_tiny",
    h5_rough = TRUE
  )

  # Get a subset of genes
  all_genes <- rownames(raw_data$response_matrix)
  gene_subset <- head(all_genes, min(5, length(all_genes)))

  # Process with gene_list
  pilot_data <- reference_data_processing(
    response_matrix = raw_data$response_matrix,
    read_umi_table = raw_data$read_umi_table,
    mapping_efficiency = raw_data$mapping_efficiency,
    gene_list = gene_subset,
    TPM_thres = 0,  # Use low threshold to ensure genes pass
    downsample_ratio = 0.6,
    D2_rough = 0.4,
    h5_only = FALSE,
    n_threads = 1
  )

  # Test that only genes in gene_list are included
  expect_true(all(pilot_data$baseline_expression_stats$response_id %in% gene_subset))

  # Test that mapping efficiency was adjusted
  expect_type(pilot_data$mapping_efficiency, "double")
})

test_that("reference_data_processing respects TPM_thres parameter", {
  # Set seed for reproducibility
  set.seed(123)

  # Get raw data
  raw_data <- reference_data_preprocessing_10x(
    path_to_top_level_output = extdata_path,
    path_to_run_level_output = "cellranger_tiny",
    h5_rough = TRUE
  )

  # Process with high TPM threshold
  pilot_data_high <- reference_data_processing(
    response_matrix = raw_data$response_matrix,
    read_umi_table = raw_data$read_umi_table,
    mapping_efficiency = raw_data$mapping_efficiency,
    TPM_thres = 100,  # High threshold
    downsample_ratio = 0.6,
    h5_only = FALSE,
    n_threads = 1
  )

  # Process with low TPM threshold
  pilot_data_low <- reference_data_processing(
    response_matrix = raw_data$response_matrix,
    read_umi_table = raw_data$read_umi_table,
    mapping_efficiency = raw_data$mapping_efficiency,
    TPM_thres = 0.01,  # Low threshold
    downsample_ratio = 0.6,
    h5_only = FALSE,
    n_threads = 1
  )

  # High threshold should give fewer genes than low threshold
  expect_lte(
    nrow(pilot_data_high$baseline_expression_stats),
    nrow(pilot_data_low$baseline_expression_stats)
  )
})
