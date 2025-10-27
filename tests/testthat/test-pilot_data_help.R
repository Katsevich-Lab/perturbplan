# Test for pilot_data_help.R functions
library(testthat)

# Get path to test data
cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")

test_that("obtain_qc_response_data loads Cell Ranger matrix correctly", {
  # Load the matrix
  response_matrix <- obtain_qc_response_data(cellranger_path)

  # Test that it returns a sparse matrix
  expect_s4_class(response_matrix, "CsparseMatrix")

  # Test that dimensions are positive
  expect_gt(nrow(response_matrix), 0)
  expect_gt(ncol(response_matrix), 0)

  # Test that row and column names exist
  expect_false(is.null(rownames(response_matrix)))
  expect_false(is.null(colnames(response_matrix)))

  # Test that gene IDs are non-empty and unique
  expect_true(all(nzchar(rownames(response_matrix))))
  expect_equal(length(unique(rownames(response_matrix))), nrow(response_matrix))

  # Test that barcodes are non-empty and unique
  expect_true(all(nzchar(colnames(response_matrix))))
  expect_equal(length(unique(colnames(response_matrix))), ncol(response_matrix))

  # Test that matrix contains non-negative values
  expect_true(all(response_matrix@x >= 0))
})

test_that("obtain_qc_response_data handles invalid paths", {
  # Test with non-existent path
  expect_error(
    obtain_qc_response_data("/path/to/nowhere"),
    "cannot open"
  )
})

test_that("obtain_qc_read_umi_table extracts UMI information correctly", {
  # Load UMI table
  read_umi_table <- obtain_qc_read_umi_table(cellranger_path)

  # Test that it returns a data frame
  expect_s3_class(read_umi_table, "data.frame")

  # Test that required columns exist
  required_cols <- c("num_reads", "UMI_id", "cell_id", "response_id")
  expect_true(all(required_cols %in% colnames(read_umi_table)))

  # Test that num_reads is positive
  expect_true(all(read_umi_table$num_reads > 0))

  # Test that there are rows of data
  expect_gt(nrow(read_umi_table), 0)

  # Test that cell_id and response_id are character/factor
  expect_true(is.character(read_umi_table$cell_id) || is.factor(read_umi_table$cell_id))
  expect_true(is.character(read_umi_table$response_id) || is.factor(read_umi_table$response_id))
})

test_that("obtain_qc_read_umi_table handles invalid paths", {
  # Test with non-existent path
  expect_error(
    obtain_qc_read_umi_table("/path/to/nowhere"),
    "does not exist"
  )
})

test_that("obtain_mapping_efficiency calculates efficiency correctly", {
  # Load UMI table
  read_umi_table <- obtain_qc_read_umi_table(cellranger_path)

  # Calculate mapping efficiency
  mapping_eff <- obtain_mapping_efficiency(read_umi_table, cellranger_path)

  # Test that it returns a numeric value
  expect_type(mapping_eff, "double")
  expect_length(mapping_eff, 1)

  # Test that efficiency is between 0 and 1
  expect_gte(mapping_eff, 0)
  expect_lte(mapping_eff, 1)

  # Test that it's not NA
  expect_false(is.na(mapping_eff))
})

test_that("obtain_expression_information fits NB model correctly", {
  # Load response matrix
  response_matrix <- obtain_qc_response_data(cellranger_path)

  # Extract expression information
  expr_info <- obtain_expression_information(
    response_matrix = response_matrix,
    TPM_thres = 0.1,
    rough = TRUE,
    n_threads = 1
  )

  # Test that it returns a data frame
  expect_s3_class(expr_info, "data.frame")

  # Test that required columns exist
  required_cols <- c("response_id", "relative_expression", "expression_size")
  expect_true(all(required_cols %in% colnames(expr_info)))

  # Test that there are rows (some genes should pass TPM threshold)
  expect_gt(nrow(expr_info), 0)

  # Test that relative_expression is positive and sums to approximately 1
  expect_true(all(expr_info$relative_expression > 0))
  expect_equal(sum(expr_info$relative_expression), 1, tolerance = 1e-6)

  # Test that expression_size is positive
  expect_true(all(expr_info$expression_size > 0))

  # Test that response_id is non-empty
  expect_true(all(nzchar(expr_info$response_id)))
})

test_that("obtain_expression_information respects TPM threshold", {
  # Load response matrix
  response_matrix <- obtain_qc_response_data(cellranger_path)

  # Test with high threshold (should return fewer genes)
  expr_info_high <- obtain_expression_information(
    response_matrix = response_matrix,
    TPM_thres = 100,  # High threshold
    rough = TRUE,
    n_threads = 1
  )

  # Test with low threshold (should return more genes)
  expr_info_low <- obtain_expression_information(
    response_matrix = response_matrix,
    TPM_thres = 0.01,  # Low threshold
    rough = TRUE,
    n_threads = 1
  )

  # High threshold should give fewer genes than low threshold
  expect_lte(nrow(expr_info_high), nrow(expr_info_low))
})

test_that("obtain_expression_information handles edge cases", {
  # Load response matrix
  response_matrix <- obtain_qc_response_data(cellranger_path)

  # Test with TPM_thres = 0 (should include all genes with non-zero expression)
  expr_info_zero <- obtain_expression_information(
    response_matrix = response_matrix,
    TPM_thres = 0,
    rough = TRUE,
    n_threads = 1
  )

  expect_gt(nrow(expr_info_zero), 0)
  expect_true(all(expr_info_zero$relative_expression > 0))
})
