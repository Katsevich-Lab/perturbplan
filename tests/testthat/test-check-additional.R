# Additional tests for R/check.R to improve coverage
# These tests focus on edge cases and error conditions

# ============================================================================
# Tests for input_check_library_computation
# ============================================================================

test_that("input_check_library_computation validates QC_data structure", {
  # Test with NULL
  expect_error(
    input_check_library_computation(QC_data = NULL, downsample_ratio = 0.7, D2_rough = 0.3),
    "QC_data.*must be a specified data frame"
  )

  # Test with non-data.frame
  expect_error(
    input_check_library_computation(QC_data = list(), downsample_ratio = 0.7, D2_rough = 0.3),
    "QC_data.*must be a specified data frame"
  )

  # Test with missing columns
  bad_df <- data.frame(wrong_col = 1:10)
  expect_error(
    input_check_library_computation(QC_data = bad_df, downsample_ratio = 0.7, D2_rough = 0.3),
    "num_reads"
  )
})

test_that("input_check_library_computation validates downsample_ratio", {
  qc_data <- data.frame(
    num_reads = 1:10,
    UMI_id = 1:10,
    cell_id = rep("cell1", 10),
    response_id = rep("gene1", 10)
  )

  # NULL downsample_ratio
  expect_error(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = NULL, D2_rough = 0.3),
    "downsample_ratio"
  )

  # Non-numeric downsample_ratio
  expect_error(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = "0.7", D2_rough = 0.3),
    "downsample_ratio"
  )

  # downsample_ratio = 0
  expect_error(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = 0, D2_rough = 0.3),
    "downsample_ratio"
  )

  # downsample_ratio > 1
  expect_error(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = 1.5, D2_rough = 0.3),
    "downsample_ratio"
  )

  # downsample_ratio < 0
  expect_error(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = -0.5, D2_rough = 0.3),
    "downsample_ratio"
  )
})

test_that("input_check_library_computation validates D2_rough", {
  qc_data <- data.frame(
    num_reads = 1:10,
    UMI_id = 1:10,
    cell_id = rep("cell1", 10),
    response_id = rep("gene1", 10)
  )

  # NULL D2_rough
  expect_error(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = 0.7, D2_rough = NULL),
    "D2_rough"
  )

  # Non-numeric D2_rough
  expect_error(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = 0.7, D2_rough = "0.3"),
    "D2_rough"
  )

  # D2_rough < 0
  expect_error(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = 0.7, D2_rough = -0.1),
    "D2_rough"
  )

  # D2_rough > 1
  expect_error(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = 0.7, D2_rough = 1.5),
    "D2_rough"
  )
})

test_that("input_check_library_computation passes with valid inputs", {
  qc_data <- data.frame(
    num_reads = 1:10,
    UMI_id = 1:10,
    cell_id = rep("cell1", 10),
    response_id = rep("gene1", 10)
  )

  # Should pass without error
  expect_silent(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = 0.7, D2_rough = 0.3)
  )

  # Boundary values
  expect_silent(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = 0.001, D2_rough = 0)
  )

  expect_silent(
    input_check_library_computation(QC_data = qc_data, downsample_ratio = 1.0, D2_rough = 1.0)
  )
})

test_that("input_check_library_computation validates empty QC_data", {
  empty_df <- data.frame(
    num_reads = numeric(0),
    UMI_id = integer(0),
    cell_id = character(0),
    response_id = character(0)
  )

  expect_error(
    input_check_library_computation(QC_data = empty_df, downsample_ratio = 0.7, D2_rough = 0.3),
    "cannot be empty"
  )
})

# ============================================================================
# Tests for input_check_power_function
# ============================================================================

test_that("input_check_power_function validates recovery_rate", {
  # NULL recovery_rate
  expect_error(
    input_check_power_function(
      recovery_rate = NULL,
      num_total_reads = 1e9,
      mapping_efficiency = 0.7
    ),
    "recovery_rate"
  )

  # Non-numeric
  expect_error(
    input_check_power_function(
      recovery_rate = "0.5",
      num_total_reads = 1e9,
      mapping_efficiency = 0.7
    ),
    "recovery_rate"
  )

  # Negative value
  expect_error(
    input_check_power_function(
      recovery_rate = -0.1,
      num_total_reads = 1e9,
      mapping_efficiency = 0.7
    ),
    "recovery_rate"
  )

  # Greater than 1
  expect_error(
    input_check_power_function(
      recovery_rate = 1.5,
      num_total_reads = 1e9,
      mapping_efficiency = 0.7
    ),
    "recovery_rate"
  )
})

test_that("input_check_power_function validates num_total_reads", {
  # NULL
  expect_error(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = NULL,
      mapping_efficiency = 0.7
    ),
    "num_total_reads"
  )

  # Non-numeric
  expect_error(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = "1e9",
      mapping_efficiency = 0.7
    ),
    "num_total_reads"
  )

  # Negative
  expect_error(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = -1000,
      mapping_efficiency = 0.7
    ),
    "num_total_reads"
  )

  # Zero
  expect_error(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = 0,
      mapping_efficiency = 0.7
    ),
    "num_total_reads"
  )
})

test_that("input_check_power_function validates mapping_efficiency", {
  # NULL
  expect_error(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = 1e9,
      mapping_efficiency = NULL
    ),
    "mapping_efficiency"
  )

  # Non-numeric
  expect_error(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = 1e9,
      mapping_efficiency = "0.7"
    ),
    "mapping_efficiency"
  )

  # Negative
  expect_error(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = 1e9,
      mapping_efficiency = -0.1
    ),
    "mapping_efficiency"
  )

  # Greater than 1
  expect_error(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = 1e9,
      mapping_efficiency = 1.5
    ),
    "mapping_efficiency"
  )

  # Zero
  expect_error(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = 1e9,
      mapping_efficiency = 0
    ),
    "mapping_efficiency"
  )
})

test_that("input_check_power_function passes with valid inputs", {
  # Create minimal valid data frames
  cells_df <- data.frame(grna_target = "gene1", num_cells = 100)
  baseline_df <- data.frame(response_id = "gene1", relative_expression = 1e-5, expression_size = 1.0)
  discovery_df <- data.frame(grna_target = "gene1", response_id = "gene1")

  # Should pass without error with all required parameters
  expect_silent(
    input_check_power_function(
      recovery_rate = 0.5,
      num_total_reads = 1e9,
      mapping_efficiency = 0.7,
      cells_per_grna = cells_df,
      baseline_relative_expression_stats = baseline_df,
      num_planned_cells = 10000,
      UMI_per_cell = 15000,
      variation = 0.3,
      fold_change_mean = 0.5,
      fold_change_sd = 0.1,
      control_group = "complement",
      side = "left",
      multiple_testing_method = "BH",
      multiple_testing_alpha = 0.05,
      cutoff = NULL,
      discovery_pairs = discovery_df,
      n_nonzero_trt_thresh = 7,
      n_nonzero_cntrl_thresh = 7
    )
  )

  # Boundary values
  expect_silent(
    input_check_power_function(
      recovery_rate = 0.001,
      num_total_reads = 1,
      mapping_efficiency = 0.001,
      cells_per_grna = cells_df,
      baseline_relative_expression_stats = baseline_df,
      num_planned_cells = 100,
      UMI_per_cell = 1000,
      variation = 0,
      fold_change_mean = 0.1,
      fold_change_sd = 0.01,
      control_group = "nt_cells",
      side = "right",
      multiple_testing_method = "BH",
      multiple_testing_alpha = 0.1,
      cutoff = NULL,
      discovery_pairs = discovery_df,
      n_nonzero_trt_thresh = 1,
      n_nonzero_cntrl_thresh = 1
    )
  )
})

# ============================================================================
# Additional edge case tests for input_check_compute_power_plan
# ============================================================================

test_that("input_check_compute_power_plan validates grid_size", {
  baseline <- data.frame(
    response_id = "ENSG001",
    relative_expression = 1e-5,
    expression_size = 1.0
  )
  library_params <- list(UMI_per_cell = 10000, variation = 0.3)

  # grid_size = 0 should fail
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 1.0,
      minimum_fold_change = 0.8,
      cells_per_target = 100,
      sequenced_reads_per_cell = 20000,
      baseline_expression_stats = baseline,
      library_parameters = library_params,
      grid_size = 0
    ),
    "grid_size"
  )

  # Negative grid_size should fail
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 1.0,
      minimum_fold_change = 0.8,
      cells_per_target = 100,
      sequenced_reads_per_cell = 20000,
      baseline_expression_stats = baseline,
      library_parameters = library_params,
      grid_size = -1
    ),
    "grid_size"
  )

  # Non-integer grid_size should fail
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 1.0,
      minimum_fold_change = 0.8,
      cells_per_target = 100,
      sequenced_reads_per_cell = 20000,
      baseline_expression_stats = baseline,
      library_parameters = library_params,
      grid_size = 5.5
    ),
    "grid_size"
  )
})

test_that("input_check_compute_power_plan validates power thresholds", {
  baseline <- data.frame(
    response_id = "ENSG001",
    relative_expression = 1e-5,
    expression_size = 1.0
  )
  library_params <- list(UMI_per_cell = 10000, variation = 0.3)

  # min_power_threshold >= max_power_threshold
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 1.0,
      minimum_fold_change = 0.8,
      cells_per_target = 100,
      sequenced_reads_per_cell = 20000,
      baseline_expression_stats = baseline,
      library_parameters = library_params,
      min_power_threshold = 0.8,
      max_power_threshold = 0.5
    ),
    "min_power_threshold.*max_power_threshold"
  )

  # min_power_threshold = max_power_threshold
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 1.0,
      minimum_fold_change = 0.8,
      cells_per_target = 100,
      sequenced_reads_per_cell = 20000,
      baseline_expression_stats = baseline,
      library_parameters = library_params,
      min_power_threshold = 0.5,
      max_power_threshold = 0.5
    ),
    "min_power_threshold.*max_power_threshold"
  )
})

test_that("input_check_compute_power_plan validates mapping_efficiency boundaries", {
  baseline <- data.frame(
    response_id = "ENSG001",
    relative_expression = 1e-5,
    expression_size = 1.0
  )
  library_params <- list(UMI_per_cell = 10000, variation = 0.3)

  # mapping_efficiency > 1
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 1.0,
      minimum_fold_change = 0.8,
      cells_per_target = 100,
      sequenced_reads_per_cell = 20000,
      baseline_expression_stats = baseline,
      library_parameters = library_params,
      mapping_efficiency = 1.5
    ),
    "mapping_efficiency"
  )

  # mapping_efficiency = 0
  expect_error(
    input_check_compute_power_plan(
      TPM_threshold = 1.0,
      minimum_fold_change = 0.8,
      cells_per_target = 100,
      sequenced_reads_per_cell = 20000,
      baseline_expression_stats = baseline,
      library_parameters = library_params,
      mapping_efficiency = 0
    ),
    "mapping_efficiency"
  )
})
