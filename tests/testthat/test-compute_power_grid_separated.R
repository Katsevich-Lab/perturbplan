test_that("compute_power_grid_separated works correctly", {

  # Create a small test dataset
  cells_reads_df <- data.frame(
    num_total_cells = c(50000, 100000),
    reads_per_cell = c(1000, 2000)
  )

  # Set common parameters for both functions
  test_params <- list(
    num_targets = 10,
    gRNAs_per_target = 2,
    non_targeting_gRNAs = 5,
    num_pairs = 100,
    tpm_threshold = 10,
    fdr_target = 0.05,
    fc_mean = 0.85,
    fc_sd = 0.15,
    prop_non_null = 0.1,
    MOI = 1,  # Reduced MOI to keep treatment cells reasonable
    biological_system = "K562",
    experimental_platform = "10x Chromium v3",
    side = "left",
    control_group = "complement",
    B = 50  # Small B for faster testing
  )

  # Set seed for reproducibility
  set.seed(123)

  # Run separated approach
  result_separated <- do.call(compute_power_grid_separated, c(list(cells_reads_df = cells_reads_df), test_params))

  # Check that result has the expected structure
  expect_true(is.data.frame(result_separated))
  expect_true("overall_power" %in% names(result_separated))
  expect_true("num_total_cells" %in% names(result_separated))
  expect_true("reads_per_cell" %in% names(result_separated))

  # Check that power values are reasonable (between 0 and 1)
  expect_true(all(result_separated$overall_power >= 0))
  expect_true(all(result_separated$overall_power <= 1))

  # Check that input values are preserved
  expect_identical(result_separated$num_total_cells, cells_reads_df$num_total_cells)
  expect_identical(result_separated$reads_per_cell, cells_reads_df$reads_per_cell)

  # Check that power curves are included
  expect_true("power_by_fc" %in% names(result_separated))
  expect_true("power_by_expr" %in% names(result_separated))
  
  # Check that power curves have the expected structure
  expect_true(is.list(result_separated$power_by_fc))
  expect_true(is.list(result_separated$power_by_expr))
})
