test_that("compute_power_grid_separated produces identical output to compute_power_grid", {

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

  # Run original function
  result_original <- do.call(compute_power_grid, c(list(cells_reads_df = cells_reads_df), test_params))

  # Set same seed for C++ version
  set.seed(123)

  # Run C++ accelerated version
  result_separated <- do.call(compute_power_grid_separated, c(list(cells_reads_df = cells_reads_df), test_params))

  # Check that both results have the same structure
  expect_identical(names(result_original), names(result_separated))
  expect_identical(nrow(result_original), nrow(result_separated))
  expect_identical(ncol(result_original), ncol(result_separated))

  # Check that overall power values are identical (or very close)
  expect_equal(result_original$overall_power, result_separated$overall_power, tolerance = 1e-6)

  # Check that num_total_cells and reads_per_cell are identical
  expect_identical(result_original$num_total_cells, result_separated$num_total_cells)
  expect_identical(result_original$reads_per_cell, result_separated$reads_per_cell)

  # Test that power functions return identical values for the same inputs
  for (i in 1:nrow(result_original)) {
    # Test power_by_fc function
    power_fc_separated <- result_separated$power_by_fc[[i]]
    power_fc_orig <- sapply(power_fc_separated$fold_change, function(fold_change){
      result_original$power_by_fc[[i]](fold_change)}
      )
    expect_equal(power_fc_orig, power_fc_separated$power, tolerance = 1e-6)

    # Test power_by_pi function
    power_expr_separated <- result_separated$power_by_expr[[i]]
    power_expr_orig <- sapply(power_expr_separated$relative_expression, function(pi){
      result_original$power_by_expr[[i]](TPM = pi * 1e6)}
    )
    expect_equal(power_expr_orig, power_expr_separated$power, tolerance = 1e-6)
  }
})
