library(dplyr)
library(testthat)

# Load test data
source("helper-compute_power_posthoc.R")

test_that("compute_power_posthoc_fixed_fc works correctly with BH method", {
  
  # Test parameters
  alpha <- 0.1
  
  # Test with BH method
  result_fixed_fc <- compute_power_posthoc_fixed_fc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7,
    side = "left",
    multiple_testing_method = "BH",
    multiple_testing_alpha = alpha
  )
  
  # Test that result structure is correct
  expect_true(is.list(result_fixed_fc))
  expect_true("individual_power" %in% names(result_fixed_fc))
  expect_true("expected_num_discoveries" %in% names(result_fixed_fc))
  expect_true(is.numeric(result_fixed_fc$expected_num_discoveries))
  expect_true(is.data.frame(result_fixed_fc$individual_power))
  
  # Test that power values are reasonable (between 0 and 1)
  expect_true(all(result_fixed_fc$individual_power$power >= 0))
  expect_true(all(result_fixed_fc$individual_power$power <= 1))
  
  # Test that expected discoveries is reasonable
  expect_true(result_fixed_fc$expected_num_discoveries >= 0)
  expect_true(result_fixed_fc$expected_num_discoveries <= nrow(discovery_pairs))
})

test_that("compute_power_posthoc_fixed_fc works with different test sides", {
  
  alpha <- 0.1
  
  # Test with right-sided test
  result_fixed_fc_right <- compute_power_posthoc_fixed_fc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    side = "right",
    multiple_testing_alpha = alpha
  )
  
  expect_true(is.list(result_fixed_fc_right))
  expect_true(result_fixed_fc_right$expected_num_discoveries >= 0)
  
  # Test with two-sided test
  result_fixed_fc_both <- compute_power_posthoc_fixed_fc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    side = "both",
    multiple_testing_alpha = alpha
  )
  
  expect_true(is.list(result_fixed_fc_both))
  expect_true(result_fixed_fc_both$expected_num_discoveries >= 0)
  
  # Test that different sides give different results (generally)
  result_left <- compute_power_posthoc_fixed_fc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    side = "left",
    multiple_testing_alpha = alpha
  )
  
  # Left and right sided tests should generally give different results
  expect_false(identical(result_left$expected_num_discoveries, result_fixed_fc_right$expected_num_discoveries))
})

test_that("compute_power_posthoc_fixed_fc works with different multiple testing methods", {
  
  alpha <- 0.1
  
  # Test with bonferroni method
  result_fixed_fc_bonf <- compute_power_posthoc_fixed_fc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    multiple_testing_method = "bonferroni",
    multiple_testing_alpha = alpha
  )
  
  expect_true(is.list(result_fixed_fc_bonf))
  expect_true(result_fixed_fc_bonf$expected_num_discoveries >= 0)
  
  # Test with BH method for comparison
  result_fixed_fc_bh <- compute_power_posthoc_fixed_fc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    multiple_testing_method = "BH",
    multiple_testing_alpha = alpha
  )
  
  # Bonferroni should generally be more conservative (fewer discoveries)
  expect_true(result_fixed_fc_bonf$expected_num_discoveries <= result_fixed_fc_bh$expected_num_discoveries)
})

test_that("compute_power_posthoc_fixed_fc works with nt_cells control group", {
  
  alpha <- 0.1
  
  # Add non-targeting cells to test data
  cells_per_grna_with_nt <- cells_per_grna |>
    dplyr::bind_rows(
      data.frame(
        grna_id = "nt_1",
        grna_target = "non-targeting",
        num_cells = 1000
      )
    )
  
  result_fixed_fc_nt <- compute_power_posthoc_fixed_fc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna_with_nt,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "nt_cells",
    fold_change = effect_size_mean,
    multiple_testing_alpha = alpha
  )
  
  # Note: compute_power_posthoc_cpp has a bug when random_assignment=FALSE with nt_cells
  # It tries to use mean_num_cells on the original cells_per_grna but only renames in grna_gene
  # So we skip the comparison for nt_cells case and just test that our function works
  expect_true(is.list(result_fixed_fc_nt))
  expect_true("individual_power" %in% names(result_fixed_fc_nt))
  expect_true("expected_num_discoveries" %in% names(result_fixed_fc_nt))
  expect_true(is.numeric(result_fixed_fc_nt$expected_num_discoveries))
})

test_that("compute_power_posthoc_fixed_fc works with custom cutoff", {
  
  cutoff_value <- 0.05
  
  result_fixed_fc_cutoff <- compute_power_posthoc_fixed_fc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    cutoff = cutoff_value
  )
  
  expect_true(is.list(result_fixed_fc_cutoff))
  expect_true(result_fixed_fc_cutoff$expected_num_discoveries >= 0)
  
  # Test with default cutoff for comparison
  result_fixed_fc_default <- compute_power_posthoc_fixed_fc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell
  )
  
  # Results should generally be different when using different cutoffs
  expect_false(identical(result_fixed_fc_cutoff$expected_num_discoveries, result_fixed_fc_default$expected_num_discoveries))
})