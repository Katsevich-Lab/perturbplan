library(dplyr)
library(testthat)

# Load test data
source("helper-compute_power_posthoc.R")

test_that("compute_power_posthoc_fixed_fc produces identical results to compute_power_posthoc_cpp with random_assignment=FALSE", {
  
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
  
  result_cpp_fixed <- compute_power_posthoc_cpp(
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
    multiple_testing_alpha = alpha,
    random_assignment = FALSE
  )
  
  # Test that results are identical
  expect_identical(result_fixed_fc, result_cpp_fixed)
  
  # Test individual components
  expect_equal(result_fixed_fc$expected_num_discoveries, result_cpp_fixed$expected_num_discoveries)
  expect_identical(result_fixed_fc$individual_power, result_cpp_fixed$individual_power)
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
  
  result_cpp_right <- compute_power_posthoc_cpp(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    side = "right",
    multiple_testing_alpha = alpha,
    random_assignment = FALSE
  )
  
  expect_identical(result_fixed_fc_right, result_cpp_right)
  
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
  
  result_cpp_both <- compute_power_posthoc_cpp(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    side = "both",
    multiple_testing_alpha = alpha,
    random_assignment = FALSE
  )
  
  expect_identical(result_fixed_fc_both, result_cpp_both)
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
  
  result_cpp_bonf <- compute_power_posthoc_cpp(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    multiple_testing_method = "bonferroni",
    multiple_testing_alpha = alpha,
    random_assignment = FALSE
  )
  
  expect_identical(result_fixed_fc_bonf, result_cpp_bonf)
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
  
  result_cpp_cutoff <- compute_power_posthoc_cpp(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    cutoff = cutoff_value,
    random_assignment = FALSE
  )
  
  expect_identical(result_fixed_fc_cutoff, result_cpp_cutoff)
})