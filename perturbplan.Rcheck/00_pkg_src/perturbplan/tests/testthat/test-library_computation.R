# This is a Rscript testing the input_check_posthoc function
library(testthat)

# test when there is no signal
test_that("input_check_posthoc works correctly!", {

  # Test valid input (expecting no error)
  expect_silent({
    compute_power_posthoc(
      discovery_pairs = discovery_pairs,
      cells_per_grna = cells_per_grna,
      baseline_expression_stats = baseline_expression_stats,
      control_group = "complement",
      fold_change_mean = effect_size_mean,
      fold_change_sd = effect_size_sd,
      num_total_cells = num_total_cell
    )
  })

  # Test invalid inputs (expecting an error)
  # grna_id is missing
  expect_error({
    compute_power_posthoc(
      discovery_pairs = discovery_pairs |> dplyr::select(grna_target),
      cells_per_grna = cells_per_grna,
      baseline_expression_stats = baseline_expression_stats,
      control_group = "complement",
      fold_change_mean = effect_size_mean,
      fold_change_sd = effect_size_sd,
      num_total_cells = num_total_cell
    )
  })
  # grna_target does not match those in cells_per_grna
  expect_error({
    compute_power_posthoc(
      discovery_pairs = discovery_pairs,
      cells_per_grna = cells_per_grna |> dplyr::filter(grna_target == "enh1"),
      baseline_expression_stats = baseline_expression_stats,
      control_group = "complement",
      fold_change_mean = effect_size_mean,
      fold_change_sd = effect_size_sd,
      num_total_cells = num_total_cell
    )
  })
  # gene set does not match in baseline_expression_stats and discovery_pairs
  expect_error({
    compute_power_posthoc(
      discovery_pairs = discovery_pairs |> dplyr::filter(gene == gene_expression$ensembl[1]),
      cells_per_grna = cells_per_grna,
      baseline_expression_stats = baseline_expression_stats,
      control_group = "complement",
      fold_change_mean = effect_size_mean,
      fold_change_sd = effect_size_sd,
      num_total_cells = num_total_cell
    )
  })
  # data types in fold_change_mean and fold_change_sd do not match
  expect_error({
    compute_power_posthoc(
      discovery_pairs = discovery_pairs |> dplyr::filter(gene == gene_expression$ensembl[1]),
      cells_per_grna = cells_per_grna,
      baseline_expression_stats = baseline_expression_stats,
      control_group = "complement",
      fold_change_mean = as.matrix(effect_size_mean),
      fold_change_sd = effect_size_sd,
      num_total_cells = num_total_cell
    )
  })
  # data structure in fold_change_mean and fold_change_sd do not match
  expect_error({
    compute_power_posthoc(
      discovery_pairs = discovery_pairs,
      cells_per_grna = cells_per_grna,
      baseline_expression_stats = baseline_expression_stats,
      control_group = "complement",
      fold_change_mean = as.data.frame(effect_size_mean),
      fold_change_sd = effect_size_sd,
      num_total_cells = num_total_cell
    )
  })
  # for high-MOI, num_total_cells is not specified
  expect_error({
    compute_power_posthoc(
      discovery_pairs = discovery_pairs,
      cells_per_grna = cells_per_grna,
      baseline_expression_stats = baseline_expression_stats,
      control_group = "complement",
      fold_change_mean = as.data.frame(effect_size_mean),
      fold_change_sd = effect_size_sd
    )
  })
  # for low-MOI, no nt_cells in cells_per_grna
  expect_error({
    compute_power_posthoc(
      discovery_pairs = discovery_pairs,
      cells_per_grna = cells_per_grna,
      baseline_expression_stats = baseline_expression_stats,
      control_group = "nt_cells",
      fold_change_mean = effect_size_mean,
      fold_change_sd = effect_size_sd
    )
  })

})
