library(dplyr)
library(testthat)
library(tidyr)

# load simulation results
file_path_cutoff <- system.file("extdata", "simulation_power_cutoff.rds", package = "perturbplan")
file_path_BH <- system.file("extdata", "simulation_power_BH.rds", package = "perturbplan")
simulation_power_cutoff <- readRDS(file_path_cutoff)
simulation_power_BH <- readRDS(file_path_BH)

# test when there is no signal
test_that("BH Cpp implementatio works!", {

  # specify alpha
  alpha <- unique(simulation_power_BH$alpha)

  # modify the cells_per_grna
  cells_per_grna <- cells_per_grna |>
    dplyr::rename(mean_num_cells = num_cells) |>
    dplyr::mutate(sd_num_cells = 0)

  # compute analytic power with distribution computation in R
  computation_bisection <- system.time({
    compute_power_bisection <- compute_power_posthoc_fixed_es(
      discovery_pairs = discovery_pairs,
      cells_per_grna = cells_per_grna,
      baseline_expression_stats = baseline_expression_stats,
      control_group = "complement",
      fold_change = effect_size_mean,
      num_total_cells = num_total_cell,
      n_nonzero_trt_thresh = 7,
      n_nonzero_cntrl_thresh = 7,
      side = "left",
      multiple_testing_method = "BH_bisection",
      multiple_testing_alpha = alpha,
      random_assignment = TRUE
    )
  })["elapsed"]

  # compute analytic power with distribution computation in Cpp
  computation_cpp <- system.time({
    compute_power_cpp <- compute_power_posthoc_cpp(
      discovery_pairs = discovery_pairs,
      cells_per_grna = cells_per_grna,
      baseline_expression_stats = baseline_expression_stats,
      control_group = "complement",
      fold_change = effect_size_mean,
      num_total_cells = num_total_cell,
      n_nonzero_trt_thresh = 7,
      n_nonzero_cntrl_thresh = 7,
      side = "left",
      multiple_testing_method = "BH_bisection",
      multiple_testing_alpha = alpha,
      random_assignment = TRUE
    )
  })["elapsed"]

  # compute the difference
  num_total_discovery_cpp <- compute_power_cpp$expected_num_discoveries
  num_total_discovery_bisection <- compute_power_bisection$expected_num_discoveries
  discovery_diff <- num_total_discovery_cpp - num_total_discovery_bisection
  expect_lt(abs(discovery_diff), 1e-10)
})
