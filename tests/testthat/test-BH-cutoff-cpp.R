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

  # compute analytic power with BH cutoff search in R
  compute_power_R_BH <- compute_power_posthoc_fixed_es(
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

  # compute analytic power with BH cutoff search in R
  compute_power_Cpp_BH <- compute_power_posthoc_fixed_es(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change = effect_size_mean,
    num_total_cells = num_total_cell,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7,
    side = "left",
    multiple_testing_method = "BH_cpp",
    multiple_testing_alpha = alpha,
    random_assignment = FALSE
  )

  # compute the difference
  discovery_diff <- compute_power_Cpp_BH$expected_num_discoveries - compute_power_R_BH$expected_num_discoveries
  expect_equal(discovery_diff, 0)
})
