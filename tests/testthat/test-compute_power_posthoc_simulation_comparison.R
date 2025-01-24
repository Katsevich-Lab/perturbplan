library(dplyr)
library(testthat)
library(tidyr)

# load simulation results
file_path_cutoff <- system.file("extdata", "simulation_power_cutoff.rds", package = "perturbplan")
file_path_BH <- system.file("extdata", "simulation_power_BH.rds", package = "perturbplan")
simulation_power_cutoff <- readRDS(file_path_cutoff)
simulation_power_BH <- readRDS(file_path_BH)

# test when there is no signal
test_that("compute_power_posthoc aligns well with the simulation results (no cutoff)!", {

  # specify alpha
  alpha <- unique(simulation_power_BH$alpha)

  # compute analytic power
  compute_power_posthoc_results <- compute_power_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change_mean = effect_size_mean,
    fold_change_sd = effect_size_sd,
    num_total_cells = num_total_cell,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7,
    side = "left",
    multiple_testing_alpha = alpha
  )

  # merge two dfs
  merged_df <- compute_power_posthoc_results$individual_power |>
    dplyr::mutate(pair = paste(grna_target, response_id, sep = "_")) |>
    dplyr::left_join(simulation_power_BH, "pair")
  # compute the absolute difference
  abs_diff <- merged_df |>
    dplyr::mutate(abs_diff = abs(power - sceptre_power)) |>
    dplyr::select(abs_diff) |> dplyr::pull() |> unique()
  expect_lt(max(abs_diff),
            0.12)
})

# # plot scatter plot
# merged_df |>
#   ggplot2::ggplot(ggplot2::aes(x = sceptre_power, y = power)) +
#   ggplot2::geom_point() + ggplot2::geom_abline() +
#   ggplot2::theme_bw() +
#   ggplot2::labs(title = "With QC")


# test when there is no signal
test_that("compute_power_posthoc aligns well with the simulation results (with cutoff)!", {

  # specify alpha
  cutoff <- unique(simulation_power_cutoff$cutoff)

  # compute analytic power
  compute_power_posthoc_results <- compute_power_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change_mean = effect_size_mean,
    fold_change_sd = effect_size_sd,
    num_total_cells = num_total_cell,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7,
    side = "left",
    cutoff = cutoff
  )

  # merge two dfs
  merged_df <- compute_power_posthoc_results$individual_power |>
    dplyr::mutate(pair = paste(grna_target, response_id, sep = "_")) |>
    dplyr::left_join(simulation_power_cutoff, "pair")
  # compute the absolute difference
  abs_diff <- merged_df |>
    dplyr::mutate(abs_diff = abs(power - sceptre_power)) |>
    dplyr::select(abs_diff) |> dplyr::pull() |> unique()
  expect_lt(max(abs_diff),
            0.09)
})

