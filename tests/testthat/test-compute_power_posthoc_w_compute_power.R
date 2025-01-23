library(testthat)

# This is a Rscript testing if the compte_power_posthoc has the same results with compute_power
test_that("Results from compute_power_posthoc aligns with those from compute_power!", {


  # significance level
  alpha <- 0.01

  # return the intermediate outcome
  compute_power_results <- compute_power(
    control_cell_vec = control_cell_vec,
    target_cell_df = target_cell_df,
    baseline_expression = baseline_expression,
    size_parameter = size_parameter,
    effect_size_mean = effect_size_mean,
    effect_size_sd = effect_size_sd,
    side = "left",
    multiple_testing_method = "BH",
    multiple_testing_alpha = alpha,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7,
    intermediate_outcome = FALSE
  )

  compute_power_posthoc_results <- compute_power_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change_mean = effect_size_mean,
    fold_change_sd = effect_size_sd,
    num_total_cells = num_total_cell,
    multiple_testing_method = "BH",
    multiple_testing_alpha = alpha,
    side = "left"
  )

  # join two dfs and compute the absolute difference
  abs_diff <- compute_power_posthoc_results$individual_power |>
    dplyr::mutate(pair = paste(grna_target, response_id, sep = "_")) |>
    dplyr::left_join(compute_power_results$individual_power, "pair") |>
    dplyr::mutate(abs_diff = abs(power - adjusted_power)) |>
    dplyr::select(abs_diff) |> dplyr::pull() |> unique()
  expect_lt(max(abs_diff),
            1e-5)
})


# test when there is no signal
test_that("compute_power_posthoc works when there is no signal!", {

  # specify cutoff
  cutoff <- 0.1

  compute_power_posthoc_results <- compute_power_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change_mean = 1,
    fold_change_sd = 0,
    num_total_cells = num_total_cell,
    n_nonzero_trt_thresh = 0,
    n_nonzero_cntrl_thresh = 0,
    side = "left",
    cutoff = cutoff
  )

  # join two dfs and compute the absolute difference
  abs_diff <- compute_power_posthoc_results$individual_power |>
    dplyr::mutate(pair = paste(grna_target, response_id, sep = "_")) |>
    dplyr::mutate(abs_diff = abs(power - cutoff)) |>
    dplyr::select(abs_diff) |> dplyr::pull() |> unique()
  expect_lt(max(abs_diff),
            1e-5)
})

# test if the mean and sd vectors computed from compute_power and compute_power_posthoc align
test_that("Results from compute_power_posthoc aligns with those from compute_power!", {


  # significance level
  alpha <- 0.1

  # return the intermediate outcome
  compute_power_results <- compute_power(
    control_cell_vec = control_cell_vec,
    target_cell_df = target_cell_df,
    baseline_expression = baseline_expression,
    size_parameter = size_parameter,
    effect_size_mean = effect_size_mean,
    effect_size_sd = effect_size_sd,
    intermediate_outcome = TRUE
  )

  ############################# rename variables ###############################
  control_group = "complement"
  fold_change_mean = effect_size_mean
  fold_change_sd = effect_size_sd
  num_total_cells = num_total_cell
  n_nonzero_trt_thresh = 0
  n_nonzero_cntrl_thresh = 0
  side = "left"
  ############################# create enhancer_gene df ########################
  enhancer_gene <- discovery_pairs |>
    # join grna df
    dplyr::left_join(
      cells_per_grna |>
        dplyr::filter(grna_target != "non-targeting") |>
        dplyr::group_by(grna_target) |>
        dplyr::summarize(num_trt_cells = sum(num_cells),
                         num_trt_cells_sq = sum(num_cells^2)) |>
        dplyr::ungroup(),
      "grna_target"
    ) |>
    # join gene expression df
    dplyr::left_join(
      baseline_expression_stats, "response_id"
    )

  ############################# obtain number of control cells #################
  if (control_group == "nt_cells") {
    num_cntrl_cells <- cells_per_grna |>
      dplyr::filter(grna_target == "non-targeting") |>
      dplyr::summarize(sum(num_cells)) |>
      dplyr::pull()
    enhancer_gene <- enhancer_gene |> dplyr::mutate(
      num_cntrl_cells = num_cntrl_cells
    )
  } else { # control_group == "complement"
    enhancer_gene <- enhancer_gene |> dplyr::mutate(
      num_cntrl_cells = num_total_cells - num_trt_cells
    )
  }

  ################## transform the scalar-valued effect size mean/sd ###########
  if(is.numeric(fold_change_mean)){
    # create the effect size matrices
    enhancer_gene <- enhancer_gene |> dplyr::mutate(
      fold_change_mean = fold_change_mean,
      fold_change_sd = fold_change_sd
    )
  }else{
    # join the enhancer_gene df and effect size mean(sd) data frames
    enhancer_gene <- enhancer_gene |>
      dplyr::left_join(
        fold_change_mean,
        c("grna_target", "response_id")
      ) |>
      dplyr::left_join(
        fold_change_sd,
        c("grna_target", "response_id")
      )
  }

  ########################### prepare for multiple testing #####################
  enhancer_gene <- enhancer_gene |>
    dplyr::group_by(grna_target, response_id) |>
    dplyr::mutate(
      # compute mean and sd of the test statistic for each pair
      test_stat_distribution = compute_distribution_teststat(
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        num_trt_cells_sq = num_trt_cells_sq,
        expression_mean = expression_mean,
        expression_size = expression_size,
        fold_change_mean = fold_change_mean,
        fold_change_sd = fold_change_sd
      ),
      # extract mean and sd from test_stat_distribution
      mean_test_stat = unlist(test_stat_distribution)["mean"],
      sd_test_stat = unlist(test_stat_distribution)["sd"],
      # compute QC probability
      QC_prob = compute_QC(
        fold_change_mean = fold_change_mean,
        expression_mean = expression_mean,
        expression_size = expression_size,
        num_cntrl_cells = num_cntrl_cells,
        num_trt_cells = num_trt_cells,
        n_nonzero_trt_thresh = n_nonzero_trt_thresh,
        n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh)
    ) |>
    dplyr::select(-test_stat_distribution) |>
    dplyr::ungroup()

  # join two dfs and compute the absolute difference
  abs_diff <- enhancer_gene |>
    dplyr::mutate(pair = paste(grna_target, response_id, sep = "_")) |>
    dplyr::left_join(compute_power_results, "pair") |>
    dplyr::mutate(abs_diff_sd = abs(sd_test_stat - asy_sd)) |>
    dplyr::select(abs_diff_sd) |> dplyr::pull()

  # check if the difference is small enough
  expect_lt(max(abs_diff), 1e-4)
})


# test when there is no signal
test_that("compute_power_posthoc works when there is no signal!", {

  # specify cutoff
  cutoff <- 0.1

  compute_power_posthoc_results <- compute_power_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change_mean = 1,
    fold_change_sd = 0,
    num_total_cells = num_total_cell,
    n_nonzero_trt_thresh = 0,
    n_nonzero_cntrl_thresh = 0,
    side = "left",
    cutoff = cutoff
  )

  # join two dfs and compute the absolute difference
  abs_diff <- compute_power_posthoc_results$individual_power |>
    dplyr::mutate(pair = paste(grna_target, response_id, sep = "_")) |>
    dplyr::mutate(abs_diff = abs(power - cutoff)) |>
    dplyr::select(abs_diff) |> dplyr::pull() |> unique()
  expect_lt(max(abs_diff),
            1e-5)
})
