# This is a Rscript computing the power function using score test

#' Compute power for each perturbation-gene pair
#'
#' @param discovery_pairs A data frame specifying which element-gene pairs to consider, with columns `grna_target` and `response_id`
#' @param cells_per_grna A data frame specifying how many cells contain each gRNA, with columns `grna_id`, `grna_target`, and `num_cells`
#' @param baseline_expression_stats A data frame specifying the baseline expression statistics for each gene, with columns `response_id`, `expression_mean`, and `expression_size`
#' @param control_group A character string specifying the control group, either "complement" or "nt_cells"
#' @param fold_change_mean A numeric value to use for mean effect size for all element-gene pairs
#' @param fold_change_sd A numeric value to use for standard deviation of effect size for all element-gene pairs
#' @param num_total_cells (Required only if control_group == "complement") A positive integer specifying the total number of cells in the experiment
#' @param cutoff (Optional) A numeric value between 0 and 1 to use as the p-value cutoff
#' @param n_nonzero_trt_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param n_nonzero_cntrl_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param side (Optional) A character string specifying the side of the test, either "left", "right", or "both"; defaults to "both"
#' @param multiple_testing_method (Optional) A character string specifying the multiple testing correction method to use, either "BH" or "bonferroni"; defaults to "BH"
#' @param multiple_testing_alpha (Optional) A numeric value between 0 and 1 specifying the alpha level for multiple testing correction; defaults to 0.1
#'
#' @return A list with two elements: `individual_power` (a data frame with columns `grna_target`, `response_id`, and `power`) and `expected_num_discoveries` (a numeric value)
#' @export
compute_power_posthoc <- function(
    discovery_pairs,
    cells_per_grna,
    baseline_expression_stats,
    control_group,
    fold_change_mean,
    fold_change_sd,
    num_total_cells = NULL,
    cutoff = NULL,
    n_nonzero_trt_thresh = 7L,
    n_nonzero_cntrl_thresh = 7L,
    side = "both",
    multiple_testing_method = "BH",
    multiple_testing_alpha = 0.1) {

  ############################# perform input checks ###########################
  input_check_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = control_group,
    fold_change_mean = fold_change_mean,
    fold_change_sd = fold_change_sd,
    num_total_cells = num_total_cells,
    cutoff = cutoff,
    n_nonzero_trt_thresh = n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
    side = side,
    multiple_testing_method = multiple_testing_method,
    multiple_testing_alpha = multiple_testing_alpha
  )

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

  ########################### correct multiplicity #############################
  # compute cutoff if it is NULL
  if(is.null(cutoff)){
    cutoff <- enhancer_gene |>
      dplyr::summarize(
        cutoff = adjusted_cutoff(mean_list = mean_test_stat,
                                 sd_list = sd_test_stat,
                                 multiple_testing_alpha = multiple_testing_alpha,
                                 multiple_testing_method = multiple_testing_method,
                                 side = side, QC_prob = QC_prob)
      ) |> dplyr::select(cutoff) |> dplyr::pull()
  }
  # compute the adjusted power
  enhancer_gene <- enhancer_gene |>
    dplyr::mutate(
      cutoff = cutoff,
      power = rejection_computation(mean_list = mean_test_stat,
                                    sd_list = sd_test_stat,
                                    side = side,
                                    cutoff = cutoff) * (1 - QC_prob)
    )

  # store individual power and rejection size as the output
  output <- list(
    individual_power = enhancer_gene |> dplyr::select(grna_target, response_id, power),
    expected_num_discoveries = sum(enhancer_gene$power)
  )

  return(output)
}



power_function <- function(

    ###################### specify experimental design parameters ##############
    recovery_rate = NULL,                                              # Library prep parameters
    num_total_reads = NULL, mapping_efficiency = NULL,                 # Sequencing parameters

    ######################## specify the power-determining parameters ##########
    cell_per_grna, baseline_relative_expression_stats,
    fold_change_mean, fold_change_sd, num_planned_cells, control_group,

    ###################### specify test-related parameters #####################
    side = "both", multiple_testing_method = "BH", multiple_testing_alpha = 0.1,
    cutoff = NULL, discovery_pairs,

    ######################## specify QC-related parameters ################
    n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L
){

  ########## compute library size with other power-determing parameters ########
  # compute the number of total cells (singletons) surviving from library preparation
  num_total_cells <- num_planned_cells * recovery_rate

  # compute the reads per cell
  reads_per_cell <- read_per_cell_after_QC(num_total_reads = num_total_reads,
                                           mapping_efficiency = mapping_efficiency,
                                           num_total_cells = num_total_cells)

  # compute the averaged library size with read per cell
  avg_library_size <- library_computation(reads_per_cell = reads_per_cell,
                                          h5_path = h5_path)

  ####### perform power calculation with power-determining parameters ##########
  # compute the baseline_expression
  baseline_expression_stats <- baseline_relative_expression_stats |>
    dplyr::mutate(expression_mean = avg_library_size * relative_expression) |>
    dplyr::select(-relative_expression)

  # compute power using function compute_power
  power_result <- compute_power_posthoc()

  # return the power_result
  return(power_result)
}
