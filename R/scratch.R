# This is a Rscript computing the power function using score test

#' Compute power for each perturbation-gene pair
#'
#' @param discovery_pairs A data frame specifying which element-gene pairs to consider, with columns `grna_target` and `response_id`; it can also have `grna_id` as a column but this is optional
#' @param cells_per_grna A data frame specifying how many cells contain each gRNA, with columns `grna_id`, `grna_target`, and `num_cells`
#' @param baseline_expression_stats A data frame specifying the baseline expression statistics for each gene, with columns `response_id`, `expression_mean`, and `expression_size`
#' @param control_group A character string specifying the control group, either "complement" or "nt_cells"
#' @param fold_change A data frame to use for fixed effect size for all gRNA-gene pairs
#' @param num_total_cells (Required only if control_group == "complement") A positive integer specifying the total number of cells in the experiment
#' @param cutoff (Optional) A numeric value between 0 and 1 to use as the p-value cutoff
#' @param n_nonzero_trt_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param n_nonzero_cntrl_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param side (Optional) A character string specifying the side of the test, either "left", "right", or "both"; defaults to "both"
#' @param multiple_testing_method (Optional) A character string specifying the multiple testing correction method to use, either "BH" or "bonferroni"; defaults to "BH"
#' @param multiple_testing_alpha (Optional) A numeric value between 0 and 1 specifying the alpha level for multiple testing correction; defaults to 0.1
#' @param random_assignment (Optional) A logical value specifying if random gRNA assignment is used or not
#'
#' @return A list with two elements: `individual_power` (a data frame with columns `grna_target`, `response_id`, and `power`) and `expected_num_discoveries` (a numeric value)
#' @export
compute_power_posthoc_fixed_es <- function(
    discovery_pairs,
    cells_per_grna,
    baseline_expression_stats,
    control_group,
    fold_change,
    num_total_cells = NULL,
    cutoff = NULL,
    n_nonzero_trt_thresh = 7L,
    n_nonzero_cntrl_thresh = 7L,
    side = "both",
    multiple_testing_method = "BH",
    multiple_testing_alpha = 0.1,
    random_assignment = FALSE) {

  # ############################# perform input checks ###########################
  # input_check_posthoc_fixed_es(
  #   discovery_pairs = discovery_pairs,
  #   cells_per_grna = cells_per_grna,
  #   baseline_expression_stats = baseline_expression_stats,
  #   control_group = control_group,
  #   num_total_cells = num_total_cells,
  #   cutoff = cutoff,
  #   n_nonzero_trt_thresh = n_nonzero_trt_thresh,
  #   n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
  #   side = side,
  #   multiple_testing_method = multiple_testing_method,
  #   multiple_testing_alpha = multiple_testing_alpha
  # )

  ############################# create grna_gene df ############################
  grna_gene <-  cells_per_grna |>
    # exclude the non-targeting gRNAs
    dplyr::filter(grna_target != "non-targeting") |>
    # join discovery pairs
    dplyr::left_join(discovery_pairs,
                     unlist(ifelse("grna_id" %in% colnames(discovery_pairs), list(c("grna_id", "grna_target")), "grna_target"))) |>
    # join gene expression df
    dplyr::left_join(baseline_expression_stats, "response_id")

  ################### obtain number of treatment and control cells #############
  # compute the number of treatment cells by grouping grna_target and response_id
  if(!random_assignment){
    grna_gene <- grna_gene |>
      dplyr::rename(mean_num_cells = num_cells) |>
      dplyr::group_by(grna_target, response_id) |>
      dplyr::mutate(num_trt_cells = sum(mean_num_cells),
                    sd_num_cells = 0) |>
      dplyr::ungroup()
  }else{
    grna_gene <- grna_gene |>
      dplyr::group_by(grna_target, response_id) |>
      dplyr::mutate(num_trt_cells = sum(mean_num_cells)) |>
      dplyr::ungroup()
  }

  # define the control cells based on control_group
  if (control_group == "nt_cells") {
    # compute the number of control cells using cells receiving non-targeting gRNAs
    num_cntrl_cells <- cells_per_grna |>
      dplyr::filter(grna_target == "non-targeting") |>
      dplyr::summarize(sum(mean_num_cells)) |>
      dplyr::pull()

    grna_gene <- grna_gene |> dplyr::mutate(num_cntrl_cells = num_cntrl_cells)
  } else { # control_group == "complement"
    grna_gene <- grna_gene |> dplyr::mutate(num_cntrl_cells = num_total_cells - num_trt_cells)
  }

  ################## transform the scalar-valued effect size mean/sd ###########
  if(is.numeric(fold_change)){
    # append the scalar to form a new column in grna_gene
    grna_gene <- grna_gene |> dplyr::mutate(fold_change = fold_change)
  }else{
    # join the grna_gene df and effect size data frames
    grna_gene <- grna_gene |> dplyr::left_join(fold_change, c("grna_id", "response_id"))
  }

  ########################### prepare for multiple testing #####################
  enhancer_gene <- grna_gene |>
    dplyr::group_by(grna_target, response_id) |>
    dplyr::summarise(
      # compute mean and sd of the test statistic for each pair
      test_stat_distribution = list({
        if(random_assignment){
          compute_distribution_teststat_fixed_es_random_assignment(
            mean_num_cells = mean_num_cells,
            sd_num_cells = sd_num_cells,
            num_cntrl_cells = num_cntrl_cells,
            expression_mean = expression_mean,
            expression_size = expression_size,
            fold_change = fold_change
          )
        }else{
          compute_distribution_teststat_fixed_es(
            num_trt_cells = num_trt_cells,
            num_cntrl_cells = num_cntrl_cells,
            num_cells = mean_num_cells,
            expression_mean = expression_mean,
            expression_size = expression_size,
            fold_change = fold_change
          )
        }
      }),

      # extract mean and sd from test_stat_distribution
      mean_test_stat = unlist(test_stat_distribution)["mean"],
      sd_test_stat = unlist(test_stat_distribution)["sd"],

      # compute QC probability
      QC_prob = compute_QC_fixed_es(
        fold_change = fold_change,
        expression_mean = expression_mean,
        expression_size = expression_size,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = mean_num_cells,
        n_nonzero_trt_thresh = n_nonzero_trt_thresh,
        n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-test_stat_distribution)

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

# This is a Rscript computing the power function using score test

#' Compute power for each perturbation-gene pair
#'
#' @param discovery_pairs A data frame specifying which element-gene pairs to consider, with columns `grna_target` and `response_id`; it can also have `grna_id` as a column but this is optional
#' @param cells_per_grna A data frame specifying how many cells contain each gRNA, with columns `grna_id`, `grna_target`, and `num_cells`
#' @param baseline_expression_stats A data frame specifying the baseline expression statistics for each gene, with columns `response_id`, `expression_mean`, and `expression_size`
#' @param control_group A character string specifying the control group, either "complement" or "nt_cells"
#' @param fold_change A data frame to use for fixed effect size for all gRNA-gene pairs
#' @param num_total_cells (Required only if control_group == "complement") A positive integer specifying the total number of cells in the experiment
#' @param cutoff (Optional) A numeric value between 0 and 1 to use as the p-value cutoff
#' @param n_nonzero_trt_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param n_nonzero_cntrl_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param side (Optional) A character string specifying the side of the test, either "left", "right", or "both"; defaults to "both"
#' @param multiple_testing_method (Optional) A character string specifying the multiple testing correction method to use, either "BH" or "bonferroni"; defaults to "BH"
#' @param multiple_testing_alpha (Optional) A numeric value between 0 and 1 specifying the alpha level for multiple testing correction; defaults to 0.1
#' @param random_assignment (Optional) A logical value specifying if random gRNA assignment is used or not
#'
#' @return A list with two elements: `individual_power` (a data frame with columns `grna_target`, `response_id`, and `power`) and `expected_num_discoveries` (a numeric value)
#' @export
compute_power_posthoc_cpp <- function(
    discovery_pairs,
    cells_per_grna,
    baseline_expression_stats,
    control_group,
    fold_change,
    num_total_cells = NULL,
    cutoff = NULL,
    n_nonzero_trt_thresh = 7L,
    n_nonzero_cntrl_thresh = 7L,
    side = "both",
    multiple_testing_method = "BH",
    multiple_testing_alpha = 0.1,
    random_assignment = FALSE) {

  # ############################# perform input checks ###########################
  # input_check_posthoc_fixed_es(
  #   discovery_pairs = discovery_pairs,
  #   cells_per_grna = cells_per_grna,
  #   baseline_expression_stats = baseline_expression_stats,
  #   control_group = control_group,
  #   num_total_cells = num_total_cells,
  #   cutoff = cutoff,
  #   n_nonzero_trt_thresh = n_nonzero_trt_thresh,
  #   n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
  #   side = side,
  #   multiple_testing_method = multiple_testing_method,
  #   multiple_testing_alpha = multiple_testing_alpha
  # )

  ############################# create grna_gene df ############################
  grna_gene <-  cells_per_grna |>
    # exclude the non-targeting gRNAs
    dplyr::filter(grna_target != "non-targeting") |>
    # join discovery pairs
    dplyr::left_join(discovery_pairs,
                     unlist(ifelse("grna_id" %in% colnames(discovery_pairs), list(c("grna_id", "grna_target")), "grna_target"))) |>
    # join gene expression df
    dplyr::left_join(baseline_expression_stats, "response_id")

  ################### obtain number of treatment and control cells #############
  # compute the number of treatment cells by grouping grna_target and response_id
  if(!random_assignment){
    grna_gene <- grna_gene |>
      dplyr::rename(mean_num_cells = num_cells) |>
      dplyr::group_by(grna_target, response_id) |>
      dplyr::mutate(num_trt_cells = sum(mean_num_cells),
                    sd_num_cells = 0) |>
      dplyr::ungroup()
  }else{
    grna_gene <- grna_gene |>
      dplyr::group_by(grna_target, response_id) |>
      dplyr::mutate(num_trt_cells = sum(mean_num_cells)) |>
      dplyr::ungroup()
  }

  # define the control cells based on control_group
  if (control_group == "nt_cells") {
    # compute the number of control cells using cells receiving non-targeting gRNAs
    num_cntrl_cells <- cells_per_grna |>
      dplyr::filter(grna_target == "non-targeting") |>
      dplyr::summarize(sum(mean_num_cells)) |>
      dplyr::pull()

    grna_gene <- grna_gene |> dplyr::mutate(num_cntrl_cells = num_cntrl_cells)
  } else { # control_group == "complement"
    grna_gene <- grna_gene |> dplyr::mutate(num_cntrl_cells = num_total_cells - num_trt_cells)
  }

  ################## transform the scalar-valued effect size mean/sd ###########
  if(is.numeric(fold_change)){
    # append the scalar to form a new column in grna_gene
    grna_gene <- grna_gene |> dplyr::mutate(fold_change = fold_change)
  }else{
    # join the grna_gene df and effect size data frames
    grna_gene <- grna_gene |> dplyr::left_join(fold_change, c("grna_id", "response_id"))
  }

  ########################### prepare for multiple testing #####################
  # 1) Collect the list-column version
  enhancer_gene_dt <- grna_gene |>
    dtplyr::lazy_dt() |>
    dplyr::group_by(grna_target, response_id) |>
    dplyr::summarise(
      # compute mean and sd of the test statistic for each pair
      test_stat_distribution = list({
        if(random_assignment){
          compute_distribution_teststat_fixed_es_random_assignment_cpp(
            mean_num_cells = mean_num_cells,
            sd_num_cells = sd_num_cells,
            num_cntrl_cells = num_cntrl_cells,
            expression_mean = expression_mean,
            expression_size = expression_size,
            fold_change = fold_change
          )
        }else{
          compute_distribution_teststat_fixed_es_cpp(
            num_trt_cells = num_trt_cells,
            num_cntrl_cells = num_cntrl_cells,
            num_cells = mean_num_cells,
            expression_mean = expression_mean,
            expression_size = expression_size,
            fold_change = fold_change
          )
        }
      }),

      # compute QC probability
      QC_prob = compute_QC_fixed_es_cpp(
        fold_change = fold_change,
        expression_mean = expression_mean,
        expression_size = expression_size,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = mean_num_cells,
        n_nonzero_trt_thresh = n_nonzero_trt_thresh,
        n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh)
    ) |>
    dplyr::ungroup() |>
    dplyr::collect()

  # 2) Expand the list‐column
  data.table::setDT(enhancer_gene_dt)
  dist_dt <- data.table::rbindlist(enhancer_gene_dt$test_stat_distribution)

  # 3) Rename those two cols
  data.table::setnames(
    dist_dt,
    old = c("mean", "sd"),
    new = c("mean_test_stat", "sd_test_stat")
  )

  # 4) Stick them back on
  enhancer_gene <- cbind(enhancer_gene_dt[, !"test_stat_distribution"], dist_dt) |> tibble::as_tibble()

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
