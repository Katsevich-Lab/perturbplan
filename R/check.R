#' Input checking function for compute_power_posthoc
#'
#' @inheritParams compute_power_posthoc
#'
#' @return NULL
input_check_posthoc <- function(
    discovery_pairs = NULL,
    cells_per_grna = NULL,
    baseline_expression_stats = NULL,
    control_group = NULL,
    fold_change_mean = NULL,
    fold_change_sd = NULL,
    num_total_cells = NULL,
    cutoff = NULL,
    n_nonzero_trt_thresh = NULL,
    n_nonzero_cntrl_thresh = NULL,
    side = NULL,
    multiple_testing_method = NULL,
    multiple_testing_alpha = NULL
){
  ############################## discovery_pairs ##############################
  if (is.null(discovery_pairs)) {
    stop("`discovery_pairs` must be specified!")
  }
  if (!is.data.frame(discovery_pairs)) {
    stop("`discovery_pairs` must be a data frame!")
  }
  needed_cols_dp <- c("grna_target", "response_id")
  if (!all(needed_cols_dp %in% colnames(discovery_pairs))) {
    stop("`discovery_pairs` must contain columns: ",
         paste(needed_cols_dp, collapse = ", "), "!")
  }
  distinct_discovery_pairs <- discovery_pairs |> dplyr::distinct()
  if(nrow(distinct_discovery_pairs) != nrow(discovery_pairs)){
    stop("There are duplicate rows in `discovery_pairs`!")
  }
  discovered_targets <- unique(discovery_pairs$grna_target)
  discovered_genes   <- unique(discovery_pairs$response_id)

  ############################## cells_per_grna ###############################
  if (is.null(cells_per_grna)) {
    stop("`cells_per_grna` must be specified!")
  }
  if (!is.data.frame(cells_per_grna)) {
    stop("`cells_per_grna` must be a data frame!")
  }
  needed_cols_cpg <- c("grna_id", "grna_target", "num_cells")
  if (!all(needed_cols_cpg %in% colnames(cells_per_grna))) {
    stop("`cells_per_grna` must contain columns: ",
         paste(needed_cols_cpg, collapse = ", "), "!")
  }
  # Check num_cells is integer
  if (!all(cells_per_grna$num_cells %% 1 == 0)) {
    stop("Column `num_cells` in `cells_per_grna` should be integer!")
  }
  # Cross-check that all discovery_pairs$grna_target are in cells_per_grna
  cpg_targets <- unique(cells_per_grna$grna_target)
  if (!all(discovered_targets %in% cpg_targets)) {
    stop("All `grna_target` values in `discovery_pairs` must appear in `cells_per_grna`!")
  }

  ############################ baseline_expression_stats ######################
  if (is.null(baseline_expression_stats)) {
    stop("`baseline_expression_stats` must be specified!")
  }
  if (!is.data.frame(baseline_expression_stats)) {
    stop("`baseline_expression_stats` must be a data frame!")
  }
  needed_cols_be <- c("response_id", "expression_mean", "expression_size")
  if (!all(needed_cols_be %in% colnames(baseline_expression_stats))) {
    stop("`baseline_expression_stats` must contain columns: ",
         paste(needed_cols_be, collapse = ", "), "!")
  }
  # Cross-check that all discovery_pairs$response_id are in baseline_expression_stats
  be_genes <- unique(baseline_expression_stats$response_id)
  if (!all(discovered_genes %in% be_genes)) {
    stop("All `response_id` values in `discovery_pairs` must appear in `baseline_expression_stats`!")
  }

  ############################## control_group ###############################
  valid_control_groups <- c("complement", "nt_cells")
  if (!control_group %in% valid_control_groups) {
    stop("`control_group` must be either 'complement' or 'nt_cells'!")
  }else{
    if((control_group == "nt_cells") && !"nt_cells" %in% cells_per_grna$grna_target){
      stop("Non-targeting cells should appear in the grna_target of `cells_per_grna`!")
    }
  }

  ############################## fold_change_mean ############################
  if (is.null(fold_change_mean)) {
    stop("`fold_change_mean` must be specified (scalar or data frame)!")
  }

  # Helper to check for minimal columns
  check_pairs_cols <- function(df, df_name) {
    req_cols <- c("grna_target", "response_id")
    if (!all(req_cols %in% colnames(df))) {
      stop("Data frame `", df_name, "` must contain columns: ",
           paste(req_cols, collapse = ", "), "!")
    }
  }

  # Check if numeric or data frame
  if (is.numeric(fold_change_mean)) {
    if (length(fold_change_mean) != 1) {
      stop("If `fold_change_mean` is numeric, it should be a single scalar!")
    }
    if (fold_change_mean <= 0) {
      stop("`fold_change_mean` scalar must be > 0!")
    }
  } else if (is.data.frame(fold_change_mean)) {
    check_pairs_cols(fold_change_mean, "fold_change_mean")
    if (!("fold_change_mean" %in% colnames(fold_change_mean))) {
      stop("Data frame `fold_change_mean` must have a column named `fold_change_mean`!")
    }
    if (any(fold_change_mean$fold_change_mean <= 0, na.rm = TRUE)) {
      stop("All `fold_change_mean` values in the data frame must be positive!")
    }
    # Ensure all (grna_target, response_id) in discovery_pairs exist in fold_change_mean
    fcmean_pairs <- unique(fold_change_mean[, c("grna_target", "response_id")])
    missing_fcmean <- setdiff(
      paste(discovery_pairs$grna_target, discovery_pairs$response_id, sep = "_"),
      paste(fcmean_pairs$grna_target, fcmean_pairs$response_id, sep = "_")
    )
    if (length(missing_fcmean) > 0) {
      stop("Some (grna_target, response_id) in `discovery_pairs` are missing in `fold_change_mean` data frame!")
    }
  } else {
    stop("`fold_change_mean` must be either a scalar numeric or a data frame!")
  }

  ############################## fold_change_sd ##############################
  if (is.null(fold_change_sd)) {
    stop("`fold_change_sd` must be specified (scalar or data frame)!")
  }
  # Must match data type of fold_change_mean
  if (is.numeric(fold_change_mean) && !is.numeric(fold_change_sd)) {
    stop("If `fold_change_mean` is numeric, `fold_change_sd` must also be numeric!")
  }
  if (is.data.frame(fold_change_mean) && !is.data.frame(fold_change_sd)) {
    stop("If `fold_change_mean` is a data frame, `fold_change_sd` must also be a data frame!")
  }

  # Check constraints
  if (is.numeric(fold_change_sd)) {
    if (length(fold_change_sd) != 1) {
      stop("If `fold_change_sd` is numeric, it should be a single scalar!")
    }
    if (fold_change_sd < 0) {
      stop("`fold_change_sd` scalar must be >= 0!")
    }
  } else if (is.data.frame(fold_change_sd)) {
    check_pairs_cols(fold_change_sd, "fold_change_sd")
    if (!("fold_change_sd" %in% colnames(fold_change_sd))) {
      stop("Data frame `fold_change_sd` must have a column named `fold_change_sd`!")
    }
    if (any(fold_change_sd$fold_change_sd < 0, na.rm = TRUE)) {
      stop("All `fold_change_sd` values in the data frame must be nonnegative!")
    }
    # Ensure all (grna_target, response_id) in discovery_pairs exist in fold_change_sd
    fcsd_pairs <- unique(fold_change_sd[, c("grna_target", "response_id")])
    missing_fcsd <- setdiff(
      paste(discovery_pairs$grna_target, discovery_pairs$response_id, sep = "_"),
      paste(fcsd_pairs$grna_target, fcsd_pairs$response_id, sep = "_")
    )
    if (length(missing_fcsd) > 0) {
      stop("Some (grna_target, response_id) in `discovery_pairs` are missing in `fold_change_sd` data frame!")
    }
  }

  ############################## num_total_cells #############################
  # Only relevant if control_group == "complement"
  if (control_group == "complement") {
    if (is.null(num_total_cells)) {
      stop("`num_total_cells` must be specified when `control_group = 'complement'`!")
    }
    if (!is.numeric(num_total_cells) || length(num_total_cells) != 1) {
      stop("`num_total_cells` must be a single numeric value!")
    }
    if (num_total_cells <= 0) {
      stop("`num_total_cells` must be strictly positive!")
    }
  }

  ############################## cutoff / MTP ###############################
  if (is.null(cutoff)) {
    # If cutoff is not specified, we require multiple_testing_method & multiple_testing_alpha
    valid_methods <- c("BH", "bonferroni")
    if (!multiple_testing_method %in% valid_methods) {
      stop("`multiple_testing_method` must be either 'BH' or 'bonferroni' if `cutoff` is NULL!")
    }
    if (!is.numeric(multiple_testing_alpha) || length(multiple_testing_alpha) != 1 ||
        multiple_testing_alpha <= 0 || multiple_testing_alpha >= 1) {
      stop("`multiple_testing_alpha` must be a single numeric in (0,1)!")
    }
  } else {
    if (!is.numeric(cutoff) || length(cutoff) != 1) {
      stop("`cutoff` must be a single numeric value!")
    }
    if (cutoff <= 0 || cutoff >= 1) {
      stop("`cutoff` must be strictly between 0 and 1!")
    }
  }

  ############################# side parameter ###############################
  valid_sides <- c("left", "right", "both")
  if (!side %in% valid_sides) {
    stop("`side` must be one of 'left', 'right', or 'both'!")
  }

  ############################# QC thresholds ################################
  if (!is.numeric(n_nonzero_trt_thresh) ||
      length(n_nonzero_trt_thresh) != 1 ||
      n_nonzero_trt_thresh < 0 ||
      n_nonzero_trt_thresh %% 1 != 0) {
    stop("`n_nonzero_trt_thresh` must be a nonnegative integer!")
  }
  if (!is.numeric(n_nonzero_cntrl_thresh) ||
      length(n_nonzero_cntrl_thresh) != 1 ||
      n_nonzero_cntrl_thresh < 0 ||
      n_nonzero_cntrl_thresh %% 1 != 0) {
    stop("`n_nonzero_cntrl_thresh` must be a nonnegative integer!")
  }

  # If all checks pass, return NULL (invisibly)
  invisible(NULL)
}



#' Input checking function for power_function
#'
#' @inheritParams power_function
#'
#' @return NULL
input_check_power_function <- function(
    recovery_rate = NULL,
    num_total_reads = NULL,
    mapping_efficiency = NULL,
    cells_per_grna = NULL,
    baseline_relative_expression_stats = NULL,
    fold_change_mean = NULL,
    fold_change_sd = NULL,
    num_planned_cells = NULL,
    control_group = NULL,
    umi_per_cell = NULL,
    umi_variation = NULL,
    side = NULL,
    multiple_testing_method = NULL,
    multiple_testing_alpha = NULL,
    cutoff = NULL,
    discovery_pairs = NULL,
    n_nonzero_trt_thresh = NULL,
    n_nonzero_cntrl_thresh = NULL
) {
  ####################### Experimental Design Checks ##########################
  if (is.null(recovery_rate) || !is.numeric(recovery_rate) || recovery_rate <= 0 || recovery_rate > 1) {
    stop("`recovery_rate` must be a numeric value in (0,1].")
  }
  if (is.null(num_total_reads) || !is.numeric(num_total_reads) || num_total_reads <= 0) {
    stop("`num_total_reads` must be a positive numeric value.")
  }
  if (is.null(mapping_efficiency) || !is.numeric(mapping_efficiency) || mapping_efficiency <= 0 || mapping_efficiency > 1) {
    stop("`mapping_efficiency` must be a numeric value in (0,1].")
  }

  ###################### Power-Determining Parameter Checks ######################
  if (is.null(cells_per_grna) || !is.data.frame(cells_per_grna)) {
    stop("`cells_per_grna` must be a specified data frame.")
  }
  if (is.null(baseline_relative_expression_stats) || !is.data.frame(baseline_relative_expression_stats)) {
    stop("`baseline_relative_expression_stats` must be a specified data frame.")
  }
  if (is.null(num_planned_cells) || !is.numeric(num_planned_cells) || num_planned_cells <= 0) {
    stop("`num_planned_cells` must be a positive numeric value.")
  }
  if (is.null(umi_per_cell) || !is.numeric(umi_per_cell) || umi_per_cell <= 0) {
    stop("`umi_per_cell` must be a positive numeric value.")
  }
  if (is.null(umi_variation) || !is.numeric(umi_variation) || umi_variation < 0 || umi_variation > 1) {
    stop("`umi_variation` must be a numeric value in [0,1].")
  }

  ###################### Effect Sizes Checks ######################
  if (is.null(fold_change_mean)) {
    stop("`fold_change_mean` must be provided (numeric or data frame).")
  }
  if (is.null(fold_change_sd)) {
    stop("`fold_change_sd` must be provided (numeric or data frame).")
  }

  ###################### Control Group Checks ######################
  valid_control_groups <- c("complement", "nt_cells")
  if (is.null(control_group) || !control_group %in% valid_control_groups) {
    stop("`control_group` must be either 'complement' or 'nt_cells'.")
  }

  ###################### Multiple Testing & Side Checks ######################
  valid_methods <- c("BH", "bonferroni")
  if (!is.null(multiple_testing_method) && !multiple_testing_method %in% valid_methods) {
    stop("`multiple_testing_method` must be either 'BH' or 'bonferroni'.")
  }
  if (!is.null(multiple_testing_alpha) && (!is.numeric(multiple_testing_alpha) || multiple_testing_alpha <= 0 || multiple_testing_alpha >= 1)) {
    stop("`multiple_testing_alpha` must be a numeric in (0,1).")
  }
  valid_sides <- c("left", "right", "both")
  if (is.null(side) || !side %in% valid_sides) {
    stop("`side` must be one of 'left', 'right', or 'both'.")
  }

  ###################### QC Thresholds ######################
  if (!is.numeric(n_nonzero_trt_thresh) || length(n_nonzero_trt_thresh) != 1 || n_nonzero_trt_thresh < 0 || n_nonzero_trt_thresh %% 1 != 0) {
    stop("`n_nonzero_trt_thresh` must be a nonnegative integer.")
  }
  if (!is.numeric(n_nonzero_cntrl_thresh) || length(n_nonzero_cntrl_thresh) != 1 || n_nonzero_cntrl_thresh < 0 || n_nonzero_cntrl_thresh %% 1 != 0) {
    stop("`n_nonzero_cntrl_thresh` must be a nonnegative integer.")
  }

  ###################### Discovery Pairs Check ######################
  if (is.null(discovery_pairs) || !is.data.frame(discovery_pairs)) {
    stop("`discovery_pairs` must be a specified data frame.")
  }

  invisible(NULL)
}


#' Input checking function for library_computation
#'
#' @inheritParams library_computation
#'
#' @return NULL
input_check_library_computation <- function(
    QC_data = NULL,
    downsample_ratio = NULL,
    D2_rough = NULL
) {
  ###################### QC_data Check ######################
  if (is.null(QC_data) || !is.data.frame(QC_data)) {
    stop("`QC_data` must be a specified data frame.")
  }
  if (!"num_reads" %in% colnames(QC_data)) {
    stop("`QC_data` must contain a `num_reads` column.")
  }
  if (!"UMI_id" %in% colnames(QC_data)) {
    stop("`QC_data` must contain a `UMI_id` column.")
  }
  if (!"cell_id" %in% colnames(QC_data)) {
    stop("`QC_data` must contain a `cell_id` column.")
  }
  if (!"response_id" %in% colnames(QC_data)) {
    stop("`QC_data` must contain a `response_id` column.")
  }
  if (nrow(QC_data) == 0) {
    stop("`QC_data` cannot be empty.")
  }

  ###################### downsample_ratio Check ######################
  if (is.null(downsample_ratio) || !is.numeric(downsample_ratio) || downsample_ratio <= 0 || downsample_ratio > 1) {
    stop("`downsample_ratio` must be a numeric value in (0,1].")
  }

  ###################### D2_rough Check ######################
  if (is.null(D2_rough) || !is.numeric(D2_rough) || D2_rough < 0 || D2_rough > 1) {
    stop("`D2_rough` must be a numeric value in [0,1].")
  }

  invisible(NULL)
}

