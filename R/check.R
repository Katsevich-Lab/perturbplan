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
    UMI_per_cell = NULL,
    variation = NULL,
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
  if (is.null(UMI_per_cell) || !is.numeric(UMI_per_cell) || UMI_per_cell <= 0) {
    stop("`UMI_per_cell` must be a positive numeric value.")
  }
  if (is.null(variation) || !is.numeric(variation) || variation < 0 || variation > 1) {
    stop("`variation` must be a numeric value in [0,1].")
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

#' Input checking function for compute_power_plan_full_grid
#'
#' @inheritParams compute_power_plan_full_grid
#'
#' @return NULL
input_check_compute_power_plan_full_grid <- function(
    tpm_threshold, minimum_fold_change, cells_per_target, reads_per_cell,
    MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4, gRNA_variability = 0.13,
    control_group = "complement", side = "left", multiple_testing_alpha = 0.05, prop_non_null = 0.1,
    baseline_expression_stats, library_parameters,
    grid_size = 10, min_power_threshold = 0.01, max_power_threshold = 0.8,
    mapping_efficiency = 0.72
) {
  
  ############################ tpm_threshold #############################
  if (missing(tpm_threshold)) {
    stop("`tpm_threshold` must be specified!")
  }
  if (is.numeric(tpm_threshold)) {
    if (any(tpm_threshold < 0)) {
      stop("`tpm_threshold` values must be non-negative!")
    }
    if (any(tpm_threshold > 1000000)) {
      stop("`tpm_threshold` values seem unreasonably large (>1,000,000 TPM)!")
    }
  } else if (is.character(tpm_threshold)) {
    if (length(tpm_threshold) != 1 || tpm_threshold != "varying") {
      stop("`tpm_threshold` must be numeric or the string 'varying'!")
    }
  } else {
    stop("`tpm_threshold` must be numeric or the string 'varying'!")
  }
  
  ######################## minimum_fold_change ############################
  if (missing(minimum_fold_change)) {
    stop("`minimum_fold_change` must be specified!")
  }
  if (is.numeric(minimum_fold_change)) {
    if (any(minimum_fold_change <= 0)) {
      stop("`minimum_fold_change` values must be positive!")
    }
  } else if (is.character(minimum_fold_change)) {
    if (length(minimum_fold_change) != 1 || minimum_fold_change != "varying") {
      stop("`minimum_fold_change` must be numeric or the string 'varying'!")
    }
  } else {
    stop("`minimum_fold_change` must be numeric or the string 'varying'!")
  }
  
  ########################## cells_per_target ##############################
  if (missing(cells_per_target)) {
    stop("`cells_per_target` must be specified!")
  }
  if (is.numeric(cells_per_target)) {
    if (any(cells_per_target <= 0)) {
      stop("`cells_per_target` values must be positive!")
    }
    if (any(cells_per_target != round(cells_per_target))) {
      stop("`cells_per_target` values must be integers!")
    }
  } else if (is.character(cells_per_target)) {
    if (length(cells_per_target) != 1 || cells_per_target != "varying") {
      stop("`cells_per_target` must be numeric or the string 'varying'!")
    }
  } else {
    stop("`cells_per_target` must be numeric or the string 'varying'!")
  }
  
  ########################## reads_per_cell #################################
  if (missing(reads_per_cell)) {
    stop("`reads_per_cell` must be specified!")
  }
  if (is.numeric(reads_per_cell)) {
    if (any(reads_per_cell <= 0)) {
      stop("`reads_per_cell` values must be positive!")
    }
    if (any(reads_per_cell != round(reads_per_cell))) {
      stop("`reads_per_cell` values must be integers!")
    }
  } else if (is.character(reads_per_cell)) {
    if (length(reads_per_cell) != 1 || reads_per_cell != "varying") {
      stop("`reads_per_cell` must be numeric or the string 'varying'!")
    }
  } else {
    stop("`reads_per_cell` must be numeric or the string 'varying'!")
  }
  
  ###################### Experimental parameters ############################
  if (!is.numeric(MOI) || length(MOI) != 1 || MOI <= 0) {
    stop("`MOI` must be a positive numeric value!")
  }
  
  if (!is.numeric(num_targets) || length(num_targets) != 1 || num_targets <= 0 || num_targets != round(num_targets)) {
    stop("`num_targets` must be a positive integer!")
  }
  
  if (!is.numeric(non_targeting_gRNAs) || length(non_targeting_gRNAs) != 1 || non_targeting_gRNAs < 0 || non_targeting_gRNAs != round(non_targeting_gRNAs)) {
    stop("`non_targeting_gRNAs` must be a non-negative integer!")
  }
  
  if (!is.numeric(gRNAs_per_target) || length(gRNAs_per_target) != 1 || gRNAs_per_target <= 0 || gRNAs_per_target != round(gRNAs_per_target)) {
    stop("`gRNAs_per_target` must be a positive integer!")
  }
  
  if (!is.numeric(gRNA_variability) || length(gRNA_variability) != 1 || gRNA_variability < 0) {
    stop("`gRNA_variability` must be a non-negative numeric value!")
  }
  
  ###################### Analysis parameters ###############################
  if (!is.character(control_group) || length(control_group) != 1 || !control_group %in% c("complement", "nt_cells")) {
    stop("`control_group` must be either 'complement' or 'nt_cells'!")
  }
  
  if (!is.character(side) || length(side) != 1 || !side %in% c("left", "right", "both")) {
    stop("`side` must be 'left', 'right', or 'both'!")
  }
  
  if (!is.numeric(multiple_testing_alpha) || length(multiple_testing_alpha) != 1 || multiple_testing_alpha <= 0 || multiple_testing_alpha >= 1) {
    stop("`multiple_testing_alpha` must be a numeric value in (0,1)!")
  }
  
  if (!is.numeric(prop_non_null) || length(prop_non_null) != 1 || prop_non_null <= 0 || prop_non_null > 1) {
    stop("`prop_non_null` must be a numeric value in (0,1]!")
  }
  
  ###################### Data inputs #####################################
  if (missing(baseline_expression_stats)) {
    stop("`baseline_expression_stats` must be specified!")
  }
  if (!is.data.frame(baseline_expression_stats)) {
    stop("`baseline_expression_stats` must be a data frame!")
  }
  required_cols_baseline <- c("response_id", "relative_expression", "expression_size")
  if (!all(required_cols_baseline %in% colnames(baseline_expression_stats))) {
    stop("`baseline_expression_stats` must contain columns: ",
         paste(required_cols_baseline, collapse = ", "), "!")
  }
  if (nrow(baseline_expression_stats) == 0) {
    stop("`baseline_expression_stats` cannot be empty!")
  }
  if (any(baseline_expression_stats$relative_expression <= 0, na.rm = TRUE)) {
    stop("`relative_expression` values must be positive!")
  }
  if (any(baseline_expression_stats$expression_size <= 0, na.rm = TRUE)) {
    stop("`expression_size` values must be positive!")
  }
  
  if (missing(library_parameters)) {
    stop("`library_parameters` must be specified!")
  }
  if (!is.list(library_parameters)) {
    stop("`library_parameters` must be a list!")
  }
  required_params_library <- c("UMI_per_cell", "variation")
  if (!all(required_params_library %in% names(library_parameters))) {
    stop("`library_parameters` must contain elements: ",
         paste(required_params_library, collapse = ", "), "!")
  }
  if (!is.numeric(library_parameters$UMI_per_cell) || library_parameters$UMI_per_cell <= 0) {
    stop("`library_parameters$UMI_per_cell` must be a positive numeric value!")
  }
  if (!is.numeric(library_parameters$variation) || library_parameters$variation <= 0) {
    stop("`library_parameters$variation` must be a positive numeric value!")
  }
  
  ###################### Grid parameters ################################
  if (!is.numeric(grid_size) || length(grid_size) != 1 || grid_size <= 0 || grid_size != round(grid_size)) {
    stop("`grid_size` must be a positive integer!")
  }
  
  if (!is.numeric(min_power_threshold) || length(min_power_threshold) != 1 || min_power_threshold <= 0 || min_power_threshold >= 1) {
    stop("`min_power_threshold` must be a numeric value in (0,1)!")
  }
  
  if (!is.numeric(max_power_threshold) || length(max_power_threshold) != 1 || max_power_threshold <= 0 || max_power_threshold >= 1) {
    stop("`max_power_threshold` must be a numeric value in (0,1)!")
  }
  
  if (min_power_threshold >= max_power_threshold) {
    stop("`min_power_threshold` must be less than `max_power_threshold`!")
  }
  
  if (!is.numeric(mapping_efficiency) || length(mapping_efficiency) != 1 || mapping_efficiency <= 0 || mapping_efficiency > 1) {
    stop("`mapping_efficiency` must be a numeric value in (0,1]!")
  }
  
  invisible(NULL)
}

#' Input checking function for cost_power_computation
#'
#' @inheritParams cost_power_computation
#'
#' @return NULL
input_check_cost_power_computation <- function(
    minimizing_variable = "tpm_threshold", fixed_variable = list(minimum_fold_change = 0.8),
    MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4, gRNA_variability = 0.13,
    control_group = "complement", side = "left", multiple_testing_alpha = 0.05, prop_non_null = 0.1,
    baseline_expression_stats, library_parameters,
    grid_size = 20, power_target = 0.8, power_precision = 0.01, min_power = 0.05, max_power = 0.95,
    cost_precision = 0.9,
    cost_per_captured_cell = 0.086, cost_per_million_reads = 0.374, cost_constraint = NULL,
    mapping_efficiency = 0.72
) {
  
  ######################## minimizing_variable ##############################
  if (!is.character(minimizing_variable) || length(minimizing_variable) != 1) {
    stop("`minimizing_variable` must be a single character string!")
  }
  valid_minimizing_vars <- c("tpm_threshold", "minimum_fold_change")
  if (!minimizing_variable %in% valid_minimizing_vars) {
    stop("`minimizing_variable` must be one of: ",
         paste(valid_minimizing_vars, collapse = ", "), "!")
  }
  
  ########################## fixed_variable #################################
  if (!is.list(fixed_variable)) {
    stop("`fixed_variable` must be a list!")
  }
  
  # Check required fixed variables based on minimizing variable
  if (minimizing_variable == "tpm_threshold") {
    if (!"minimum_fold_change" %in% names(fixed_variable)) {
      stop("When minimizing tpm_threshold, `fixed_variable` must contain 'minimum_fold_change'!")
    }
    if (!is.numeric(fixed_variable$minimum_fold_change) || fixed_variable$minimum_fold_change <= 0) {
      stop("`fixed_variable$minimum_fold_change` must be a positive numeric value!")
    }
  } else if (minimizing_variable == "minimum_fold_change") {
    if (!"tpm_threshold" %in% names(fixed_variable)) {
      stop("When minimizing minimum_fold_change, `fixed_variable` must contain 'tpm_threshold'!")
    }
    if (!is.numeric(fixed_variable$tpm_threshold) || fixed_variable$tpm_threshold < 0) {
      stop("`fixed_variable$tpm_threshold` must be a positive numeric value!")
    }
  }
  
  # Check optional fixed variables
  if ("cells_per_target" %in% names(fixed_variable)) {
    if (!is.numeric(fixed_variable$cells_per_target) || fixed_variable$cells_per_target <= 0 || 
        fixed_variable$cells_per_target != round(fixed_variable$cells_per_target)) {
      stop("`fixed_variable$cells_per_target` must be a positive integer!")
    }
  }
  
  if ("reads_per_cell" %in% names(fixed_variable)) {
    if (!is.numeric(fixed_variable$reads_per_cell) || fixed_variable$reads_per_cell <= 0 || 
        fixed_variable$reads_per_cell != round(fixed_variable$reads_per_cell)) {
      stop("`fixed_variable$reads_per_cell` must be a positive integer!")
    }
  }
  
  ###################### Experimental parameters ############################
  if (!is.numeric(MOI) || length(MOI) != 1 || MOI <= 0) {
    stop("`MOI` must be a positive numeric value!")
  }
  
  if (!is.numeric(num_targets) || length(num_targets) != 1 || num_targets <= 0 || num_targets != round(num_targets)) {
    stop("`num_targets` must be a positive integer!")
  }
  
  if (!is.numeric(non_targeting_gRNAs) || length(non_targeting_gRNAs) != 1 || non_targeting_gRNAs < 0 || non_targeting_gRNAs != round(non_targeting_gRNAs)) {
    stop("`non_targeting_gRNAs` must be a non-negative integer!")
  }
  
  if (!is.numeric(gRNAs_per_target) || length(gRNAs_per_target) != 1 || gRNAs_per_target <= 0 || gRNAs_per_target != round(gRNAs_per_target)) {
    stop("`gRNAs_per_target` must be a positive integer!")
  }
  
  if (!is.numeric(gRNA_variability) || length(gRNA_variability) != 1 || gRNA_variability < 0) {
    stop("`gRNA_variability` must be a non-negative numeric value!")
  }
  
  ###################### Analysis parameters ###############################
  if (!is.character(control_group) || length(control_group) != 1 || !control_group %in% c("complement", "nt_cells")) {
    stop("`control_group` must be either 'complement' or 'nt_cells'!")
  }
  
  if (!is.character(side) || length(side) != 1 || !side %in% c("left", "right", "both")) {
    stop("`side` must be 'left', 'right', or 'both'!")
  }
  
  if (!is.numeric(multiple_testing_alpha) || length(multiple_testing_alpha) != 1 || multiple_testing_alpha <= 0 || multiple_testing_alpha >= 1) {
    stop("`multiple_testing_alpha` must be a numeric value in (0,1)!")
  }
  
  if (!is.numeric(prop_non_null) || length(prop_non_null) != 1 || prop_non_null <= 0 || prop_non_null > 1) {
    stop("`prop_non_null` must be a numeric value in (0,1]!")
  }
  
  ###################### Data inputs #####################################
  if (missing(baseline_expression_stats)) {
    stop("`baseline_expression_stats` must be specified!")
  }
  if (!is.data.frame(baseline_expression_stats)) {
    stop("`baseline_expression_stats` must be a data frame!")
  }
  required_cols_baseline <- c("response_id", "relative_expression", "expression_size")
  if (!all(required_cols_baseline %in% colnames(baseline_expression_stats))) {
    stop("`baseline_expression_stats` must contain columns: ",
         paste(required_cols_baseline, collapse = ", "), "!")
  }
  if (nrow(baseline_expression_stats) == 0) {
    stop("`baseline_expression_stats` cannot be empty!")
  }
  if (any(baseline_expression_stats$relative_expression <= 0, na.rm = TRUE)) {
    stop("`relative_expression` values must be positive!")
  }
  if (any(baseline_expression_stats$expression_size <= 0, na.rm = TRUE)) {
    stop("`expression_size` values must be positive!")
  }
  
  if (missing(library_parameters)) {
    stop("`library_parameters` must be specified!")
  }
  if (!is.list(library_parameters)) {
    stop("`library_parameters` must be a list!")
  }
  required_params_library <- c("UMI_per_cell", "variation")
  if (!all(required_params_library %in% names(library_parameters))) {
    stop("`library_parameters` must contain elements: ",
         paste(required_params_library, collapse = ", "), "!")
  }
  if (!is.numeric(library_parameters$UMI_per_cell) || library_parameters$UMI_per_cell <= 0) {
    stop("`library_parameters$UMI_per_cell` must be a positive numeric value!")
  }
  if (!is.numeric(library_parameters$variation) || library_parameters$variation <= 0) {
    stop("`library_parameters$variation` must be a positive numeric value!")
  }
  
  ###################### Power optimization parameters ###################
  if (!is.numeric(grid_size) || length(grid_size) != 1 || grid_size <= 0 || grid_size != round(grid_size)) {
    stop("`grid_size` must be a positive integer!")
  }
  
  if (!is.numeric(power_target) || length(power_target) != 1 || power_target <= 0 || power_target >= 1) {
    stop("`power_target` must be a numeric value in (0,1)!")
  }
  
  if (!is.numeric(power_precision) || length(power_precision) != 1 || power_precision <= 0 || power_precision >= 1) {
    stop("`power_precision` must be a numeric value in (0,1)!")
  }
  
  if (!is.numeric(min_power) || length(min_power) != 1 || min_power <= 0 || min_power >= 1) {
    stop("`min_power` must be a numeric value in (0,1)!")
  }
  
  if (!is.numeric(max_power) || length(max_power) != 1 || max_power <= 0 || max_power >= 1) {
    stop("`max_power` must be a numeric value in (0,1)!")
  }
  
  if (min_power >= max_power) {
    stop("`min_power` must be less than `max_power`!")
  }
  
  if (power_target <= min_power || power_target >= max_power) {
    stop("`power_target` must be between `min_power` and `max_power`!")
  }
  
  ###################### Cost parameters ################################
  if (!is.numeric(cost_precision) || length(cost_precision) != 1 || cost_precision <= 0 || cost_precision > 1) {
    stop("`cost_precision` must be a numeric value in (0,1]!")
  }
  
  if (!is.numeric(cost_per_captured_cell) || length(cost_per_captured_cell) != 1 || cost_per_captured_cell < 0) {
    stop("`cost_per_captured_cell` must be a non-negative numeric value!")
  }
  
  if (!is.numeric(cost_per_million_reads) || length(cost_per_million_reads) != 1 || cost_per_million_reads < 0) {
    stop("`cost_per_million_reads` must be a non-negative numeric value!")
  }
  
  if (!is.null(cost_constraint)) {
    if (!is.numeric(cost_constraint) || length(cost_constraint) != 1 || cost_constraint <= 0) {
      stop("`cost_constraint` must be NULL or a positive numeric value!")
    }
  }
  
  if (!is.numeric(mapping_efficiency) || length(mapping_efficiency) != 1 || mapping_efficiency <= 0 || mapping_efficiency > 1) {
    stop("`mapping_efficiency` must be a numeric value in (0,1]!")
  }
  
  invisible(NULL)
}

