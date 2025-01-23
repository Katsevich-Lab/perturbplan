#' Input checking function (currently only available for input_check)
#'
#' @inheritParams compute_power
#'
#' @return NULL
input_check <- function(
    ######################## power-determining parameters ######################
    target_cell_df = NULL, control_cell_vec = NULL,
    size_parameter = NULL, baseline_expression = NULL,
    effect_size_mean = NULL, effect_size_sd = NULL,

    ######################## QC-related parameters #############################
    n_nonzero_trt_thresh = NULL, n_nonzero_cntrl_thresh = NULL,

    ####################### test-related parameters ############################
    side = NULL, multiple_testing_method = NULL, multiple_testing_alpha = NULL, cutoff = NULL,

    ###################### output parameter ####################################
    intermediate_outcome = NULL
){
  ########################## target_cell_df ####################################
  # argument missingness checking
  if(is.null(target_cell_df)){
    stop("target_cell_df should be specified!")
  }
  # data frame checking
  if(!is.data.frame(target_cell_df)){
    stop("target_cell_df should be a data frame!")
  }
  # colnames checking
  if(!dplyr::setequal(colnames(target_cell_df), c("gRNA_id", "gRNA_target", "num_cells"))){
    stop("Column names of target_cell_df should be exactly: gRNA_id, gRNA_target and num_cells!")
  }
  # data type checking
  if(!all(target_cell_df$num_cells %% 1 == 0)){
    stop("Column num_cells in target_cell_df should be integer!")
  }

  ########################## control_cell_vec ##################################
  # argument missingness checking
  if(is.null(control_cell_vec)){
    stop("control_cell_vec should be specified!")
  }
  # vector checking
  if(!is.vector(control_cell_vec)){
    stop("control_cell_vec should be a vector!")
  }
  # names missingness checking
  if(any(is.null(names(control_cell_vec)))){
    stop("Every position in control_cell_vec should be named with a particular genome element of interest!")
  }
  # unique names checking
  if(length(control_cell_vec) > length(unique(names(control_cell_vec)))){
    stop("Repeated control cell size for the same genome element is not allowed!")
  }
  ## interaction checking with target_cell_df
  # set of gRNA_target should be same for control and treatment cell
  if(!dplyr::setequal(names(control_cell_vec), unique(target_cell_df$gRNA_target))){
    stop("gRNA_target in target_cell_df should match the names of the contro_cell_vec if the same set of genome elements are of interest!")
  }

  ########################## size_parameter ####################################
  # argument missingness checking
  if(is.null(size_parameter)){
    stop("Size parameter for genes of interest should be specified!")
  }
  # vector checking
  if(!is.vector(size_parameter)){
    stop("size_parameter should be a vector!")
  }
  # names missingness checking
  if(any(is.null(names(size_parameter)))){
    stop("Every position in size_parameter should be named with a particular gene of interest!")
  }
  # unique names checking
  if(length(size_parameter) > length(unique(names(size_parameter)))){
    stop("Repeated size parameter for the same gene is not allowed!")
  }
  # positivity checking
  if(any(size_parameter <= 0)){
    stop("elements in size_parameter should all be positive!")
  }

  ########################## baseline_expression ###############################
  # argument missingness checking
  if(is.null(baseline_expression)){
    stop("baseline_expression should be specified!")
  }
  # vector checking
  if(!is.vector(baseline_expression)){
    stop("baseline_expression should be a vector!")
  }
  # names missingness checking
  if(any(is.null(names(baseline_expression)))){
    stop("Every position in baseline_expression should be named with a particular gene of interest!")
  }
  # unique names checking
  if(length(baseline_expression) > length(unique(names(baseline_expression)))){
    stop("Repeated baseline expression for the same gene is not allowed!")
  }
  # positivity checking
  if(any(baseline_expression <= 0)){
    stop("Elements in baseline_expression should all be positive!")
  }
  ## interaction checking with size_parameter
  # gene set should be same for size_parameter and baseline_expression
  if(!dplyr::setequal(names(baseline_expression), names(size_parameter))){
    stop("Set of genes provided in baseline_expression and size_parameter should be same!")
  }

  ########################## effect_size_mean ##################################
  # argument missingness checking
  if(is.null(effect_size_mean)){
    stop("Mean effect size matrix/scalar should be specified!")
  }
  # data type checking
  if(is.matrix(effect_size_mean)){
    # if matrix, rownames missingness checking
    if(any(is.null(rownames(effect_size_mean)))){
      stop("Rownames in matrix effect_size_mean should be specified!")
    }
    # colnames missingness checking
    if(any(is.null(colnames(effect_size_mean)))){
      stop("Colnames in matrix effect_size_mean should be specified!")
    }
    ## interaction checking with control_cell_vec
    if(!dplyr::setequal(rownames(effect_size_mean), names(control_cell_vec))){
      stop("Rownames of effect_size_mean should equal to the names of control_cell_vec containing all the genome elements of interest!")
    }
    ## interaction checking with baseline_expression
    if(!dplyr::setequal(colnames(effect_size_mean), names(baseline_expression))){
      stop("Colnames of effect_size_mean should equal to the names of baseline_expression!")
    }
  }else{
    # test if a scalar is inputted
    if(any(!is.numeric(effect_size_mean), (length(effect_size_mean) > 1))){
      stop("Either a scalar or a matrix should be inputted for effect_size_mean!")
    }
  }
  # value should be positive
  if(min(effect_size_mean < 0)){
    stop("The fold changes should be a positive number!")
  }

  ########################## effect_size_sd ####################################
  # argument missingness checking
  if(is.null(effect_size_sd)){
    stop("Sd effect size matrix/scalar should be specified!")
  }
  ## interaction checking with effect_size_mean
  if(!all(class(effect_size_sd) == class(effect_size_mean))){
    stop("Data type of effect_size_sd and effect_size_mean should be same!")
  }
  if(is.matrix(effect_size_sd)){
    # test if effect_size_sd and effect_size_mean have the same set of rownames and colnames
    if(dplyr::setequal(rownames(effect_size_sd), rownames(effect_size_mean))){
      if(!dplyr::setequal(colnames(effect_size_sd), colnames(effect_size_mean))){
        stop("Genes in effect_size_sd should be same as those in effect_size_mean!")
      }
    }else{
      stop("Genome elements in effect_size_sd should be same as those in effect_size_mean!")
    }
  }else{
    # test if both effect_size_sd and effect_size_mean have the same length
    if(length(effect_size_sd) != length(effect_size_mean)){
      stop("Length of effect_size_sd should be same as that of effect_size_mean!")
    }
  }

  #################### test-related parameters #################################
  # argument missingness checking
  if(is.null(cutoff)){
    if(any(is.null(c(multiple_testing_method, multiple_testing_alpha)))){
      stop("Either the cutoff should be specified, or the multiple testing method should be specified!")
    }else{
      # multple testing method should be BH or Bonferroni
      if(!(multiple_testing_method %in% c("BH", "Bonferroni"))){
        stop("Multple testing method should be either BH or bonferroni!")
      }
    }
  }else{
    # cutoff should be between 0 and 1
    if(any(cutoff >= 1, cutoff <= 0)){
      stop("Input cutoff should be strictly between 0 and 1!")
    }
  }

  return(NULL)
}


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
