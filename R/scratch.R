
# This is a Rscript computing the power function using score test

#' Compute power for each perturbation-gene pair with fixed fold change
#'
#' This function computes power for perturb-seq experiments with fixed (non-random) gRNA assignment.
#' It uses C++ implementations for computational efficiency.
#'
#' @param discovery_pairs A data frame specifying which element-gene pairs to consider, with columns `grna_target` and `response_id`; it can also have `grna_id` as a column but this is optional
#' @param cells_per_grna A data frame specifying how many cells contain each gRNA, with columns `grna_id`, `grna_target`, and `num_cells`
#' @param baseline_expression_stats A data frame specifying the baseline expression statistics for each gene, with columns `response_id`, `expression_mean`, and `expression_size`
#' @param control_group A character string specifying the control group, either "complement" or "nt_cells"
#' @param fold_change A numeric value or data frame to use for fixed effect size for all gRNA-gene pairs
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
compute_power_posthoc_fixed_fc <- function(
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
    multiple_testing_alpha = 0.1) {

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
  grna_gene <- grna_gene |>
    dplyr::group_by(grna_target, response_id) |>
    dplyr::mutate(num_trt_cells = sum(num_cells)) |>
    dplyr::ungroup()

  # define the control cells based on control_group
  if (control_group == "nt_cells") {
    # compute the number of control cells using cells receiving non-targeting gRNAs
    num_cntrl_cells <- cells_per_grna |>
      dplyr::filter(grna_target == "non-targeting") |>
      dplyr::summarize(sum(num_cells)) |>
      dplyr::pull()

    grna_gene <- grna_gene |> dplyr::mutate(num_cntrl_cells = num_cntrl_cells)
  } else { # control_group == "complement"
    grna_gene <- grna_gene |> dplyr::mutate(num_cntrl_cells = num_total_cells - num_trt_cells)
  }

  ################## transform the scalar-valued effect size ###########
  if(is.numeric(fold_change)){
    # append the scalar to form a new column in grna_gene
    grna_gene <- grna_gene |> dplyr::mutate(fold_change = fold_change)
  }else{
    # join the grna_gene df and effect size data frames
    grna_gene <- grna_gene |> dplyr::left_join(fold_change, c("grna_id", "response_id"))
  }

  ########################### prepare for multiple testing #####################
  # Use data.table for efficiency
  enhancer_gene_dt <- grna_gene |>
    dtplyr::lazy_dt() |>
    dplyr::group_by(grna_target, response_id) |>
    dplyr::summarise(
      # compute mean and sd of the test statistic for each pair
      test_stat_distribution = list(
        compute_distribution_teststat_fixed_es_cpp(
          num_trt_cells = num_trt_cells,
          num_cntrl_cells = num_cntrl_cells,
          num_cells = num_cells,
          expression_mean = expression_mean,
          expression_size = expression_size,
          fold_change = fold_change
        )
      ),

      # compute QC probability
      QC_prob = compute_QC_fixed_es_cpp(
        fold_change = fold_change,
        expression_mean = expression_mean,
        expression_size = expression_size,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = num_cells,
        n_nonzero_trt_thresh = n_nonzero_trt_thresh,
        n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh)
    ) |>
    dplyr::ungroup() |>
    dplyr::collect()

  # Expand the list‐column
  data.table::setDT(enhancer_gene_dt)
  dist_dt <- data.table::rbindlist(enhancer_gene_dt$test_stat_distribution)

  # Rename columns
  data.table::setnames(
    dist_dt,
    old = c("mean", "sd"),
    new = c("mean_test_stat", "sd_test_stat")
  )

  # Combine results
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
  power_values <- rejection_computation_cpp(mean_list = enhancer_gene$mean_test_stat,
                                            sd_list = enhancer_gene$sd_test_stat,
                                            side = side,
                                            cutoff = cutoff) * (1 - enhancer_gene$QC_prob)
  
  enhancer_gene <- enhancer_gene |>
    dplyr::mutate(
      cutoff = cutoff,
      power = power_values
    )

  # store individual power and rejection size as the output
  output <- list(
    individual_power = enhancer_gene |> dplyr::select(grna_target, response_id, power),
    expected_num_discoveries = sum(enhancer_gene$power)
  )

  return(output)
}


# This is a Rscript computing the power function using score test



#' Compute power in underspecified setup
#'
#' @inheritParams calculate_power_grid
#' @param cells_reads_df A dataframe including columns num_total_cells and reads_per_cell
#' @param side (Optional) A character string specifying the side of the test, either "left", "right", or "both"; defaults to "both"
#' @param control_group A character string specifying the control group, either "complement" or "nt_cells"
#' @param B (Optional) A scalar specifying how many Monte Carlo samples used
#'
#' @return A dataframe concatenating the overall_power, power_by_fc and power_by_expr
compute_power_grid <- function(
    cells_reads_df,
    num_targets = 100,
    gRNAs_per_target = 4,
    non_targeting_gRNAs = 10,
    num_pairs = 1000,
    tpm_threshold = 10,
    fdr_target = 0.05,
    fc_mean = 0.85,
    fc_sd = 0.15,
    prop_non_null = 0.1,
    MOI = 10,
    biological_system = "K562",
    experimental_platform = "10x Chromium v3",
    side = "left",
    control_group = "complement",
    B = 500
){

  ############### extract the baseline expression information ##################
  fc_expression_info <- extract_fc_expression_info(fold_change_mean = fc_mean, fold_change_sd = fc_sd,
                                                   biological_system = biological_system, B = B)

  ################## compute the library size information ######################
  # TBD

  ############### compute the power for the cells-reads grid ###################
  power_df <- cells_reads_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      power_output = list(
        compute_underspecified_power(
          # experimental information
          num_total_cells = num_total_cells, library_size = reads_per_cell, MOI = MOI, num_targets = num_targets, gRNAs_per_target = gRNAs_per_target, non_targeting_gRNAs = non_targeting_gRNAs,
          # analysis information
          multiple_testing_alpha = fdr_target, multiple_testing_method = "BH", control_group = control_group, side = side, num_pairs = num_pairs,
          # effect size and baseline expression information
          fc_expression_df = fc_expression_info$fc_expression_df, expression_dispersion_curve = fc_expression_info$expression_dispersion_curve, prop_non_null = prop_non_null)
      ),
      overall_power = power_output$overall_power,
      power_by_fc = list(power_output$power_by_fc),
      power_by_expr = list(power_output$power_by_expr)
    ) |>
    dplyr::select(-power_output)

  # return the output dataframe
  return(power_df)
}


#' Helper function that returns named vector directly (cleaner than original)
#'
#' @param num_trt_cells Number of treatment cells
#' @param num_cntrl_cells Number of control cells
#' @param expression_mean Mean expression level
#' @param expression_size Size parameter for negative binomial
#' @param fold_change_mean Mean fold change
#' @return Named vector with mean and sd of test statistic
compute_test_stat_clean <- function(num_trt_cells, num_cntrl_cells, expression_mean, expression_size, fold_change_mean) {
  # This is the same logic as compute_distribution_teststat but returns named vector directly
  num_trt_cells_sq <- num_trt_cells^2
  fold_change_sd <- 0

  # compute treatment/control cells proportion
  num_test_cells <- num_trt_cells + num_cntrl_cells
  trt_test_prop <- num_trt_cells / num_test_cells
  cntrl_test_prop <- 1 - trt_test_prop

  # define treatment/control/pooled mean expression
  trt_expression_mean <- expression_mean * fold_change_mean
  cntrl_expression_mean <- expression_mean
  pooled_expression_mean <- trt_expression_mean * trt_test_prop + cntrl_expression_mean * cntrl_test_prop

  # compute the square of the denominator in the score statistic
  pooled_var <- var_nb(mean = pooled_expression_mean, size = expression_size)
  denominator_sq <- pooled_var * (1 / num_cntrl_cells + 1 / num_trt_cells)

  ########################## compute the asymptotic sd of test stat ############
  # compute the control group variance
  cntrl_var <- var_nb(mean = cntrl_expression_mean, size = expression_size) / num_cntrl_cells

  # compute the treatment group variance
  trt_var_within <- (trt_expression_mean + (cntrl_expression_mean^2 * (fold_change_sd^2 + fold_change_mean^2) / expression_size)) / num_trt_cells
  trt_var_across <- cntrl_expression_mean^2 * fold_change_sd^2 * num_trt_cells_sq / num_trt_cells^2
  trt_var <- trt_var_within + trt_var_across

  # compute the asymptotic sd
  sd <- sqrt((cntrl_var + trt_var) / denominator_sq)

  ################# compute the asymptotic mean of test stat ###################
  mean <- cntrl_expression_mean * (fold_change_mean - 1) / sqrt(denominator_sq)

  # Return named vector directly (no list wrapper!)
  return(c(mean = mean, sd = sd))
}


#' Compute power grid using separated Monte Carlo approach
#'
#' This function separates Monte Carlo integration accuracy from output grid resolution,
#' providing an efficient and flexible approach to power analysis.
#'
#' @param cells_reads_df Data frame with columns num_total_cells and reads_per_cell
#' @param num_targets Number of targets to test
#' @param gRNAs_per_target Number of gRNAs per target
#' @param non_targeting_gRNAs Number of non-targeting gRNAs
#' @param num_pairs Number of pairs for multiple testing
#' @param tpm_threshold TPM threshold (currently unused)
#' @param fdr_target Target false discovery rate
#' @param fc_mean Mean fold change for effect size distribution
#' @param fc_sd Standard deviation of fold change distribution
#' @param prop_non_null Proportion of non-null hypotheses
#' @param MOI Multiplicity of infection
#' @param biological_system Biological system for baseline expression
#' @param experimental_platform Experimental platform
#' @param side Test sidedness ("left", "right", "both")
#' @param control_group Control group type ("complement" or "nt_cells")
#' @param B Number of Monte Carlo samples for integration
#' @param fc_curve_points Number of points for fold change curve
#' @param expr_curve_points Number of points for expression curve
#' @return Data frame with power analysis results
#' @export
compute_power_grid_separated <- function(
    cells_reads_df,
    num_targets = 100,
    gRNAs_per_target = 4,
    non_targeting_gRNAs = 10,
    num_pairs = 1000,
    tpm_threshold = 10,
    fdr_target = 0.05,
    fc_mean = 0.85,
    fc_sd = 0.15,
    prop_non_null = 0.1,
    MOI = 10,
    biological_system = "K562",
    experimental_platform = "10x Chromium v3",
    side = "left",
    control_group = "complement",
    B = 500,
    fc_curve_points = 10,
    expr_curve_points = 10
){

  ############### extract Monte Carlo samples for integration ##################
  set.seed(1)  # Reproducible results
  fc_expression_info <- extract_fc_expression_info(
    fold_change_mean = fc_mean, fold_change_sd = fc_sd,
    biological_system = biological_system, B = B
  )

  # Define systematic output grids
  fc_range <- range(fc_expression_info$fc_expression_df$fold_change)
  expr_range <- range(fc_expression_info$fc_expression_df$relative_expression)

  # Use directional fold change grid based on test sidedness
  if (side == "left") {
    # For left-sided tests (knockdown), use FC < 1
    fc_output_grid <- seq(min(fc_range[1], 0.5), 1, length.out = fc_curve_points)
  } else if (side == "right") {
    # For right-sided tests (overexpression), use FC > 1
    fc_output_grid <- seq(1, max(fc_range[2], 2), length.out = fc_curve_points)
  } else {
    # For two-sided tests, use full range
    fc_output_grid <- seq(fc_range[1], fc_range[2], length.out = fc_curve_points)
  }

  # Use log-spaced points for gene expression evaluation grid
  expr_output_grid <- 10^seq(log10(expr_range[1]), log10(expr_range[2]), length.out = expr_curve_points)

  ############### compute the power for the cells-reads grid ###################
  power_df <- cells_reads_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      power_output = list(
        .compute_underspecified_power_separated(
          # experimental information
          num_total_cells = num_total_cells, library_size = reads_per_cell, MOI = MOI,
          num_targets = num_targets, gRNAs_per_target = gRNAs_per_target,
          non_targeting_gRNAs = non_targeting_gRNAs,
          # analysis information
          multiple_testing_alpha = fdr_target, multiple_testing_method = "BH",
          control_group = control_group, side = side, num_pairs = num_pairs,
          # separated approach information
          fc_expression_df = fc_expression_info$fc_expression_df,
          expression_dispersion_curve = fc_expression_info$expression_dispersion_curve,
          fc_output_grid = fc_output_grid,
          expr_output_grid = expr_output_grid,
          prop_non_null = prop_non_null
        )
      ),
      overall_power = power_output$overall_power,
      power_by_fc = list(power_output$power_by_fc),
      power_by_expr = list(power_output$power_by_expr)
    ) |>
    dplyr::select(-power_output)

  # return the output dataframe
  return(power_df)
}


#' Compute power grid using efficient C++ test statistic computation
#'
#' This function is identical to compute_power_grid_separated() except it uses
#' the C++ implementation compute_distribution_teststat_fixed_es_cpp() instead of
#' the R function compute_test_stat_clean() for test statistic computation.
#'
#' @inheritParams compute_power_grid_separated
#' @return Data frame with power analysis results
#' @export
compute_power_grid_efficient <- function(
    cells_reads_df,
    num_targets = 100,
    gRNAs_per_target = 4,
    non_targeting_gRNAs = 10,
    num_pairs = 1000,
    tpm_threshold = 10,
    fdr_target = 0.05,
    fc_mean = 0.85,
    fc_sd = 0.15,
    prop_non_null = 0.1,
    MOI = 10,
    biological_system = "K562",
    experimental_platform = "10x Chromium v3",
    side = "left",
    control_group = "complement",
    B = 500,
    fc_curve_points = 10,
    expr_curve_points = 10
){

  ############### extract Monte Carlo samples for integration ##################
  set.seed(1)  # Reproducible results
  fc_expression_info <- extract_fc_expression_info(
    fold_change_mean = fc_mean, fold_change_sd = fc_sd,
    biological_system = biological_system, B = B
  )

  # Define systematic output grids
  fc_range <- range(fc_expression_info$fc_expression_df$fold_change)
  expr_range <- range(fc_expression_info$fc_expression_df$relative_expression)

  # Use directional fold change grid based on test sidedness
  if (side == "left") {
    # For left-sided tests (knockdown), use FC < 1
    fc_output_grid <- seq(min(fc_range[1], 0.5), 1, length.out = fc_curve_points)
  } else if (side == "right") {
    # For right-sided tests (overexpression), use FC > 1
    fc_output_grid <- seq(1, max(fc_range[2], 2), length.out = fc_curve_points)
  } else {
    # For two-sided tests, use full range
    fc_output_grid <- seq(fc_range[1], fc_range[2], length.out = fc_curve_points)
  }

  # Use log-spaced points for gene expression evaluation grid
  expr_output_grid <- 10^seq(log10(expr_range[1]), log10(expr_range[2]), length.out = expr_curve_points)

  ############### compute the power for the cells-reads grid ###################
  power_df <- cells_reads_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      power_output = list(
        .compute_underspecified_power_efficient(
          # experimental information
          num_total_cells = num_total_cells, library_size = reads_per_cell, MOI = MOI,
          num_targets = num_targets, gRNAs_per_target = gRNAs_per_target,
          non_targeting_gRNAs = non_targeting_gRNAs,
          # analysis information
          multiple_testing_alpha = fdr_target, multiple_testing_method = "BH",
          control_group = control_group, side = side, num_pairs = num_pairs,
          # separated approach information
          fc_expression_df = fc_expression_info$fc_expression_df,
          expression_dispersion_curve = fc_expression_info$expression_dispersion_curve,
          fc_output_grid = fc_output_grid,
          expr_output_grid = expr_output_grid,
          prop_non_null = prop_non_null
        )
      ),
      overall_power = power_output$overall_power,
      power_by_fc = list(power_output$power_by_fc),
      power_by_expr = list(power_output$power_by_expr)
    ) |>
    dplyr::select(-power_output)

  # return the output dataframe
  return(power_df)
}


#' Internal function for efficient separated power computation using C++
#'
#' @param num_total_cells Total number of cells
#' @param library_size Library size (reads per cell)
#' @param MOI Multiplicity of infection
#' @param num_targets Number of targets
#' @param gRNAs_per_target Number of gRNAs per target
#' @param non_targeting_gRNAs Number of non-targeting gRNAs
#' @param multiple_testing_alpha Alpha level for multiple testing
#' @param multiple_testing_method Multiple testing method
#' @param control_group Control group type
#' @param side Test sidedness
#' @param num_pairs Number of pairs
#' @param fc_expression_df Data frame with fold change and expression info
#' @param expression_dispersion_curve Function for expression-size relationship
#' @param fc_output_grid Grid points for fold change curve
#' @param expr_output_grid Grid points for expression curve
#' @param prop_non_null Proportion of non-null hypotheses
#' @return List with overall power and power curves
.compute_underspecified_power_efficient <- function(
  # experimental information
  num_total_cells, library_size, MOI = 10, num_targets = 100, gRNAs_per_target = 4, non_targeting_gRNAs = 10,
  # analysis information
  multiple_testing_alpha = 0.05, multiple_testing_method = "BH", control_group = "complement", side = "left", num_pairs = 1000,
  # separated approach information
  fc_expression_df, expression_dispersion_curve, fc_output_grid, expr_output_grid, prop_non_null = 0.1){

  ################ compute the treatment and control cells #####################
  num_trt_cells <- gRNAs_per_target * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
  num_cntrl_cells <- round(switch(control_group,
                                  complement = {
                                    num_total_cells - num_trt_cells
                                  },
                                  nt_cells = {
                                    non_targeting_gRNAs * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
                                  }))

  ##############################################################################
  #################### compute Monte Carlo integration ########################

  # Compute test statistics for B Monte Carlo samples (for integration)
  mc_expression_means <- library_size * fc_expression_df$relative_expression

  mc_test_stats <- vector("list", nrow(fc_expression_df))
  for (i in 1:nrow(fc_expression_df)) {
    # Use C++ function instead of R function
    mc_test_stats[[i]] <- compute_distribution_teststat_fixed_es_cpp(
      fold_change = fc_expression_df$fold_change[i],
      expression_mean = mc_expression_means[i],
      expression_size = fc_expression_df$expression_size[i],
      num_trt_cells = num_trt_cells,
      num_cntrl_cells = num_cntrl_cells,
      num_cells = num_trt_cells  # For single gRNA case, all cells have same fold change
    )
  }

  # Extract vectors for Monte Carlo samples
  mc_means <- sapply(mc_test_stats, function(x) x[["mean"]])
  mc_sds <- sapply(mc_test_stats, function(x) x[["sd"]])

  ########################## compute the cutoff ################################
  sig_cutoff <- switch(multiple_testing_method,
                       BH = {
                         compute_underspecified_BH_cutoff_cpp(
                           mean_list = mc_means,
                           sd_list = mc_sds,
                           side = side,
                           multiple_testing_alpha = multiple_testing_alpha,
                           prop_non_null = prop_non_null,
                           num_pairs = num_pairs
                         )
                       })

  ####################### compute overall power ################################
  mc_powers <- rejection_computation_cpp(mean_list = mc_means,
                                         sd_list = mc_sds,
                                         side = side,
                                         cutoff = sig_cutoff)
  overall_power <- mean(mc_powers)

  ##############################################################################
  #################### compute output curves ###################################

  ##################### compute FC curve at specific points ###################
  power_by_fc_list <- vector("list", length(fc_output_grid))
  for (i in 1:length(fc_output_grid)) {
    fc_val <- fc_output_grid[i]

    # Compute power for this FC across all Monte Carlo expression samples
    fc_test_stats <- vector("list", nrow(fc_expression_df))
    for (j in 1:nrow(fc_expression_df)) {
      # Use C++ function instead of R function
      fc_test_stats[[j]] <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = fc_val,  # Use systematic FC grid point
        expression_mean = mc_expression_means[j],  # Use MC expression samples
        expression_size = fc_expression_df$expression_size[j],
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = num_trt_cells  # For single gRNA case
      )
    }

    fc_means <- sapply(fc_test_stats, function(x) x[["mean"]])
    fc_sds <- sapply(fc_test_stats, function(x) x[["sd"]])

    fc_powers <- rejection_computation_cpp(mean_list = fc_means,
                                           sd_list = fc_sds,
                                           side = side,
                                           cutoff = sig_cutoff)
    power_by_fc_list[[i]] <- mean(fc_powers)
  }

  power_by_fc <- data.frame(
    fold_change = fc_output_grid,
    power = unlist(power_by_fc_list)
  )

  ################## compute expression curve at specific points ###############
  power_by_expr_list <- vector("list", length(expr_output_grid))
  for (i in 1:length(expr_output_grid)) {
    expr_val <- expr_output_grid[i]
    expr_mean <- library_size * expr_val

    # Use expression_dispersion_curve to get size parameter
    expr_size <- expression_dispersion_curve(expr_val)

    # Compute power for this expression across all Monte Carlo FC samples
    expr_test_stats <- vector("list", nrow(fc_expression_df))
    for (j in 1:nrow(fc_expression_df)) {
      # Use C++ function instead of R function
      expr_test_stats[[j]] <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = fc_expression_df$fold_change[j],  # Use MC FC samples
        expression_mean = expr_mean,  # Use systematic expression grid point
        expression_size = expr_size,
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = num_trt_cells  # For single gRNA case
      )
    }

    expr_means <- sapply(expr_test_stats, function(x) x[["mean"]])
    expr_sds <- sapply(expr_test_stats, function(x) x[["sd"]])

    expr_powers <- rejection_computation_cpp(mean_list = expr_means,
                                             sd_list = expr_sds,
                                             side = side,
                                             cutoff = sig_cutoff)
    power_by_expr_list[[i]] <- mean(expr_powers)
  }

  power_by_expr <- data.frame(
    relative_expression = expr_output_grid,
    power = unlist(power_by_expr_list)
  )

  # return the output
  return(list(
    overall_power = overall_power,
    power_by_fc = power_by_fc,
    power_by_expr = power_by_expr,
    grid_summary = list(
      mc_samples = nrow(fc_expression_df),
      fc_output_points = length(fc_output_grid),
      expr_output_points = length(expr_output_grid),
      total_computations = nrow(fc_expression_df) * (1 + length(fc_output_grid) + length(expr_output_grid)),
      cutoff = sig_cutoff
    )
  ))
}


#' Internal function for separated power computation
#'
#' @param num_total_cells Total number of cells
#' @param library_size Library size (reads per cell)
#' @param MOI Multiplicity of infection
#' @param num_targets Number of targets
#' @param gRNAs_per_target Number of gRNAs per target
#' @param non_targeting_gRNAs Number of non-targeting gRNAs
#' @param multiple_testing_alpha Alpha level for multiple testing
#' @param multiple_testing_method Multiple testing method
#' @param control_group Control group type
#' @param side Test sidedness
#' @param num_pairs Number of pairs
#' @param fc_expression_df Data frame with fold change and expression info
#' @param expression_dispersion_curve Function for expression-size relationship
#' @param fc_output_grid Grid points for fold change curve
#' @param expr_output_grid Grid points for expression curve
#' @param prop_non_null Proportion of non-null hypotheses
#' @return List with overall power and power curves
.compute_underspecified_power_separated <- function(
  # experimental information
  num_total_cells, library_size, MOI = 10, num_targets = 100, gRNAs_per_target = 4, non_targeting_gRNAs = 10,
  # analysis information
  multiple_testing_alpha = 0.05, multiple_testing_method = "BH", control_group = "complement", side = "left", num_pairs = 1000,
  # separated approach information
  fc_expression_df, expression_dispersion_curve, fc_output_grid, expr_output_grid, prop_non_null = 0.1){

  ################ compute the treatment and control cells #####################
  num_trt_cells <- gRNAs_per_target * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
  num_cntrl_cells <- round(switch(control_group,
                                  complement = {
                                    num_total_cells - num_trt_cells
                                  },
                                  nt_cells = {
                                    non_targeting_gRNAs * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
                                  }))

  ##############################################################################
  #################### compute Monte Carlo integration ########################

  # Compute test statistics for B Monte Carlo samples (for integration)
  mc_expression_means <- library_size * fc_expression_df$relative_expression

  mc_test_stats <- vector("list", nrow(fc_expression_df))
  for (i in 1:nrow(fc_expression_df)) {
    mc_test_stats[[i]] <- compute_test_stat_clean(
      num_trt_cells = num_trt_cells,
      num_cntrl_cells = num_cntrl_cells,
      expression_mean = mc_expression_means[i],
      expression_size = fc_expression_df$expression_size[i],
      fold_change_mean = fc_expression_df$fold_change[i]
    )
  }

  # Extract vectors for Monte Carlo samples
  mc_means <- sapply(mc_test_stats, function(x) x["mean"])
  mc_sds <- sapply(mc_test_stats, function(x) x["sd"])

  ########################## compute the cutoff ################################
  sig_cutoff <- switch(multiple_testing_method,
                       BH = {
                         compute_underspecified_BH_cutoff_cpp(
                           mean_list = mc_means,
                           sd_list = mc_sds,
                           side = side,
                           multiple_testing_alpha = multiple_testing_alpha,
                           prop_non_null = prop_non_null,
                           num_pairs = num_pairs
                         )
                       })

  ####################### compute overall power ################################
  mc_powers <- rejection_computation_cpp(mean_list = mc_means,
                                         sd_list = mc_sds,
                                         side = side,
                                         cutoff = sig_cutoff)
  overall_power <- mean(mc_powers)

  ##############################################################################
  #################### compute output curves ###################################

  ##################### compute FC curve at specific points ###################
  power_by_fc_list <- vector("list", length(fc_output_grid))
  for (i in 1:length(fc_output_grid)) {
    fc_val <- fc_output_grid[i]

    # Compute power for this FC across all Monte Carlo expression samples
    fc_test_stats <- vector("list", nrow(fc_expression_df))
    for (j in 1:nrow(fc_expression_df)) {
      fc_test_stats[[j]] <- compute_test_stat_clean(
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        expression_mean = mc_expression_means[j],  # Use MC expression samples
        expression_size = fc_expression_df$expression_size[j],
        fold_change_mean = fc_val  # Use systematic FC grid point
      )
    }

    fc_means <- sapply(fc_test_stats, function(x) x["mean"])
    fc_sds <- sapply(fc_test_stats, function(x) x["sd"])

    fc_powers <- rejection_computation_cpp(mean_list = fc_means,
                                           sd_list = fc_sds,
                                           side = side,
                                           cutoff = sig_cutoff)
    power_by_fc_list[[i]] <- mean(fc_powers)
  }

  power_by_fc <- data.frame(
    fold_change = fc_output_grid,
    power = unlist(power_by_fc_list)
  )

  ################## compute expression curve at specific points ###############
  power_by_expr_list <- vector("list", length(expr_output_grid))
  for (i in 1:length(expr_output_grid)) {
    expr_val <- expr_output_grid[i]
    expr_mean <- library_size * expr_val

    # Use expression_dispersion_curve to get size parameter
    expr_size <- expression_dispersion_curve(expr_val)

    # Compute power for this expression across all Monte Carlo FC samples
    expr_test_stats <- vector("list", nrow(fc_expression_df))
    for (j in 1:nrow(fc_expression_df)) {
      expr_test_stats[[j]] <- compute_test_stat_clean(
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        expression_mean = expr_mean,  # Use systematic expression grid point
        expression_size = expr_size,
        fold_change_mean = fc_expression_df$fold_change[j]  # Use MC FC samples
      )
    }

    expr_means <- sapply(expr_test_stats, function(x) x["mean"])
    expr_sds <- sapply(expr_test_stats, function(x) x["sd"])

    expr_powers <- rejection_computation_cpp(mean_list = expr_means,
                                             sd_list = expr_sds,
                                             side = side,
                                             cutoff = sig_cutoff)
    power_by_expr_list[[i]] <- mean(expr_powers)
  }

  power_by_expr <- data.frame(
    relative_expression = expr_output_grid,
    power = unlist(power_by_expr_list)
  )

  # return the output
  return(list(
    overall_power = overall_power,
    power_by_fc = power_by_fc,
    power_by_expr = power_by_expr,
    grid_summary = list(
      mc_samples = nrow(fc_expression_df),
      fc_output_points = length(fc_output_grid),
      expr_output_points = length(expr_output_grid),
      total_computations = nrow(fc_expression_df) * (1 + length(fc_output_grid) + length(expr_output_grid)),
      cutoff = sig_cutoff
    )
  ))
}


#' Example usage of the separated Monte Carlo power analysis
#'
#' This function demonstrates how to use the optimized power analysis approach
#' with reasonable default parameters.
#'
#' @param num_cells Vector of cell numbers to test
#' @param reads_per_cell Vector of read depths to test
#' @param B Number of Monte Carlo samples (default: 100)
#' @param curve_points Number of points for power curves (default: 10)
#' @return Data frame with power analysis results
#' @export
#' @examples
#' # Quick power analysis for different experimental conditions
#' result <- example_power_analysis(
#'   num_cells = c(10000, 20000),
#'   reads_per_cell = c(500, 1000)
#' )
#' print(result$overall_power)
example_power_analysis <- function(
  num_cells = c(10000, 20000),
  reads_per_cell = c(500, 1000),
  B = 100,
  curve_points = 10
) {

  # Create experimental design grid
  cells_reads_df <- data.frame(
    num_total_cells = num_cells,
    reads_per_cell = reads_per_cell
  )

  # Run separated power analysis
  result <- compute_power_grid_separated(
    cells_reads_df = cells_reads_df,
    B = B,
    num_targets = 100,
    fc_mean = 0.85,
    fc_sd = 0.15,
    prop_non_null = 0.1,
    fc_curve_points = curve_points,
    expr_curve_points = curve_points
  )

  return(result)
}
