# Power calculation function using optimized approach for Shiny app integration

#' Calculate power grid for app heatmap visualization
#'
#' This function provides power analysis functionality for the Shiny application.
#' It creates a grid of cell/read combinations and computes power for each combination.
#'
#' @param num_targets Number of targets
#' @param gRNAs_per_target Number of gRNAs per target
#' @param non_targeting_gRNAs Number of non-targeting gRNAs
#' @param num_pairs Number of pairs analyzed
#' @param tpm_threshold Minimum TPM threshold
#' @param fdr_target FDR target level
#' @param fc_mean Fold-change mean
#' @param fc_sd Fold-change SD
#' @param prop_non_null Proportion of non-null pairs
#' @param MOI Multiplicity of infection
#' @param biological_system Biological system
#' @param experimental_platform Experimental platform
#' @param side Test sidedness ("left", "right", "both")
#' @param control_group Control group type ("complement" or "nt_cells")
#'
#' @return List with power grid, cell/read sequences, and parameters
#' @export
calculate_power_grid <- function(
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
  control_group = "complement"
) {

  # Create grid for heatmap visualization
  cells_seq <- round(seq(5000, 50000, length.out = 20))
  reads_seq <- round(seq(2000, 50000, length.out = 20))

  # Create cells-reads data frame for power calculation
  cells_reads_df <- expand.grid(
    num_total_cells = cells_seq,
    reads_per_cell = reads_seq
  )

  # Call the optimized efficient Monte Carlo power function
  power_results <- compute_power_grid_efficient(
    cells_reads_df = cells_reads_df,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    num_pairs = num_pairs,
    tpm_threshold = tpm_threshold,
    fdr_target = fdr_target,
    fc_mean = fc_mean,
    fc_sd = fc_sd,
    prop_non_null = prop_non_null,
    MOI = MOI,
    biological_system = biological_system,
    experimental_platform = experimental_platform,
    side = side,
    control_group = control_group,
    B = 100,  # Monte Carlo samples for good accuracy vs speed balance
    fc_curve_points = 10,  # Sufficient resolution for curves
    expr_curve_points = 10
  )

  # Transform to expected format for heatmap
  power_grid <- data.frame(
    cells = power_results$num_total_cells,
    reads = power_results$reads_per_cell,
    power = power_results$overall_power
  )

  # Calculate expected discoveries
  expected_discoveries <- sum(power_grid$power) * prop_non_null

  # Return structure compatible with Shiny app
  list(
    power_grid = power_grid,
    cells_seq = cells_seq,
    reads_seq = reads_seq,
    expected_discoveries = expected_discoveries,
    parameters = list(
      num_targets = num_targets,
      gRNAs_per_target = gRNAs_per_target,
      non_targeting_gRNAs = non_targeting_gRNAs,
      num_pairs = num_pairs,
      tpm_threshold = tpm_threshold,
      fdr_target = fdr_target,
      fc_mean = fc_mean,
      fc_sd = fc_sd,
      prop_non_null = prop_non_null,
      MOI = MOI,
      biological_system = biological_system,
      experimental_platform = experimental_platform
    ),
    # Additional data for potential future use
    power_curves = list(
      fc_curves = power_results$power_by_fc,
      expr_curves = power_results$power_by_expr
    )
  )
}





#' Compute power grid using efficient C++ test statistic computation
#'
#' This function uses C++ implementations for computational efficiency in power analysis
#' for perturb-seq experiments across different experimental conditions.
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

  ########################## compute the library size ##########################
  read_UMI_info <- extract_library_info(biological_system = biological_system)
  UMI_per_cell <- read_UMI_info$UMI_per_cell
  variation <- read_UMI_info$variation

  # compute the library size
  cells_reads_df <- cells_reads_df |>
    dplyr::mutate(
      library_size = fit_read_UMI_curve(reads_per_cell = reads_per_cell, UMI_per_cell = !!UMI_per_cell, variation = !!variation)
    )

  ############### compute the power for the cells-reads grid ###################
  power_df <- cells_reads_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      power_output = list(
        .compute_power_plan_efficient(
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



#' Internal function for efficient separated power computation using C++ Monte Carlo
#'
#' This function replaces the R Monte Carlo for loop with C++ implementation
#' for improved performance.
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
.compute_power_plan_efficient <- function(
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

  # Use C++ function to compute test statistics for all Monte Carlo samples
  mc_results <- compute_monte_carlo_teststat_cpp(
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells
  )

  # Extract vectors for Monte Carlo samples
  mc_means <- mc_results$means
  mc_sds <- mc_results$sds

  ########################## compute the cutoff ################################
  sig_cutoff <- switch(multiple_testing_method,
                       BH = {
                         compute_BH_plan(
                           mean_list = mc_means,
                           sd_list = mc_sds,
                           side = side,
                           multiple_testing_alpha = multiple_testing_alpha,
                           prop_non_null = prop_non_null
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

  # Compute fold change curve using C++ function for performance
  power_by_fc <- compute_fc_curve_cpp(
    fc_output_grid = fc_output_grid,
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    side = side,
    cutoff = sig_cutoff
  )

  # Compute expression curve using C++ function for performance
  power_by_expr <- compute_expression_curve_cpp(
    expr_output_grid = expr_output_grid,
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    expression_dispersion_curve = expression_dispersion_curve,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    side = side,
    cutoff = sig_cutoff
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
    reads_per_cell = c(1e4, 2e4),
    B = 100,
    curve_points = 10
) {

  # Create experimental design grid
  cells_reads_df <- data.frame(
    num_total_cells = num_cells,
    reads_per_cell = reads_per_cell
  )

  # Run separated power analysis
  result <- compute_power_grid_efficient(
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
