# Core power calculation utilities

#' Compute overall power for power analysis (core utility function)
#'
#' This function computes overall power and BH cutoff for a single experimental design.
#' This is a core utility function used by grid-based power analysis functions.
#'
#' @param num_trt_cells Number of treatment cells
#' @param num_cntrl_cells Number of control cells
#' @param library_size Library size (effective UMIs per cell)
#' @param multiple_testing_alpha Alpha level for multiple testing (default: 0.05)
#' @param multiple_testing_method Multiple testing method (default: "BH")
#' @param side Test sidedness ("left", "right", "both", default: "left")
#' @param fc_expression_df Data frame with fold change and expression info
#' @param prop_non_null Proportion of non-null hypotheses (default: 0.1)
#' @param return_full_results If TRUE, return list with all intermediate results; if FALSE, return only overall power
#' 
#' @return Overall power value (scalar) or list with full results depending on return_full_results
#' 
#' @details
#' This function serves as the core power calculation utility that:
#' \itemize{
#'   \item Computes test statistic distributions for each gene
#'   \item Applies Benjamini-Hochberg multiple testing correction
#'   \item Calculates overall statistical power
#' }
#' 
#' The function delegates to the optimized C++ implementation for performance.
#' 
#' @export
compute_power_plan_overall <- function(
    # experimental information
  num_trt_cells, num_cntrl_cells, library_size,
  # analysis information
  multiple_testing_alpha = 0.05, multiple_testing_method = "BH", side = "left",
  # separated approach information
  fc_expression_df, prop_non_null = 0.1, return_full_results = FALSE){

  # Call the C++ implementation for improved performance
  return(compute_power_plan_overall_cpp(
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    multiple_testing_alpha = multiple_testing_alpha,
    multiple_testing_method = multiple_testing_method,
    side = side,
    prop_non_null = prop_non_null,
    return_full_results = return_full_results
  ))
}

#' Compute power analysis for experimental design grid
#'
#' @description
#' Simplified version of identify_cell_read_range that returns a clean dataframe
#' with experimental design combinations and their corresponding power values.
#'
#' @param fc_expression_df Data frame with fold change and expression information.
#' @param library_parameters List containing UMI_per_cell and variation parameters.
#' @param grid_size Integer. Number of points in each dimension of the grid (default: 10).
#' @param min_power_threshold Numeric. Minimum power threshold for cell range determination (default: 0.01).
#' @param max_power_threshold Numeric. Maximum power threshold for cell range determination (default: 0.8).
#' @param MOI Numeric. Multiplicity of infection for cell allocation calculations (default: 10).
#' @param num_targets Integer. Number of targets for cell allocation calculations (default: 100).
#' @param gRNAs_per_target Integer. Number of gRNAs per target (default: 4).
#' @param non_targeting_gRNAs Integer. Number of non-targeting gRNAs (default: 10).
#' @param control_group String. Control group type: "complement" or "nt_cells" (default: "complement").
#' @param multiple_testing_alpha Numeric. Alpha level for multiple testing (default: 0.05).
#' @param side String. Test sidedness: "left", "right", or "both" (default: "left").
#' @param prop_non_null Numeric. Proportion of non-null hypotheses (default: 0.1).
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{cells_per_target}{Numeric. Number of treatment cells per target}
#'   \item{num_total_cells}{Numeric. Total number of cells in the experiment}
#'   \item{reads_per_cell}{Numeric. Sequencing reads per cell}
#'   \item{library_size}{Numeric. Effective library size (UMIs)}
#'   \item{overall_power}{Numeric. Statistical power for this experimental design}
#' }
#'
#' @details
#' This function simplifies the experimental design process by:
#' \enumerate{
#'   \item Determining optimal reads per cell range using library size curves
#'   \item Determining optimal cell count range based on power thresholds
#'   \item Creating logarithmically-spaced grids for both dimensions
#'   \item Computing power for all combinations
#'   \item Returning a clean dataframe ready for analysis
#' }
#'
#' @export
compute_power_plan_per_grid <- function(
  fc_expression_df,
  library_parameters,
  grid_size = 10,
  min_power_threshold = 0.01,
  max_power_threshold = 0.8,
  MOI = 10,
  num_targets = 100,
  gRNAs_per_target = 4,
  non_targeting_gRNAs = 10,
  control_group = "complement",
  multiple_testing_alpha = 0.05,
  side = "left",
  prop_non_null = 0.1
) {

  # Extract needed data
  UMI_per_cell <- library_parameters$UMI_per_cell
  variation <- library_parameters$variation

  # Step 1: Determine reads per cell range using library size curves
  reads_range <- identify_reads_range_cpp(
    UMI_per_cell = UMI_per_cell,
    variation = variation
  )

  min_reads_per_cell <- reads_range$min_reads_per_cell
  max_reads_per_cell <- reads_range$max_reads_per_cell

  # Step 2: Determine cell count range based on power thresholds
  cell_range <- identify_cell_range_cpp(
    min_reads_per_cell = min_reads_per_cell,
    max_reads_per_cell = max_reads_per_cell,
    fc_expression_df = fc_expression_df,
    UMI_per_cell = UMI_per_cell,
    variation = variation,
    MOI = MOI,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    control_group = control_group,
    multiple_testing_alpha = multiple_testing_alpha,
    side = side,
    prop_non_null = prop_non_null,
    min_power_threshold = min_power_threshold,
    max_power_threshold = max_power_threshold
  )

  min_total_cells <- cell_range$min_cells
  max_total_cells <- cell_range$max_cells

  # Step 3: Create logarithmically-spaced grids
  # Logarithmic spacing for cells (total cells)
  cells_seq <- exp(seq(log(min_total_cells), log(max_total_cells), length.out = grid_size))

  # Logarithmic spacing for reads per cell
  reads_seq <- exp(seq(log(min_reads_per_cell), log(max_reads_per_cell), length.out = grid_size))

  # Step 4: Calculate treatment cell counts and control cell counts
  # Convert total cells to treatment cells using experimental design parameters
  total_gRNAs <- num_targets * gRNAs_per_target + non_targeting_gRNAs
  cells_per_target_seq <- (cells_seq * gRNAs_per_target * MOI) / total_gRNAs

  # Calculate control cells based on control group type
  if (control_group == "complement") {
    # Control cells = total cells - treatment cells
    num_cntrl_cells_seq <- cells_seq - cells_per_target_seq
  } else {
    # nt_cells: control cells = non-targeting gRNA cells
    num_cntrl_cells_seq <- (cells_seq * non_targeting_gRNAs * MOI) / total_gRNAs
  }

  # Step 5: Create all combinations and compute power
  design_grid <- expand.grid(
    cells_idx = 1:grid_size,
    reads_idx = 1:grid_size
  ) |>
    dplyr::mutate(
      cells_per_target = cells_per_target_seq[cells_idx],
      num_total_cells = cells_seq[cells_idx],
      reads_per_cell = reads_seq[reads_idx],
      num_cntrl_cells = num_cntrl_cells_seq[cells_idx]
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      library_size = fit_read_UMI_curve_cpp(
        reads_per_cell = reads_per_cell,
        UMI_per_cell = UMI_per_cell,
        variation = variation
      ),
      overall_power = compute_single_power_cpp(
        num_cells = num_total_cells,
        reads_per_cell = reads_per_cell,
        fc_expression_df = fc_expression_df,
        UMI_per_cell = UMI_per_cell,
        variation = variation,
        MOI = MOI,
        num_targets = num_targets,
        gRNAs_per_target = gRNAs_per_target,
        non_targeting_gRNAs = non_targeting_gRNAs,
        control_group = control_group,
        multiple_testing_alpha = multiple_testing_alpha,
        side = side,
        prop_non_null = prop_non_null
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(cells_per_target, num_total_cells, reads_per_cell, library_size, overall_power)

  return(design_grid)
}

#' Compute power analysis for full parameter grid
#'
#' @description
#' This function integrates compute_power_plan_per_grid() to create a comprehensive
#' power analysis across multiple parameter combinations (TPM thresholds, fold changes).
#'
#' @param tpm_threshold Numeric or character. TPM threshold value or "varying" for auto-selection.
#' @param minimum_fold_change Numeric or character. Minimum fold change value or "varying" for auto-selection.
#' @param MOI Numeric. Multiplicity of infection (default: 10).
#' @param num_targets Integer. Number of targets (default: 100).
#' @param non_targeting_gRNAs Integer. Number of non-targeting gRNAs (default: 10).
#' @param gRNAs_per_target Integer. Number of gRNAs per target (default: 4).
#' @param gRNA_variability Numeric. Standard deviation for gRNA effect variation (default: 0.13).
#' @param control_group String. Control group type (default: "complement").
#' @param side String. Test sidedness (default: "left").
#' @param multiple_testing_alpha Numeric. FDR level (default: 0.05).
#' @param prop_non_null Numeric. Proportion of non-null hypotheses (default: 0.1).
#' @param baseline_expression_stats Data frame. Baseline expression statistics.
#' @param library_parameters List. Library parameters with UMI_per_cell and variation.
#' @param grid_size Integer. Grid size for each dimension (default: 10).
#' @param min_power_threshold Numeric. Minimum power threshold (default: 0.01).
#' @param max_power_threshold Numeric. Maximum power threshold (default: 0.8).
#'
#' @return Data frame with comprehensive power analysis results across parameter combinations.
#'
#' @details
#' This function provides comprehensive power analysis by:
#' \enumerate{
#'   \item Expanding parameter combinations (TPM thresholds, fold changes)
#'   \item Creating fold change expression data for each combination
#'   \item Running compute_power_plan_per_grid() for each parameter set
#'   \item Combining results into a flat dataframe for analysis
#' }
#'
#' @export
compute_power_plan_full_grid <- function(
    # power-determining parameters
    tpm_threshold, minimum_fold_change,
    # experimental parameters
    MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4, gRNA_variability = 0.13,
    # analysis parameters
    control_group = "complement", side = "left", multiple_testing_alpha = 0.05, prop_non_null = 0.1,
    # data inputs
    baseline_expression_stats, library_parameters,
    # grid parameters
    grid_size = 10, min_power_threshold = 0.01, max_power_threshold = 0.8
){

  ####################### construct the tpm_threshold ##########################
  if(tpm_threshold == "varying"){
    tpm_threshold_list <- unname(round(quantile(baseline_expression_stats$relative_expression,
                                                probs = seq(0.1, 0.9, length.out = 5)) * 1e6))
  }else{
    tpm_threshold_list <- tpm_threshold
  }

  ####################### construct the minimum_fold_change ####################
  if(minimum_fold_change == "varying"){
    minimum_fold_change_list <- switch (side,
      both = {c(seq(0.5, 0.95, length.out = 5), seq(1.05, 1.5, length.out = 5))},
      left = {seq(0.5, 0.95, length.out = 10)},
      right = {seq(1.05, 1.5, length.out = 10)}
    )
  }else{
    minimum_fold_change_list <- minimum_fold_change
  }

  # construct parameter grid by expanding the grid
  parameter_grid <- expand.grid(
    minimum_fold_change = minimum_fold_change_list,
    tpm_threshold = tpm_threshold_list
  )

  ###################### Generate complete power analysis grid #####################
  parameter_grid <- parameter_grid |>
    dplyr::rowwise() |>
    dplyr::mutate(
      # Create fold change expression dataframe for this combination
      fc_expression_df = list({
        # Filter expression data by TPM threshold
        filtered_expression_df <- baseline_expression_stats |>
          dplyr::filter(relative_expression >= tpm_threshold / 1e6)

        # Add fold change parameters
        filtered_expression_df |>
          dplyr::rowwise() |>
          dplyr::mutate(
            grna_effects = list(stats::rnorm(n = gRNAs_per_target,
                                            mean = minimum_fold_change,
                                            sd = gRNA_variability) |> pmax(0)),
            avg_fold_change = mean(grna_effects),
            avg_fold_change_sq = mean(grna_effects^2)
          ) |>
          dplyr::select(-grna_effects) |>  # Remove intermediate column
          dplyr::ungroup()
      }),
      # Generate experimental design grid with power calculations
      power_grid = list(
        compute_power_plan_per_grid(
          fc_expression_df = fc_expression_df,
          library_parameters = library_parameters,
          grid_size = grid_size,
          min_power_threshold = min_power_threshold,
          max_power_threshold = max_power_threshold,
          MOI = MOI,
          num_targets = num_targets,
          gRNAs_per_target = gRNAs_per_target,
          non_targeting_gRNAs = non_targeting_gRNAs,
          control_group = control_group,
          multiple_testing_alpha = multiple_testing_alpha,
          side = side,
          prop_non_null = prop_non_null
        )
      )
    ) |>
    dplyr::ungroup()

  # Unnest the power grids to create a flat dataframe
  result <- parameter_grid |>
    dplyr::select(-fc_expression_df) |>  # Remove intermediate data
    tidyr::unnest(power_grid) |>
    dplyr::select(minimum_fold_change, tpm_threshold,
                  cells_per_target, num_total_cells, reads_per_cell,
                  library_size, overall_power)

  return(result)
}