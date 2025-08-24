# Core power calculation utilities

# Suppress R CMD check notes for NSE (non-standard evaluation) variables
utils::globalVariables(c("total_cost", "library_cost", "sequencing_cost", ".data",
                         "minimum_cost", "min_cells", "max_cells", "min_reads",
                         "max_reads", "cost_precision", "cost_of_interest", "cost_grid"))

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
#' @param cells_per_target Numeric, numeric vector, or character. Number of cells per target, custom sequence, or "varying" for auto-generated grid.
#' @param reads_per_cell Numeric, numeric vector, or character. Reads per cell, custom sequence, or "varying" for auto-generated grid.
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
#' @param mapping_efficiency Numeric. Mapping efficiency for raw reads to usable reads (default: 0.72).
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{cells_per_target}{Numeric. Number of treatment cells per target}
#'   \item{reads_per_cell}{Numeric. Sequencing reads per cell}
#'   \item{library_size}{Numeric. Effective library size (UMIs)}
#'   \item{overall_power}{Numeric. Statistical power for this experimental design}
#'   \item{num_captured_cells}{Numeric. Number of captured cells}
#'   \item{raw_reads_per_cell}{Numeric. Raw reads per cell before mapping}
#' }
#'
#' @details
#' This function simplifies the experimental design process by:
#' \enumerate{
#'   \item Determining parameter sequences based on input type:
#'     - Numeric value: Single fixed value (length = 1)
#'     - Numeric vector: Custom sequence (length = vector length)
#'     - "varying": Auto-generated sequence using optimal ranges (length = grid_size)
#'   \item Creating dynamic grids based on sequence lengths:
#'     - Grid size = length(cells_seq) × length(reads_seq)
#'     - Examples: 1×1 (both fixed), 5×1 (custom cells, fixed reads), 10×10 (both varying)
#'   \item Computing power for all parameter combinations
#'   \item Returning a clean dataframe ready for analysis
#' }
#'
#' @export
compute_power_plan_per_grid <- function(
  cells_per_target,
  reads_per_cell,
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
  prop_non_null = 0.1,
  mapping_efficiency = 0.72
) {

  # Extract needed data
  UMI_per_cell <- library_parameters$UMI_per_cell
  variation <- library_parameters$variation

  # Step 1: Determine reads per cell sequence
  if (is.numeric(reads_per_cell)) {
    # Both single value and custom sequence
    reads_seq <- reads_per_cell
  } else if (reads_per_cell == "varying") {
    # Auto-generated sequence using library size curves
    reads_range <- identify_reads_range_cpp(
      UMI_per_cell = UMI_per_cell,
      variation = variation
    )
    min_reads_per_cell <- reads_range$min_reads_per_cell
    max_reads_per_cell <- reads_range$max_reads_per_cell
    reads_seq <- exp(seq(log(min_reads_per_cell), log(max_reads_per_cell), length.out = grid_size))
  }

  # Step 2: Determine cell sequence
  if (is.numeric(cells_per_target)) {
    # Both single value and custom sequence - convert to total cells
    total_gRNAs <- num_targets * gRNAs_per_target + non_targeting_gRNAs
    cells_seq <- (cells_per_target * total_gRNAs) / (gRNAs_per_target * MOI)
  } else if (cells_per_target == "varying") {
    # Auto-generated sequence using power thresholds
    # Need reads range for cell range calculation
    if (is.numeric(reads_per_cell)) {
      min_reads_per_cell <- min(reads_per_cell)
      max_reads_per_cell <- max(reads_per_cell)
    } else {
      # Already calculated above when reads_per_cell == "varying"
    }

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
    cells_seq <- exp(seq(log(min_total_cells), log(max_total_cells), length.out = grid_size))
  }

  # Step 3: Determine grid sizes based on sequence lengths
  cells_grid_size <- length(cells_seq)
  reads_grid_size <- length(reads_seq)

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

  # Step 5: Create dynamic grid based on which parameters are varying
  design_grid <- expand.grid(
    cells_idx = 1:cells_grid_size,
    reads_idx = 1:reads_grid_size
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
    dplyr::select(cells_per_target, num_total_cells,
                  reads_per_cell, library_size, overall_power) |>
    dplyr::mutate(num_captured_cells = num_total_cells,
                  raw_reads_per_cell = reads_per_cell / mapping_efficiency) |>
    dplyr::select(-num_total_cells)

  return(design_grid)
}

#' Compute power analysis for full parameter grid
#'
#' @description
#' This function integrates compute_power_plan_per_grid() to create a comprehensive
#' power analysis across multiple parameter combinations (TPM thresholds, fold changes).
#'
#' @param TPM_threshold Numeric, numeric vector, or character. TPM threshold value, custom sequence, or "varying" for auto-selection.
#' @param minimum_fold_change Numeric, numeric vector, or character. Minimum fold change value, custom sequence, or "varying" for auto-selection.
#' @param cells_per_target Numeric, numeric vector, or character. Number of cells per target, custom sequence, or "varying" for auto-generated grid.
#' @param reads_per_cell Numeric, numeric vector, or character. Reads per cell, custom sequence, or "varying" for auto-generated grid.
#' @param MOI Numeric. Multiplicity of infection (default: 10).
#' @param num_targets Integer. Number of targets (default: 100).
#' @param non_targeting_gRNAs Integer. Number of non-targeting gRNAs (default: 10).
#' @param gRNAs_per_target Integer. Number of gRNAs per target (default: 4).
#' @param gRNA_variability Numeric. Standard deviation for gRNA effect variation (default: 0.13).
#' @param assay String. Assay type: "Perturb-seq" or "TAP-seq" (default: "Perturb-seq").
#' @param control_group String. Control group type (default: "complement").
#' @param side String. Test sidedness (default: "left").
#' @param multiple_testing_alpha Numeric. FDR level (default: 0.05).
#' @param prop_non_null Numeric. Proportion of non-null hypotheses (default: 0.1).
#' @param baseline_expression_stats Data frame. Baseline expression statistics.
#' @param library_parameters List. Library parameters with UMI_per_cell and variation.
#' @param grid_size Integer. Grid size for each dimension (default: 10).
#' @param min_power_threshold Numeric. Minimum power threshold (default: 0.01).
#' @param max_power_threshold Numeric. Maximum power threshold (default: 0.8).
#' @param mapping_efficiency Numeric. Mapping efficiency for raw reads to usable reads (default: 0.72).
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
#' @importFrom stats quantile
#' @export
compute_power_plan_full_grid <- function(
    # power-determining parameters
    TPM_threshold, minimum_fold_change, cells_per_target, reads_per_cell,
    # experimental parameters
    MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4, gRNA_variability = 0.13, assay = "Perturb-seq",
    # analysis parameters
    control_group = "complement", side = "left", multiple_testing_alpha = 0.05, prop_non_null = 0.1,
    # data inputs
    baseline_expression_stats, library_parameters,
    # grid parameters
    grid_size = 10, min_power_threshold = 0.01, max_power_threshold = 0.8,
    # efficiency of library preparation and sequencing platform
    mapping_efficiency = 0.72
){

  ####################### Input validation ####################################
  input_check_compute_power_plan_full_grid(
    TPM_threshold = TPM_threshold, minimum_fold_change = minimum_fold_change,
    cells_per_target = cells_per_target, reads_per_cell = reads_per_cell,
    MOI = MOI, num_targets = num_targets, non_targeting_gRNAs = non_targeting_gRNAs,
    gRNAs_per_target = gRNAs_per_target, gRNA_variability = gRNA_variability,
    control_group = control_group, side = side, multiple_testing_alpha = multiple_testing_alpha,
    prop_non_null = prop_non_null, baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters, grid_size = grid_size,
    min_power_threshold = min_power_threshold, max_power_threshold = max_power_threshold,
    mapping_efficiency = mapping_efficiency
  )

  ####################### construct the TPM_threshold ##########################
  if (is.numeric(TPM_threshold)) {
    TPM_threshold_list <- TPM_threshold  # Works for both single values and vectors
  } else if (TPM_threshold == "varying") {
    TPM_threshold_list <- unname(round(quantile(baseline_expression_stats$relative_expression,
                                                probs = seq(0.1, 0.9, length.out = 5)) * 1e6))
  }

  ####################### construct the minimum_fold_change ####################
  if (is.numeric(minimum_fold_change)) {
    minimum_fold_change_list <- minimum_fold_change  # Works for both single values and vectors
  } else if (minimum_fold_change == "varying") {
    minimum_fold_change_list <- switch (side,
      both = {c(seq(0.5, 0.95, length.out = 5), seq(1.05, 1.5, length.out = 5))},
      left = {seq(0.5, 0.95, length.out = 10)},
      right = {seq(1.05, 1.5, length.out = 10)}
    )
  }

  # construct parameter grid by expanding the grid
  parameter_grid <- expand.grid(
    minimum_fold_change = minimum_fold_change_list,
    TPM_threshold = TPM_threshold_list
  )

  ###################### Generate complete power analysis grid #####################
  parameter_grid <- parameter_grid |>
    dplyr::rowwise() |>
    dplyr::mutate(
      # Create fold change expression dataframe for this combination
      fc_expression_df = list({
        # Filter expression data by TPM threshold
        if (assay == "perturb-seq") {
          filtered_expression_df <- baseline_expression_stats |>
            dplyr::filter(relative_expression >= TPM_threshold / 1e6)
        } else if (assay == "TAP-seq") {
          filtered_expression_df <- baseline_expression_stats |>
            dplyr::filter(relative_expression_main >= TPM_threshold / 1e6)
        }

        # Add fold change parameters
        filtered_expression_df |>
          dplyr::rowwise() |>
          dplyr::mutate(
            grna_effects = list(stats::rnorm(n = gRNAs_per_target,
                                            mean = minimum_fold_change,
                                            sd = gRNA_variability) |> pmax(.Machine$double.eps)),
            avg_fold_change = mean(grna_effects),
            avg_fold_change_sq = mean(grna_effects^2)
          ) |>
          dplyr::select(-grna_effects) |>  # Remove intermediate column
          dplyr::ungroup()
      }),
      # Generate experimental design grid with power calculations
      power_grid = list(
        compute_power_plan_per_grid(
          cells_per_target = cells_per_target,
          reads_per_cell = reads_per_cell,
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
          prop_non_null = prop_non_null,
          mapping_efficiency = mapping_efficiency
        )
      )
    ) |>
    dplyr::ungroup()

  # Unnest the power grids to create a flat dataframe
  result <- parameter_grid |>
    dplyr::select(-fc_expression_df) |>  # Remove intermediate data
    tidyr::unnest(power_grid) |>
    dplyr::select(minimum_fold_change, TPM_threshold,
                  cells_per_target, num_captured_cells, raw_reads_per_cell,
                  library_size, overall_power)

  return(result)
}


#' Cost-Power Optimization for Perturb-seq Experimental Design
#'
#' Performs cost-constrained optimization to find the minimal parameter value that achieves
#' target statistical power within a specified budget for perturb-seq experiments.
#'
#' @param minimizing_variable Character. The parameter to minimize during optimization.
#'   Options: "TPM_threshold", "minimum_fold_change", "cells_per_target", "reads_per_cell", or "cost". Default: "TPM_threshold".
#' @param fixed_variable List. Fixed values for other analysis parameters. Can include:
#'   \itemize{
#'     \item \code{minimum_fold_change}: Fixed fold change threshold (when optimizing TPM_threshold, or required for cells_per_target/reads_per_cell/cost)
#'     \item \code{TPM_threshold}: Fixed TPM threshold (when optimizing minimum_fold_change, or required for cells_per_target/reads_per_cell/cost)
#'     \item \code{cells_per_target}: Fixed cells per target (otherwise uses "varying")
#'     \item \code{reads_per_cell}: Fixed reads per cell (otherwise uses "varying")
#'   }
#' @param MOI Numeric. Multiplicity of infection (default: 10).
#' @param num_targets Integer. Number of targets (default: 100).
#' @param non_targeting_gRNAs Integer. Number of non-targeting gRNAs (default: 10).
#' @param gRNAs_per_target Integer. Number of gRNAs per target (default: 4).
#' @param gRNA_variability Numeric. gRNA variability parameter (default: 0.13).
#' @param control_group Character. Control group type: "complement" or "non_targeting" (default: "complement").
#' @param side Character. Test side: "left", "right", or "both" (default: "left").
#' @param multiple_testing_alpha Numeric. Multiple testing significance level (default: 0.05).
#' @param prop_non_null Numeric. Proportion of non-null hypotheses (default: 0.1).
#' @param baseline_expression_stats Data frame. Baseline expression statistics with columns:
#'   \code{response_id}, \code{relative_expression}, \code{expression_size}.
#' @param library_parameters List. Library parameters containing \code{UMI_per_cell} and \code{variation}.
#' @param grid_size Integer. Grid size for coarse search (default: 20).
#' @param power_target Numeric. Target statistical power (default: 0.8).
#' @param power_precision Numeric. Acceptable precision around power target (default: 0.01).
#' @param cost_per_captured_cell Numeric. Cost per captured cell in dollars (default: 0.086).
#' @param cost_per_million_reads Numeric. Cost per million sequencing reads in dollars (default: 0.374).
#' @param cost_constraint Numeric. Maximum budget constraint in dollars (default: 1e4).
#'   Set to \code{NULL} to disable cost constraints.
#' @param mapping_efficiency Numeric. Sequencing mapping efficiency (default: 0.72).
#'
#' @return Data frame with power analysis results including:
#'   \itemize{
#'     \item Analysis parameters (TPM_threshold, minimum_fold_change, etc.)
#'     \item Experimental design (cells_per_target, num_captured_cells, raw_reads_per_cell)
#'     \item Power metrics (overall_power)
#'     \item Cost breakdown (library_cost, sequencing_cost, total_cost) - when cost_constraint is specified
#'     \item Power threshold indicator (meets_threshold) - when cost_constraint is NULL
#'   }
#'
#' @details
#' This function implements a two-stage optimization algorithm:
#'
#' \strong{Stage 1: Coarse Grid Search}
#' \enumerate{
#'   \item Creates parameter grid for the minimizing variable:
#'     \itemize{
#'       \item TPM_threshold: 20 log-spaced values from 1 to 1000 TPM
#'       \item minimum_fold_change: 20 values based on test side (left: 0.01-1.0, right: 1.0-10.0, both: combined)
#'     }
#'   \item Runs power analysis across experimental design space
#'   \item Applies dual filtering:
#'     \itemize{
#'       \item Power filter: \code{power_target ± power_precision}
#'       \item Cost filter: \code{total_cost ≤ budget_precision × cost_constraint}
#'     }
#'   \item Identifies minimum parameter value meeting both constraints
#' }
#'
#' \strong{Stage 2: Fine Grid Search}
#' \enumerate{
#'   \item Re-runs power analysis with optimal parameter and higher resolution (100 grid points)
#'   \item Combines coarse and fine search results for comprehensive output
#' }
#'
#' \strong{Cost Model:}
#' \deqn{Total Cost = Library Cost + Sequencing Cost}
#' \deqn{Library Cost = cost\_per\_captured\_cell \times num\_captured\_cells}
#' \deqn{Sequencing Cost = cost\_per\_million\_reads \times raw\_reads\_per\_cell \times num\_captured\_cells / 10^6}
#'
#' @examples
#' \dontrun{
#' # Load pilot data
#' pilot_data <- get_pilot_data_from_package("K562")
#'
#' # Optimize TPM threshold with fixed fold change
#' result1 <- cost_power_computation(
#'   minimizing_variable = "TPM_threshold",
#'   fixed_variable = list(minimum_fold_change = 0.8),
#'   baseline_expression_stats = pilot_data$baseline_expression_stats,
#'   library_parameters = pilot_data$library_parameters,
#'   power_target = 0.8,
#'   cost_constraint = 15000,
#'   cost_precision = 0.9
#' )
#'
#' # Optimize fold change with fixed TPM threshold
#' result2 <- cost_power_computation(
#'   minimizing_variable = "minimum_fold_change",
#'   fixed_variable = list(TPM_threshold = 10),
#'   baseline_expression_stats = pilot_data$baseline_expression_stats,
#'   library_parameters = pilot_data$library_parameters,
#'   power_target = 0.9,
#'   cost_constraint = NULL  # No cost constraint
#' )
#' }
#'
#' @seealso
#' \code{\link{compute_power_plan_full_grid}} for the underlying power analysis
#' \code{\link{find_optimal_cost_design}} for cost optimization

#' Cost-Constrained Power Analysis for Perturb-seq Experiments
#'
#' Performs comprehensive power analysis across experimental design space with optional
#' cost constraints for perturb-seq experiments. Computes power across parameter grids
#' and applies filtering based on power targets and budget constraints.
#'
#' @param minimizing_variable Character. The parameter to vary during analysis.
#'   Options: "TPM_threshold" or "minimum_fold_change". Default: "TPM_threshold".
#' @param fixed_variable List. Fixed values for other analysis parameters. Can include:
#'   \itemize{
#'     \item \code{minimum_fold_change}: Fixed fold change threshold (when varying TPM_threshold)
#'     \item \code{TPM_threshold}: Fixed TPM threshold (when varying minimum_fold_change)
#'     \item \code{cells_per_target}: Fixed cells per target (otherwise uses "varying")
#'     \item \code{reads_per_cell}: Fixed reads per cell (otherwise uses "varying")
#'   }
#' @param MOI Numeric. Multiplicity of infection (default: 10).
#' @param num_targets Integer. Number of targets (default: 100).
#' @param non_targeting_gRNAs Integer. Number of non-targeting gRNAs (default: 10).
#' @param gRNAs_per_target Integer. Number of gRNAs per target (default: 4).
#' @param gRNA_variability Numeric. gRNA variability parameter (default: 0.13).
#' @param assay String. Assay type: "perturb-seq" or "TAP-seq" (default: "perturb-seq").
#' @param control_group Character. Control group type: "complement" or "non_targeting" (default: "complement").
#' @param side Character. Test side: "left", "right", or "both" (default: "left").
#' @param multiple_testing_alpha Numeric. Multiple testing significance level (default: 0.05).
#' @param prop_non_null Numeric. Proportion of non-null hypotheses (default: 0.1).
#' @param baseline_expression_stats Data frame. Baseline expression statistics with columns:
#'   \code{response_id}, \code{relative_expression}, \code{expression_size}.
#' @param library_parameters List. Library parameters containing \code{UMI_per_cell} and \code{variation}.
#' @param grid_size Integer. Grid size for parameter search (default: 20).
#' @param power_target Numeric. Target statistical power (default: 0.8).
#' @param power_precision Numeric. Acceptable precision around power target (default: 0.01).
#' @param min_power Numeric. Minimum power threshold for grid search (default: 0.05).
#' @param max_power Numeric. Maximum power threshold for grid search (default: 0.95).
#' @param cost_precision Numeric. Cost utilization factor (default: 0.9).
#'   Filters designs with total cost ≤ cost_precision × cost_constraint.
#' @param cost_per_captured_cell Numeric. Cost per captured cell in dollars (default: 0.086).
#' @param cost_per_million_reads Numeric. Cost per million sequencing reads in dollars (default: 0.374).
#' @param cost_constraint Numeric. Maximum budget constraint in dollars (default: NULL).
#'   Set to \code{NULL} to disable cost constraints.
#' @param mapping_efficiency Numeric. Sequencing mapping efficiency (default: 0.72).
#'
#' @return Data frame with power analysis results including:
#'   \itemize{
#'     \item Analysis parameters (TPM_threshold, minimum_fold_change, etc.)
#'     \item Experimental design (cells_per_target, num_captured_cells, raw_reads_per_cell)
#'     \item Power metrics (overall_power)
#'     \item Cost breakdown (library_cost, sequencing_cost, total_cost)
#'   }
#'
#' @details
#' This function performs comprehensive power analysis by:
#' \enumerate{
#'   \item Setting up parameter grids based on the minimizing variable
#'   \item Computing power across experimental design space
#'   \item Calculating costs for each design
#'   \item Applying validation checks via \code{check_power_results()}
#' }
#'
#' Parameter grid generation:
#' \itemize{
#'   \item \code{TPM_threshold}: Uses quantiles of baseline expression (10th to 99th percentile)
#'   \item \code{minimum_fold_change}: Uses ranges based on test side (left: 0.5-0.9, right: 1-10, both: combined)
#' }
#'
#' @examples
#' \dontrun{
#' # Load pilot data
#' pilot_data <- get_pilot_data_from_package("K562")
#'
#' # Compute power across TPM_threshold range
#' result1 <- cost_power_computation(
#'   minimizing_variable = "TPM_threshold",
#'   fixed_variable = list(minimum_fold_change = 0.8),
#'   baseline_expression_stats = pilot_data$baseline_expression_stats,
#'   library_parameters = pilot_data$library_parameters,
#'   power_target = 0.8,
#'   cost_constraint = 15000
#' )
#'
#' # Compute power across fold change range
#' result2 <- cost_power_computation(
#'   minimizing_variable = "minimum_fold_change",
#'   fixed_variable = list(TPM_threshold = 50),
#'   baseline_expression_stats = pilot_data$baseline_expression_stats,
#'   library_parameters = pilot_data$library_parameters,
#'   power_target = 0.8,
#'   cost_constraint = NULL
#' )
#'
#' # Optimize cost across all experimental designs
#' result3 <- cost_power_computation(
#'   minimizing_variable = "cost",
#'   fixed_variable = list(TPM_threshold = 50, minimum_fold_change = 0.8),
#'   baseline_expression_stats = pilot_data$baseline_expression_stats,
#'   library_parameters = pilot_data$library_parameters,
#'   power_target = 0.8,
#'   cost_constraint = NULL
#' )
#'
#' # Optimize cells per target with fixed detection parameters
#' result4 <- cost_power_computation(
#'   minimizing_variable = "cells_per_target",
#'   fixed_variable = list(TPM_threshold = 50, minimum_fold_change = 0.8),
#'   baseline_expression_stats = pilot_data$baseline_expression_stats,
#'   library_parameters = pilot_data$library_parameters,
#'   power_target = 0.8,
#'   cost_constraint = 10000
#' )
#' }
#'
#' @export
cost_power_computation <- function(minimizing_variable = "TPM_threshold", fixed_variable = list(minimum_fold_change = 0.8),
                                   # experimental parameters
                                   MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4, gRNA_variability = 0.13, assay = "perturb-seq",
                                   # analysis parameters
                                   control_group = "complement", side = "left", multiple_testing_alpha = 0.05, prop_non_null = 0.1,
                                   # data inputs
                                   baseline_expression_stats, library_parameters,
                                   # grid parameters for power
                                   grid_size = 20, power_target = 0.8, power_precision = 0.01, min_power = 0.05, max_power = 0.95,
                                   # grid parameter for budget
                                   cost_precision = 0.9,
                                   # cost parameters
                                   cost_per_captured_cell = 0.086, cost_per_million_reads = 0.374, cost_constraint = NULL,
                                   # efficiency of library preparation and sequencing platform
                                   mapping_efficiency = 0.72){

  ####################### Input validation ####################################
  input_check_cost_power_computation(
    minimizing_variable = minimizing_variable, fixed_variable = fixed_variable,
    MOI = MOI, num_targets = num_targets, non_targeting_gRNAs = non_targeting_gRNAs,
    gRNAs_per_target = gRNAs_per_target, gRNA_variability = gRNA_variability,
    control_group = control_group, side = side, multiple_testing_alpha = multiple_testing_alpha,
    prop_non_null = prop_non_null, baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters, grid_size = grid_size, power_target = power_target,
    power_precision = power_precision, min_power = min_power, max_power = max_power,
    cost_precision = cost_precision, cost_per_captured_cell = cost_per_captured_cell,
    cost_per_million_reads = cost_per_million_reads, cost_constraint = cost_constraint,
    mapping_efficiency = mapping_efficiency
  )

  # set the seed for computation
  set.seed(1)

  ################# Power-determining parameters setup #########################
  # specify the parameter grid for TPM threshold and minimum_fold_change
  switch (minimizing_variable,
    TPM_threshold = {
      # specify the TPM_threshold based on if TPM_threshold is given in the fixed_variable list or not
      if(is.null(fixed_variable$TPM_threshold)){
        TPM_threshold <- unname(quantile(baseline_expression_stats$relative_expression, probs = seq(0.1, .99, length.out = 20))) * 1e6
      }else{
        TPM_threshold <- fixed_variable$TPM_threshold
      }
      minimum_fold_change <- fixed_variable$minimum_fold_change
    },
    minimum_fold_change = {
      # specify the minimum_fold_change based on if minimum_fold_change is given in the fixed_variable list or not
      if(is.null(fixed_variable$minimum_fold_change)){
        minimum_fold_change <- switch(side,
                                      left = {seq(0.5, 0.9, length.out = 20)},
                                      right = {seq(1.1, 1.5, length.out = 20)},
                                      both = {c(seq(0.5, 0.9, length.out = 10), seq(1.1, 1.5, length.out = 10))})
      }else{
        minimum_fold_change <- fixed_variable$minimum_fold_change
      }
      TPM_threshold <- fixed_variable$TPM_threshold
    },
    cells_per_target = ,
    reads_per_cell = ,
    cost = {
      # When minimizing over cells_per_target or reads_per_cell, get fixed values from fixed_variable
      TPM_threshold <- fixed_variable$TPM_threshold
      minimum_fold_change <- fixed_variable$minimum_fold_change
    }
  )

  # specify cells per target and reads per cell
  cells_per_target <- ifelse(is.null(fixed_variable$cells_per_target), "varying", fixed_variable$cells_per_target)
  reads_per_cell <- ifelse(is.null(fixed_variable$reads_per_cell), "varying", fixed_variable$reads_per_cell)

  ########################## Perform grid power search #########################
  # perform power calculation
  cost_power_df <- compute_power_plan_full_grid(
    # power-determining parameters
    TPM_threshold = TPM_threshold, minimum_fold_change = minimum_fold_change, cells_per_target = cells_per_target, reads_per_cell = reads_per_cell,
    # experimental parameters
    MOI = MOI, num_targets = num_targets, non_targeting_gRNAs = non_targeting_gRNAs, gRNAs_per_target = gRNAs_per_target, gRNA_variability = gRNA_variability, assay = assay,
    # analysis parameters
    control_group = control_group, side = side, multiple_testing_alpha = multiple_testing_alpha, prop_non_null = prop_non_null,
    # data inputs
    baseline_expression_stats = baseline_expression_stats, library_parameters = library_parameters,
    # grid parameters
    grid_size = grid_size, min_power_threshold = max(power_target - 0.3, min_power), max_power_threshold = min(power_target + 0.3, max_power),
    # efficiency of library preparation and sequencing platform
    mapping_efficiency = mapping_efficiency
  ) |>
    # compute the cost
    dplyr::mutate(
      library_cost = cost_per_captured_cell * num_captured_cells,
      sequencing_cost = cost_per_million_reads * raw_reads_per_cell * num_captured_cells / 1e6,
      total_cost = library_cost + sequencing_cost
    )

  ########################## Check the computation #############################
  final_results <- check_power_results(power_df = cost_power_df,
                                        cost_constraint = cost_constraint, cost_precision = cost_precision,
                                        power_target = power_target, power_precision = power_precision)

  # return the final results
  return(final_results)
}

#' Validate Power Analysis Results
#'
#' Validates that power analysis results would have sufficient data after applying
#' cost and power filters. Returns the original unfiltered data if validation passes.
#'
#' @param power_df Data frame. Power analysis results with columns including
#'   \code{overall_power}, \code{total_cost}.
#' @param cost_constraint Numeric. Maximum budget constraint in dollars.
#'   Set to \code{NULL} to disable cost filtering.
#' @param cost_precision Numeric. Cost utilization factor.
#'   Filters designs with total cost ≤ cost_precision × cost_constraint.
#' @param power_target Numeric. Target statistical power.
#' @param power_precision Numeric. Acceptable precision around power target.
#'   Filters designs with power ≥ power_target - power_precision.
#'
#' @return Data frame. The original unfiltered power analysis results (if validation passes).
#'
#' @details
#' This function validates filtering viability without actually filtering:
#' \enumerate{
#'   \item Cost check: Verifies that \code{total_cost ≤ cost_precision × cost_constraint}
#'     would leave at least one row (only when cost_constraint is not NULL)
#'   \item Power check: Verifies that \code{overall_power ≥ power_target - power_precision}
#'     would leave at least one row
#' }
#'
#' Throws informative errors when validation checks fail, otherwise returns the original data.
#'
#' @examples
#' \dontrun{
#' # Validate filtering viability and return original data
#' validated_results <- check_power_results(
#'   power_df = power_analysis_results,
#'   cost_constraint = 15000,
#'   cost_precision = 0.9,
#'   power_target = 0.8,
#'   power_precision = 0.01
#' )
#' }
#'
check_power_results <- function(power_df,
                                cost_constraint, cost_precision,
                                power_target, power_precision){
  ####################### check cost filtering viability #########################
  if(!is.null(cost_constraint)){
    # check if cost filter would leave any results
    cost_check_df <- power_df |> dplyr::filter(total_cost <= cost_precision * cost_constraint)

    # Return error message when filtering fails!
    if(nrow(cost_check_df) == 0){
      stop("The cost optimization failed! Consider increasing cost budget!")
    }
  }

  ####################### check power filtering viability ##########################
  # check if power filter would leave any results
  power_check_df <- power_df |> dplyr::filter(overall_power >= power_target - power_precision)

  # Return the error if the filtering fails!
  if(nrow(power_check_df) == 0){
    stop("The power optimization failed! Consider adjusting power target!")
  }

  # Return the original unfiltered data
  return(power_df)
}

#' Find Optimal Cost-Efficient Experimental Designs
#'
#' @description
#' Identifies cost-optimal experimental designs that achieve target statistical power
#' within specified precision bounds. This function processes cost-power analysis results
#' to find minimal-cost designs for each parameter level and generates detailed cost
#' grids for design optimization.
#'
#' @param cost_power_df Data frame. Output from \code{\link{cost_power_computation}}
#'   containing power analysis results with cost calculations. Must include columns:
#'   \code{overall_power}, \code{total_cost}, \code{cells_per_target},
#'   \code{raw_reads_per_cell}, plus the specified minimizing variable (except for cost optimization).
#' @param minimizing_variable Character. The parameter being optimized. Must be one of:
#'   \itemize{
#'     \item "TPM_threshold": TPM expression threshold optimization
#'     \item "minimum_fold_change": Minimum fold change threshold optimization
#'     \item "cost": Total cost optimization across all experimental designs
#'   }
#' @param power_target Numeric. Target statistical power level (typically 0.8 for 80% power).
#'   Must be between 0 and 1.
#' @param power_precision Numeric. Acceptable precision around power target. Designs with
#'   power within \code{power_target ± power_precision} are considered acceptable.
#'   Must be between 0 and 1.
#' @param MOI Numeric. Multiplicity of infection parameter for experimental design
#'   calculations (default: 10). Used to compute number of captured cells.
#' @param num_targets Integer. Number of target genes in the experiment (default: 100).
#'   Used for cost calculations.
#' @param non_targeting_gRNAs Integer. Number of non-targeting gRNAs in the experiment
#'   (default: 10). Used to calculate total library size and captured cell requirements.
#' @param gRNAs_per_target Integer. Number of gRNAs per target gene (default: 4).
#'   Used to calculate total gRNAs and experimental design parameters.
#' @param cost_per_captured_cell Numeric. Cost per captured cell in dollars
#'   (default: 0.086). Used for library preparation cost calculations.
#' @param cost_per_million_reads Numeric. Cost per million sequencing reads in dollars
#'   (default: 0.374). Used for sequencing cost calculations.
#' @param cost_grid_size Integer. Number of grid points for cost optimization grid
#'   (default: 200). Higher values provide finer resolution but longer computation time.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{optimal_cost_power_df}{Data frame with optimal power-cost combinations,
#'     including columns from input plus minimum cost information and cost precision.}
#'   \item{optimal_cost_grid}{Data frame with nested cost grids for each parameter level,
#'     containing detailed design alternatives within cost precision bounds.}
#' }
#'
#' @details
#' This function implements a three-stage cost optimization process:
#'
#' \strong{Stage 1: Power Filtering}
#' \enumerate{
#'   \item Filters input data to designs achieving power within target ± precision
#'   \item Ensures only viable designs (meeting power requirements) are considered
#' }
#'
#' \strong{Stage 2: Cost Optimization}
#' \enumerate{
#'   \item Groups designs by minimizing variable (e.g., TPM_threshold levels)
#'   \item Identifies minimum cost for each parameter level
#'   \item Computes cost precision (1% of minimum cost) for grid generation
#'   \item Records parameter ranges (min/max cells and reads per cell) for each level
#' }
#'
#' \strong{Stage 3: Design Grid Generation}
#' \enumerate{
#'   \item Creates log-spaced grids within parameter ranges for each level
#'   \item Computes detailed cost components (library + sequencing costs)
#'   \item Filters to designs within cost precision bounds (±1% of minimum cost)
#'   \item Applies sampling to reduce redundant designs while preserving diversity
#' }
#'
#' \strong{Cost Model:}
#' \deqn{Total Cost = Library Cost + Sequencing Cost}
#' \deqn{Library Cost = cost\_per\_captured\_cell \times num\_captured\_cells}
#' \deqn{Sequencing Cost = cost\_per\_million\_reads \times reads\_per\_cell \times num\_captured\_cells / 10^6}
#' \deqn{num\_captured\_cells = \frac{(gRNAs\_per\_target \times num\_targets + non\_targeting\_gRNAs) \times cells\_per\_target}{gRNAs\_per\_target \times MOI}}
#'
#' The function is designed to work with output from \code{cost_power_computation()}
#' and provides fine-grained cost optimization for experimental design selection.
#'
#' @examples
#' \dontrun{
#' # Load pilot data and perform cost-power analysis
#' pilot_data <- get_pilot_data_from_package("K562")
#' cost_results <- cost_power_computation(
#'   minimizing_variable = "TPM_threshold",
#'   fixed_variable = list(minimum_fold_change = 0.8),
#'   baseline_expression_stats = pilot_data$baseline_expression_stats,
#'   library_parameters = pilot_data$library_parameters,
#'   power_target = 0.8,
#'   cost_constraint = 15000
#' )
#'
#' # Find optimal cost-efficient designs
#' optimal_designs <- find_optimal_cost_design(
#'   cost_power_df = cost_results,
#'   minimizing_variable = "TPM_threshold",
#'   power_target = 0.8,
#'   power_precision = 0.02,
#'   num_targets = 100,
#'   non_targeting_gRNAs = 10,
#'   gRNAs_per_target = 4,
#'   cost_grid_size = 100
#' )
#'
#' # Examine optimal designs
#' head(optimal_designs$optimal_cost_power_df)
#'
#' # Examine detailed cost grids
#' optimal_designs$optimal_cost_grid$cost_grid[[1]]
#'
#' # Find globally optimal cost design across all parameters
#' cost_optimal <- find_optimal_cost_design(
#'   cost_power_df = cost_results,
#'   minimizing_variable = "cost",
#'   power_target = 0.8,
#'   power_precision = 0.02,
#'   cost_grid_size = 50
#' )
#' }
#'
#' @seealso
#' \code{\link{cost_power_computation}} for the underlying cost-power analysis
#'
#' @export
find_optimal_cost_design <- function(cost_power_df, minimizing_variable,
                                     power_target, power_precision,
                                     MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4,
                                     cost_per_captured_cell = 0.086, cost_per_million_reads = 0.374,
                                     cost_grid_size = 200){

  # Input validation
  input_check_find_optimal_cost_design(
    cost_power_df = cost_power_df,
    minimizing_variable = minimizing_variable,
    power_target = power_target,
    power_precision = power_precision,
    MOI = MOI,
    num_targets = num_targets,
    non_targeting_gRNAs = non_targeting_gRNAs,
    gRNAs_per_target = gRNAs_per_target,
    cost_per_captured_cell = cost_per_captured_cell,
    cost_per_million_reads = cost_per_million_reads,
    cost_grid_size = cost_grid_size
  )

  # filter the dataframe based on power_target
  cost_power_df_filtered <- cost_power_df |>
    dplyr::filter((overall_power > power_target - power_precision) & (overall_power < power_target + power_precision))

  # Return error if cost_power_df_filtered has dimension 0
  if(nrow(cost_power_df_filtered) == 0){
    stop("0 row preserved after applying power filtering! Try adjusting power_precision!")
  }

  # depend on the minimizing_variable
  switch (minimizing_variable,
    "cost" = {
      # compute the optimal cost dataframe
      optimal_design_df <- cost_power_df_filtered |>
        dplyr::summarise(minimum_cost = min(total_cost),
                         min_reads = min(raw_reads_per_cell),
                         max_reads = max(raw_reads_per_cell),
                         min_cells = min(cells_per_target),
                         max_cells = max(cells_per_target)) |>
        dplyr::ungroup() |>
        dplyr::mutate(cost_precision = unname(minimum_cost / 100))
    },
    "TPM_threshold" = ,
    "minimum_fold_change" = {
      # compute the optimal cost dataframe
      optimal_design_df <- cost_power_df_filtered |>
        dplyr::group_by(.data[[minimizing_variable]]) |>
        dplyr::summarise(minimum_cost = min(total_cost),
                         min_reads = min(raw_reads_per_cell),
                         max_reads = max(raw_reads_per_cell),
                         min_cells = min(cells_per_target),
                         max_cells = max(cells_per_target)) |>
        dplyr::ungroup() |>
        dplyr::mutate(cost_precision = unname(minimum_cost / 100))
    }
  )

  # obtain cost grid
  cost_grid_df <- optimal_design_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      cost_grid = list(
        # expand the grid
        expand.grid(
          cells_per_target = 10^seq(log10(min_cells), log10(max_cells), length.out = cost_grid_size),
          raw_reads_per_cell = 10^seq(log10(min_reads), log10(max_reads), length.out = cost_grid_size)
        ) |>
          dplyr::mutate(
            # compute number of captured cells
            num_captured_cells = ((gRNAs_per_target * num_targets + non_targeting_gRNAs) * cells_per_target / gRNAs_per_target) / MOI,
            # compute the library cost, sequencing cost and total cost
            library_cost = cost_per_captured_cell * num_captured_cells,
            sequencing_cost = cost_per_million_reads * raw_reads_per_cell * num_captured_cells / 1e6,
            total_cost = library_cost + sequencing_cost
          ) |>
          dplyr::filter(total_cost > minimum_cost - cost_precision & total_cost < minimum_cost + cost_precision) |>
          dplyr::mutate(cost_of_interest = round(minimum_cost)) |>
          dplyr::group_by(num_captured_cells, cost_of_interest) |>
          dplyr::slice_sample(n = 1) |>
          dplyr::ungroup() |>
          dplyr::group_by(raw_reads_per_cell, cost_of_interest) |>
          dplyr::slice_sample(n = 1) |>
          dplyr::ungroup()
      )
    ) |>
    dplyr::ungroup()

  # merge optimal cost and power dataframe
  if(minimizing_variable == "cost"){
    optimal_cost_power_df <- cost_power_df_filtered |>
      tidyr::crossing(optimal_design_df)
  }else{
    optimal_cost_power_df <- cost_power_df_filtered |>
      dplyr::left_join(optimal_design_df, by = minimizing_variable)
  }

  # unnest the tibble
  optimal_cost_grid <- cost_grid_df |>
    tidyr::unnest(cost_grid)

  # return the results
  return(list(
    optimal_cost_power_df = optimal_cost_power_df,
    optimal_cost_grid = optimal_cost_grid
  ))
}

#' Calculate optimal experimental design parameters under cost constraints
#'
#' @description
#' This function determines the optimal experimental design parameters (cells per target
#' and reads per cell) given a total cost constraint. It handles two optimization scenarios:
#'
#' \itemize{
#'   \item When \code{reads_per_cell} is NULL: Calculates the maximum reads per cell
#'         achievable given a fixed \code{cells_per_target} within the cost constraint
#'   \item When \code{cells_per_target} is NULL: Calculates the maximum cells per target
#'         achievable given a fixed \code{reads_per_cell} within the cost constraint
#' }
#'
#' @param cost_per_captured_cell Numeric. Cost per captured cell in dollars (default: 0.086).
#'   Used for library preparation cost calculations.
#' @param cost_per_million_reads Numeric. Cost per million sequencing reads in dollars
#'   (default: 0.374). Used for sequencing cost calculations.
#' @param cost_constraint Numeric. Total budget constraint in dollars. Must be positive.
#' @param MOI Numeric. Multiplicity of infection parameter (default: 10). Used to compute
#'   the number of captured cells from cells per target.
#' @param num_targets Integer. Number of target genes in the experiment (default: 100).
#' @param non_targeting_gRNAs Integer. Number of non-targeting gRNAs (default: 10).
#' @param gRNAs_per_target Integer. Number of gRNAs per target gene (default: 4).
#' @param reads_per_cell Numeric or NULL. If provided, this parameter is fixed and
#'   \code{cells_per_target} will be optimized. If NULL, this parameter will be optimized.
#' @param cells_per_target Numeric or NULL. If provided, this parameter is fixed and
#'   \code{reads_per_cell} will be optimized. If NULL, this parameter will be optimized.
#' @param mapping_efficiency Numeric. Mapping efficiency of sequencing platform
#'   (default: 0.72). Used to convert between raw reads and mapped reads.
#'
#' @return A list containing:
#' \describe{
#'   \item{cells_per_target}{Numeric. Optimized or provided cells per target value}
#'   \item{reads_per_cell}{Numeric. Optimized or provided reads per cell value}
#' }
#'
#' @details
#' The function uses a cost model where total cost = cell preparation cost + sequencing cost:
#' \itemize{
#'   \item Cell preparation cost = \code{cost_per_captured_cell * num_captured_cells}
#'   \item Sequencing cost = \code{cost_per_million_reads * total_reads / 1e6}
#'   \item \code{num_captured_cells = ((gRNAs_per_target * num_targets + non_targeting_gRNAs) * cells_per_target / gRNAs_per_target) / MOI}
#'   \item \code{total_reads = num_captured_cells * reads_per_cell / mapping_efficiency}
#' }
#'
#' **Optimization Logic:**
#'
#' **Scenario 1 (reads_per_cell = NULL):** Given fixed \code{cells_per_target}, maximizes
#' \code{reads_per_cell} within cost constraint by allocating remaining budget after cell
#' preparation costs to sequencing.
#'
#' **Scenario 2 (cells_per_target = NULL):** Given fixed \code{reads_per_cell}, maximizes
#' \code{cells_per_target} (equivalently, captured cells) within cost constraint using
#' the total cost per captured cell (including both preparation and sequencing costs).
#'
#' @section Error Handling:
#' The function validates that the cost constraint is sufficient for meaningful experiments:
#' \itemize{
#'   \item Ensures at least 10 reads per cell can be achieved
#'   \item Ensures at least 10 captured cells can be achieved
#'   \item Provides clear error messages when constraints are too tight
#' }
#'
#' @examples
#' \dontrun{
#' # Optimize reads per cell given fixed cells per target
#' result1 <- obtain_fixed_variable_constraining_cost(
#'   cost_constraint = 1000,
#'   cells_per_target = 100,
#'   reads_per_cell = NULL
#' )
#'
#' # Optimize cells per target given fixed reads per cell
#' result2 <- obtain_fixed_variable_constraining_cost(
#'   cost_constraint = 1000,
#'   cells_per_target = NULL,
#'   reads_per_cell = 5000
#' )
#' }
#'
#' @seealso
#' \code{\link{cost_power_computation}} for cost-constrained power analysis that uses this function
#'
#' @export
obtain_fixed_variable_constraining_cost <- function(
    cost_per_captured_cell = 0.086, cost_per_million_reads = 0.374, cost_constraint,
    MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4,
    reads_per_cell, cells_per_target, mapping_efficiency = 0.72){

  # Input validation
  if (missing(cost_constraint) || !is.numeric(cost_constraint) || length(cost_constraint) != 1 || cost_constraint <= 0) {
    stop("`cost_constraint` must be a positive numeric value!")
  }

  if (!is.numeric(cost_per_captured_cell) || cost_per_captured_cell <= 0) {
    stop("`cost_per_captured_cell` must be positive!")
  }

  if (!is.numeric(cost_per_million_reads) || cost_per_million_reads <= 0) {
    stop("`cost_per_million_reads` must be positive!")
  }

  if (!is.numeric(mapping_efficiency) || mapping_efficiency <= 0 || mapping_efficiency > 1) {
    stop("`mapping_efficiency` must be between 0 and 1!")
  }

  # Check that exactly one parameter is NULL for optimization
  if (is.null(reads_per_cell) && is.null(cells_per_target)) {
    stop("Exactly one of `reads_per_cell` or `cells_per_target` must be NULL for optimization!")
  }

  if (!is.null(reads_per_cell) && !is.null(cells_per_target)) {
    stop("Exactly one of `reads_per_cell` or `cells_per_target` must be NULL for optimization!")
  }

  # Validate non-NULL parameters
  if (!is.null(reads_per_cell) && (!is.numeric(reads_per_cell) || reads_per_cell <= 0)) {
    stop("`reads_per_cell` must be positive when provided!")
  }

  if (!is.null(cells_per_target) && (!is.numeric(cells_per_target) || cells_per_target <= 0)) {
    stop("`cells_per_target` must be positive when provided!")
  }

  # Optimization logic: determine which parameter to optimize
  if(is.null(reads_per_cell)){
    ################## Scenario 1: Optimize reads_per_cell given fixed cells_per_target ##################

    # Step 1: Calculate number of captured cells needed for the given cells_per_target
    num_captured_cells <- ((gRNAs_per_target * num_targets + non_targeting_gRNAs) * cells_per_target / gRNAs_per_target) / MOI

    # Step 2: Calculate cell preparation costs
    cells_budget <- cost_per_captured_cell * num_captured_cells

    # Step 3: Calculate remaining budget for sequencing
    reads_budget <- cost_constraint - cells_budget
    if(reads_budget < 0){
      stop("Cost constraint is too tight to get any reads per cell!")
    }else{
      # Step 4: Calculate maximum reads per cell achievable with remaining budget
      # Formula: reads_budget = (total_reads / 1e6) * cost_per_million_reads
      # where total_reads = num_captured_cells * reads_per_cell / mapping_efficiency
      reads_per_cell <- (1e6 * mapping_efficiency * reads_budget / cost_per_million_reads) / num_captured_cells

      # Step 5: Validate minimum meaningful sequencing depth
      if(reads_per_cell < 10){
        stop("Cost constraint is too tight to get more than 10 reads per cell!")
      }
    }
  }else{
    ################## Scenario 2: Optimize cells_per_target given fixed reads_per_cell ##################

    # Step 1: Calculate raw reads per cell (before mapping efficiency)
    raw_reads_per_cell <- reads_per_cell / mapping_efficiency

    # Step 2: Calculate total cost per captured cell (preparation + sequencing)
    cost_per_captured_cell_total <- raw_reads_per_cell * cost_per_million_reads / 1e6 + cost_per_captured_cell

    # Step 3: Calculate maximum number of captured cells within budget
    num_captured_cells <- round(cost_constraint / cost_per_captured_cell_total)

    # Step 4: Convert captured cells back to cells per target
    cells_per_target <- num_captured_cells * MOI * gRNAs_per_target / (gRNAs_per_target * num_targets + non_targeting_gRNAs)

    # Step 5: Validate minimum meaningful cell count
    if(num_captured_cells < 10){
      stop("Cost constraint is too tight to get more than 10 captured cells!")
    }
  }

  # Return optimized experimental design parameters
  return(list(
    cells_per_target = cells_per_target,
    reads_per_cell = reads_per_cell
  ))
}
