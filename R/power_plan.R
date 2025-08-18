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
#' @param cell_recovery_rate Numeric. Rate of cell recovery from loading to capture (default: 0.5).
#' @param QC_rate Numeric. Quality control rate from captured to usable cells (default: 1.0).
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{cells_per_target}{Numeric. Number of treatment cells per target}
#'   \item{reads_per_cell}{Numeric. Sequencing reads per cell}
#'   \item{library_size}{Numeric. Effective library size (UMIs)}
#'   \item{overall_power}{Numeric. Statistical power for this experimental design}
#'   \item{num_usable_cells}{Numeric. Number of usable cells after QC}
#'   \item{num_captured_cells}{Numeric. Number of captured cells}
#'   \item{num_loaded_cells}{Numeric. Number of loaded cells}
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
  mapping_efficiency = 0.72,
  cell_recovery_rate = 0.5,
  QC_rate = 1
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
    dplyr::mutate(num_usable_cells = num_total_cells,
                  num_captured_cells = num_usable_cells / QC_rate,
                  num_loaded_cells = num_captured_cells / cell_recovery_rate,
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
#' @param tpm_threshold Numeric, numeric vector, or character. TPM threshold value, custom sequence, or "varying" for auto-selection.
#' @param minimum_fold_change Numeric, numeric vector, or character. Minimum fold change value, custom sequence, or "varying" for auto-selection.
#' @param cells_per_target Numeric, numeric vector, or character. Number of cells per target, custom sequence, or "varying" for auto-generated grid.
#' @param reads_per_cell Numeric, numeric vector, or character. Reads per cell, custom sequence, or "varying" for auto-generated grid.
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
#' @param mapping_efficiency Numeric. Mapping efficiency for raw reads to usable reads (default: 0.72).
#' @param cell_recovery_rate Numeric. Rate of cell recovery from loading to capture (default: 0.5).
#' @param QC_rate Numeric. Quality control rate from captured to usable cells (default: 1.0).
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
    tpm_threshold, minimum_fold_change, cells_per_target, reads_per_cell,
    # experimental parameters
    MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4, gRNA_variability = 0.13,
    # analysis parameters
    control_group = "complement", side = "left", multiple_testing_alpha = 0.05, prop_non_null = 0.1,
    # data inputs
    baseline_expression_stats, library_parameters,
    # grid parameters
    grid_size = 10, min_power_threshold = 0.01, max_power_threshold = 0.8,
    # efficiency of library preparation and sequencing platform
    mapping_efficiency = 0.72, cell_recovery_rate = 0.5, QC_rate = 1
){

  ####################### construct the tpm_threshold ##########################
  if (is.numeric(tpm_threshold)) {
    tpm_threshold_list <- tpm_threshold  # Works for both single values and vectors
  } else if (tpm_threshold == "varying") {
    tpm_threshold_list <- unname(round(quantile(baseline_expression_stats$relative_expression,
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
          mapping_efficiency = mapping_efficiency,
          cell_recovery_rate = cell_recovery_rate,
          QC_rate = QC_rate
        )
      )
    ) |>
    dplyr::ungroup()

  # Unnest the power grids to create a flat dataframe
  result <- parameter_grid |>
    dplyr::select(-fc_expression_df) |>  # Remove intermediate data
    tidyr::unnest(power_grid) |>
    dplyr::select(minimum_fold_change, tpm_threshold,
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
#'   Options: "tpm_threshold" or "minimum_fold_change". Default: "tpm_threshold".
#' @param fixed_variable List. Fixed values for other analysis parameters. Can include:
#'   \itemize{
#'     \item \code{minimum_fold_change}: Fixed fold change threshold (when optimizing tpm_threshold)
#'     \item \code{tpm_threshold}: Fixed TPM threshold (when optimizing minimum_fold_change)
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
#' @param budget_precision Numeric. Budget utilization factor (default: 0.9).
#'   Filters designs with total cost ≤ budget_precision × cost_constraint.
#' @param cost_per_captured_cell Numeric. Cost per captured cell in dollars (default: 0.086).
#' @param cost_per_million_reads Numeric. Cost per million sequencing reads in dollars (default: 0.374).
#' @param cost_constraint Numeric. Maximum budget constraint in dollars (default: 1e4).
#'   Set to \code{NULL} to disable cost constraints.
#' @param mapping_efficiency Numeric. Sequencing mapping efficiency (default: 0.72).
#' @param cell_recovery_rate Numeric. Cell recovery rate (default: 0.5).
#' @param QC_rate Numeric. Quality control pass rate (default: 1).
#'
#' @return Data frame with power analysis results including:
#'   \itemize{
#'     \item Analysis parameters (tpm_threshold, minimum_fold_change, etc.)
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
#'       \item tpm_threshold: 20 log-spaced values from 1 to 1000 TPM
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
#' result1 <- cost_power_optimization(
#'   minimizing_variable = "tpm_threshold",
#'   fixed_variable = list(minimum_fold_change = 0.8),
#'   baseline_expression_stats = pilot_data$baseline_expression_stats,
#'   library_parameters = pilot_data$library_parameters,
#'   power_target = 0.8,
#'   cost_constraint = 15000,
#'   budget_precision = 0.9
#' )
#'
#' # Optimize fold change with fixed TPM threshold
#' result2 <- cost_power_optimization(
#'   minimizing_variable = "minimum_fold_change",
#'   fixed_variable = list(tpm_threshold = 10),
#'   baseline_expression_stats = pilot_data$baseline_expression_stats,
#'   library_parameters = pilot_data$library_parameters,
#'   power_target = 0.9,
#'   cost_constraint = NULL  # No cost constraint
#' )
#' }
#'
#' @seealso
#' \code{\link{compute_power_plan_full_grid}} for the underlying power analysis
#' \code{\link{cost_computation}} for cost calculation details
#'
#' @export
cost_power_optimization <- function(minimizing_variable = "tpm_threshold", fixed_variable = list(minimum_fold_change = 0.8),
                                    # experimental parameters
                                    MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4, gRNA_variability = 0.13,
                                    # analysis parameters
                                    control_group = "complement", side = "left", multiple_testing_alpha = 0.05, prop_non_null = 0.1,
                                    # data inputs
                                    baseline_expression_stats, library_parameters,
                                    # grid parameters for power
                                    grid_size = 20, power_target = 0.8, power_precision = 0.01,
                                    # grid parameters for budge
                                    budget_precision = 0.9,
                                    # cost parameters
                                    cost_per_captured_cell = 0.086, cost_per_million_reads = 0.374, cost_constraint = NULL,
                                    # efficiency of library preparation and sequencing platform
                                    mapping_efficiency = 0.72, cell_recovery_rate = 0.5, QC_rate = 1){

  # set the seed for computation
  set.seed(1)

  ################# Power-determining parameters setup #########################
  # specify the parameter grid for TPM threshold and minimum_fold_change
  switch (minimizing_variable,
    tpm_threshold = {
      tpm_threshold <- unname(quantile(baseline_expression_stats$relative_expression, probs = seq(0.1, .99, length.out = 20))) * 1e6
      minimum_fold_change <- fixed_variable$minimum_fold_change
    },
    minimum_fold_change = {
      minimum_fold_change <- switch(side,
                                    left = {seq(0.5, 0.9, length.out = 20)},
                                    right = {10^{seq(0, 1, length.out = 20)}},
                                    both = {c(seq(0.5, 0.9, length.out = 10), 10^{seq(0, 1, length.out = 10)})})
      tpm_threshold <- fixed_variable$tpm_threshold
    }
  )

  # specify cells per target and reads per cell
  cells_per_target <- ifelse(is.null(fixed_variable$cells_per_target), "varying", fixed_variable$cells_per_target)
  reads_per_cell <- ifelse(is.null(fixed_variable$reads_per_cell), "varying", fixed_variable$reads_per_cell)

  ########################## Perform grid power search #########################
  # perform power calculation
  power_search <- compute_power_plan_full_grid(
    # power-determining parameters
    tpm_threshold = tpm_threshold, minimum_fold_change = minimum_fold_change, cells_per_target = cells_per_target, reads_per_cell = reads_per_cell,
    # experimental parameters
    MOI = MOI, num_targets = num_targets, non_targeting_gRNAs = non_targeting_gRNAs, gRNAs_per_target = gRNAs_per_target, gRNA_variability = gRNA_variability,
    # analysis parameters
    control_group = control_group, side = side, multiple_testing_alpha = multiple_testing_alpha, prop_non_null = prop_non_null,
    # data inputs
    baseline_expression_stats = baseline_expression_stats, library_parameters = library_parameters,
    # grid parameters
    grid_size = grid_size, min_power_threshold = max(power_target - 0.1, 0.05), max_power_threshold = min(power_target + 0.1, 0.95),
    # efficiency of library preparation and sequencing platform
    mapping_efficiency = mapping_efficiency, cell_recovery_rate = cell_recovery_rate, QC_rate = QC_rate
  )

  # extract the optimal value
  if(!is.null(cost_constraint)){
    # filter the grid based on cost_constraint
    cost_power_rough_search <- power_search |>
      # compute the cost
      dplyr::mutate(
        library_cost = cost_per_captured_cell * num_captured_cells,
        sequencing_cost = cost_per_million_reads * raw_reads_per_cell * num_captured_cells / 1e6,
        total_cost = library_cost + sequencing_cost
      ) |>
      # filter power
      dplyr::filter(
        (overall_power > power_target - power_precision) & (overall_power < power_target + power_precision)
      ) |>
      # constraining the cost
      dplyr::filter(total_cost <= budget_precision * cost_constraint)

    # Return error message when filtering fails!
    if(nrow(cost_power_rough_search) == 0){
      stop("The cost-power optimization failed! Consider adjusting power target or cost budget!")
    }else{
      minimizing_variable_sequence <- cost_power_rough_search |> dplyr::pull(minimizing_variable)
      optimal_value <- switch (minimizing_variable,
        tpm_threshold = {min(minimizing_variable_sequence)},
        minimum_fold_change = {
          minimizing_variable_sequence[which.min(abs(minimizing_variable_sequence - 1))]
        }
      )
    }

    # perform a finer grid search given the optimal value
    switch(minimizing_variable,
           tpm_threshold = {tpm_threshold <- optimal_value},
           minimum_fold_change = {minimum_fold_change <- optimal_value})
    power_results <- compute_power_plan_full_grid(
      # power-determining parameters
      tpm_threshold = tpm_threshold, minimum_fold_change = minimum_fold_change, cells_per_target = cells_per_target, reads_per_cell = reads_per_cell,
      # experimental parameters
      MOI = MOI, num_targets = num_targets, non_targeting_gRNAs = non_targeting_gRNAs, gRNAs_per_target = gRNAs_per_target, gRNA_variability = gRNA_variability,
      # analysis parameters
      control_group = control_group, side = side, multiple_testing_alpha = multiple_testing_alpha, prop_non_null = prop_non_null,
      # data inputs
      baseline_expression_stats = baseline_expression_stats, library_parameters = library_parameters,
      # grid parameters
      grid_size = 100, min_power_threshold = max(power_target - 0.1, 0.05), max_power_threshold = min(power_target + 0.1, 0.95),
      # efficiency of library preparation and sequencing platform
      mapping_efficiency = mapping_efficiency, cell_recovery_rate = cell_recovery_rate, QC_rate = QC_rate
    ) |>
      # compute the cost
      dplyr::mutate(
        library_cost = cost_per_captured_cell * num_captured_cells,
        sequencing_cost = cost_per_million_reads * raw_reads_per_cell * num_captured_cells / 1e6,
        total_cost = library_cost + sequencing_cost
      ) |>
      # select filter the power
      dplyr::filter(
        (overall_power > power_target - power_precision) & (overall_power < power_target + power_precision)
      )

    # return error when cost-power analysis fails!
    if(min(power_results$total_cost) > cost_constraint){
      stop("The cost-power optimization failed! Consider adjusting power target or cost budget!")
    }

  }else{
    # tag each column based on power_target and power_precision
    power_results <- power_search |>
      dplyr::mutate(
        meets_threshold = ifelse(overall_power >= power_target, TRUE, FALSE)
      ) |>
      # select filter the power
      dplyr::filter(
        (overall_power > power_target - power_precision) & (overall_power < power_target + power_precision)
      )
  }

  # return the final results
  return(power_results)
}




