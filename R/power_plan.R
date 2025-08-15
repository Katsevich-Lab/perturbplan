# Power calculation function using optimized approach for Shiny app integration

# Suppress R CMD check warnings for variables used in dplyr contexts
utils::globalVariables(c("library_size", "reads_per_cell", "num_trt_cells", "num_cntrl_cells", "cells_idx", "reads_idx"))

#' Calculate power grid for app heatmap visualization (lightweight)
#'
#' This function provides lightweight power analysis for the Shiny application heatmap.
#' It only computes overall power for each cell/read combination, not detailed curves.
#'
#' @param fc_expression_info List from extract_fc_expression_info() containing fc_expression_df and expression_dispersion_curve
#' @param cells_reads_df Data frame with pre-computed experimental design containing
#'   reads_per_cell, num_trt_cells, num_cntrl_cells, and library_size columns.
#'   Must be provided (from identify_cell_read_range() -> convert_design_to_dataframe()).
#' @param fdr_target FDR target level
#' @param prop_non_null Proportion of non-null pairs
#' @param side Test sidedness ("left", "right", "both")
#'
#' @return List with power grid, cell/read sequences
#' @export
calculate_power_grid <- function(
  fc_expression_info,
  cells_reads_df,
  fdr_target = 0.05,
  prop_non_null = 0.1,
  side = "left"
) {

  # Validate that cells_reads_df is provided
  if (missing(cells_reads_df) || is.null(cells_reads_df)) {
    stop("cells_reads_df is required and must be provided from identify_cell_read_range() -> convert_design_to_dataframe()")
  }

  # Validate provided experimental design
  required_cols <- c("reads_per_cell", "num_trt_cells", "num_cntrl_cells", "library_size")
  missing_cols <- setdiff(required_cols, colnames(cells_reads_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in cells_reads_df: ", paste(missing_cols, collapse = ", "))
  }

  # Call the lightweight power function (overall power only, no curves)
  power_results <- compute_power_grid_overall(
    cells_reads_df = cells_reads_df,
    fc_expression_info = fc_expression_info,
    fdr_target = fdr_target,
    prop_non_null = prop_non_null,
    side = side
  )

  # Transform to expected format for heatmap
  power_grid <- data.frame(
    cells = round(power_results$num_trt_cells),  # Display treatment cells as integers
    reads = power_results$reads_per_cell,
    power = power_results$overall_power,
    num_cntrl_cells = round(power_results$num_cntrl_cells),  # Pre-computed control cells
    library_size = power_results$library_size  # Pre-computed library sizes
  )

  # Return structure compatible with Shiny app
  list(
    power_grid = power_grid,
    cells_seq = sort(unique(power_grid$cells)),  # Treatment cell sequence for axis
    reads_seq = sort(unique(power_results$reads_per_cell))
  )
}

#' Compute power grid for overall power analysis (no curves)
#'
#' This function computes only overall power for each cell/read combination
#' without the expensive curve calculations. Used for heatmap generation.
#' Experimental design parameters (MOI, control group, etc.) should be pre-computed
#' into treatment and control cell counts in the cells_reads_df.
#'
#' @param cells_reads_df Data frame with columns reads_per_cell, library_size, num_trt_cells, and num_cntrl_cells
#' @param fc_expression_info List from extract_fc_expression_info() containing fc_expression_df and expression_dispersion_curve
#' @param fdr_target Target false discovery rate
#' @param prop_non_null Proportion of non-null hypotheses
#' @param side Test sidedness ("left", "right", "both")
#' @return Data frame with power analysis results (overall power only)
#' @export
compute_power_grid_overall <- function(
    cells_reads_df,
    fc_expression_info,
    fdr_target = 0.05,
    prop_non_null = 0.1,
    side = "left"
){

  ########################## use pre-computed library size ##########################
  # cells_reads_df should always contain pre-computed library_size from identify_cell_read_range()
  if (!"library_size" %in% colnames(cells_reads_df)) {
    stop("library_size column missing from cells_reads_df. Expected pre-computed values from identify_cell_read_range().")
  }

  ############### compute the power for the cells-reads grid ###################
  power_df <- cells_reads_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      overall_power = compute_power_plan_overall(
        # experimental information
        num_trt_cells = num_trt_cells, num_cntrl_cells = num_cntrl_cells, library_size = library_size,
        # analysis information
        multiple_testing_alpha = fdr_target, multiple_testing_method = "BH",
        side = side,
        # separated approach information
        fc_expression_df = fc_expression_info$fc_expression_df,
        prop_non_null = prop_non_null
      )
    )

  # return the output dataframe
  return(power_df)
}

#' Compute overall power for power analysis (lightweight computation)
#'
#' This function computes overall power and BH cutoff without expensive curve calculations.
#' Used for heatmap generation and as a base for detailed curve computation.
#' Treatment and control cell counts should be pre-computed.
#'
#' @param num_trt_cells Number of treatment cells (pre-computed)
#' @param num_cntrl_cells Number of control cells (pre-computed)
#' @param library_size Library size (reads per cell)
#' @param multiple_testing_alpha Alpha level for multiple testing
#' @param multiple_testing_method Multiple testing method
#' @param side Test sidedness
#' @param fc_expression_df Data frame with fold change and expression info
#' @param prop_non_null Proportion of non-null hypotheses
#' @param return_full_results If TRUE, return list with all intermediate results; if FALSE, return only overall power
#' @return Overall power value (scalar) or list with full results depending on return_full_results
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