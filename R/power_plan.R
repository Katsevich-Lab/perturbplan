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