// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

/*------------------------------------------------------------ *
 *  Cell Range Determination for Power Analysis (C++)         *
 *------------------------------------------------------------ */

// Forward declaration for function from other files
double compute_single_power_cpp(
  double num_cells,
  double reads_per_cell,
  DataFrame fc_expression_df,
  double UMI_per_cell,
  double variation,
  double MOI,
  int num_targets,
  int gRNAs_per_target,
  int non_targeting_gRNAs,
  std::string control_group,
  double multiple_testing_alpha,
  std::string side,
  double prop_non_null);

//' Identify optimal cell count range based on power thresholds (C++)
//'
//' @description
//' Determines minimum and maximum cell counts for power analysis using binary search.
//' Uses a cross-search strategy to guarantee logical ordering (min_cells <= max_cells):
//' finds min_cells with best-case reads and max_cells with worst-case reads.
//'
//' @param min_reads_per_cell Numeric. Minimum reads per cell from library size range
//' @param max_reads_per_cell Numeric. Maximum reads per cell from library size range
//' @param fc_expression_df DataFrame with fold change and expression info
//' @param UMI_per_cell Numeric. Maximum UMI per cell parameter from S-M curve
//' @param variation Numeric. Variation parameter from S-M curve
//' @param MOI Numeric. Multiplicity of infection (default 10)
//' @param num_targets Integer. Number of targets (default 100)
//' @param gRNAs_per_target Integer. gRNAs per target (default 4)
//' @param non_targeting_gRNAs Integer. Non-targeting gRNAs (default 10)
//' @param control_group String. Control group type ("complement" or "nt_cells", default "complement")
//' @param multiple_testing_alpha Numeric. FDR target level (default 0.05)
//' @param side String. Test sidedness ("left", "right", "both", default "left")
//' @param prop_non_null Numeric. Proportion of non-null hypotheses (default 0.1)
//' @param min_power_threshold Numeric. Minimum power threshold for min cells (default 0.01)
//' @param max_power_threshold Numeric. Target power threshold for max cells (default 0.8)
//' @param cell_lower_bound Numeric. Lower bound for total cell search (default 100)
//' @param cell_upper_bound Numeric. Upper bound for total cell search (default 1e7)
//'
//' @return List with min_cells, max_cells, and achieved power values
//'
//' @details
//' This function performs two binary searches using a cross-search strategy:
//' \itemize{
//'   \item Find minimum cells: Where power >= min_power_threshold using max_reads_per_cell (best-case)
//'   \item Find maximum cells: Where power >= max_power_threshold using min_reads_per_cell (worst-case)
//' }
//' 
//' This cross-search strategy ensures min_cells <= max_cells and provides robust
//' experimental design ranges from minimally acceptable to well-powered studies.
//'
//' @export
// [[Rcpp::export]]
List identify_cell_range_cpp(
  double min_reads_per_cell,
  double max_reads_per_cell,
  DataFrame fc_expression_df,
  double UMI_per_cell,
  double variation,
  double MOI = 10.0,
  int num_targets = 100,
  int gRNAs_per_target = 4,
  int non_targeting_gRNAs = 10,
  std::string control_group = "complement",
  double multiple_testing_alpha = 0.05,
  std::string side = "left",
  double prop_non_null = 0.1,
  double min_power_threshold = 0.01,
  double max_power_threshold = 0.8,
  double cell_lower_bound = 100.0,
  double cell_upper_bound = 1e7) {
  
  // Input validation
  if (min_reads_per_cell <= 0) {
    stop("min_reads_per_cell must be positive");
  }
  if (max_reads_per_cell <= 0) {
    stop("max_reads_per_cell must be positive");
  }
  if (min_reads_per_cell > max_reads_per_cell) {
    stop("min_reads_per_cell must be <= max_reads_per_cell");
  }
  if (UMI_per_cell <= 0) {
    stop("UMI_per_cell must be positive");
  }
  if (variation < 0) {
    stop("variation must be non-negative");
  }
  if (min_power_threshold <= 0 || min_power_threshold >= 1) {
    stop("min_power_threshold must be between 0 and 1");
  }
  if (max_power_threshold <= 0 || max_power_threshold >= 1) {
    stop("max_power_threshold must be between 0 and 1");
  }
  if (min_power_threshold >= max_power_threshold) {
    stop("min_power_threshold must be < max_power_threshold");
  }
  if (cell_lower_bound <= 0) {
    stop("cell_lower_bound must be positive");
  }
  if (cell_upper_bound <= cell_lower_bound) {
    stop("cell_upper_bound must be > cell_lower_bound");
  }
  
  double tolerance = 1.0;  // 1 cell tolerance for convergence
  
  // Helper function to convert treatment cells to total cells
  auto trt_to_total_cells = [&](double num_trt_cells) -> double {
    // Total gRNAs in the pool
    double total_gRNAs = num_targets * gRNAs_per_target + non_targeting_gRNAs;
    
    // Calculate total cells based on treatment cells and experimental design
    // num_trt_cells = (gRNAs_per_target * total_cells * MOI) / total_gRNAs
    // Therefore: total_cells = (num_trt_cells * total_gRNAs) / (gRNAs_per_target * MOI)
    return (num_trt_cells * total_gRNAs) / (gRNAs_per_target * MOI);
  };
  
  //============================================================================
  // Step 1: Binary search for minimum treatment cells (power >= min_threshold)
  // Use max_reads_per_cell to find where power first reaches min_power_threshold
  // (Best-case sequencing: easier to achieve min power with high read depth)
  //============================================================================
  
  double lower_min_trt = cell_lower_bound;  // Treatment cells
  double upper_min_trt = cell_upper_bound;  // Treatment cells
  double min_trt_cells = cell_upper_bound;  // Default if not found
  double actual_min_power = 0.0;
  bool min_found = false;
  
  // Check if min_power_threshold is achievable within bounds
  double max_total_cells = trt_to_total_cells(cell_upper_bound);
  double power_at_max = compute_single_power_cpp(
    max_total_cells, max_reads_per_cell, fc_expression_df,
    UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
    non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
  );
  
  if (power_at_max < min_power_threshold) {
    Rcout << "Warning: Minimum power threshold (" << min_power_threshold 
          << ") not achievable within treatment cell bounds. Using upper bound." << std::endl;
    min_trt_cells = cell_upper_bound;
    actual_min_power = power_at_max;
  } else {
    // Binary search for minimum treatment cells
    while (upper_min_trt - lower_min_trt > tolerance) {
      double mid_trt_cells = (lower_min_trt + upper_min_trt) / 2.0;
      double mid_total_cells = trt_to_total_cells(mid_trt_cells);
      
      double power = compute_single_power_cpp(
        mid_total_cells, max_reads_per_cell, fc_expression_df,
        UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
        non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
      );
      
      if (power >= min_power_threshold) {
        upper_min_trt = mid_trt_cells;  // Found sufficient power, try fewer treatment cells
        min_found = true;
      } else {
        lower_min_trt = mid_trt_cells;  // Need more treatment cells
      }
    }
    
    min_trt_cells = upper_min_trt;
    double min_total_cells = trt_to_total_cells(min_trt_cells);
    actual_min_power = compute_single_power_cpp(
      min_total_cells, max_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  //============================================================================
  // Step 2: Binary search for maximum treatment cells (power >= max_threshold)
  // Use min_reads_per_cell to find where power reaches max_power_threshold
  // (Cross-search strategy: guarantees min_cells <= max_cells logical ordering)
  //============================================================================
  
  double lower_max_trt = cell_lower_bound;  // Treatment cells
  double upper_max_trt = cell_upper_bound;  // Treatment cells
  double max_trt_cells = cell_upper_bound;  // Default if not found
  double actual_max_power = 0.0;
  bool max_found = false;
  
  // Check if max_power_threshold is achievable within bounds (reuse max_total_cells)
  power_at_max = compute_single_power_cpp(
    max_total_cells, min_reads_per_cell, fc_expression_df,
    UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
    non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
  );
  
  if (power_at_max < max_power_threshold) {
    Rcout << "Warning: Maximum power threshold (" << max_power_threshold 
          << ") not achievable within treatment cell bounds. Using upper bound." << std::endl;
    max_trt_cells = cell_upper_bound;
    actual_max_power = power_at_max;
  } else {
    // Binary search for maximum treatment cells
    while (upper_max_trt - lower_max_trt > tolerance) {
      double mid_trt_cells = (lower_max_trt + upper_max_trt) / 2.0;
      double mid_total_cells = trt_to_total_cells(mid_trt_cells);
      
      double power = compute_single_power_cpp(
        mid_total_cells, min_reads_per_cell, fc_expression_df,
        UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
        non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
      );
      
      if (power >= max_power_threshold) {
        upper_max_trt = mid_trt_cells;  // Found target power, try fewer treatment cells
        max_found = true;
      } else {
        lower_max_trt = mid_trt_cells;  // Need more treatment cells
      }
    }
    
    max_trt_cells = upper_max_trt;
    double max_total_cells_final = trt_to_total_cells(max_trt_cells);
    actual_max_power = compute_single_power_cpp(
      max_total_cells_final, min_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  //============================================================================
  // Step 3: Apply practical thresholds and ensure logical consistency
  //============================================================================
  
  // Apply practical limits to treatment cells (no need to cap at 100 since that's our minimum)
  if (min_trt_cells < 100.0) {
    Rcout << "Info: min_trt_cells (" << min_trt_cells << ") < 100. Setting min_trt_cells to 100." << std::endl;
    min_trt_cells = 100.0;
    // Recalculate actual power at capped value
    double min_total_cells_capped = trt_to_total_cells(min_trt_cells);
    actual_min_power = compute_single_power_cpp(
      min_total_cells_capped, min_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  if (max_trt_cells > cell_upper_bound) {
    Rcout << "Info: max_trt_cells (" << max_trt_cells << ") > " << cell_upper_bound << ". Setting max_trt_cells to " << cell_upper_bound << "." << std::endl;
    max_trt_cells = cell_upper_bound;
    // Recalculate actual power at capped value
    double max_total_cells_capped = trt_to_total_cells(max_trt_cells);
    actual_max_power = compute_single_power_cpp(
      max_total_cells_capped, max_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  // Ensure min_trt_cells <= max_trt_cells after capping
  if (min_trt_cells > max_trt_cells) {
    Rcout << "Warning: After capping, min_trt_cells (" << min_trt_cells << ") > max_trt_cells (" << max_trt_cells 
          << "). Setting both to 100." << std::endl;
    min_trt_cells = 100.0;
    max_trt_cells = 100.0;
    // Recalculate powers at common value
    double common_total_cells = trt_to_total_cells(100.0);
    actual_min_power = compute_single_power_cpp(
      common_total_cells, min_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
    actual_max_power = compute_single_power_cpp(
      common_total_cells, max_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  // Calculate final total cell ranges for return
  double min_total_cells_final = trt_to_total_cells(min_trt_cells);
  double max_total_cells_final = trt_to_total_cells(max_trt_cells);
  
  // Check for degenerate cases and provide user-friendly error messages
  if (std::abs(min_total_cells_final - max_total_cells_final) < 1.0) {
    stop("Power threshold gap too small: min_power_threshold (" + std::to_string(min_power_threshold) + 
         ") and max_power_threshold (" + std::to_string(max_power_threshold) + 
         ") result in the same optimal cell count (" + std::to_string(min_total_cells_final) + 
         " cells). Please increase the gap between power thresholds for meaningful experimental design ranges.");
  }
  
  // Note: Removed check for equal min/max reads per cell - this is valid when user specifies fixed reads per cell
  
  // Return comprehensive results
  return List::create(
    Named("min_cells") = min_total_cells_final,           // Total cells (for backward compatibility)
    Named("max_cells") = max_total_cells_final,           // Total cells (for backward compatibility)
    Named("min_trt_cells") = min_trt_cells,               // Treatment cells (new)
    Named("max_trt_cells") = max_trt_cells,               // Treatment cells (new)
    Named("min_reads_per_cell") = min_reads_per_cell,
    Named("max_reads_per_cell") = max_reads_per_cell,
    Named("min_power_achieved") = actual_min_power,
    Named("max_power_achieved") = actual_max_power,
    Named("min_power_threshold") = min_power_threshold,
    Named("max_power_threshold") = max_power_threshold,
    Named("min_threshold_found") = min_found,
    Named("max_threshold_found") = max_found,
    Named("cell_search_bounds") = NumericVector::create(cell_lower_bound, cell_upper_bound)
  );
}