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
//' Uses minimum reads per cell to find where power first reaches 1%, and maximum reads 
//' per cell to find where power reaches 80%.
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
//' @param cell_lower_bound Numeric. Lower bound for cell search (default 100)
//' @param cell_upper_bound Numeric. Upper bound for cell search (default 1e6)
//'
//' @return List with min_cells, max_cells, and achieved power values
//'
//' @details
//' This function performs two binary searches:
//' \itemize{
//'   \item Find minimum cells: Where power ≥ 1% using minimum reads per cell
//'   \item Find maximum cells: Where power ≥ 80% using maximum reads per cell
//' }
//' 
//' The resulting cell range spans from barely useful (1% power) to highly powered 
//' (80% power) experiments, providing guidance for experimental design.
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
  double cell_upper_bound = 1e6) {
  
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
  if (variation <= 0) {
    stop("variation must be positive");
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
  
  //============================================================================
  // Step 1: Binary search for minimum cells (power >= min_threshold)
  // Use min_reads_per_cell to find where power first reaches min_power_threshold
  //============================================================================
  
  double lower_min = cell_lower_bound;
  double upper_min = cell_upper_bound;
  double min_cells = cell_upper_bound;  // Default if not found
  double actual_min_power = 0.0;
  bool min_found = false;
  
  // Check if min_power_threshold is achievable within bounds
  double power_at_max = compute_single_power_cpp(
    cell_upper_bound, min_reads_per_cell, fc_expression_df,
    UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
    non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
  );
  
  if (power_at_max < min_power_threshold) {
    Rcout << "Warning: Minimum power threshold (" << min_power_threshold 
          << ") not achievable within cell bounds. Using upper bound." << std::endl;
    min_cells = cell_upper_bound;
    actual_min_power = power_at_max;
  } else {
    // Binary search for minimum cells
    while (upper_min - lower_min > tolerance) {
      double mid_cells = (lower_min + upper_min) / 2.0;
      
      double power = compute_single_power_cpp(
        mid_cells, min_reads_per_cell, fc_expression_df,
        UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
        non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
      );
      
      if (power >= min_power_threshold) {
        upper_min = mid_cells;  // Found sufficient power, try fewer cells
        min_found = true;
      } else {
        lower_min = mid_cells;  // Need more cells
      }
    }
    
    min_cells = upper_min;
    actual_min_power = compute_single_power_cpp(
      min_cells, min_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  //============================================================================
  // Step 2: Binary search for maximum cells (power >= max_threshold)
  // Use max_reads_per_cell to find where power reaches max_power_threshold
  //============================================================================
  
  double lower_max = cell_lower_bound;
  double upper_max = cell_upper_bound;
  double max_cells = cell_upper_bound;  // Default if not found
  double actual_max_power = 0.0;
  bool max_found = false;
  
  // Check if max_power_threshold is achievable within bounds
  power_at_max = compute_single_power_cpp(
    cell_upper_bound, max_reads_per_cell, fc_expression_df,
    UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
    non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
  );
  
  if (power_at_max < max_power_threshold) {
    Rcout << "Warning: Maximum power threshold (" << max_power_threshold 
          << ") not achievable within cell bounds. Using upper bound." << std::endl;
    max_cells = cell_upper_bound;
    actual_max_power = power_at_max;
  } else {
    // Binary search for maximum cells
    while (upper_max - lower_max > tolerance) {
      double mid_cells = (lower_max + upper_max) / 2.0;
      
      double power = compute_single_power_cpp(
        mid_cells, max_reads_per_cell, fc_expression_df,
        UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
        non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
      );
      
      if (power >= max_power_threshold) {
        upper_max = mid_cells;  // Found target power, try fewer cells
        max_found = true;
      } else {
        lower_max = mid_cells;  // Need more cells
      }
    }
    
    max_cells = upper_max;
    actual_max_power = compute_single_power_cpp(
      max_cells, max_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  //============================================================================
  // Step 3: Apply practical thresholds and ensure logical consistency
  //============================================================================
  
  // Apply practical limits
  if (min_cells > 1000.0) {
    Rcout << "Info: min_cells (" << min_cells << ") > 1000. Setting min_cells to 1000." << std::endl;
    min_cells = 1000.0;
    // Recalculate actual power at capped value
    actual_min_power = compute_single_power_cpp(
      min_cells, min_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  if (max_cells > 1e6) {
    Rcout << "Info: max_cells (" << max_cells << ") > 1e6. Setting max_cells to 1e6." << std::endl;
    max_cells = 1e6;
    // Recalculate actual power at capped value
    actual_max_power = compute_single_power_cpp(
      max_cells, max_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  // Ensure min_cells <= max_cells after capping
  if (min_cells > max_cells) {
    Rcout << "Warning: After capping, min_cells (" << min_cells << ") > max_cells (" << max_cells 
          << "). Setting both to 1000." << std::endl;
    min_cells = 1000.0;
    max_cells = 1000.0;
    // Recalculate powers at common value
    actual_min_power = compute_single_power_cpp(
      min_cells, min_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
    actual_max_power = compute_single_power_cpp(
      max_cells, max_reads_per_cell, fc_expression_df,
      UMI_per_cell, variation, MOI, num_targets, gRNAs_per_target,
      non_targeting_gRNAs, control_group, multiple_testing_alpha, side, prop_non_null
    );
  }
  
  // Return comprehensive results
  return List::create(
    Named("min_cells") = min_cells,
    Named("max_cells") = max_cells,
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