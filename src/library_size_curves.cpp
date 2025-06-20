// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

/*------------------------------------------------------------ *
 *  S-M Curve Implementation (C++)                             *
 *------------------------------------------------------------ */

//' Compute effective library size from read depth using UMI saturation curve (C++)
//'
//' @description
//' C++ implementation of the saturation-magnitude (S-M) curve that relates 
//' sequencing reads to unique UMI counts, accounting for PCR amplification 
//' variability and UMI saturation.
//'
//' @param reads_per_cell Numeric vector. Total reads per cell.
//' @param UMI_per_cell Numeric. Maximum UMI per cell parameter from S-M curve fit.
//' @param variation Numeric. Variation parameter characterizing PCR bias from S-M curve fit.
//'
//' @return Numeric vector. Effective library size in UMIs for each read depth.
//'
//' @details
//' This C++ implementation provides significant performance improvements over the R version
//' for large-scale power analysis computations. The S-M curve formula:
//' \deqn{effective\_UMI = UMI\_per\_cell \times (1 - exp(-reads\_per\_cell / UMI\_per\_cell) \times (1 + variation \times reads\_per\_cell^2 / (2 \times UMI\_per\_cell^2)))}
//'
//' @seealso \code{\link{fit_read_UMI_curve}} for R version
//' @export
// [[Rcpp::export]]
NumericVector fit_read_UMI_curve_cpp(NumericVector reads_per_cell, 
                                     double UMI_per_cell, 
                                     double variation) {
  
  // Input validation
  if (UMI_per_cell <= 0) {
    stop("UMI_per_cell must be positive");
  }
  if (variation <= 0) {
    stop("variation must be positive");
  }
  
  int n = reads_per_cell.size();
  NumericVector effective_UMI(n);
  
  // Precompute constants
  double inv_UMI = 1.0 / UMI_per_cell;
  double var_factor = variation / (2.0 * UMI_per_cell * UMI_per_cell);
  
  // Vectorized computation
  for (int i = 0; i < n; i++) {
    double reads = reads_per_cell[i];
    
    // Input validation for each read value
    if (reads < 0) {
      stop("reads_per_cell values must be non-negative");
    }
    
    // S-M curve formula with optimized computation
    double exp_term = exp(-reads * inv_UMI);
    double variation_term = 1.0 + var_factor * reads * reads;
    effective_UMI[i] = UMI_per_cell * (1.0 - exp_term * variation_term);
  }
  
  return effective_UMI;
}

/*------------------------------------------------------------ *
 *  Library Size Range Identification (C++)                   *
 *------------------------------------------------------------ */

//' Identify optimal reads per cell range for power analysis grid (C++)
//'
//' @description
//' C++ implementation that determines the minimum and maximum reads per cell values 
//' for power analysis grid generation using binary search on the S-M curve.
//'
//' @param experimental_platform String. Experimental platform identifier.
//' @param UMI_per_cell Numeric. Maximum UMI per cell parameter.
//' @param variation Numeric. Variation parameter for S-M curve.
//'
//' @return List with min_reads_per_cell and max_reads_per_cell elements.
//'
//' @details
//' This C++ implementation uses efficient binary search to find the reads per cell
//' that achieves approximately 80% UMI saturation. Platform-specific minimums:
//' - "10x Chromium v3": 500 reads/cell
//' - "Other": 1000 reads/cell
//' - Default: 500 reads/cell
//'
//' @seealso \code{\link{identify_library_size_range}} for R version
//' @export
// [[Rcpp::export]]
List identify_library_size_range_cpp(std::string experimental_platform,
                                     double UMI_per_cell,
                                     double variation) {
  
  // Input validation
  if (UMI_per_cell <= 0) {
    stop("UMI_per_cell must be positive");
  }
  if (variation <= 0) {
    stop("variation must be positive");
  }
  
  // Step 1: Determine minimum reads per cell based on experimental platform
  int min_reads_per_cell;
  if (experimental_platform == "10x Chromium v3") {
    min_reads_per_cell = 500;
  } else if (experimental_platform == "Other") {
    min_reads_per_cell = 1000;
  } else {
    min_reads_per_cell = 500;  // Fallback default
  }
  
  // Step 2: Determine maximum reads per cell for ~80% UMI saturation
  double target_UMI = 0.8 * UMI_per_cell;
  double upper_bound = 10.0 * UMI_per_cell;  // Generous upper limit
  
  // Step 3: Check corner case - can we even reach 80% saturation?
  NumericVector upper_reads = NumericVector::create(upper_bound);
  NumericVector upper_bound_UMI_vec = fit_read_UMI_curve_cpp(upper_reads, UMI_per_cell, variation);
  double upper_bound_UMI = upper_bound_UMI_vec[0];
  
  int max_reads_per_cell;
  
  if (upper_bound_UMI < target_UMI) {
    // Corner case: Even generous upper bound doesn't reach 80% saturation
    max_reads_per_cell = static_cast<int>(round(upper_bound));
    
    // Calculate actual saturation percentage for informative message
    double actual_saturation = round(100.0 * upper_bound_UMI / UMI_per_cell * 10.0) / 10.0;
    
    Rcout << "Note: 80% UMI saturation not achievable with practical read depths. "
          << "Using maximum practical depth (" << max_reads_per_cell << " reads/cell) "
          << "which achieves " << actual_saturation << "% saturation." << std::endl;
    
  } else {
    // Step 4: Normal case - use binary search to find 80% saturation point
    double lower_bound = static_cast<double>(min_reads_per_cell);
    double tolerance = 0.01 * UMI_per_cell;  // 1% tolerance for convergence
    
    // Binary search loop
    while (upper_bound - lower_bound > 1.0) {
      double mid_point = (lower_bound + upper_bound) / 2.0;
      
      NumericVector mid_reads = NumericVector::create(mid_point);
      NumericVector current_UMI_vec = fit_read_UMI_curve_cpp(mid_reads, UMI_per_cell, variation);
      double current_UMI = current_UMI_vec[0];
      
      // Check if we've found the target within tolerance
      if (std::abs(current_UMI - target_UMI) < tolerance) {
        max_reads_per_cell = static_cast<int>(round(mid_point));
        break;
      }
      
      // Update bounds for next iteration
      if (current_UMI < target_UMI) {
        lower_bound = mid_point;
      } else {
        upper_bound = mid_point;
      }
    }
    
    // If loop ended without breaking, use the upper bound
    if (upper_bound - lower_bound <= 1.0) {
      max_reads_per_cell = static_cast<int>(round(upper_bound));
    }
  }
  
  // Step 5: Ensure minimum < maximum (sanity check)
  if (min_reads_per_cell >= max_reads_per_cell) {
    Rcout << "Warning: Minimum reads per cell (" << min_reads_per_cell 
          << ") >= maximum (" << max_reads_per_cell 
          << "). Adjusting minimum to ensure valid range." << std::endl;
    min_reads_per_cell = std::max(100, max_reads_per_cell - 1000);
  }
  
  // Return the range as a list
  return List::create(
    Named("min_reads_per_cell") = min_reads_per_cell,
    Named("max_reads_per_cell") = max_reads_per_cell
  );
}

/*------------------------------------------------------------ *
 *  Vectorized S-M Curve for Grid Generation (C++)            *
 *------------------------------------------------------------ */

//' Generate reads per cell grid using S-M curve analysis (C++)
//'
//' @description
//' Convenience function that combines range identification with grid generation
//' for power analysis heatmaps.
//'
//' @param experimental_platform String. Experimental platform identifier.
//' @param UMI_per_cell Numeric. Maximum UMI per cell parameter.
//' @param variation Numeric. Variation parameter for S-M curve.
//' @param grid_size Integer. Number of points in the grid (default: 10).
//'
//' @return NumericVector. Sequence of reads per cell values for grid.
//'
//' @export
// [[Rcpp::export]]
NumericVector generate_reads_grid_cpp(std::string experimental_platform,
                                      double UMI_per_cell,
                                      double variation,
                                      int grid_size = 10) {
  
  // Get the range using our C++ function
  List range_result = identify_library_size_range_cpp(experimental_platform, 
                                                      UMI_per_cell, 
                                                      variation);
  
  int min_reads = range_result["min_reads_per_cell"];
  int max_reads = range_result["max_reads_per_cell"];
  
  // Generate grid
  NumericVector grid(grid_size);
  double step = static_cast<double>(max_reads - min_reads) / (grid_size - 1);
  
  for (int i = 0; i < grid_size; i++) {
    grid[i] = round(min_reads + i * step);
  }
  
  return grid;
}