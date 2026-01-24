// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <complex>
#include <algorithm>
using namespace Rcpp;

/*------------------------------------------------------------ *
 *  S-M Curve Implementation (C++)                             *
 *------------------------------------------------------------ */

// Internal helper for simple ZTNB prediction (backward compatibility within this file)
NumericVector fit_read_UMI_curve_simple(NumericVector reads_per_cell,
                                       double UMI_per_cell,
                                       double variation) {
  int n = reads_per_cell.size();
  NumericVector effective_UMI(n);

  // Precompute constants
  double inv_variation = -1.0 / variation;
  double var_over_UMI = variation / UMI_per_cell;

  // Vectorized computation using simple ZTNB formula
  for (int i = 0; i < n; i++) {
    double reads = reads_per_cell[i];
    if (reads < 0) {
      stop("reads_per_cell values must be non-negative");
    }

    double base = 1.0 + var_over_UMI * reads;
    double power_term = pow(base, inv_variation);
    effective_UMI[i] = UMI_per_cell * (1.0 - power_term);
  }

  return effective_UMI;
}

//' Compute effective library size from read depth using preseqR saturation curve (C++)
//'
//' @description
//' C++ implementation of the preseqR-based saturation curve that relates
//' sequencing reads to unique UMI counts. Supports both ZTNB and RFA methods.
//'
//' @param reads_per_cell Numeric vector. Total reads per cell.
//' @param rSAC_fn_wrapper List. Parameters from library_estimation containing:
//'   \itemize{
//'     \item method_used: "ZTNB" or "RFA"
//'     \item reads_norm: Normalization constant
//'     \item n_cells: Number of cells
//'     \item For ZTNB: L, size, mu
//'     \item For RFA: valid_estimator, coefs_real, coefs_imag, poles_real, poles_imag, or constant_value
//'   }
//'
//' @return Numeric vector. Effective library size in UMIs for each read depth.
//'
//' @details
//' This C++ implementation provides significant performance improvements over the R version
//' for large-scale power analysis computations. Supports two methods:
//' \itemize{
//'   \item ZTNB: Uses L * P(X > 0 | size, mu * t)
//'   \item RFA: Uses rational function approximation with complex arithmetic
//' }
//'
//' @seealso \code{\link{fit_read_UMI_curve}} for R wrapper
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector fit_read_UMI_curve_cpp(NumericVector reads_per_cell,
                                     List rSAC_fn_wrapper) {

  // Validate inputs
  if (reads_per_cell.size() == 0) {
    stop("reads_per_cell has length 0");
  }

  // Extract common parameters
  double reads_norm = as<double>(rSAC_fn_wrapper["reads_norm"]);
  double n_cells = as<double>(rSAC_fn_wrapper["n_cells"]);
  std::string method_used = as<std::string>(rSAC_fn_wrapper["method_used"]);

  // Normalize reads_per_cell to the scale used during estimation
  int n = reads_per_cell.size();
  NumericVector t(n);
  for (int i = 0; i < n; i++) {
    t[i] = reads_per_cell[i] / reads_norm;
  }

  NumericVector predictions(n);

  if (method_used == "ZTNB") {
    // Use ZTNB closed-form formula: L * pnbinom(0, size = size, mu = mu * t, lower.tail = FALSE)
    double L = as<double>(rSAC_fn_wrapper["L"]);
    double size = as<double>(rSAC_fn_wrapper["size"]);
    double mu = as<double>(rSAC_fn_wrapper["mu"]);

    for (int i = 0; i < n; i++) {
      // P(X > 0) = 1 - P(X <= 0) = pnbinom(0, size, mu*t, lower.tail=FALSE)
      // Use pnbinom_mu for mu parameterization
      double mu_scaled = mu * t[i];
      double prob = R::pnbinom_mu(0.0, size, mu_scaled, 0, 0); // lower_tail=0 (gives P(X>0)), log=0
      predictions[i] = (L * prob) / n_cells;
    }

  } else if (method_used == "RFA") {
    // Use RFA (ds.rSAC) formula
    bool valid_estimator = as<bool>(rSAC_fn_wrapper["valid_estimator"]);

    if (!valid_estimator) {
      // Invalid estimator - return constant
      double constant_value = as<double>(rSAC_fn_wrapper["constant_value"]);
      for (int i = 0; i < n; i++) {
        predictions[i] = constant_value / n_cells;
      }
    } else {
      // Valid RFA estimator - use formula: Re(coefs %*% (x/(x - poles)))
      NumericVector coefs_real = as<NumericVector>(rSAC_fn_wrapper["coefs_real"]);
      NumericVector coefs_imag = as<NumericVector>(rSAC_fn_wrapper["coefs_imag"]);
      NumericVector poles_real = as<NumericVector>(rSAC_fn_wrapper["poles_real"]);
      NumericVector poles_imag = as<NumericVector>(rSAC_fn_wrapper["poles_imag"]);

      int num_terms = coefs_real.size();

      for (int i = 0; i < n; i++) {
        double x = t[i];
        std::complex<double> sum(0.0, 0.0);

        for (int j = 0; j < num_terms; j++) {
          std::complex<double> coef(coefs_real[j], coefs_imag[j]);
          std::complex<double> pole(poles_real[j], poles_imag[j]);

          // Compute x / (x - pole)
          std::complex<double> ratio = x / (x - pole);

          sum += coef * ratio;
        }

        predictions[i] = sum.real() / n_cells;
      }
    }
  } else {
    stop("Unknown method_used: " + method_used);
  }

  // Handle edge cases: ensure no negative values, NaN, or Inf
  for (int i = 0; i < n; i++) {
    if (std::isnan(predictions[i]) || std::isinf(predictions[i])) {
      predictions[i] = 0.0;
    }
    if (predictions[i] < 0.0) {
      predictions[i] = 0.0;
    }
  }

  return predictions;
}

/*------------------------------------------------------------ *
 *  Library Size Range Identification (C++)                   *
 *------------------------------------------------------------ */

//' Identify optimal reads per cell range for power analysis grid (C++)
//'
//' @description
//' C++ implementation that determines the minimum and maximum reads per cell values
//' for power analysis grid generation using binary search on the S-M curve.
//' Uses saturation-based thresholds (10% and 98%) instead of platform-specific minimums.
//'
//' @param experimental_platform String. Experimental platform identifier (kept for compatibility, not used).
//' @param UMI_per_cell Numeric. Maximum UMI per cell parameter.
//' @param variation Numeric. Variation parameter for S-M curve.
//'
//' @return List with min_reads_per_cell and max_reads_per_cell elements.
//'
//' @details
//' This C++ implementation uses efficient binary search to find the reads per cell
//' range for power analysis. Uses saturation-based thresholds:
//' - Minimum reads: 10% UMI saturation (dynamic based on UMI_per_cell)
//' - Maximum reads: 98% UMI saturation (diminishing returns beyond this point)
//'
//' @seealso \code{\link{identify_library_size_range}} for R version
//' @keywords internal
//' @export
// [[Rcpp::export]]
List identify_library_size_range_cpp(std::string experimental_platform,
                                     double UMI_per_cell,
                                     double variation) {
  
  // Input validation
  if (UMI_per_cell <= 0) {
    stop("UMI_per_cell must be positive");
  }
  if (variation < 0) {
    stop("variation must be non-negative");
  }
  
  // Step 1: Determine minimum reads per cell based on 10% UMI saturation
  double target_min_UMI = 0.1 * UMI_per_cell;  // 10% saturation target
  double lower_bound = 100.0;   // Start search from 100 reads
  double upper_bound_search = 5.0 * UMI_per_cell;  // Reasonable upper bound for search
  
  // Binary search to find reads that achieve 10% saturation
  double tolerance = 1.0;  // 1 read tolerance
  int min_reads_per_cell = static_cast<int>(lower_bound);  // Default fallback
  
  while (upper_bound_search - lower_bound > tolerance) {
    double mid_reads = (lower_bound + upper_bound_search) / 2.0;
    NumericVector mid_reads_vec = NumericVector::create(mid_reads);
    NumericVector current_UMI_vec = fit_read_UMI_curve_simple(mid_reads_vec, UMI_per_cell, variation);
    double current_UMI = current_UMI_vec[0];
    
    if (current_UMI >= target_min_UMI) {
      upper_bound_search = mid_reads;  // Found target, try fewer reads
    } else {
      lower_bound = mid_reads;  // Need more reads
    }
  }
  
  min_reads_per_cell = static_cast<int>(std::ceil(upper_bound_search));

  // Step 2: Determine maximum reads per cell for ~98% UMI saturation
  double target_UMI = 0.98 * UMI_per_cell;
  double upper_bound = 10.0 * UMI_per_cell;  // Generous upper limit

  // Step 3: Check corner case - can we even reach 98% saturation?
  NumericVector upper_reads = NumericVector::create(upper_bound);
  NumericVector upper_bound_UMI_vec = fit_read_UMI_curve_simple(upper_reads, UMI_per_cell, variation);
  double upper_bound_UMI = upper_bound_UMI_vec[0];
  
  int max_reads_per_cell;
  
  if (upper_bound_UMI < target_UMI) {
    // Corner case: Even generous upper bound doesn't reach 98% saturation
    max_reads_per_cell = static_cast<int>(round(upper_bound));

    // Calculate actual saturation percentage for informative message
    double actual_saturation = round(100.0 * upper_bound_UMI / UMI_per_cell * 10.0) / 10.0;

    Rcout << "Note: 98% UMI saturation not achievable with practical read depths. "
          << "Using maximum practical depth (" << max_reads_per_cell << " reads/cell) "
          << "which achieves " << actual_saturation << "% saturation." << std::endl;
    
  } else {
    // Step 4: Normal case - use binary search to find 98% saturation point
    double lower_bound = static_cast<double>(min_reads_per_cell);
    double tolerance = 0.01 * UMI_per_cell;  // 1% tolerance for convergence
    
    // Binary search loop
    while (upper_bound - lower_bound > 1.0) {
      double mid_point = (lower_bound + upper_bound) / 2.0;
      
      NumericVector mid_reads = NumericVector::create(mid_point);
      NumericVector current_UMI_vec = fit_read_UMI_curve_simple(mid_reads, UMI_per_cell, variation);
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
//' @keywords internal
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

/*------------------------------------------------------------ *
 *  Streamlined Library Size Range Functions (Clean API)     *
 *------------------------------------------------------------ */

//' Identify optimal reads per cell range (streamlined version)
//'
//' @description
//' Streamlined C++ implementation that determines the minimum and maximum reads per cell values
//' for power analysis grid generation using binary search on the S-M curve.
//' Uses saturation-based thresholds (10% and 95%) with a clean API.
//'
//' @param UMI_per_cell Numeric. Maximum UMI per cell parameter.
//' @param variation Numeric. Variation parameter for S-M curve.
//'
//' @return List with min_reads_per_cell and max_reads_per_cell elements.
//'
//' @details
//' This streamlined version removes the unused experimental_platform parameter.
//' Uses efficient binary search to find the reads per cell range for power analysis:
//' - Minimum reads: 10% UMI saturation (dynamic based on UMI_per_cell)
//' - Maximum reads: 98% UMI saturation (diminishing returns beyond this point)
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
List identify_reads_range_cpp(double UMI_per_cell, double variation) {
  
  // Input validation
  if (UMI_per_cell <= 0) {
    stop("UMI_per_cell must be positive");
  }
  if (variation < 0) {
    stop("variation must be non-negative");
  }
  
  // Step 1: Determine minimum reads per cell based on 10% UMI saturation
  double target_min_UMI = 0.1 * UMI_per_cell;  // 10% saturation target
  double lower_bound = 100.0;   // Start search from 100 reads
  double upper_bound_search = 5.0 * UMI_per_cell;  // Reasonable upper bound for search
  
  // Binary search to find reads that achieve 10% saturation
  double tolerance = 1.0;  // 1 read tolerance
  int min_reads_per_cell = static_cast<int>(lower_bound);  // Default fallback
  
  while (upper_bound_search - lower_bound > tolerance) {
    double mid_reads = (lower_bound + upper_bound_search) / 2.0;
    NumericVector mid_reads_vec = NumericVector::create(mid_reads);
    NumericVector current_UMI_vec = fit_read_UMI_curve_simple(mid_reads_vec, UMI_per_cell, variation);
    double current_UMI = current_UMI_vec[0];
    
    if (current_UMI >= target_min_UMI) {
      upper_bound_search = mid_reads;  // Found target, try fewer reads
    } else {
      lower_bound = mid_reads;  // Need more reads
    }
  }
  
  min_reads_per_cell = static_cast<int>(std::ceil(upper_bound_search));

  // Step 2: Determine maximum reads per cell for ~98% UMI saturation
  double target_UMI = 0.98 * UMI_per_cell;
  double upper_bound = 10.0 * UMI_per_cell;  // Generous upper limit

  // Step 3: Check corner case - can we even reach 98% saturation?
  NumericVector upper_reads = NumericVector::create(upper_bound);
  NumericVector upper_bound_UMI_vec = fit_read_UMI_curve_simple(upper_reads, UMI_per_cell, variation);
  double upper_bound_UMI = upper_bound_UMI_vec[0];
  
  int max_reads_per_cell;
  
  if (upper_bound_UMI < target_UMI) {
    // Corner case: Even generous upper bound doesn't reach 98% saturation
    max_reads_per_cell = static_cast<int>(round(upper_bound));

    // Calculate actual saturation percentage for informative message
    double actual_saturation = round(100.0 * upper_bound_UMI / UMI_per_cell * 10.0) / 10.0;

    Rcout << "Note: 98% UMI saturation not achievable with practical read depths. "
          << "Using maximum practical depth (" << max_reads_per_cell << " reads/cell) "
          << "which achieves " << actual_saturation << "% saturation." << std::endl;
    
  } else {
    // Step 4: Normal case - use binary search to find 98% saturation point
    double lower_bound = static_cast<double>(min_reads_per_cell);
    double tolerance = 0.01 * UMI_per_cell;  // 1% tolerance for convergence
    
    // Binary search loop
    while (upper_bound - lower_bound > 1.0) {
      double mid_point = (lower_bound + upper_bound) / 2.0;
      
      NumericVector mid_reads = NumericVector::create(mid_point);
      NumericVector current_UMI_vec = fit_read_UMI_curve_simple(mid_reads, UMI_per_cell, variation);
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

//' Generate reads per cell grid (streamlined version)
//'
//' @description
//' Streamlined convenience function that combines range identification with grid generation
//' for power analysis heatmaps.
//'
//' @param UMI_per_cell Numeric. Maximum UMI per cell parameter.
//' @param variation Numeric. Variation parameter for S-M curve.
//' @param grid_size Integer. Number of points in the grid (default: 10).
//'
//' @return NumericVector. Sequence of reads per cell values for grid.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector generate_reads_grid_streamlined_cpp(double UMI_per_cell,
                                                  double variation,
                                                  int grid_size = 10) {
  
  // Get the range using our streamlined C++ function
  List range_result = identify_reads_range_cpp(UMI_per_cell, variation);
  
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

