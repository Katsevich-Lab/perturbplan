// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

/*------------------------------------------------------------ *
 *  Overall Power Computation (C++)                            *
 *------------------------------------------------------------ */

// Forward declarations for functions from other files
List compute_monte_carlo_teststat_cpp(DataFrame fc_expression_df,
                                     double library_size,
                                     double num_trt_cells,
                                     double num_cntrl_cells);

double compute_BH_plan(const NumericVector &mean_list,
                      const NumericVector &sd_list,
                      const std::string   &side,
                      double               multiple_testing_alpha,
                      double               prop_non_null);

NumericVector rejection_computation_cpp(const NumericVector &mean_list,
                                       const NumericVector &sd_list,
                                       const std::string   &side,
                                       double               cutoff);

//' Compute overall power for power analysis (C++)
//'
//' @description
//' C++ implementation of compute_power_plan_overall that provides significant 
//' performance improvements for power analysis computations.
//'
//' @param fc_expression_df DataFrame with fold change and expression info
//' @param library_size Numeric. Effective library size
//' @param num_trt_cells Numeric. Number of treatment cells  
//' @param num_cntrl_cells Numeric. Number of control cells
//' @param multiple_testing_alpha Numeric. FDR target level (default 0.05)
//' @param multiple_testing_method String. Method ("BH" only supported)
//' @param side String. Test sidedness ("left", "right", "both")
//' @param prop_non_null Numeric. Proportion of non-null hypotheses (default 0.1)
//' @param return_full_results Logical. Return full results or just overall power
//'
//' @return Numeric overall power (if return_full_results=FALSE) or List with full results
//'
//' @details
//' This C++ implementation orchestrates existing optimized C++ functions:
//' \itemize{
//'   \item compute_monte_carlo_teststat_cpp(): Monte Carlo test statistics
//'   \item compute_BH_plan(): Benjamini-Hochberg significance cutoff
//'   \item rejection_computation_cpp(): Power calculations
//' }
//' 
//' The function provides identical results to the R version while offering
//' significant performance improvements by eliminating R function call overhead.
//'
//' @seealso \code{\link{compute_power_plan_overall}} for R version
//' @export
// [[Rcpp::export]]
SEXP compute_power_plan_overall_cpp(DataFrame fc_expression_df,
                                   double library_size,
                                   double num_trt_cells, 
                                   double num_cntrl_cells,
                                   double multiple_testing_alpha = 0.05,
                                   std::string multiple_testing_method = "BH",
                                   std::string side = "left",
                                   double prop_non_null = 0.1,
                                   bool return_full_results = false) {
  
  // Input validation
  if (library_size <= 0) {
    stop("library_size must be positive");
  }
  if (num_trt_cells <= 0) {
    stop("num_trt_cells must be positive");
  }
  if (num_cntrl_cells <= 0) {
    stop("num_cntrl_cells must be positive");
  }
  if (multiple_testing_alpha <= 0 || multiple_testing_alpha >= 1) {
    stop("multiple_testing_alpha must be between 0 and 1");
  }
  if (prop_non_null < 0 || prop_non_null > 1) {
    stop("prop_non_null must be between 0 and 1");
  }
  
  // Validate multiple testing method
  if (multiple_testing_method != "BH") {
    stop("Only 'BH' multiple testing method is currently supported");
  }
  
  // Validate side
  if (side != "left" && side != "right" && side != "both") {
    stop("side must be 'left', 'right', or 'both'");
  }
  
  // Step 1: Compute Monte Carlo test statistics (call existing C++ function)
  List mc_results = compute_monte_carlo_teststat_cpp(fc_expression_df, 
                                                    library_size,
                                                    num_trt_cells, 
                                                    num_cntrl_cells);
  
  // Extract vectors for Monte Carlo samples
  NumericVector mc_means = mc_results["means"];
  NumericVector mc_sds = mc_results["sds"];
  
  // Step 2: Compute significance cutoff (call existing C++ function)
  double sig_cutoff = compute_BH_plan(mc_means, mc_sds, side, 
                                     multiple_testing_alpha, prop_non_null);
  
  // Step 3: Compute overall power (call existing C++ function)
  NumericVector mc_powers = rejection_computation_cpp(mc_means, mc_sds, 
                                                      side, sig_cutoff);
  double overall_power = mean(mc_powers);
  
  // Step 4: Return results
  if (return_full_results) {
    return List::create(
      Named("overall_power") = overall_power,
      Named("sig_cutoff") = sig_cutoff,
      Named("num_trt_cells") = num_trt_cells,
      Named("num_cntrl_cells") = num_cntrl_cells
    );
  } else {
    return wrap(overall_power);
  }
}