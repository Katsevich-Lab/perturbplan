// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

/*------------------------------------------------------------ *
 *  Overall Power Computation (C++)                            *
 *------------------------------------------------------------ */

// Forward declarations for functions from other files
double compute_BH_plan(const NumericVector &mean_list,
                      const NumericVector &sd_list,
                      const std::string   &side,
                      double               multiple_testing_alpha,
                      double               prop_non_null);

List compute_distribution_teststat_random_es_cpp(double num_trt_cell, 
                                                  double num_cntrl_cell, 
                                                  double expression_mean, 
                                                  double expression_size, 
                                                  double avg_fold_change, 
                                                  double avg_fold_change_sq);

//' Compute Monte Carlo test statistics for power analysis with random effect sizes
//' 
//' @param fc_expression_df Data frame with fold change and expression information
//' @param library_size Library size parameter
//' @param num_trt_cells Number of treatment cells
//' @param num_cntrl_cells Number of control cells
//' @return List with Monte Carlo mean and standard deviation vectors
//'
//' @export
// [[Rcpp::export]]
List compute_monte_carlo_teststat_cpp(DataFrame fc_expression_df,
                                         double library_size,
                                         double num_trt_cells,
                                         double num_cntrl_cells) {
  
  const int n_samples = fc_expression_df.nrows();
  
  // Extract vectors from data frame
  NumericVector relative_expression = fc_expression_df["relative_expression"];
  NumericVector expression_size = fc_expression_df["expression_size"];
  NumericVector avg_fold_change = fc_expression_df["avg_fold_change"];
  NumericVector avg_fold_change_sq = fc_expression_df["avg_fold_change_sq"];
  
  // Pre-compute expression means
  NumericVector mc_expression_means = library_size * relative_expression;
  
  // Initialize output vectors
  NumericVector mc_means(n_samples);
  NumericVector mc_sds(n_samples);
  
  // Loop through Monte Carlo samples
  for (int i = 0; i < n_samples; i++) {
    // Call the random effect size C++ function
    List test_stats = compute_distribution_teststat_random_es_cpp(
      num_trt_cells,
      num_cntrl_cells,
      mc_expression_means[i],
      expression_size[i],
      avg_fold_change[i],
      avg_fold_change_sq[i]
    );
    
    // Extract results
    mc_means[i] = as<double>(test_stats["mean"]);
    mc_sds[i] = as<double>(test_stats["sd"]);
  }
  
  // Return as list
  return List::create(
    Named("means") = mc_means,
    Named("sds") = mc_sds
  );
}

NumericVector rejection_computation_cpp(const NumericVector &mean_list,
                                       const NumericVector &sd_list,
                                       const std::string   &side,
                                       double               cutoff);

NumericVector fit_read_UMI_curve_cpp(NumericVector reads_per_cell, 
                                     double UMI_per_cell, 
                                     double variation);

//' Compute overall power for power analysis (C++)
//'
//' @description
//' C++ implementation of compute_power_plan_overall that provides significant 
//' performance improvements for power analysis computations. Uses random effect sizes
//' format with avg_fold_change and avg_fold_change_sq columns.
//'
//' @param fc_expression_df DataFrame with fold change and expression info. Must contain:
//'   \itemize{
//'     \item relative_expression: Relative expression levels
//'     \item expression_size: Size parameters for negative binomial distribution
//'     \item avg_fold_change: Average fold change across perturbations
//'     \item avg_fold_change_sq: Average of squared fold changes (second moment)
//'   }
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
//' This C++ implementation uses optimized C++ functions for random effect sizes:
//' \itemize{
//'   \item compute_monte_carlo_teststat_cpp(): Monte Carlo test statistics for random effect sizes
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
  
  // Step 1: Compute Monte Carlo test statistics for random effect sizes
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

/*------------------------------------------------------------ *
 *  Single Power Computation for Cell Range Determination     *
 *------------------------------------------------------------ */

//' Compute power for a single experimental design point (C++)
//'
//' @description
//' Lightweight C++ function that computes power for a single cell count and 
//' read depth combination. Optimized for use in binary search algorithms for 
//' cell range determination.
//'
//' @param num_cells Numeric. Total number of cells in the experiment
//' @param reads_per_cell Numeric. Sequencing reads per cell
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
//'
//' @return Numeric. Overall power for the specified experimental design point
//'
//' @details
//' This function efficiently computes power for a single experimental condition by:
//' \itemize{
//'   \item Converting read depth to library size using S-M curve (fit_read_UMI_curve_cpp)
//'   \item Calculating treatment and control cell counts based on experimental design
//'   \item Computing overall power using compute_power_plan_overall_cpp
//' }
//' 
//' The function is designed for use in binary search algorithms that determine 
//' optimal cell count ranges based on power thresholds.
//'
//' @seealso \code{\link{compute_power_plan_overall_cpp}} for full power analysis
//' @export
// [[Rcpp::export]]
double compute_single_power_cpp(
  double num_cells,
  double reads_per_cell,
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
  double prop_non_null = 0.1) {
  
  // Input validation
  if (num_cells <= 0) {
    stop("num_cells must be positive");
  }
  if (reads_per_cell <= 0) {
    stop("reads_per_cell must be positive");
  }
  if (UMI_per_cell <= 0) {
    stop("UMI_per_cell must be positive");
  }
  if (variation < 0) {
    stop("variation must be non-negative");
  }
  if (MOI <= 0) {
    stop("MOI must be positive");
  }
  if (num_targets <= 0) {
    stop("num_targets must be positive");
  }
  if (gRNAs_per_target <= 0) {
    stop("gRNAs_per_target must be positive");
  }
  if (non_targeting_gRNAs < 0) {
    stop("non_targeting_gRNAs must be non-negative");
  }
  if (control_group != "complement" && control_group != "nt_cells") {
    stop("control_group must be 'complement' or 'nt_cells'");
  }
  
  // Step 1: Convert reads per cell to library size using S-M curve
  NumericVector reads_vec = NumericVector::create(reads_per_cell);
  NumericVector library_size_vec = fit_read_UMI_curve_cpp(reads_vec, UMI_per_cell, variation);
  double library_size = library_size_vec[0];
  
  // Step 2: Calculate treatment and control cell counts based on experimental design
  double total_gRNAs = num_targets * gRNAs_per_target + non_targeting_gRNAs;
  double num_trt_cells = (gRNAs_per_target * num_cells * MOI) / total_gRNAs;
  
  double num_cntrl_cells;
  if (control_group == "complement") {
    num_cntrl_cells = num_cells - num_trt_cells;
  } else { // "nt_cells"
    num_cntrl_cells = (non_targeting_gRNAs * num_cells * MOI) / total_gRNAs;
  }
  
  // Ensure positive cell counts
  if (num_trt_cells <= 0) {
    stop("Calculated num_trt_cells is not positive. Check experimental design parameters.");
  }
  if (num_cntrl_cells <= 0) {
    stop("Calculated num_cntrl_cells is not positive. Check experimental design parameters.");
  }
  
  // Step 3: Compute overall power using existing function
  SEXP power_result = compute_power_plan_overall_cpp(
    fc_expression_df,
    library_size,
    num_trt_cells,
    num_cntrl_cells,
    multiple_testing_alpha,
    "BH",  // Always use BH method
    side,
    prop_non_null,
    false  // return_full_results = false, we only need overall power
  );
  
  // Extract power value and return
  double overall_power = as<double>(power_result);
  return overall_power;
}