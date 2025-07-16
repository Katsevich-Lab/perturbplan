// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// Forward declarations for functions from other files
NumericVector rejection_computation_cpp(const NumericVector &mean_list,
                                       const NumericVector &sd_list,
                                       const std::string   &side,
                                       double               cutoff);

List compute_distribution_teststat_fixed_es_cpp(
    NumericVector fold_change,
    NumericVector expression_mean,
    NumericVector expression_size,
    NumericVector num_trt_cells,
    NumericVector num_cntrl_cells,
    NumericVector num_cells);

List compute_distribution_teststat_random_es_cpp(
    double num_trt_cell,
    double num_cntrl_cell,
    double expression_mean,
    double expression_size,
    double avg_fold_change,
    double avg_fold_change_sq);

/*------------------------------------------------------------ *
 *  Compute power curve across fold change values (C++)        *
 *------------------------------------------------------------ */
// [[Rcpp::export]]
DataFrame compute_fc_curve_cpp(NumericVector fc_output_grid,
                              DataFrame fc_expression_df,
                              double library_size,
                              double num_trt_cells,
                              double num_cntrl_cells,
                              std::string side,
                              double cutoff) {
  
  // Power curve functionality is temporarily disabled
  stop("Power curve functionality is temporarily disabled. Use heatmap and slice visualization instead.");
  
  const int n_grid = fc_output_grid.size();
  const int n_samples = fc_expression_df.nrows();
  
  // Extract vectors from data frame
  NumericVector fc_df_relative_expr = fc_expression_df["relative_expression"];
  NumericVector fc_df_expr_size = fc_expression_df["expression_size"];
  
  // Pre-compute expression means for efficiency
  NumericVector mc_expression_means = library_size * fc_df_relative_expr;
  
  NumericVector power_values(n_grid);
  
  // Loop through fold change grid points
  for (int i = 0; i < n_grid; i++) {
    double fc_val = fc_output_grid[i];
    
    NumericVector fc_means(n_samples);
    NumericVector fc_sds(n_samples);
    
    // Compute test statistics for this FC across all Monte Carlo expression samples
    for (int j = 0; j < n_samples; j++) {
      // Prepare inputs for C++ function (need vectors, not scalars)
      NumericVector fc_vec(1, fc_val);
      NumericVector expr_mean_vec(1, mc_expression_means[j]);
      NumericVector expr_size_vec(1, fc_df_expr_size[j]);
      NumericVector num_trt_vec(1, num_trt_cells);
      NumericVector num_cntrl_vec(1, num_cntrl_cells);
      NumericVector num_cells_vec(1, num_trt_cells);  // For single gRNA case
      
      List test_stats = compute_distribution_teststat_fixed_es_cpp(
        fc_vec, expr_mean_vec, expr_size_vec,
        num_trt_vec, num_cntrl_vec, num_cells_vec
      );
      
      fc_means[j] = as<double>(test_stats["mean"]);
      fc_sds[j] = as<double>(test_stats["sd"]);
    }
    
    // Compute rejection probabilities and average power
    NumericVector fc_powers = rejection_computation_cpp(fc_means, fc_sds, side, cutoff);
    power_values[i] = mean(fc_powers);
  }
  
  // Return as DataFrame
  return DataFrame::create(
    Named("fold_change") = fc_output_grid,
    Named("power") = power_values
  );
}

/*------------------------------------------------------------ *
 *  Compute power curve across expression values (C++)         *
 *------------------------------------------------------------ */
// [[Rcpp::export]]
DataFrame compute_expression_curve_cpp(NumericVector expr_output_grid,
                                      DataFrame fc_expression_df,
                                      double library_size,
                                      Function expression_dispersion_curve,
                                      double num_trt_cells,
                                      double num_cntrl_cells,
                                      std::string side,
                                      double cutoff) {
  
  // Power curve functionality is temporarily disabled
  stop("Power curve functionality is temporarily disabled. Use heatmap and slice visualization instead.");
  
  const int n_grid = expr_output_grid.size();
  const int n_samples = fc_expression_df.nrows();
  
  // Extract fold change vector from data frame
  NumericVector fc_df_fold_change = fc_expression_df["fold_change"];
  
  NumericVector power_values(n_grid);
  
  // Loop through expression grid points
  for (int i = 0; i < n_grid; i++) {
    double expr_val = expr_output_grid[i];
    double expr_mean = library_size * expr_val;
    
    // Use expression_dispersion_curve to get size parameter
    NumericVector expr_size_result = expression_dispersion_curve(expr_val);
    double expr_size = expr_size_result[0];
    
    NumericVector expr_means(n_samples);
    NumericVector expr_sds(n_samples);
    
    // Compute test statistics for this expression across all Monte Carlo FC samples
    for (int j = 0; j < n_samples; j++) {
      // Prepare inputs for C++ function (need vectors, not scalars)
      NumericVector fc_vec(1, fc_df_fold_change[j]);
      NumericVector expr_mean_vec(1, expr_mean);
      NumericVector expr_size_vec(1, expr_size);
      NumericVector num_trt_vec(1, num_trt_cells);
      NumericVector num_cntrl_vec(1, num_cntrl_cells);
      NumericVector num_cells_vec(1, num_trt_cells);  // For single gRNA case
      
      List test_stats = compute_distribution_teststat_fixed_es_cpp(
        fc_vec, expr_mean_vec, expr_size_vec,
        num_trt_vec, num_cntrl_vec, num_cells_vec
      );
      
      expr_means[j] = as<double>(test_stats["mean"]);
      expr_sds[j] = as<double>(test_stats["sd"]);
    }
    
    // Compute rejection probabilities and average power
    NumericVector expr_powers = rejection_computation_cpp(expr_means, expr_sds, side, cutoff);
    power_values[i] = mean(expr_powers);
  }
  
  // Return as DataFrame
  return DataFrame::create(
    Named("relative_expression") = expr_output_grid,
    Named("power") = power_values
  );
}


/*------------------------------------------------------------ *
 *  Compute Monte Carlo test statistics for random ES (C++)    *
 *------------------------------------------------------------ */
//' Compute Monte Carlo Test Statistics for Random Effect Sizes
//' 
//' @description
//' Computes Monte Carlo test statistics for random effect sizes across multiple
//' expression samples. This function is the counterpart to compute_monte_carlo_teststat_cpp
//' but uses random effect sizes characterized by avg_fold_change and avg_fold_change_sq
//' instead of fixed fold changes.
//' 
//' @param fc_expression_df DataFrame containing Monte Carlo expression samples with columns:
//'   \itemize{
//'     \item \code{relative_expression}: Relative expression levels
//'     \item \code{expression_size}: Size parameters for negative binomial distribution
//'     \item \code{avg_fold_change}: Average fold change across perturbations
//'     \item \code{avg_fold_change_sq}: Average of squared fold changes (second moment)
//'   }
//' @param library_size Numeric. Library size for scaling expression means
//' @param num_trt_cells Numeric. Number of treatment cells
//' @param num_cntrl_cells Numeric. Number of control cells
//' 
//' @return A list containing:
//' \describe{
//'   \item{means}{NumericVector. Asymptotic means of test statistics}
//'   \item{sds}{NumericVector. Asymptotic standard deviations of test statistics}
//' }
//' 
//' @details
//' This function processes Monte Carlo samples where each sample has random effect sizes
//' characterized by their first and second moments (avg_fold_change and avg_fold_change_sq).
//' It calls compute_distribution_teststat_random_es_cpp for each sample to compute the
//' asymptotic distribution of the test statistic.
//' 
//' The key difference from compute_monte_carlo_teststat_cpp is that it handles random
//' effect sizes rather than fixed fold changes, making it suitable for scenarios where
//' perturbation effects vary across cells or conditions.
//' 
//' @seealso \code{\link{compute_monte_carlo_teststat_cpp}} for fixed effect sizes
//' @seealso \code{\link{compute_distribution_teststat_random_es_cpp}} for single sample computation
//' 
//' @export
// [[Rcpp::export]]
List compute_monte_carlo_teststat_new_cpp(DataFrame fc_expression_df,
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