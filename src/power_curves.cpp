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