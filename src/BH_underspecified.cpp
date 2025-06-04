// BH_underspecified.cpp - Underspecified BH cutoff search with prop_non_null

#include <Rcpp.h>
using namespace Rcpp;

// Forward declare rejection computation function
NumericVector rejection_computation_cpp(const NumericVector &mean_list,
                                       const NumericVector &sd_list,
                                       const std::string   &side,
                                       double               cutoff);

// Helper function to compute power using vectorized rejection computation
double compute_marginal_power_cpp(double cutoff,
                                  const NumericVector &mean_list,
                                  const NumericVector &sd_list,
                                  const std::string   &side) {
  NumericVector rejection_probs = rejection_computation_cpp(mean_list, sd_list, side, cutoff);
  return mean(rejection_probs);
}

// [[Rcpp::export]]
double compute_underspecified_BH_cutoff_cpp(const NumericVector &mean_list,
                                            const NumericVector &sd_list,
                                            const std::string   &side,
                                            double               multiple_testing_alpha,
                                            double               prop_non_null,
                                            int                  num_pairs) {
  
  const int n = mean_list.size();
  if (n < 1) stop("mean_list must have at least one element.");
  if (sd_list.size() != n) stop("mean_list and sd_list must have identical length.");
  
  // Helper function: FDP(t) - multiple_testing_alpha
  // FDP(t) = t / (1 - prop_non_null + prop_non_null * power(t))
  auto f = [&](double t) -> double {
    double power_t = compute_marginal_power_cpp(t, mean_list, sd_list, side);
    double denominator = 1.0 - prop_non_null + prop_non_null * power_t;
    return t / denominator - multiple_testing_alpha;
  };
  
  // Define bounds
  double lower = 1.0 / num_pairs;
  double upper = multiple_testing_alpha;
  double f_low = f(lower);
  double f_high = f(upper);
  
  // Tolerance
  const double tol = 1.0 / num_pairs;
  
  // Check if solution exists at lower bound
  if (f_low >= 0) {
    return lower;
  }
  
  // Bisection search with iteration limit
  for (int iter = 0; iter < 100; ++iter) {
    double mid = 0.5 * (lower + upper);
    double f_mid = f(mid);
    
    // Check convergence
    if (f_mid < 0 && std::abs(f_mid) < tol) {
      return mid;
    }
    
    // Update bounds
    if (f_mid > 0) {
      upper = mid;
    } else {
      lower = mid;
    }
    
    // Check bracket convergence
    if ((upper - lower) < tol) {
      break;
    }
  }
  
  return 0.5 * (lower + upper);
}