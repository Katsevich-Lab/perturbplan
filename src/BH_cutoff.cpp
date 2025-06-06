// BH_cutoff.cpp - Merged BH cutoff implementations
#include <Rcpp.h>
#include <functional>
using namespace Rcpp;

/*------------------------------------------------------------ *
 *  Vectorised rejection probability                           *
 *------------------------------------------------------------ */
// [[Rcpp::export]]
NumericVector rejection_computation_cpp(const NumericVector &mean_list,
                                        const NumericVector &sd_list,
                                        const std::string   &side,
                                        double               cutoff)
{
  const int n = mean_list.size();
  if (sd_list.size() != n)
    stop("mean_list and sd_list must have identical length.");

  NumericVector prob(n);

  if (side == "left") {
    const double thr = R::qnorm(cutoff, 0.0, 1.0, /*lower_tail=*/1, /*log=*/0);
    for (int i = 0; i < n; ++i)
      prob[i] = R::pnorm(thr, mean_list[i], sd_list[i], /*lower_tail=*/1, 0);

  } else if (side == "right") {
    const double thr = R::qnorm(1.0 - cutoff, 0.0, 1.0, 1, 0);
    for (int i = 0; i < n; ++i)
      prob[i] = R::pnorm(thr, mean_list[i], sd_list[i], /*lower_tail=*/0, 0);

  } else if (side == "both" || side == "two.sided") {
    const double thr_hi = R::qnorm(1.0 - cutoff / 2.0, 0.0, 1.0, 1, 0);
    const double thr_lo = R::qnorm(cutoff / 2.0,        0.0, 1.0, 1, 0);
    for (int i = 0; i < n; ++i)
      prob[i] = R::pnorm(thr_hi, mean_list[i], sd_list[i], 0, 0) +
        R::pnorm(thr_lo, mean_list[i], sd_list[i], 1, 0);

  } else {
    stop("side must be 'left', 'right', 'both', or 'two.sided'.");
  }
  return prob;
}

/*------------------------------------------------------------ *
 *  FDP estimate (vector QC_prob)                              *
 *------------------------------------------------------------ */
// [[Rcpp::export]]
double compute_FDP_posthoc(const NumericVector &mean_list,
                           const NumericVector &sd_list,
                           const std::string   &side,
                           double               cutoff,
                           const NumericVector &QC_prob)
{
  const int n = mean_list.size();
  if (sd_list.size() != n || QC_prob.size() != n)
    stop("mean_list, sd_list, and QC_prob must have identical length.");

  NumericVector rej = rejection_computation_cpp(mean_list, sd_list, side, cutoff);

  double num_hypo_adj  = 0.0;
  double rejection_sz  = 0.0;

  for (int i = 0; i < n; ++i) {
    const double w = 1.0 - QC_prob[i];
    num_hypo_adj += w;
    rejection_sz += rej[i] * w;
  }
  return num_hypo_adj * cutoff / rejection_sz;
}


/*------------------------------------------------------------ *
 *  Unified bisection search for BH cutoff                     *
 *------------------------------------------------------------ */
double bisection_search(double lower, double upper, double tolerance, double target_alpha,
                         std::function<double(double)> fdp_function) {
  
  // Check early termination conditions
  double f_low = fdp_function(lower) - target_alpha;
  
  // If even the lower bound has FDP > α, no solution exists
  if (f_low > 0) {
    return 0.0;
  }
  
  // If lower bound satisfies FDP ≤ α exactly (within tolerance), return it
  if (f_low >= 0 || std::abs(f_low) < tolerance) {
    return lower;
  }
  
  // Bisection search with iteration limit
  for (int iter = 0; iter < 100; ++iter) {
    double mid = 0.5 * (lower + upper);
    double f_mid = fdp_function(mid) - target_alpha;
    
    // Check for convergence at midpoint
    if (std::abs(f_mid) < tolerance) {
      return mid;
    }
    
    // Update bounds
    if (f_mid > 0) {
      upper = mid;
    } else {
      lower = mid;
    }
    
    // Check bracket convergence
    if ((upper - lower) < tolerance) {
      break;
    }
  }
  
  return 0.5 * (lower + upper);
}

/*------------------------------------------------------------ *
 *  BH cutoff with bisection search (for standard QC analysis) *
 *------------------------------------------------------------ */
// [[Rcpp::export]]
double compute_BH_posthoc(const NumericVector &mean_list,
                          const NumericVector &sd_list,
                          const std::string   &side,
                          double               multiple_testing_alpha,
                          const NumericVector &QC_prob) {
  int n = mean_list.size();
  if (n < 1) stop("mean_list must have at least one element.");

  // recycle scalar QC_prob
  NumericVector qc = (QC_prob.size() == 1
                        ? NumericVector(n, QC_prob[0])
                          : QC_prob);
  if (qc.size() != n)
    stop("QC_prob must have length 1 or length(mean_list).");

  // Setup bounds and tolerance for post-hoc analysis
  double alpha = multiple_testing_alpha;
  double lower = alpha / n;
  double upper = alpha;
  double tolerance = alpha / n;
  
  // Define FDP function for post-hoc analysis
  auto fdp_function = [&](double t) -> double {
    return compute_FDP_posthoc(mean_list, sd_list, side, t, qc);
  };
  
  // Use unified bisection search
  return bisection_search(lower, upper, tolerance, alpha, fdp_function);
}

/*------------------------------------------------------------ *
 *  FDP estimate for planning analysis with prop_non_null      *
 *------------------------------------------------------------ */
// [[Rcpp::export]]
double compute_FDP_plan(const NumericVector &mean_list,
                        const NumericVector &sd_list,
                        const std::string   &side,
                        double               cutoff,
                        double               prop_non_null) {
  
  const int n = mean_list.size();
  if (n < 1) stop("mean_list must have at least one element.");
  if (sd_list.size() != n) stop("mean_list and sd_list must have identical length.");
  
  // Compute marginal power: mean of rejection probabilities across all test statistics
  NumericVector rejection_probs = rejection_computation_cpp(mean_list, sd_list, side, cutoff);
  double power_t = mean(rejection_probs);
  
  // FDP formula: FDP(t) = t / (1 - prop_non_null + prop_non_null * power(t))
  double denominator = 1.0 - prop_non_null + prop_non_null * power_t;
  
  return cutoff / denominator;
}

/*------------------------------------------------------------ *
 *  BH cutoff for underspecified analysis with prop_non_null   *
 *------------------------------------------------------------ */
// [[Rcpp::export]]
double compute_BH_plan(const NumericVector &mean_list,
                       const NumericVector &sd_list,
                       const std::string   &side,
                       double               multiple_testing_alpha,
                       double               prop_non_null,
                       int                  num_pairs) {
  
  const int n = mean_list.size();
  if (n < 1) stop("mean_list must have at least one element.");
  if (sd_list.size() != n) stop("mean_list and sd_list must have identical length.");
  
  // Setup bounds and tolerance for planning analysis
  double lower = 1.0 / num_pairs;
  double upper = multiple_testing_alpha;
  double tolerance = 1.0 / num_pairs;
  
  // Define FDP function for planning analysis
  auto fdp_function = [&](double t) -> double {
    return compute_FDP_plan(mean_list, sd_list, side, t, prop_non_null);
  };
  
  // Use unified bisection search
  return bisection_search(lower, upper, tolerance, multiple_testing_alpha, fdp_function);
}
