// BH_cutoff_search.cpp

#include <Rcpp.h>
using namespace Rcpp;

// forward‐declare the helper you already have in another .cpp
double fdp_hat_cpp(double t,
                   const NumericVector &mean_list,
                   const NumericVector &sd_list,
                   const std::string   &side,
                   const NumericVector &QC_prob);

// [[Rcpp::export]]
double BH_cutoff_bi(const NumericVector &mean_list,
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

  double alpha = multiple_testing_alpha;
  double t_low  = alpha / n;
  double t_high = alpha;
  // if even the max cutoff has FDP > α, give up
  if (fdp_hat_cpp(t_low, mean_list, sd_list, side, qc) > alpha)
    return 0.0;

  // tolerance = α / n
  const double tol = alpha / n;

  for (int iter = 0; iter < 100; ++iter) {
    double t_mid = 0.5 * (t_low + t_high);
    double f_mid = fdp_hat_cpp(t_mid, mean_list, sd_list, side, qc) - alpha;
    if (f_mid > 0) {
      t_high = t_mid;
    } else {
      t_low = t_mid;
    }
    if ((t_high - t_low) < tol)
      break;
  }
  return 0.5 * (t_low + t_high);
}
