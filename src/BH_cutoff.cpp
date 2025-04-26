// BH_cutoff.cpp          <<—— replace the old file with this
#include <Rcpp.h>
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
double FDP_estimate_cpp(const NumericVector &mean_list,
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
 *  Wrapper: FDP at a given t                                  *
 *------------------------------------------------------------ */
// [[Rcpp::export]]
double fdp_hat_cpp(double t,
                   const NumericVector &mean_list,
                   const NumericVector &sd_list,
                   const std::string   &side,
                   const NumericVector &QC_prob)
{
  return FDP_estimate_cpp(mean_list, sd_list, side, t, QC_prob);
}

/*------------------------------------------------------------ *
 *  Main grid-search (no batch_size)                           *
 *------------------------------------------------------------ */
// [[Rcpp::export]]
double BH_cutoff_cpp(const NumericVector &mean_list,
                     const NumericVector &sd_list,
                     const std::string   &side,
                     double               multiple_testing_alpha,
                     const NumericVector &QC_prob)
{
  const int n = mean_list.size();
  if (n < 1) stop("mean_list must have at least one element.");

  // accept scalar QC_prob by recycling
  NumericVector qc = QC_prob.size() == 1 ?
  NumericVector(n, QC_prob[0]) : QC_prob;
  if (qc.size() != n) stop("QC_prob must have length 1 or length(mean_list).");

  const double t_max = multiple_testing_alpha;
  const double t_min = multiple_testing_alpha / n;
  const double step  = (t_max - t_min) / (n - 1);

  for (int i = 0; i < n; ++i) {
    double t_i = t_max - i * step;
    if (fdp_hat_cpp(t_i, mean_list, sd_list, side, qc) <= multiple_testing_alpha)
      return t_i;            // first admissible cutoff
  }
  return 0.0;                // none passed
}
