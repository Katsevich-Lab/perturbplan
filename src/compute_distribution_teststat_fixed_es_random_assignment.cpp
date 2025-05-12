// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "var_nb.h"
using namespace Rcpp;

// [[Rcpp::export]]
List compute_distribution_teststat_fixed_es_random_assignment_cpp(
    NumericVector fold_change,
    NumericVector expression_mean,
    NumericVector expression_size,
    NumericVector num_cntrl_cells,
    NumericVector mean_num_cells,
    NumericVector sd_num_cells,
    int B = 1000)
{
  /* ---- constants & checks ---- */
  const int    K              = fold_change.size();          // # gRNAs
  const double mu             = expression_mean[0];
  const double size           = expression_size[0];
  const double n_cntrl        = num_cntrl_cells[0];
  const double mean_n_cells   = mean_num_cells[0];
  const double sd_n_cells     = sd_num_cells[0];

  /* ---- totals & proportions ---- */
  const double total_n_trt = mean_n_cells * K;
  const double n_total     = total_n_trt + n_cntrl;
  const double trt_prop    = total_n_trt / n_total;
  const double cntrl_prop  = 1.0 - trt_prop;

  /* ---- means ---- */
  NumericVector trt_mu_per_guide = mu * fold_change;
  const double  trt_expression_mean   = mean(trt_mu_per_guide);
  const double  cntrl_expression_mean = mu;
  const double  pooled_expression_mean = trt_prop * trt_expression_mean + cntrl_prop * cntrl_expression_mean;

  /* ---- denominator ---- */
  const double pooled_var = var_nb(pooled_expression_mean, size);
  const double denom_sq   = pooled_var * (total_n_trt / n_cntrl + 1.0);

  /* ---- Monte-Carlo for across-gRNA variance ---- */
  NumericVector sim = Rcpp::rnorm(B, total_n_trt, std::sqrt((double)K) * sd_n_cells);

  // clamp negatives to zero (in-place)
  for (double &v : sim)
    if (v < 0.0) v = 0.0;

    const double mean_sim           = mean(sim);
    const double mean_sqrt_sim      = mean(sqrt(sim));
    const double var_sqrt_trt_assign = mean_sim - mean_sqrt_sim * mean_sqrt_sim;

    const double across_across_var = var_sqrt_trt_assign * std::pow(cntrl_expression_mean - trt_expression_mean, 2);

    const double across_within_var = total_n_trt * (mean(pow(trt_mu_per_guide, 2)) - std::pow(mean(trt_mu_per_guide), 2)) / 2.0;

    const double across_var = across_across_var + across_within_var;

    /* ---- within-gRNA variance ---- */
    double within_var = 0.0;
    for (int i = 0; i < K; ++i)
      within_var += var_nb(trt_mu_per_guide[i], size);
    within_var /= K;

    /* ---- sd and mean ---- */
    const double sd   = std::sqrt((across_var + within_var) / denom_sq);
    const double mean = std::sqrt(total_n_trt) * (trt_expression_mean - cntrl_expression_mean) / std::sqrt(denom_sq);

    return List::create(Named("mean") = mean,
                        Named("sd")   = sd);
}
