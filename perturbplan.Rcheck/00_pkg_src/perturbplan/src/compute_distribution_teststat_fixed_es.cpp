// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "var_nb.h"          // <-- our shared helper
using namespace Rcpp;

// [[Rcpp::export]]
List compute_distribution_teststat_fixed_es_cpp(
    NumericVector fold_change,
    NumericVector expression_mean,
    NumericVector expression_size,
    NumericVector num_trt_cells,
    NumericVector num_cntrl_cells,
    NumericVector num_cells)
{
  /* ---- unique scalars ---- */
  const double mu      = expression_mean[0];
  const double size    = expression_size[0];
  const double n_trt   = num_trt_cells[0];
  const double n_cntrl = num_cntrl_cells[0];

  const int K = fold_change.size();
  if (num_cells.size() != K)
    stop("`fold_change` and `num_cells` must have the same length.");

  /* ---- treatment mean per guide ---- */
  NumericVector trt_mu_per_guide = mu * fold_change;       // vectorised
  double trt_expression_mean = sum(trt_mu_per_guide * num_cells) / n_trt;

  /* ---- pooled mean & denominator ---- */
  const double n_total   = n_trt + n_cntrl;
  const double trt_prop  = n_trt   / n_total;
  const double cntrl_prop = 1.0 - trt_prop;

  const double cntrl_expression_mean = mu;
  const double pooled_expression_mean = trt_prop * trt_expression_mean + cntrl_prop * cntrl_expression_mean;

  const double pooled_var = var_nb(pooled_expression_mean, size);
  const double denom_sq   = pooled_var * (1.0 / n_cntrl + 1.0 / n_trt);

  /* ---- sd components ---- */
  const double cntrl_var = var_nb(cntrl_expression_mean, size) / n_cntrl;

  double trt_var = 0.0;
  for (int i = 0; i < K; ++i)
    trt_var += var_nb(trt_mu_per_guide[i], size) * num_cells[i];
  trt_var /= (n_trt * n_trt);

  const double sd   = std::sqrt((cntrl_var + trt_var) / denom_sq);
  const double mean = (trt_expression_mean - cntrl_expression_mean) / std::sqrt(denom_sq);

  return List::create(Named("mean") = mean,
                      Named("sd")   = sd);
}
