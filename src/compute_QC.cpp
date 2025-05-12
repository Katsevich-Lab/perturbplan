// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "zero_prob.h"
using namespace Rcpp;

/* helper: replicate a scalar probability `p` exactly `n` times -------------*/
static NumericVector replicate_prob(double p, int n) {
  NumericVector out(n);
  std::fill(out.begin(), out.end(), p);
  return out;
}

// -------------------------------------------------------------------------
// [[Rcpp::export]]
double compute_QC_fixed_es_cpp(
    NumericVector fold_change,           // length K (per-gRNA FC)
    NumericVector expression_mean,       // vector; treat as scalar(unique)
    NumericVector expression_size,       // idem
    NumericVector num_cntrl_cells,       // idem
    NumericVector num_cells,             // ***scalar*** mean_num_cells here
    int           n_nonzero_trt_thresh = 7,
    int           n_nonzero_cntrl_thresh = 7)
{
  /* ---- pull unique scalars ------------------------------------------- */
  const double mu       = expression_mean[0];
  const double size     = expression_size[0];
  const int    n_cntrl  = static_cast<int>(num_cntrl_cells[0]);
  const int    mean_nc  = static_cast<int>(num_cells[0]);     // scalar

  const int K = fold_change.size();

  /* ---- 1. non-zero probabilities ------------------------------------ */
  const double p_cntrl_nonzero =
  1.0 - zero_prob(1.0, mu, size);                         // scalar

  NumericVector p_trt_nonzero(K);
  for (int i = 0; i < K; ++i)
    p_trt_nonzero[i] = 1.0 - zero_prob(fold_change[i], mu, size);

  /* ---- 2. control arm: Binomial tail ------------------------------- */
  const double p_cntrl_kept = R::pbinom(
    n_nonzero_cntrl_thresh - 1,
    n_cntrl,
    p_cntrl_nonzero,
    /*lower_tail=*/0, /*log_p=*/0);

  /* ---- 3. treatment arm: Poisson-binomial via R -------------------- */
  // mean_nc = expected number of treatment cells PER gRNA
  const int total_nc = mean_nc * K;
  NumericVector probs = replicate_prob(mean(p_trt_nonzero), total_nc);

  Function ppbinom("ppbinom", Environment::namespace_env("PoissonBinomial"));
  double p_trt_kept = as<double>(ppbinom(
    _["x"]         = n_nonzero_trt_thresh - 1,
    _["probs"]     = probs,
    _["lower.tail"]= false));

  /* ---- 4. final QC probability ------------------------------------- */
  const double kept  = p_cntrl_kept * p_trt_kept;
  return 1.0 - kept;                     // probability that QC filters pair
}
