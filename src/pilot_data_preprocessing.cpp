// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef HAS_OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::SparseVector;
using Rcpp::NumericVector;

/*--------------------------------------------------------------
 *  Rough moment estimator of the size parameter
 *--------------------------------------------------------------*/
double theta_rough(const Eigen::SparseVector<double>& y,
                   const Eigen::VectorXd&             mu)
{
  const int n  = mu.size();
  const int nz = y.nonZeros();

  double denom = 0.0;

  // ----- 1. non-zero cells (y != 0) --------------------------------
  // Iterate over raw pointers instead of InnerIterator
  const double* y_val = y.valuePtr();
  const auto* y_idx = y.innerIndexPtr();
  for (int k = 0; k < nz; ++k, ++y_val, ++y_idx) {
    const int    j   = *y_idx;          // column index
    const double ratio = (*y_val) * (1.0 / mu[j]);   // one multiplication
    const double diff  = ratio - 1.0;
    denom += diff * diff;
  }

  // ----- 2. zero cells (y == 0) ------------------------------------
  // (0 / mu − 1)^2 == 1  → simply add their count
  denom += static_cast<double>(n - nz);

  // ----- 3. return rough estimate ----------------------------------
  return static_cast<double>(n) / denom;
}

/*--------------------------------------------------------------
 *  Refined Newton iteration (moment-matching)
 *--------------------------------------------------------------*/
double theta_refined(double                          t0,
                     const Eigen::SparseVector<double>& y,
                     const Eigen::VectorXd&            mu,
                     double                           dfr,
                     int                              limit,
                     double                           eps)
{
  const int n  = mu.size();
  const int nz = y.nonZeros();

  // Pre-allocate arrays of reciprocals for the current t0
  // (updated each Newton step)
  Eigen::VectorXd mu_recip = mu.cwiseInverse();            // 1 / mu
  Eigen::VectorXd mu_sq    = mu.array().square();          // mu^2

  for (int iter = 0; iter < limit; ++iter) {
    t0 = std::fabs(t0);

    double num = 0.0;   // numerator   Σ (diff^2 / denom)
    double den = 0.0;   // denominator Σ (diff^2 / (mu + t0)^2)

    // ---------- 1. non-zero cells only ------------------------
    const double* y_val = y.valuePtr();
    const auto* y_idx = y.innerIndexPtr();
    for (int k = 0; k < nz; ++k, ++y_val, ++y_idx) {
      const int    j   = *y_idx;
      const double diff = (*y_val) - mu[j];

      // denom = mu + mu^2 / t0  -> use reciprocal to save division
      const double denom   = mu[j] + mu_sq[j] * (1.0 / t0);
      const double denom2  = std::pow(mu[j] + t0, 2.0);

      num += (diff * diff) / denom;
      den += (diff * diff) / denom2;
    }

    // ---------- 2. zero cells (y == 0)  -----------------------
    // diff = -mu   => diff^2 = mu^2
    const int n_zero = n - nz;
    if (n_zero > 0) {
      // Use the mean of mu and mu^2 to approximate the aggregate
      const double mean_mu   = mu.mean();
      const double mean_mu2  = mu_sq.mean();

      const double diff2  = mean_mu2;              // (0 - mu)^2
      const double denomZ = mean_mu + mean_mu2 / t0;
      const double denom2Z = std::pow(mean_mu + t0, 2.0);

      num += n_zero * (diff2 / denomZ);
      den += n_zero * (diff2 / denom2Z);
    }

    const double delta = (num - dfr) / den;
    if (std::fabs(delta) < eps) break;
    t0 -= delta;
  }
  return t0;
}

/*--------------------------------------------------------------
 *  Vectorised entry point
 *--------------------------------------------------------------*/
// [[Rcpp::export]]
Rcpp::NumericVector theta_batch_cpp(const Eigen::MappedSparseMatrix<double>& Y,
                                    const NumericVector&             library_size,
                                    const NumericVector&             rel_expr,
                                    bool                              rough      = false,
                                    int                               n_threads  = 0)
{
  const int G = Y.rows();                // genes
  const int C = Y.cols();                // cells
  const double eps = std::pow(DBL_EPSILON, 0.25);

  Rcpp::NumericVector theta(G);

#ifdef HAS_OPENMP
  if (n_threads > 0)
    omp_set_num_threads(n_threads);
#pragma omp parallel for schedule(static)
#endif
  for (int g = 0; g < G; ++g) {

    // Sparse row for current gene
    SparseVector<double, Eigen::ColMajor> y = Y.row(g);

    // Expected counts for this gene
    VectorXd mu = Eigen::Map<const VectorXd>(library_size.begin(), C) *
                  rel_expr[g];

    // Rough estimate
    double t0 = theta_rough(y, mu);

    double t_est = t0;
    if (!rough)
      t_est = theta_refined(t0, y, mu, C - 1, 50, eps);

    // Clip to stable range
    t_est = std::max(0.01, std::min(t_est, 1e3));
    theta[g] = t_est;
  }
  return theta;
}