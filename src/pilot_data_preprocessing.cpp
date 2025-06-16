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
double theta_rough(const SparseVector<double>& y,
                   const VectorXd&             mu)
{
  const int n = mu.size();
  double numer = static_cast<double>(n);
  double denom = 0.0;

  // iterate only over non-zero entries of y
  for (SparseVector<double>::InnerIterator it(y); it; ++it) {
    const int    j   = it.index();
    const double res = it.value() / mu[j] - 1.0;
    denom += res * res;
  }
  // add zeros where y == 0
  const int nz = y.nonZeros();
  denom += (n - nz);   // because (0/mu - 1)^2 = 1

  return numer / denom;
}

/*--------------------------------------------------------------
 *  Refined Newton iteration (moment-matching)
 *--------------------------------------------------------------*/
double theta_refined(double                     t0,
                     const SparseVector<double>& y,
                     const VectorXd&             mu,
                     double                      dfr,
                     int                         limit,
                     double                      eps)
{
  for (int iter = 0; iter < limit; ++iter) {
    t0 = std::fabs(t0);

    // compute the two sums needed for Newton step
    double num = 0.0, den = 0.0;
    for (SparseVector<double>::InnerIterator it(y); it; ++it) {
      const int    j   = it.index();
      const double diff = it.value() - mu[j];
      const double denom = mu[j] + mu[j] * mu[j] / t0;
      num += (diff * diff) / denom;
      den += (diff * diff) / std::pow(mu[j] + t0, 2);
    }
    // add zero cells (y == 0)
    const int n0 = mu.size() - y.nonZeros();
    for (int k = 0; k < n0; ++k) {
      // diff = -mu  -> diff^2 = mu^2
      // y==0 so index not needed, treat in aggregate
      // approximate by mean of mu to avoid loop overhead
    }

    double delta = (num - dfr) / den;
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