// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef HAS_OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using Eigen::VectorXd;

// ─────────────────────────────────────────────────────────────────────────────
//  Rough dispersion estimate for one sparse row (no temporary SparseVector)
// ─────────────────────────────────────────────────────────────────────────────
inline double theta_rough_row(const Eigen::MappedSparseMatrix<double>& Y,
                              int                               row,
                              const VectorXd&                   mu) {
  const int n   = mu.size();
  double     ss = 0.0;                // sum of squared residuals
  int        nz = 0;                  // number of non-zero entries

  for (Eigen::MappedSparseMatrix<double>::InnerIterator it(Y, row); it; ++it) {
    const int    j   = it.index();        // column index
    const double res = it.value() / mu[j] - 1.0;
    ss += res * res;
    ++nz;
  }
  ss += static_cast<double>(n - nz);       // contribution of zeros
  return static_cast<double>(n) / ss;      // rough theta
}

// ─────────────────────────────────────────────────────────────────────────────
//  Refined Newton iteration (one row)
// ─────────────────────────────────────────────────────────────────────────────
inline double theta_refined_row(double                           t0,
                                const Eigen::MappedSparseMatrix<double>& Y,
                                int                               row,
                                const VectorXd&                   mu,
                                double                            dfr,
                                int                               limit,
                                double                            eps) {
  const int n   = mu.size();
  const auto* val_ptr = Y.valuePtr();   // direct raw pointers
  const auto* idx_ptr = Y.innerIndexPtr();

  for (int iter = 0; iter < limit; ++iter) {
    t0 = std::fabs(t0);
    double num = 0.0, den = 0.0;
    int    nz  = 0;

    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(Y, row); it; ++it) {
      const int    j    = it.index();
      const double diff = it.value() - mu[j];
      const double mu_j = mu[j];

      num += (diff * diff) / (mu_j + mu_j * mu_j / t0);
      den += (diff * diff) / std::pow(mu_j + t0, 2);
      ++nz;
    }

    // aggregate contribution from zero cells
    const int n0 = n - nz;
    if (n0) {
      const double mu_avg   = mu.mean();
      const double diff0sq  = mu_avg * mu_avg;
      num += n0 * (diff0sq / (mu_avg + mu_avg * mu_avg / t0));
      den += n0 * (diff0sq / std::pow(mu_avg + t0, 2));
    }

    const double delta = (num - dfr) / den;
    if (std::fabs(delta) < eps) break;
    t0 -= delta;
  }
  return t0;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Vectorised front-end exposed to R
// ─────────────────────────────────────────────────────────────────────────────
// [[Rcpp::export]]
Rcpp::NumericVector theta_batch_eigen_optimized(
        const Eigen::MappedSparseMatrix<double>& Y,
        const Rcpp::NumericVector&        library_size,
        const Rcpp::NumericVector&        rel_expr,
        bool                              rough      = false,
        int                               n_threads  = 0) {

  const int G   = Y.rows();                   // genes
  const int C   = Y.cols();                   // cells
  const double eps = std::pow(DBL_EPSILON, 0.25);

  Rcpp::NumericVector theta(G);

  // map once, reuse in every loop
  const VectorXd lib_vec = Eigen::Map<const VectorXd>(library_size.begin(), C);

#ifdef HAS_OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
#pragma omp parallel
#endif
  {
    VectorXd mu(C);                 // thread-local working buffer

#ifdef HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int g = 0; g < G; ++g) {
      // compute mu for this gene (dense, one vector multiply)
      mu = lib_vec * rel_expr[g];

      // rough estimate
      double t_est = theta_rough_row(Y, g, mu);

      // optional refined Newton iteration
      if (!rough)
        t_est = theta_refined_row(t_est, Y, g, mu,
                                  /*dfr=*/C - 1, /*limit=*/50, eps);

      // clip to stable range
      theta[g] = std::max(0.01, std::min(t_est, 1e3));
    }
  }
  return theta;
}