// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef HAS_OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::SparseVector;

// Rough estimate of dispersion (theta)
double theta_rough(const SparseVector<double>& y,
                   const VectorXd&             mu) {
  const int n = mu.size();
  const int nnz = y.nonZeros();
  const double* y_val = y.valuePtr();
  const int* y_idx = y.innerIndexPtr();

  double sum_sq = 0.0;

#pragma omp simd reduction(+:sum_sq)
  for (int i = 0; i < nnz; ++i) {
    const int j = y_idx[i];
    const double res = y_val[i] / mu[j] - 1.0;
    sum_sq += res * res;
  }

  // Add (n - nnz) terms for zero entries, where (0/mu - 1)^2 = 1
  sum_sq += static_cast<double>(n - nnz);

  return static_cast<double>(n) / sum_sq;
}

// Refined moment-based Newton iteration
double theta_refined(double                     t0,
                     const SparseVector<double>& y,
                     const VectorXd&             mu,
                     double                      dfr,
                     int                         limit,
                     double                      eps) {
  const int n = mu.size();
  const int nnz = y.nonZeros();
  const double* y_val = y.valuePtr();
  const int* y_idx = y.innerIndexPtr();

  for (int iter = 0; iter < limit; ++iter) {
    t0 = std::fabs(t0);
    double num = 0.0, den = 0.0;

#pragma omp simd reduction(+:num,den)
    for (int i = 0; i < nnz; ++i) {
      const int j = y_idx[i];
      const double diff = y_val[i] - mu[j];
      num += (diff * diff) / (mu[j] + mu[j] * mu[j] / t0);
      den += (diff * diff) / std::pow(mu[j] + t0, 2);
    }

    // Approximate zero entries by assuming their mu is average
    double mu_avg = mu.sum() / static_cast<double>(n);
    const int n0 = n - nnz;
    if (n0 > 0) {
      double diff0_sq = mu_avg * mu_avg;
      num += n0 * (diff0_sq / (mu_avg + mu_avg * mu_avg / t0));
      den += n0 * (diff0_sq / std::pow(mu_avg + t0, 2));
    }

    double delta = (num - dfr) / den;
    if (std::fabs(delta) < eps) break;
    t0 -= delta;
  }

  return t0;
}

// [[Rcpp::export]]
NumericVector theta_batch_cpp(const Eigen::MappedSparseMatrix<double>& Y,
                              const NumericVector& library_size,
                              const NumericVector& rel_expr,
                              bool rough = false,
                              int n_threads = 0) {
  const int G = Y.rows(); // number of genes
  const int C = Y.cols(); // number of cells
  const double eps = std::pow(DBL_EPSILON, 0.25);
  NumericVector theta(G);

#ifdef HAS_OPENMP
  if (n_threads > 0)
    omp_set_num_threads(n_threads);
#pragma omp parallel for schedule(static)
#endif
  for (int g = 0; g < G; ++g) {
    // Sparse vector for gene g
    SparseVector<double, Eigen::ColMajor> y = Y.row(g);

    // Expected counts for gene g across cells
    const VectorXd mu = Eigen::Map<const VectorXd>(library_size.begin(), C) * rel_expr[g];

    // Estimate dispersion
    double t0 = theta_rough(y, mu);
    double t_est = rough ? t0 : theta_refined(t0, y, mu, C - 1, 50, eps);

    // Clip result for numerical stability
    theta[g] = std::max(0.01, std::min(t_est, 1e3));
  }

  return theta;
}