// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef HAS_OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using Eigen::VectorXd;

// Pilot estimate of NB dispersion (theta)
double compute_theta_rough(NumericVector y, NumericVector mu) {
  double n = static_cast<double>(y.size());
  double denom = Rcpp::sum((y / mu - 1) * (y / mu - 1));
  return n / denom;
}

// Refined moment-based estimate of theta
std::vector<double> nb_theta_mm(double t0, NumericVector y, NumericVector mu, double dfr, int limit, double eps) {
  int it = 0;
  double del = 1.0;
  while (++it < limit && std::fabs(del) > eps) {
    t0 = std::fabs(t0);
    del = (
      Rcpp::sum(Rcpp::pow(y - mu, 2) /
                (mu + Rcpp::pow(mu, 2) / t0)) - dfr
    ) / Rcpp::sum(Rcpp::pow(y - mu, 2) /
                  Rcpp::pow(mu + t0, 2));
    t0 -= del;
  }

  double warning = 0.0;
  if (t0 < 0 || it == limit || !std::isfinite(t0)) warning = 1.0;
  return { t0, warning };
}

// [[Rcpp::export]]
double compute_theta_cpp(NumericVector y, NumericVector mu, double dfr, int limit, double eps, bool rough) {
  double t0 = compute_theta_rough(y, mu);
  double estimate = t0;

  if (!rough) {
    try {
      std::vector<double> out = nb_theta_mm(t0, y, mu, dfr, limit, eps);
      estimate = out[0];
    } catch (...) {
      // silent failback to rough
    }
  }

  return estimate;
}

// [[Rcpp::export]]
Rcpp::NumericVector theta_batch_cpp(const Eigen::MappedSparseMatrix<double> &Y,
                                    const Rcpp::NumericVector &library_size,
                                    const Rcpp::NumericVector &rel_expr,
                                    bool rough = false,
                                    int n_threads = 0) {
  const int G = Y.rows(); // number of genes
  const int C = Y.cols(); // number of cells
  const double eps = std::pow(DBL_EPSILON, 0.25);

  Rcpp::NumericVector theta(G);

#ifdef HAS_OPENMP
  if (n_threads > 0)
    omp_set_num_threads(n_threads);

#pragma omp parallel for schedule(static)
#endif
  for (int g = 0; g < G; ++g) {
    // 1. Extract sparse row vector y for gene g
    Eigen::SparseVector<double, Eigen::ColMajor> y = Y.row(g);

    // 2. Convert y to NumericVector
    Rcpp::NumericVector y_vec(y.nonZeros());
    std::copy(y.valuePtr(), y.valuePtr() + y.nonZeros(), y_vec.begin());

    // 3. Compute mu = library_size * rel_expr[g]
    VectorXd mu_eig = Eigen::Map<const VectorXd>(library_size.begin(), C) * rel_expr[g];
    Rcpp::NumericVector mu_vec(mu_eig.data(), mu_eig.data() + C);

    // 4. Estimate theta using compute_theta_cpp
    double th = compute_theta_cpp(y_vec, mu_vec, C - 1, 50, eps, rough);

    // 5. Clip theta to range [0.01, 1000]
    theta[g] = std::min(std::max(th, 0.01), 1e3);
  }

  return theta;
}