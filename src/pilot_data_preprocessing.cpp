// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <cmath>
#include <algorithm>     // for std::clamp

#ifdef HAS_OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using Eigen::VectorXd;

// ─────────────────────────────────────────────────────────────────────────────
//  Rough dispersion estimate for one sparse row (no temporary SparseVector)
// ─────────────────────────────────────────────────────────────────────────────
#include <RcppEigen.h>
using Eigen::VectorXd;

inline double theta_rough_row(const Eigen::MappedSparseMatrix<double>& Y,
                              int                               row,
                              const Eigen::VectorXd&            mu) {
  const int n   = mu.size();
  double     ss = 0.0;  // Sum of squared residuals
  int        nz = 0;    // Number of non-zero entries in the row

  // Loop through non-zero elements in the specified row
  for (Eigen::MappedSparseMatrix<double>::InnerIterator it(Y, row); it; ++it) {
    const int j = it.index();          // Column index
    const double mu_j = mu[j];

    // Safety check: mu[j] must be finite and non-zero
    if (!std::isfinite(mu_j) || mu_j == 0.0) {
      return -99;
    }

    const double res = it.value() / mu_j - 1.0;

    // Safety check: residual must be finite
    if (!std::isfinite(res)) {
      return -99;
    }

    ss += res * res;
    ++nz;
  }

  // Add contribution from zero entries (assumed residual = -1)
  ss += static_cast<double>(n - nz);

  // Safety check before final division
  if (ss <= 0.0 || !std::isfinite(ss)) {
    return -99;
  }

  double theta = static_cast<double>(n) / ss;

  // Final check on computed theta
  if (!std::isfinite(theta)) {
    return -99;
  }

  return theta;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Refined Newton iteration (one row)  — with step damping
// ─────────────────────────────────────────────────────────────────────────────
inline double theta_refined_row(double                                  t0,
                                const Eigen::MappedSparseMatrix<double>& Y,
                                int                                      row,
                                const VectorXd&                          mu,
                                int                                      limit = 50,
                                double                                   eps   = 1e-8)
{
  const int n = mu.size();

  for (int iter = 0; iter < limit; ++iter) {
    t0 = std::fabs(t0);                       // keep θ positive

    double num = 0.0, den = 0.0;
    int    nz  = 0;

    // ---------- non-zero counts --------------------------------------------
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(Y, row); it; ++it) {
      const int    j   = it.index();
      const double diff = it.value() - mu[j];
      const double mj   = mu[j];

      num += (diff * diff) / (mj + mj * mj / t0);
      den += (diff * diff) / std::pow(mj + t0, 2);
      ++nz;
    }

    // ---------- zero counts -------------------------------------------------
    const int n0 = n - nz;
    if (n0) {
      const double mu_avg  = mu.mean();
      const double diff0sq = mu_avg * mu_avg;

      num += n0 * (diff0sq / (mu_avg + mu_avg * mu_avg / t0));
      den += n0 * (diff0sq / std::pow(mu_avg + t0, 2));
    }

    // ---------- Newton step with damping -----------------------------------
    if (std::fabs(den) < 1e-12) break;          // avoid division blow-up
    double step = (num - static_cast<double>(n - 1)) / den;

    // step size limitation
    step = std::clamp(step, -5.0, 5.0);

    if (std::fabs(step) < eps * (1.0 + t0))     // convergence
      break;

    t0 -= step;
  }

  // guard against NaN / Inf
  if (!std::isfinite(t0)) t0 = -99;
  return t0;
}


// ─────────────────────────────────────────────────────────────────────────────
//  Newton–Raphson MLE of dispersion (single sparse row, Boost digamma)
// ─────────────────────────────────────────────────────────────────────────────
inline double theta_mle_row(double                                  t0,
                            const Eigen::MappedSparseMatrix<double>&Y,
                            int                                     row,
                            const VectorXd&                         mu,
                            int                                     limit  = 100,
                            double                                  eps    = 1e-8)
{
  const int n = mu.size();
  using boost::math::digamma;
  using boost::math::polygamma;          // polygamma(1,x) = trigamma(x)

  for (int iter = 0; iter < limit; ++iter) {
    t0 = std::fabs(t0);                  // keep positive

    double score = 0.0, info = 0.0;
    int    nz    = 0;

    // ----- non-zero counts ---------------------------------------------------
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(Y, row); it; ++it) {
      const int    j  = it.index();
      const double y  = it.value();
      const double mj = mu[j];

      score += digamma(y + t0) - digamma(t0)
             + std::log(t0) + 1.0
             - std::log(t0 + mj)
             - (y + t0) / (t0 + mj);

      info  += polygamma(1, y + t0) - polygamma(1, t0)
             + 1.0 / t0
             - 2.0 / (t0 + mj)
             + (y + t0) / std::pow(t0 + mj, 2);

      ++nz;
    }

    // ----- zero counts -------------------------------------------------------
    const int n0 = n - nz;
    if (n0) {
      const double mu_avg = mu.mean();

      score += n0 * ( std::log(t0) + 1.0
                    - std::log(t0 + mu_avg)
                    - t0 / (t0 + mu_avg) );

      info  += n0 * ( 1.0 / t0
                    - 2.0 / (t0 + mu_avg)
                    + t0 / std::pow(t0 + mu_avg, 2) );
    }

    // ----- Newton step with damping -----------------------------------------
    if (std::fabs(info) < 1e-9) break;         // singular → give up
    double step = score / info;
    step = std::clamp(step, -5.0, 5.0);        // cap huge jumps

    if (std::fabs(step) < eps * (1.0 + t0))    // converged
      break;

    t0 -= step;
  }

  // final guard
  if (!std::isfinite(t0)) t0 = 1.0;
  return t0;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Vectorised front-end exposed to R
// ─────────────────────────────────────────────────────────────────────────────
// [[Rcpp::export]]
Rcpp::NumericVector theta_batch_cpp(
        const Eigen::MappedSparseMatrix<double>& Y,
        const Rcpp::NumericVector& library_size,
        const Rcpp::NumericVector& rel_expr,
        bool rough = false,
        int n_threads = 0) {

  const int G = Y.rows();  // number of genes
  const int C = Y.cols();  // number of cells
  const double eps = std::pow(DBL_EPSILON, 0.25);

  if (library_size.size() != C || rel_expr.size() != G) {
    Rcpp::stop("Dimension mismatch: rel_expr or library_size does not match Y.");
  }

  Rcpp::Rcout << "[theta_batch_cpp] Starting computation" << std::endl;
  Rcpp::Rcout << " - G (genes) = " << G << std::endl;
  Rcpp::Rcout << " - C (cells) = " << C << std::endl;
  Rcpp::Rcout << " - rel_expr size = " << rel_expr.size() << std::endl;
  Rcpp::Rcout << " - library_size size = " << library_size.size() << std::endl;

  Rcpp::NumericVector theta(G);
  const VectorXd lib_vec = Eigen::Map<const VectorXd>(library_size.begin(), C);

#ifdef HAS_OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
#pragma omp parallel
#endif
  {
    VectorXd mu(C);  // thread-local buffer

#ifdef HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int g = 0; g < G; ++g) {

      if (g < 5) {
        Rcpp::Rcout << "[gene " << g << "] rel_expr = " << rel_expr[g] << std::endl;
      }

      mu = lib_vec * rel_expr[g];

      double t_0 = theta_rough_row(Y, g, mu);
      double t_est = t_0;

      if (!std::isfinite(t_0) || t_0 <= 0.0) {
        Rcpp::Rcout << "[gene " << g << "] Invalid t_0 from rough estimation, using -99" << std::endl;
        theta[g] = -99;
        continue;
      }

      if (g < 5) {
        Rcpp::Rcout << "[gene " << g << "] t_0 = " << t_0 << std::endl;
      }

      try {
        if (rough) {
          t_est = theta_refined_row(t_0, Y, g, mu);
          if (g < 5) Rcpp::Rcout << "[gene " << g << "] theta_refined_row result = " << t_est << std::endl;
        } else {
          t_est = theta_mle_row(t_0, Y, g, mu);
          if (g < 5) Rcpp::Rcout << "[gene " << g << "] theta_mle_row result = " << t_est << std::endl;
        }

        if (!std::isfinite(t_est) || t_est <= 0.0) {
          Rcpp::Rcout << "[gene " << g << "] Invalid t_est after main method, using -99" << std::endl;
          t_est = -99;
        }

      } catch (...) {
        Rcpp::Rcout << "[gene " << g << "] Exception caught during estimation, using -99" << std::endl;
        t_est = -99;
      }

      theta[g] = t_est;
    }
  }

  Rcpp::Rcout << "[theta_batch_cpp] Done." << std::endl;
  return theta;
}