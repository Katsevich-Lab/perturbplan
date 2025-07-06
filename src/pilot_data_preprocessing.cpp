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
  if (!std::isfinite(t0)) t0 = 1.0;
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
      double t_0 = theta_rough_row(Y, g, mu);
      double t_est = t_0;

      try {
        if (rough) {
          // Attempt refined Newton iteration (fast)
          t_est = theta_refined_row(t_0, Y, g, mu);
        } else {
          // Attempt full MLE with digamma
          t_est = theta_mle_row(t_0, Y, g, mu);
        }

        // Check if result is valid
        if (!std::isfinite(t_est) || t_est <= 0.0) {
          // Fall back to rough estimate
          t_est = t_0;
        }

        // If rough fallback also fails, use a conservative default
        if (!std::isfinite(t_est) || t_est <= 0.0) {
          t_est = 1.0;
        }

      } catch (...) {
        // Catch numerical errors from digamma/polygamma/log, etc.
        // Fallback to rough estimate
        t_est = t_0;

        if (!std::isfinite(t_est) || t_est <= 0.0) {
          t_est = 1.0;
        }
      }

      theta[g] = t_est;
    }
  }

  return theta;
}