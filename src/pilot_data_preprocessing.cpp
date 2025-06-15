// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

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
              Rcpp::pow(mu + t0, 2));;
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
    }
  }


  return estimate;
}