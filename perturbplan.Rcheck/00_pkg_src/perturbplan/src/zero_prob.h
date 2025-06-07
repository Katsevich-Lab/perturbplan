#ifndef ZERO_PROB_H
#define ZERO_PROB_H

/*  P(Y = 0) for NB(size = r, mean = μ)  =  (r / (μ + r))^r              */
inline double zero_prob(double fold_change_mean,
                        double expression_mean,
                        double expression_size)
{
  const double trt_mu = expression_mean * fold_change_mean;
  const double ratio  = expression_size / (trt_mu + expression_size);
  return std::pow(ratio, expression_size);
}
#endif
