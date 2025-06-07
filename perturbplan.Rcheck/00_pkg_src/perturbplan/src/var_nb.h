#ifndef VAR_NB_H            // include-guard
#define VAR_NB_H

/* Negative-binomial variance: Var(Y) = μ + μ² / size                    */
inline double var_nb(double mean, double size) {
  return mean + (mean * mean) / size;
}

#endif
