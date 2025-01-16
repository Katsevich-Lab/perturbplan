# This is a Rscript including the QC probability computation
#' QC probability computation
#'
#' @param baseline_expression Baseline gene expression (vector; length J)
#' @param num_control Number of control cells (vector; length L)
#' @param num_trt Number of treatment cells (vector; length L)
#' @inheritParams power_function
#'
#' @return QC probability
#' @importFrom stats pbinom
#' @export

QC_prob <- function(effect_size_mean, baseline_expression, size_parameter,
                    num_control, num_trt, n_nonzero_trt, n_nonzero_ctl){

  # compute the nonzero probability for control group (of legnth num_gene)
  prob_vec_ctl <- 1 - zero_prob(effect_size_mean = effect_size_mean * 0,
                                baseline_expression = baseline_expression,
                                size_parameter = size_parameter)

  # compute the nonzero probability for treatment group
  prob_vec_trt <- 1 -  zero_prob(effect_size_mean = effect_size_mean,
                                 baseline_expression = baseline_expression,
                                 size_parameter = size_parameter)

  # construct size matrix for pbinom
  num_gene <- length(baseline_expression)
  size_ctl_mat <- matrix(rep(num_control, num_gene), ncol = num_gene)
  size_trt_mat <- matrix(rep(num_trt, num_gene), ncol = num_gene)

  # compute the control group being kept
  p_ctl_kept <- pbinom(q = n_nonzero_ctl - 1, size = size_ctl_mat, prob = prob_vec_ctl, lower.tail = FALSE)

  # compute the treatment group being kept
  p_trt_kept <- pbinom(q = n_nonzero_trt - 1, size = size_trt_mat, prob = prob_vec_trt, lower.tail = FALSE)

  # compute the kept probability
  kept_prob <- as.vector(p_ctl_kept * p_trt_kept)

  # return the QC probability
  return(1 - kept_prob)
}


#' zero probability for each element and gene pair
#'
#' @inheritParams QC_prob
#'
#' @return The zero probability of length L x J
#' @export

zero_prob <- function(effect_size_mean, baseline_expression, size_parameter){

  # construct the matrix with baseline_expression and vector with size parameter
  num_element <- nrow(effect_size_mean)
  baseline_expression_mat <- matrix(rep(baseline_expression, times = num_element),
                                    nrow = num_element, byrow = TRUE)
  size_parameter_mat <- matrix(rep(size_parameter, times = num_element),
                               nrow = num_element, byrow = TRUE)

  # compute mean_expression matrix
  mean_expression_mat <- baseline_expression_mat * effect_size_mean

  # return the zero probability matrix
  return((size_parameter_mat / (mean_expression_mat + size_parameter_mat))^size_parameter_mat)
}
