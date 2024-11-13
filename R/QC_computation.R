# This is a Rscript including the QC probability computation
#' QC probability computation
#'
#' @param effect_size Effect size matrix with L row and J column (element, gene)
#' @param baseline_expression Baseline gene expression vector of length J
#' @param size_parameter Size parameter for gene expression of length J
#' @param num_control Number of control cells
#' @param num_trt Number of treatment cells (receiving perturbation)
#' @param n_thresh Cut off QC size
#'
#' @return QC probability
#' @importFrom stats pbinom
#' @export

QC_prob <- function(effect_size, baseline_expression, size_parameter, 
                    num_control, num_trt, n_thresh){
  
  # compute the nonzero probability for control group (of legnth num_gene)
  prob_vec_ctl <- 1 - zero_prob(effect_size = effect_size * 0, 
                                baseline_expression = baseline_expression, 
                                size_parameter = size_parameter)
  
  # compute the nonzero probability for treatment group 
  prob_vec_trt <- 1 -  zero_prob(effect_size = effect_size, 
                                 baseline_expression = baseline_expression, 
                                 size_parameter = size_parameter)
  
  # compute the control group being kept
  p_ctl_kept <- pbinom(q = n_thresh, size = num_control, prob = prob_vec_ctl, lower.tail = FALSE)
  
  # compute the treatment group being kept
  p_trt_kept <- pbinom(q = n_thresh, size = num_trt, prob = prob_vec_trt, lower.tail = FALSE)
  
  # compute the kept probability
  kept_prob <- p_ctl_kept * p_trt_kept
  
  # return the QC probability
  return(1 - kept_prob)
}


#' zero probability for each element and gene pair
#'
#' @param effect_size Effect size matrix with L row and J column (element, gene)
#' @param baseline_expression Baseline gene expression vector of length J
#' @param size_parameter Size parameter for gene expression of length J
#'
#' @return The zero probability of length L * J
#' @export 

zero_prob <- function(effect_size, baseline_expression, size_parameter){
  
  # construct the matrix with baseline_expression and vector with size parameter
  num_element <- nrow(effect_size)
  baseline_expression_mat <- matrix(rep(baseline_expression, times = num_element),
                                    nrow = num_element, byrow = TRUE)
  size_parameter_vec <- as.vector(matrix(rep(size_parameter, times = num_element),
                                         nrow = num_element, byrow = TRUE))
  
  # compute mean_expression and transform the matrix to vector
  mean_expression_vec <- as.vector(baseline_expression_mat * exp(effect_size))
  
  # return the zero probability
  return((size_parameter_vec / (mean_expression_vec + size_parameter_vec))^size_parameter_vec)
}