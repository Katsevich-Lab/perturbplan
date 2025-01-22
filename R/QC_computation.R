# This is a Rscript including the QC probability computation
#' QC probability computation
#'
#' @inheritParams compute_power
#' @param target_cell_vec Number of cell per element (vector; length J).
#'
#' @return QC probability
#' @importFrom stats pbinom
#' @export

QC_prob <- function(effect_size_mean, baseline_expression, size_parameter,
                    control_cell_vec, target_cell_vec, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh){

  # obtain the number of genes and elements
  num_gene <- length(baseline_expression)
  num_element <- length(control_cell_vec)

  # unchanged fold signal
  no_signal <- matrix(1, nrow = nrow(effect_size_mean), ncol = ncol(effect_size_mean),
                      dimnames = list(
                        gRNA_target = names(control_cell_vec),
                        gene = names(baseline_expression)
                      ))

  # compute the nonzero probability for control group (of legnth num_gene)
  prob_vec_ctl <- 1 - zero_prob(effect_size_mean = no_signal,
                                baseline_expression = baseline_expression,
                                size_parameter = size_parameter)

  # compute the nonzero probability for treatment group
  prob_vec_trt <- 1 -  zero_prob(effect_size_mean = effect_size_mean,
                                 baseline_expression = baseline_expression,
                                 size_parameter = size_parameter)

  # construct size matrix for pbinom
  trial_ctl_mat <- matrix(rep(control_cell_vec, num_gene),
                          nrow = num_element, ncol = num_gene,
                          dimnames = list(
                            gRNA_target = names(control_cell_vec),
                            gene = names(baseline_expression)
                          ))
  trial_trt_mat <- matrix(rep(target_cell_vec, num_gene),
                          nrow = num_element, ncol = num_gene,
                          dimnames = list(
                            gRNA_target = names(control_cell_vec),
                            gene = names(baseline_expression)
                          ))

  # compute the control group being kept
  p_ctl_kept <- pbinom(q = n_nonzero_cntrl_thresh - 1, size = trial_ctl_mat, prob = prob_vec_ctl, lower.tail = FALSE)

  # compute the treatment group being kept
  p_trt_kept <- pbinom(q = n_nonzero_trt_thresh - 1, size = trial_trt_mat, prob = prob_vec_trt, lower.tail = FALSE)

  # compute the kept probability
  kept_prob <- p_ctl_kept * p_trt_kept

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
                                    nrow = num_element, byrow = TRUE,
                                    dimnames = list(
                                      gRNA_target = rownames(effect_size_mean),
                                      gene = colnames(effect_size_mean)
                                    ))
  size_parameter_mat <- matrix(rep(size_parameter, times = num_element),
                               nrow = num_element, byrow = TRUE,
                               dimnames = list(
                                 gRNA_target = rownames(effect_size_mean),
                                 gene = colnames(effect_size_mean)
                               ))

  # compute mean_expression matrix
  mean_expression_mat <- baseline_expression_mat * effect_size_mean

  # return the zero probability matrix
  return((size_parameter_mat / (mean_expression_mat + size_parameter_mat))^size_parameter_mat)
}

#' Compute QC probability for each enhancer-gene pair
#'
#' @inheritParams compute_distribution_teststat
#' @param n_nonzero_trt_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param n_nonzero_cntrl_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#'
#' @importFrom stats pbinom
#' @return Probability of a enhancer-gene pair being filtered due to QC
compute_QC <- function(fold_change_mean, expression_mean, expression_size,
                       num_cntrl_cells, num_trt_cells,
                       n_nonzero_trt_thresh, n_nonzero_cntrl_thresh){

  # compute the nonzero probability for control group (of legnth num_gene)
  cntrl_nonzero_prob <- 1 - compute_zero_prob(fold_change_mean = 1,
                                              expression_mean = expression_mean,
                                              expression_size = expression_size)

  # compute the nonzero probability for treatment group
  trt_nonzero_prob <- 1 - compute_zero_prob(fold_change_mean = fold_change_mean,
                                            expression_mean = expression_mean,
                                            expression_size = expression_size)

  # compute the control group being kept
  cntrl_kept_prob <- pbinom(q = n_nonzero_cntrl_thresh - 1,
                            size = num_cntrl_cells,
                            prob = cntrl_nonzero_prob, lower.tail = FALSE)

  # compute the treatment group being kept
  trt_kept_prob <- pbinom(q = n_nonzero_trt_thresh - 1,
                          size = num_trt_cells,
                          prob = trt_nonzero_prob, lower.tail = FALSE)

  # compute the kept probability
  kept_prob <- cntrl_kept_prob * trt_kept_prob

  # return the QC probability
  return(1 - kept_prob)
}



#' Compute probability mass of NB distribution at zero
#'
#' @inheritParams compute_QC
#'
#' @return Probability of a NB variable being 0
compute_zero_prob <- function(fold_change_mean, expression_mean, expression_size){

  # compute treatment mean expression
  trt_expression_mean <- expression_mean * fold_change_mean

  # compute probability mass of NB distribution being zero
  zero_prob <- (expression_size / (trt_expression_mean + expression_size))^expression_size

  # return the zero probability matrix
  return(zero_prob)
}
