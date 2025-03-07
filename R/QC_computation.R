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


#' Compute QC probability when fixed effect size is used
#' @inheritParams compute_distribution_teststat_fixed_es
#' @param n_nonzero_trt_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param n_nonzero_cntrl_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @importFrom PoissonBinomial ppbinom
#' @return Probability of a enhancer-gene pair being filtered due to QC
compute_QC_fixed_es <- function(
    fold_change,
    expression_mean, expression_size,
    num_cntrl_cells, num_cells,
    n_nonzero_trt_thresh = 7, n_nonzero_cntrl_thresh = 7){

  # unique the value
  n_nonzero_trt_thresh <- unique(n_nonzero_trt_thresh)
  n_nonzero_cntrl_thresh <- unique(n_nonzero_trt_thresh)
  num_cntrl_cells <- unique(num_cntrl_cells)

  # compute the nonzero probability for control group (of legnth num_gene)
  cntrl_nonzero_prob <- 1 - compute_zero_prob(fold_change_mean = 1,
                                              expression_mean = expression_mean,
                                              expression_size = expression_size)

  # compute the nonzero probability for treatment group
  trt_nonzero_prob <- 1 - compute_zero_prob(fold_change_mean = fold_change,
                                            expression_mean = expression_mean,
                                            expression_size = expression_size)

  # compute the control group being kept
  cntrl_kept_prob <- stats::pbinom(q = n_nonzero_cntrl_thresh - 1,
                                   size = num_cntrl_cells,
                                   prob = unique(cntrl_nonzero_prob),
                                   lower.tail = FALSE)

  # compute the treatment group being kept
  trt_kept_prob <- PoissonBinomial::ppbinom(x = n_nonzero_trt_thresh - 1,
                                            probs = rep(trt_nonzero_prob, num_cells),
                                            lower.tail = FALSE)

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
