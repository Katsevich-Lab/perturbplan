# This is a Rscript implementing the help function for power/cost analysis

#' Variance of NB distribution
#'
#' @param mean mean gene expression.
#' @param size size parameter.
#'
#' @return variance of NB distribution.

var_nb <- function(mean, size){

  # compute the variance
  mean + mean^2 / size
}

#' Compute adjusted significance cutoff for multiple testing correction
#'
#' @description
#' This function computes the adjusted significance level (cutoff) for multiple
#' hypothesis testing using either Benjamini-Hochberg (BH) or Bonferroni correction,
#' accounting for quality control failures.
#'
#' @param mean_list Numeric vector. Mean values of test statistics for each hypothesis.
#' @param sd_list Numeric vector. Standard deviation values of test statistics for each hypothesis.
#' @param multiple_testing_alpha Numeric. Target false discovery rate or family-wise error rate.
#' @param multiple_testing_method Character. Multiple testing method, either "BH" or "bonferroni".
#' @param side Character. Test sidedness: "left", "right", or "both".
#' @param QC_prob Numeric vector. The probability of failing QC for each hypothesis.
#'
#' @return Numeric. The adjusted significance level (cutoff threshold).
#'
#' @details
#' The function implements:
#' \itemize{
#'   \item \strong{BH method}: Uses bisection search to find the appropriate cutoff that
#'     controls the false discovery rate at the specified level
#'   \item \strong{Bonferroni method}: Divides the target alpha by the effective number of
#'     tests (accounting for QC failures)
#' }
#'
#' @seealso 
#' \code{\link{BH_cutoff_bisection}} for the BH-specific implementation
#' \code{\link{compute_power_posthoc}} for the main power analysis function
#' @export
adjusted_cutoff <- function(mean_list, sd_list, multiple_testing_alpha, multiple_testing_method, side, QC_prob){

  # compute the adjusted cutoff/significance level
  adjusted_sig_level <- switch (multiple_testing_method,
                                BH = {
                                  BH_cutoff_bisection(mean_list = mean_list,
                                                      sd_list = sd_list,
                                                      multiple_testing_alpha = multiple_testing_alpha,
                                                      side = side,
                                                      QC_prob = QC_prob)
                                },
                                bonferroni = {
                                  num_hypo_adjusted <- sum(1 - QC_prob)
                                  multiple_testing_alpha / num_hypo_adjusted
                                }
  )

  # return the adjusted cutoff/significance level
  return(adjusted_sig_level)
}

#' @useDynLib perturbplan, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Benjamini–Hochberg cutoff with bisection search (C++ back-end)
#'
#' Thin wrapper that validates inputs and forwards to the compiled routine.
#'
#' @inheritParams adjusted_cutoff
#'
#' @return Adjusted cutoff/significance level.
BH_cutoff_bisection <- function(mean_list, sd_list, side, multiple_testing_alpha, QC_prob)
{
  compute_BH_posthoc(mean_list, sd_list, side, multiple_testing_alpha, QC_prob)
}

#' Estimate false discovery proportion (FDP) based on rejection probabilities
#'
#' @description
#' This function estimates the false discovery proportion (FDP) by computing
#' the expected number of false discoveries divided by the expected total
#' number of discoveries, accounting for quality control failures.
#'
#' @param mean_list Numeric vector. Mean values of test statistics for each hypothesis.
#' @param sd_list Numeric vector. Standard deviation values of test statistics for each hypothesis.
#' @param side Character. Test sidedness: "left", "right", or "both".
#' @param cutoff Numeric. Significance threshold for rejecting hypotheses.
#' @param QC_prob Numeric vector. The probability of failing QC for each hypothesis.
#'
#' @return Numeric. The estimated false discovery proportion (FDP).
#'
#' @details
#' The FDP is computed as:
#' \deqn{FDP = \frac{E[\text{False Discoveries}]}{E[\text{Total Discoveries}]}}
#' 
#' where:
#' \itemize{
#'   \item False discoveries are assumed to follow the null distribution
#'   \item Total discoveries include both true and false positives
#'   \item QC failure probabilities are incorporated into the calculations
#' }
#'
#' @seealso 
#' \code{\link{adjusted_cutoff}} for computing appropriate cutoffs
#' \code{\link{rejection_computation}} for computing rejection probabilities
#' @export
FDP_estimate <- function(mean_list, sd_list, side, cutoff, QC_prob){

  # adjust the number of hypothesis by taking QC probability into consideration
  num_hypo_adjusted <- sum(1 - QC_prob)

  # define the function with cutoff
  rejection_size <- sum(rejection_computation(mean_list = mean_list,
                                              sd_list = sd_list,
                                              side = side,
                                              cutoff = cutoff) * (1 - QC_prob))

  # return the FDP estimate
  return(num_hypo_adjusted * cutoff / rejection_size)
}

#' Compute the rejection probability.
#'
#' @inheritParams adjusted_cutoff
#' @inheritParams compute_power_posthoc
#'
#' @return The rejection probablity.
#' @export

rejection_computation <- function(mean_list, sd_list, side, cutoff){

  # compute different rejection probability based on sideness of the test
  rejection_prob <- switch (side,
                            left = {
                              stats::pnorm(stats::qnorm(cutoff), mean = mean_list, sd = sd_list)
                            },
                            right = {
                              stats::pnorm(stats::qnorm(1 - cutoff),
                                           mean = mean_list, sd = sd_list,
                                           lower.tail = FALSE)
                            },
                            both = {
                              stats::pnorm(stats::qnorm(1 - cutoff / 2), mean = mean_list, sd = sd_list,
                                           lower.tail = FALSE) + stats::pnorm(stats::qnorm(cutoff / 2),
                                                                              mean = mean_list, sd = sd_list)
                            }
  )

  # return the rejection probability list
  return(rejection_prob)
}


# #' Compute the score test statistic.
# #'
# #' @param X Treatment/control indicator.
# #' @param Y Outcome for two groups.
# #' @param size_parameter Size parameter.
# #'
# #' @return Score test statistic.
# #' @export

# score_test <- function(X, Y, size_parameter){
#
#   # compute the number of treat and number of control group
#   n_trt <- sum(X == 1)
#   n_ctl <- sum(X == 0)
#
#   # compute the sample-mean
#   trt_mean <- mean(Y[X == 1])
#   ctl_mean <- mean(Y[X == 0])
#
#   # compute the pooled sample mean
#   pooled_mean <- mean(Y)
#
#   # compute the numerator and denominator of the test statistic
#   test_stat_numerator <- trt_mean - ctl_mean
#   test_stat_denominator <- sqrt(pooled_mean * (1 + pooled_mean / size_parameter)) * sqrt(1 / n_trt + 1 / n_ctl)
#
#   # compute the test statistic
#   test_stat <- test_stat_numerator / test_stat_denominator
#
#   # return the score test statistic
#   return(test_stat)
# }

#' Compute asymptotic distribution of score test statistic
#'
#' @description
#' This function computes the asymptotic mean and standard deviation of the
#' score test statistic used for differential expression analysis in single-cell
#' perturbation experiments.
#'
#' @param num_trt_cells Integer. Number of treatment cells in the score test.
#' @param num_cntrl_cells Integer. Number of control cells in the score test.
#' @param num_trt_cells_sq Numeric. Squared number of treatment cells (used for variance calculations).
#' @param expression_mean Numeric. Mean gene expression level.
#' @param expression_size Numeric. Size parameter in the negative binomial distribution.
#' @param fold_change_mean Numeric. Mean fold change effect size.
#' @param fold_change_sd Numeric. Standard deviation of fold change effect size.
#'
#' @return A list with elements:
#' \describe{
#'   \item{mean}{Asymptotic mean of the test statistic}
#'   \item{sd}{Asymptotic standard deviation of the test statistic}
#' }
#'
#' @details
#' The function computes the asymptotic distribution parameters for a score test
#' statistic under the negative binomial model. The calculations account for:
#' \itemize{
#'   \item Random effect sizes following a normal distribution
#'   \item Negative binomial distribution for gene expression
#'   \item Unequal sample sizes between treatment and control groups
#'   \item Pooled variance estimation for the denominator
#' }
#'
#' The score test statistic follows an asymptotically normal distribution
#' under both null and alternative hypotheses.
#'
#' @seealso 
#' \code{\link{var_nb}} for negative binomial variance calculation
#' \code{\link{rejection_computation}} for computing power from these distributions
#' @export
compute_distribution_teststat <- function(num_trt_cells, num_cntrl_cells, num_trt_cells_sq,
                                          expression_mean, expression_size,
                                          fold_change_mean, fold_change_sd){

  # compute treatment/control cells proportion
  num_test_cells <- num_trt_cells + num_cntrl_cells
  trt_test_prop <- num_trt_cells / num_test_cells
  cntrl_test_prop <- 1 - trt_test_prop

  # define treatment/control/pooled mean expression
  trt_expression_mean <- expression_mean * fold_change_mean
  cntrl_expression_mean <- expression_mean
  pooled_expression_mean <- trt_expression_mean * trt_test_prop + cntrl_expression_mean * cntrl_test_prop

  # compute the square of the denominator in the score statistic
  pooled_var <- var_nb(mean = pooled_expression_mean, size = expression_size)
  denominator_sq <- pooled_var * (1 / num_cntrl_cells + 1 / num_trt_cells)

  ########################## compute the asymptotic sd of test stat ############
  # compute the control group variance
  cntrl_var <- var_nb(mean = cntrl_expression_mean, size = expression_size) / num_cntrl_cells

  # compute the treatment group variance
  trt_var_within <- (trt_expression_mean + (cntrl_expression_mean^2 * (fold_change_sd^2 + fold_change_mean^2) / expression_size)) / num_trt_cells
  trt_var_across <- cntrl_expression_mean^2 * fold_change_sd^2 * num_trt_cells_sq / num_trt_cells^2
  trt_var <- trt_var_within + trt_var_across

  # compute the asymptotic sd
  sd <- sqrt((cntrl_var + trt_var) / denominator_sq)

  ################# compute the asymptotic mean of test stat ###################
  mean <- cntrl_expression_mean * (fold_change_mean - 1) / sqrt(denominator_sq)

  # return the mean and sd vector
  return(
    list(stats::setNames(c(mean, sd), c("mean", "sd")))
  )
}


