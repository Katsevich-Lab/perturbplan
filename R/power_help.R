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

#' Compute the adjusted significance level with either BH or Bonferroni procedure.
#'
#' @param mean_list Asymptotic mean of test statistic
#' @param sd_list Asymptotic sd of test statistic
#' @param QC_prob The probability of failing QC
#' @inheritParams compute_power_posthoc
#'
#' @return The adjusted significance level.

adjusted_cutoff <- function(mean_list, sd_list, multiple_testing_alpha, multiple_testing_method, side, QC_prob){

  # compute the adjusted cutoff/significance level
  adjusted_sig_level <- switch (multiple_testing_method,
                                BH = {
                                  BH_cutoff(mean_list = mean_list,
                                            sd_list = sd_list,
                                            multiple_testing_alpha = multiple_testing_alpha,
                                            side = side,
                                            QC_prob = QC_prob)
                                },
                                bonferroni = {
                                  num_hypo_adjusted <- sum(1 - QC_prob)
                                  multiple_testing_alpha / num_hypo_adjusted
                                },
                                BH_cpp = {
                                  BH_cutoff_efficient(mean_list = mean_list,
                                                      sd_list = sd_list,
                                                      multiple_testing_alpha = multiple_testing_alpha,
                                                      side = side,
                                                      QC_prob = QC_prob)
                                },
                                BH_bisection = {
                                  BH_cutoff_bisection(mean_list = mean_list,
                                                      sd_list = sd_list,
                                                      multiple_testing_alpha = multiple_testing_alpha,
                                                      side = side,
                                                      QC_prob = QC_prob)
                                }
  )

  # return the adjusted cutoff/significance level
  return(adjusted_sig_level)
}

#' Compute the adjusted cutoff/significance level applying BH procedure.
#'
#' @inheritParams adjusted_cutoff
#'
#' @return Adjusted cutoff/significance level.

BH_cutoff <- function(mean_list, sd_list, side, multiple_testing_alpha, QC_prob){

  # compute the FDP estimate with the given significance level
  FDP <- function(t){FDP_estimate(mean_list = mean_list,
                                  sd_list = sd_list,
                                  side = side,
                                  cutoff = t, QC_prob = QC_prob)}

  # do a grid search to obtain the adjusted cutoff
  num_hypotheses <- length(mean_list)
  t_vals <- seq(multiple_testing_alpha, multiple_testing_alpha / num_hypotheses,  length.out = num_hypotheses)

  # if the grid search line is too long, we split it
  if(num_hypotheses > 1e4){

    ## check in batches with 1e4 size
    t_hat <- NULL
    tolerance <-  1 / num_hypotheses
    batch <- 0
    splitted_tvals <- split(t_vals, ceiling(seq_along(t_vals) / 1e4))
    while(is.null(t_hat)){

      # add another batch
      batch <- batch + 1

      # extract the current t_vals
      cur_t_vals <- splitted_tvals[[batch]]

      # grid search in next batch
      fdp_hat_vals <- sapply(cur_t_vals, FDP)

      # if the above condition is not true then do the following
      if(all(fdp_hat_vals > multiple_testing_alpha)){
        next
      } else{
        t_hat <- cur_t_vals[min(which(fdp_hat_vals <= multiple_testing_alpha))]
      }

      # if this is the last batch but sill no t_hat searched, specify zero
      if(is.null(t_hat) & (batch == length(splitted_tvals))){
        t_hat <- 0
      }

    }
  }else{

    # do not split
    fdp_hat_vals <- sapply(t_vals, FDP)
    if(all(fdp_hat_vals > multiple_testing_alpha)){
      t_hat <- 0
    }else{
      t_hat <- t_vals[min(which(fdp_hat_vals <= multiple_testing_alpha))]
    }
  }

  # return the adjusted cutoff/significance level
  return(t_hat)
}

#' Benjamini–Hochberg cutoff (C++ back-end)
#'
#' Thin wrapper that validates inputs and forwards to the compiled routine.
#'
#' @inheritParams adjusted_cutoff
#'
#' @return Adjusted cutoff/significance level.
BH_cutoff_efficient <- function(mean_list, sd_list, side, multiple_testing_alpha, QC_prob)
{
  BH_cutoff_cpp(mean_list, sd_list, side, multiple_testing_alpha, QC_prob)
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
  BH_cutoff_bi(mean_list, sd_list, side, multiple_testing_alpha, QC_prob)
}

#' FDP estimate based on rejection probability.
#'
#' @inheritParams adjusted_cutoff
#' @inheritParams compute_power_posthoc
#'
#' @return FDP estimate.

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

#' Compute mean and sd of the score test statistic
#'
#' @inheritParams compute_power_posthoc
#' @param num_trt_cells Number of treatment cells in score test
#' @param num_cntrl_cells Number of control cells in score test
#' @param num_trt_cells_sq Squared number of control cells in score test
#' @param expression_mean Mean gene expression
#' @param expression_size Size parameter in NB distribution
#'
#' @return A list including mean and sd of the test statistic
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


#' Compute the mean and sd of the score test statistic with fixed effect size and random gRNA assignment
#'
#' @inheritParams compute_distribution_teststat_fixed_es
#' @param mean_num_cells Mean of number of treatment cells
#' @param sd_num_cells Standard deviation of number of treatment cells
#' @param B Number of Monte Carlo replications to estiamte variance
#'
#' @return A list including mean and sd of the test statistic
compute_distribution_teststat_fixed_es_random_assignment <- function(
    fold_change,
    expression_mean, expression_size,
    num_cntrl_cells, mean_num_cells, sd_num_cells, B = 1000
){

  # obtian the number of gRNAs
  num_grna <- length(fold_change)

  # obtain the unique value for variables except for fold_change and num_cells
  mean_num_cells <- unique(mean_num_cells)
  expression_mean <- unique(expression_mean)
  expression_size <- unique(expression_size)
  total_num_trt_cells <- mean_num_cells * num_grna
  num_cntrl_cells <- unique(num_cntrl_cells)

  # compute treatment/control cells proportion
  num_test_cells <- total_num_trt_cells + num_cntrl_cells
  trt_test_prop <- total_num_trt_cells / num_test_cells
  cntrl_test_prop <- 1 - trt_test_prop

  # define treatment/control/pooled mean expression
  trt_expression_mean_per_guide <- expression_mean * fold_change
  trt_expression_mean <- mean(trt_expression_mean_per_guide)
  cntrl_expression_mean <- expression_mean
  pooled_expression_mean <- trt_test_prop * trt_expression_mean +  cntrl_test_prop * cntrl_expression_mean

  # compute the square of the denominator in the score statistic
  pooled_var <- var_nb(mean = pooled_expression_mean, size = expression_size)
  denominator_sq <- pooled_var * (total_num_trt_cells / num_cntrl_cells + 1)

  ########################## compute the asymptotic sd of test stat ############
  # compute the across gRNA variance
  simulated_trt_assignment <- stats::rnorm(n = B, mean = total_num_trt_cells, sd = sqrt(num_grna) * sd_num_cells)
  simulated_trt_assignment[simulated_trt_assignment < 0] <- 0
  var_sqrt_trt_assignment <- mean(simulated_trt_assignment) - (mean(sqrt(simulated_trt_assignment)))^2
  across_across_var <- var_sqrt_trt_assignment * (cntrl_expression_mean - trt_expression_mean)^2
  across_within_var <- total_num_trt_cells * (mean(trt_expression_mean_per_guide^2) - mean(trt_expression_mean_per_guide)^2) / 2
  across_var <- across_across_var + across_within_var

  # compute the within gRNA variance
  within_var <- mean(var_nb(mean = trt_expression_mean_per_guide, size = rep(expression_size, num_grna)))

  # compute the asymptotic sd
  sd <- sqrt((across_var + within_var) / denominator_sq)

  ################# compute the asymptotic mean of test stat ###################
  mean <- sqrt(total_num_trt_cells) * (trt_expression_mean - cntrl_expression_mean) / sqrt(denominator_sq)

  # return the mean and sd vector
  return(
    list(stats::setNames(c(mean, sd), c("mean", "sd")))
  )
}
