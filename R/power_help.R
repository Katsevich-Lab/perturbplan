# This is a Rscript implementing the help function for power/cost analysis

#' Compute the averaged library size per singlet.
#'
#' @param UMI_s Averaged UMI count per singlet (estimate from pilot data given cell-type).
#' @param read_c Averaged read count per cell.
#' @param doublet_rate The doublet fraction among cells.
#' @param doublet_factor Ratio of averaged UMI count per doublet to UMI per singlet.
#'
#' @return Averaged library size per singlet.

library_computation <- function(UMI_s, read_c, doublet_rate, doublet_factor){

  # compute the averaged UMI count per cell
  UMI_c <- UMI_s * (1 + doublet_rate * (doublet_factor - 1))

  # apply the Geometric random variable formula
  library_s <- UMI_s * (1 - exp(- read_c / UMI_c))

  # return the library size value per singlet
  return(library_s)
}

#' Compute averaged number of sequencing reads for each cell after mapping.
#'
#' @param planned_read Number of reads planned before the experiment.
#' @param mapping_efficiency Mapping efficiency for sequenced reads.
#' @param planned_cell Number of cells planned before the experiment.
#' @param recovery_rate Fraction of cells surviving after the library preparation.
#'
#' @return Averaged number of sequencing reads for each cell.

read_per_cell <- function(planned_read, mapping_efficiency,
                          planned_cell, recovery_rate){

  # compute the number of cell surviving the library preparation step
  N_survive <-  planned_cell * recovery_rate

  # compute the number of mapped reads
  mapped_reads <- planned_read * mapping_efficiency

  # return the reads per cell after accounting for mapping and library preparation
  return(mapped_reads / N_survive)
}

#' Compute the gene expression related quantity in the power function.
#'
#' @inheritParams distribution_teststat
#'
#' @return The gene expression related part in the power function.

mean_expression_computation <- function(baseline_expression,
                                        effect_size_mean, control_cell_vec, target_cell_vec){

  # extract information
  num_element <- length(control_cell_vec)
  num_gene <- length(baseline_expression)
  num_cell <- target_cell_vec + control_cell_vec

  # compute the gene expression related part in power formula
  baseline_mat <- matrix(rep(baseline_expression, num_element),
                         nrow = num_element, ncol = num_gene, byrow = TRUE,
                         dimnames = list(
                           gRNA_target = names(control_cell_vec),
                           gene = names(baseline_expression)
                         ))
  trt_mat <- baseline_mat * effect_size_mean
  ctl_mat <- baseline_mat
  pooled_mat <- trt_mat * (target_cell_vec / num_cell) + ctl_mat * (control_cell_vec / num_cell)

  # return the gene expression related part
  return(list(
    baseline_mat = baseline_mat,
    trt_mat = trt_mat,
    ctl_mat = ctl_mat,
    pooled_mat = pooled_mat
  ))
}


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

#' Adjusted power based on adjusted significance level.
#'
#' @param mean_list Asymptotic mean of test statistic (vector; length L x J).
#' @param sd_list Asymptotic sd of test statistic (vector; length L x J).
#' @param QC_prob QC probability for each enhancer-gene pair (vector; length L x J).
#' @inheritParams compute_power
#'
#' @return Adjusted power list including adjusted power and discovery size estimate.

adjusted_power <- function(mean_list, sd_list,
                           multiple_testing_alpha = NULL, multiple_testing_method = NULL, side,
                           QC_prob, cutoff = NULL){

  if(is.null(cutoff)){

    # compute the adjusted cutoff
    cutoff <- adjusted_cutoff(mean_list = mean_list,
                              sd_list = sd_list,
                              multiple_testing_alpha = multiple_testing_alpha,
                              multiple_testing_method = multiple_testing_method,
                              side = side, QC_prob = QC_prob)
  }

  # compute the adjusted power
  adjusted_power <- rejection_computation(mean_list = mean_list,
                                          sd_list = sd_list,
                                          side = side,
                                          cutoff = cutoff)

  # compute the discovery set
  discovery_size <- sum(adjusted_power * (1 - QC_prob))

  # compute the rejection probability with the adjusted cutoff
  return(output = list(
    adjusted_power = adjusted_power,
    discovery_size_estimate = discovery_size
  ))
}

#' Compute the adjusted significance level with either BH or Bonferroni procedure.
#'
#' @inheritParams adjusted_power
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
#' @importFrom dplyr if_else

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

#' FDP estimate based on rejection probability.
#'
#' @inheritParams adjusted_power
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
#' @inheritParams adjusted_power
#'
#' @return The rejection probablity.
#' @importFrom stats qnorm pnorm
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


#' Compute the score test statistic.
#'
#' @param X Treatment/control indicator.
#' @param Y Outcome for two groups.
#' @param size_parameter Size parameter.
#'
#' @return Score test statistic.
#' @export

score_test <- function(X, Y, size_parameter){

  # compute the number of treat and number of control group
  n_trt <- sum(X == 1)
  n_ctl <- sum(X == 0)

  # compute the sample-mean
  trt_mean <- mean(Y[X == 1])
  ctl_mean <- mean(Y[X == 0])

  # compute the pooled sample mean
  pooled_mean <- mean(Y)

  # compute the numerator and denominator of the test statistic
  test_stat_numerator <- trt_mean - ctl_mean
  test_stat_denominator <- sqrt(pooled_mean * (1 + pooled_mean / size_parameter)) * sqrt(1 / n_trt + 1 / n_ctl)

  # compute the test statistic
  test_stat <- test_stat_numerator / test_stat_denominator

  # return the score test statistic
  return(test_stat)
}

#' Compute the mean and sd approximation of the test statistic.
#'
#' @inheritParams compute_power
#' @param target_cell_vec Number of cell per element (vector; length J).
#' @param target_cell_vec_sq Squared number of cell per element (vector; length J).
#'
#' @return Mean and sd vector of length (L x J).
distribution_teststat <- function(control_cell_vec, target_cell_vec, target_cell_vec_sq,
                                  baseline_expression, size_parameter,
                                  effect_size_mean, effect_size_sd){

  # compute the mean gene expression appearing in the power formula
  mean_expression <- mean_expression_computation(baseline_expression = baseline_expression,
                                                 effect_size_mean = effect_size_mean,
                                                 control_cell_vec = control_cell_vec,
                                                 target_cell_vec = target_cell_vec)

  # compute the different variance quantities
  pooled_mean <- mean_expression$pooled_mat
  trt_mean <- mean_expression$trt_mat
  ctl_mean <- mean_expression$ctl_mat

  # compute the number of total elements and total genes
  num_element <- length(control_cell_vec)
  num_gene <- length(size_parameter)

  # fill the size parameter matrix
  size_mat <- matrix(rep(size_parameter, num_element),
                     nrow = num_element, ncol = num_gene, byrow = TRUE,
                     dimnames = list(
                       gRNA_target = names(control_cell_vec),
                       gene = names(size_parameter)
                     ))

  # compute the square of the denominator in the score statistic
  pooled_var <- var_nb(mean = pooled_mean, size = size_mat)
  denominator_sq <- sweep(pooled_var, 1, (1 / control_cell_vec + 1 / target_cell_vec), "*")

  ########################## compute the asymptotic sd of test stat ############
  # compute the control group variance
  ctl_var <- var_nb(mean = ctl_mean, size = size_mat)

  # compute the treatment group variance
  trt_var_part1 <- var_nb(mean = trt_mean, size = size_mat)
  trt_var_part2 <- trt_mean^2 * effect_size_sd^2 / size_mat + sweep(trt_mean^2 * effect_size_sd^2, 1,
                                                                    target_cell_vec_sq / target_cell_vec, "*")
  trt_var <- trt_var_part1 + trt_var_part2

  # compute the asymptotic sd
  asy_var_1 <- sweep(ctl_var, 1, 1 / control_cell_vec, "*")
  asy_var_2 <- sweep(trt_var, 1, 1 / target_cell_vec, "*")
  asy_var <- (asy_var_1 + asy_var_2) / denominator_sq
  asy_sd <- sqrt(asy_var)

  ################# compute the asymptotic mean of test stat ###################
  asy_mean <- ctl_mean * (effect_size_mean - 1) / sqrt(denominator_sq)

  # return the mean and sd vector
  return(list(
    mean = asy_mean,
    sd = asy_sd
  ))
}
