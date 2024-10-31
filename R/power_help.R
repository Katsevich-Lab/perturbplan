# This is a Rscript implementing the help function for power/cost analysis

#' Compute the averaged library size per singlet
#'
#' @param UMI_s Averaged UMI count per singlet (estimate from pilot data given cell-type)
#' @param read_c Averaged read count per cell
#' @param doublet_rate The doublet fraction among cells
#' @param doublet_factor Ratio of averaged UMI count per doublet to UMI per singlet
#'
#' @return Averaged library size per singlet
#' @export

library_computation <- function(UMI_s, read_c, doublet_rate, doublet_factor){
  
  # compute the averaged UMI count per cell
  UMI_c <- UMI_s * (1 + doublet_rate * (doublet_factor - 1))
  
  # apply the Geometric random variable formula
  library_s <- UMI_s * (1 - exp(- read_c / UMI_c))
  
  # return the library size value per singlet
  return(library_s)
}

#' Compute the weighted gRNA efficiency by the cell size
#'
#' @param pi_mat gRNA efficiency matrix of dim L by K (row: genomic elements, column: different gRNA)
#' @param cell_mat Cell size matrix of dim L by K (row: genomic elements, column: different gRNA)
#'
#' @return A vector of weighted gRNA efficiency vector by cell size (of length L)
#' @export

efficiency_computation <- function(pi_mat, cell_mat){
  
  # convert the cell_mat to a weight matrix
  weight_mat <- cell_mat / rowSums(cell_mat)
  
  # compute the weighted efficiency
  weighted_pi_mat <- weight_mat * pi_mat
  
  # return the efficiency vector
  return(rowSums(weighted_pi_mat))
}

#' Compute averaged number of sequencing reads for each cell after mapping
#'
#' @param planned_read Number of reads planned before the experiment
#' @param mapping_efficiency Mapping efficiency for sequenced reads
#' @param planned_cell Number of cells planned before the experiment
#' @param recovery_rate Fraction of cells surviving after the library preparation
#'
#' @return Averaged number of sequencing reads for each cell
#' @export

read_per_cell <- function(planned_read, mapping_efficiency,
                          planned_cell, recovery_rate){
  
  # compute the number of cell surviving the library preparation step
  N_survive <-  planned_cell * recovery_rate
  
  # compute the number of mapped reads
  mapped_reads <- planned_read * mapping_efficiency
  
  # return the reads per cell after accounting for mapping and library preparation
  return(mapped_reads / N_survive)
}

#' Compute the gene expression related quantity in the power function
#'
#' @param expression_level_list List of relative gene expression level
#' @param dispersion_list Dispersion parameter for each gene
#' @param library_size Library size for each cell considered
#'
#' @return The gene expression related part in the power function
#' @export

gene_part_computation <- function(expression_level_list, dispersion_list, library_size){
  
  # compute the gene expression related part in power formula
  baseline_expression <- library_size * expression_level_list
  gene_part <- sqrt(baseline_expression / (1 + baseline_expression / dispersion_list))
  
  # return the gene expression related part
  return(gene_part)
}

#' Adjusted power based on adjusted significance level
#'
#' @param mean_list List of mean under local alternative
#' @param sig_level Significance level imposed by users
#' @param correction Either BH or Bonferroni 
#' @param sideness Sideness of the testing procedure
#'
#' @return Adjusted power list including adjusted power and discovery size estimate
#' @export

adjusted_power <- function(mean_list, sig_level, correction, sideness){
  
  # compute the adjusted cutoff
  adjusted_cutoff <- adjusted_cutoff(mean_list = mean_list, 
                                     sig_level = sig_level, 
                                     correction = correction, 
                                     sideness = sideness)
  
  # compute the adjusted power
  adjusted_power <- rejection_computation(mean_list = mean_list, 
                                          sideness = sideness,
                                          sig_level = adjusted_cutoff)
  
  # compute the discovery set
  discovery_size <- sum(adjusted_power)
    
  # compute the rejection probability with the adjusted cutoff
  return(output = list(
    adjusted_power = adjusted_power,
    discovery_size_estimate = discovery_size
  ))
}

#' Compute the adjusted significance level with either BH or Bonferroni procedure
#'
#' @param mean_list List of mean under local alternative
#' @param sig_level Significance level imposed by users
#' @param correction Either BH or Bonferroni 
#' @param sideness Sideness of the testing procedure
#'
#' @return The adjusted significance level
#' @export

adjusted_cutoff <- function(mean_list, sig_level, correction, sideness){
  
  # extract the number of hypotheses from the mean_list
  num_hypotheses <- length(mean_list)
  
  # compute the adjusted cutoff/significance level
  adjusted_sig_level <- switch (correction,
                                BH = {
                                  BH_cutoff(mean_list = mean_list, 
                                            sig_level = sig_level, 
                                            sideness = sideness)
                                },
                                Bonferroni = {
                                  sig_level / num_hypotheses
                                }
  )
  
  # return the adjusted cutoff/significance level
  return(adjusted_sig_level)
}

#' Compute the adjusted cutoff/significance level applying BH procedure
#'
#' @param mean_list List of mean under local alternative
#' @param sideness Sideness of the testing procedure
#' @param sig_level Significance level imposed by users
#'
#' @return Adjusted cutoff/significance level
#' @importFrom stats uniroot
#' @export

BH_cutoff <- function(mean_list, sideness, sig_level){
  
  # compute the FDP estimate with the given significance level
  FDP <- function(t){FDP_estimate(mean_list = mean_list,
                                  sideness = sideness,
                                  sig_level = t)}
  
  # do a grid search to obtain the adjusted cutoff
  num_hypotheses <- length(mean_list)
  t_vals <- seq(sig_level / num_hypotheses, sig_level, length.out = num_hypotheses)
  fdp_hat_vals <- sapply(t_vals, FDP)
  if(all(fdp_hat_vals > sig_level)){
    t_hat <- 0
  } else{
    t_hat <- t_vals[max(which(fdp_hat_vals <= sig_level))]
  }
  
  # return the adjusted cutoff/significance level
  return(t_hat)
}

#' FDP estimate based on rejection probability
#'
#' @param mean_list List of mean under local alternative
#' @param sideness Sideness of the testing procedure
#' @param sig_level Significance level imposed by users
#'
#' @return FDP estimate
#' @export

FDP_estimate <- function(mean_list, sideness, sig_level){
  
  # extract the number of hypotheses considered
  num_hypotheses <- length(mean_list)
  
  # define the function with cutoff
  rejection_size <- sum(rejection_computation(mean_list = mean_list,
                                              sideness = sideness,
                                              sig_level = sig_level))
  
  # return the FDP estimate 
  return(num_hypotheses * sig_level / rejection_size)
}

#' Compute the rejection probability
#'
#' @param mean_list List of mean under local alternative
#' @param sideness Sideness of the testing procedure
#' @param sig_level Significance level imposed by users
#'
#' @return The rejection probablity
#' @importFrom stats qnorm pnorm
#' @export

rejection_computation <- function(mean_list, sideness, sig_level){
  
  # compute different rejection probability based on sideness 
  rejection_prob <- switch (sideness,
                            left = {
                              stats::pnorm(stats::qnorm(sig_level), mean = mean_list)
                            },
                            right = {
                              stats::pnorm(stats::qnorm(1 - sig_level), mean = mean_list, lower.tail = FALSE)
                            },
                            both = {
                              stats::pnorm(stats::qnorm(1 - sig_level / 2), mean = mean_list,  
                                           lower.tail = FALSE) + stats::pnorm(stats::qnorm(sig_level / 2), mean = mean_list)
                            }
  )
  
  # return the rejection probability list
  return(rejection_prob)
}