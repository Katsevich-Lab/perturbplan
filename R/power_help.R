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
