#' Compute the averaged library size per singlet.
#'
#' @param reads_per_cell Averaged read count per cell.
#' @param h5_path H5 file path to be used
#'
#' @return Averaged library size per singleton.

library_computation <- function(reads_per_cell, h5_path){

  # compute the averaged UMI count per cell
  UMIs_per_cell <-

  # return the library size value per singlet
  return(UMIs_per_cell)
}

#' Compute averaged number of reads after QC step.
#'
#' @param num_total_reads Number of reads planned before the experiment.
#' @param mapping_efficiency Mapping efficiency for sequenced reads.
#' @param num_total_cells Number of cells planned after the library preparation.
#'
#' @return Averaged number of reads for each cell (singleton).

read_per_cell_after_QC <- function(num_total_reads, mapping_efficiency,
                                   num_total_cells){

  # compute the number of mapped reads
  mapped_reads <- num_total_reads * mapping_efficiency

  # return the reads per cell after accounting for mapping and library preparation
  return(mapped_reads / num_total_cells)
}
