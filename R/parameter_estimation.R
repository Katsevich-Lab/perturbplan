#' Compute the average total UMI per cell and UMI variation parameters.
#'
#' @param QC_data QC'd data coming from the function obtain_qc_data in parameter_estimation_help.R
#' @param library_model S-M curve model fitted by the function model_fit in helper-S-M-curve.R
#'
#' @return Named vector of average total UMI per cell and UMI variation.

library_estimation <- function(QC_data, library_model){

  umi_per_cell <- stats::coef(library_model)["total_UMIs"] / length(unique(QC_data$cell_id))
  umi_variation <- stats::coef(library_model)["D2"] / length(unique(QC_data$cell_id))

  return(
    stats::setNames(c(umi_per_cell, umi_variation), c("umi_per_cell", "umi_variation"))
  )
}

