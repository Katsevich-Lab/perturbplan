#' Compute the average total UMI per cell and UMI variation parameters.
#'
#' @inheritParams library_computation
#'
#' @return Named vector of average total UMI per cell and UMI variation.
#' @export

library_estimation <- function(QC_data, downsample_ratio=0.7, D2_rough=0.3){

  library_model <- library_computation(QC_data, downsample_ratio, D2_rough)
  total_UMIs <- stats::coef(library_model)["total_UMIs"]
  umi_variation <- stats::coef(library_model)["D2"]

  return(
    stats::setNames(c(total_UMIs, umi_variation), c("umi_per_cell", "umi_variation"))
  )
}

