#' Obtain a data frame with information of all QC'd reads.
#'
#' @param path_to_outs_folder The path to the outs folder of the cellranger output.
#'
#' @return A data frame with columns `num_reads`, `UMI_id`, `cell_id`, and `response_id`
obtain_qc_data <- function(path_to_outs_folder){

  # specify the particular h5 file of interest
  raw_count_file_path <- sprintf("%s/molecule_info.h5", path_to_outs_folder)
  qc_info_file_path <- sprintf("%s/filtered_feature_bc_matrix.h5", path_to_outs_folder)

  ###################### construct the raw data frame ##########################
  raw_count_file <- rhdf5::h5read(raw_count_file_path, "count")
  umi_idx <- rhdf5::h5read(raw_count_file_path, "umi")

  # obtain cell index
  barcode_idx <- rhdf5::h5read(raw_count_file_path, "barcode_idx")
  cell_barcodes <- rhdf5::h5read(raw_count_file_path, "barcodes")
  cell_idx <- cell_barcodes[barcode_idx + 1]

  # append the gem_group
  gem_group <- rhdf5::h5read(raw_count_file_path, "gem_group")
  cell_id_with_gem <- paste(cell_idx, gem_group, sep = "-")

  # obtain gene index for each RNA
  RNA_idx <- rhdf5::h5read(raw_count_file_path, "feature_idx")
  gene_reference <- rhdf5::h5read(raw_count_file_path, "features")
  gene_idx <- gene_reference$id[RNA_idx + 1]

  # store the data frame
  raw_data_frame <- data.frame(
    num_reads = raw_count_file,
    UMI_id = umi_idx + 1,
    cell_id = cell_id_with_gem,
    response_id = gene_idx
  )

  ############################ QC the raw data #################################
  qc_cell <- rhdf5::h5read(qc_info_file_path, "matrix/barcodes")

  # QC the raw data
  qc_df <- raw_data_frame |> dplyr::filter(cell_id %in% qc_cell)

  # return the reads vector
  return(qc_df)
}


#' Obtain the summary statistics of the QC'd data.
#'
#' @param QC_data The QC'd data coming from the function obtain_qc_data
#'
#' @return A named vector with the number of cells and average reads per cell.
summary_data <- function(QC_data){

  # extract the number of total cells
  num_cells <- length(unique(QC_data$cell_id))

  # extract the number of reads per cell
  total_reads <- sum(QC_data$num_reads)
  num_reads_per_cell <- total_reads / num_cells

  # output the summary statistics
  return(
    setNames(c(num_cells, num_reads_per_cell), c("num_cells", "avg_reads"))
  )
}


#' Fit the S-M curve between # mapped reads and # ovserved UMIs.
#'
#' @param QC_data The QC'd data coming from the function obtain_qc_data.
#' @param downsample_ratio The ratio of the size of the downsampled dataset to the one of the original dataset.
#' @param D2_rough The rough estimate of D2 in the S-M curve model.
#'
#' @return A fitted S-M curve model in the form of a nlsLM object.
library_computation <- function(QC_data, downsample_ratio=0.7, D2_rough=0.3){

  ########################### downsample the data ##############################
  # obtain the observed reads vector
  num_observed_reads <- sum(QC_data$num_reads)
  reads_vec <- rep(seq_along(QC_data$num_reads), QC_data$num_reads)

  # compute the number of observed UMIs (before downsampling)
  num_observed_umis <- length(unique(reads_vec))

  # perform downsampling and append the results together with observed reads-UMIs
  num_downsampled_reads <- round(num_observed_reads * downsample_ratio)
  num_downsampled_UMIs <- sapply(num_downsampled_reads, function(reads) length(unique(sample(reads_vec, reads))))
  down_sample_added <- data.frame(num_reads = num_downsampled_reads, num_UMIs = num_downsampled_UMIs)
  down_sample_df <- down_sample_added |>
    dplyr::bind_rows(data.frame(num_reads = num_observed_reads, num_UMIs = num_observed_umis)) |>
    dplyr::arrange(num_reads, num_UMIs)

  ####################### fit nonlinear model ##################################
  delicate_initial <- (1+D2_rough) * num_observed_reads^2 / (2 * (num_observed_reads - num_observed_umis))
  rough_initial <- num_observed_umis
  inital_num_UMIs_vec <- stats::setNames(c(delicate_initial, rough_initial), c("delicate", "rough"))

  # fit model with different initial values on total UMIs
  fitted_output <- lapply(inital_num_UMIs_vec, function(initial_UMIs){

    # do the model fitting
    nlm_fitting <- minpack.lm::nlsLM(
      num_UMIs ~ total_UMIs * (1 - exp(-num_reads / total_UMIs) * (1 + D2 * num_reads^2 / (2 * num_reads^2))),
      data = down_sample_df,
      start = list(total_UMIs = initial_UMIs, D2 = D2_rough),
      upper = c(Inf, 1),
      lower = c(0, 0)
    )

    # return the model and in-sample relative loss
    relative_loss <- sum((stats::predict(nlm_fitting) / QC_data$umi_UMIs - 1)^2)
    output_list <- list(nlm_fitting, relative_loss)
    names(output_list) <- c("fitted_model", "relative_error")
    return(output_list)
  })

  # choose the model with lower relative error
  if(fitted_output$delicate$relative_error > fitted_output$rough$relative_error){
    final_model <- fitted_output$rough$fitted_model
  }else{
    final_model <- fitted_output$rough$fitted_model
  }
  return(final_model)
}

