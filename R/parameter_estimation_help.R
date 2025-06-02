#' Obtain gene-level expression and filtering information based on TPM.
#'
#' @param TPM_thres TPM threshold for gene filtering (default = 10).
#' @param response_matrix Optional expression matrix (genes × cells). If NULL, it will be read from file.
#' @param file_type Either "mtx" (default) or "odm", indicating the file format.
#' @param path_to_gene_expression Path to the folder containing expression matrix files.
#' @param cell_covariates A data frame of cell-level covariates (used in precomputation).
#'
#' @return A data frame with genes passing the TPM threshold, with relative expression and size parameter.
obtain_expression_information <- function(TPM_thres = 10, response_matrix = NULL,
                                          file_type = "mtx", path_to_gene_expression,
                                          cell_covariates = NULL) {
  if (is.null(response_matrix)) {
    if (file_type == "odm") {
      gene_odm <- ondisc::read_odm(odm_fp = paste0(path_to_gene_expression, "/matrix.odm"),
                                   metadata_fp = paste0(path_to_gene_expression, "/metadata.rds"))
      ok_cells_gene <- (ondisc::get_cell_covariates(gene_odm)$n_umis) != 0L
      response_matrix <- gene_odm[, ok_cells_gene]
      rownames(response_matrix) <- ondisc::get_feature_ids(gene_odm)
    }

    if (file_type == "mtx") {
      response_matrix <- Matrix::readMM(file.path(path_to_gene_expression, "/matrix.mtx.gz"))
      genes <- data.table::fread(file.path(path_to_gene_expression, "/features.tsv.gz"), header = FALSE)
      barcodes <- data.table::fread(file.path(path_to_gene_expression, "/barcodes.tsv.gz"), header = FALSE)
      rownames(response_matrix) <- genes$V2
      colnames(response_matrix) <- barcodes$V1
    }
  }

  cell_id <- 1:ncol(response_matrix)
  chunk_list <- split(cell_id, cut(seq_along(cell_id), breaks = 10, labels = FALSE))
  sum_expression <- setNames(numeric(nrow(response_matrix)), rownames(response_matrix))

  for (chunk in chunk_list) {
    sum_expression <- Matrix::rowSums(response_matrix[, chunk]) + sum_expression
  }

  relative_expression <- sum_expression / sum(sum_expression)
  TPM <- relative_expression * 1e6

  tpm_filtered_genes <- rownames(response_matrix)[TPM >= TPM_thres]

  gene_info <- data.frame(
    response_id = tpm_filtered_genes,
    stringsAsFactors = FALSE
  ) |> dplyr::mutate(
    relative_expression_at_scale = relative_expression[response_id],
    expression_size_at_scale = unlist(sapply(
      rownames(response_matrix),
      function(response_id) {
        sceptre:::perform_response_precomputation(response_matrix[response_id, ],
                                                  covariate_matrix = cell_covariates)$theta
      }
    ))[response_id]
  )

  return(gene_info)
}

#' Generate random gRNA–gene discovery pairs for control or simulation.
#'
#' @param num_targets Integer. Number of pseudo gRNAs to simulate.
#' @param pairs_per_target Integer. Number of genes to assign per gRNA.
#' @param gene_info Data frame. Must contain a column named `response_id` with gene names.
#'
#' @return A data frame with columns `response_id` and `grna_id`, each row representing a pseudo discovery pair.
#' @export
obtain_random_pairs <- function(num_targets, pairs_per_target, gene_info) {
  # Check input validity
  if (!("response_id" %in% colnames(gene_info))) {
    stop("`gene_info` must contain a column named `response_id`.")
  }

  tpm_filtered_genes <- gene_info$response_id

  if (length(tpm_filtered_genes) < pairs_per_target) {
    stop("Not enough genes to sample for each gRNA. Reduce `pairs_per_target`.")
  }

  discovery_pairs <- NULL

  for (i in seq_len(num_targets)) {
    # randomly sample 'pairs_per_target' genes
    target_sample <- sample(tpm_filtered_genes, pairs_per_target)

    # generate the corresponding gRNA IDs
    grna_sample <- rep(paste0("grna_target_", i), pairs_per_target)

    # append to output
    discovery_pairs <- rbind(
      discovery_pairs,
      data.frame(
        response_id = target_sample,
        grna_id = grna_sample,
        stringsAsFactors = FALSE
      )
    )
  }

  return(discovery_pairs)
}



#' Obtain a data frame with information of all QC'd reads for library estimation.
#'
#' @param path_to_h5_file The path to the outs folder of the cellranger output.
#'
#' @return A data frame with columns `num_reads`, `UMI_id`, `cell_id`, and `response_id`
obtain_qc_h5_data <- function(path_to_h5_file){

  # specify the particular h5 file of interest
  raw_count_file_path <- sprintf("%s/molecule_info.h5", path_to_h5_file)
  qc_info_file_path <- sprintf("%s/filtered_feature_bc_matrix.h5", path_to_h5_file)

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


#' Compute mapping efficiency from QC'd molecule info and Cell Ranger metrics.
#'
#' @param QC_data A data frame from `obtain_qc_h5_data()`, must contain a `num_reads` column.
#' @param path_to_metrics_summary Path to the folder containing Cell Ranger's `metrics_summary.csv`.
#'
#' @return A numeric value representing the mapping efficiency (mapped_reads / total_reads).
#' @export
obtain_mapping_efficiency <- function(QC_data, path_to_metrics_summary) {
  # Check input
  if (!"num_reads" %in% colnames(QC_data)) {
    stop("QC_data must contain a column named `num_reads`.")
  }

  # Read metrics_summary.csv and extract total read count
  metrics_summary <- read.csv(file.path(path_to_metrics_summary, "metrics_summary.csv"),
                              check.names = FALSE)

  if (!"Number of Reads" %in% colnames(metrics_summary)) {
    stop("`metrics_summary.csv` must contain a column named 'Number of Reads'.")
  }

  # Remove commas and convert to numeric
  total_reads <- as.numeric(gsub(",", "", metrics_summary$`Number of Reads`))

  # Compute total mapped reads (after QC)
  mapped_reads <- sum(QC_data$num_reads)

  # Compute mapping efficiency
  mapping_efficiency <- mapped_reads / total_reads

  return(mapping_efficiency)
}


#' Obtain the summary statistics of the QC'd data.
#'
#' @param QC_data The QC'd data coming from the function obtain_qc_h5_data
#'
#' @return A named vector with the number of cells and average reads per cell.
summary_h5_data <- function(QC_data){

  # extract the number of total cells
  num_cells <- length(unique(QC_data$cell_id))

  # extract the number of reads per cell
  total_reads <- sum(QC_data$num_reads)
  num_reads_per_cell <- total_reads / num_cells

  # output the summary statistics
  return(
    stats::setNames(c(num_cells, num_reads_per_cell), c("num_cells", "avg_reads"))
  )
}


#' Fit the S-M curve between # mapped reads and # ovserved UMIs.
#'
#' @param QC_data The QC'd data coming from the function obtain_qc_h5_data.
#' @param downsample_ratio The ratio of the size of the downsampled dataset to the one of the original dataset.
#' @param D2_rough The rough estimate of D2 in the S-M curve model.
#'
#' @return A fitted S-M curve model in the form of a nlsLM object.
library_computation <- function(QC_data, downsample_ratio = 0.7, D2_rough = 0.3){

  ########################### downsample the data ##############################
  # obtain the observed reads vector
  cell_num <- length(unique(QC_data$cell_id))
  num_observed_reads <- sum(QC_data$num_reads)
  reads_vec <- rep(seq_along(QC_data$num_reads), QC_data$num_reads)

  # compute the number of observed UMIs (before downsampling)
  num_observed_umis <- length(unique(reads_vec))

  # perform downsampling and append the results together with observed reads-UMIs
  num_downsampled_reads <- round(num_observed_reads * downsample_ratio)
  num_downsampled_UMIs <- sapply(num_downsampled_reads, function(reads) length(unique(sample(reads_vec, reads))))
  down_sample_added <- data.frame(num_reads = num_downsampled_reads / cell_num, num_UMIs = num_downsampled_UMIs / cell_num)
  down_sample_df <- down_sample_added |>
    dplyr::bind_rows(data.frame(num_reads = num_observed_reads / cell_num, num_UMIs = num_observed_umis / cell_num)) |>
    dplyr::arrange(num_reads, num_UMIs)

  ####################### fit nonlinear model ##################################
  delicate_initial <- (1 + D2_rough) * (num_observed_reads / cell_num)^2 / (2 * (num_observed_reads - num_observed_umis) / cell_num)
  rough_initial <- num_observed_umis / cell_num
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
    relative_loss <- sum((stats::predict(nlm_fitting) / (down_sample_df$num_UMIs) - 1)^2)
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

# nb_estimation <- function(discovery_pairs, response_matrix, covariate_matrix,
#                           grna_matrix = NULL) ) {
#   library(dplyr)
#   library(sceptre)
#   library(ondisc)
#
#   # calculate relative expression
#   cell_id <- 1:ncol(response_matrix)
#   chunk_list <- split(cell_id, cut(seq_along(cell_id), breaks = 10, labels = FALSE))
#   sum_expression <- setNames(numeric(nrow(response_matrix)), rownames(response_matrix))
#   for (chunk in 1:length(chunk_list)) {
#     sum_expression <- apply(response_matrix[, chunk_list[[chunk]]], 1, sum) + sum_expression
#   }
#   relative_expression <- sum_expression / sum(sum_expression)
#
#   # get all parameters for response data
#   gene_summary_info <- data.frame(
#     response_id = rownames(response_matrix)
#   ) |> dplyr::mutate(
#     relative_expression_at_scale = relative_expression[response_id],
#     expression_size_at_scale = unlist(sapply(rownames(response_matrix),
#                                              function(response_id) sceptre:::perform_response_precomputation(response_matrix[response_id, ],
#                                                                               covariate_matrix = cell_covariates)$theta))[response_id])
#   )
#
#   # calculate cells per guide and effect size (if possible)
#   if (!is.null(grna_matrix)) {
#     sceptre_object <- sceptre::import_data(response_matrix = response_matrix,
#                                            grna_matrix = grna_matrix,
#                                            grna_target_data_frame = discovery_pairs,
#                                            moi = "high")
#
#     # set analysis parameters and run qc
#     sceptre_object <- sceptre::set_analysis_parameters(
#       sceptre_object = sceptre_object,
#       discovery_pairs = discovery_pairs,
#       grna_integration_strategy = "singleton",
#       control_group = "complement",
#       resampling_mechanism = "permutations",
#       side = "left",
#       formula_object = formula(~ .)) |>
#       sceptre::assign_grnas(method = "thresholding", threshold = 1) |>
#       sceptre::run_qc(
#         n_nonzero_trt_thresh = 7,
#         n_nonzero_cntrl_thresh = 7,
#         response_n_umis_range = c(0, 1),
#         response_n_nonzero_range = c(0, 1),
#         p_mito_threshold = 1,
#       )
#
#     # run_discovery_analysis for sceptre_object
#     sceptre_object <- sceptre::run_discovery_analysis(sceptre_object)
#
#     # select the variables
#     sceptre_discovery_result <- sceptre_object@discovery_result |>
#       dplyr::filter(!is.na(fold_change)) |>
#       dplyr::select(grna_id, response_id, fold_change, se_fold_change) |>
#       dplyr::rename(fold_change_at_scale = fold_change, se_fold_change_at_scale = se_fold_change) |>
#       dplyr::distinct()
#
#
#     # join the dataframes
#     discovery_pairs <- discovery_pairs |> dplyr::mutate(grna_gene = sprintf("%s_%s", grna_id, response_id))
#     joined_fold_change <- discovery_pairs |>
#       dplyr::left_join(gene_summary_info, "response_id") |>
#       dplyr::mutate(grna_gene = sprintf("%s_%s", grna_id, response_id)) |>
#       dplyr::filter(grna_id != NA, grna_gene %in% discovery_pairs$grna_gene) |>
#       dplyr::left_join(discovery_pairs |> dplyr::select(grna_gene, grna_target, target_type,
#                                                            num_oracle_cells, num_total_plan_cells), by = "grna_gene")
#   }
#   return(joined_fold_change)
# }
