#' @importFrom methods as
#' @importFrom stats setNames
#' @importFrom utils read.csv
NULL

#' Load and QC gene expression matrix from Cell Ranger .mtx format
#'
#' @description
#' This function reads a sparse expression matrix from a Cell Ranger output directory
#' (in `.mtx` format) and performs quality control by removing genes with missing,
#' empty, or duplicated names.
#'
#' @param path_to_gene_expression Character. Path to the folder containing
#' `matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz`.
#'
#' @return A sparse gene-by-cell expression matrix of class `dgCMatrix`, where
#' only genes with valid and unique names are retained. Row names are set to gene symbols.
#'
#' @export
obtain_qc_response_data <- function(path_to_gene_expression) {
  # Read sparse matrix (.mtx) and convert to efficient format
  response_matrix <- as(Matrix::readMM(file.path(path_to_gene_expression, "matrix.mtx.gz")), "dgCMatrix")

  # Read features (gene names)
  genes <- data.table::fread(file.path(path_to_gene_expression, "features.tsv.gz"), header = FALSE)
  gene_names <- genes$V2

  # Apply QC: remove empty, NA, or duplicated gene names
  valid_idx <- which(!is.na(gene_names) & gene_names != "" & !duplicated(gene_names))

  # Subset response_matrix and gene_names
  gene_names <- gene_names[valid_idx]
  response_matrix <- response_matrix[valid_idx, , drop = FALSE]

  # Assign cleaned gene names to matrix
  rownames(response_matrix) <- gene_names

  return(response_matrix)
}

#' Estimate Gene Expression Metrics from a Sparse Expression Matrix
#'
#' This function processes a gene-by-cell sparse expression matrix (e.g., read from a `.mtx` file),
#' computes relative expression (TPM), filters genes by a TPM threshold, and estimates the
#' expression-size (theta) for each retained gene using a fast parallel method.
#'
#' @param response_matrix A \code{CsparseMatrix} object representing gene-by-cell expression values.
#'        Usually read with \code{Matrix::readMM()} and coerced via \code{as(., "CsparseMatrix")}.
#' @param TPM_thres A numeric threshold for filtering genes by TPM (Transcripts Per Million).
#'        Genes with TPM below this value will be excluded. Default is 0.
#' @param rough Logical; if \code{TRUE}, use a faster but rougher estimation method for theta.
#'        Default is \code{FALSE}.
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{response_id}{Character vector of gene names passing the TPM threshold.}
#'   \item{relative_expression}{Numeric vector of each gene's relative expression (as a proportion).}
#'   \item{expression_size}{Estimated expression-size parameter \eqn{\theta} for each gene.}
#' }
#'
#' @details
#' The theta parameter is clipped to the interval [0.01, 1e3] to avoid numerical instability.
#'
#' @importMethodsFrom Matrix [
#' @importFrom Matrix rowSums colSums
#' @importClassesFrom Matrix CsparseMatrix dgCMatrix
#' @export
obtain_expression_information <- function(response_matrix,
                                          TPM_thres = 1,
                                          rough     = FALSE) {

  ## 1. library size per cell
  print("Start relative expression")
  library_size <- Matrix::colSums(response_matrix)

  ## 2. gene-level totals & TPM
  gene_sum <- Matrix::rowSums(response_matrix)
  rel_expr <- gene_sum / sum(gene_sum)
  TPM      <- rel_expr * 1e6

  keep_gene <- names(TPM)[TPM >= TPM_thres]
  print("Finish relative expression")
  if (length(keep_gene) == 0)
    stop("No genes pass TPM threshold")

  ## 3. parallel estimation of theta
  print("Start dispersion estimation")
  n_threads <- as.integer(Sys.getenv("NSLOTS", unset = "1"))
  theta_vec <- theta_batch_cpp(
    response_matrix[keep_gene, , drop = FALSE],
    library_size,
    rel_expr[keep_gene],
    rough = rough,
    n_threads = n_threads
  )
  print("Finish dispersion estimation")

  ## 4. assemble result
  data.frame(
    response_id         = keep_gene,
    relative_expression = rel_expr[keep_gene],
    expression_size     = theta_vec,
    stringsAsFactors    = FALSE
  )
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



#' Obtain a data frame with information of all QC'd reads for library estimation
#'
#' @description
#' This function reads HDF5 files from Cell Ranger output to extract molecule-level
#' information and applies quality control filtering to obtain a comprehensive
#' data frame for library size estimation.
#'
#' @param path_to_h5_file Character. The path to the outs folder of the Cell Ranger output
#' containing `molecule_info.h5` and `filtered_feature_bc_matrix.h5` files.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{num_reads}{Number of reads per molecule}
#'   \item{UMI_id}{Unique molecular identifier}
#'   \item{cell_id}{Cell barcode identifier with gem group}
#'   \item{response_id}{Gene identifier}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reads raw molecule information from `molecule_info.h5`
#'   \item Extracts cell barcodes, UMI identifiers, and gene features
#'   \item Filters molecules to only include those from QC-passed cells
#'   \item Returns a cleaned data frame for downstream analysis
#' }
#'
#' @export
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
  metrics_summary <- read.csv(file.path(path_to_metrics_summary, "/metrics_summary.csv"),
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


#' Obtain summary statistics of QC'd molecular data
#'
#' @description
#' This function computes basic summary statistics from quality-controlled
#' molecular data, providing key metrics for library size estimation.
#'
#' @param QC_data Data frame. The QC'd data from \code{\link{obtain_qc_h5_data}}
#' containing columns \code{num_reads}, \code{UMI_id}, \code{cell_id}, and \code{response_id}.
#'
#' @return A named numeric vector with elements:
#' \describe{
#'   \item{num_cells}{Total number of unique cells}
#'   \item{avg_reads}{Average number of reads per cell}
#' }
#'
#' @details
#' The function calculates:
#' \itemize{
#'   \item Total number of unique cell barcodes
#'   \item Total reads summed across all molecules
#'   \item Average reads per cell (total reads divided by number of cells)
#' }
#'
#' @seealso \code{\link{obtain_qc_h5_data}} for generating the input data
#' @export
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


#' Fit saturation-magnitude (S-M) curve between mapped reads and observed UMIs
#'
#' @description
#' This function fits a nonlinear saturation curve model to estimate the relationship
#' between the number of mapped reads per cell and the number of observed UMIs per cell.
#' The model is essential for library size estimation in single-cell RNA sequencing.
#'
#' @param QC_data Data frame. The QC'd molecular data from \code{\link{obtain_qc_h5_data}}
#' containing columns \code{num_reads}, \code{UMI_id}, \code{cell_id}, and \code{response_id}.
#' @param downsample_ratio Numeric. The ratio for downsampling the dataset (default: 0.7).
#' Must be between 0 and 1.
#' @param D2_rough Numeric. Rough estimate of the D2 parameter in the S-M curve model (default: 0.3).
#' Represents the variation parameter in the saturation curve.
#'
#' @return A fitted S-M curve model object of class \code{nlsLM} from the \code{minpack.lm} package.
#' The model parameters include:
#' \describe{
#'   \item{total_UMIs}{Total number of UMIs per cell}
#'   \item{variation}{Variation parameter characterizing PCR bias and saturation}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Downsamples the read data to create multiple observation points
#'   \item Computes UMI counts for different read depths
#'   \item Fits a nonlinear saturation curve: \code{UMI = total_UMIs * (1 - exp(-reads/total_UMIs) * (1 + variation * reads^2/(2*total_UMIs^2)))}
#'   \item Returns the best-fitting model from multiple initial parameter values
#' }
#'
#' @seealso
#' \code{\link{obtain_qc_h5_data}} for input data preparation
#' \code{\link{fit_read_UMI_curve}} for using the fitted parameters
#' @export
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

