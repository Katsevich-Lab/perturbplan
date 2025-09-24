#' @importFrom methods as
#' @importFrom stats setNames
#' @importFrom utils read.csv
NULL

#' @importFrom methods as
#' @importFrom stats setNames coef predict
#' @importFrom utils read.csv
#' @importFrom Matrix readMM rowSums colSums
#' @importClassesFrom Matrix CsparseMatrix dgCMatrix
#' @importFrom data.table fread
#' @importFrom dplyr filter bind_rows arrange
#' @importFrom minpack.lm nlsLM
NULL

#' Load and QC gene expression matrix from Cell Ranger output
#'
#' @description
#' Reads a sparse gene-by-cell matrix (\code{matrix.mtx.gz}) together with
#' feature annotations from a Cell Ranger run folder. The required files are
#' located in the sub-directory \code{outs/filtered_feature_bc_matrix/}.
#'
#' @param path_to_cellranger_output Character. Path to an SRR-named folder
#'   (e.g. \code{"SRR12345678"}). This folder must contain:
#'   \itemize{
#'     \item \code{outs/filtered_feature_bc_matrix/matrix.mtx.gz}
#'     \item \code{outs/filtered_feature_bc_matrix/features.tsv.gz}
#'     \item \code{outs/filtered_feature_bc_matrix/barcodes.tsv.gz}
#'   }
#'
#' @examples
#' # Load example Cell Ranger output
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' @return A \code{dgCMatrix} (genes × cells) with cleaned, unique row names and column names.
#' @export
obtain_qc_response_data <- function(path_to_cellranger_output) {
  # Construct path to filtered matrix directory
  mat_dir <- file.path(path_to_cellranger_output, "outs", "filtered_feature_bc_matrix")

  # Read sparse matrix in Matrix Market format
  response_matrix <- as(Matrix::readMM(file.path(mat_dir, "matrix.mtx.gz")), "CsparseMatrix")

  # Read features.tsv.gz: typically contains gene_id (V1), gene_name (V2), and feature_type (V3)
  genes <- data.table::fread(file.path(mat_dir, "features.tsv.gz"), header = FALSE)
  gene_ids <- genes$V1

  # Keep only unique, non-empty gene IDs
  valid_gene <- which(!is.na(gene_ids) & nzchar(gene_ids) & !duplicated(gene_ids))
  response_matrix <- response_matrix[valid_gene, , drop = FALSE]
  rownames(response_matrix) <- gene_ids[valid_gene]

  # Read cell barcodes
  barcodes <- readLines(file.path(mat_dir, "barcodes.tsv.gz"))

  # Set as column names of the matrix
  valid_barcode <- which(!is.na(barcodes) & nzchar(barcodes) & !duplicated(barcodes))
  response_matrix <- response_matrix[, valid_barcode, drop = FALSE]
  colnames(response_matrix) <- barcodes[valid_barcode]

  return(response_matrix)
}

#' Estimate gene-level dispersion (theta)
#'
#' @param response_matrix \code{CsparseMatrix} (genes × cells).
#' @param TPM_thres Numeric. Filter threshold on TPM. Default \code{1}.
#' @param rough Logical. If \code{TRUE}, use rough C++ estimator; otherwise use
#'   refined/MLE. Default \code{FALSE}.
#' @param n_threads Integer controlling parallelism:
#'   \itemize{
#'     \item \code{NULL} – auto-detect (prefer \env{NSLOTS}, else local cores)
#'     \item \code{NA}   – force use of \env{NSLOTS}
#'     \item positive integer – user-supplied core count
#'   }
#' @examples
#' # Get response matrix from Cell Ranger output
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' response_matrix <- obtain_qc_response_data(cellranger_path)
#'
#' # Extract expression information
#' expr_info <- obtain_expression_information(
#'   response_matrix = response_matrix,
#'   TPM_thres = 0.1,
#'   rough = TRUE
#' )
#'
#' # Examine results
#' head(expr_info)
#' dim(expr_info)
#'
#' @return \describe{
#'   \item{response_id}{Gene symbol passing the TPM filter}
#'   \item{relative_expression}{Proportion of total counts}
#'   \item{expression_size}{Estimated dispersion \eqn{\theta}}
#' }
#'
#' @export
obtain_expression_information <- function(response_matrix,
                                          TPM_thres = 0.1,
                                          rough     = FALSE,
                                          n_threads = NULL) {
  # ensure response_matrix is CsparseMatrix
  if (!inherits(response_matrix, "CsparseMatrix")){
    response_matrix <- as(response_matrix, "dgCMatrix")
  }
  # --- decide #threads ------------------------------------------------------
  if (is.null(n_threads)) {
    ns <- Sys.getenv("NSLOTS", unset = "")
    n_threads <- if (nzchar(ns)) as.integer(ns) else parallel::detectCores()
  } else if (is.na(n_threads)) {
    n_threads <- as.integer(Sys.getenv("NSLOTS", unset = "1"))
  } else {
    n_threads <- as.integer(n_threads)
  }
  if (n_threads < 1L) n_threads <- 1L

  # --- library size & TPM ---------------------------------------------------
  message("Start relative expression calculation @ ", Sys.time())
  lib_size <- Matrix::colSums(response_matrix)
  rel_expr <- Matrix::rowSums(response_matrix) / sum(response_matrix)
  names(rel_expr) <- rownames(response_matrix)
  TPM <- rel_expr * 1e6
  keep_gene <- names(TPM)[TPM >= TPM_thres]
  message("Finish relative expression calculation @ ", Sys.time())
  if (!length(keep_gene)) stop("No genes pass TPM threshold")
  # print the number of genes
  message("Number of genes passing TPM threshold: ", length(keep_gene))


  # --- dispersion estimation ------------------------------------------------
  message("Start dispersion estimation (", n_threads, " thread(s)) @ ", Sys.time())
  # check all elements in keep_gene are in rownames(response_matrix)
  if (!all(keep_gene %in% rownames(response_matrix))) {
    # show how many genes are not in it
    missing_genes <- setdiff(keep_gene, rownames(response_matrix))
    message("Missing genes: ", length(missing_genes), " (", paste(utils::head(missing_genes, 10), collapse = ", "), "...)")
    stop("Some genes in keep_gene are not present in response_matrix rownames")
  }
  theta_vec <- theta_batch_cpp(
    Matrix::t(response_matrix[keep_gene, , drop = FALSE]),
    lib_size,
    rel_expr[keep_gene],
    rough      = rough,
    n_threads  = n_threads
  )
  theta_vec[theta_vec == -99] <- NA_real_
  message("Finish dispersion estimation @ ", Sys.time())

  data.frame(
    response_id         = keep_gene,
    relative_expression = rel_expr[keep_gene],
    expression_size     = theta_vec,
    stringsAsFactors    = FALSE
  )
}




#' Extract QC-filtered molecule information from Cell Ranger HDF5 files
#'
#' @param path_to_cellranger_output Character. Folder containing
#'   \code{outs/molecule_info.h5} and \code{outs/filtered_feature_bc_matrix.h5}.
#'
#' @examples
#' # Extract read/UMI information from Cell Ranger output
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' # Examine the data
#' head(qc_table)
#' dim(qc_table)
#' summary(qc_table)
#'
#' @return Data frame with columns \code{num_reads}, \code{UMI_id},
#'   \code{cell_id}, \code{response_id}.
#' @export
obtain_qc_read_umi_table <- function(path_to_cellranger_output) {
  raw_path <- file.path(path_to_cellranger_output, "outs", "molecule_info.h5")
  qc_path  <- file.path(path_to_cellranger_output, "outs", "filtered_feature_bc_matrix.h5")

  count   <- rhdf5::h5read(raw_path, "count")
  umi_idx <- rhdf5::h5read(raw_path, "umi")
  barcodes <- rhdf5::h5read(raw_path, "barcodes")
  barcode_idx <- rhdf5::h5read(raw_path, "barcode_idx")
  gem_group <- rhdf5::h5read(raw_path, "gem_group")
  cell_id <- paste(barcodes[barcode_idx + 1], gem_group, sep = "-")
  RNA_idx <- rhdf5::h5read(raw_path, "feature_idx")
  gene_id <- rhdf5::h5read(raw_path, "features")$id[RNA_idx + 1]

  raw_df <- data.frame(num_reads = count,
                       UMI_id    = umi_idx + 1,
                       cell_id   = cell_id,
                       response_id = gene_id)

  qc_cells <- rhdf5::h5read(qc_path, "matrix/barcodes")
  dplyr::filter(raw_df, cell_id %in% qc_cells)
}

#' Mapping efficiency of a Cell Ranger run
#'
#' @param QC_data Output of \code{obtain_qc_read_umi_table()}.
#' @param path_to_cellranger_output Folder containing `outs/metrics_summary.csv`.
#'
#' @examples
#' # Get mapping efficiency from Cell Ranger output
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' qc_data <- obtain_qc_read_umi_table(cellranger_path)
#' @return Numeric proportion: mapped / total reads.
#' @export
obtain_mapping_efficiency <- function(QC_data, path_to_cellranger_output) {
  if (!"num_reads" %in% names(QC_data))
    stop("QC_data must contain `num_reads`.")
  csv <- file.path(path_to_cellranger_output, "outs", "metrics_summary.csv")
  metrics <- read.csv(csv, check.names = FALSE)
  if (!"Number of Reads" %in% names(metrics))
    stop("Missing 'Number of Reads' in metrics_summary.csv")
  total_reads  <- as.numeric(gsub(",", "", metrics$`Number of Reads`))
  mapped_reads <- sum(QC_data$num_reads)
  mapped_reads / total_reads
}




#' Obtain summary statistics of QC'd molecular data
#'
#' @description
#' This function computes basic summary statistics from quality-controlled
#' molecular data, providing key metrics for library size estimation.
#'
#' @param QC_data Data frame. The QC'd data from \code{\link{obtain_qc_read_umi_table}}
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
#' @examples
#' # Get QC data and summarize
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' qc_data <- obtain_qc_read_umi_table(cellranger_path)
#' @seealso \code{\link{obtain_qc_read_umi_table}} for generating the input data
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
#' @return A list with elements:
#' \describe{
#'   \item{UMI_per_cell}{Total UMI per cell parameter}
#'   \item{variation}{Variation parameter characterizing PCR bias}
#' }
#' @keywords internal

library_estimation <- function(QC_data, downsample_ratio=0.7, D2_rough=0.3){
  library_model <- library_computation(QC_data, downsample_ratio, D2_rough)
  total_UMIs <- stats::coef(library_model)["total_UMIs"]
  umi_variation <- stats::coef(library_model)["D2"]

  return(list(
    UMI_per_cell = unname(as.numeric(total_UMIs)),
    variation = unname(as.numeric(umi_variation))
  ))
}


#' Fit saturation-magnitude (S-M) curve between mapped reads and observed UMIs
#'
#' @description
#' This function fits a nonlinear saturation curve model to estimate the relationship
#' between the number of mapped reads per cell and the number of observed UMIs per cell.
#' The model is essential for library size estimation in single-cell RNA sequencing.
#'
#' @param QC_data Data frame. The QC'd molecular data from \code{\link{obtain_qc_read_umi_table}}
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
#'   \item Issues a warning if the relative error of the fitted model exceeds 7%
#' }
#'
#' @examples
#' # Get QC data and compute library parameters
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' qc_data <- obtain_qc_read_umi_table(cellranger_path)
#'
#' # Compute library size parameters
#' lib_params <- library_computation(
#'   QC_data = qc_data,
#'   downsample_ratio = 0.7,
#'   D2_rough = 0.3
#' )
#'
#' # View parameters
#' print(lib_params)
#'
#' @seealso
#' \code{\link{obtain_qc_read_umi_table}} for input data preparation
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
      num_UMIs ~ total_UMIs * (1 - exp(-num_reads / total_UMIs) * (1 + D2 * num_reads^2 / (2 * total_UMIs^2))),
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
    final_model <- fitted_output$delicate$fitted_model
  }
  # add a warning about the relative error
  if (!is.null(final_model$relative_error) &&
      !is.na(final_model$relative_error) &&
      final_model$relative_error > 0.05) {
    perc_error <- round(100 * final_model$relative_error, 2)
    warning(
      sprintf("The relative error of the fitted model is %.2f%%. Consider adjusting downsample_ratio or D2_rough.", perc_error)
    )
  }
  return(final_model)
}
