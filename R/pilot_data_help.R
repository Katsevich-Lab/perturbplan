#' @importFrom methods as
#' @importFrom stats setNames
#' @importFrom utils read.csv
NULL

#' @importFrom methods as
#' @importFrom stats setNames coef predict
#' @importFrom utils read.csv
#' @importFrom Matrix readMM rowSums colSums
#' @importClassesFrom Matrix CsparseMatrix CsparseMatrix
#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#' @importFrom dplyr filter bind_rows arrange
#' @importFrom minpack.lm nlsLM
NULL

#' Load and QC Gene Expression Matrix from Cell Ranger Output
#'
#' @description
#' Reads a sparse gene-by-cell expression matrix from Cell Ranger output and performs
#' quality control checks. This function is used internally by
#' \code{\link{reference_data_preprocessing_10x}}.
#'
#' @param path_to_cellranger_output Character. Path to a Cell Ranger run folder
#'   (e.g., \code{"SRR12345678"}). This folder must contain:
#'   \itemize{
#'     \item \code{outs/filtered_feature_bc_matrix/matrix.mtx.gz}
#'     \item \code{outs/filtered_feature_bc_matrix/features.tsv.gz}
#'     \item \code{outs/filtered_feature_bc_matrix/barcodes.tsv.gz}
#'   }
#'
#' @return A sparse \code{CsparseMatrix} (genes as rows, cells as columns) with:
#'   \itemize{
#'     \item Unique, non-empty gene IDs as row names
#'     \item Unique, non-empty cell barcodes as column names
#'     \item Duplicate genes and barcodes removed (keeping first occurrence)
#'   }
#'
#' @details
#' In some cases, the subfolder \code{filtered_feature_bc_matrix/} may need to be
#' produced by unzipping the \code{filtered_feature_bc_matrix.tar.gz} file from
#' Cell Ranger output.
#'
#' The function:
#' \enumerate{
#'   \item Reads the sparse matrix in Matrix Market format
#'   \item Converts to column-compressed sparse format (CsparseMatrix)
#'   \item Reads gene annotations from features.tsv.gz
#'   \item Removes duplicate or empty gene IDs
#'   \item Reads cell barcodes from barcodes.tsv.gz
#'   \item Removes duplicate or empty barcodes
#' }
#'
#' @examples
#' # Load example Cell Ranger output
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' response_matrix <- obtain_qc_response_data(cellranger_path)
#'
#' # Inspect the matrix
#' dim(response_matrix)
#' class(response_matrix)
#'
#' @seealso \code{\link{reference_data_preprocessing_10x}} for aggregating data from
#'   multiple Cell Ranger runs
#' @keywords internal
#' @export
obtain_qc_response_data <- function(path_to_cellranger_output) {
  # Construct path to filtered matrix directory
  mat_dir <- file.path(path_to_cellranger_output, "outs", "filtered_feature_bc_matrix")

  # Read sparse matrix in Matrix Market format
  m <- Matrix::readMM(file.path(mat_dir, "matrix.mtx.gz"))

  # If the matrix is in pattern triplet format (ngTMatrix),
  # convert it to numeric triplet (dgTMatrix), then to column-compressed (CsparseMatrix)
  if (inherits(m, "nMatrix")) {
    m <- methods::as(m, "dMatrix")
  }
  if (!inherits(m, "CsparseMatrix")) {
    m <- methods::as(m, "CsparseMatrix")
  }

  response_matrix <- m

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

#' Fit Negative Binomial Model to Estimate Gene Expression Parameters
#'
#' @description
#' Fits a negative binomial model to gene expression data to estimate relative expression
#' levels and dispersion parameters for each gene. This function is used internally by
#' \code{\link{reference_data_processing}}.
#'
#' @param response_matrix Sparse matrix (genes as rows, cells as columns). Typically a
#'   \code{CsparseMatrix} from \code{\link{obtain_qc_response_data}}.
#' @param TPM_thres Numeric. Expression threshold in TPM (Transcripts Per Million) for
#'   filtering low-expression genes. Genes with TPM below this threshold are excluded.
#'   Default: 0.1.
#' @param rough Logical. If TRUE, uses fast C++ estimator for dispersion. If FALSE,
#'   uses refined maximum likelihood estimation. Default: FALSE.
#' @param n_threads Integer or NULL controlling parallelism:
#'   \itemize{
#'     \item \code{NULL} – auto-detect (prefer \env{NSLOTS} environment variable, else
#'       use \code{parallel::detectCores()})
#'     \item \code{NA} – force use of \env{NSLOTS} only
#'     \item positive integer – user-specified thread count
#'   }
#'
#' @return Data frame with three columns:
#' \describe{
#'   \item{response_id}{Gene identifier (e.g., Ensembl ID) for genes passing TPM threshold}
#'   \item{relative_expression}{Estimated relative expression proportion (sums to 1
#'     across all genes)}
#'   \item{expression_size}{Estimated dispersion parameter \eqn{\theta} from negative
#'     binomial model. Small values indicate high biological variability.}
#' }
#'
#' @details
#' ## Negative Binomial Model
#'
#' For each gene, the model is:
#'
#' \deqn{\text{gene_expression} \sim \text{NB}(\text{mean} = \text{library_size} \times \text{relative_expression}, \text{size} = \text{expression_size})}
#'
#' where \code{library_size} is the total UMI count per cell and \code{relative_expression}
#' and \code{expression_size} are the fitted parameters.
#'
#' ## Processing Steps
#'
#' \enumerate{
#'   \item Calculates library sizes (total UMIs per cell)
#'   \item Computes relative expression (gene counts / total counts)
#'   \item Converts to TPM scale and filters genes below threshold
#'   \item Estimates dispersion parameters using C++ implementation
#'   \item Returns data frame with fitted parameters
#' }
#'
#' @examples
#' # Get response matrix from Cell Ranger output
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' response_matrix <- obtain_qc_response_data(cellranger_path)
#'
#' # Extract expression information
#' expr_info <- obtain_expression_information(
#'   response_matrix = response_matrix,
#'   TPM_thres = 0.1,
#'   rough = TRUE,
#'   n_threads = 1
#' )
#'
#' # Examine results
#' head(expr_info)
#' dim(expr_info)
#' summary(expr_info$expression_size)
#'
#' @seealso \code{\link{reference_data_processing}} for the complete pilot data
#'   preprocessing workflow
#' @keywords internal
#' @export
obtain_expression_information <- function(response_matrix,
                                          TPM_thres = 0.1,
                                          rough     = FALSE,
                                          n_threads = NULL) {
  # ensure response_matrix is CsparseMatrix (specifically CsparseMatrix for C++ compatibility)
  if (!inherits(response_matrix, "CsparseMatrix")){
    response_matrix <- as(response_matrix, "CsparseMatrix")
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
  # Ensure transposed matrix is also CsparseMatrix for C++ compatibility
  t_matrix <- Matrix::t(response_matrix[keep_gene, , drop = FALSE])
  if (!inherits(t_matrix, "CsparseMatrix")) {
    t_matrix <- as(t_matrix, "CsparseMatrix")
  }

  theta_vec <- theta_batch_cpp(
    t_matrix,
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




#' Extract UMI-Level Molecule Information from Cell Ranger HDF5 Files
#'
#' @description
#' Extracts QC-filtered UMI-level molecule information from Cell Ranger HDF5 files.
#' This function is used internally by \code{\link{reference_data_preprocessing_10x}}.
#'
#' @param path_to_cellranger_output Character. Path to Cell Ranger run folder
#'   containing:
#'   \itemize{
#'     \item \code{outs/molecule_info.h5} – Raw molecule information
#'     \item \code{outs/filtered_feature_bc_matrix.h5} – QC-filtered cell barcodes
#'   }
#'
#' @return Data frame with UMI-level molecule information containing columns:
#' \describe{
#'   \item{num_reads}{Number of reads supporting this UMI-cell combination}
#'   \item{UMI_id}{UMI index (1-based)}
#'   \item{cell_id}{Cell barcode with GEM group suffix (e.g., "ACGTACGT-1")}
#'   \item{response_id}{Gene identifier (e.g., Ensembl ID)}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Reads raw molecule information from \code{molecule_info.h5}
#'   \item Reads QC-filtered cell barcodes from \code{filtered_feature_bc_matrix.h5}
#'   \item Filters molecule data to retain only QC-passed cells
#'   \item Constructs cell IDs with GEM group suffixes
#'   \item Returns data frame with read counts per UMI per cell
#' }
#'
#' This data is used for fitting the library saturation (S-M) curve in
#' \code{\link{library_computation}}.
#'
#' @examples
#' # Extract read/UMI information from Cell Ranger output
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' qc_table <- obtain_qc_read_umi_table(cellranger_path)
#'
#' # Examine the data
#' head(qc_table)
#' dim(qc_table)
#' summary(qc_table$num_reads)
#'
#' @seealso
#' \code{\link{reference_data_preprocessing_10x}} for aggregating data from multiple runs.
#'
#' \code{\link{library_computation}} for fitting saturation curves using this data.
#' @keywords internal
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

  # Read feature information
  features <- rhdf5::h5read(raw_path, "features")
  gene_id <- features$id[RNA_idx + 1]
  feature_type <- features$feature_type[RNA_idx + 1]

  raw_df <- data.frame(num_reads = count,
                       UMI_id    = umi_idx + 1,
                       cell_id   = cell_id,
                       response_id = gene_id,
                       feature_type = feature_type)

  qc_cells <- rhdf5::h5read(qc_path, "matrix/barcodes")

  # Filter for QC-passed cells and Gene Expression features only and exclude feature_type
  raw_df |> dplyr::filter(cell_id %in% qc_cells & feature_type == "Gene Expression") |>
    dplyr::select(-feature_type)
}

#' Calculate Naive Mapping Efficiency from Cell Ranger Metrics
#'
#' @description
#' Computes the naive mapping efficiency as the proportion of total reads that map
#' to the transcriptome. This function is used internally by
#' \code{\link{reference_data_preprocessing_10x}}.
#'
#' @param QC_data Data frame. Output of \code{\link{obtain_qc_read_umi_table}} containing
#'   a \code{num_reads} column with read counts per UMI.
#' @param path_to_cellranger_output Character. Path to Cell Ranger run folder containing
#'   \code{outs/metrics_summary.csv} with a "Number of Reads" column.
#'
#' @return Numeric value between 0 and 1 representing the proportion of total reads
#'   that successfully mapped to the transcriptome.
#'
#' @details
#' The function calculates:
#'
#' \deqn{\text{mapping_efficiency} = \frac{\text{mapped_reads}}{\text{total_reads}}}
#'
#' where:
#' \itemize{
#'   \item \code{mapped_reads} = sum of \code{num_reads} from QC_data
#'   \item \code{total_reads} = "Number of Reads" from metrics_summary.csv
#' }
#'
#' ## Important Notes
#'
#' \itemize{
#'   \item The \code{metrics_summary.csv} file must contain a column named "Number of Reads"
#'   \item This column may need to be added or edited manually when Cell Ranger is run
#'     with multiple libraries or samples
#'   \item The function removes commas from the "Number of Reads" field before conversion
#'   \item This gives a "naive" estimate that will be adjusted in
#'     \code{\link{reference_data_processing}} when a gene list is specified
#' }
#'
#' @examples
#' # Get mapping efficiency from Cell Ranger output
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' qc_data <- obtain_qc_read_umi_table(cellranger_path)
#' mapping_eff <- obtain_mapping_efficiency(qc_data, cellranger_path)
#'
#' # View result
#' print(mapping_eff)
#'
#' @seealso \code{\link{reference_data_preprocessing_10x}} for the complete aggregation workflow
#' @keywords internal
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
#' @seealso \code{\link{obtain_qc_read_umi_table}} for generating the input data
#' @keywords internal
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


#' Fit Saturation-Magnitude (S-M) Curve Between Reads and UMIs
#'
#' @description
#' Fits a nonlinear saturation curve model to estimate the relationship between mapped
#' reads per cell and observed UMIs per cell. The model accounts for both UMI saturation
#' at high read depths and PCR amplification variability. This function is used internally
#' by \code{\link{reference_data_processing}}.
#'
#' @param QC_data Data frame. UMI-level molecule information from
#'   \code{\link{obtain_qc_read_umi_table}} containing columns \code{num_reads},
#'   \code{UMI_id}, \code{cell_id}, and \code{response_id}.
#' @param downsample_ratio Numeric or numeric vector. Proportion(s) for downsampling
#'   the dataset to create additional observation points. Must be between 0 and 1.
#'   Can be a vector for multiple downsampling levels, but one level is often sufficient.
#'   Default: 0.7.
#' @param D2_rough Numeric. Rough prior estimate for the variation parameter (D2) in
#'   the S-M curve model. Represents PCR amplification bias. Typically 0.3 for perturb-seq,
#'   higher (e.g., 0.8) for TAP-seq. Default: 0.3.
#'
#' @return A fitted S-M curve model object of class \code{nlsLM} from the
#'   \code{minpack.lm} package. The model has two fitted parameters accessible via
#'   \code{coef()}:
#' \describe{
#'   \item{total_UMIs}{Maximum UMI count per cell at sequencing saturation}
#'   \item{D2}{Variation parameter characterizing PCR amplification bias (0 to 1)}
#' }
#'
#' @details
#' ## Saturation Model
#'
#' The S-M curve model is:
#'
#' \deqn{\text{UMI} = \text{total_UMIs} \times \left(1 - \exp\left(-\frac{\text{reads}}{\text{total_UMIs}}\right) \times \left(1 + D2 \times \frac{\text{reads}^2}{2 \times \text{total_UMIs}^2}\right)\right)}
#'
#' where:
#' \itemize{
#'   \item \code{reads}: Number of mapped reads per cell (independent variable)
#'   \item \code{UMI}: Number of observed UMIs per cell (dependent variable)
#'   \item \code{total_UMIs}: Maximum UMI per cell at saturation (fitted parameter)
#'   \item \code{D2}: Variation parameter for PCR bias, between 0 and 1 (fitted parameter)
#' }
#'
#' ## Fitting Procedure
#'
#' \enumerate{
#'   \item Expands read data by replicating UMI indices according to read counts
#'   \item Downsamples the read data at specified ratio(s) to create multiple observation points
#'   \item Counts unique UMIs at each downsampled read depth
#'   \item Fits nonlinear model using two different initial parameter sets:
#'     \itemize{
#'       \item "Delicate": Uses prior D2_rough and derives initial total_UMIs
#'       \item "Rough": Uses observed UMI count as initial total_UMIs
#'     }
#'   \item Selects model with lower relative prediction error
#'   \item Warns if relative error exceeds 5\%
#' }
#'
#' ## Important Notes
#'
#' \itemize{
#'   \item The toy example data has very few reads, so fitted parameters may be sensitive
#'     to random seed and prior specification
#'   \item In practice with real data, the function demonstrates robustness to both random
#'     seed choice and moderate prior misspecification
#'   \item Multiple downsampling ratios can be provided as a vector for more observation
#'     points, but typically one ratio suffices
#' }
#'
#' @examples
#' # Set seed for reproducibility (required for small toy datasets)
#' set.seed(123)
#'
#' # Get QC data and compute library parameters
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' qc_data <- obtain_qc_read_umi_table(cellranger_path)
#'
#' # Fit saturation curve
#' lib_model <- library_computation(
#'   QC_data = qc_data,
#'   downsample_ratio = 0.7,
#'   D2_rough = 0.3
#' )
#'
#' # View fitted parameters
#' coef(lib_model)
#'
#' # Extract specific parameters
#' total_umis <- coef(lib_model)["total_UMIs"]
#' variation <- coef(lib_model)["D2"]
#'
#' @seealso
#' \code{\link{obtain_qc_read_umi_table}} for input data preparation.
#'
#' \code{\link{reference_data_processing}} for the complete preprocessing workflow.
#'
#' \code{\link{library_estimation}} for extracting parameters from the fitted model.
#' @keywords internal
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
