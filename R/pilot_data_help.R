#' @importFrom methods as
#' @importFrom stats setNames
#' @importFrom utils read.csv
NULL

#' @importFrom methods as
#' @importFrom stats setNames coef predict dnbinom
#' @importFrom utils read.csv
#' @importFrom Matrix readMM rowSums colSums
#' @importClassesFrom Matrix CsparseMatrix CsparseMatrix
#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#' @importFrom dplyr filter bind_rows arrange
#' @importFrom preseqR preseqR.ztnb.em
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
#' Only molecules with \code{feature_type == "Gene Expression"} are retained; other
#' feature types (e.g., "Antibody Capture", "CRISPR Guide Capture") are filtered out.
#' This function is used internally by \code{\link{reference_data_preprocessing_10x}}.
#'
#' @param path_to_cellranger_output Character. Path to Cell Ranger run folder
#'   containing:
#'   \itemize{
#'     \item \code{outs/molecule_info.h5} – Raw molecule information with required datasets:
#'       \itemize{
#'         \item \code{count}: Number of reads per molecule
#'         \item \code{umi}: UMI indices
#'         \item \code{barcodes}: Cell barcodes (without GEM group suffix)
#'         \item \code{barcode_idx}: Barcode indices (0-based)
#'         \item \code{gem_group}: GEM group identifiers
#'         \item \code{feature_idx}: Feature indices (0-based)
#'         \item \code{features/id}: Gene identifiers
#'         \item \code{features/feature_type}: Feature types (e.g., "Gene Expression")
#'       }
#'     \item \code{outs/filtered_feature_bc_matrix.h5} – QC-filtered cell barcodes with:
#'       \itemize{
#'         \item \code{matrix/barcodes}: Cell barcodes passing QC filters
#'       }
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
#'   \item Filters to retain only molecules with \code{feature_type == "Gene Expression"}
#'   \item Constructs cell IDs with GEM group suffixes
#'   \item Returns data frame with read counts per UMI per cell for Gene Expression features only
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
#' \strong{Note:} This function only supports Cell Ranger count output format,
#' not Cell Ranger multi.
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
#'   \item \strong{Only Cell Ranger count format is supported.} Cell Ranger multi uses a different
#'     \code{metrics_summary.csv} format (row-based with "Library Type" and "Metric Name" columns)
#'     and is not compatible with this function
#'   \item The \code{metrics_summary.csv} file must contain a column named "Number of Reads"
#'     (Cell Ranger count format where metric names are column headers)
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
  # Call library_computation which now directly returns the parameter list
  return(library_computation(QC_data, downsample_ratio, D2_rough))
}


#' Fit Saturation-Magnitude (S-M) Curve Between Reads and UMIs Using PreseqR
#'
#' @description
#' Fits a zero-truncated negative binomial model using preseqR to estimate the relationship
#' between mapped reads per cell and observed UMIs per cell. The model accounts for both
#' UMI saturation at high read depths and PCR amplification variability. This function is
#' used internally by \code{\link{reference_data_processing}}.
#'
#' @param QC_data Data frame. UMI-level molecule information from
#'   \code{\link{obtain_qc_read_umi_table}} containing columns \code{num_reads},
#'   \code{UMI_id}, \code{cell_id}, and \code{response_id}.
#' @param downsample_ratio Numeric. Not used in preseqR method but kept for
#'   API compatibility. Default: 0.7.
#' @param D2_rough Numeric. Not used in preseqR method but kept for
#'   API compatibility. Default: 0.3.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{UMI_per_cell}{Maximum UMI count per cell at sequencing saturation}
#'   \item{variation}{Variation parameter characterizing UMI richness (1/size from ZTNB model)}
#' }
#'
#' @details
#' ## PreseqR Model
#'
#' The function uses preseqR's zero-truncated negative binomial (ZTNB) model to fit
#' the read-UMI distribution. The saturation curve is:
#'
#' \deqn{\text{UMI} = \text{saturation_UMIs} \times \left(1 - \left(1 + \text{variation} \times \frac{\text{reads}}{\text{saturation_UMIs}}\right)^{-1/\text{variation}}\right)}
#'
#' where:
#' \itemize{
#'   \item \code{reads}: Number of mapped reads per cell (independent variable)
#'   \item \code{UMI}: Number of observed UMIs per cell (dependent variable)
#'   \item \code{saturation_UMIs}: Maximum UMI per cell at saturation
#'   \item \code{variation}: UMI richness variation (1/size parameter from ZTNB fit)
#' }
#'
#' ## Fitting Procedure
#'
#' \enumerate{
#'   \item Creates read-UMI frequency table from QC data
#'   \item Fits ZTNB model using \code{preseqR.ztnb.em()}
#'   \item Extracts size and mu parameters
#'   \item Calculates saturation UMI count per cell
#'   \item Computes UMI richness variation as 1/size
#' }
#'
#' ## Important Notes
#'
#' \itemize{
#'   \item This method replaces the previous downsampling-based approach
#'   \item Parameters downsample_ratio and D2_rough are ignored but retained for compatibility
#'   \item The preseqR method is more robust and doesn't require manual tuning
#'   \item For large datasets, the ZTNB fitting may take several seconds
#' }
#'
#' @examples
#' # Get QC data and compute library parameters
#' cellranger_path <- system.file("extdata/cellranger_tiny", package = "perturbplan")
#' qc_data <- obtain_qc_read_umi_table(cellranger_path)
#'
#' # Fit saturation curve using preseqR
#' lib_params <- library_computation(QC_data = qc_data)
#'
#' # View fitted parameters
#' lib_params$UMI_per_cell
#' lib_params$variation
#'
#' @seealso
#' \code{\link{obtain_qc_read_umi_table}} for input data preparation.
#'
#' \code{\link{reference_data_processing}} for the complete preprocessing workflow.
#'
#' \code{\link{library_estimation}} for the wrapper function.
#' @keywords internal
#' @export
library_computation <- function(QC_data, downsample_ratio = 0.7, D2_rough = 0.3){

  # Create read-UMI frequency table
  read_umi_summary <- QC_data$num_reads |> table()

  # Extract read counts (names) and frequencies (values)
  preseq_input <- cbind(as.integer(names(read_umi_summary)), as.vector(read_umi_summary))

  # Fit ZTNB model using preseqR
  preseq_output <- preseqR::preseqR.ztnb.em(preseq_input)

  # Extract parameters from ZTNB fit
  size <- as.numeric(preseq_output$size)
  mu   <- as.numeric(preseq_output$mu)

  # Calculate summary statistics
  S0 <- sum(preseq_input[, 2])                        # initial distinct UMIs
  R0 <- sum(preseq_input[, 1] * preseq_input[, 2])   # initial total reads
  cell_num <- length(unique(QC_data$cell_id))
  umis_per_cell  <- S0 / cell_num
  reads_per_cell <- R0 / cell_num

  # Calculate saturation parameters
  p0 <- 1 - dnbinom(0, size = size, mu = mu)   # P(seen at least once) at baseline
  saturation_UMIs_per_cell <- umis_per_cell / p0

  # Return parameters in expected format
  return(list(
    UMI_per_cell = saturation_UMIs_per_cell,
    variation = 1 / size
  ))
}
