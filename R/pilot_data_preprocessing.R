#' Aggregate Expression and QC Data from Multiple SRR-Level Cell Ranger Outputs
#'
#' @description
#' This function aggregates gene expression matrices and h5 molecule-level data
#' from multiple SRR-level Cell Ranger output directories. It aligns all matrices
#' to a common set of genes, performs basic consistency checks, and optionally
#' restricts to a user-defined subset of run-level directories.
#'
#' @param path_to_top_level_output Character. Path to the top-level directory
#'   containing Cell Ranger run-level subdirectories.
#' @param path_to_run_level_output Optional character vector. A subset of run-level
#'   directory names (not full paths). These should match the basename of folders
#'   inside \code{path_to_top_level_output}. Unmatched entries will trigger a warning.
#' @param h5_rough Logical. If TRUE (default), h5 data will only be extracted from
#'   the first SRR folder. If FALSE, data from all SRRs will be combined.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{response_matrix}{A matrix of gene expression values (common genes only),
#'     combined across SRR directories.}
#'   \item{read_umi_table}{A data frame of molecule-level QC data from one or more SRRs,
#'     including the SRR label.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Lists all SRR directories under the given top-level folder.
#'   \item Optionally filters to a subset of run-level names.
#'   \item Reads response matrices and retains only shared genes across SRRs.
#'   \item Optionally reads h5 QC data from one or all SRRs.
#' }
#'
#' @seealso \code{\link{obtain_qc_response_data}}, \code{\link{obtain_qc_read_umi_table}}
#' @export
reference_data_preprocessing_10x <- function(path_to_top_level_output,
                                             path_to_run_level_output = NULL,
                                             h5_rough = TRUE) {
  run_dirs <- list.dirs(path_to_top_level_output, recursive = FALSE, full.names = TRUE)
  run_dir_names <- basename(run_dirs)

  if (!is.null(path_to_run_level_output)) {
    missing_dirs <- setdiff(path_to_run_level_output, run_dir_names)
    if (length(missing_dirs) > 0) {
      warning("The following directories were not found in the run directories: ",
              paste(missing_dirs, collapse = ", "))
    }
    run_dirs <- run_dirs[run_dir_names %in% path_to_run_level_output]
    run_dir_names <- basename(run_dirs)  # update names after filtering
  }

  if (length(run_dirs) == 0) {
    stop("No valid run directories found.")
  }

  # Read response matrices
  mats <- lapply(run_dirs, perturbplan::obtain_qc_response_data)
  gene_lists <- lapply(mats, rownames)
  common_genes <- Reduce(intersect, gene_lists)

  if (length(common_genes) == 0) {
    stop("No common genes found across SRRs.")
  }

  # Warn if some genes are missing
  missing_genes <- lapply(mats, function(mat) setdiff(rownames(mat), common_genes))
  if (any(sapply(missing_genes, length) > 0)) {
    warning("Some genes are missing in some matrices: ",
            paste(sapply(missing_genes, function(x) paste(x, collapse = ", ")), collapse = "; "))
  }

  # Subset matrices
  mats <- lapply(mats, function(mat) mat[common_genes, , drop = FALSE])
  response_matrix <- do.call(cbind, mats)

  # Read h5 data
  if (h5_rough) {
    read_umi_table <- perturbplan::obtain_qc_read_umi_table(run_dirs[1]) |>
      dplyr::mutate(srr_idx = run_dir_names[1])
  } else {
    read_umi_table <- list()
    for (i in seq_along(run_dirs)) {
      current <- perturbplan::obtain_qc_read_umi_table(run_dirs[i]) |>
        dplyr::mutate(srr_idx = run_dir_names[i])
      read_umi_table[[i]] <- current
    }
    read_umi_table <- dplyr::bind_rows(read_umi_table)
  }

  return(list(
    response_matrix = response_matrix,
    read_umi_table = read_umi_table
  ))
}




#' Pilot Dataset Preprocessing for Power Analysis
#'
#' @description
#' Further process sequencing data to extract gene-level expression and library
#' statistics required by the PerturbPlan simulation framework. Outputs are
#' compatible with built-in pilot examples.
#'
#' @param response_matrix Matrix or NULL. Gene-by-cell matrix of normalized expression
#'   responses, typically from \code{\link{reference_data_preprocessing_10x}}. If \code{h5_only = TRUE},
#'   this can be NULL.
#' @param read_umi_table Data frame. QC information from molecule_info.h5 or filtered_feature_bc_matrix.h5,
#'   as obtained via \code{\link{obtain_qc_read_umi_table}}.
#' @param n_threads Integer. Number of threads used for parallel processing. Default: NULL (single-threaded).
#' @param downsample_ratio Numeric. Proportion of downsampling used for library size estimation. Default: 0.7.
#' @param D2_rough Numeric. Rough prior value for library variation parameter. Default: 0.3.
#' @param h5_only Logical. If TRUE, skips baseline expression and dispersion estimation
#'   steps (only processes read_umi_table). Default: FALSE.
#'
#' @return A list containing:
#' \describe{
#'   \item{baseline_expression}{List with:
#'     \itemize{
#'       \item \code{baseline_expression}: Data frame with gene-level expression statistics.
#'       \item \code{expression_dispersion_curve}: Function mapping expression to dispersion.
#'     }}
#'   \item{library_parameters}{List with:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Estimated UMI/cell count.
#'       \item \code{variation}: Estimated variation parameter for PCR bias.
#'     }}
#' }
#'
#' @details
#' This function executes the core steps in the pilot data setup for PerturbPlan:
#' \enumerate{
#'   \item Computes gene expression means and variances from response matrix.
#'   \item Fits a nonparametric dispersion curve.
#'   \item Extracts library-level statistics from HDF5 molecule info.
#'   \item Outputs a list structure matching internal pilot examples.
#' }
#'
#' @seealso
#' \code{\link{obtain_expression_information}},
#' \code{\link{obtain_expression_dispersion_curve}},
#' \code{\link{obtain_qc_read_umi_table}},
#' \code{\link{library_estimation}}
#'
#' @export
reference_data_preprocessing <- function(response_matrix = NULL,
                                         read_umi_table,
                                         n_threads = NULL,
                                         downsample_ratio = 0.7,
                                         D2_rough = 0.3,
                                         h5_only = FALSE
                                        ) {

  message("Starting pilot data preprocessing @ ", Sys.time())
  if (!h5_only){
  message("Step 1: Computing gene expression information...")
  baseline_expression_df <- obtain_expression_information(
    response_matrix = response_matrix,
    n_threads = n_threads
  )

  message("Step 2: Computing expression-dispersion curve...")
  expression_dispersion_curve <- obtain_expression_dispersion_curve(baseline_expression_df)
  } else {
    message("Skipping Step 1 and Step 2 as h5_only is TRUE")
    baseline_expression_df <- NULL
    expression_dispersion_curve <- NULL
  }

  message("Step 3: Estimating library parameters...")
  library_params <- library_estimation(
    QC_data = read_umi_table,
    downsample_ratio = downsample_ratio,
    D2_rough = D2_rough
  )

  # Construct the final output structure to match pilot_example.rds format
  result <- list(
    baseline_expression = list(
      baseline_expression = baseline_expression_df,
      expression_dispersion_curve = expression_dispersion_curve
    ),
    library_parameters = library_params  # Already has UMI_per_cell and variation
  )

  message("Completed pilot data preprocessing @ ", Sys.time())
  message("Processed ", nrow(baseline_expression_df), " genes")
  message("Library parameters: UMI_per_cell = ", round(library_params$UMI_per_cell),
          ", variation = ", signif(library_params$variation, 3))

  return(result)
}
