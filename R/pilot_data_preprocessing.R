#' Reference Data Aggregation from Multi-SRR Cell Ranger Outputs
#'
#' @description
#' Aggregates gene expression response matrices and molecular information
#' from multiple SRR-level Cell Ranger output folders, extracting common genes
#' and aligning them into a unified response matrix. Optionally processes only a
#' subset of run-level folders, and summarizes gene consistency across SRRs.
#'
#' @param path_to_top_level_output Character. Path to the top-level directory
#'   containing Cell Ranger run-level subdirectories, each with multiple SRR subdirectories
#' @param path_to_run_level_output Optional character vector. Subset of specific
#'   run-level directories to include. If not found within the top-level path,
#'   unmatched entries will trigger a warning.
#' @param h5_rough Logical. Whether to extract h5 molecular information from only
#'   the first SRR (TRUE) or from all SRRs and combine them (FALSE). Default is TRUE.
#'
#' @return A list with:
#' \describe{
#'   \item{response_matrix}{Combined gene expression matrix with consistent genes across SRRs.}
#'   \item{h5_data}{QC and molecule-level data extracted from HDF5 files.}
#' }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Identifies all run-level and SRR-level directories.
#'   \item Extracts gene expression matrices from each SRR via \code{obtain_qc_response_data}.
#'   \item Intersects gene sets to ensure consistency across all SRRs.
#'   \item Extracts h5 QC data from either one or all SRRs via \code{obtain_qc_h5_data}.
#' }
#'
#' @seealso
#' \code{\link{obtain_qc_response_data}}, \code{\link{obtain_qc_h5_data}}
#'
#' @export
reference_data_preprocessing_10x <- function(path_to_top_level_output,
                                             path_to_run_level_output = NULL,
                                             h5_rough = TRUE
) {
  # Obtain the srr directories from the top-level output directory
  run_dirs <- list.dirs(parent_dir, recursive = FALSE, full.names = TRUE)
  if (!is.null(path_to_run_level_output)) {
    # check if every element in the vector path_to_run_level_output belongs to run_dirs
    # if not, tell the users which elements are unmatched and keep the matched ones in run_dirs
    missing_dirs <- setdiff(path_to_run_level_output, run_dirs)
    if (length(missing_dirs) > 0) {
      warning("The following directories were not found in the run directories: ",
              paste(missing_dirs, collapse = ", "))
    }
    run_dirs <- intersect(run_dirs, path_to_run_level_output)
  }

  srr_dirs <- unlist(lapply(run_dirs, function(rd) {
    list.dirs(rd, recursive = FALSE, full.names = TRUE)
  }), use.names = FALSE)

  # Read the response matrices from each srr directory
  mats <- list()
  for (i in seq_along(srr_dirs)) {
    mats[[i]] <- perturbplan::obtain_qc_response_data(srr_dirs[i])
  }

  all_genes_list <- list()
  for (i in seq_along(mats)) {
    all_genes_list[[i]] <- rownames(mats[[i]])
  }

  common_genes <- Reduce(intersect, all_genes_list)
  # record the genes present in some mats but not in others
  missing_genes <- lapply(mats, function(mat) setdiff(rownames(mat), common_genes))
  if (any(sapply(missing_genes, length) > 0)) {
    warning("Some genes are missing in some matrices: ",
            paste(sapply(missing_genes, function(x) paste(x, collapse = ", ")), collapse = "; "))
  }

  for (i in seq_along(mats)) {
    mats[[i]] <- mats[[i]][common_genes, , drop = FALSE]
  }
  response_matrix <- do.call(cbind, mats)

  # Read the h5 data from each srr directory
  h5_data <- NULL
  if (h5_rough) {
    h5_data <- perturbplan::obtain_qc_h5_data(srr_dirs[1])|>
      mutate(
        srr_idx = srr_dirs[1]
      )
  } else {
    current_data <- perturbplan::obtain_qc_h5_data(srr_dirs[i])|>
      mutate(
        srr_idx = srr_dirs[i]
      )
    h5_data <- dplyr::bind_rows(h5_data, current_data)
  }

  return(list(response_matrix = response_matrix, h5_data = h5_data))
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
#' @param h5_data Data frame. QC information from molecule_info.h5 or filtered_feature_bc_matrix.h5,
#'   as obtained via \code{\link{obtain_qc_h5_data}}.
#' @param n_threads Integer. Number of threads used for parallel processing. Default: NULL (single-threaded).
#' @param downsample_ratio Numeric. Proportion of downsampling used for library size estimation. Default: 0.7.
#' @param D2_rough Numeric. Rough prior value for library variation parameter. Default: 0.3.
#' @param h5_only Logical. If TRUE, skips baseline expression and dispersion estimation
#'   steps (only processes h5_data). Default: FALSE.
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
#' \code{\link{obtain_qc_h5_data}},
#' \code{\link{library_estimation}}
#'
#' @export
reference_data_preprocessing <- function(response_matrix = NULL,
                                         h5_data,
                                         n_threads = NULL,
                                         downsample_ratio = 0.7,
                                         D2_rough = 0.3,
                                         h5_only = FALSE
                                        ) {

  # Input validation
  if (!dir.exists(path_to_cellranger_output)) {
    stop("Cell Ranger output directory does not exist: ", path_to_cellranger_output)
  }

  # Check for required subdirectories and files
  required_paths <- c(
    file.path(path_to_cellranger_output, "outs", "filtered_feature_bc_matrix"),
    file.path(path_to_cellranger_output, "outs", "molecule_info.h5"),
    file.path(path_to_cellranger_output, "outs", "filtered_feature_bc_matrix.h5")
  )

  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths) > 0) {
    stop("Missing required Cell Ranger output files/directories: ",
         paste(missing_paths, collapse = ", "))
  }

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
    QC_data = h5_data,
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
