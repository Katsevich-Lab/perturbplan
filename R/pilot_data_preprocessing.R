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
#'    \item{mapping_efficiency}{A numeric value representing the mapping efficiency}
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

  mapping_efficiency <- NULL
  # Read h5 data to get read umi table and mapping efficiency
  if (h5_rough) {
    read_umi_table <- perturbplan::obtain_qc_read_umi_table(run_dirs[1]) |>
      dplyr::mutate(srr_idx = run_dir_names[1])
    mapping_efficiency <- obtain_mapping_efficiency(read_umi_table,run_dirs[1])
  } else {
    read_umi_table <- list()
    for (i in seq_along(run_dirs)) {
      current <- perturbplan::obtain_qc_read_umi_table(run_dirs[i]) |>
        dplyr::mutate(srr_idx = run_dir_names[i])
      mapping_efficiency[i] <- obtain_mapping_efficiency(current,run_dirs[i])
      read_umi_table[[i]] <- current
    }
    read_umi_table <- dplyr::bind_rows(read_umi_table)
  
    if(length(mapping_efficiency)>1){
      message("Multiple SRR runs detected. Using median mapping efficiency across runs: ",
              paste(round(mapping_efficiency,3), collapse = ", "))
    }
    mapping_efficiency <- median(mapping_efficiency,na.rm=TRUE)
  }

  return(list(
    response_matrix = response_matrix,
    read_umi_table = read_umi_table,
    mapping_efficiency = mapping_efficiency
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
#' @param TPM_thres Numeric. Threshold for filtering low-expression genes during preprocessing.
#' @param downsample_ratio Numeric. Proportion of downsampling used for library size estimation. Default: 0.7.
#' @param D2_rough Numeric. Rough prior value for library variation parameter. Default: 0.3.
#' @param h5_only Logical. If TRUE, skips baseline expression estimation step (only processes read_umi_table). Default: FALSE.
#'  @param mapping_efficiency Numeric. Estimated mapping efficiency from
#'  \code{obtain_mapping_efficiency}.
#'  @param TAP Logical. If TRUE, applies Targeted Analysis Pipeline (TAP)
#'  @param primer_threshold Numeric. Threshold for primer efficiency banding in TAP. Default: 0.2.
#'  @param TAP_target_list Character vector. List of target genes for TAP filtering. If NULL, no filtering is applied.
#'  @param TAP_origin_reference List. Reference data containing baseline expression statistics for TAP.  
#'
#' @return A list containing:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene-level expression statistics including
#'     response_id, relative_expression, and expression_size columns.}
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
#'   \item Extracts library-level statistics from HDF5 molecule info.
#'   \item Outputs a simplified list structure for power analysis.
#' }
#'
#' @seealso
#' \code{\link{obtain_expression_information}},
#' \code{\link{obtain_qc_read_umi_table}},
#' \code{\link{library_estimation}}
#'
#' @export
reference_data_preprocessing <- function(response_matrix = NULL, read_umi_table, mapping_efficiency = NULL,
                                         n_threads = NULL, 
                                         TPM_thres = 0.1, downsample_ratio = 0.7, D2_rough = 0.3, h5_only = FALSE,
                                         TAP = FALSE, primer_threshold = 0.2, TAP_target_list = NULL, TAP_origin_reference = NULL
                                        ) {

  message("Starting pilot data preprocessing @ ", Sys.time())
  if (!h5_only){
  message("Step 1: Computing gene expression information...")
  baseline_expression_df <- obtain_expression_information(
    response_matrix = response_matrix,
    TPM_thres = TPM_thres,  # No filtering during preprocessing
    n_threads = n_threads
  )
  } else {
    message("Skipping Step 1 as h5_only is TRUE")
    baseline_expression_df <- NULL
  }

  message("Step 2: Estimating library parameters...")
  library_params <- library_estimation(
    QC_data = read_umi_table,
    downsample_ratio = downsample_ratio,
    D2_rough = D2_rough
  )
  
  if (TAP) {
    # Correct for primer efficiency and target gene representation in baseline_expression_df
    if (is.null(TAP_target_list)) {
      warning("TAP_target_list must be provided to filter for target genes.")
    } else {
      baseline_expression_df <- baseline_expression_df |>
        dplyr::filter(response_id %in% TAP_target_list)
    }
    if (is.null(TAP_origin_reference)) {
      warning("TAP_origin_reference must be provided to correct for primer efficiency.")
    } else {
      representation_target_genes <- TAP_origin_reference$baseline_expression_stats |> 
                                            dplyr::filter(!is.na(relative_expression)) |>
                                            dplyr::filter(response_id %in% unique(baseline_expression_df$response_id)) |> 
                                            dplyr::summarise(
                                              representation = sum(relative_expression)
                                            ) |> 
                                            dplyr::pull()
      
      relative_expression_df <- baseline_expression_df |> 
        dplyr::distinct() |> 
        dplyr::left_join(
          TAP_origin_reference$baseline_expression_stats |> 
            dplyr:::mutate(Perturb_tpm = relative_expression * 1e6) |>
            dplyr::select(Perturb_tpm, response_id) |>
            dplyr::distinct(), 
          by = "response_id"
        ) |>
        dplyr::mutate(Tap_tpm = relative_expression * representation_target_genes * 1e6) 
      
      # extract the final discovery pairs
      baseline_expression_df <- relative_expression_df |>
        mutate(in_band = between(Tap_tpm, Perturb_tpm * primer_threshold, Perturb_tpm / primer_threshold))
       
      # if any genes are out of the band, print a warning
      if (sum(!baseline_expression_df$in_band)>0){
        message("The following genes are out of the primer efficiency band:")
        cat(baseline_expression_df |> 
                dplyr::filter(!in_band) |>
                dplyr::mutate(expression_status = ifelse(Tap_tpm > Perturb_tpm, "too high", "too low")) |>
                dplyr::select(response_id, Tap_tpm, Perturb_tpm, expression_status))
      }
      # filter the genes that are in the band
      baseline_expression_df <- baseline_expression_df |>  
        dplyr::filter(in_band) |>
        dplyr::select(response_id, relative_expression, expression_size)
    }
    # Normalize relative expression to sum to 1
    baseline_expression_df <- baseline_expression_df |>
      dplyr::mutate(relative_expression = relative_expression / sum(relative_expression)) |>
      dplyr::distinct()
    
    # Correct library parameters
    final_representation_target_genes <- TAP_origin_reference$baseline_expression_stats |> 
      dplyr::filter(!is.na(relative_expression)) |>
      dplyr::filter(response_id %in% unique(baseline_expression_df$response_id)) |> 
      dplyr::summarise(
          representation = sum(relative_expression)
        ) |> 
      dplyr::pull()
    library_params$UMI_per_cell <- library_params$UMI_per_cell * final_representation_target_genes
    
    # Correct mapping efficiency
    mapping_efficiency <- mapping_efficiency * sum(baseline_expression_df$relative_expression)
  }

  # Construct the final output structure with simplified baseline expression
  result <- list(
    baseline_expression_stats = baseline_expression_df,
    library_parameters = library_params,  # Already has UMI_per_cell and variation
    mapping_efficiency = mapping_efficiency
  )

  message("Completed pilot data preprocessing @ ", Sys.time())
  message("Processed ", nrow(baseline_expression_df), " genes")
  message("Library parameters: UMI_per_cell = ", round(library_params$UMI_per_cell),
          ", variation = ", signif(library_params$variation, 3))

  return(result)
}
