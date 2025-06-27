#' Comprehensive pilot data preprocessing from Cell Ranger output
#'
#' @description
#' This function performs all pilot data preprocessing steps to extract baseline 
#' expression information and library parameters from Cell Ranger output. The output
#' structure matches the format expected by the perturbplan power analysis pipeline.
#'
#' @param path_to_cellranger_output Character. Path to Cell Ranger output folder
#'   containing the required subdirectories and files including:
#'   \itemize{
#'     \item \code{outs/filtered_feature_bc_matrix/} (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
#'     \item \code{outs/molecule_info.h5}
#'     \item \code{outs/filtered_feature_bc_matrix.h5}
#'   }
#' @param rough Logical. Whether to use rough dispersion estimation (default: FALSE).
#' @param n_threads Integer. Number of threads for parallel processing (default: NULL).
#' @param downsample_ratio Numeric. Downsampling ratio for library estimation (default: 0.7).
#' @param D2_rough Numeric. Initial D2 parameter estimate (default: 0.3).
#'
#' @return A list with two main elements matching the pilot_example.rds structure:
#' \describe{
#'   \item{baseline_expression}{List containing:
#'     \itemize{
#'       \item \code{baseline_expression}: Data frame with columns \code{response_id}, 
#'             \code{relative_expression}, and \code{expression_size}
#'       \item \code{expression_dispersion_curve}: Function for predicting dispersion 
#'             from expression level
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Total UMI per cell parameter
#'       \item \code{variation}: Variation parameter characterizing PCR bias
#'     }
#'   }
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Load and QC gene expression matrix using \code{\link{obtain_qc_response_data}}
#'   \item Calculate gene expression information and dispersion using \code{\link{obtain_expression_information}}
#'   \item Generate expression-dispersion curve using \code{\link{obtain_expression_dispersion_curve}}
#'   \item Extract molecular information from HDF5 files using \code{\link{obtain_qc_h5_data}}
#'   \item Estimate library parameters using \code{\link{library_estimation}}
#' }
#'
#' The output structure is designed to be compatible with existing perturbplan functions
#' and matches the format of the built-in pilot data examples.
#'
#' @seealso
#' \code{\link{obtain_qc_response_data}}, \code{\link{obtain_expression_information}},
#' \code{\link{obtain_expression_dispersion_curve}}, \code{\link{obtain_qc_h5_data}}, 
#' \code{\link{library_estimation}}
#'
#' @export
pilot_data_preprocessing <- function(path_to_cellranger_output,
                                   rough = FALSE,
                                   n_threads = NULL,
                                   downsample_ratio = 0.7,
                                   D2_rough = 0.3) {

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

  # Step 1: Load and QC gene expression matrix
  message("Step 1: Loading gene expression matrix...")
  response_matrix <- obtain_qc_response_data(path_to_cellranger_output)

  # Step 2: Obtain gene expression information and dispersion
  message("Step 2: Computing gene expression information...")
  baseline_expression_df <- obtain_expression_information(
    response_matrix = response_matrix,
    TPM_thres = 0,  # No filtering during preprocessing
    rough = rough,
    n_threads = n_threads
  )

  # Step 3: Generate expression-dispersion curve
  message("Step 3: Computing expression-dispersion curve...")
  expression_dispersion_curve <- obtain_expression_dispersion_curve(baseline_expression_df)

  # Step 4: Extract QC molecular data
  message("Step 4: Extracting molecular information...")
  QC_data <- obtain_qc_h5_data(path_to_cellranger_output)

  # Step 5: Estimate library parameters
  message("Step 5: Estimating library parameters...")
  library_params <- library_estimation(
    QC_data = QC_data,
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