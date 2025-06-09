
# Suppress R CMD check warnings for variables used in dplyr contexts
utils::globalVariables(c("response_id"))

#' Extract fold change and expression information for power analysis
#'
#' @description
#' This function combines fold change effect size sampling with baseline expression
#' data to create a comprehensive dataset for Monte Carlo power analysis simulations.
#' It can handle both user-specified genes and random sampling scenarios.
#'
#' @param fold_change_mean Numeric. Mean of the fold change distribution.
#' @param fold_change_sd Numeric. Standard deviation of the fold change distribution.
#' @param biological_system Character. Biological system for baseline expression (default: "K562").
#' @param B Integer. Number of Monte Carlo samples to generate when gene_list is NULL (default: 200).
#'   Ignored when gene_list is provided.
#' @param gene_list Character vector. Optional list of Ensembl gene IDs to use for analysis.
#'   If provided, expression parameters will be extracted for ALL specified genes (no sampling).
#'   If NULL (default), B genes are randomly sampled from baseline data.
#'
#' @return A list with elements:
#' \describe{
#'   \item{fc_expression_df}{Data frame with sampled fold changes and expression parameters}
#'   \item{expression_dispersion_curve}{Function relating expression mean to dispersion}
#' }
#'
#' @details
#' The function operates in two modes:
#' \itemize{
#'   \item \strong{Gene-specific mode} (gene_list provided): Uses ALL specified genes, no sampling
#'   \item \strong{Random sampling mode} (gene_list = NULL): Randomly samples B genes from baseline
#' }
#' 
#' In both modes:
#' \itemize{
#'   \item Sets a random seed for reproducibility
#'   \item Generates fold change values from a normal distribution (one per gene)
#'   \item Returns combined data for Monte Carlo integration
#' }
#'
#' @seealso \code{\link{extract_baseline_expression}} for baseline data extraction
#' @export
extract_fc_expression_info <- function(fold_change_mean, fold_change_sd, biological_system =  "K562", B = 200, gene_list = NULL){

  # set the random seed
  set.seed(1)

  ############## combine expression and effect size information ################
  baseline_expression_stats <- extract_baseline_expression(biological_system = biological_system)
  
  # Handle gene-specific vs random sampling scenarios
  if (!is.null(gene_list)) {
    # User provided specific genes - use ALL specified genes (no sampling)
    baseline_df <- baseline_expression_stats$baseline_expression
    
    # Check if baseline data has gene identifiers (assuming 'response_id' column exists)
    if ("response_id" %in% colnames(baseline_df)) {
      # Filter for specified genes
      specified_genes_df <- baseline_df |> 
        dplyr::filter(response_id %in% gene_list)
      
      # Check if we found any matching genes
      if (nrow(specified_genes_df) == 0) {
        stop("No matching genes found in baseline expression data. Please check gene IDs.")
      } else {
        # Use ALL specified genes (no sampling)
        expression_df <- specified_genes_df
        n_genes <- nrow(expression_df)
      }
    } else {
      stop("Baseline expression data does not contain response_id column. Cannot use specified gene list.")
    }
  } else {
    # No specific genes provided - use random sampling of B genes
    expression_df <- baseline_expression_stats$baseline_expression |> dplyr::slice_sample(n = B)
    n_genes <- B
  }
  
  # Combine fold changes with expression parameters
  # Generate fold changes to match the number of genes (either B or length of gene_list)
  fc_expression_df <- data.frame(
    fold_change = stats::rnorm(n = n_genes, mean = fold_change_mean, sd = fold_change_sd)
  ) |>
    dplyr::bind_cols(expression_df)

  ################## extract the expression-dispersion curve ###################
  expression_dispersion_curve <- baseline_expression_stats$expression_dispersion_curve

  # return the data frame
  return(list(
    fc_expression_df = fc_expression_df,
    expression_dispersion_curve = expression_dispersion_curve
  ))
}

#' Extract baseline expression data by biological system
#'
#' @description
#' This function loads pre-computed baseline expression statistics and
#' expression-dispersion relationships for specified biological systems.
#'
#' @param biological_system Character. Biological system identifier (default: "K562").
#' Currently supports "K562" cells.
#'
#' @return A list with elements:
#' \describe{
#'   \item{baseline_expression}{Data frame with expression parameters for genes}
#'   \item{expression_dispersion_curve}{Function mapping expression mean to dispersion}
#' }
#'
#' @details
#' The function loads system-specific data files containing:
#' \itemize{
#'   \item Gene expression means and dispersion parameters
#'   \item Fitted dispersion curves for the negative binomial model
#'   \item System-specific calibration parameters
#' }
#'
#' @seealso \code{\link{extract_fc_expression_info}} for combining with effect sizes
extract_baseline_expression <- function(biological_system = "K562"){

  # sample baseline expression based on biological system
  switch(biological_system,
         K562 = {

           # load the Gasperini baseline expression list
           rds_path <- system.file("extdata/baseline_expression", "Gasperini_expression.rds", package = "perturbplan", mustWork = TRUE)
           baseline_expression_list <- readRDS(rds_path)

         })

  # return the data frame with the susbampled rows
  return(list(
    baseline_expression = baseline_expression_list$baseline_expression,
    expression_dispersion_curve = baseline_expression_list$expression_dispersion_curve
  ))
}

#' Extract library size parameters by biological system
#'
#' @description
#' This function loads pre-computed library size parameters including UMI saturation
#' curve parameters for specified biological systems.
#'
#' @param biological_system Character. Biological system identifier (default: "K562").
#' Currently supports "K562" cells.
#'
#' @return A list with elements:
#' \describe{
#'   \item{UMI_per_cell}{Total UMI per cell parameter}
#'   \item{variation}{Variation parameter characterizing PCR bias}
#' }
#'
#' @details
#' The function loads system-specific library parameters fitted from
#' single-cell RNA-seq data using the saturation-magnitude (S-M) curve model.
#' These parameters are used in \code{\link{fit_read_UMI_curve}} to convert
#' read depth to effective library size.
#'
#' @seealso 
#' \code{\link{fit_read_UMI_curve}} for using these parameters
#' \code{\link{library_computation}} for fitting these parameters from data
extract_library_info <- function(biological_system = "K562"){

  # sample baseline expression based on biological system
  switch(biological_system,
         K562 = {

           # load the Gasperini baseline expression list
           rds_path <- system.file("extdata/library_info", "Gasperini_library.rds", package = "perturbplan", mustWork = TRUE)
           library_info <- readRDS(rds_path)

         })

  # return the data frame with the susbampled rows
  return(list(
    UMI_per_cell = library_info$S_M_curve_params["UMI_per_cell"],
    variation = library_info$S_M_curve_params["variation"]
  ))
}

#' Compute effective library size from read depth using UMI saturation curve
#'
#' @description
#' This function implements the saturation-magnitude (S-M) curve to convert
#' read depth per cell into effective library size (observed UMIs), accounting
#' for PCR bias and UMI saturation effects.
#'
#' @param reads_per_cell Numeric. Number of reads per cell.
#' @param UMI_per_cell Numeric. Total UMI per cell parameter (asymptotic maximum).
#' @param variation Numeric. Variation parameter characterizing PCR bias and saturation.
#'
#' @return Numeric. Effective library size (number of observed UMIs per cell).
#'
#' @details
#' The function implements the S-M curve model:
#' \deqn{UMI_{obs} = UMI_{total} \times \left(1 - \exp\left(-\frac{reads}{UMI_{total}}\right) \times \left(1 + variation \times \frac{reads^2}{2 \times UMI_{total}^2}\right)\right)}
#'
#' where:
#' \itemize{
#'   \item \code{UMI_obs} is the number of observed UMIs (library size)
#'   \item \code{UMI_total} is the total number of UMIs per cell
#'   \item \code{reads} is the number of reads per cell
#'   \item \code{variation} accounts for PCR bias and technical variation
#' }
#'
#' @seealso 
#' \code{\link{extract_library_info}} for obtaining the parameters
#' \code{\link{library_computation}} for fitting parameters from data
#' @export
fit_read_UMI_curve <- function(reads_per_cell, UMI_per_cell, variation){
  UMI_per_cell * (1 - exp(-reads_per_cell / UMI_per_cell) * (1 + variation * reads_per_cell^2 / (2 * UMI_per_cell^2)))
}
