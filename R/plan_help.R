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
#' @param tpm_threshold Numeric. Minimum TPM threshold (default: 10). Genes with expression 
#'   levels below tpm_threshold/1e6 are filtered out before power calculation.
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
#'   \item Filters genes below TPM threshold (relative_expression < tpm_threshold/1e6)
#'   \item Generates fold change values from a normal distribution (one per gene)
#'   \item Returns combined data for Monte Carlo integration
#' }
#'
#' @seealso \code{\link{extract_baseline_expression}} for baseline data extraction
#' @export
extract_fc_expression_info <- function(fold_change_mean, fold_change_sd, biological_system =  "K562", B = 200, gene_list = NULL, tpm_threshold = 10){

  # set the random seed
  set.seed(1)

  ############## combine expression and effect size information ################
  baseline_expression_stats <- extract_baseline_expression(biological_system = biological_system)
  baseline_df <- baseline_expression_stats$baseline_expression
  
  #################### apply TPM threshold filtering FIRST ###################
  # Convert TPM threshold to relative expression scale (TPM / 1e6)
  tpm_threshold_relative <- tpm_threshold / 1e6
  
  # Filter the full baseline dataset by TPM threshold first
  if ("relative_expression" %in% colnames(baseline_df)) {
    pre_filter_n <- nrow(baseline_df)
    filtered_baseline_df <- baseline_df |> 
      dplyr::filter(relative_expression >= tpm_threshold_relative)
    post_filter_n <- nrow(filtered_baseline_df)
    
    # Check if we have any genes left after filtering
    if (post_filter_n == 0) {
      stop("No genes remain after TPM threshold filtering. Consider lowering tpm_threshold.")
    }
    
    # Print filtering summary
    cat("TPM filtering: Kept", post_filter_n, "out of", pre_filter_n, "genes (threshold:", tpm_threshold, "TPM)\n")
  } else {
    warning("No relative_expression column found for TPM filtering. Using unfiltered data.")
    filtered_baseline_df <- baseline_df
    post_filter_n <- nrow(filtered_baseline_df)
  }
  
  ################# sample/subset genes from filtered pool ###################
  # Handle gene-specific vs random sampling scenarios from the filtered pool
  if (!is.null(gene_list)) {
    # User provided specific genes - extract from filtered pool
    if ("response_id" %in% colnames(filtered_baseline_df)) {
      # Filter for specified genes from the TPM-filtered pool
      specified_genes_df <- filtered_baseline_df |> 
        dplyr::filter(response_id %in% gene_list)
      
      # Check if we found any matching genes in the filtered pool
      if (nrow(specified_genes_df) == 0) {
        stop("No matching genes found in TPM-filtered baseline expression data. Genes may be below TPM threshold or not in dataset.")
      } else {
        # Use ALL specified genes that passed TPM filtering
        expression_df <- specified_genes_df
        n_genes <- nrow(expression_df)
        cat("Gene-specific mode: Using", n_genes, "specified genes that passed TPM filtering\n")
      }
    } else {
      stop("Baseline expression data does not contain response_id column. Cannot use specified gene list.")
    }
  } else {
    # No specific genes provided - sample B genes from filtered pool
    if (post_filter_n < B) {
      warning("Fewer genes (", post_filter_n, ") available after TPM filtering than requested (", B, "). Using all available genes.")
      expression_df <- filtered_baseline_df
      n_genes <- post_filter_n
    } else {
      expression_df <- filtered_baseline_df |> dplyr::slice_sample(n = B)
      n_genes <- B
    }
  }
  
  # Combine fold changes with expression parameters
  # Generate fold changes to match the number of genes (after filtering and sampling)
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
#' @export
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
#' This function computes the effective library size (in UMIs) from sequencing read depth
#' using fitted saturation curves that account for PCR amplification bias and UMI saturation.
#'
#' @param reads_per_cell Numeric. Total reads per cell.
#' @param UMI_per_cell Numeric. Maximum UMI per cell parameter from S-M curve fit.
#' @param variation Numeric. Variation parameter characterizing PCR bias from S-M curve fit.
#'
#' @return Numeric. Effective library size in UMIs.
#'
#' @details
#' The saturation-magnitude (S-M) curve model relates sequencing reads to unique UMI counts
#' accounting for:
#' \itemize{
#'   \item PCR amplification variability
#'   \item UMI saturation at high read depths
#'   \item Platform-specific technical biases
#' }
#'
#' The relationship follows: UMI = UMI_per_cell * (1 - exp(-reads_per_cell / (UMI_per_cell * variation)))
#'
#' @seealso \code{\link{extract_library_info}} for obtaining curve parameters
#' @export
fit_read_UMI_curve <- function(reads_per_cell, UMI_per_cell, variation){

  # Apply the saturation curve transformation
  effective_UMI <- UMI_per_cell * (1 - exp(-reads_per_cell / (UMI_per_cell * variation)))
  
  return(effective_UMI)
}