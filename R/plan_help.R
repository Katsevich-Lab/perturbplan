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
#' @param custom_baseline_data List. Optional custom baseline expression data with same structure
#'   as extract_baseline_expression() output. If provided, this data is used instead of the
#'   default biological_system data. Must contain baseline_expression data frame and
#'   expression_dispersion_curve function.
#' @param gRNAs_per_target Integer. Number of gRNAs per target (default: 4).
#'   Each target will have gRNAs_per_target individual gRNA effect sizes drawn from the
#'   specified fold change distribution. avg_fold_change and avg_fold_change_sq are
#'   calculated as the mean and mean-of-squares of these gRNA effect sizes.
#'
#' @return A list with elements:
#' \describe{
#'   \item{fc_expression_df}{Data frame with avg_fold_change, avg_fold_change_sq, and expression parameters}
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
#'   \item Generates gRNAs_per_target effect sizes per target from a normal distribution
#'   \item Calculates avg_fold_change and avg_fold_change_sq from the gRNA effect sizes
#'   \item Returns combined data for Monte Carlo integration with random effect sizes
#' }
#'
#' @seealso \code{\link{extract_baseline_expression}} for baseline data extraction
#' @export
extract_fc_expression_info <- function(fold_change_mean, fold_change_sd, biological_system =  "K562", B = 200, gene_list = NULL, tpm_threshold = 10, custom_baseline_data = NULL, gRNAs_per_target = 4){

  # set the random seed
  set.seed(1)

  ############## combine expression and effect size information ################
  # Use custom baseline data if provided, otherwise use default biological system data
  if (!is.null(custom_baseline_data)) {
    baseline_expression_stats <- custom_baseline_data
  } else {
    baseline_expression_stats <- extract_baseline_expression(biological_system = biological_system)
  }
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

  ################# sample B genes from filtered pool with importance weights ###################
  # Always sample B genes, but use importance sampling when gene_list is provided
  if (!is.null(gene_list)) {
    # User provided specific genes - use importance sampling based on gene weights
    if ("response_id" %in% colnames(filtered_baseline_df)) {

      # Calculate gene weights from the original gene list
      gene_weights <- table(gene_list)
      unique_genes <- names(gene_weights)

      # Filter baseline data for unique genes that passed TPM filtering
      unique_genes_df <- filtered_baseline_df |>
        dplyr::filter(response_id %in% unique_genes)

      # Check if we found any matching genes in the filtered pool
      if (nrow(unique_genes_df) == 0) {
        stop("No matching genes found in TPM-filtered baseline expression data. Genes may be below TPM threshold or not in dataset.")
      }

      # Check which requested genes were not found
      found_genes <- unique_genes_df$response_id
      missing_genes <- setdiff(unique_genes, found_genes)
      if (length(missing_genes) > 0) {
        warning("Some requested genes were filtered out: ", paste(missing_genes, collapse = ", "))
      }

      # Use importance sampling: sample B genes based on weights from found genes
      weights <- as.numeric(gene_weights[found_genes])

      # Sample B genes with replacement according to importance weights
      sampled_indices <- sample(seq_len(nrow(unique_genes_df)),
                              size = B,
                              prob = weights,
                              replace = TRUE)

      expression_df <- unique_genes_df[sampled_indices, ]
      n_genes <- B
      cat("Gene-specific mode with importance sampling: Using", length(found_genes), "unique genes, sampled", B, "times with weights\n")

    } else {
      stop("Baseline expression data does not contain response_id column. Cannot use specified gene list.")
    }
  } else {
    # No specific genes provided - sample B genes uniformly from filtered pool with replacement
    sampled_indices <- sample(seq_len(nrow(filtered_baseline_df)),
                             size = B,
                             replace = TRUE)
    expression_df <- filtered_baseline_df[sampled_indices, ]
    n_genes <- B
    cat("Random mode with replacement: Sampled", B, "genes from", post_filter_n, "available genes\n")
  }

  # Combine fold changes with expression parameters
  # Generate gRNA effect sizes for each target (K gRNAs per target)
  # Each target will have gRNAs_per_target individual gRNA effect sizes
  avg_fold_change <- numeric(n_genes)
  avg_fold_change_sq <- numeric(n_genes)
  
  for (i in 1:n_genes) {
    # Generate gRNAs_per_target effect sizes for this target
    grna_effects <- stats::rnorm(n = gRNAs_per_target, 
                                mean = fold_change_mean, 
                                sd = fold_change_sd)
    
    # Calculate moments: mean and mean of squares
    avg_fold_change[i] <- mean(grna_effects)
    avg_fold_change_sq[i] <- mean(grna_effects^2)
  }
  
  # Create DataFrame with new random effect size format
  fc_expression_df <- data.frame(
    avg_fold_change = avg_fold_change,
    avg_fold_change_sq = avg_fold_change_sq
  ) |>
    dplyr::bind_cols(expression_df)

  ################## extract the expression-dispersion curve ###################
  expression_dispersion_curve <- baseline_expression_stats$expression_dispersion_curve

  # return the data frame
  return(list(
    fc_expression_df = fc_expression_df,
    expression_dispersion_curve = expression_dispersion_curve,
    fold_change_mean = fold_change_mean  # Include for adaptive grid generation
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
  baseline_expression_list <- switch(biological_system,
         K562 = {

           # load the Gasperini baseline expression list
           rds_path <- system.file("extdata/baseline_expression", "Gasperini_expression.rds", package = "perturbplan", mustWork = TRUE)
           readRDS(rds_path)

         },
         # Default case for unsupported biological systems
         {
           stop("Unsupported biological system: '", biological_system, "'. ", 
                "Supported systems: 'K562'. For other systems, use custom baseline data.")
         })

  # return the data frame with the susbampled rows
  return(baseline_expression_list)
}

#' Validate custom baseline expression data
#'
#' @description
#' This function validates that custom baseline expression data has the required
#' structure and data types for use in power analysis.
#'
#' @param data Data frame containing custom baseline expression data
#' @param file_path Character. Optional file path for error messages (default: "uploaded file")
#'
#' @return List with validation results:
#' \describe{
#'   \item{valid}{Logical. TRUE if data passes all validation checks}
#'   \item{data}{Data frame. Cleaned and validated data (if valid=TRUE)}
#'   \item{errors}{Character vector. Error messages (if valid=FALSE)}
#'   \item{warnings}{Character vector. Warning messages}
#'   \item{summary}{Character. Summary statistics for display}
#' }
#'
#' @details
#' Required columns:
#' - response_id: Character vector of gene IDs (preferably Ensembl IDs)
#' - relative_expression: Numeric vector of expression values (TPM/1e6 scale)
#' - expression_size: Numeric vector of dispersion parameters (positive values)
#'
#' @export
validate_custom_baseline <- function(data, file_path = "uploaded file") {
  
  errors <- character(0)
  warnings <- character(0)
  
  # Check if data is a data frame
  if (!is.data.frame(data)) {
    errors <- c(errors, "Data must be a data frame")
    return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
  }
  
  # Check required columns exist
  required_cols <- c("response_id", "relative_expression", "expression_size")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # If missing required columns, return early
  if (length(errors) > 0) {
    return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
  }
  
  # Check data types and values
  if (!is.character(data$response_id) && !is.factor(data$response_id)) {
    errors <- c(errors, "response_id column must be character or factor")
  }
  
  if (!is.numeric(data$relative_expression)) {
    errors <- c(errors, "relative_expression column must be numeric")
  } else {
    # Check for negative values
    if (any(data$relative_expression < 0, na.rm = TRUE)) {
      errors <- c(errors, "relative_expression values must be non-negative")
    }
    # Check for extremely large values (potential TPM instead of relative scale)
    if (any(data$relative_expression > 1, na.rm = TRUE)) {
      max_val <- max(data$relative_expression, na.rm = TRUE)
      warnings <- c(warnings, paste0("Some relative_expression values > 1 (max: ", round(max_val, 4), 
                                    "). Ensure values are on TPM/1e6 scale, not raw TPM."))
    }
  }
  
  if (!is.numeric(data$expression_size)) {
    errors <- c(errors, "expression_size column must be numeric")
  } else {
    # Check for non-positive values
    if (any(data$expression_size <= 0, na.rm = TRUE)) {
      errors <- c(errors, "expression_size values must be positive")
    }
  }
  
  # Check for missing values
  for (col in required_cols) {
    if (any(is.na(data[[col]]))) {
      na_count <- sum(is.na(data[[col]]))
      errors <- c(errors, paste0("Column '", col, "' contains ", na_count, " missing values"))
    }
  }
  
  # Check for duplicate gene IDs
  if (any(duplicated(data$response_id))) {
    dup_count <- sum(duplicated(data$response_id))
    warnings <- c(warnings, paste0(dup_count, " duplicate gene IDs found. Only first occurrence will be used."))
    data <- data[!duplicated(data$response_id), ]
  }
  
  # Check for reasonable number of genes
  n_genes <- nrow(data)
  if (n_genes < 100) {
    warnings <- c(warnings, paste0("Only ", n_genes, " genes provided. Consider using more genes for robust analysis."))
  }
  
  # Validate Ensembl gene ID format (warn if not)
  ensembl_pattern <- "^ENSG[0-9]{11}$"
  non_ensembl <- !grepl(ensembl_pattern, data$response_id)
  if (any(non_ensembl)) {
    non_ensembl_count <- sum(non_ensembl)
    warnings <- c(warnings, paste0(non_ensembl_count, " gene IDs do not match Ensembl format (ENSGXXXXXXXXXXX)"))
  }
  
  # Create summary statistics
  if (length(errors) == 0) {
    # Convert relative expression back to TPM for display
    tpm_values <- data$relative_expression * 1e6
    summary_text <- paste0(
      "Loaded custom baseline expression (", formatC(n_genes, format = "d", big.mark = ","), " genes)<br/>",
      "Average TPM: ", round(mean(tpm_values, na.rm = TRUE), 1), 
      ", Range: ", round(min(tpm_values, na.rm = TRUE), 1), " - ", round(max(tpm_values, na.rm = TRUE), 1), "<br/>",
      "Expression size range: ", round(min(data$expression_size, na.rm = TRUE), 2), 
      " - ", round(max(data$expression_size, na.rm = TRUE), 2)
    )
  } else {
    summary_text <- ""
  }
  
  # Return validation results
  valid <- length(errors) == 0
  return(list(
    valid = valid,
    data = if (valid) data else NULL,
    errors = errors,
    warnings = warnings,
    summary = summary_text
  ))
}

#' Validate custom baseline expression RDS data
#'
#' @description
#' This function validates that custom baseline expression RDS data has the required
#' structure matching the output of extract_baseline_expression().
#'
#' @param data List object from RDS file containing baseline expression data
#' @param file_path Character. Optional file path for error messages (default: "uploaded file")
#'
#' @return List with validation results:
#' \describe{
#'   \item{valid}{Logical. TRUE if data passes all validation checks}
#'   \item{data}{List. Validated data structure (if valid=TRUE)}
#'   \item{errors}{Character vector. Error messages (if valid=FALSE)}
#'   \item{warnings}{Character vector. Warning messages}
#'   \item{summary}{Character. Summary statistics for display}
#' }
#'
#' @details
#' Required structure:
#' - List with two elements: 'baseline_expression' and 'expression_dispersion_curve'
#' - baseline_expression: Data frame with columns 'response_id', 'relative_expression', 'expression_size'
#' - expression_dispersion_curve: Function that takes a numeric vector and returns dispersion values
#'
#' @export
validate_custom_baseline_rds <- function(data, file_path = "uploaded file") {
  
  errors <- character(0)
  warnings <- character(0)
  
  # Check if data is a list
  if (!is.list(data)) {
    if (is.data.frame(data)) {
      errors <- c(errors, "Uploaded file contains a data frame, but should be a list. Please use the data frame as 'baseline_expression' element within a list structure.")
    } else {
      errors <- c(errors, paste("RDS data must be a list, but received:", class(data)[1], ". Please check the file structure."))
    }
    return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
  }
  
  # Check required list elements exist
  required_elements <- c("baseline_expression", "expression_dispersion_curve")
  missing_elements <- setdiff(required_elements, names(data))
  if (length(missing_elements) > 0) {
    errors <- c(errors, paste("Missing required list elements:", paste(missing_elements, collapse = ", ")))
  }
  
  # If missing required elements, return early
  if (length(errors) > 0) {
    return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
  }
  
  # Validate baseline_expression data frame
  baseline_df <- data$baseline_expression
  if (!is.data.frame(baseline_df)) {
    errors <- c(errors, "baseline_expression must be a data frame")
  } else {
    # Check required columns exist
    required_cols <- c("response_id", "relative_expression", "expression_size")
    missing_cols <- setdiff(required_cols, colnames(baseline_df))
    if (length(missing_cols) > 0) {
      errors <- c(errors, paste("baseline_expression missing columns:", paste(missing_cols, collapse = ", ")))
    } else {
      # Validate column data types and values
      if (!is.character(baseline_df$response_id) && !is.factor(baseline_df$response_id)) {
        errors <- c(errors, "response_id column must be character or factor")
      }
      
      if (!is.numeric(baseline_df$relative_expression)) {
        errors <- c(errors, "relative_expression column must be numeric")
      } else {
        if (any(baseline_df$relative_expression < 0, na.rm = TRUE)) {
          errors <- c(errors, "relative_expression values must be non-negative")
        }
        if (any(baseline_df$relative_expression > 1, na.rm = TRUE)) {
          max_val <- max(baseline_df$relative_expression, na.rm = TRUE)
          warnings <- c(warnings, paste0("Some relative_expression values > 1 (max: ", round(max_val, 4), 
                                        "). Ensure values are on TPM/1e6 scale, not raw TPM."))
        }
      }
      
      if (!is.numeric(baseline_df$expression_size)) {
        errors <- c(errors, "expression_size column must be numeric")
      } else {
        if (any(baseline_df$expression_size <= 0, na.rm = TRUE)) {
          errors <- c(errors, "expression_size values must be positive")
        }
      }
      
      # Check for missing values
      for (col in required_cols) {
        if (any(is.na(baseline_df[[col]]))) {
          na_count <- sum(is.na(baseline_df[[col]]))
          errors <- c(errors, paste0("baseline_expression column '", col, "' contains ", na_count, " missing values"))
        }
      }
      
      # Check for duplicate gene IDs
      if (any(duplicated(baseline_df$response_id))) {
        dup_count <- sum(duplicated(baseline_df$response_id))
        warnings <- c(warnings, paste0(dup_count, " duplicate gene IDs found. Only first occurrence will be used."))
        baseline_df <- baseline_df[!duplicated(baseline_df$response_id), ]
        data$baseline_expression <- baseline_df
      }
    }
  }
  
  # Validate expression_dispersion_curve function
  if (!is.function(data$expression_dispersion_curve)) {
    errors <- c(errors, "expression_dispersion_curve must be a function")
  } else {
    # Test if the function works with a numeric input
    tryCatch({
      test_result <- data$expression_dispersion_curve(c(1, 10, 100))
      if (!is.numeric(test_result) || length(test_result) != 3) {
        errors <- c(errors, "expression_dispersion_curve must return numeric values of same length as input")
      }
    }, error = function(e) {
      errors <- c(errors, paste("expression_dispersion_curve function error:", e$message))
    })
  }
  
  # Create summary statistics
  if (length(errors) == 0 && is.data.frame(baseline_df)) {
    n_genes <- nrow(baseline_df)
    
    # Check for reasonable number of genes
    if (n_genes < 100) {
      warnings <- c(warnings, paste0("Only ", n_genes, " genes provided. Consider using more genes for robust analysis."))
    }
    
    # Validate Ensembl gene ID format (warn if not)
    ensembl_pattern <- "^ENSG[0-9]{11}$"
    non_ensembl <- !grepl(ensembl_pattern, baseline_df$response_id)
    if (any(non_ensembl)) {
      non_ensembl_count <- sum(non_ensembl)
      warnings <- c(warnings, paste0(non_ensembl_count, " gene IDs do not match Ensembl format (ENSGXXXXXXXXXXX)"))
    }
    
    # Convert relative expression back to TPM for display
    tpm_values <- baseline_df$relative_expression * 1e6
    summary_text <- paste0(
      "Loaded custom baseline expression (", formatC(n_genes, format = "d", big.mark = ","), " genes)<br/>",
      "Average TPM: ", round(mean(tpm_values, na.rm = TRUE), 1)
    )
  } else {
    summary_text <- ""
  }
  
  # Return validation results
  valid <- length(errors) == 0
  return(list(
    valid = valid,
    data = if (valid) data else NULL,
    errors = errors,
    warnings = warnings,
    summary = summary_text
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
#' @param custom_library_data List. Optional custom library data with UMI_per_cell and variation parameters.
#' If provided, this data will be used instead of the biological_system data.
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
#' When custom_library_data is provided, it should be a list with:
#' \itemize{
#'   \item UMI_per_cell: Maximum UMI per cell parameter (positive numeric)
#'   \item variation: Variation parameter for PCR bias (positive numeric)
#' }
#'
#' @seealso
#' \code{\link{fit_read_UMI_curve}} for using these parameters
#' \code{\link{library_computation}} for fitting these parameters from data
#' \code{\link{validate_custom_library_rds}} for custom data validation
#' @export
extract_library_info <- function(biological_system = "K562", custom_library_data = NULL){

  # Use custom library data if provided, otherwise use default biological system data
  if (!is.null(custom_library_data)) {
    return(custom_library_data)
  }

  # sample baseline expression based on biological system
  switch(biological_system,
         K562 = {

           # load the Gasperini baseline expression list
           rds_path <- system.file("extdata/library_info", "Gasperini_library.rds", package = "perturbplan", mustWork = TRUE)
           library_info <- readRDS(rds_path)

         })

  # Extract scalar values from named vector
  params <- library_info$S_M_curve_params

  # return the data frame with the susbampled rows
  return(list(
    UMI_per_cell = unname(as.numeric(params[["UMI_per_cell"]])),
    variation = unname(as.numeric(params[["variation"]]))
  ))
}

#' Validate custom library RDS file structure and content
#'
#' @description
#' This function validates that an uploaded RDS file contains valid library parameters 
#' with the correct structure and value ranges for power analysis.
#'
#' @param data The loaded RDS data to validate
#' @param filename Character. The original filename for error reporting (optional)
#'
#' @return A list with validation results:
#' \describe{
#'   \item{valid}{Logical. TRUE if validation passed}
#'   \item{data}{The validated data if valid, NULL otherwise}
#'   \item{errors}{Character vector of error messages}
#'   \item{warnings}{Character vector of warning messages}
#'   \item{summary}{Character. HTML summary for display}
#' }
#'
#' @details
#' Expected RDS structure:
#' \code{
#' list(
#'   UMI_per_cell = numeric_value,  # Positive number
#'   variation = numeric_value      # Positive number between 0 and 1
#' )
#' }
#'
#' @export
validate_custom_library_rds <- function(data, filename = "uploaded file") {
  errors <- character(0)
  warnings <- character(0)
  
  # Check if data is a list
  if (!is.list(data)) {
    errors <- c(errors, "RDS file must contain a list object")
    return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
  }
  
  # Check required elements
  required_elements <- c("UMI_per_cell", "variation")
  missing_elements <- setdiff(required_elements, names(data))
  if (length(missing_elements) > 0) {
    errors <- c(errors, paste("Missing required elements:", paste(missing_elements, collapse = ", ")))
  }
  
  # Check for unexpected elements
  unexpected_elements <- setdiff(names(data), required_elements)
  if (length(unexpected_elements) > 0) {
    warnings <- c(warnings, paste("Unexpected elements will be ignored:", paste(unexpected_elements, collapse = ", ")))
  }
  
  # Validate UMI_per_cell
  if ("UMI_per_cell" %in% names(data)) {
    umi_val <- data$UMI_per_cell
    if (!is.numeric(umi_val) || length(umi_val) != 1) {
      errors <- c(errors, "UMI_per_cell must be a single numeric value")
    } else if (is.na(umi_val) || !is.finite(umi_val)) {
      errors <- c(errors, "UMI_per_cell cannot be NA or infinite")
    } else if (umi_val <= 0) {
      errors <- c(errors, "UMI_per_cell must be positive")
    } else if (umi_val < 1000) {
      warnings <- c(warnings, "UMI_per_cell is unusually low (< 1000)")
    } else if (umi_val > 50000) {
      warnings <- c(warnings, "UMI_per_cell is unusually high (> 50000)")
    }
  }
  
  # Validate variation
  if ("variation" %in% names(data)) {
    var_val <- data$variation
    if (!is.numeric(var_val) || length(var_val) != 1) {
      errors <- c(errors, "variation must be a single numeric value")
    } else if (is.na(var_val) || !is.finite(var_val)) {
      errors <- c(errors, "variation cannot be NA or infinite")
    } else if (var_val <= 0) {
      errors <- c(errors, "variation must be positive")
    } else if (var_val > 1) {
      warnings <- c(warnings, "variation parameter is unusually high (> 1)")
    }
  }
  
  # Generate summary if validation passed
  if (length(errors) == 0) {
    summary_text <- paste0(
      "Loaded custom library parameters<br/>",
      "UMI per cell: ", formatC(data$UMI_per_cell, format = "d", big.mark = ","), "<br/>",
      "Variation: ", formatC(data$variation, format = "e")
    )
  } else {
    summary_text <- ""
  }
  
  # Return validation results
  valid <- length(errors) == 0
  return(list(
    valid = valid,
    data = if (valid) data else NULL,
    errors = errors,
    warnings = warnings,
    summary = summary_text
  ))
}

#' Validate combined pilot data RDS file structure
#'
#' @description
#' This function validates the structure of a combined pilot data RDS file that contains
#' both baseline expression data and library parameters in a single list.
#'
#' @param data A list object loaded from an RDS file, expected to contain:
#'   \itemize{
#'     \item baseline_expression: A list with baseline expression data and dispersion curve
#'     \item library_parameters: A list with UMI_per_cell and variation parameters
#'   }
#' @param file_path Character. Path or description of the uploaded file for error messages.
#'
#' @return A list with elements:
#' \describe{
#'   \item{valid}{Logical indicating if validation passed}
#'   \item{data}{The validated data (if valid) or NULL}
#'   \item{errors}{Character vector of error messages}
#'   \item{warnings}{Character vector of warning messages}
#'   \item{summary}{HTML-formatted summary text for display}
#' }
#'
#' @details
#' This function validates the overall structure of the combined pilot data, then
#' delegates validation of individual components to \code{\link{validate_custom_baseline_rds}}
#' and \code{\link{validate_custom_library_rds}}.
#'
#' Expected structure:
#' \preformatted{
#' list(
#'   baseline_expression = list(
#'     baseline_expression = data.frame(...),
#'     expression_dispersion_curve = function(v) {...}
#'   ),
#'   library_parameters = list(
#'     UMI_per_cell = numeric_value,
#'     variation = numeric_value
#'   )
#' )
#' }
#'
#' @seealso
#' \code{\link{validate_custom_baseline_rds}} for baseline expression validation
#' \code{\link{validate_custom_library_rds}} for library parameter validation
#' @export
validate_combined_pilot_data <- function(data, file_path = "uploaded file") {
  
  errors <- character(0)
  warnings <- character(0)
  summary_parts <- character(0)
  
  # Check if data is a list
  if (!is.list(data)) {
    errors <- c(errors, paste("RDS data must be a list, but received:", class(data)[1]))
    return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
  }
  
  # Check required top-level elements
  required_elements <- c("baseline_expression", "library_parameters")
  missing_elements <- setdiff(required_elements, names(data))
  if (length(missing_elements) > 0) {
    errors <- c(errors, paste("Missing required list elements:", paste(missing_elements, collapse = ", ")))
    errors <- c(errors, "Expected structure: list(baseline_expression = ..., library_parameters = ...)")
    return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
  }
  
  # Check for unexpected top-level elements
  unexpected_elements <- setdiff(names(data), required_elements)
  if (length(unexpected_elements) > 0) {
    warnings <- c(warnings, paste("Unexpected top-level elements will be ignored:", paste(unexpected_elements, collapse = ", ")))
  }
  
  # Validate baseline_expression component
  baseline_result <- validate_custom_baseline_rds(data$baseline_expression, 
                                                  paste0(file_path, " > baseline_expression"))
  
  # Validate library_parameters component  
  library_result <- validate_custom_library_rds(data$library_parameters,
                                                paste0(file_path, " > library_parameters"))
  
  # Combine validation results
  all_errors <- c(errors, baseline_result$errors, library_result$errors)
  all_warnings <- c(warnings, baseline_result$warnings, library_result$warnings)
  
  # Create combined summary
  if (baseline_result$valid && library_result$valid) {
    summary_parts <- c(baseline_result$summary, library_result$summary)
    combined_summary <- paste(summary_parts, collapse = "<br/>")
  } else {
    combined_summary <- ""
  }
  
  # Overall validation status
  valid <- length(all_errors) == 0
  
  # Prepare validated data
  validated_data <- NULL
  if (valid) {
    validated_data <- list(
      baseline_expression = baseline_result$data,
      library_parameters = library_result$data
    )
  }
  
  return(list(
    valid = valid,
    data = validated_data,
    errors = all_errors,
    warnings = all_warnings,
    summary = combined_summary
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
#'
#' @seealso \code{\link{extract_library_info}} for obtaining curve parameters
#' @export
fit_read_UMI_curve <- function(reads_per_cell, UMI_per_cell, variation){
  
  # Wrapper function that calls the optimized C++ implementation
  return(fit_read_UMI_curve_cpp(reads_per_cell, UMI_per_cell, variation))
}

#' Identify optimal reads per cell range for power analysis grid
#'
#' @description
#' This function determines the minimum and maximum reads per cell values for 
#' power analysis grid generation based on experimental platform capabilities
#' and S-M curve saturation analysis.
#'
#' @param experimental_platform Character. Experimental platform identifier
#'   (e.g., "10x Chromium v3", "Other").
#' @param library_info List. Output from extract_library_info() containing
#'   UMI_per_cell and variation parameters for S-M curve analysis.
#'
#' @return List with elements:
#' \describe{
#'   \item{min_reads_per_cell}{Minimum reads per cell based on platform}
#'   \item{max_reads_per_cell}{Maximum reads per cell for ~80% UMI saturation}
#' }
#'
#' @details
#' This function is a wrapper around the optimized C++ implementation 
#' \code{identify_library_size_range_cpp} which provides significant 
#' performance improvements for power analysis computations.
#' 
#' The function operates in two phases:
#' 
#' **Minimum determination**: Platform-specific minimum sequencing depth
#' based on typical experimental capabilities and quality thresholds.
#' 
#' **Maximum determination**: Uses binary search on the S-M curve to find
#' the reads per cell that achieves approximately 80% UMI saturation. If
#' 80% saturation is not achievable within practical limits (10x UMI_per_cell),
#' returns the practical upper bound.
#'
#' @seealso 
#' \code{\link{fit_read_UMI_curve}} for S-M curve evaluation
#' \code{\link{extract_library_info}} for obtaining library parameters
#' \code{\link{identify_library_size_range_cpp}} for C++ implementation
#' @export
identify_library_size_range <- function(experimental_platform, library_info) {
  
  # Input validation for library_info structure
  if (!is.list(library_info) || !all(c("UMI_per_cell", "variation") %in% names(library_info))) {
    stop("library_info must be a list with UMI_per_cell and variation elements")
  }
  
  # Extract parameters and call optimized C++ implementation
  UMI_per_cell <- library_info$UMI_per_cell
  variation <- library_info$variation
  
  # Wrapper around the C++ implementation
  return(identify_library_size_range_cpp(experimental_platform, UMI_per_cell, variation))
}

#' Identify optimal cell and read range for power analysis with logarithmic grid generation
#'
#' @description
#' This function combines library size range determination with cell range optimization
#' to generate logarithmically-spaced experimental design grids for power analysis.
#' It integrates identify_library_size_range_cpp() and identify_cell_range_cpp() with
#' treatment/control cell allocation calculations.
#'
#' @param experimental_platform String. Experimental platform identifier for read range determination.
#' @param fc_expression_info List from extract_fc_expression_info() containing fc_expression_df and expression_dispersion_curve.
#' @param library_info List from extract_library_info() containing UMI_per_cell and variation parameters.
#' @param grid_size Integer. Number of points in each dimension of the grid (default: 10).
#' @param min_power_threshold Numeric. Minimum power threshold for cell range determination (default: 0.01).
#' @param max_power_threshold Numeric. Maximum power threshold for cell range determination (default: 0.8).
#' @param MOI Numeric. Multiplicity of infection for cell allocation calculations (default: 10).
#' @param num_targets Integer. Number of targets for cell allocation calculations (default: 100).
#' @param gRNAs_per_target Integer. Number of gRNAs per target (default: 4).
#' @param non_targeting_gRNAs Integer. Number of non-targeting gRNAs (default: 10).
#' @param control_group String. Control group type: "complement" or "nt_cells" (default: "complement").
#' @param multiple_testing_alpha Numeric. Alpha level for multiple testing (default: 0.05).
#' @param side String. Test sidedness: "left", "right", or "both" (default: "left").
#' @param prop_non_null Numeric. Proportion of non-null hypotheses (default: 0.1).
#'
#' @return A list containing:
#' \describe{
#'   \item{cells_seq}{Numeric vector of logarithmically-spaced total cell counts}
#'   \item{reads_seq}{Numeric vector of logarithmically-spaced reads per cell values}
#'   \item{library_size_seq}{Numeric vector of library sizes corresponding to reads_seq}
#'   \item{num_trt_cells_seq}{Numeric vector of treatment cell counts for each cell count}
#'   \item{num_cntrl_cells_seq}{Numeric vector of control cell counts for each cell count}
#'   \item{cell_range}{List with detailed cell range determination results}
#'   \item{reads_range}{List with reads per cell range determination results}
#'   \item{grid_size}{Integer grid size used}
#' }
#'
#' @details
#' The function operates in three main steps:
#' \enumerate{
#'   \item \strong{Reads Range}: Uses identify_library_size_range_cpp() to determine optimal reads per cell range based on UMI saturation curves
#'   \item \strong{Cell Range}: Uses identify_cell_range_cpp() to determine optimal cell count range based on power thresholds
#'   \item \strong{Grid Generation}: Creates logarithmically-spaced sequences and pre-computes treatment/control cell allocations
#' }
#'
#' Both cell and read sequences use logarithmic spacing to provide good coverage across multiple orders of magnitude.
#' Treatment and control cell counts are calculated based on MOI, gRNA ratios, and control group specification.
#'
#' @examples
#' \dontrun{
#' # Extract required info
#' fc_info <- extract_fc_expression_info(0.8, 0.1, "K562", B = 100)
#' lib_info <- extract_library_info("K562")
#' 
#' # Generate experimental design
#' design <- identify_cell_read_range(
#'   experimental_platform = "10x Chromium v3",
#'   fc_expression_info = fc_info,
#'   library_info = lib_info,
#'   grid_size = 10
#' )
#' 
#' # Examine the design
#' print(design$cells_seq)      # Total cells: 100, 279, 777, ...
#' print(design$reads_seq)      # Reads/cell: 500, 1753, 6145, ...
#' print(design$num_trt_cells_seq)   # Treatment cells for each total
#' }
#'
#' @seealso 
#' \code{\link{identify_library_size_range_cpp}} for reads per cell range determination
#' \code{\link{identify_cell_range_cpp}} for cell count range determination
#' @export
identify_cell_read_range <- function(
  experimental_platform,
  fc_expression_info, 
  library_info,
  grid_size = 10,
  min_power_threshold = 0.01,
  max_power_threshold = 0.8,
  MOI = 10,
  num_targets = 100,
  gRNAs_per_target = 4,
  non_targeting_gRNAs = 10,
  control_group = "complement",
  multiple_testing_alpha = 0.05,
  side = "left",
  prop_non_null = 0.1
) {
  
  # Input validation
  if (grid_size < 2) {
    stop("grid_size must be at least 2")
  }
  if (min_power_threshold <= 0 || min_power_threshold >= 1) {
    stop("min_power_threshold must be between 0 and 1")
  }
  if (max_power_threshold <= 0 || max_power_threshold >= 1) {
    stop("max_power_threshold must be between 0 and 1")
  }
  if (min_power_threshold >= max_power_threshold) {
    stop("min_power_threshold must be < max_power_threshold")
  }
  
  # Step 1: Determine optimal reads per cell range using S-M curve analysis
  reads_range <- identify_library_size_range(
    experimental_platform = experimental_platform,
    library_info = library_info
  )
  
  # Step 2: Determine optimal cell count range based on power thresholds
  cell_range <- identify_cell_range_cpp(
    min_reads_per_cell = reads_range$min_reads_per_cell,
    max_reads_per_cell = reads_range$max_reads_per_cell,
    fc_expression_df = fc_expression_info$fc_expression_df,
    UMI_per_cell = library_info$UMI_per_cell,
    variation = library_info$variation,
    MOI = MOI,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    control_group = control_group,
    multiple_testing_alpha = multiple_testing_alpha,
    side = side,
    prop_non_null = prop_non_null,
    min_power_threshold = min_power_threshold,
    max_power_threshold = max_power_threshold,
    cell_lower_bound = 100.0,   # Start from 100 treatment cells for practical experimental design
    cell_upper_bound = 1e5      # Cap at 1e5 treatment cells for more reasonable ranges
  )
  
  # Step 3: Generate logarithmically-spaced sequences
  
  # Cell sequence (logarithmic spacing)
  log_min_cells <- log10(cell_range$min_cells)
  log_max_cells <- log10(cell_range$max_cells)
  log_cells_seq <- seq(log_min_cells, log_max_cells, length.out = grid_size)
  cells_seq <- round(10^log_cells_seq)
  
  # Reads sequence (logarithmic spacing)
  log_min_reads <- log10(reads_range$min_reads_per_cell)
  log_max_reads <- log10(reads_range$max_reads_per_cell)
  log_reads_seq <- seq(log_min_reads, log_max_reads, length.out = grid_size)
  reads_seq <- round(10^log_reads_seq)
  
  # Step 4: Pre-compute treatment and control cell allocations
  total_gRNAs <- num_targets * gRNAs_per_target + non_targeting_gRNAs
  num_trt_cells_seq <- (gRNAs_per_target * cells_seq * MOI) / total_gRNAs
  
  if (control_group == "complement") {
    num_cntrl_cells_seq <- cells_seq - num_trt_cells_seq
  } else {  # "nt_cells"
    num_cntrl_cells_seq <- (non_targeting_gRNAs * cells_seq * MOI) / total_gRNAs
  }
  
  # Round to whole cells
  num_trt_cells_seq <- round(num_trt_cells_seq)
  num_cntrl_cells_seq <- round(num_cntrl_cells_seq)
  
  # Calculate library sizes corresponding to reads_seq
  library_size_seq <- fit_read_UMI_curve(
    reads_per_cell = reads_seq,
    UMI_per_cell = library_info$UMI_per_cell,
    variation = library_info$variation
  )
  
  # Return comprehensive experimental design
  return(list(
    cells_seq = cells_seq,
    reads_seq = reads_seq,
    library_size_seq = library_size_seq,
    num_trt_cells_seq = num_trt_cells_seq,
    num_cntrl_cells_seq = num_cntrl_cells_seq,
    cell_range = cell_range,
    reads_range = reads_range,
    grid_size = grid_size
  ))
}

#' Convert identify_cell_read_range output to data frame for calculate_power_grid
#'
#' @description
#' This helper function converts the output from identify_cell_read_range() into
#' the data frame format expected by calculate_power_grid(). It creates the full
#' experimental design grid by expanding all combinations of cells and reads.
#'
#' @param design_output List output from identify_cell_read_range() containing
#'   cells_seq, reads_seq, library_size_seq, num_trt_cells_seq, and num_cntrl_cells_seq.
#'
#' @return Data frame with columns: num_total_cells, reads_per_cell, library_size,
#'   num_trt_cells, num_cntrl_cells. Each row represents one experimental 
#'   condition (cell count × reads per cell combination).
#'
#' @details
#' The function creates all combinations of cell counts and reads per cell values,
#' then assigns the appropriate treatment and control cell counts and library sizes 
#' based on the total cell count and reads per cell for each row.
#'
#' @examples
#' \dontrun{
#' # Generate experimental design
#' fc_info <- extract_fc_expression_info(0.8, 0.1, "K562", B = 100)
#' lib_info <- extract_library_info("K562")
#' design <- identify_cell_read_range("10x Chromium v3", fc_info, lib_info)
#' 
#' # Convert to data frame for power analysis
#' cells_reads_df <- convert_design_to_dataframe(design)
#' print(head(cells_reads_df))
#' }
#'
#' @seealso \code{\link{identify_cell_read_range}} for experimental design generation
#' @export
convert_design_to_dataframe <- function(design_output) {
  
  # Input validation
  required_elements <- c("cells_seq", "reads_seq", "num_trt_cells_seq", "num_cntrl_cells_seq")
  missing_elements <- setdiff(required_elements, names(design_output))
  if (length(missing_elements) > 0) {
    stop("Missing required elements in design_output: ", paste(missing_elements, collapse = ", "))
  }
  
  # Create all combinations of cells and reads
  cells_reads_df <- expand.grid(
    cells_idx = 1:length(design_output$cells_seq),
    reads_idx = 1:length(design_output$reads_seq)
  ) |>
    dplyr::mutate(
      reads_per_cell = design_output$reads_seq[reads_idx],
      library_size = if("library_size_seq" %in% names(design_output)) design_output$library_size_seq[reads_idx] else NA,
      num_trt_cells = design_output$num_trt_cells_seq[cells_idx],
      num_cntrl_cells = design_output$num_cntrl_cells_seq[cells_idx]
    ) |>
    dplyr::select(reads_per_cell, library_size, num_trt_cells, num_cntrl_cells)
  
  return(cells_reads_df)
}
