# Variable bindings are handled in R/perturbplan.R

#' Get pilot data from package data directory
#'
#' @description
#' Internal function to load baseline expression and library parameters from the
#' pilot datasets stored in the package data/ directory.
#'
#' @param biological_system Character. The biological system name (e.g., "K562", "A549", "THP-1", "T_CD8", "iPSC")
#'
#' @return A list containing:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene expression data or list with baseline_expression_stats data frame}
#'   \item{library_parameters}{List with UMI_per_cell and variation parameters}
#' }
#'
#' @examples
#' # Load pilot data for K562 cells
#' k562_data <- get_pilot_data_from_package("K562")
#'
#' # View the structure
#' str(k562_data)
#'
#' # Access baseline expression data
#' head(k562_data$baseline_expression_stats)
#'
#' # Access library parameters including mapping efficiency
#' k562_data$library_parameters
#' cat("Mapping efficiency:", k562_data$mapping_efficiency)
#'
#' # The mapping efficiency affects power calculations by determining
#' # what fraction of sequencing reads contribute to gene expression
#' # Higher mapping efficiency means more effective sequencing depth
#'
#' # Compare mapping efficiency across cell types
#' a549_data <- get_pilot_data_from_package("A549")
#' cat("K562 mapping efficiency:", k562_data$mapping_efficiency)
#' cat("A549 mapping efficiency:", a549_data$mapping_efficiency)
#' @importFrom utils data
#' @keywords internal
#' @export
get_pilot_data_from_package <- function(biological_system) {
  # Map biological system names to data file names
  data_mapping <- list(
    "K562" = "K562_10x",
    "A549" = "A549_10x",
    "THP-1" = "THP_1_10x",
    "T_CD8" = "T_CD8_10x",
    "iPSC" = "iPSC_10x",
    "iPSC_neuron" = "iPSC_neuron_10x"
  )

  # Check if biological system is supported
  if (!biological_system %in% names(data_mapping)) {
    # Load reference datasets to show available options
    env <- environment()
    data("reference_expression_datasets", package = "perturbplan", envir = env)
    ref_datasets <- get("reference_expression_datasets", envir = env)
    available_systems <- paste(ref_datasets$cell_type, collapse = ", ")
    stop("Unsupported biological system: '", biological_system, "'. Available options: ", available_systems)
  }

  # Get the data object name
  data_name <- data_mapping[[biological_system]]

  # Load the data from the package data/ directory
  data(list = data_name, package = "perturbplan", envir = environment())

  # Return the loaded data object
  return(get(data_name, envir = environment()))
}

#' Extract fold change and expression information for power analysis
#'
#' @description
#' This function combines fold change effect size sampling with baseline expression
#' data to create a comprehensive dataset for Monte Carlo power analysis simulations.
#' It can handle both user-specified genes and random sampling scenarios.
#'
#' @param minimum_fold_change Numeric. Minimum expected fold change effect (mean of gRNA effect distribution).
#' @param gRNA_variability Numeric. Standard deviation of gRNA effect sizes, representing variability between gRNAs targeting the same gene.
#' @param biological_system Character. Biological system for baseline expression. Available options:
#'   "K562", "A549", "THP-1", "T_CD8", "iPSC" (default: "K562").
#' @param B Integer. Number of Monte Carlo samples to generate when gene_list is NULL (default: 200).
#'   Ignored when gene_list is provided.
#' @param gene_list Character vector. Optional list of Ensembl gene IDs to use for analysis.
#'   If provided, expression parameters will be extracted for ALL specified genes (no sampling).
#'   If NULL (default), B genes are randomly sampled from baseline data.
#' @param TPM_threshold Numeric. Minimum TPM threshold (default: 10). Genes with expression
#'   levels below TPM_threshold/1e6 are filtered out before power calculation.
#' @param custom_pilot_data List. Optional custom pilot data. If provided,
#'   this data is used instead of the default biological_system data. Must contain
#'   baseline_expression_stats (data frame with relative_expression and expression_size columns)
#'   and library_parameters (with UMI_per_cell and variation). See
#'   \code{\link{reference_data_preprocessing_10x}} for processing 10x Cell Ranger output and
#'   \code{\link{reference_data_processing}} for further pilot data processing.
#' @param gRNAs_per_target Integer. Number of gRNAs per target (default: 4).
#'   Each target will have gRNAs_per_target individual gRNA effect sizes drawn from the
#'   specified fold change distribution. avg_fold_change and avg_fold_change_sq are
#'   calculated as the mean and mean-of-squares of these gRNA effect sizes.
#'
#' @return A list with elements:
#' \describe{
#'   \item{fc_expression_df}{Data frame with avg_fold_change, avg_fold_change_sq, and expression parameters}
#'   \item{minimum_fold_change}{Numeric. The input minimum fold change for reference}
#'   \item{pilot_data}{List. Complete pilot data object for further use}
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
#'   \item Filters genes below TPM threshold (relative_expression < TPM_threshold/1e6)
#'   \item Generates gRNAs_per_target effect sizes per target from a normal distribution
#'   \item Calculates avg_fold_change and avg_fold_change_sq from the gRNA effect sizes
#'   \item Returns combined data for Monte Carlo integration with random effect sizes
#' }
#'
#' @examples
#' # Extract fold change and expression information
#' fc_expr_data <- extract_fc_expression_info(
#'   minimum_fold_change = 0.8,
#'   gRNA_variability = 0.13,
#'   biological_system = "K562",
#'   B = 200,
#'   TPM_threshold = 10,
#'   gRNAs_per_target = 4
#' )
#' @seealso \code{\link{get_pilot_data_from_package}} for direct pilot data access
#' @keywords internal
#' @export
extract_fc_expression_info <- function(minimum_fold_change, gRNA_variability, biological_system =  "K562", B = 200, gene_list = NULL, TPM_threshold = 10, custom_pilot_data = NULL, gRNAs_per_target = 4){

  # set the random seed
  set.seed(1)

  ############## combine expression and effect size information ################
  # Use custom pilot data if provided, otherwise load from data/ directory
  if (!is.null(custom_pilot_data)) {
    pilot_data <- custom_pilot_data
    # Try new key first, fall back to old key for backward compatibility
    baseline_expression_stats <- if (!is.null(custom_pilot_data$baseline_expression_stats)) {
      custom_pilot_data$baseline_expression_stats
    } else {
      custom_pilot_data$baseline_expression
    }
  } else {
    # Load complete pilot data from data/ directory based on biological_system
    pilot_data <- get_pilot_data_from_package(biological_system)
    # Try new key first, fall back to old key for backward compatibility
    baseline_expression_stats <- if (!is.null(pilot_data$baseline_expression_stats)) {
      pilot_data$baseline_expression_stats
    } else {
      pilot_data$baseline_expression
    }
  }

  # Handle both old nested structure and new simplified structure
  if (is.data.frame(baseline_expression_stats)) {
    # New simplified structure: baseline_expression_stats is directly a data frame
    baseline_df <- baseline_expression_stats
  } else {
    # Old nested structure: baseline_expression contains baseline_expression data frame
    baseline_df <- baseline_expression_stats$baseline_expression
  }

  #################### apply TPM threshold filtering FIRST ###################
  # Convert TPM threshold to relative expression scale (TPM / 1e6)
  TPM_threshold_relative <- TPM_threshold / 1e6

  # Filter the full baseline dataset by TPM threshold first
  if ("relative_expression" %in% colnames(baseline_df)) {
    pre_filter_n <- nrow(baseline_df)
    filtered_baseline_df <- baseline_df |>
      dplyr::filter(relative_expression >= TPM_threshold_relative)
    post_filter_n <- nrow(filtered_baseline_df)

    # Check if we have any genes left after filtering
    if (post_filter_n == 0) {
      stop("No genes remain after TPM threshold filtering. Consider lowering TPM_threshold.")
    }

    # Print filtering summary
    cat("TPM filtering: Kept", post_filter_n, "out of", pre_filter_n, "genes (threshold:", TPM_threshold, "TPM)\n")
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
                                 mean = minimum_fold_change,
                                 sd = gRNA_variability) |>
      pmax(.Machine$double.eps)  # Ensure no negative effects

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

  # return the data frame with complete pilot data access
  result <- list(
    fc_expression_df = fc_expression_df,
    minimum_fold_change = minimum_fold_change,  # Include for adaptive grid generation
    pilot_data = pilot_data  # Always include pilot data (either built-in or custom)
  )

  return(result)
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
#' @keywords internal
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
#' structure with baseline expression data and dispersion curve function.
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
#' @keywords internal
validate_custom_baseline_rds <- function(data, file_path = "uploaded file") {

  errors <- character(0)
  warnings <- character(0)

  # Handle both data frame (new format) and list (old format) inputs
  if (is.data.frame(data)) {
    # New simplified format: direct data frame
    baseline_df <- data
  } else if (is.list(data)) {
    # Check if it's new structure (baseline_expression_stats) or old structure (baseline_expression)
    if ("baseline_expression_stats" %in% names(data)) {
      # New structure - validate it has the baseline_expression_stats data frame
      if (!is.data.frame(data$baseline_expression_stats)) {
        errors <- c(errors, "baseline_expression_stats must be a data frame")
        return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
      }
      baseline_df <- data$baseline_expression_stats
    } else if ("baseline_expression" %in% names(data)) {
      # Old nested structure - validate it has the baseline_expression data frame
      if (!is.data.frame(data$baseline_expression)) {
        errors <- c(errors, "baseline_expression must be a data frame")
        return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
      }
      baseline_df <- data$baseline_expression

      # Warn about deprecated key name but don't fail
      warnings <- c(warnings, "Key 'baseline_expression' is deprecated, use 'baseline_expression_stats' instead")

      # Warn about deprecated structure but don't fail
      if ("expression_dispersion_curve" %in% names(data)) {
        warnings <- c(warnings, "expression_dispersion_curve is deprecated and will be ignored")
      }
    } else {
      errors <- c(errors, "List must contain 'baseline_expression_stats' (preferred) or 'baseline_expression' element")
      return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
    }
  } else {
    errors <- c(errors, paste("RDS data must be a data frame or list, but received:", class(data)[1]))
    return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
  }

  # Now validate the baseline_df we extracted
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

  # Prepare validated data in simplified format
  validated_data <- NULL
  if (valid) {
    if (is.data.frame(data)) {
      # Already in new format
      validated_data <- data
    } else {
      # Convert from old format to new format
      validated_data <- baseline_df
    }
  }

  return(list(
    valid = valid,
    data = validated_data,
    errors = errors,
    warnings = warnings,
    summary = summary_text
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
#' @keywords internal
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
#' Also, it checks whether the value of parameters makes sense.
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
#' @examples
#' # set seed for reproducibility
#' set.seed(123)
#' # First create pilot data using the preprocessing pipeline
#' extdata_path <- system.file("extdata", package = "perturbplan")
#' # Get raw data from 10x output
#' raw_data <- reference_data_preprocessing_10x(
#'   path_to_top_level_output = extdata_path,
#'   path_to_run_level_output = "cellranger_tiny",
#'   h5_rough = TRUE
#' )
#' # Process into final pilot data format
#' pilot_data <- reference_data_processing(
#'   response_matrix = raw_data$response_matrix,
#'   read_umi_table = raw_data$read_umi_table,
#'   mapping_efficiency = raw_data$mapping_efficiency,
#'   TPM_thres = 0.1,
#'   h5_only = FALSE
#' )
#'
#' # Validate the processed pilot data
#' validation_result <- validate_combined_pilot_data(pilot_data)
#' if (validation_result$valid) {
#'   cat("Validation passed:", validation_result$summary)
#' } else {
#'   cat("Validation failed:", validation_result$errors)
#' }
#'
#' pilot_data_corrected <- pilot_data
#' pilot_data_corrected$library_parameters$variation <- 0.3
#' validation_result2 <- validate_combined_pilot_data(pilot_data_corrected)
#' if (validation_result2$valid) {
#'   cat("Validation passed:", validation_result2$summary)
#' } else {
#'   cat("Validation failed:", validation_result2$errors)
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

  # Check for required elements - prefer new naming but support backward compatibility
  has_new_baseline <- "baseline_expression_stats" %in% names(data)
  has_old_baseline <- "baseline_expression" %in% names(data)
  has_library <- "library_parameters" %in% names(data)

  # Determine which baseline key to use
  baseline_key <- NULL
  if (has_new_baseline) {
    baseline_key <- "baseline_expression_stats"
  } else if (has_old_baseline) {
    baseline_key <- "baseline_expression"
    warnings <- c(warnings, "Using deprecated 'baseline_expression' key. Consider updating to 'baseline_expression_stats'")
  }

  # Check required elements
  if (is.null(baseline_key)) {
    errors <- c(errors, "Missing baseline expression data. Expected 'baseline_expression_stats' or 'baseline_expression'")
  }
  if (!has_library) {
    errors <- c(errors, "Missing required 'library_parameters' element")
  }

  if (length(errors) > 0) {
    errors <- c(errors, "Expected structure: list(baseline_expression_stats = data.frame(...), library_parameters = list(...))")
    return(list(valid = FALSE, data = NULL, errors = errors, warnings = warnings, summary = ""))
  }

  # Check for unexpected top-level elements
  expected_elements <- c("baseline_expression_stats", "baseline_expression", "library_parameters")
  unexpected_elements <- setdiff(names(data), expected_elements)
  if (length(unexpected_elements) > 0) {
    warnings <- c(warnings, paste("Unexpected top-level elements will be ignored:", paste(unexpected_elements, collapse = ", ")))
  }

  # Validate baseline expression component
  baseline_result <- validate_custom_baseline_rds(data[[baseline_key]],
                                                  paste0(file_path, " > ", baseline_key))

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

  # Prepare validated data - always use new naming in output
  validated_data <- NULL
  if (valid) {
    validated_data <- list(
      baseline_expression_stats = baseline_result$data,
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
#' @examples
#' # Get library parameters from pilot data
#' pilot_data <- get_pilot_data_from_package("K562")
#' library_params <- pilot_data$library_parameters
#'
#' # Define read depths to test
#' read_depths <- c(10000, 25000, 50000, 100000)
#'
#' # Calculate effective library sizes
#' effective_umis <- fit_read_UMI_curve(
#'   reads_per_cell = read_depths,
#'   UMI_per_cell = library_params$UMI_per_cell,
#'   variation = library_params$variation
#' )
#'
#' # View the results
#' data.frame(
#'   reads_per_cell = read_depths,
#'   effective_UMI = effective_umis,
#'   saturation_pct = round(100 * effective_umis / library_params$UMI_per_cell, 1)
#' )
#'
#' @examples
#' # Get library parameters from pilot data
#' pilot_data <- get_pilot_data_from_package("K562")
#' library_params <- pilot_data$library_parameters
#'
#' # Define read depths to test
#' read_depths <- c(10000, 25000, 50000, 100000)
#'
#' # Calculate effective library sizes
#' effective_umis <- fit_read_UMI_curve(
#'   reads_per_cell = read_depths,
#'   UMI_per_cell = library_params$UMI_per_cell,
#'   variation = library_params$variation
#' )
#'
#' # View the results
#' data.frame(
#'   reads_per_cell = read_depths,
#'   effective_UMI = effective_umis,
#'   saturation_pct = round(100 * effective_umis / library_params$UMI_per_cell, 1)
#' )
#'
#' @seealso \code{\link{get_pilot_data_from_package}} for obtaining curve parameters
#' @keywords internal
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
#' @param library_parameters List. Library parameters containing
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
#' \code{\link{get_pilot_data_from_package}} for obtaining library parameters
#' \code{\link{identify_library_size_range_cpp}} for C++ implementation
#' @keywords internal
identify_library_size_range <- function(experimental_platform, library_parameters) {

  # Input validation for library_parameters structure
  if (!is.list(library_parameters) || !all(c("UMI_per_cell", "variation") %in% names(library_parameters))) {
    stop("library_parameters must be a list with UMI_per_cell and variation elements")
  }

  # Extract parameters and call optimized C++ implementation
  UMI_per_cell <- library_parameters$UMI_per_cell
  variation <- library_parameters$variation

  # Wrapper around the C++ implementation
  return(identify_library_size_range_cpp(experimental_platform, UMI_per_cell, variation))
}

#' Extract baseline expression information without fold change augmentation
#'
#' @description
#' This function extracts and processes baseline expression data for power analysis
#' without adding fold change parameters. It handles pilot data loading, TPM filtering,
#' and gene sampling. This is a modularized version of the first part of extract_fc_expression_info.
#'
#' @param biological_system Character. Biological system for baseline expression. Available options:
#'   "K562", "A549", "THP-1", "T_CD8", "iPSC" (default: "K562").
#' @param B Integer. Number of Monte Carlo samples to generate when gene_list is NULL (default: 200).
#'   Ignored when gene_list is provided.
#' @param gene_list Character vector. Optional list of Ensembl gene IDs to use for analysis.
#'   If provided, expression parameters will be extracted for ALL specified genes (no sampling).
#'   If NULL (default), B genes are randomly sampled from baseline data.
#' @param TPM_threshold Numeric. Minimum TPM threshold (default: 10). Genes with expression
#'   levels below TPM_threshold/1e6 are filtered out before power calculation.
#' @param custom_pilot_data List. Optional custom pilot data. If provided,
#'   this data is used instead of the default biological_system data. Must contain
#'   baseline_expression_stats (data frame with gene expression data)
#'   and library_parameters (with UMI_per_cell and variation).
#'
#' @return A list with elements:
#' \describe{
#'   \item{expression_df}{Data frame with baseline expression parameters (response_id, relative_expression, expression_size)}
#'   \item{pilot_data}{Complete pilot data object for further use}
#'   \item{n_genes}{Integer number of genes in the processed dataset}
#' }
#'
#' @details
#' The function operates in two modes:
#' \itemize{
#'   \item \strong{Gene-specific mode} (gene_list provided): Uses ALL specified genes with importance sampling
#'   \item \strong{Random sampling mode} (gene_list = NULL): Randomly samples B genes from baseline
#' }
#'
#' Processing steps:
#' \enumerate{
#'   \item Load pilot data (custom or from package)
#'   \item Apply TPM threshold filtering
#'   \item Sample genes according to specified mode
#'   \item Return baseline expression data ready for fold change augmentation
#' }
#'
#' @examples
#' # Extract expression info from K562 system
#' expr_info <- extract_expression_info(
#'   biological_system = "K562",
#'   B = 100,
#'   TPM_threshold = 5
#' )
#'
#' # Examine the results
#' head(expr_info$expression_df)
#' print(paste("Extracted", expr_info$n_genes, "genes"))
#'
#' # Use with custom gene list
#' gene_list <- c("ENSG00000141510", "ENSG00000157764", "ENSG00000175899")
#' custom_expr <- extract_expression_info(
#'   biological_system = "K562",
#'   gene_list = gene_list,
#'   TPM_threshold = 1
#' )
#'
#' # Check expression levels
#' summary(custom_expr$expression_df$relative_expression)
#'
#' @keywords internal
extract_expression_info <- function(biological_system = "K562", B = 200, gene_list = NULL, TPM_threshold = 10, custom_pilot_data = NULL) {

  # set the random seed for reproducibility
  set.seed(1)

  ############## Load pilot data ################
  # Use custom pilot data if provided, otherwise load from data/ directory
  if (!is.null(custom_pilot_data)) {
    pilot_data <- custom_pilot_data
    # Try new key first, fall back to old key for backward compatibility
    baseline_expression_stats <- if (!is.null(custom_pilot_data$baseline_expression_stats)) {
      custom_pilot_data$baseline_expression_stats
    } else {
      custom_pilot_data$baseline_expression
    }
  } else {
    # Load complete pilot data from data/ directory based on biological_system
    pilot_data <- get_pilot_data_from_package(biological_system)
    # Try new key first, fall back to old key for backward compatibility
    baseline_expression_stats <- if (!is.null(pilot_data$baseline_expression_stats)) {
      pilot_data$baseline_expression_stats
    } else {
      pilot_data$baseline_expression
    }
  }

  # Handle both old nested structure and new simplified structure
  if (is.data.frame(baseline_expression_stats)) {
    # New simplified structure: baseline_expression_stats is directly a data frame
    baseline_df <- baseline_expression_stats
  } else {
    # Old nested structure: baseline_expression contains baseline_expression data frame
    baseline_df <- baseline_expression_stats$baseline_expression
  }

  #################### apply TPM threshold filtering FIRST ###################
  # Convert TPM threshold to relative expression scale (TPM / 1e6)
  TPM_threshold_relative <- TPM_threshold / 1e6

  # Filter the full baseline dataset by TPM threshold first
  if ("relative_expression" %in% colnames(baseline_df)) {
    pre_filter_n <- nrow(baseline_df)
    filtered_baseline_df <- baseline_df |>
      dplyr::filter(relative_expression >= TPM_threshold_relative)
    post_filter_n <- nrow(filtered_baseline_df)

    # Check if we have any genes left after filtering
    if (post_filter_n == 0) {
      stop("No genes remain after TPM threshold filtering. Consider lowering TPM_threshold.")
    }

    # Print filtering summary
    cat("TPM filtering: Kept", post_filter_n, "out of", pre_filter_n, "genes (threshold:", TPM_threshold, "TPM)\n")
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

  # return the baseline expression data without fold change augmentation
  result <- list(
    expression_df = expression_df,
    pilot_data = pilot_data,
    n_genes = n_genes
  )

  return(result)
}

#' Compute experimental cost for perturb-seq experiments
#'
#' @description
#' This function calculates the total cost of a perturb-seq experiment by combining
#' library preparation costs and sequencing costs based on the experimental platform,
#' sequencing platform, number of captured cells, and raw reads per cell.
#'
#' @param experimental_platform Character. The experimental platform used for single-cell capture.
#'   Currently supported: "10x Chromium v3" (default).
#' @param sequencing_platform Character. The sequencing platform used for RNA-seq.
#'   Currently supported: "NovaSeq X 25B" (default).
#' @param num_captured_cells Numeric. Number of captured cells in the experiment.
#' @param raw_reads_per_cell Numeric. Number of raw sequencing reads per cell.
#'
#' @return Numeric. Total experimental cost in USD combining library preparation and sequencing costs.
#'
#' @details
#' The cost calculation includes two main components:
#' \itemize{
#'   \item **Library preparation cost**: Based on cost per captured cell for the experimental platform
#'   \item **Sequencing cost**: Based on cost per million reads for the sequencing platform
#' }
#'
#' Current cost parameters:
#' \itemize{
#'   \item 10x Chromium v3: $0.086 per captured cell
#'   \item NovaSeq X 25B: $0.374 per million reads
#' }
#'
#' Total cost = (cost_per_captured_cell × num_captured_cells) +
#'              (cost_per_million_reads × raw_reads_per_cell × num_captured_cells / 1e6)
#'
#' @examples
#' # Calculate cost for a typical experiment
#' cost_result <- cost_computation(
#'   experimental_platform = "10x Chromium v3",
#'   sequencing_platform = "NovaSeq X 25B",
#'   num_captured_cells = 10000,
#'   raw_reads_per_cell = 50000
#' )
#' @keywords internal
cost_computation <- function(experimental_platform = "10x Chromium v3",
                             sequencing_platform = "NovaSeq X 25B",
                             num_captured_cells, raw_reads_per_cell){

  # Obtain cost per captured cell
  cost_per_captured_cell <- switch(experimental_platform,
    `10x Chromium v3` = 0.086
  )

  # Obtain cost per million reads
  cost_per_million_reads <- switch (sequencing_platform,
    `NovaSeq X 25B` = 0.374
  )

  # compute the total cost for library preparation
  library_preparation_cost <- cost_per_captured_cell * num_captured_cells

  # compute the sequencing cost
  sequencing_cost <- cost_per_million_reads * raw_reads_per_cell * num_captured_cells / 1e6

  # compute the total cost
  return(library_preparation_cost + sequencing_cost)

}

