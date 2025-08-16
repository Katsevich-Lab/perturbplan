# Script to create custom pilot data RDS file for PerturbPlan
# This example shows how to create a properly formatted RDS file for custom pilot data upload
# The output can be used as input for custom_pilot_data parameter in extract_expression_info() 
# and extract_fc_expression_info() functions
#
# Usage: Run this script from the perturbplan package directory
# R -e "source('inst/extdata/create_custom_pilot_data.R')"

#' Create custom pilot data from K562_10x dataset
#' 
#' This function demonstrates how to construct a custom pilot data object
#' that follows the standardized structure expected by PerturbPlan functions.
#' 
#' @param subset_size Integer. Number of genes to include in the subset (default: 1000)
#' @param seed Integer. Random seed for reproducible gene selection (default: 123)
#' 
#' @return List with baseline_expression_stats and library_parameters
create_custom_pilot_data <- function(subset_size = 1000, seed = 123) {
  
  cat("Creating custom pilot data from K562_10x dataset...\n")
  
  # Load the K562_10x data (now uses updated structure)
  data("K562_10x", package = "perturbplan")
  
  # Extract baseline expression stats and library parameters
  # Should now use the new structure after package reinstallation
  baseline_stats <- K562_10x$baseline_expression_stats
  library_params <- K562_10x$library_parameters
  
  cat("Original K562 data contains", nrow(baseline_stats), "genes\n")
  
  # Set seed for reproducible results
  set.seed(seed)
  
  # Create subset of genes (can be random or by specific criteria)
  if (subset_size >= nrow(baseline_stats)) {
    cat("Subset size larger than available genes, using all", nrow(baseline_stats), "genes\n")
    subset_indices <- seq_len(nrow(baseline_stats))
  } else {
    # Random subset for demonstration
    subset_indices <- sample(seq_len(nrow(baseline_stats)), subset_size)
    cat("Randomly selected", subset_size, "genes from", nrow(baseline_stats), "available\n")
  }
  
  # Create subset of baseline expression data
  baseline_subset <- baseline_stats[subset_indices, ]
  
  # Optionally modify the data (examples):
  # 1. Scale expression levels
  # baseline_subset$relative_expression <- baseline_subset$relative_expression * 1.2
  
  # 2. Filter by expression level (keep genes with TPM > 5)
  # tpm_values <- baseline_subset$relative_expression * 1e6
  # keep_high_expr <- tpm_values > 5
  # baseline_subset <- baseline_subset[keep_high_expr, ]
  # cat("After TPM filtering:", nrow(baseline_subset), "genes remain\n")
  
  # 3. Adjust library parameters if needed
  # library_params$UMI_per_cell <- library_params$UMI_per_cell * 0.8
  # library_params$variation <- library_params$variation * 1.1
  
  # Construct the custom pilot data in the expected format
  custom_pilot_data <- list(
    baseline_expression_stats = baseline_subset,
    library_parameters = library_params
  )
  
  # Validate the structure
  cat("Validating custom pilot data structure...\n")
  
  # Check required components
  required_components <- c("baseline_expression_stats", "library_parameters")
  has_components <- all(required_components %in% names(custom_pilot_data))
  
  if (!has_components) {
    stop("Custom pilot data missing required components: ", 
         paste(setdiff(required_components, names(custom_pilot_data)), collapse = ", "))
  }
  
  # Check baseline_expression_stats structure
  baseline_df <- custom_pilot_data$baseline_expression_stats
  required_cols <- c("response_id", "relative_expression", "expression_size")
  has_cols <- all(required_cols %in% colnames(baseline_df))
  
  if (!has_cols) {
    stop("baseline_expression_stats missing required columns: ",
         paste(setdiff(required_cols, colnames(baseline_df)), collapse = ", "))
  }
  
  # Check library_parameters structure
  library_df <- custom_pilot_data$library_parameters
  required_params <- c("UMI_per_cell", "variation")
  has_params <- all(required_params %in% names(library_df))
  
  if (!has_params) {
    stop("library_parameters missing required elements: ",
         paste(setdiff(required_params, names(library_df)), collapse = ", "))
  }
  
  # Report final structure
  cat("✓ Custom pilot data created successfully!\n")
  cat("✓ Structure: baseline_expression_stats (", nrow(baseline_df), " genes) + library_parameters\n")
  cat("✓ Library parameters: UMI_per_cell =", library_df$UMI_per_cell, 
      ", variation =", round(library_df$variation, 4), "\n")
  
  # Calculate some summary statistics
  tpm_values <- baseline_df$relative_expression * 1e6
  cat("✓ Expression range: TPM", round(min(tpm_values), 2), "-", round(max(tpm_values), 2), "\n")
  cat("✓ Median TPM:", round(median(tpm_values), 2), "\n")
  
  return(custom_pilot_data)
}

# Create different variants of custom pilot data
cat("=== Creating Custom Pilot Data Examples ===\n\n")

# Example 1: Small subset for testing (100 genes)
cat("1. Creating small test dataset (100 genes)...\n")
small_pilot_data <- create_custom_pilot_data(subset_size = 100, seed = 123)
saveRDS(small_pilot_data, "inst/extdata/custom_pilot_data_small.rds")
cat("   Saved to: inst/extdata/custom_pilot_data_small.rds\n\n")

# Example 2: Medium subset for analysis (1000 genes)
cat("2. Creating medium dataset (1000 genes)...\n")
medium_pilot_data <- create_custom_pilot_data(subset_size = 1000, seed = 456)
saveRDS(medium_pilot_data, "inst/extdata/custom_pilot_data_medium.rds")
cat("   Saved to: inst/extdata/custom_pilot_data_medium.rds\n\n")

# Example 3: Large subset for comprehensive analysis (5000 genes)
cat("3. Creating large dataset (5000 genes)...\n")
large_pilot_data <- create_custom_pilot_data(subset_size = 5000, seed = 789)
saveRDS(large_pilot_data, "inst/extdata/custom_pilot_data_large.rds")
cat("   Saved to: inst/extdata/custom_pilot_data_large.rds\n\n")

# Test that the custom pilot data works with PerturbPlan functions
cat("=== Testing Custom Pilot Data with PerturbPlan Functions ===\n\n")

cat("Testing extract_expression_info() with custom pilot data...\n")
result1 <- extract_expression_info(
  biological_system = "K562", 
  B = 5, 
  custom_pilot_data = small_pilot_data
)
cat("✓ extract_expression_info() succeeded - returned", nrow(result1$expression_df), "genes\n")

cat("Testing extract_fc_expression_info() with custom pilot data...\n")
result2 <- extract_fc_expression_info(
  minimum_fold_change = 1.5, 
  gRNA_variability = 0.1, 
  biological_system = "K562", 
  B = 5, 
  custom_pilot_data = small_pilot_data
)
cat("✓ extract_fc_expression_info() succeeded - returned", nrow(result2$fc_expression_df), "genes\n")

cat("\n=== Summary ===\n")
cat("Custom pilot data files created successfully!\n")
cat("Files available:\n")
cat("- custom_pilot_data_small.rds (100 genes) - for quick testing\n")
cat("- custom_pilot_data_medium.rds (1000 genes) - for typical analysis\n")
cat("- custom_pilot_data_large.rds (5000 genes) - for comprehensive analysis\n")
cat("\nUsage example:\n")
cat("custom_data <- readRDS('inst/extdata/custom_pilot_data_medium.rds')\n")
cat("result <- extract_expression_info(custom_pilot_data = custom_data, B = 200)\n")

cat("\nStructure of custom pilot data:\n")
cat("list(\n")
cat("  baseline_expression_stats = data.frame(\n")
cat("    response_id = c('ENSG00000...', ...),\n")
cat("    relative_expression = c(1.23e-05, ...),  # TPM/1e6 scale\n")
cat("    expression_size = c(0.45, ...)           # Dispersion parameters\n")
cat("  ),\n")
cat("  library_parameters = list(\n")
cat("    UMI_per_cell = 59158,     # From K562 data\n")
cat("    variation = 0.397         # From K562 data\n")
cat("  )\n")
cat(")\n")