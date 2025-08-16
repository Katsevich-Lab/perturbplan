# Script to create custom baseline expression RDS file for PerturbPlan
# This example shows how to create a properly formatted RDS file for custom baseline upload

library(perturbplan)

# Method 1: Create from scratch (minimal example)
create_minimal_baseline <- function() {
  # Example with 500 genes
  n_genes <- 500
  
  # Generate example gene IDs (Ensembl format)
  gene_ids <- paste0("ENSG", sprintf("%011d", 1:n_genes))
  
  # Generate realistic expression values
  # TPM values typically range from 0.1 to 10000, convert to relative scale (TPM/1e6)
  set.seed(123)
  tpm_values <- exp(rnorm(n_genes, mean = log(100), sd = 2))  # Log-normal distribution
  relative_expression <- tpm_values / 1e6
  
  # Generate expression size parameters (dispersion, typically 0.1 to 10)
  expression_size <- runif(n_genes, min = 0.1, max = 5.0)
  
  # Create baseline expression data frame
  baseline_df <- data.frame(
    response_id = gene_ids,
    relative_expression = relative_expression,
    expression_size = expression_size,
    stringsAsFactors = FALSE
  )
  
  # Create expression dispersion curve function
  # This function models the relationship between expression level and dispersion
  dispersion_curve <- function(v) {
    # Typical mean-variance relationship in scRNA-seq data
    pmax(0.01, 0.1 + 0.5 / sqrt(v))
  }
  
  # Create the simplified structure (new format)
  custom_baseline <- baseline_df
  
  return(custom_baseline)
}

# Method 2: Modify existing data (recommended)
create_from_existing <- function() {
  # Load default baseline data
  default_data <- get_pilot_data_from_package("K562")
  
  # Extract baseline expression stats (handles both old and new key names)
  baseline_stats <- if (!is.null(default_data$baseline_expression_stats)) {
    default_data$baseline_expression_stats
  } else {
    default_data$baseline_expression$baseline_expression
  }
  
  # Subset to specific genes of interest (example: first 1000 genes)
  subset_indices <- 1:1000
  subset_baseline <- baseline_stats[subset_indices, ]
  
  # You can modify expression values if needed
  # Example: Increase expression levels by 20%
  # subset_baseline$relative_expression <- subset_baseline$relative_expression * 1.2
  
  # Create modified baseline structure (new simplified format)
  custom_baseline <- subset_baseline
  
  return(custom_baseline)
}

# Create examples
if (interactive()) {
  # Generate minimal example
  minimal_baseline <- create_minimal_baseline()
  saveRDS(minimal_baseline, "minimal_custom_baseline.rds")
  
  # Generate example from existing data
  existing_baseline <- create_from_existing()
  saveRDS(existing_baseline, "example_custom_baseline.rds")
  
  # Validate both files
  cat("Minimal baseline validation:\n")
  result1 <- validate_custom_baseline_rds(minimal_baseline)
  cat("Valid:", result1$valid, "\n")
  
  cat("\nExisting baseline validation:\n")
  result2 <- validate_custom_baseline_rds(existing_baseline)
  cat("Valid:", result2$valid, "\n")
  
  cat("\nFiles created successfully!\n")
  cat("- minimal_custom_baseline.rds: Simple example with 500 genes\n")
  cat("- example_custom_baseline.rds: Subset of default K562 data\n")
}