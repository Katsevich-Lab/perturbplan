# Create example combined pilot data file
# This script demonstrates how to create a combined pilot data RDS file
# that includes both baseline expression and library parameters

# Load existing K562 data directly from RDS files
baseline_rds_path <- system.file("extdata/baseline_expression", "Gasperini_expression.rds", package = "perturbplan", mustWork = TRUE)
library_rds_path <- system.file("extdata/library_info", "Gasperini_library.rds", package = "perturbplan", mustWork = TRUE)

# If running from package source directory, use relative paths
if (!file.exists(baseline_rds_path)) {
  baseline_rds_path <- "inst/extdata/baseline_expression/Gasperini_expression.rds"
  library_rds_path <- "inst/extdata/library_info/Gasperini_library.rds"
}

# Load the data directly
baseline_expression_raw <- readRDS(baseline_rds_path)
library_info_raw <- readRDS(library_rds_path)

# Format baseline expression data (simplified structure)
baseline_data <- baseline_expression_raw$baseline_expression

# Format library data (same as extract_library_info)
params <- library_info_raw$S_M_curve_params
library_data <- list(
  UMI_per_cell = unname(as.numeric(params[["UMI_per_cell"]])),
  variation = unname(as.numeric(params[["variation"]]))
)

# Combine into the expected structure with informative names
combined_pilot_data <- list(
  baseline_expression = baseline_data,
  library_parameters = library_data
)

# Simple validation of structure
if (is.list(combined_pilot_data) && 
    "baseline_expression" %in% names(combined_pilot_data) &&
    "library_parameters" %in% names(combined_pilot_data)) {
  cat("✓ Combined pilot data structure looks correct\n")
  cat("Components:", paste(names(combined_pilot_data), collapse = ", "), "\n")
  
  # Show sample info
  n_genes <- nrow(combined_pilot_data$baseline_expression)
  cat("Baseline expression: ", n_genes, " genes\n")
  cat("Library parameters: UMI_per_cell =", combined_pilot_data$library_parameters$UMI_per_cell, 
      ", variation =", combined_pilot_data$library_parameters$variation, "\n")
} else {
  cat("✗ Structure validation failed\n")
}

# Save the example file
output_path <- "inst/extdata/example_combined_pilot_data.rds"
saveRDS(combined_pilot_data, output_path)
cat("Saved example combined pilot data to:", output_path, "\n")

# Example of how users would create their own custom combined data:
cat("\n# Example user code to create custom combined pilot data:\n")
cat("combined_pilot_data <- list(\n")
cat("  baseline_expression = data.frame(\n")
cat("    response_id = c('ENSG00000141510', 'ENSG00000157764', ...),\n")
cat("    relative_expression = c(1.23e-05, 4.56e-06, ...),  # TPM/1e6 scale\n")
cat("    expression_size = c(0.45, 1.23, ...)              # Dispersion parameters\n")
cat("  ),\n")
cat("  library_parameters = list(\n")
cat("    UMI_per_cell = 15000,    # Maximum UMI per cell\n")
cat("    variation = 0.25         # PCR bias variation\n")
cat("  )\n")
cat(")\n")
cat("saveRDS(combined_pilot_data, 'my_custom_pilot_data.rds')\n")