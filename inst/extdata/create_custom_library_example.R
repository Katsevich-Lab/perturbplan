# Example script for creating custom library parameters for PerturbPlan
# This script demonstrates how to create a custom library RDS file
# for use with the custom library functionality in the Shiny app

library(perturbplan)

# Example 1: Creating custom library data from scratch
# =====================================================

# Define custom library parameters
# UMI_per_cell: Maximum UMI per cell parameter (typically 1000-50000)
# variation: Variation parameter for PCR bias (typically 0.1-1.0)

custom_library_example <- list(
  UMI_per_cell = 15000,  # Example: 15K UMI per cell maximum
  variation = 0.25       # Example: moderate variation
)

# Validate the structure
validation_result <- validate_custom_library_rds(custom_library_example)
if (validation_result$valid) {
  cat("Custom library data is valid!\n")
  cat("Summary:", validation_result$summary, "\n")
} else {
  cat("Validation errors:", paste(validation_result$errors, collapse = ", "), "\n")
}

# Save as RDS file
saveRDS(custom_library_example, "my_custom_library.rds")
cat("Custom library saved to 'my_custom_library.rds'\n")

# Example 2: Modifying default library parameters
# ===============================================

# Extract default library info for K562 (recommended approach)
fc_result <- extract_fc_expression_info(minimum_fold_change = -1, gRNA_variability = 0.25, biological_system = "K562", B = 10)
default_library <- fc_result$pilot_data$library_parameters
cat("Default K562 library parameters:\n")
cat("  UMI_per_cell:", default_library$UMI_per_cell, "\n")
cat("  variation:", default_library$variation, "\n")

# Create modified version
modified_library <- list(
  UMI_per_cell = default_library$UMI_per_cell * 1.2,  # 20% higher UMI capacity
  variation = default_library$variation * 0.8         # 20% lower variation
)

# Save modified version
saveRDS(modified_library, "modified_k562_library.rds")
cat("Modified library saved to 'modified_k562_library.rds'\n")

# Example 3: Different cell types with varying parameters
# =======================================================

# Example for a high-efficiency cell type
high_efficiency_library <- list(
  UMI_per_cell = 25000,  # Higher capture efficiency
  variation = 0.15       # Lower PCR bias
)

# Example for a low-efficiency cell type  
low_efficiency_library <- list(
  UMI_per_cell = 8000,   # Lower capture efficiency
  variation = 0.4        # Higher PCR bias
)

# Save different cell type examples
saveRDS(high_efficiency_library, "high_efficiency_library.rds")
saveRDS(low_efficiency_library, "low_efficiency_library.rds")

cat("High efficiency library saved to 'high_efficiency_library.rds'\n")
cat("Low efficiency library saved to 'low_efficiency_library.rds'\n")

# Usage instructions
cat("\n=== Usage Instructions ===\n")
cat("1. Upload any of these RDS files in the Shiny app under 'Library size parameters'\n")
cat("2. Select 'Custom' and choose your RDS file\n")
cat("3. The app will validate and display the parameters\n")
cat("4. Proceed with power analysis using your custom library parameters\n")

# Validation tips
cat("\n=== Validation Tips ===\n")
cat("- UMI_per_cell should be positive (typically 1000-50000)\n")
cat("- variation should be positive (typically 0.1-1.0)\n")
cat("- Both parameters must be single numeric values\n")
cat("- The RDS file must contain a list with exactly these two elements\n")