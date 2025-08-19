#!/usr/bin/env Rscript

# Simple test script for cost_power_computation function
cat("Testing cost_power_computation function...\n")

library(perturbplan)

# Load built-in data and create test parameters
cat("Loading test data...\n")
data("K562_10x", package = "perturbplan")
baseline_expression_stats <- extract_expression_info(B = 1e3, tpm_threshold = 0,
                                                     biological_system = "K562")$expression_df
library_parameters <- list(UMI_per_cell = 15000, variation = 0.25)

cat("Baseline expression data loaded:", nrow(baseline_expression_stats), "genes\n")

# Test 1: Basic tpm_threshold computation with cost constraint
cat("\nTest 1: Basic tmp_threshold computation with cost constraint...\n")
tryCatch({
  result1 <- cost_power_computation(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.8),
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters,
    num_targets = 50,
    power_target = 0.7,
    cost_constraint = 50000,
    grid_size = 3
  )
  cat("✓ Test 1 PASSED:", nrow(result1), "rows returned\n")
  cat("  Sample columns:", paste(names(result1)[1:5], collapse = ", "), "...\n")
  cat("  Cost range: $", round(min(result1$total_cost)), " - $", round(max(result1$total_cost)), "\n")
  cat("  Power range:", round(min(result1$overall_power), 3), " - ", round(max(result1$overall_power), 3), "\n")
}, error = function(e) {
  cat("✗ Test 1 FAILED:", e$message, "\n")
})

# Test 2: minimum_fold_change computation without cost constraint
cat("\nTest 2: minimum_fold_change computation without cost constraint...\n")
tryCatch({
  result2 <- cost_power_computation(
    minimizing_variable = "minimum_fold_change",
    fixed_variable = list(tpm_threshold = 100),
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters,
    num_targets = 50,
    power_target = 0.6,
    cost_constraint = NULL,
    grid_size = 3
  )
  cat("✓ Test 2 PASSED:", nrow(result2), "rows returned\n")
  cat("  Sample columns:", paste(names(result2)[1:5], collapse = ", "), "...\n")
  cat("  Power range:", round(min(result2$overall_power), 3), " - ", round(max(result2$overall_power), 3), "\n")
}, error = function(e) {
  cat("✗ Test 2 FAILED:", e$message, "\n")
})

# Test 3: Fixed experimental design
cat("\nTest 3: Fixed experimental design...\n")
tryCatch({
  result3 <- cost_power_computation(
    minimizing_variable = "minimum_fold_change",
    fixed_variable = list(
      tmp_threshold = 50,
      cells_per_target = 500,
      reads_per_cell = 8000
    ),
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = library_parameters,
    num_targets = 30,
    power_target = 0.6,
    cost_constraint = 30000,
    grid_size = 3
  )
  cat("✓ Test 3 PASSED:", nrow(result3), "rows returned\n")
  cat("  Fixed design - cells per target:", unique(result3$cells_per_target), "\n")
  cat("  Fixed design - reads per cell:", unique(result3$reads_per_cell), "\n")
}, error = function(e) {
  cat("✗ Test 3 FAILED:", e$message, "\n")
})

cat(paste0("\n", strrep("=", 60), "\n"))
cat("Testing completed! Run as: Rscript test_cost_computation_simple.R\n")