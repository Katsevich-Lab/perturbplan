#!/usr/bin/env Rscript

# Test script for cost_power_optimization function
# Run this script to test the cost optimization functionality

library(perturbplan)

# Load pilot data
cat("Loading pilot data...\n")
pilot_data <- get_pilot_data_from_package("K562")
baseline_expression_stats <- extract_expression_info(B = 1e3, tpm_threshold = 0,
                                                     biological_system = "K562")$expression_df

# Test 1: Optimize tpm_threshold with cost constraint
cat("Test 1: Optimizing tpm_threshold with cost constraint...\n")
tryCatch({
  result1 <- cost_power_optimization(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.85),
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = pilot_data$library_parameters,
    num_targets = 3000,
    MOI = 10,
    power_target = 0.8,
    cost_constraint = 60000,
    grid_size = 50
  )
  cat("✓ Test 1 PASSED: ", nrow(result1), " rows returned\n")
  cat("  Columns:", paste(names(result1), collapse = ", "), "\n")
  print(head(result1, 2))
}, error = function(e) {
  cat("✗ Test 1 FAILED:", e$message, "\n")
})

cat("\n" %+% strrep("-", 60) %+% "\n")

# Test 2: Optimize minimum_fold_change with cost constraint
cat("Test 2: Optimizing minimum_fold_change with cost constraint...\n")
tryCatch({
  result2 <- cost_power_optimization(
    minimizing_variable = "minimum_fold_change",
    fixed_variable = list(tpm_threshold = 50,
                          cells_per_target = 400,
                          reads_per_cell = 1e4),
    baseline_expression_stats = baseline_expression_stats,
    library_parameters = pilot_data$library_parameters,
    num_targets = 3000,
    MOI = 10,
    power_target = 0.8,
    cost_constraint = NULL,
    grid_size = 50
  )
  cat("✓ Test 2 PASSED: ", nrow(result2), " rows returned\n")
  cat("  Columns:", paste(names(result2), collapse = ", "), "\n")
  print(head(result2, 2))
}, error = function(e) {
  cat("✗ Test 2 FAILED:", e$message, "\n")
})

cat("\n" %+% strrep("-", 60) %+% "\n")

# Test 3: Optimize without cost constraint
cat("Test 3: Optimizing tmp_threshold without cost constraint...\n")
tryCatch({
  result3 <- cost_power_optimization(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(minimum_fold_change = 0.8),
    baseline_expression_stats = pilot_data$baseline_expression_stats,
    library_parameters = pilot_data$library_parameters,
    power_target = 0.8,
    cost_constraint = NULL,
    grid_size = 5
  )
  cat("✓ Test 3 PASSED: ", nrow(result3), " rows returned\n")
  cat("  Columns:", paste(names(result3), collapse = ", "), "\n")
  print(head(result3, 2))
}, error = function(e) {
  cat("✗ Test 3 FAILED:", e$message, "\n")
})

cat("\n" %+% strrep("-", 60) %+% "\n")

# Test 4: Fixed experimental design (no varying parameters)
cat("Test 4: Fixed experimental design optimization...\n")
tryCatch({
  result4 <- cost_power_optimization(
    minimizing_variable = "tpm_threshold",
    fixed_variable = list(
      minimum_fold_change = 0.8,
      cells_per_target = 1000,
      reads_per_cell = 10000
    ),
    baseline_expression_stats = pilot_data$baseline_expression_stats,
    library_parameters = pilot_data$library_parameters,
    power_target = 0.8,
    cost_constraint = 15000,
    grid_size = 3
  )
  cat("✓ Test 4 PASSED: ", nrow(result4), " rows returned\n")
  cat("  Columns:", paste(names(result4), collapse = ", "), "\n")
  print(head(result4, 2))
}, error = function(e) {
  cat("✗ Test 4 FAILED:", e$message, "\n")
})

cat("\n" %+% strrep("-", 60) %+% "\n")

# Test 5: Lower power target to avoid threshold issues
cat("Test 5: Lower power target optimization...\n")
tryCatch({
  result5 <- cost_power_optimization(
    minimizing_variable = "minimum_fold_change",
    fixed_variable = list(tpm_threshold = 100),
    baseline_expression_stats = pilot_data$baseline_expression_stats,
    library_parameters = pilot_data$library_parameters,
    power_target = 0.6,
    cost_constraint = 20000,
    grid_size = 3
  )
  cat("✓ Test 5 PASSED: ", nrow(result5), " rows returned\n")
  cat("  Columns:", paste(names(result5), collapse = ", "), "\n")
  print(head(result5, 2))
}, error = function(e) {
  cat("✗ Test 5 FAILED:", e$message, "\n")
})

cat("\n" %+% strrep("=", 60) %+% "\n")
cat("Test completed! Run the script as: Rscript test_cost_optimization.R\n")
cat("or source it in R: source('test_cost_optimization.R')\n")
