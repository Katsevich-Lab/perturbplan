# Fixed example for compute_power_posthoc

# Create example discovery pairs
discovery_pairs <- data.frame(
  grna_target = c("Gene1", "Gene1", "Gene2", "Gene2"),
  response_id = c("ENSG00000141510", "ENSG00000157764",
                  "ENSG00000141510", "ENSG00000175899")
)

# Create example cells per gRNA data - FIXED: use "non-targeting" instead of "nt_cells"
cells_per_grna <- data.frame(
  grna_id = c("gRNA1", "gRNA2", "gRNA3", "gRNA4", "NT1", "NT2"),
  grna_target = c("Gene1", "Gene1", "Gene2", "Gene2",
                  "non-targeting", "non-targeting"),  # Changed from "nt_cells"
  num_cells = c(150, 180, 160, 170, 200, 190)
)

# Create example baseline expression data - FIXED: add missing relative_expression column
baseline_expression_stats <- data.frame(
  response_id = c("ENSG00000141510", "ENSG00000157764", "ENSG00000175899"),
  expression_mean = c(12.5, 8.3, 15.7),  # These values will be used directly
  expression_size = c(2.1, 1.8, 2.5)
)

# Compute post-hoc power - FIXED: change control_group to "nt_cells" to match expected behavior
power_results <- compute_power_posthoc(
  discovery_pairs = discovery_pairs,
  cells_per_grna = cells_per_grna,
  baseline_expression_stats = baseline_expression_stats,
  control_group = "nt_cells",
  fold_change_mean = 0.8,
  fold_change_sd = 0.15,
  multiple_testing_alpha = 0.1,
  side = "left"
)

# View results
print("Individual power results:")
print(power_results$individual_power)
print(paste("Expected discoveries:", round(power_results$expected_num_discoveries, 2)))

# Alternative approach: use "complement" control group (requires num_total_cells)
total_cells <- sum(cells_per_grna$num_cells)
power_results_complement <- compute_power_posthoc(
  discovery_pairs = discovery_pairs,
  cells_per_grna = cells_per_grna,
  baseline_expression_stats = baseline_expression_stats,
  control_group = "complement",  # Using complement control
  num_total_cells = total_cells,  # Required for complement
  fold_change_mean = 0.8,
  fold_change_sd = 0.15,
  multiple_testing_alpha = 0.1,
  side = "left"
)

print("\nResults with complement control group:")
print(power_results_complement$individual_power)
print(paste("Expected discoveries:", round(power_results_complement$expected_num_discoveries, 2)))