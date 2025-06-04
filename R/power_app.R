# Power calculation function using optimized separated Monte Carlo approach
# Integrates the new compute_power_grid_separated() function with Shiny app interface

#' Calculate power grid for app heatmap visualization
#' 
#' This function now uses the optimized separated Monte Carlo approach for
#' accurate and efficient power analysis while maintaining the interface
#' expected by the Shiny application.
#' 
#' @param num_targets Number of targets
#' @param gRNAs_per_target Number of gRNAs per target
#' @param non_targeting_gRNAs Number of non-targeting gRNAs
#' @param num_pairs Number of pairs analyzed
#' @param tpm_threshold Minimum TPM threshold
#' @param fdr_target FDR target level
#' @param fc_mean Fold-change mean
#' @param fc_sd Fold-change SD
#' @param prop_non_null Proportion of non-null pairs
#' @param MOI Multiplicity of infection
#' @param biological_system Biological system
#' @param experimental_platform Experimental platform
#' 
#' @return List with power grid, cell/read sequences, and parameters
#' @export
calculate_power_grid <- function(
  num_targets = 100,
  gRNAs_per_target = 4,
  non_targeting_gRNAs = 10,
  num_pairs = 1000,
  tpm_threshold = 10,
  fdr_target = 0.05,
  fc_mean = 0.85,
  fc_sd = 0.15,
  prop_non_null = 0.1,
  MOI = 10,
  biological_system = "K562",
  experimental_platform = "10x Chromium v3"
) {
  
  # Create grid for heatmap visualization
  cells_seq <- round(seq(5000, 50000, length.out = 20))
  reads_seq <- round(seq(2000, 50000, length.out = 20))
  
  # Create cells-reads data frame for power calculation
  cells_reads_df <- expand.grid(
    num_total_cells = cells_seq,
    reads_per_cell = reads_seq
  )
  
  # Call the optimized separated Monte Carlo power function
  power_results <- compute_power_grid_separated(
    cells_reads_df = cells_reads_df,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    num_pairs = num_pairs,
    tpm_threshold = tpm_threshold,
    fdr_target = fdr_target,
    fc_mean = fc_mean,
    fc_sd = fc_sd,
    prop_non_null = prop_non_null,
    MOI = MOI,
    biological_system = biological_system,
    experimental_platform = experimental_platform,
    side = "left",  # Standard for CRISPRi knockdown experiments
    control_group = "complement",  # Use complement cells as control
    B = 100,  # Monte Carlo samples for good accuracy vs speed balance
    fc_curve_points = 10,  # Sufficient resolution for curves
    expr_curve_points = 10
  )
  
  # Transform to expected format for heatmap
  power_grid <- data.frame(
    cells = power_results$num_total_cells,
    reads = power_results$reads_per_cell,
    power = power_results$overall_power
  )
  
  # Calculate expected discoveries
  expected_discoveries <- sum(power_grid$power) * prop_non_null
  
  # Return structure compatible with Shiny app
  list(
    power_grid = power_grid,
    cells_seq = cells_seq,
    reads_seq = reads_seq,
    expected_discoveries = expected_discoveries,
    parameters = list(
      num_targets = num_targets,
      gRNAs_per_target = gRNAs_per_target,
      non_targeting_gRNAs = non_targeting_gRNAs,
      num_pairs = num_pairs,
      tpm_threshold = tpm_threshold,
      fdr_target = fdr_target,
      fc_mean = fc_mean,
      fc_sd = fc_sd,
      prop_non_null = prop_non_null,
      MOI = MOI,
      biological_system = biological_system,
      experimental_platform = experimental_platform
    ),
    # Additional data for potential future use
    power_curves = list(
      fc_curves = power_results$power_by_fc,
      expr_curves = power_results$power_by_expr
    )
  )
}