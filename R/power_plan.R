# Power calculation function using optimized approach for Shiny app integration

# Suppress R CMD check warnings for variables used in dplyr contexts
utils::globalVariables(c("library_size", "num_total_cells", "reads_per_cell", "num_trt_cells", "num_cntrl_cells"))

#' Calculate power grid for app heatmap visualization (lightweight)
#'
#' This function provides lightweight power analysis for the Shiny application heatmap.
#' It only computes overall power for each cell/read combination, not detailed curves.
#' Use calculate_power_curves() for detailed curve computation on selected tiles.
#'
#' @param num_targets Number of targets
#' @param gRNAs_per_target Number of gRNAs per target
#' @param non_targeting_gRNAs Number of non-targeting gRNAs
#' @param fdr_target FDR target level
#' @param prop_non_null Proportion of non-null pairs
#' @param MOI Multiplicity of infection
#' @param side Test sidedness ("left", "right", "both")
#' @param control_group Control group type ("complement" or "nt_cells")
#' @param fc_expression_info List from extract_fc_expression_info() containing fc_expression_df and expression_dispersion_curve
#' @param library_info List from extract_library_info() containing UMI_per_cell and variation parameters
#' @param cells_reads_df Optional data frame with pre-computed experimental design containing
#'   num_total_cells, reads_per_cell, num_trt_cells, num_cntrl_cells columns. 
#'   If NULL (default), uses hardcoded grid for backward compatibility.
#'
#' @return List with power grid, cell/read sequences, and parameters
#' @export
calculate_power_grid <- function(
  fc_expression_info,
  library_info,
  num_targets = 100,
  gRNAs_per_target = 4,
  non_targeting_gRNAs = 10,
  fdr_target = 0.05,
  prop_non_null = 0.1,
  MOI = 10,
  side = "left",
  control_group = "complement",
  cells_reads_df = NULL
) {

  # Store original parameter to know if it was provided
  cells_reads_df_provided <- !is.null(cells_reads_df)

  # Use provided experimental design or create default grid
  if (is.null(cells_reads_df)) {
    # Create default grid for heatmap visualization (backward compatibility)
    cells_seq <- round(seq(5000, 50000, length.out = 10))
    reads_seq <- round(seq(2000, 50000, length.out = 10))

    # Create cells-reads data frame for power calculation
    cells_reads_df <- expand.grid(
      num_total_cells = cells_seq,
      reads_per_cell = reads_seq
    ) |>
      dplyr::mutate(
        num_trt_cells = gRNAs_per_target * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs),
        num_cntrl_cells = switch(control_group,
          complement = num_total_cells - num_trt_cells,
          nt_cells = non_targeting_gRNAs * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
        )
      )
  } else {
    # Validate provided experimental design
    required_cols <- c("num_total_cells", "reads_per_cell", "num_trt_cells", "num_cntrl_cells")
    missing_cols <- setdiff(required_cols, colnames(cells_reads_df))
    if (length(missing_cols) > 0) {
      stop("Missing required columns in cells_reads_df: ", paste(missing_cols, collapse = ", "))
    }
  }

  # Call the lightweight power function (overall power only, no curves)
  power_results <- compute_power_grid_overall(
    cells_reads_df = cells_reads_df,
    fc_expression_info = fc_expression_info,
    library_info = library_info,
    fdr_target = fdr_target,
    prop_non_null = prop_non_null,
    side = side
  )

  # Transform to expected format for heatmap
  power_grid <- data.frame(
    cells = round(power_results$num_trt_cells),  # Display treatment cells as integers
    reads = power_results$reads_per_cell,
    power = power_results$overall_power
  )

  # Return structure compatible with Shiny app
  list(
    power_grid = power_grid,
    cells_seq = sort(unique(power_grid$cells)),  # Treatment cell sequence for axis
    reads_seq = if (!cells_reads_df_provided) reads_seq else sort(unique(power_results$reads_per_cell)),
    parameters = list(
      num_targets = num_targets,
      gRNAs_per_target = gRNAs_per_target,
      non_targeting_gRNAs = non_targeting_gRNAs,
      fdr_target = fdr_target,
      prop_non_null = prop_non_null,
      MOI = MOI
    )
  )
}

#' Calculate power curves for specific cell/read combinations
#'
#' This function computes detailed power curves (TPM and fold-change) for specific
#' selected cell/read combinations. More efficient than computing all curves upfront.
#'
#' @param selected_tiles Data frame with columns 'cells' and 'reads' for selected tiles
#' @param fc_expression_info List from extract_fc_expression_info() containing fc_expression_df and expression_dispersion_curve
#' @param library_info List from extract_library_info() containing UMI_per_cell and variation parameters
#' @param num_targets Number of targets
#' @param gRNAs_per_target Number of gRNAs per target
#' @param non_targeting_gRNAs Number of non-targeting gRNAs
#' @param tpm_threshold Minimum TPM threshold
#' @param fdr_target FDR target level
#' @param prop_non_null Proportion of non-null pairs
#' @param MOI Multiplicity of infection
#' @param side Test sidedness ("left", "right", "both")
#' @param control_group Control group type ("complement" or "nt_cells")
#'
#' @return List with power curves for selected tiles
#' @export
calculate_power_curves <- function(
  selected_tiles,
  fc_expression_info,
  library_info,
  num_targets = 100,
  gRNAs_per_target = 4,
  non_targeting_gRNAs = 10,
  tpm_threshold = 10,
  fdr_target = 0.05,
  prop_non_null = 0.1,
  MOI = 10,
  side = "left",
  control_group = "complement"
) {

  # Create cells_reads_df for selected tiles only
  cells_reads_df <- data.frame(
    num_total_cells = selected_tiles$cells,
    reads_per_cell = selected_tiles$reads
  ) |>
    dplyr::mutate(
      num_trt_cells = gRNAs_per_target * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs),
      num_cntrl_cells = switch(control_group,
        complement = num_total_cells - num_trt_cells,
        nt_cells = non_targeting_gRNAs * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
      )
    )

  # Call the detailed power function for selected tiles only
  power_results <- compute_power_grid_full(
    cells_reads_df = cells_reads_df,
    fc_expression_info = fc_expression_info,
    library_info = library_info,
    tpm_threshold = tpm_threshold,
    fdr_target = fdr_target,
    prop_non_null = prop_non_null,
    side = side,
    fc_curve_points = 10,  # Sufficient resolution for curves
    expr_curve_points = 10
  )

  # Return power curves with cell/read information
  list(
    power_curves = list(
      fc_curves = power_results$power_by_fc,
      expr_curves = power_results$power_by_expr
    ),
    tiles_info = selected_tiles
  )
}

#' Compute power grid for overall power analysis (no curves)
#'
#' This function computes only overall power for each cell/read combination
#' without the expensive curve calculations. Used for heatmap generation.
#' Experimental design parameters (MOI, control group, etc.) should be pre-computed
#' into treatment and control cell counts in the cells_reads_df.
#'
#' @param cells_reads_df Data frame with columns num_total_cells, reads_per_cell, num_trt_cells, and num_cntrl_cells
#' @param fc_expression_info List from extract_fc_expression_info() containing fc_expression_df and expression_dispersion_curve
#' @param library_info List from extract_library_info() containing UMI_per_cell and variation parameters
#' @param fdr_target Target false discovery rate
#' @param prop_non_null Proportion of non-null hypotheses
#' @param side Test sidedness ("left", "right", "both")
#' @return Data frame with power analysis results (overall power only)
#' @export
compute_power_grid_overall <- function(
    cells_reads_df,
    fc_expression_info,
    library_info,
    fdr_target = 0.05,
    prop_non_null = 0.1,
    side = "left"
){

  ########################## compute the library size ##########################
  # Check if library_size is already provided (from identify_cell_read_range)
  if (!"library_size" %in% colnames(cells_reads_df)) {
    # Backward compatibility: compute library size if not provided
    UMI_per_cell <- library_info$UMI_per_cell
    variation <- library_info$variation
    
    cells_reads_df <- cells_reads_df |>
      dplyr::mutate(
        library_size = fit_read_UMI_curve(reads_per_cell = reads_per_cell, UMI_per_cell = !!UMI_per_cell, variation = !!variation)
      )
  }
  # If library_size already exists, use pre-computed values (more efficient!)

  ############### compute the power for the cells-reads grid ###################
  power_df <- cells_reads_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      overall_power = compute_power_plan_overall(
        # experimental information
        num_trt_cells = num_trt_cells, num_cntrl_cells = num_cntrl_cells, library_size = library_size,
        # analysis information
        multiple_testing_alpha = fdr_target, multiple_testing_method = "BH",
        side = side,
        # separated approach information
        fc_expression_df = fc_expression_info$fc_expression_df,
        prop_non_null = prop_non_null
      )
    )

  # return the output dataframe
  return(power_df)
}





#' Compute full power grid with curves using C++ test statistic computation
#'
#' This function uses C++ implementations for computational efficiency in power analysis
#' for perturb-seq experiments across different experimental conditions. Experimental
#' design parameters should be pre-computed into treatment and control cell counts.
#'
#' @param cells_reads_df Data frame with columns num_total_cells, reads_per_cell, num_trt_cells, and num_cntrl_cells
#' @param fc_expression_info List from extract_fc_expression_info() containing fc_expression_df and expression_dispersion_curve
#' @param library_info List from extract_library_info() containing UMI_per_cell and variation parameters
#' @param tpm_threshold TPM threshold for expression curve generation
#' @param fdr_target Target false discovery rate
#' @param prop_non_null Proportion of non-null hypotheses
#' @param side Test sidedness ("left", "right", "both")
#' @param fc_curve_points Number of points for fold change curve
#' @param expr_curve_points Number of points for expression curve
#' @return Data frame with power analysis results
#' @export
compute_power_grid_full <- function(
    cells_reads_df,
    fc_expression_info,
    library_info,
    tpm_threshold = 10,
    fdr_target = 0.05,
    prop_non_null = 0.1,
    side = "left",
    fc_curve_points = 10,
    expr_curve_points = 10
){

  ############### extract components from fc_expression_info ##################
  set.seed(1)  # Reproducible results
  fc_expression_df <- fc_expression_info$fc_expression_df
  expression_dispersion_curve <- fc_expression_info$expression_dispersion_curve
  fold_change_mean <- fc_expression_info$fold_change_mean

  # Define systematic output grids
  fc_range <- range(fc_expression_df$fold_change)
  expr_range <- range(fc_expression_df$relative_expression)

  # Use adaptive fold change grid based on fold_change_mean and test sidedness
  if (side == "left") {
    # For left-sided tests (knockdown), focus on FC < 1
    if (fold_change_mean < 1) {
      # Mean suggests knockdown effects, focus on range below 1
      fc_output_grid <- seq(min(fc_range[1], 0.5), 1, length.out = fc_curve_points)
    } else {
      # Mean > 1 but testing left side, still focus below 1 but start from mean
      fc_output_grid <- seq(min(fc_range[1], fold_change_mean - 0.5), 1, length.out = fc_curve_points)
    }
  } else if (side == "right") {
    # For right-sided tests (overexpression), focus on FC > 1
    if (fold_change_mean > 1) {
      # Mean suggests overexpression effects, focus on range above 1
      fc_output_grid <- seq(1, max(fc_range[2], fold_change_mean + 0.5), length.out = fc_curve_points)
    } else {
      # Mean < 1 but testing right side, still focus above 1
      fc_output_grid <- seq(1, max(fc_range[2], 2), length.out = fc_curve_points)
    }
  } else {
    # For both-sided tests, use same grid as either left or right based on fold_change_mean
    if (fold_change_mean < 1) {
      # Use left-side grid (knockdown range)
      fc_output_grid <- seq(min(fc_range[1], 0.5), 1, length.out = fc_curve_points)
    } else {
      # Use right-side grid (overexpression range)
      fc_output_grid <- seq(1, max(fc_range[2], fold_change_mean + 0.5), length.out = fc_curve_points)
    }
  }

  # Use log-spaced points for gene expression evaluation grid, starting from TPM threshold
  # Convert TPM threshold to relative expression scale (TPM / 1e6)
  tpm_threshold_relative <- tpm_threshold / 1e6
  expr_min <- max(expr_range[1], tpm_threshold_relative)  # Start from TPM threshold or data minimum, whichever is higher
  expr_output_grid <- 10^seq(log10(expr_min), log10(expr_range[2]), length.out = expr_curve_points)

  ########################## compute the library size ##########################
  # Check if library_size is already provided (from identify_cell_read_range)
  if (!"library_size" %in% colnames(cells_reads_df)) {
    # Backward compatibility: compute library size if not provided
    UMI_per_cell <- library_info$UMI_per_cell
    variation <- library_info$variation
    
    cells_reads_df <- cells_reads_df |>
      dplyr::mutate(
        library_size = fit_read_UMI_curve(reads_per_cell = reads_per_cell, UMI_per_cell = !!UMI_per_cell, variation = !!variation)
      )
  }
  # If library_size already exists, use pre-computed values (more efficient!)

  ############### compute the power for the cells-reads grid ###################
  power_df <- cells_reads_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      power_output = list(
        compute_power_plan_full(
          # experimental information
          num_trt_cells = num_trt_cells, num_cntrl_cells = num_cntrl_cells, library_size = library_size,
          # analysis information
          multiple_testing_alpha = fdr_target, multiple_testing_method = "BH",
          side = side,
          # separated approach information
          fc_expression_df = fc_expression_df,
          expression_dispersion_curve = expression_dispersion_curve,
          fc_output_grid = fc_output_grid,
          expr_output_grid = expr_output_grid,
          prop_non_null = prop_non_null
        )
      ),
      overall_power = power_output$overall_power,
      power_by_fc = list(power_output$power_by_fc),
      power_by_expr = list(power_output$power_by_expr)
    ) |>
    dplyr::select(-power_output)

  # return the output dataframe
  return(power_df)
}



#' Compute full power analysis with curves using C++ Monte Carlo
#'
#' This function computes detailed power curves using C++ implementations for 
#' improved performance. Treatment and control cell counts should be pre-computed.
#'
#' @param num_trt_cells Number of treatment cells (pre-computed)
#' @param num_cntrl_cells Number of control cells (pre-computed)
#' @param library_size Library size (reads per cell)
#' @param multiple_testing_alpha Alpha level for multiple testing
#' @param multiple_testing_method Multiple testing method
#' @param side Test sidedness
#' @param fc_expression_df Data frame with fold change and expression info
#' @param expression_dispersion_curve Function for expression-size relationship
#' @param fc_output_grid Grid points for fold change curve
#' @param expr_output_grid Grid points for expression curve
#' @param prop_non_null Proportion of non-null hypotheses
#' @return List with overall power and power curves
#' @export
compute_power_plan_full <- function(
    # experimental information
  num_trt_cells, num_cntrl_cells, library_size,
  # analysis information
  multiple_testing_alpha = 0.05, multiple_testing_method = "BH", side = "left",
  # separated approach information
  fc_expression_df, expression_dispersion_curve, fc_output_grid, expr_output_grid, prop_non_null = 0.1){

  ################ compute shared results (overall power, cutoff, cell counts) ################
  lightweight_results <- compute_power_plan_overall(
    # experimental information
    num_trt_cells = num_trt_cells, num_cntrl_cells = num_cntrl_cells, library_size = library_size,
    # analysis information
    multiple_testing_alpha = multiple_testing_alpha, multiple_testing_method = multiple_testing_method,
    side = side,
    # separated approach information
    fc_expression_df = fc_expression_df, prop_non_null = prop_non_null,
    return_full_results = TRUE
  )
  
  # Extract results from lightweight computation
  overall_power <- lightweight_results$overall_power
  sig_cutoff <- lightweight_results$sig_cutoff
  num_trt_cells <- lightweight_results$num_trt_cells
  num_cntrl_cells <- lightweight_results$num_cntrl_cells

  ##############################################################################
  #################### compute output curves ###################################

  # Compute fold change curve using C++ function for performance
  power_by_fc <- compute_fc_curve_cpp(
    fc_output_grid = fc_output_grid,
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    side = side,
    cutoff = sig_cutoff
  )

  # Compute expression curve using C++ function for performance
  power_by_expr <- compute_expression_curve_cpp(
    expr_output_grid = expr_output_grid,
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    expression_dispersion_curve = expression_dispersion_curve,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    side = side,
    cutoff = sig_cutoff
  )

  # return the output
  return(list(
    overall_power = overall_power,
    power_by_fc = power_by_fc,
    power_by_expr = power_by_expr,
    grid_summary = list(
      mc_samples = nrow(fc_expression_df),
      fc_output_points = length(fc_output_grid),
      expr_output_points = length(expr_output_grid),
      total_computations = nrow(fc_expression_df) * (1 + length(fc_output_grid) + length(expr_output_grid)),
      cutoff = sig_cutoff
    )
  ))
}

#' Compute overall power for power analysis (lightweight computation)
#'
#' This function computes overall power and BH cutoff without expensive curve calculations.
#' Used for heatmap generation and as a base for detailed curve computation.
#' Treatment and control cell counts should be pre-computed.
#'
#' @param num_trt_cells Number of treatment cells (pre-computed)
#' @param num_cntrl_cells Number of control cells (pre-computed)
#' @param library_size Library size (reads per cell)
#' @param multiple_testing_alpha Alpha level for multiple testing
#' @param multiple_testing_method Multiple testing method
#' @param side Test sidedness
#' @param fc_expression_df Data frame with fold change and expression info
#' @param prop_non_null Proportion of non-null hypotheses
#' @param return_full_results If TRUE, return list with all intermediate results; if FALSE, return only overall power
#' @return Overall power value (scalar) or list with full results depending on return_full_results
#' @export
compute_power_plan_overall <- function(
    # experimental information
  num_trt_cells, num_cntrl_cells, library_size,
  # analysis information
  multiple_testing_alpha = 0.05, multiple_testing_method = "BH", side = "left",
  # separated approach information
  fc_expression_df, prop_non_null = 0.1, return_full_results = FALSE){

  # Call the C++ implementation for improved performance
  return(compute_power_plan_overall_cpp(
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    num_trt_cells = num_trt_cells,
    num_cntrl_cells = num_cntrl_cells,
    multiple_testing_alpha = multiple_testing_alpha,
    multiple_testing_method = multiple_testing_method,
    side = side,
    prop_non_null = prop_non_null,
    return_full_results = return_full_results
  ))
}

