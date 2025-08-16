########################## New workflow ########################################
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
#' @param tpm_threshold Numeric. Minimum TPM threshold (default: 10). Genes with expression
#'   levels below tpm_threshold/1e6 are filtered out before power calculation.
#' @param custom_pilot_data List. Optional custom pilot data. If provided,
#'   this data is used instead of the default biological_system data. Must contain
#'   baseline_expression (with baseline_expression data frame and expression_dispersion_curve function)
#'   and library_parameters (with UMI_per_cell and variation).
#'
#' @return A list with elements:
#' \describe{
#'   \item{expression_df}{Data frame with baseline expression parameters (response_id, relative_expression, expression_size)}
#'   \item{expression_dispersion_curve}{Function relating expression mean to dispersion}
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
#' @export
extract_expression_info <- function(biological_system = "K562", B = 200, gene_list = NULL, tpm_threshold = 10, custom_pilot_data = NULL) {

  # set the random seed for reproducibility
  set.seed(1)

  ############## Load pilot data ################
  # Use custom pilot data if provided, otherwise load from data/ directory
  if (!is.null(custom_pilot_data)) {
    pilot_data <- custom_pilot_data
    baseline_expression_stats <- custom_pilot_data$baseline_expression
  } else {
    # Load complete pilot data from data/ directory based on biological_system
    pilot_data <- perturbplan:::get_pilot_data_from_package(biological_system)
    baseline_expression_stats <- pilot_data$baseline_expression
  }
  baseline_df <- baseline_expression_stats$baseline_expression

  #################### apply TPM threshold filtering FIRST ###################
  # Convert TPM threshold to relative expression scale (TPM / 1e6)
  tpm_threshold_relative <- tpm_threshold / 1e6

  # Filter the full baseline dataset by TPM threshold first
  if ("relative_expression" %in% colnames(baseline_df)) {
    pre_filter_n <- nrow(baseline_df)
    filtered_baseline_df <- baseline_df |>
      dplyr::filter(relative_expression >= tpm_threshold_relative)
    post_filter_n <- nrow(filtered_baseline_df)

    # Check if we have any genes left after filtering
    if (post_filter_n == 0) {
      stop("No genes remain after TPM threshold filtering. Consider lowering tpm_threshold.")
    }

    # Print filtering summary
    cat("TPM filtering: Kept", post_filter_n, "out of", pre_filter_n, "genes (threshold:", tpm_threshold, "TPM)\n")
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

  ################## extract the expression-dispersion curve ###################
  expression_dispersion_curve <- baseline_expression_stats$expression_dispersion_curve

  # return the baseline expression data without fold change augmentation
  result <- list(
    expression_df = expression_df,
    expression_dispersion_curve = expression_dispersion_curve,
    pilot_data = pilot_data,
    n_genes = n_genes
  )

  return(result)
}

#' Compute power analysis for experimental design grid
#'
#' @description
#' Simplified version of identify_cell_read_range that returns a clean dataframe
#' with experimental design combinations and their corresponding power values.
#'
#' @param fc_expression_df Data frame with fold change and expression information.
#' @param library_info List containing UMI_per_cell and variation parameters.
#' @param grid_size Integer. Number of points in each dimension of the grid (default: 10).
#' @param min_power_threshold Numeric. Minimum power threshold for cell range determination (default: 0.01).
#' @param max_power_threshold Numeric. Maximum power threshold for cell range determination (default: 0.8).
#' @param MOI Numeric. Multiplicity of infection for cell allocation calculations (default: 10).
#' @param num_targets Integer. Number of targets for cell allocation calculations (default: 100).
#' @param gRNAs_per_target Integer. Number of gRNAs per target (default: 4).
#' @param non_targeting_gRNAs Integer. Number of non-targeting gRNAs (default: 10).
#' @param control_group String. Control group type: "complement" or "nt_cells" (default: "complement").
#' @param multiple_testing_alpha Numeric. Alpha level for multiple testing (default: 0.05).
#' @param side String. Test sidedness: "left", "right", or "both" (default: "left").
#' @param prop_non_null Numeric. Proportion of non-null hypotheses (default: 0.1).
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{cells_per_target}{Numeric. Number of treatment cells per target}
#'   \item{num_total_cells}{Numeric. Total number of cells in the experiment}
#'   \item{reads_per_cell}{Numeric. Sequencing reads per cell}
#'   \item{library_size}{Numeric. Effective library size (UMIs)}
#'   \item{overall_power}{Numeric. Statistical power for this experimental design}
#' }
#'
#' @details
#' This function simplifies the experimental design process by:
#' \enumerate{
#'   \item Determining optimal reads per cell range using library size curves
#'   \item Determining optimal cell count range based on power thresholds
#'   \item Creating logarithmically-spaced grids for both dimensions
#'   \item Computing power for all combinations
#'   \item Returning a clean dataframe ready for analysis
#' }
#'
#' @export
compute_power_plan_per_grid <- function(
  fc_expression_df,
  library_info,
  grid_size = 10,
  min_power_threshold = 0.01,
  max_power_threshold = 0.8,
  MOI = 10,
  num_targets = 100,
  gRNAs_per_target = 4,
  non_targeting_gRNAs = 10,
  control_group = "complement",
  multiple_testing_alpha = 0.05,
  side = "left",
  prop_non_null = 0.1
) {

  # Extract needed data
  UMI_per_cell <- library_info$UMI_per_cell
  variation <- library_info$variation

  # Step 1: Determine reads per cell range using library size curves
  reads_range <- identify_reads_range_cpp(
    UMI_per_cell = UMI_per_cell,
    variation = variation
  )

  min_reads_per_cell <- reads_range$min_reads_per_cell
  max_reads_per_cell <- reads_range$max_reads_per_cell

  # Step 2: Determine cell count range based on power thresholds
  cell_range <- identify_cell_range_cpp(
    min_reads_per_cell = min_reads_per_cell,
    max_reads_per_cell = max_reads_per_cell,
    fc_expression_df = fc_expression_df,
    UMI_per_cell = UMI_per_cell,
    variation = variation,
    MOI = MOI,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    control_group = control_group,
    multiple_testing_alpha = multiple_testing_alpha,
    side = side,
    prop_non_null = prop_non_null,
    min_power_threshold = min_power_threshold,
    max_power_threshold = max_power_threshold
  )

  min_total_cells <- cell_range$min_cells
  max_total_cells <- cell_range$max_cells

  # Step 3: Create logarithmically-spaced grids
  # Logarithmic spacing for cells (total cells)
  cells_seq <- exp(seq(log(min_total_cells), log(max_total_cells), length.out = grid_size))

  # Logarithmic spacing for reads per cell
  reads_seq <- exp(seq(log(min_reads_per_cell), log(max_reads_per_cell), length.out = grid_size))

  # Step 4: Calculate treatment cell counts and control cell counts
  # Convert total cells to treatment cells using experimental design parameters
  total_gRNAs <- num_targets * gRNAs_per_target + non_targeting_gRNAs
  cells_per_target_seq <- (cells_seq * gRNAs_per_target * MOI) / total_gRNAs

  # Calculate control cells based on control group type
  if (control_group == "complement") {
    # Control cells = total cells - treatment cells
    num_cntrl_cells_seq <- cells_seq - cells_per_target_seq
  } else {
    # nt_cells: control cells = non-targeting gRNA cells
    num_cntrl_cells_seq <- (cells_seq * non_targeting_gRNAs * MOI) / total_gRNAs
  }

  # Step 5: Create all combinations and compute power
  design_grid <- expand.grid(
    cells_idx = 1:grid_size,
    reads_idx = 1:grid_size
  ) |>
    dplyr::mutate(
      cells_per_target = cells_per_target_seq[cells_idx],
      num_total_cells = cells_seq[cells_idx],
      reads_per_cell = reads_seq[reads_idx],
      num_cntrl_cells = num_cntrl_cells_seq[cells_idx]
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      library_size = fit_read_UMI_curve_cpp(
        reads_per_cell = reads_per_cell,
        UMI_per_cell = UMI_per_cell,
        variation = variation
      ),
      overall_power = compute_single_power_cpp(
        num_cells = num_total_cells,
        reads_per_cell = reads_per_cell,
        fc_expression_df = fc_expression_df,
        UMI_per_cell = UMI_per_cell,
        variation = variation,
        MOI = MOI,
        num_targets = num_targets,
        gRNAs_per_target = gRNAs_per_target,
        non_targeting_gRNAs = non_targeting_gRNAs,
        control_group = control_group,
        multiple_testing_alpha = multiple_testing_alpha,
        side = side,
        prop_non_null = prop_non_null
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(cells_per_target, num_total_cells, reads_per_cell, library_size, overall_power)

  return(design_grid)
}

#' Compute power analysis for full parameter grid
#'
#' @description
#' This function integrates compute_power_plan_per_grid() to create a comprehensive
#' power analysis across multiple parameter combinations (TPM thresholds, fold changes).
#'
#' @export
compute_power_plan_full_grid <- function(
    # power-determining parameters
    tpm_threshold, minimum_fold_change,
    # experimental parameters
    MOI = 10, num_targets = 100, non_targeting_gRNAs = 10, gRNAs_per_target = 4, gRNA_variability = 0.13,
    # analysis parameters
    control_group = "complement", side = "left", multiple_testing_alpha = 0.05, prop_non_null = 0.1,
    # data inputs
    baseline_expression_stats, library_info,
    # grid parameters
    grid_size = 10, min_power_threshold = 0.01, max_power_threshold = 0.8
){

  ####################### construct the tpm_threshold ##########################
  if(tpm_threshold == "varying"){
    tpm_threshold_list <- unname(round(quantile(baseline_expression_stats$relative_expression,
                                                probs = seq(0.1, 0.9, length.out = 5)) * 1e6))
  }else{
    tpm_threshold_list <- tpm_threshold
  }

  ####################### construct the minimum_fold_change ####################
  if(minimum_fold_change == "varying"){
    minimum_fold_change_list <- switch (side,
      both = {c(seq(0.5, 0.95, length.out = 5), seq(1.05, 1.5, length.out = 5))},
      left = {seq(0.5, 0.95, length.out = 10)},
      right = {seq(1.05, 1.5, length.out = 10)}
    )
  }else{
    minimum_fold_change_list <- minimum_fold_change
  }

  # construct parameter grid by expanding the grid
  parameter_grid <- expand.grid(
    minimum_fold_change = minimum_fold_change_list,
    tpm_threshold = tpm_threshold_list
  )

  ###################### Generate complete power analysis grid #####################
  parameter_grid <- parameter_grid |>
    dplyr::rowwise() |>
    dplyr::mutate(
      # Create fold change expression dataframe for this combination
      fc_expression_df = list({
        # Filter expression data by TPM threshold
        filtered_expression_df <- baseline_expression_stats |>
          dplyr::filter(relative_expression >= tpm_threshold / 1e6)

        # Add fold change parameters
        filtered_expression_df |>
          dplyr::rowwise() |>
          dplyr::mutate(
            grna_effects = list(stats::rnorm(n = gRNAs_per_target,
                                            mean = minimum_fold_change,
                                            sd = gRNA_variability) |> pmax(0)),
            avg_fold_change = mean(grna_effects),
            avg_fold_change_sq = mean(grna_effects^2)
          ) |>
          dplyr::select(-grna_effects) |>  # Remove intermediate column
          dplyr::ungroup()
      }),
      # Generate experimental design grid with power calculations
      power_grid = list(
        compute_power_plan_per_grid(
          fc_expression_df = fc_expression_df,
          library_info = library_info,
          grid_size = grid_size,
          min_power_threshold = min_power_threshold,
          max_power_threshold = max_power_threshold,
          MOI = MOI,
          num_targets = num_targets,
          gRNAs_per_target = gRNAs_per_target,
          non_targeting_gRNAs = non_targeting_gRNAs,
          control_group = control_group,
          multiple_testing_alpha = multiple_testing_alpha,
          side = side,
          prop_non_null = prop_non_null
        )
      )
    ) |>
    dplyr::ungroup()

  # Unnest the power grids to create a flat dataframe
  result <- parameter_grid |>
    dplyr::select(-fc_expression_df) |>  # Remove intermediate data
    tidyr::unnest(power_grid) |>
    dplyr::select(minimum_fold_change, tpm_threshold,
                  cells_per_target, num_total_cells, reads_per_cell,
                  library_size, overall_power)

  return(result)
}
