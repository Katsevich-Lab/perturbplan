# This is a scratch script
# library size computation
extract_library_parameters <- function(gene_list = "full", cell_type = "K562", assay_type = "Perturb_seq"){

  # load pilot data
  pilot_data <- perturbplan:::get_pilot_data_from_package(cell_type)

  # extract baseline expression
  baseline_expression <- pilot_data$baseline_expression$baseline_expression

  # extract library parameters
  library_parameters <- pilot_data$library_parameters

  # define the UMI_per_cell
  UMI_per_cell <- switch (assay_type,
                          Perturb_seq = {
                            library_parameters$UMI_per_cell
                          },
                          TAP_seq = {
                            # extract the baseline expression for targeted genes
                            library_parameters$UMI_per_cell * sum(baseline_expression[gene_list, "relative_expression"])
                          }
  )

  # return the library parameters
  return(list(
    UMI_per_cell = UMI_per_cell, variation = library_parameters$variation
  ))
}

# power calculation
power_calculation <- function(
    # Experimental information
  cell_type = "K562", assay_type = "Perturb_seq", gene_list = "full",
  # Sequencing information
  reads_per_cell = 1e4, target_mapping_efficiency = 1,
  # Target information
  gRNAs_per_target = 4, cells_per_target = 500, minimum_fold_change = 0.9, gRNA_variability = 0.13,
  # Analysis information
  discovery_pairs = NULL, prop_non_null = 0.1, multiple_testing_alpha = 0.1, tpm_threshold = 10, side = "left",
  # Pilot data for TAP-seq
  custom_pilot_data = NULL, library_parameters = NULL){

  # extract the whole Tx data
  whole_tx_pilot <- perturbplan:::get_pilot_data_from_package(cell_type)

  ######################## construct fc_expression_info ########################
  # extract the TPM, FC and gRNAs per target grid
  parameter_df <- suppressWarnings(
    expand.grid(
      minimum_fold_change = minimum_fold_change,
      tpm_threshold = tpm_threshold,
      gRNAs_per_target = gRNAs_per_target
    ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        subsetted_genes = list(
          whole_tx_pilot$baseline_expression$baseline_expression |>
            dplyr::filter(relative_expression >= tpm_threshold / 1e6) |>
            dplyr::select(response_id) |>
            dplyr::pull()
        ),
        fc_expression_df = list(
          extract_fc_expression_info(
            minimum_fold_change = minimum_fold_change,
            gRNA_variability = gRNA_variability,
            biological_system = cell_type,
            B = 1000,
            gene_list = discovery_pairs$response_id[discovery_pairs$response_id %in% subsetted_genes],
            tpm_threshold = 0,
            gRNAs_per_target = gRNAs_per_target,
            custom_pilot_data = custom_pilot_data
          )$fc_expression_df
        )
      )
  )

  ################### factorize in cells and reads #############################
  # create cells_reads_df
  parameter_df <- parameter_df |>
    tidyr::crossing(
      expand.grid(
        reads_per_cell = reads_per_cell,
        num_trt_cells = cells_per_target
      )
    ) |>
    dplyr::mutate(num_cntrl_cells = 1e5 - num_trt_cells)

  ######################### compute library size ###############################
  # extract library parameters
  if(is.null(library_parameters)){
    library_parameters <- extract_library_parameters(gene_list = gene_list, cell_type = cell_type,
                                                     assay_type = assay_type)
  }

  parameter_df <- parameter_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      library_size = fit_read_UMI_curve_cpp(reads_per_cell = reads_per_cell * target_mapping_efficiency,
                                            UMI_per_cell = library_parameters$UMI_per_cell,
                                            variation = library_parameters$variation)
    ) |>
    dplyr::ungroup()

  ########################### compute overall power ############################
  parameter_df <- parameter_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      overall_power = compute_power_plan_overall(
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        library_size = library_size,
        multiple_testing_alpha = multiple_testing_alpha,
        side = side,
        fc_expression_df = fc_expression_df,
        prop_non_null = prop_non_null
      )
    ) |>
    dplyr::ungroup()

  # return the parameter_grid output
  return(parameter_df |>
           dplyr::rename(cells_per_target = num_trt_cells) |>
           dplyr::select(cells_per_target, reads_per_cell,
                         minimum_fold_change, tpm_threshold, gRNAs_per_target,
                         library_size, overall_power))
}
