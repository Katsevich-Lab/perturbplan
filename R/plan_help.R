

extract_fc_expression_info <- function(fold_change_mean, fold_change_sd, biological_system =  "K562", B = 200){

  # set the random seed
  set.seed(1)

  ############## combine expression and effect size information ################
  baseline_expression_stats <- extract_baseline_expression(biological_system = biological_system)
  fc_expression_df <- data.frame(
    fold_change = stats::rnorm(n = B, mean = fold_change_mean, sd = fold_change_sd)
  ) |>
    dplyr::bind_cols(baseline_expression_stats$baseline_expression |> dplyr::slice_sample(n = B))

  ################## extract the expression-dispersion curve ###################
  expression_dispersion_curve <- baseline_expression_stats$expression_dispersion_curve

  # return the data frame
  return(list(
    fc_expression_df = fc_expression_df,
    expression_dispersion_curve = expression_dispersion_curve
  ))
}

extract_baseline_expression <- function(biological_system = "K562"){

  # sample baseline expression based on biological system
  switch(biological_system,
         K562 = {

           # load the Gasperini baseline expression list
           rds_path <- system.file("extdata/baseline_expression", "Gasperini.rds", package = "perturbplan", mustWork = TRUE)
           baseline_expression_list <- readRDS(rds_path)

         })

  # return the data frame with the susbampled rows
  return(list(
    baseline_expression = baseline_expression_list$baseline_expression,
    expression_dispersion_curve = baseline_expression_list$expression_dispersion_curve
  ))
}

#' Compute power curve across fold change values
#'
#' @param fc_output_grid Vector of fold change values to evaluate
#' @param fc_expression_df Data frame with Monte Carlo samples (fold_change, relative_expression, expression_size)
#' @param library_size Library size (reads per cell)
#' @param num_trt_cells Number of treatment cells
#' @param num_cntrl_cells Number of control cells
#' @param side Test sidedness ("left", "right", "both")
#' @param cutoff Significance cutoff for power computation
#' @return Data frame with fold_change and power columns
compute_fc_curve <- function(fc_output_grid, fc_expression_df, library_size, 
                            num_trt_cells, num_cntrl_cells, side, cutoff) {
  
  # Pre-compute expression means for efficiency
  mc_expression_means <- library_size * fc_expression_df$relative_expression
  
  power_by_fc_list <- vector("list", length(fc_output_grid))
  for (i in 1:length(fc_output_grid)) {
    fc_val <- fc_output_grid[i]
    
    # Compute power for this FC across all Monte Carlo expression samples
    fc_test_stats <- vector("list", nrow(fc_expression_df))
    for (j in 1:nrow(fc_expression_df)) {
      # Use C++ function for test statistic computation
      fc_test_stats[[j]] <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = fc_val,  # Use systematic FC grid point
        expression_mean = mc_expression_means[j],  # Use MC expression samples
        expression_size = fc_expression_df$expression_size[j],
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = num_trt_cells  # For single gRNA case
      )
    }
    
    fc_means <- sapply(fc_test_stats, function(x) x[["mean"]])
    fc_sds <- sapply(fc_test_stats, function(x) x[["sd"]])
    
    fc_powers <- rejection_computation_cpp(mean_list = fc_means,
                                          sd_list = fc_sds,
                                          side = side,
                                          cutoff = cutoff)
    power_by_fc_list[[i]] <- mean(fc_powers)
  }
  
  return(data.frame(
    fold_change = fc_output_grid,
    power = unlist(power_by_fc_list)
  ))
}

#' Compute power curve across expression values
#'
#' @param expr_output_grid Vector of relative expression values to evaluate
#' @param fc_expression_df Data frame with Monte Carlo samples (fold_change, relative_expression, expression_size)
#' @param library_size Library size (reads per cell)
#' @param expression_dispersion_curve Function mapping relative expression to size parameter
#' @param num_trt_cells Number of treatment cells
#' @param num_cntrl_cells Number of control cells
#' @param side Test sidedness ("left", "right", "both")
#' @param cutoff Significance cutoff for power computation
#' @return Data frame with relative_expression and power columns
compute_expression_curve <- function(expr_output_grid, fc_expression_df, library_size, 
                                    expression_dispersion_curve, num_trt_cells, 
                                    num_cntrl_cells, side, cutoff) {
  
  power_by_expr_list <- vector("list", length(expr_output_grid))
  for (i in 1:length(expr_output_grid)) {
    expr_val <- expr_output_grid[i]
    expr_mean <- library_size * expr_val
    
    # Use expression_dispersion_curve to get size parameter
    expr_size <- expression_dispersion_curve(expr_val)
    
    # Compute power for this expression across all Monte Carlo FC samples
    expr_test_stats <- vector("list", nrow(fc_expression_df))
    for (j in 1:nrow(fc_expression_df)) {
      # Use C++ function for test statistic computation
      expr_test_stats[[j]] <- compute_distribution_teststat_fixed_es_cpp(
        fold_change = fc_expression_df$fold_change[j],  # Use MC FC samples
        expression_mean = expr_mean,  # Use systematic expression grid point
        expression_size = expr_size,
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        num_cells = num_trt_cells  # For single gRNA case
      )
    }
    
    expr_means <- sapply(expr_test_stats, function(x) x[["mean"]])
    expr_sds <- sapply(expr_test_stats, function(x) x[["sd"]])
    
    expr_powers <- rejection_computation_cpp(mean_list = expr_means,
                                            sd_list = expr_sds,
                                            side = side,
                                            cutoff = cutoff)
    power_by_expr_list[[i]] <- mean(expr_powers)
  }
  
  return(data.frame(
    relative_expression = expr_output_grid,
    power = unlist(power_by_expr_list)
  ))
}

