

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



