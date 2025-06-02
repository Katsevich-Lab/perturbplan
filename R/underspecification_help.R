
compute_underspecified_power <- function(
  # experimental information
  num_total_cells, library_size, MOI = 10, num_targets = 100, gRNAs_per_target = 4, non_targeting_gRNAs = 10,
  # analysis information
  multiple_testing_alpha = 0.05, multiple_testing_method = "BH", control_group = "complement", side = "left", num_pairs = 1000,
  # effect size and baseline expression information
  fc_expression_df, expression_dispersion_curve, prop_non_null = 0.1){

  ################ compute the treatment and control cells #####################
  num_trt_cells <- gRNAs_per_target * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
  num_cntrl_cells <- round(switch(control_group,
                                  complement = {
                                    num_total_cells - num_trt_cells
                                  },
                                  nt_cells = {
                                    non_targeting_gRNAs * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
                                  }))

  ##############################################################################
  ########################## compute overall power #############################

  ############## compute the sufficient statistics dataframe ###################
  power_df <- compute_power_function(fc_expression_df = fc_expression_df, library_size = library_size,
                                     num_trt_cells = num_trt_cells, num_cntrl_cells = num_cntrl_cells,
                                     side = side)

  ########################## compute the cutoff ################################
  sig_cutoff <- switch(multiple_testing_method,
                       BH = {
                         compute_underspecified_BH_cutoff(power_f = power_df$marginal_power_function,
                                                          multiple_testing_alpha = multiple_testing_alpha,
                                                          prop_non_null = prop_non_null, num_pairs = num_pairs)
                       })

  ####################### compute the overall power ############################
  overall_power <- power_df$marginal_power_function(cutoff = sig_cutoff)

  ##############################################################################
  #################### compute intra-experiment curves #########################

  ##################### compute the power on fold change #######################
  power_on_fc <- function(fold_change){
    power_df$power_on_fc(cutoff = sig_cutoff, fold_change = fold_change)
  }

  ################## compute the power on relative expression ##################
  power_on_pi <- function(TPM){

    # extract relative expression and dispersion
    relative_expression <- TPM / 1e6
    expression_size <- expression_dispersion_curve(relative_expression)

    # return the power curve
    return(mean(sapply(fc_expression_df$fold_change, function(fold_change){
      power_df$conditional_power_function(relative_expression = relative_expression,
                                          expression_size = expression_size,
                                          cutoff = sig_cutoff,
                                          fold_change = fold_change)
    })))
  }

  # return the output
  return(list(
    overall_power = overall_power,
    power_on_fc = power_on_fc,
    power_on_pi = power_on_pi
  ))
}

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

compute_power_function <- function(fc_expression_df, library_size, num_trt_cells, num_cntrl_cells, side){

  # compute the conditional power function
  conditional_power_function <- function(relative_expression, expression_size, cutoff, fold_change, side_use = side){

    # extract mean and sd of test statistic
    mean_sd_list <- compute_distribution_teststat(num_trt_cells = num_trt_cells, num_cntrl_cells = num_cntrl_cells,
                                                  num_trt_cells_sq = num_trt_cells^2, expression_mean = library_size * relative_expression,
                                                  expression_size = expression_size, fold_change_mean = fold_change, fold_change_sd = 0)

    # compute the rejection rate
    rejection_computation(mean_list = mean_sd_list[[1]]["mean"],
                          sd_list = mean_sd_list[[1]]["sd"],
                          side = side_use,
                          cutoff = cutoff)
  }

  # compute the marginal power function
  marginal_power_function <- function(cutoff){

    # obtain conditional_power_df
    conditional_power_df <- fc_expression_df |>
      dplyr::rowwise() |>
      dplyr::mutate(conditional_power = conditional_power_function(
        relative_expression = relative_expression,
        expression_size = expression_size,
        cutoff = cutoff,
        fold_change = fold_change,
        side_use = side
      ))
    return(mean(conditional_power_df$conditional_power))
  }

  # compute the power against fold change
  power_on_fc <- function(cutoff, fold_change){

    # obtain conditional_power_df
    conditional_power_df <- fc_expression_df |>
      dplyr::select(relative_expression, expression_size) |>
      dplyr::rowwise() |>
      dplyr::mutate(conditional_power = conditional_power_function(
        relative_expression = relative_expression,
        expression_size = expression_size,
        cutoff = cutoff,
        fold_change = fold_change,
        side_use = side
      ))
    return(mean(conditional_power_df$conditional_power))
  }

  # return the function
  return(list(
    conditional_power_function = conditional_power_function,
    marginal_power_function = marginal_power_function,
    power_on_fc = power_on_fc
  ))
}

compute_underspecified_BH_cutoff <- function(power_f,
                                             multiple_testing_alpha,
                                             prop_non_null,
                                             num_pairs){

  ##################### Bisection to find cutoff ###############################
  # Helper: FDP(t) – multiple_testing_alpha  (monotone in t)
  f <- function(t) {
    t / (1 - prop_non_null + prop_non_null * power_f(cutoff = t)) - multiple_testing_alpha
  }

  # define lower and upper end
  lower <- 1 / num_pairs
  upper <- multiple_testing_alpha
  f_low  <- f(lower)
  f_high <- f(upper)

  # compute the function in starting points
  tol <- 1 / num_pairs
  if (f_low >= 0){
    return(lower)
  }else{
    iter <- 0
    repeat {

      # update iter and compute the mean of lower and upper
      iter <- iter + 1
      mid  <- (lower + upper) / 2
      fmid <- f(mid)

      # tell if the convergence is met
      if (fmid < 0 & abs(fmid) < tol)
        return(mid)

      # Shrink bracket, preserving sign change
      if (fmid > 0) {
        upper <- mid
      } else {
        lower <- mid
      }
    }
  }
}

