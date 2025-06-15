# Test for compute_power_plan_overall function
# Validates the main power analysis pipeline integration

library(testthat)

test_that("compute_power_plan_overall matches manual step-by-step calculation", {
  
  # Simple, controlled test case for precise verification
  fc_expression_df <- data.frame(
    fold_change = c(0.5, 1.0, 2.0),
    relative_expression = c(1e-5, 5e-5, 1e-4),
    expression_size = c(0.5, 1.0, 2.0)
  )
  
  # Fixed experimental parameters
  num_total_cells <- 5000
  library_size <- 10000
  MOI <- 10
  num_targets <- 20
  gRNAs_per_target <- 3
  non_targeting_gRNAs <- 15
  multiple_testing_alpha <- 0.05
  multiple_testing_method <- "BH"
  control_group <- "complement"
  side <- "both"
  prop_non_null <- 0.3
  
  # Get integrated result
  integrated_result <- compute_power_plan_overall(
    num_total_cells = num_total_cells,
    library_size = library_size,
    MOI = MOI,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    multiple_testing_alpha = multiple_testing_alpha,
    multiple_testing_method = multiple_testing_method,
    control_group = control_group,
    side = side,
    fc_expression_df = fc_expression_df,
    prop_non_null = prop_non_null,
    return_full_results = TRUE
  )
  
  # Verify output structure
  expect_type(integrated_result, "list")
  expect_true(all(c("overall_power", "sig_cutoff", "num_trt_cells", "num_cntrl_cells") %in% names(integrated_result)))
  expect_true(integrated_result$overall_power >= 0 && integrated_result$overall_power <= 1)
  expect_true(is.finite(integrated_result$overall_power))
  expect_true(integrated_result$sig_cutoff > 0 && integrated_result$sig_cutoff <= multiple_testing_alpha)
  expect_true(is.finite(integrated_result$sig_cutoff))
  
  # Manual calculation step 1: Cell count allocation
  # Replicate exact logic from the function
  manual_num_trt_cells <- gRNAs_per_target * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
  manual_num_cntrl_cells <- switch(control_group,
                                  complement = {
                                    num_total_cells - manual_num_trt_cells
                                  },
                                  nt_cells = {
                                    non_targeting_gRNAs * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
                                  })
  manual_num_cntrl_cells <- round(manual_num_cntrl_cells)
  
  # Compare cell count calculations
  expect_equal(integrated_result$num_trt_cells, manual_num_trt_cells, tolerance = 1e-12)
  expect_equal(integrated_result$num_cntrl_cells, manual_num_cntrl_cells, tolerance = 1e-12)
  
  # Manual calculation step 2: Monte Carlo test statistics
  manual_mc_results <- compute_monte_carlo_teststat_cpp(
    fc_expression_df = fc_expression_df,
    library_size = library_size,
    num_trt_cells = manual_num_trt_cells,
    num_cntrl_cells = manual_num_cntrl_cells
  )
  
  # Verify Monte Carlo results structure
  expect_equal(length(manual_mc_results$means), nrow(fc_expression_df))
  expect_equal(length(manual_mc_results$sds), nrow(fc_expression_df))
  expect_true(all(is.finite(manual_mc_results$means)))
  expect_true(all(is.finite(manual_mc_results$sds)))
  expect_true(all(manual_mc_results$sds > 0))
  
  # Manual calculation step 3: Multiple testing cutoff
  manual_cutoff <- switch(multiple_testing_method,
                         BH = {
                           compute_BH_plan(
                             mean_list = manual_mc_results$means,
                             sd_list = manual_mc_results$sds,
                             side = side,
                             multiple_testing_alpha = multiple_testing_alpha,
                             prop_non_null = prop_non_null
                           )
                         })
  
  # Compare cutoff calculations
  expect_equal(integrated_result$sig_cutoff, manual_cutoff, tolerance = 1e-12)
  
  # Manual calculation step 4: Overall power computation
  manual_powers <- rejection_computation_cpp(
    mean_list = manual_mc_results$means,
    sd_list = manual_mc_results$sds,
    side = side,
    cutoff = manual_cutoff
  )
  manual_overall_power <- mean(manual_powers)
  
  # Compare final power calculation
  expect_equal(integrated_result$overall_power, manual_overall_power, tolerance = 1e-12)
  
  # Test simple return mode consistency
  simple_result <- compute_power_plan_overall(
    num_total_cells = num_total_cells,
    library_size = library_size,
    MOI = MOI,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    multiple_testing_alpha = multiple_testing_alpha,
    multiple_testing_method = multiple_testing_method,
    control_group = control_group,
    side = side,
    fc_expression_df = fc_expression_df,
    prop_non_null = prop_non_null,
    return_full_results = FALSE
  )
  
  expect_type(simple_result, "double")
  expect_equal(simple_result, integrated_result$overall_power, tolerance = 1e-12)
  
  # Test different control group method
  integrated_result_nt <- compute_power_plan_overall(
    num_total_cells = num_total_cells,
    library_size = library_size,
    MOI = MOI,
    num_targets = num_targets,
    gRNAs_per_target = gRNAs_per_target,
    non_targeting_gRNAs = non_targeting_gRNAs,
    multiple_testing_alpha = multiple_testing_alpha,
    multiple_testing_method = multiple_testing_method,
    control_group = "nt_cells",
    side = side,
    fc_expression_df = fc_expression_df,
    prop_non_null = prop_non_null,
    return_full_results = TRUE
  )
  
  # Manual calculation for nt_cells control group
  manual_num_cntrl_cells_nt <- non_targeting_gRNAs * num_total_cells * MOI / (num_targets * gRNAs_per_target + non_targeting_gRNAs)
  manual_num_cntrl_cells_nt <- round(manual_num_cntrl_cells_nt)
  
  expect_equal(integrated_result_nt$num_cntrl_cells, manual_num_cntrl_cells_nt, tolerance = 1e-12)
  expect_equal(integrated_result_nt$num_trt_cells, manual_num_trt_cells, tolerance = 1e-12)  # Same treatment cells
  
  # Verify nt_cells gives different result from complement
  expect_false(isTRUE(all.equal(integrated_result$overall_power, integrated_result_nt$overall_power)))
  
  # Test different test sides
  for (test_side in c("left", "right", "both")) {
    side_result <- compute_power_plan_overall(
      num_total_cells = num_total_cells,
      library_size = library_size,
      fc_expression_df = fc_expression_df,
      side = test_side,
      return_full_results = FALSE
    )
    
    expect_true(is.finite(side_result))
    expect_true(side_result >= 0 && side_result <= 1)
  }
  
  # Print summary for manual inspection
  cat("\n--- Power Plan Overall Test Summary ---")
  cat("\nIntegrated overall power:", round(integrated_result$overall_power, 4))
  cat("\nManual overall power:", round(manual_overall_power, 4))
  cat("\nSignificance cutoff:", round(integrated_result$sig_cutoff, 6))
  cat("\nTreatment cells:", round(integrated_result$num_trt_cells, 1))
  cat("\nControl cells:", integrated_result$num_cntrl_cells)
  cat("\nMonte Carlo samples:", length(manual_mc_results$means))
  cat("\nAll manual calculations match integrated results exactly")
})