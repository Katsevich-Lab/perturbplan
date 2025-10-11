# Test input validation for find_optimal_cost_design function

library(testthat)

# Create test data setup function
setup_validation_test_data <- function() {
  # Valid cost-power data frame
  cost_power_df <- data.frame(
    TPM_threshold = rep(c(5, 10, 15), each = 3),
    overall_power = c(0.7, 0.8, 0.9, 0.75, 0.85, 0.95, 0.6, 0.7, 0.8),
    total_cost = c(1000, 1500, 2000, 1100, 1600, 2100, 900, 1400, 1900),
    cells_per_target = rep(c(500, 1000, 1500), 3),
    sequenced_reads_per_cell = rep(c(8000, 10000, 12000), 3),
    minimum_fold_change = rep(0.8, 9)
  )
  
  list(cost_power_df = cost_power_df)
}

test_that("input_check_find_optimal_cost_design validates cost_power_df", {
  test_data <- setup_validation_test_data()
  
  # Test missing cost_power_df
  expect_error(
    input_check_find_optimal_cost_design(
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "cost_power_df.*must be specified"
  )
  
  # Test non-data.frame cost_power_df
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = list(a = 1, b = 2),
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "cost_power_df.*must be a data frame"
  )
  
  # Test empty cost_power_df
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = data.frame(),
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "cost_power_df.*cannot be empty"
  )
  
  # Test missing required columns
  incomplete_df <- test_data$cost_power_df[, !names(test_data$cost_power_df) %in% "overall_power"]
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = incomplete_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "cost_power_df.*missing required columns.*overall_power"
  )
})

test_that("input_check_find_optimal_cost_design validates minimizing_variable", {
  test_data <- setup_validation_test_data()
  
  # Test non-character minimizing_variable
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = 123,
      power_target = 0.8,
      power_precision = 0.05
    ),
    "minimizing_variable.*must be a single character string"
  )
  
  # Test invalid minimizing_variable
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "invalid_var",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "minimizing_variable.*must be one of.*TPM_threshold.*minimum_fold_change"
  )
  
  # Test minimizing_variable not in data
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df[, !names(test_data$cost_power_df) %in% "TPM_threshold"],
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "cost_power_df.*must contain a column named.*TPM_threshold"
  )
})

test_that("input_check_find_optimal_cost_design validates power parameters", {
  test_data <- setup_validation_test_data()
  
  # Test missing power_target
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_precision = 0.05
    ),
    "power_target.*must be specified"
  )
  
  # Test invalid power_target
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 1.5,
      power_precision = 0.05
    ),
    "power_target.*must be between 0 and 1"
  )
  
  # Test invalid power_precision
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = -0.1
    ),
    "power_precision.*must be between 0 and 1"
  )
  
  # Test power_precision >= power_target
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.9
    ),
    "power_precision.*should be smaller than.*power_target"
  )
})

test_that("input_check_find_optimal_cost_design validates experimental parameters", {
  test_data <- setup_validation_test_data()
  
  # Test invalid MOI
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05,
      MOI = -5
    ),
    "MOI.*must be a positive numeric value"
  )
  
  # Test invalid num_targets
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05,
      num_targets = 50.5
    ),
    "num_targets.*must be a positive integer"
  )
})

test_that("input_check_find_optimal_cost_design validates cost parameters", {
  test_data <- setup_validation_test_data()
  
  # Test invalid cost_per_captured_cell
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05,
      cost_per_captured_cell = -0.1
    ),
    "cost_per_captured_cell.*must be a non-negative numeric value"
  )
  
  # Test invalid cost_per_million_reads
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05,
      cost_per_million_reads = -0.5
    ),
    "cost_per_million_reads.*must be a non-negative numeric value"
  )
  
  # Test invalid cost_grid_size
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05,
      cost_grid_size = 0
    ),
    "cost_grid_size.*must be a positive integer"
  )
  
  # Test warning for large cost_grid_size
  expect_warning(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05,
      cost_grid_size = 1500
    ),
    "cost_grid_size.*very large.*may result in long computation times"
  )
})

test_that("input_check_find_optimal_cost_design validates data columns", {
  test_data <- setup_validation_test_data()
  
  # Test non-numeric overall_power
  invalid_df <- test_data$cost_power_df
  invalid_df$overall_power <- as.character(invalid_df$overall_power)
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = invalid_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "Column.*overall_power.*must be numeric"
  )
  
  # Test missing values in required columns
  invalid_df <- test_data$cost_power_df
  invalid_df$total_cost[1] <- NA
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = invalid_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "Column.*total_cost.*cannot contain missing values"
  )
  
  # Test non-positive values in cost columns
  invalid_df <- test_data$cost_power_df
  invalid_df$cells_per_target[1] <- 0
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = invalid_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "Column.*cells_per_target.*must contain positive values"
  )
  
  # Test overall_power out of bounds
  invalid_df <- test_data$cost_power_df
  invalid_df$overall_power[1] <- 1.5
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = invalid_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    ),
    "Column.*overall_power.*must contain values between 0 and 1"
  )
})

test_that("input_check_find_optimal_cost_design validates power criteria feasibility", {
  # Create data where no designs can meet power criteria
  impossible_df <- data.frame(
    TPM_threshold = c(5, 10, 15),
    overall_power = c(0.1, 0.2, 0.3),  # All too low
    total_cost = c(1000, 1500, 2000),
    cells_per_target = c(500, 1000, 1500),
    sequenced_reads_per_cell = c(8000, 10000, 12000),
    minimum_fold_change = rep(0.8, 3)
  )
  
  expect_error(
    input_check_find_optimal_cost_design(
      cost_power_df = impossible_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,  # Target 80% but data only has 10-30%
      power_precision = 0.05
    ),
    "No designs.*can meet the power criteria.*Power range.*Target range"
  )
})

test_that("input_check_find_optimal_cost_design passes with valid inputs", {
  test_data <- setup_validation_test_data()
  
  # Should not throw any errors with valid inputs
  expect_silent(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.8,
      power_precision = 0.05
    )
  )
  
  # Test with minimum_fold_change as minimizing variable
  expect_silent(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "minimum_fold_change",
      power_target = 0.7,
      power_precision = 0.1
    )
  )
  
  # Test with different parameter values
  expect_silent(
    input_check_find_optimal_cost_design(
      cost_power_df = test_data$cost_power_df,
      minimizing_variable = "TPM_threshold",
      power_target = 0.85,
      power_precision = 0.02,
      MOI = 15,
      num_targets = 50,
      cost_per_captured_cell = 0.1,
      cost_per_million_reads = 0.5,
      cost_grid_size = 100
    )
  )
})