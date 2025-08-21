# Test for obtain_fixed_variable_constraining_cost function
# Validates cost-constrained experimental design parameter optimization

library(testthat)

# Test basic functionality: optimize reads_per_cell given fixed cells_per_target
test_that("obtain_fixed_variable_constraining_cost optimizes reads_per_cell correctly", {
  result <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 1000,
    cells_per_target = 100,
    reads_per_cell = NULL
  )
  
  expect_type(result, "list")
  expect_named(result, c("cells_per_target", "reads_per_cell"))
  expect_equal(result$cells_per_target, 100)
  expect_true(result$reads_per_cell > 0)
  expect_true(result$reads_per_cell > 10)  # Should meet minimum threshold
})

# Test basic functionality: optimize cells_per_target given fixed reads_per_cell
test_that("obtain_fixed_variable_constraining_cost optimizes cells_per_target correctly", {
  result <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 1000,
    cells_per_target = NULL,
    reads_per_cell = 5000
  )
  
  expect_type(result, "list")
  expect_named(result, c("cells_per_target", "reads_per_cell"))
  expect_equal(result$reads_per_cell, 5000)
  expect_true(result$cells_per_target > 0)
  expect_true(result$cells_per_target > 10)  # Should be meaningful
})

# Test input validation: missing cost_constraint
test_that("obtain_fixed_variable_constraining_cost validates cost_constraint", {
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cells_per_target = 100,
      reads_per_cell = NULL
    ),
    "`cost_constraint` must be a positive numeric value!"
  )
  
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = -100,
      cells_per_target = 100,
      reads_per_cell = NULL
    ),
    "`cost_constraint` must be a positive numeric value!"
  )
  
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = 0,
      cells_per_target = 100,
      reads_per_cell = NULL
    ),
    "`cost_constraint` must be a positive numeric value!"
  )
})

# Test input validation: exactly one parameter must be NULL
test_that("obtain_fixed_variable_constraining_cost validates parameter combination", {
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = 1000,
      cells_per_target = NULL,
      reads_per_cell = NULL
    ),
    "Exactly one of `reads_per_cell` or `cells_per_target` must be NULL for optimization!"
  )
  
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = 1000,
      cells_per_target = 100,
      reads_per_cell = 5000
    ),
    "Exactly one of `reads_per_cell` or `cells_per_target` must be NULL for optimization!"
  )
})

# Test input validation: positive parameter values
test_that("obtain_fixed_variable_constraining_cost validates parameter values", {
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = 1000,
      cells_per_target = -100,
      reads_per_cell = NULL
    ),
    "`cells_per_target` must be positive when provided!"
  )
  
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = 1000,
      cells_per_target = NULL,
      reads_per_cell = -5000
    ),
    "`reads_per_cell` must be positive when provided!"
  )
  
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = 1000,
      cells_per_target = 100,
      reads_per_cell = NULL,
      cost_per_captured_cell = -0.086
    ),
    "`cost_per_captured_cell` must be positive!"
  )
  
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = 1000,
      cells_per_target = 100,
      reads_per_cell = NULL,
      mapping_efficiency = 1.5
    ),
    "`mapping_efficiency` must be between 0 and 1!"
  )
})

# Test cost constraint too tight scenarios
test_that("obtain_fixed_variable_constraining_cost handles tight cost constraints", {
  # Scenario 1: Cost too tight for reads per cell
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = 10,  # Very low budget
      cells_per_target = 1000,  # High cell count
      reads_per_cell = NULL
    ),
    "Cost constraint is too tight to get any reads per cell!"
  )
  
  # Scenario 2: Cost too tight for captured cells
  expect_error(
    obtain_fixed_variable_constraining_cost(
      cost_constraint = 0.5,  # Very low budget
      cells_per_target = NULL,
      reads_per_cell = 50000  # Very high read depth
    ),
    "Cost constraint is too tight to get more than 10 captured cells!"
  )
})

# Test different experimental parameters
test_that("obtain_fixed_variable_constraining_cost works with different experimental parameters", {
  result1 <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 2000,
    cells_per_target = 200,
    reads_per_cell = NULL,
    MOI = 5,
    num_targets = 50,
    gRNAs_per_target = 3
  )
  
  result2 <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 2000,
    cells_per_target = 200,
    reads_per_cell = NULL,
    MOI = 15,
    num_targets = 200,
    gRNAs_per_target = 6
  )
  
  expect_type(result1, "list")
  expect_type(result2, "list")
  expect_equal(result1$cells_per_target, 200)
  expect_equal(result2$cells_per_target, 200)
  # Different experimental parameters should give different reads per cell
  expect_false(result1$reads_per_cell == result2$reads_per_cell)
})

# Test cost model calculations are consistent
test_that("obtain_fixed_variable_constraining_cost cost calculations are internally consistent", {
  cost_budget <- 1500
  cells_target <- 150
  
  # Get optimized reads per cell
  result <- obtain_fixed_variable_constraining_cost(
    cost_constraint = cost_budget,
    cells_per_target = cells_target,
    reads_per_cell = NULL,
    cost_per_captured_cell = 0.1,
    cost_per_million_reads = 0.4,
    mapping_efficiency = 0.8
  )
  
  # Calculate expected costs manually
  num_captured_cells <- ((4 * 100 + 10) * cells_target / 4) / 10
  cell_cost <- 0.1 * num_captured_cells
  total_reads <- num_captured_cells * result$reads_per_cell / 0.8
  read_cost <- (total_reads / 1e6) * 0.4
  total_cost <- cell_cost + read_cost
  
  # Total cost should be close to budget (within small tolerance for rounding)
  expect_true(abs(total_cost - cost_budget) < 1)
})

# Test different cost parameters
test_that("obtain_fixed_variable_constraining_cost responds correctly to cost parameters", {
  base_params <- list(
    cost_constraint = 1000,
    cells_per_target = 100,
    reads_per_cell = NULL
  )
  
  # Higher cost per cell should reduce reads per cell
  result_low_cell_cost <- do.call(obtain_fixed_variable_constraining_cost, 
                                  c(base_params, list(cost_per_captured_cell = 0.05)))
  result_high_cell_cost <- do.call(obtain_fixed_variable_constraining_cost, 
                                   c(base_params, list(cost_per_captured_cell = 0.15)))
  
  expect_true(result_low_cell_cost$reads_per_cell > result_high_cell_cost$reads_per_cell)
  
  # Higher cost per read should reduce reads per cell
  result_low_read_cost <- do.call(obtain_fixed_variable_constraining_cost, 
                                  c(base_params, list(cost_per_million_reads = 0.2)))
  result_high_read_cost <- do.call(obtain_fixed_variable_constraining_cost, 
                                   c(base_params, list(cost_per_million_reads = 0.6)))
  
  expect_true(result_low_read_cost$reads_per_cell > result_high_read_cost$reads_per_cell)
})

# Test mapping efficiency effects
test_that("obtain_fixed_variable_constraining_cost handles mapping efficiency correctly", {
  base_params <- list(
    cost_constraint = 1000,
    cells_per_target = 100,
    reads_per_cell = NULL
  )
  
  # Higher mapping efficiency should allow more reads per cell (more efficient)
  result_low_eff <- do.call(obtain_fixed_variable_constraining_cost, 
                            c(base_params, list(mapping_efficiency = 0.5)))
  result_high_eff <- do.call(obtain_fixed_variable_constraining_cost, 
                             c(base_params, list(mapping_efficiency = 0.9)))
  
  expect_true(result_high_eff$reads_per_cell > result_low_eff$reads_per_cell)
})

# Test symmetry: optimizing in both directions should be consistent
test_that("obtain_fixed_variable_constraining_cost optimization is consistent", {
  # Start with a budget and optimize reads per cell
  result1 <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 1000,
    cells_per_target = 150,
    reads_per_cell = NULL,
    cost_per_captured_cell = 0.08,
    cost_per_million_reads = 0.35
  )
  
  # Now use those reads per cell to optimize cells per target with same budget
  result2 <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 1000,
    cells_per_target = NULL,
    reads_per_cell = result1$reads_per_cell,
    cost_per_captured_cell = 0.08,
    cost_per_million_reads = 0.35
  )
  
  # Should get back approximately the same cells per target
  expect_true(abs(result2$cells_per_target - 150) < 5)  # Allow small tolerance for rounding
})

# Test edge case: minimum thresholds (additional edge cases)
test_that("obtain_fixed_variable_constraining_cost handles extreme parameter combinations", {
  # Test very low reads per cell threshold
  result <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 100,
    cells_per_target = 10,  # Small experiment
    reads_per_cell = NULL
  )
  
  expect_true(result$reads_per_cell >= 10)
  
  # Test very low cells per target scenario
  result2 <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 100,
    cells_per_target = NULL,
    reads_per_cell = 1000  # Moderate read depth
  )
  
  expect_true(result2$cells_per_target > 0)
})

# Test return value structure
test_that("obtain_fixed_variable_constraining_cost returns correct structure", {
  result <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 1000,
    cells_per_target = 100,
    reads_per_cell = NULL
  )
  
  expect_type(result, "list")
  expect_length(result, 2)
  expect_named(result, c("cells_per_target", "reads_per_cell"))
  expect_type(result$cells_per_target, "double")
  expect_type(result$reads_per_cell, "double")
  expect_true(is.finite(result$cells_per_target))
  expect_true(is.finite(result$reads_per_cell))
})

# Test with realistic experimental scenarios
test_that("obtain_fixed_variable_constraining_cost works with realistic scenarios", {
  # Realistic scenario 1: Small-scale experiment
  small_exp <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 500,
    cells_per_target = 50,
    reads_per_cell = NULL,
    num_targets = 20
  )
  
  expect_true(small_exp$reads_per_cell > 1000)
  expect_equal(small_exp$cells_per_target, 50)
  
  # Realistic scenario 2: Large-scale experiment  
  large_exp <- obtain_fixed_variable_constraining_cost(
    cost_constraint = 5000,
    cells_per_target = NULL,
    reads_per_cell = 3000,
    num_targets = 500
  )
  
  expect_true(large_exp$cells_per_target > 100)
  expect_equal(large_exp$reads_per_cell, 3000)
})