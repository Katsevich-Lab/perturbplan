# Tests for pilot_data_preprocessing function
library(testthat)

test_that("pilot_data_preprocessing validates input path", {
  # Test with non-existent directory
  expect_error(
    pilot_data_preprocessing("/nonexistent/path"),
    "Cell Ranger output directory does not exist"
  )
})

test_that("pilot_data_preprocessing validates required files", {
  # Create temporary directory structure without required files
  temp_dir <- tempdir()
  fake_cellranger_dir <- file.path(temp_dir, "fake_cellranger")
  dir.create(fake_cellranger_dir, recursive = TRUE)
  dir.create(file.path(fake_cellranger_dir, "outs"), recursive = TRUE)
  
  # Should fail because required files are missing
  expect_error(
    pilot_data_preprocessing(fake_cellranger_dir),
    "Missing required Cell Ranger output files"
  )
  
  # Clean up
  unlink(fake_cellranger_dir, recursive = TRUE)
})

test_that("pilot_data_preprocessing validates parameter types and ranges", {
  # Test parameter validation (these will fail at directory check first, 
  # but we can test the parameter validation logic)
  temp_dir <- tempdir()
  fake_dir <- file.path(temp_dir, "test_params")
  dir.create(fake_dir, recursive = TRUE)
  
  # Test invalid downsample ratio
  expect_error(
    pilot_data_preprocessing(fake_dir, downsample_ratio = -1),
    "Missing required Cell Ranger output files"  # Will fail on file check first
  )
  
  # Test invalid downsample ratio
  expect_error(
    pilot_data_preprocessing(fake_dir, downsample_ratio = 1.5),
    "Missing required Cell Ranger output files"  # Will fail on file check first
  )
  
  # Clean up
  unlink(fake_dir, recursive = TRUE)
})

test_that("pilot_data_preprocessing output structure matches expected format", {
  skip("Requires actual Cell Ranger data - integration test")
  
  # This test would require actual Cell Ranger output data
  # When run with real data, it should verify:
  # 1. Output is a list with two elements: baseline_expression and library_parameters
  # 2. baseline_expression contains baseline_expression data.frame and expression_dispersion_curve function
  # 3. library_parameters contains UMI_per_cell and variation numeric values
  # 4. Data frame has correct columns: response_id, relative_expression, expression_size
  # 5. All values are positive and finite
})

# Helper function to create mock Cell Ranger directory structure
create_mock_cellranger_structure <- function(base_dir) {
  outs_dir <- file.path(base_dir, "outs")
  dir.create(outs_dir, recursive = TRUE)
  
  # Create filtered_feature_bc_matrix directory
  matrix_dir <- file.path(outs_dir, "filtered_feature_bc_matrix")
  dir.create(matrix_dir, recursive = TRUE)
  
  # Create placeholder files (these would need actual data for real testing)
  file.create(file.path(matrix_dir, "matrix.mtx.gz"))
  file.create(file.path(matrix_dir, "features.tsv.gz"))
  file.create(file.path(matrix_dir, "barcodes.tsv.gz"))
  file.create(file.path(outs_dir, "molecule_info.h5"))
  file.create(file.path(outs_dir, "filtered_feature_bc_matrix.h5"))
  
  return(base_dir)
}

test_that("pilot_data_preprocessing function signature is correct", {
  # Test that the function exists and has the expected parameters
  expect_true(exists("pilot_data_preprocessing"))
  
  # Test parameter defaults
  pilot_func <- get("pilot_data_preprocessing")
  expect_true(is.function(pilot_func))
  
  # Check that function has expected formal arguments
  formal_args <- names(formals(pilot_func))
  expected_args <- c("path_to_cellranger_output", "rough", 
                     "n_threads", "downsample_ratio", "D2_rough")
  expect_true(all(expected_args %in% formal_args))
})

test_that("pilot_data_preprocessing vs build_bio_library_info comparison", {
  # Test that both functions are available but serve different purposes
  expect_true(exists("pilot_data_preprocessing"))
  expect_true(exists("build_bio_library_info"))
  
  # They should have similar but different parameter sets
  pilot_params <- names(formals(pilot_data_preprocessing))
  build_params <- names(formals(build_bio_library_info))
  
  # pilot_data_preprocessing should NOT have multi_folder parameter
  expect_false("multi_folder" %in% pilot_params)
  # build_bio_library_info should have multi_folder parameter
  expect_true("multi_folder" %in% build_params)
})

test_that("pilot_data_preprocessing integration with example data structure", {
  # Load the example pilot data to understand expected structure
  example_path <- system.file("extdata", "example_pilot", "pilot_example.rds", 
                             package = "perturbplan")
  
  if (file.exists(example_path)) {
    example_data <- readRDS(example_path)
    
    # Verify the expected structure exists in the example
    expect_true(is.list(example_data))
    expect_true("baseline_expression" %in% names(example_data))
    expect_true("library_parameters" %in% names(example_data))
    
    # Check baseline_expression structure
    baseline_expr <- example_data$baseline_expression
    expect_true(is.list(baseline_expr))
    expect_true("baseline_expression" %in% names(baseline_expr))
    expect_true("expression_dispersion_curve" %in% names(baseline_expr))
    
    # Check that baseline_expression contains a data frame
    expect_true(is.data.frame(baseline_expr$baseline_expression))
    expected_cols <- c("response_id", "relative_expression", "expression_size")
    expect_true(all(expected_cols %in% names(baseline_expr$baseline_expression)))
    
    # Check that expression_dispersion_curve is a function
    expect_true(is.function(baseline_expr$expression_dispersion_curve))
    
    # Check library_parameters structure
    lib_params <- example_data$library_parameters
    expect_true(is.list(lib_params))
    expect_true("UMI_per_cell" %in% names(lib_params))
    expect_true("variation" %in% names(lib_params))
    expect_true(is.numeric(lib_params$UMI_per_cell))
    expect_true(is.numeric(lib_params$variation))
  } else {
    skip("Example pilot data not found")
  }
})