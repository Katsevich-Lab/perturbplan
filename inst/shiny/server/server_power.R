# Server module for power calculations and reactive values

create_power_server <- function(input, output, session) {
  # Reactive values
  planned <- reactiveVal(FALSE)
  output$need_plan <- reactive(!planned())
  outputOptions(output, "need_plan", suspendWhenHidden = FALSE)
  observeEvent(input$plan_btn, {
    # Validate that if "Other" is selected, custom pilot data is provided
    if ((input$biological_system == "Other" || input$experimental_platform == "Other") && 
        (is.null(custom_baseline()) || is.null(custom_library()))) {
      showNotification(
        "Please upload custom pilot data when 'Other' is selected for biological system or experimental platform.",
        type = "error", duration = 10
      )
      return()
    }
    planned(TRUE)
  })
  
  # Gene list file upload handling
  gene_list <- reactiveVal(NULL)
  target_list <- reactiveVal(NULL)  # Store targets from uploaded CSV
  
  observeEvent(input$gene_list_file, {
    req(input$gene_list_file, input$gene_list_mode == "custom")
    
    tryCatch({
      # Read the uploaded CSV file
      file_ext <- tools::file_ext(input$gene_list_file$name)
      
      if (file_ext == "csv") {
        # Read CSV with headers
        uploaded_data <- read.csv(input$gene_list_file$datapath, 
                                stringsAsFactors = FALSE)
        
        # Validate required columns
        if (!all(c("grna_target", "response_id") %in% colnames(uploaded_data))) {
          showNotification("CSV must contain 'grna_target' and 'response_id' columns", type = "error")
          gene_list(NULL)
          output$gene_list_uploaded <- reactive(FALSE)
          outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
          return()
        }
        
        # Extract gene list preserving duplicates
        genes <- as.character(uploaded_data$response_id)
        targets <- as.character(uploaded_data$grna_target)
        
        # Clean the gene list and targets
        genes <- trimws(genes)  # Remove whitespace
        genes <- genes[genes != ""]  # Remove empty lines
        targets <- trimws(targets)
        targets <- targets[targets != ""]
        # Note: NOT removing duplicates - preserve multiplicity
        
        # Store the gene list and targets
        gene_list(genes)
        target_list(targets)
        
        # Calculate statistics for status
        total_pairs <- length(genes)
        unique_genes_uploaded <- length(unique(genes))
        unique_targets <- length(unique(targets))
        
        # Analyze TPM filtering and gene availability
        tryCatch({
          # Get baseline expression data to check gene availability
          if (input$biological_system == "Other") {
            # Use custom baseline data if available
            if (!is.null(custom_baseline())) {
              baseline_expression_stats <- custom_baseline()
              baseline_df <- baseline_expression_stats$baseline_expression
            } else {
              # No custom data available, skip filtering analysis with a user-friendly message
              output$gene_list_status <- renderUI({
                HTML(sprintf("Loaded %d pairs (%d unique genes, %d unique targets)<br/><em style='color:orange;'>Gene filtering analysis unavailable: Please upload custom pilot data for 'Other' biological system</em>", 
                            total_pairs, unique_genes_uploaded, unique_targets))
              })
              output$gene_list_uploaded <- reactive(TRUE)
              outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
              return()  # Exit early to avoid further processing
            }
          } else {
            # Use built-in data (using unexported function)
            baseline_expression_stats <- perturbplan:::extract_baseline_expression(biological_system = input$biological_system)
            baseline_df <- baseline_expression_stats$baseline_expression
          }
          
          # Apply TPM threshold filtering
          tpm_threshold_relative <- input$tpm_threshold / 1e6
          if ("relative_expression" %in% colnames(baseline_df)) {
            filtered_baseline_df <- baseline_df |>
              dplyr::filter(relative_expression >= tpm_threshold_relative)
          } else {
            filtered_baseline_df <- baseline_df
          }
          
          # Check which genes are found in the database before and after filtering
          unique_genes_list <- unique(genes)
          
          # Genes found in original database
          genes_in_database <- unique_genes_list[unique_genes_list %in% baseline_df$response_id]
          genes_not_in_database <- setdiff(unique_genes_list, genes_in_database)
          
          # Genes found after TPM filtering
          genes_after_filtering <- unique_genes_list[unique_genes_list %in% filtered_baseline_df$response_id]
          genes_filtered_by_tpm <- setdiff(genes_in_database, genes_after_filtering)
          
          # Calculate proportions
          prop_not_in_database <- length(genes_not_in_database) / unique_genes_uploaded
          prop_filtered_by_tpm <- length(genes_filtered_by_tpm) / unique_genes_uploaded
          prop_available <- length(genes_after_filtering) / unique_genes_uploaded
          
          # Create comprehensive status message
          basic_msg <- sprintf("Loaded %d pairs (%d unique genes, %d unique targets)", 
                               total_pairs, unique_genes_uploaded, unique_targets)
          
          filtering_msg <- sprintf("Gene filtering: %.1f%% not in database, %.1f%% filtered by TPM threshold (≥%d), %.1f%% available for analysis",
                                  prop_not_in_database * 100, prop_filtered_by_tpm * 100, input$tpm_threshold, prop_available * 100)
          
          # Combine messages for inline display
          combined_msg <- paste(basic_msg, filtering_msg, sep = "<br/>")
          
          # Update status with combined message in UI
          output$gene_list_status <- renderUI({
            HTML(combined_msg)
          })
          
        }, error = function(e) {
          # Show error for debugging and fallback to basic status
          output$gene_list_status <- renderUI({
            HTML(sprintf("Loaded %d pairs (%d unique genes, %d unique targets)<br/><em style='color:red;'>Filtering analysis failed: %s</em>", 
                        total_pairs, unique_genes_uploaded, unique_targets, e$message))
          })
        })
        
        output$gene_list_uploaded <- reactive(TRUE)
        outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
        
        # Check for target count mismatch and show warning if needed
        if (!is.null(input$num_targets) && unique_targets != input$num_targets) {
          showNotification(
            paste0("Warning: CSV file contains ", unique_targets, " unique targets, ",
                   "but experimental choice specifies ", input$num_targets, " targets. ",
                   "The analysis will use the ", unique_targets, " targets from the CSV file."),
            type = "warning",
            duration = 8  # Show warning for 8 seconds
          )
        }
        
      } else {
        showNotification("Please upload a CSV file with grna_target and response_id columns", type = "error")
        gene_list(NULL)
        target_list(NULL)
        output$gene_list_uploaded <- reactive(FALSE)
        outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
      }
      
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      gene_list(NULL)
      target_list(NULL)
      output$gene_list_uploaded <- reactive(FALSE)
      outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
    })
  })
  
  # Reset gene list when file is removed or mode changes to random
  observe({
    if (is.null(input$gene_list_file) || input$gene_list_mode == "random") {
      gene_list(NULL)
      target_list(NULL)
      output$gene_list_uploaded <- reactive(FALSE)
      outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
    }
  })
  
  # Reset gene list when switching to random mode
  observeEvent(input$gene_list_mode, {
    if (input$gene_list_mode == "random") {
      gene_list(NULL)
      target_list(NULL)
      output$gene_list_uploaded <- reactive(FALSE)
      outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
    }
  })
  
  # Check for target count mismatch when num_targets changes (for custom mode)
  observeEvent(input$num_targets, {
    # Only check if in custom mode and a file is uploaded
    if (input$gene_list_mode == "custom" && !is.null(target_list()) && length(target_list()) > 0) {
      unique_targets_uploaded <- length(unique(target_list()))
      
      if (unique_targets_uploaded != input$num_targets) {
        showNotification(
          paste0("Warning: CSV file contains ", unique_targets_uploaded, " unique targets, ",
                 "but experimental choice specifies ", input$num_targets, " targets. ",
                 "The analysis will use the ", unique_targets_uploaded, " targets from the CSV file."),
          type = "warning",
          duration = 8  # Show warning for 8 seconds
        )
      }
    }
  })

  # Combined pilot data file upload handling
  custom_baseline <- reactiveVal(NULL)
  custom_library <- reactiveVal(NULL)
  
  observeEvent(input$pilot_data_file, {
    req(input$pilot_data_file)
    # Require file if either "Other" is selected OR custom pilot data is chosen
    req(input$pilot_data_choice == "custom" || input$biological_system == "Other" || input$experimental_platform == "Other")
    
    # Check file size (limit to 50MB for RDS files)
    file_size_mb <- file.size(input$pilot_data_file$datapath) / (1024^2)
    if (file_size_mb > 50) {
      showNotification(
        paste("File size (", round(file_size_mb, 1), "MB) exceeds the 50MB limit. Please use a smaller dataset."),
        type = "error", duration = 10
      )
      custom_baseline(NULL)
      output$baseline_uploaded <- reactive(FALSE)
      outputOptions(output, "baseline_uploaded", suspendWhenHidden = FALSE)
      return()
    }
    
    tryCatch({
      # Read the uploaded RDS file
      file_ext <- tolower(tools::file_ext(input$pilot_data_file$name))
      
      if (file_ext == "rds") {
        # Read RDS file
        uploaded_data <- readRDS(input$pilot_data_file$datapath)
        
        # Validate the combined pilot data structure
        validation_result <- perturbplan::validate_combined_pilot_data(uploaded_data, input$pilot_data_file$name)
        
        if (validation_result$valid) {
          # Extract components from validated data
          custom_baseline(validation_result$data$baseline_expression)
          custom_library(validation_result$data$library_parameters)
          
          # Create success message with summary and warnings
          status_msg <- validation_result$summary
          if (length(validation_result$warnings) > 0) {
            warning_msg <- paste0("<br/><em style='color:orange;'>", 
                                paste(validation_result$warnings, collapse = "<br/>"), 
                                "</em>")
            status_msg <- paste0(status_msg, warning_msg)
          }
          
          # Update status display
          output$pilot_data_status <- renderUI({
            HTML(status_msg)
          })
          
          output$pilot_data_uploaded <- reactive(TRUE)
          outputOptions(output, "pilot_data_uploaded", suspendWhenHidden = FALSE)
          
        } else {
          # Show validation errors
          error_msg <- paste0("Validation failed:<br/>", 
                            paste(validation_result$errors, collapse = "<br/>"))
          
          output$pilot_data_status <- renderUI({
            HTML(paste0("<em style='color:red;'>", error_msg, "</em>"))
          })
          
          custom_baseline(NULL)
          custom_library(NULL)
          output$pilot_data_uploaded <- reactive(FALSE)
          outputOptions(output, "pilot_data_uploaded", suspendWhenHidden = FALSE)
        }
        
      } else {
        showNotification("Please upload an RDS file with the required pilot data structure", type = "error")
        custom_baseline(NULL)
        custom_library(NULL)
        output$pilot_data_uploaded <- reactive(FALSE)
        outputOptions(output, "pilot_data_uploaded", suspendWhenHidden = FALSE)
      }
      
    }, error = function(e) {
      # Enhanced error handling for different types of errors
      error_msg <- if (grepl("cannot open the connection", e$message, ignore.case = TRUE)) {
        "Cannot read the uploaded file. Please ensure it's a valid RDS file."
      } else if (grepl("magic number", e$message, ignore.case = TRUE) || grepl("format", e$message, ignore.case = TRUE)) {
        "File appears to be corrupted or not a valid RDS file. Please check the file format."
      } else if (grepl("version", e$message, ignore.case = TRUE)) {
        "RDS file was created with a newer version of R. Please re-save the file with your current R version."
      } else {
        paste("Error reading baseline file:", e$message)
      }
      
      showNotification(error_msg, type = "error", duration = 10)
      output$pilot_data_status <- renderUI({
        HTML(paste0("<em style='color:red;'>", error_msg, "</em>"))
      })
      custom_baseline(NULL)
      custom_library(NULL)
      output$pilot_data_uploaded <- reactive(FALSE)
      outputOptions(output, "pilot_data_uploaded", suspendWhenHidden = FALSE)
    })
  })
  
  # Reset pilot data when choice changes to default AND neither biological_system nor experimental_platform is "Other"
  observe({
    # Reset custom data if file is removed OR (pilot_data_choice is "default" AND neither system/platform is "Other")
    should_reset <- is.null(input$pilot_data_file) || 
                   (input$pilot_data_choice == "default" && 
                    input$biological_system != "Other" && 
                    input$experimental_platform != "Other")
    
    if (should_reset) {
      custom_baseline(NULL)
      custom_library(NULL)
      output$pilot_data_uploaded <- reactive(FALSE)
      outputOptions(output, "pilot_data_uploaded", suspendWhenHidden = FALSE)
    }
  })

  # Extract fold-change and expression information first
  fc_expression_info <- reactive({
    req(planned())
    # Extract expression information including gene list if provided
    perturbplan::extract_fc_expression_info(
      fold_change_mean = input$fc_mean,
      fold_change_sd = input$fc_sd,
      biological_system = input$biological_system,
      B = 1000,  # Monte Carlo samples for good accuracy
      gene_list = if(input$gene_list_mode == "custom") gene_list() else NULL,  # Use gene list only in custom mode
      tpm_threshold = input$tpm_threshold,  # Apply TPM filtering
      custom_baseline_data = if(input$pilot_data_choice == "custom" || input$biological_system == "Other" || input$experimental_platform == "Other") custom_baseline() else NULL  # Use custom baseline if provided
    )
  })

  # Extract library information (UMI saturation parameters)
  library_info <- reactive({
    req(planned())
    # Extract library parameters for biological system or use custom data from pilot data
    perturbplan::extract_library_info(
      biological_system = input$biological_system,
      custom_library_data = if(input$pilot_data_choice == "custom" || input$biological_system == "Other" || input$experimental_platform == "Other") custom_library() else NULL
    )
  })

  # Generate optimized cell and read ranges using modular approach
  cells_reads_df <- reactive({
    req(planned(), fc_expression_info(), library_info())
    
    # Call identify_cell_read_range to generate experimental design
    cell_read_range <- perturbplan::identify_cell_read_range(
      experimental_platform = input$experimental_platform,
      fc_expression_info = fc_expression_info(),
      library_info = library_info(),
      grid_size = 10,  # 10x10 grid for Shiny app
      min_power_threshold = 0.05,  # Moderate minimum for good visualization range
      max_power_threshold = 0.9,   # High maximum for wide power range
      MOI = input$MOI,
      num_targets = if(input$gene_list_mode == "custom" && !is.null(target_list())) {
        length(unique(target_list()))  # Use actual unique targets from CSV
      } else {
        input$num_targets  # Use user input for random mode
      },
      gRNAs_per_target = input$gRNAs_per_target,
      non_targeting_gRNAs = input$non_targeting_gRNAs,
      control_group = input$control_group,
      multiple_testing_alpha = input$fdr_target,
      side = input$side,
      prop_non_null = input$prop_non_null
    )
    
    # Convert to format expected by calculate_power_grid
    perturbplan::convert_design_to_dataframe(cell_read_range)
  })

  # Real power calculation using perturbplan package with pre-computed design
  power_results <- reactive({
    req(planned(), fc_expression_info(), cells_reads_df())
    # Call the package function with pre-computed experimental design
    perturbplan::calculate_power_grid(
      fc_expression_info = fc_expression_info(),
      num_targets = if(input$gene_list_mode == "custom" && !is.null(target_list())) {
        length(unique(target_list()))  # Use actual unique targets from CSV
      } else {
        input$num_targets  # Use user input for random mode
      },
      gRNAs_per_target = input$gRNAs_per_target, 
      non_targeting_gRNAs = input$non_targeting_gRNAs,
      fdr_target = input$fdr_target,
      prop_non_null = input$prop_non_null,
      MOI = input$MOI,
      side = input$side,
      control_group = input$control_group,
      cells_reads_df = cells_reads_df()  # Pass pre-computed experimental design
    )
  })

  # Grid setup - get from power calculation results
  cells_seq <- reactive({
    req(planned())
    power_results()$cells_seq
  })
  
  reads_seq <- reactive({
    req(planned())
    power_results()$reads_seq
  })
  
  dcells <- reactive({
    seq_vals <- cells_seq()
    diff(seq_vals)[1]
  })
  
  dreads <- reactive({
    seq_vals <- reads_seq()
    diff(seq_vals)[1]
  })

  gridDF <- reactive({
    req(planned())
    power_results()$power_grid
  })

  # Return all reactive values and functions for use by other modules
  return(list(
    planned = planned,
    fc_expression_info = fc_expression_info,  # Add fc_expression_info for curves server
    library_info = library_info,  # Add library_info for curves server
    power_results = power_results,
    cells_seq = cells_seq,
    reads_seq = reads_seq,
    dcells = dcells,
    dreads = dreads,
    gridDF = gridDF,
    gene_list = gene_list
  ))
}