# Server module for power calculations and reactive values

create_power_server <- function(input, output, session) {
  # Reactive values
  planned <- reactiveVal(FALSE)
  output$need_plan <- reactive(!planned())
  outputOptions(output, "need_plan", suspendWhenHidden = FALSE)
  observeEvent(input$plan_btn, planned(TRUE))
  
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
          baseline_expression_stats <- perturbplan::extract_baseline_expression(biological_system = input$biological_system)
          baseline_df <- baseline_expression_stats$baseline_expression
          
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
          
          # Create detailed status message
          status_msg <- sprintf("Loaded %d pairs (%d unique genes, %d unique targets)", 
                               total_pairs, unique_genes_uploaded, unique_targets)
          
          filtering_msg <- sprintf("Gene filtering: %.1f%% not in database, %.1f%% filtered by TPM threshold (≥%d), %.1f%% available for analysis",
                                  prop_not_in_database * 100, prop_filtered_by_tpm * 100, input$tpm_threshold, prop_available * 100)
          
          # Update status with both messages
          output$gene_list_status <- renderUI({
            HTML(paste(status_msg, filtering_msg, sep = "<br/>"))
          })
          
          # Show warnings for significant filtering
          if (prop_not_in_database > 0.1) {  # More than 10% missing
            showNotification(
              paste0("Warning: ", round(prop_not_in_database * 100, 1), "% of unique genes (", length(genes_not_in_database), 
                     " genes) were not found in the ", input$biological_system, " expression database."),
              type = "warning",
              duration = 6
            )
          }
          
          if (prop_filtered_by_tpm > 0.1) {  # More than 10% filtered by TPM
            showNotification(
              paste0("Warning: ", round(prop_filtered_by_tpm * 100, 1), "% of unique genes (", length(genes_filtered_by_tpm), 
                     " genes) were filtered out due to TPM threshold ≥ ", input$tpm_threshold, "."),
              type = "warning", 
              duration = 6
            )
          }
          
        }, error = function(e) {
          # Fallback to basic status if filtering analysis fails
          output$gene_list_status <- renderUI({
            HTML(sprintf("Loaded %d pairs (%d unique genes, %d unique targets)", total_pairs, unique_genes_uploaded, unique_targets))
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
      tpm_threshold = input$tpm_threshold  # Apply TPM filtering
    )
  })

  # Extract library information (UMI saturation parameters)
  library_info <- reactive({
    req(planned())
    # Extract library parameters for biological system
    perturbplan::extract_library_info(
      biological_system = input$biological_system
    )
  })

  # Real power calculation using perturbplan package
  power_results <- reactive({
    req(planned(), fc_expression_info(), library_info())
    # Call the package function with extracted expression and library info
    perturbplan::calculate_power_grid(
      fc_expression_info = fc_expression_info(),
      library_info = library_info(),
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
      experimental_platform = input$experimental_platform,
      side = input$side,
      control_group = input$control_group
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