# Server module for power calculations and reactive values

create_power_server <- function(input, output, session) {
  # Reactive values
  planned <- reactiveVal(FALSE)
  output$need_plan <- reactive(!planned())
  outputOptions(output, "need_plan", suspendWhenHidden = FALSE)
  observeEvent(input$plan_btn, planned(TRUE))
  
  # Gene list file upload handling
  gene_list <- reactiveVal(NULL)
  
  observeEvent(input$gene_list_file, {
    req(input$gene_list_file)
    
    tryCatch({
      # Read the uploaded file
      file_ext <- tools::file_ext(input$gene_list_file$name)
      
      if (file_ext %in% c("csv", "txt")) {
        # Try reading as CSV first, then as simple text
        if (file_ext == "csv") {
          uploaded_data <- read.csv(input$gene_list_file$datapath, 
                                  header = FALSE, 
                                  stringsAsFactors = FALSE)
          genes <- as.character(uploaded_data[,1])
        } else {
          genes <- readLines(input$gene_list_file$datapath)
        }
        
        # Clean the gene list
        genes <- trimws(genes)  # Remove whitespace
        genes <- genes[genes != ""]  # Remove empty lines
        genes <- unique(genes)  # Remove duplicates
        
        # Store the gene list
        gene_list(genes)
        
        # Update status
        output$gene_list_status <- renderText({
          sprintf("Loaded %d unique genes", length(genes))
        })
        
        output$gene_list_uploaded <- reactive(TRUE)
        outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
        
      } else {
        showNotification("Please upload a CSV or TXT file", type = "error")
        gene_list(NULL)
        output$gene_list_uploaded <- reactive(FALSE)
        outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
      }
      
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      gene_list(NULL)
      output$gene_list_uploaded <- reactive(FALSE)
      outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
    })
  })
  
  # Reset gene list when file is removed
  observe({
    if (is.null(input$gene_list_file)) {
      gene_list(NULL)
      output$gene_list_uploaded <- reactive(FALSE)
      outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
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
      gene_list = gene_list()  # Use uploaded gene list if available
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
      num_targets = input$num_targets,
      gRNAs_per_target = input$gRNAs_per_target, 
      non_targeting_gRNAs = input$non_targeting_gRNAs,
      tpm_threshold = input$tpm_threshold,
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