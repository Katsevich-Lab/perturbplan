# Main server logic
# Extracted from app-prototype.R and modularized

server <- function(input, output, session) {
  
  # Source server modules
  source("server/reactive_data.R", local = TRUE)
  source("server/selection_logic.R", local = TRUE)  
  source("server/plot_rendering.R", local = TRUE)
  
  # Core reactive values and helpers
  is_sel <- function(x) identical(sel$type, x)
  toggle  <- function(v, x) if (x %in% v) setdiff(v, x) else c(v, x)
  planned <- reactiveVal(FALSE)
  
  # UI state management
  output$need_plan <- reactive(!planned())
  outputOptions(output,  "need_plan", suspendWhenHidden = FALSE)
  observeEvent(input$plan_btn, planned(TRUE))
  
  # Selection state
  sel <- reactiveValues(type=NULL, idx=integer(0),
                        tiles=data.frame(row=integer(0), col=integer(0)))
  slice_mode <- reactiveVal(NULL)
  slice_x    <- reactiveVal(numeric(0))
  
  # Clear selections
  observeEvent(input$clear, {
    sel$type <- NULL; sel$idx <- integer(0); sel$tiles <- sel$tiles[0,]
    slice_mode(NULL); slice_x(numeric(0))
    updateTextInput(session,"overall_points","")
    updateTextInput(session,"slice_points","")
  })
  
  # Download handler for power grid data with parameters
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("perturbplan_results_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(planned())
      
      # Get power results and grid data
      results <- power_results()
      grid_data <- results$power_grid
      params <- results$parameters
      
      # Create header with parameters as comments
      header_lines <- c(
        "# PerturbPlan Power Analysis Results",
        paste("# Generated on:", Sys.time()),
        "# ",
        "# Experimental Parameters:",
        paste("# Number of targets:", params$num_targets),
        paste("# gRNAs per target:", params$gRNAs_per_target),
        paste("# Non-targeting gRNAs:", params$non_targeting_gRNAs),
        paste("# Number of pairs analyzed:", params$num_pairs),
        paste("# TPM threshold:", params$tpm_threshold),
        paste("# FDR target level:", params$fdr_target),
        paste("# Fold-change mean:", params$fc_mean),
        paste("# Fold-change SD:", params$fc_sd),
        paste("# Proportion of non-null pairs:", params$prop_non_null),
        paste("# MOI:", params$MOI),
        paste("# Biological system:", params$biological_system),
        paste("# Experimental platform:", params$experimental_platform),
        "# ",
        "# Data columns: cells, reads, power"
      )
      
      # Write header and data
      writeLines(header_lines, file)
      write.table(grid_data, file, append = TRUE, sep = "\t", 
                  row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  )
  
  # Make reactive values available to other modules
  values <- reactiveValues(
    sel = sel,
    slice_mode = slice_mode,
    slice_x = slice_x,
    planned = planned,
    is_sel = is_sel,
    toggle = toggle
  )
}