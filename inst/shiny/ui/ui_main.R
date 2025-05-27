# Main UI layout - converted to shinydashboard
# Uses dashboard structure with sidebar and main body

ui <- dashboardPage(
  # Header with title and download button
  dashboardHeader(
    title = "PerturbPlan",
    # Download button in header
    tags$li(class = "dropdown",
      conditionalPanel(
        condition = "!output.need_plan",
        downloadButton("download_data", "Download Results", 
                      class = "btn-success", 
                      style = "margin-top: 8px; margin-right: 10px;")
      )
    )
  ),
  
  # Sidebar with experimental parameters
  dashboardSidebar(
    useShinyjs(),
    width = 275,
    
    # Add custom CSS for better contrast and readability
    tags$head(tags$style(HTML("
      /* Sidebar styling */
      .main-sidebar, .left-side {
        background-color: #f8f9fa !important;
      }
      
      /* Make header logo area same width as sidebar */
      .main-header .logo {
        width: 275px !important;
      }
      
      .main-header .navbar {
        margin-left: 275px !important;
      }
      
      /* Sidebar text */
      .sidebar-menu, .sidebar-menu > li > a, .sidebar-menu > li.header {
        color: #212529 !important;
        background-color: transparent !important;
      }
      
      /* Form elements in sidebar */
      .main-sidebar .form-control {
        background-color: #ffffff !important;
        border: 1px solid #ced4da !important;
        color: #212529 !important;
      }
      
      /* Labels */
      .main-sidebar label {
        color: #495057 !important;
        font-weight: 500;
      }
      
      /* Collapse panels - light gray background for headers */
      .main-sidebar .panel-default > .panel-heading {
        background-color: #f8f9fa !important;
        border: 1px solid #dee2e6 !important;
        color: #212529 !important;
        font-weight: 600 !important;
        position: relative !important;
      }
      
      .main-sidebar .panel-default > .panel-heading a {
        color: #212529 !important;
        text-decoration: none !important;
        font-weight: 600 !important;
        display: block !important;
        padding-right: 30px !important;
      }
      
      /* Add down arrow indicators */
      .main-sidebar .panel-default > .panel-heading a:after {
        content: '▼' !important;
        position: absolute !important;
        right: 15px !important;
        top: 50% !important;
        transform: translateY(-50%) !important;
        font-size: 12px !important;
        color: #6c757d !important;
        transition: transform 0.3s ease !important;
      }
      
      .main-sidebar .panel-default > .panel-heading.collapsed a:after {
        transform: translateY(-50%) rotate(-90deg) !important;
      }
      
      .main-sidebar .panel-body {
        background-color: #ffffff !important;
        border-color: #dee2e6 !important;
        padding: 10px 15px !important;
      }
      
      /* Buttons */
      .main-sidebar .btn {
        margin-bottom: 10px;
      }
      
      /* Content area */
      .content-wrapper {
        background-color: #ffffff !important;
      }
      
      /* Reduced spacing between input boxes */
      .main-sidebar .form-group {
        margin-bottom: 8px !important;
      }
      
      /* Select inputs */
      .main-sidebar .selectize-input {
        background-color: #ffffff !important;
        border: 1px solid #ced4da !important;
        color: #212529 !important;
      }
    "))),
    
    # Experimental parameters
    bsCollapse(
      id   = "panels",
      open = "Experimental choices",

      # --- Experimental choices -----------------------------------
      bsCollapsePanel(
        "Experimental choices",
        selectInput("biological_system",   "Biological system:",   c("K562")),
        selectInput("experimental_platform","Experimental platform:", c("10x Chromium v3")),
        numericInput("MOI",               "MOI:",                 10, 1, 30, 0.5),
        numericInput("num_targets",       "Number of targets:",   100, 1, 1000),
        numericInput("gRNAs_per_target",  "Number of gRNAs per target:",     4, 1, 10),
        numericInput("non_targeting_gRNAs","Number of non-targeting gRNAs:",10, 0, 100)
      ),

      # --- Analysis choices ---------------------------------------
      bsCollapsePanel(
        "Analysis choices",
        numericInput("num_pairs",     "Number of pairs analyzed:",  1000, 10, 10000),
        numericInput("tpm_threshold", "Minimum TPM threshold:",      10, 0, 10, 0.5),
        numericInput("fdr_target",    "FDR target level:", 0.05, 0.001, 0.1, 0.001)
      ),

      # --- Assumed effect sizes -----------------------------------
      bsCollapsePanel(
        "Assumed effect sizes",
        numericInput("fc_mean", "Fold-change mean:", 0.85, 1.1, 10, 0.05),
        numericInput("fc_sd",   "Fold-change SD:",   0.15, 0.1,  5, 0.05),
        numericInput("prop_non_null", "Proportion of non-null pairs:", 0.1, 0, 1, 0.01)
      )
    ),

    hr(),
    actionButton("plan_btn", "Plan", class = "btn-success", width = "90%")
  ),
  
  # Main dashboard body
  dashboardBody(
    useShinyjs(),
    source("ui/main_tabs.R", local = TRUE)$value
  )
)