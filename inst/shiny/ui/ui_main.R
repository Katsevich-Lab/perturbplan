# UI definition using shinydashboard structure
# Updated with new parameters and features while maintaining modular structure

ui <- dashboardPage(
  # Dashboard header
  dashboardHeader(title = "Perturb-seq Power Planner"),
  
  # Dashboard sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overall power (heatmap)", tabName = "overall_heatmap", icon = icon("th")),
      menuItem("Overall power (slice)", tabName = "overall_slice", icon = icon("chart-line")),
      menuItem("Power over TPM & FC", tabName = "per_pair", icon = icon("project-diagram"))
    ),
    
    # Parameter panels
    tags$div(
      style = "padding: 15px;",
      
      # Experimental choices
      h4("Experimental choices"),
      selectInput("biological_system", "Biological system:", c("K562")),
      selectInput("experimental_platform", "Experimental platform:", c("10x Chromium v3")),
      numericInput("MOI", "MOI:", 10, 1, 30, 0.5),
      numericInput("num_targets", "Number of targets:", 100, 1, 1000),
      numericInput("gRNAs_per_target", "Number of gRNAs per target:", 4, 1, 10),
      numericInput("non_targeting_gRNAs", "Number of non-targeting gRNAs:", 10, 0, 100),
      
      br(),
      
      # Analysis choices
      h4("Analysis choices"),
      numericInput("tpm_threshold", "Minimum TPM threshold:", 10, 0, 10, 0.5),
      numericInput("fdr_target", "FDR target level:", 0.05, 0.001, 0.1, 0.001),
      selectInput("side", "Test side:", 
                  choices = c("Left (knockdown)" = "left", 
                             "Right (overexpression)" = "right"), 
                  selected = "left"),
      selectInput("control_group", "Control group:", 
                  choices = c("Complement cells" = "complement", 
                             "Non-targeting cells" = "nt_cells"), 
                  selected = "complement"),
      
      br(),
      
      # Assumed effect sizes
      h4("Assumed effect sizes"),
      numericInput("fc_mean", "Fold-change mean:", 0.85, 1.1, 10, 0.05),
      numericInput("fc_sd", "Fold-change SD:", 0.15, 0.1, 5, 0.05),
      numericInput("prop_non_null", "Proportion of non-null pairs:", 0.1, 0, 1, 0.01),
      
      br(),
      
      actionButton("plan_btn", "Plan", class = "btn-success", width = "100%")
    )
  ),
  
  # Dashboard body
  dashboardBody(
    useShinyjs(),
    
    # Include main tabs
    source("ui/main_tabs.R", local = TRUE)$value
  )
)