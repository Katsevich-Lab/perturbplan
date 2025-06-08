# UI definition using shinydashboard structure
# Updated with new parameters and features while maintaining modular structure

ui <- dashboardPage(
  # Dashboard header
  dashboardHeader(title = "Perturb-seq Power Planner"),
  
  # Dashboard sidebar
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      menuItem("Overall power (heatmap)", tabName = "overall_heatmap", icon = icon("th")),
      menuItem("Overall power (slice)", tabName = "overall_slice", icon = icon("chart-line")),
      menuItem("Power over TPM & FC", tabName = "per_pair", icon = icon("project-diagram"))
    ),
    
    # Parameter panels - make scrollable while keeping right panel fixed
    tags$div(
      style = "padding: 15px; height: calc(100vh - 100px); overflow-y: auto; position: fixed; width: 230px;",
      
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
  
  # Dashboard body - fixed position, no scrolling
  dashboardBody(
    useShinyjs(),
    
    # Custom CSS to ensure right panel is fixed
    tags$head(
      tags$style(HTML("
        .content-wrapper {
          position: fixed;
          top: 50px;
          left: 230px;
          right: 0;
          bottom: 0;
          overflow: hidden;
        }
        .content {
          height: 100%;
          overflow-y: auto;
          padding: 15px;
        }
      "))
    ),
    
    # Tab content
    tabItems(
      # Overall power heatmap tab
      tabItem(
        tabName = "overall_heatmap",
        fluidRow(
          box(
            title = "Power versus number of cells and reads per cell",
            status = "primary",
            solidHeader = TRUE,
            width = 8,
            
            # Message before planning
            conditionalPanel(
              condition = "output.need_plan",
              h4("Select your parameters and click \"Plan\"")
            ),
            
            # Heatmap (after Plan)
            conditionalPanel(
              condition = "!output.need_plan",
              plotOutput("heat", height = 460, click = "heat_click")
            )
          ),
          
          box(
            title = "Controls",
            status = "info", 
            solidHeader = TRUE,
            width = 4,
            
            conditionalPanel(
              condition = "!output.need_plan",
              radioButtons(
                "mode", "Drill down by:",
                c("Number of cells (click one or more rows)" = "cells",
                  "Reads per cell (click one or more columns)" = "reads", 
                  "Both (click one or more tiles)" = "tile"),
                selected = "cells"
              ),
              uiOutput("overall_box_ui"),
              actionButton("go_overall", "Go", class = "btn-primary", width = "100%"),
              br(), br(),
              actionButton("clear", "Clear", width = "100%")
            )
          )
        )
      ),
      
      # Overall power slice tab
      tabItem(
        tabName = "overall_slice",
        fluidRow(
          box(
            title = textOutput("slice_title"),
            status = "primary",
            solidHeader = TRUE,
            width = 8,
            plotOutput("slice_plot", height = 460, click = "slice_click")
          ),
          
          box(
            title = "Slice Controls",
            status = "info",
            solidHeader = TRUE, 
            width = 4,
            uiOutput("slice_box_ui"),
            actionButton("go_slice", "Go", class = "btn-primary", width = "100%"),
            br(), br(),
            actionButton("slice_clear", "Clear", width = "100%")
          )
        )
      ),
      
      # Power over TPM & FC tab
      tabItem(
        tabName = "per_pair",
        fluidRow(
          box(
            title = "Power over TPM & FC",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("pp_combined", height = "400px", width = "100%")
          )
        )
      )
    )
  )
)