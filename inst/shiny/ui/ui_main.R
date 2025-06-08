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
    
    # Parameter panels with collapsible sections
    tags$div(
      style = "padding: 10px; overflow-y: auto; height: calc(100vh - 100px);",
      
      # Experimental choices - collapsible
      box(
        title = "Experimental choices", 
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        selectInput("biological_system", "Biological system:", c("K562")),
        selectInput("experimental_platform", "Experimental platform:", c("10x Chromium v3")),
        fluidRow(
          column(6, numericInput("MOI", "MOI:", 10, 1, 30, 0.5)),
          column(6, numericInput("num_targets", "Targets:", 100, 1, 1000))
        ),
        fluidRow(
          column(6, numericInput("gRNAs_per_target", "gRNAs/target:", 4, 1, 10)),
          column(6, numericInput("non_targeting_gRNAs", "Non-targeting:", 10, 0, 100))
        )
      ),
      
      # Analysis choices - collapsible
      box(
        title = "Analysis choices", 
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        fluidRow(
          column(6, numericInput("tpm_threshold", "TPM threshold:", 10, 0, 10, 0.5)),
          column(6, numericInput("fdr_target", "FDR target:", 0.05, 0.001, 0.1, 0.001))
        ),
        selectInput("side", "Test side:", 
                    choices = c("Left (knockdown)" = "left", 
                               "Right (overexpression)" = "right"), 
                    selected = "left"),
        selectInput("control_group", "Control group:", 
                    choices = c("Complement cells" = "complement", 
                               "Non-targeting cells" = "nt_cells"), 
                    selected = "complement")
      ),
      
      # Effect sizes - collapsible  
      box(
        title = "Assumed effect sizes", 
        status = "warning",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        fluidRow(
          column(6, numericInput("fc_mean", "FC mean:", 0.85, 1.1, 10, 0.05)),
          column(6, numericInput("fc_sd", "FC SD:", 0.15, 0.1, 5, 0.05))
        ),
        numericInput("prop_non_null", "Proportion of non-null pairs:", 0.1, 0, 1, 0.01)
      ),
      
      # Prominent Plan button
      br(),
      actionButton("plan_btn", "Plan Analysis", 
                   class = "btn-success btn-lg", 
                   width = "100%",
                   style = "font-size: 18px; font-weight: bold;")
    )
  ),
  
  # Dashboard body
  dashboardBody(
    useShinyjs(),
    
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
              div(
                class = "alert alert-info",
                style = "margin: 20px;",
                icon("info-circle"),
                " Select your parameters in the sidebar and click \"Plan Analysis\" to generate the power heatmap."
              )
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