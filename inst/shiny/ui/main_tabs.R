# Main tab panels for the dashboard
# Contains all tab content with real power calculation integration

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