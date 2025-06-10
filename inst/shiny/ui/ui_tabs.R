# Tab panel UI modules for different views

create_heatmap_tab <- function() {
  tabPanel(
    "Overall power (heatmap)",
    value = "overall_heatmap",
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
      
      conditionalPanel(
        condition = "!output.need_plan",
        box(
          status = "info", 
          solidHeader = TRUE,
          width = 4,
          
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
  )
}

create_slice_tab <- function() {
  tabPanel(
    "Overall power (slice)",
    value = "overall_slice",
    fluidRow(
      box(
        title = textOutput("slice_title"),
        status = "primary",
        solidHeader = TRUE,
        width = 8,
        plotOutput("slice_plot", height = 460, click = "slice_click")
      ),
      
      box(
        status = "info",
        solidHeader = TRUE, 
        width = 4,
        uiOutput("slice_box_ui"),
        actionButton("go_slice", "Go", class = "btn-primary", width = "100%"),
        br(), br(),
        actionButton("slice_clear", "Clear", width = "100%")
      )
    )
  )
}

create_curves_tab <- function() {
  tabPanel(
    "Power over TPM & FC",
    value = "per_pair",
    fluidRow(
      box(
        title = "Power over TPM & FC",
        status = "primary",
        solidHeader = TRUE,
        width = 8,
        
        # Subtabs for different plot types
        tabsetPanel(
          id = "curves_subtabs",
          tabPanel(
            "Expression (TPM)",
            value = "tpm_plot",
            plotOutput("pp_tpm", height = "600px", width = "100%")
          ),
          tabPanel(
            "Fold Change",
            value = "fc_plot", 
            plotOutput("pp_fc", height = "600px", width = "100%")
          )
        )
      ),
      
      # Control box for display options
      conditionalPanel(
        condition = "output.curves_available",
        box(
          status = "info",
          solidHeader = TRUE,
          width = 4,
          
          radioButtons(
            "curves_display_mode", "Display designs:",
            c("All together" = "all_together",
              "Facet over cells" = "facet_cells", 
              "Facet over reads per cell" = "facet_reads"),
            selected = "all_together"
          )
        )
      )
    )
  )
}

create_main_tabs <- function() {
  tabsetPanel(
    id = "main_tabs",
    create_heatmap_tab(),
    create_slice_tab(),
    create_curves_tab()
  )
}