# Tab panel UI modules for different views

create_overall_power_tab <- function() {
  tabPanel(
    "Overall power",
    value = "overall_power",
    fluidRow(
      box(
        status = "primary",
        solidHeader = TRUE,
        width = 8,
        height = "570px",
        
        # Sub-tabs inside the box
        tabsetPanel(
          id = "overall_subtabs",
          
          # Heatmap sub-tab
          tabPanel(
            "Heatmap",
            value = "overall_heatmap",
            
            # Message before planning
            conditionalPanel(
              condition = "output.need_plan",
              h4("Select your parameters and click \"Plan\"")
            ),
            
            # Heatmap (after Plan)
            conditionalPanel(
              condition = "!output.need_plan",
              plotOutput("heat", height = 480, click = "heat_click")
            )
          ),
          
          # Slice sub-tab
          tabPanel(
            "Slice",
            value = "overall_slice",
            
            # Message when no slices selected
            conditionalPanel(
              condition = "output.need_plan || !output.slice_available",
              h4("Select slices to plot on the heatmap")
            ),
            
            # Slice plot (when slices are available)
            conditionalPanel(
              condition = "!output.need_plan && output.slice_available",
              plotOutput("slice_plot", height = 480, click = "slice_click")
            )
          )
        )
      ),
      
      conditionalPanel(
        condition = "!output.need_plan",
        box(
          status = "info", 
          solidHeader = TRUE,
          width = 4,
          height = "570px",
          
          # Heatmap-specific controls (only show when on heatmap tab)
          conditionalPanel(
            condition = "input.overall_subtabs == 'overall_heatmap'",
            radioButtons(
              "mode", "Drill down by:",
              c("Number of cells per target (click one or more rows)" = "cells",
                "Reads per cell (click one or more columns)" = "reads", 
                "Both (click one or more tiles)" = "tile"),
              selected = "cells"
            ),
            uiOutput("overall_box_ui"),
            actionButton("go_overall", "Go", class = "btn-primary", width = "100%"),
            br(), br(),
            actionButton("clear", "Clear", width = "100%")
          ),
          
          # Slice-specific controls (only show when on slice tab)
          conditionalPanel(
            condition = "input.overall_subtabs == 'overall_slice'",
            uiOutput("slice_box_ui"),
            actionButton("go_slice", "Go", class = "btn-primary", width = "100%"),
            br(), br(),
            actionButton("slice_clear", "Clear", width = "100%")
          )
        )
      )
    )
  )
}

create_curves_tab <- function() {
  tabPanel(
    "Drill-down power",
    value = "per_pair",
    fluidRow(
      box(
        status = "primary",
        solidHeader = TRUE,
        width = 8,
        height = "570px",
        
        # Subtabs for different plot types
        tabsetPanel(
          id = "curves_subtabs",
          tabPanel(
            "Expression",
            value = "tpm_plot",
            plotOutput("pp_tpm", height = "480px", width = "100%")
          ),
          tabPanel(
            "Fold Change",
            value = "fc_plot", 
            plotOutput("pp_fc", height = "480px", width = "100%")
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
          height = "570px",
          
          radioButtons(
            "curves_display_mode", "Display designs:",
            c("All together" = "all_together",
              "Facet over cells per target" = "facet_cells", 
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
    create_overall_power_tab(),
    create_curves_tab()
  )
}