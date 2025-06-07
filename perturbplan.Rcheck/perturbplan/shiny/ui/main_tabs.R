# Main content tabs
# Contains the three main visualization tabs

tabsetPanel(
  id = "tabs",

  # === Overall tab ==============================================
  tabPanel(
    "Overall power (heatmap)",

    # Message before planning
    conditionalPanel(
      condition = "output.need_plan",
      h4("Select your parameters and click \"Plan\"")
    ),

    # Heat-map and controls (after Plan)
    conditionalPanel(
      condition = "!output.need_plan",
      h4("Power versus number of cells and reads per cell"),
      fluidRow(
        column(
          8,
          plotOutput("heat", height = 460, click = "heat_click")
        ),
        column(
          4,
          wellPanel(
            radioButtons(
              "mode", "Drill down by:",
              c("Number of cells (click one or more rows)"   = "cells",
                "Reads per cell (click one or more columns)" = "reads",
                "Both (click one or more tiles)"             = "tile"),
              selected = "cells"
            ),
            uiOutput("overall_box_ui"),
            actionButton("go_overall", "Go",   class = "btn-primary", width = "100%"),
            actionButton("clear",      "Clear",                       width = "100%")
          )
        )
      )
    )
  ),

  # === Slice tab ===============================================
  tabPanel(
    "Overall power (slice)",
    h4(textOutput("slice_title")),
    fluidRow(
      column(
        8,
        plotOutput("slice_plot", height = 460, click = "slice_click")
      ),
      column(
        4,
        wellPanel(
          uiOutput("slice_box_ui"),
          actionButton("go_slice",    "Go",   class = "btn-primary", width = "100%"),
          actionButton("slice_clear", "Clear",                       width = "100%")
        )
      )
    )
  ),

  # === Per-pair tab ============================================
  tabPanel(
    "Per-pair power",
    h4(textOutput("pair_title")),
    fluidRow(
      column(
        8,
        plotOutput("pp_combined", height = "400px", width = "100%")
      ),
      column(
        4
        # empty spacer
      )
    )
  )
)