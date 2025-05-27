# Main UI layout
# Extracted from app-prototype.R and modularized

ui <- fluidPage(
  useShinyjs(),

  sidebarLayout(
    # ------------------------------------------------------------------
    sidebarPanel(
      width = 3,
      tags$div(
        style = "overflow-y:auto; height:92vh; position:sticky; top:0;",
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
        actionButton("plan_btn", "Plan", class = "btn-success", width = "100%")
      )
    ),

    # ------------------------------------------------------------------
    mainPanel(
      width = 9,
      source("ui/main_tabs.R", local = TRUE)$value
    )
  )
)