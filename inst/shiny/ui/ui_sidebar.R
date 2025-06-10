# Sidebar UI module for parameter inputs

create_sidebar <- function() {
  dashboardSidebar(
    # Parameter panels - make scrollable with collapsible sections
    tags$div(
      style = "padding: 15px; max-height: 85vh; overflow-y: auto;",
      
      # Experimental choices (collapsible)
      tags$div(
        style = "border: 1px solid #ddd; border-radius: 4px; margin-bottom: 10px;",
        tags$div(
          id = "exp-header",
          style = "background-color: #f5f5f5; padding: 10px 15px; cursor: pointer; border-bottom: 1px solid #ddd; color: #333;",
          onclick = "toggleSection('experimental-content', 'exp-chevron')",
          tags$i(id = "exp-chevron", class = "fa fa-chevron-down", style = "margin-right: 8px; color: #333;"),
          tags$strong("Experimental choices")
        ),
        tags$div(
          id = "experimental-content",
          style = "padding: 15px;",
          selectInput("biological_system", "Biological system:", c("K562")),
          selectInput("experimental_platform", "Experimental platform:", c("10x Chromium v3")),
          numericInput("MOI", "MOI:", 10, 1, 30, 0.5),
          numericInput("num_targets", "Number of targets:", 100, 1, 1000),
          numericInput("gRNAs_per_target", "Number of gRNAs per target:", 4, 1, 10),
          numericInput("non_targeting_gRNAs", "Number of non-targeting gRNAs:", 10, 0, 100)
        )
      ),
      
      br(),
      
      # Analysis choices (collapsible)
      tags$div(
        style = "border: 1px solid #ddd; border-radius: 4px; margin-bottom: 10px;",
        tags$div(
          id = "analysis-header",
          style = "background-color: #f5f5f5; padding: 10px 15px; cursor: pointer; border-bottom: 1px solid #ddd; color: #333;",
          onclick = "toggleSection('analysis-content', 'analysis-chevron')",
          tags$i(id = "analysis-chevron", class = "fa fa-chevron-right", style = "margin-right: 8px; color: #333;"),
          tags$strong("Analysis choices")
        ),
        tags$div(
          id = "analysis-content",
          style = "padding: 15px;",
          # Perturbation-gene pairs to analyze
          selectInput("gene_list_mode", "Perturbation-gene pairs to analyze:",
                     choices = c("Random" = "random", "Custom" = "custom"),
                     selected = "random"),
          conditionalPanel(
            condition = "input.gene_list_mode == 'custom'",
            tags$div(
              style = "background-color: #fff3cd; border: 1px solid #ffeaa7; border-radius: 3px; padding: 6px; margin: 5px 0;",
              tags$small(
                tags$i(class = "fa fa-info-circle", style = "color: #856404; margin-right: 3px;"),
                tags$strong("Format: "), "CSV with 'grna_target' and 'response_id' columns (Ensembl gene IDs)",
                style = "color: #856404; font-size: 11px;"
              )
            ),
            tags$div(
              style = "margin-bottom: 15px;",
              fileInput("gene_list_file", 
                       label = NULL,
                       accept = c(".csv"),
                       placeholder = "Choose CSV file...")
            ),
            conditionalPanel(
              condition = "output.gene_list_uploaded",
              tags$div(
                style = "background-color: #d4edda; border: 1px solid #c3e6cb; border-radius: 4px; padding: 8px; margin: 0 0 15px 0;",
                tags$i(class = "fa fa-check-circle", style = "color: #155724; margin-right: 5px;"),
                textOutput("gene_list_status", inline = TRUE)
              )
            )
          ),
          
          # Minimum TPM threshold
          numericInput("tpm_threshold", "Minimum TPM threshold:", 10, 0, 10, 0.5),
          
          # Test side
          selectInput("side", "Test side:", 
                      choices = c("Left (knockdown)" = "left", 
                                 "Right (overexpression)" = "right"), 
                      selected = "left"),
          
          # Control group
          selectInput("control_group", "Control group:", 
                      choices = c("Complement cells" = "complement", 
                                 "Non-targeting cells" = "nt_cells"), 
                      selected = "complement"),
          
          # FDR target level
          numericInput("fdr_target", "FDR target level:", 0.05, 0.001, 0.1, 0.001)
        )
      ),
      
      br(),
      
      # Assumed effect sizes (collapsible)
      tags$div(
        style = "border: 1px solid #ddd; border-radius: 4px; margin-bottom: 10px;",
        tags$div(
          id = "effects-header",
          style = "background-color: #f5f5f5; padding: 10px 15px; cursor: pointer; border-bottom: 1px solid #ddd; color: #333;",
          onclick = "toggleSection('effects-content', 'effects-chevron')",
          tags$i(id = "effects-chevron", class = "fa fa-chevron-right", style = "margin-right: 8px; color: #333;"),
          tags$strong("Assumed effect sizes")
        ),
        tags$div(
          id = "effects-content",
          style = "padding: 15px;",
          numericInput("fc_mean", "Fold-change mean:", 0.85, 1.1, 10, 0.05),
          numericInput("fc_sd", "Fold-change SD:", 0.15, 0.1, 5, 0.05),
          numericInput("prop_non_null", "Proportion of non-null pairs:", 0.1, 0, 1, 0.01)
        )
      ),
      
      br(),
      
      tags$div(
        style = "text-align: center; padding: 0 20px;",
        actionButton("plan_btn", "Plan", class = "btn-success", style = "width: 200px; max-width: 90%;")
      )
    )
  )
}