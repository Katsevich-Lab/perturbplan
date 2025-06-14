# Sidebar UI module for parameter inputs

create_sidebar <- function() {
  dashboardSidebar(
    # Parameter panels - make scrollable with collapsible sections
    tags$div(
      style = "padding: 10px; max-height: 90vh; overflow-y: auto;",
      
      # Experimental choices (collapsible)
      tags$div(
        style = "border-radius: 4px; margin-bottom: 5px;",
        tags$div(
          id = "exp-header",
          style = "padding: 10px 15px; cursor: pointer; border-radius: 4px 4px 0 0;",
          onclick = "toggleSection('experimental-content', 'exp-chevron')",
          tags$i(id = "exp-chevron", class = "fa fa-chevron-down", style = "margin-right: 8px;"),
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
      
      # Analysis choices (collapsible)
      tags$div(
        style = "border-radius: 4px; margin-bottom: 5px;",
        tags$div(
          id = "analysis-header",
          style = "padding: 10px 15px; cursor: pointer; border-radius: 4px 4px 0 0;",
          onclick = "toggleSection('analysis-content', 'analysis-chevron')",
          tags$i(id = "analysis-chevron", class = "fa fa-chevron-right", style = "margin-right: 8px;"),
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
              class = "file-upload-info",
              style = "border-radius: 3px; padding: 6px; margin: 5px 0;",
              tags$small(
                tags$i(class = "fa fa-info-circle", style = "margin-right: 3px;"),
                tags$strong("Format: "), "CSV with 'grna_target' and 'response_id' columns (Ensembl gene IDs)",
                style = "font-size: 11px;"
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
                class = "file-upload-success",
                style = "border-radius: 4px; padding: 8px; margin: 0 0 15px 0;",
                tags$i(class = "fa fa-check-circle", style = "margin-right: 5px;"),
                htmlOutput("gene_list_status", inline = TRUE)
              )
            )
          ),
          
          # Minimum TPM threshold
          numericInput("tpm_threshold", "Minimum TPM threshold:", 10, 0, 10, 0.5),
          
          # Test side
          selectInput("side", "Test side:", 
                      choices = c("Left (knockdown)" = "left", 
                                 "Right (overexpression)" = "right",
                                 "Both (two-sided)" = "both"), 
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
      
      # Assumed effect sizes (collapsible)
      tags$div(
        style = "border-radius: 4px; margin-bottom: 5px;",
        tags$div(
          id = "effects-header",
          style = "padding: 10px 15px; cursor: pointer; border-radius: 4px 4px 0 0;",
          onclick = "toggleSection('effects-content', 'effects-chevron')",
          tags$i(id = "effects-chevron", class = "fa fa-chevron-right", style = "margin-right: 8px;"),
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
      
      # Pilot data choice (collapsible)
      tags$div(
        style = "border-radius: 4px; margin-bottom: 5px;",
        tags$div(
          id = "pilot-header",
          style = "padding: 10px 15px; cursor: pointer; border-radius: 4px 4px 0 0;",
          onclick = "toggleSection('pilot-content', 'pilot-chevron')",
          tags$i(id = "pilot-chevron", class = "fa fa-chevron-right", style = "margin-right: 8px;"),
          tags$strong("Pilot data choice")
        ),
        tags$div(
          id = "pilot-content",
          style = "padding: 15px;",
          
          # Baseline expression choice
          selectInput("baseline_choice", "Baseline expression:",
                     choices = c("Default" = "default", "Custom" = "custom"),
                     selected = "default"),
          
          # Conditional panel for custom baseline upload
          conditionalPanel(
            condition = "input.baseline_choice == 'custom'",
            tags$div(
              class = "file-upload-info",
              style = "border-radius: 3px; padding: 6px; margin: 5px 0;",
              tags$small(
                tags$i(class = "fa fa-info-circle", style = "margin-right: 3px;"),
                tags$strong("Format: "), "RDS file with same structure as extract_baseline_expression() output",
                style = "font-size: 11px;"
              )
            ),
            tags$div(
              style = "margin-bottom: 15px;",
              fileInput("baseline_file", 
                       label = NULL,
                       accept = c(".rds", ".RDS"),
                       placeholder = "Choose RDS file...")
            ),
            conditionalPanel(
              condition = "output.baseline_uploaded",
              tags$div(
                class = "file-upload-success",
                style = "border-radius: 4px; padding: 8px; margin: 0 0 15px 0;",
                tags$i(class = "fa fa-check-circle", style = "margin-right: 5px;"),
                htmlOutput("baseline_status", inline = TRUE)
              )
            )
          )
        ),
        
        # Library size parameters (collapsible)
        tags$div(
          style = "border-radius: 4px; margin-bottom: 5px;",
          tags$div(
            id = "library-header",
            style = "padding: 10px 15px; cursor: pointer; border-radius: 4px 4px 0 0;",
            onclick = "toggleSection('library-content', 'library-chevron')",
            tags$i(id = "library-chevron", class = "fa fa-chevron-right", style = "margin-right: 8px;"),
            tags$strong("Library size parameters")
          ),
          tags$div(
            id = "library-content",
            style = "padding: 15px; display: none;",
            selectInput("library_choice", "Library parameters:",
                       choices = c("Default" = "default", "Custom" = "custom"),
                       selected = "default"),
            conditionalPanel(
              condition = "input.library_choice == 'custom'",
              tags$div(
                style = "font-size: 11px; color: #666; margin-bottom: 8px;",
                "Upload RDS file with UMI_per_cell and variation parameters"
              )
            ),
            conditionalPanel(
              condition = "input.library_choice == 'custom'",
              tags$div(
                style = "margin-bottom: 15px;",
                fileInput("library_file", 
                         label = NULL,
                         accept = c(".rds", ".RDS"),
                         placeholder = "Choose RDS file...")
              ),
              conditionalPanel(
                condition = "output.library_uploaded",
                tags$div(
                  class = "file-upload-success",
                  style = "border-radius: 4px; padding: 8px; margin: 0 0 15px 0;",
                  tags$i(class = "fa fa-check-circle", style = "margin-right: 5px;"),
                  htmlOutput("library_status", inline = TRUE)
                )
              )
            )
          )
        )
      ),
      
      # Horizontal separator line
      tags$hr(class = "sidebar-separator"),
      
      tags$div(
        style = "text-align: center; padding: 0 20px;",
        actionButton("plan_btn", "Plan", class = "btn-success", style = "width: 200px; max-width: 90%;")
      )
    )
  )
}