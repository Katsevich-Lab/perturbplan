# UI definition with tabs moved to main content area
# Updated with new parameters and features while maintaining modular structure

ui <- dashboardPage(
  # Dashboard header
  dashboardHeader(
    title = "PerturbPlan",
    tags$li(
      class = "dropdown",
      style = "margin: 8px 10px 0 0;",
      conditionalPanel(
        condition = "!output.need_plan",
        downloadButton(
          "download_results",
          "Download Results",
          class = "btn-primary",
          style = "background-color: white; color: #3c8dbc; border: 1px solid #3c8dbc; margin-top: 7px;"
        )
      )
    )
  ),
  
  # Dashboard sidebar - only parameters
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
          style = "padding: 15px; display: block;",
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
          style = "padding: 15px; display: none;",
          numericInput("tpm_threshold", "Minimum TPM threshold:", 10, 0, 10, 0.5),
          numericInput("fdr_target", "FDR target level:", 0.05, 0.001, 0.1, 0.001),
          selectInput("side", "Test side:", 
                      choices = c("Left (knockdown)" = "left", 
                                 "Right (overexpression)" = "right"), 
                      selected = "left"),
          selectInput("control_group", "Control group:", 
                      choices = c("Complement cells" = "complement", 
                                 "Non-targeting cells" = "nt_cells"), 
                      selected = "complement")
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
          style = "padding: 15px; display: none;",
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
  ),
  
  # Dashboard body with tabs
  dashboardBody(
    useShinyjs(),
    
    # Add custom CSS to adjust sidebar width and color
    tags$style(HTML("
      .main-sidebar {
        width: 350px !important;
        background-color: #3c8dbc !important;
      }
      
      .main-sidebar .sidebar {
        background-color: #3c8dbc !important;
      }
      
      .content-wrapper, .right-side {
        margin-left: 350px !important;
      }
      
      @media (max-width: 767px) {
        .content-wrapper, .right-side {
          margin-left: 0 !important;
        }
      }
    ")),
    
    # Add custom JavaScript for collapsible sections with auto-collapse
    tags$script(HTML("
      function toggleSection(contentId, chevronId) {
        var content = document.getElementById(contentId);
        var chevron = document.getElementById(chevronId);
        
        if (content.style.display === 'none') {
          // Close other sections first
          collapseOtherSections(contentId);
          
          content.style.display = 'block';
          chevron.className = 'fa fa-chevron-down';
        } else {
          content.style.display = 'none';
          chevron.className = 'fa fa-chevron-right';
        }
      }
      
      function collapseOtherSections(currentContentId) {
        var sections = ['experimental-content', 'analysis-content', 'effects-content'];
        var chevrons = ['exp-chevron', 'analysis-chevron', 'effects-chevron'];
        
        for (var i = 0; i < sections.length; i++) {
          if (sections[i] !== currentContentId) {
            var content = document.getElementById(sections[i]);
            var chevron = document.getElementById(chevrons[i]);
            
            if (content && content.style.display !== 'none') {
              content.style.display = 'none';
              chevron.className = 'fa fa-chevron-right';
            }
          }
        }
      }
      
      function checkSectionCompletion(sectionId) {
        var section = document.getElementById(sectionId);
        if (!section) return false;
        
        var inputs = section.querySelectorAll('input, select');
        var allFilled = true;
        
        for (var i = 0; i < inputs.length; i++) {
          var input = inputs[i];
          if (input.type === 'number' && (input.value === '' || input.value === null)) {
            allFilled = false;
            break;
          } else if (input.tagName === 'SELECT' && input.value === '') {
            allFilled = false;
            break;
          }
        }
        
        return allFilled;
      }
      
      function autoCollapseIfComplete() {
        setTimeout(function() {
          var sections = [
            {content: 'experimental-content', chevron: 'exp-chevron'},
            {content: 'analysis-content', chevron: 'analysis-chevron'},
            {content: 'effects-content', chevron: 'effects-chevron'}
          ];
          
          for (var i = 0; i < sections.length; i++) {
            var section = sections[i];
            var content = document.getElementById(section.content);
            var chevron = document.getElementById(section.chevron);
            
            if (content && content.style.display !== 'none' && 
                checkSectionCompletion(section.content)) {
              
              // Find next incomplete section to open
              var nextSectionIndex = (i + 1) % sections.length;
              var nextSection = sections[nextSectionIndex];
              
              // If next section is incomplete, switch to it
              if (!checkSectionCompletion(nextSection.content)) {
                content.style.display = 'none';
                chevron.className = 'fa fa-chevron-right';
                
                var nextContent = document.getElementById(nextSection.content);
                var nextChevron = document.getElementById(nextSection.chevron);
                if (nextContent) {
                  nextContent.style.display = 'block';
                  nextChevron.className = 'fa fa-chevron-down';
                }
                break;
              }
            }
          }
        }, 500);
      }
      
      // Monitor input changes
      $(document).ready(function() {
        $('body').on('change', 'input, select', function() {
          autoCollapseIfComplete();
        });
      });
    ")),
    
    # Tab panel in main content area
    tabsetPanel(
      id = "main_tabs",
      
      # Overall power heatmap tab
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
          
          box(
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
      ),
      
      # Power over TPM & FC tab
      tabPanel(
        "Power over TPM & FC",
        value = "per_pair",
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