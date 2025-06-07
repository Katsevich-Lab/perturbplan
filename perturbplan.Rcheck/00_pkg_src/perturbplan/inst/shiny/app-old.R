library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyBS)

ui <- dashboardPage(
  dashboardHeader(title = "PerturbPlan"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Experimental Design", tabName = "design", icon = icon("flask"))
    )
  ),

  dashboardBody(
    fluidRow(
      # Left column: Input parameters (wrapped in a parent box)
      box(
        width = 4,  # Half of the page width for inputs
        title = "Experimental Design Parameters", status = "primary", solidHeader = TRUE,

        box(
          id = "box1", title = "Biological System",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = FALSE,
          width = 12,

          # Radio buttons to ask if the user has pilot data
          radioButtons("has_pilot", "Do you have pilot data?", choices = c("Yes", "No"), inline = TRUE, selected = character(0)),

          # Conditional panel for uploading pilot data, but only once the user has made a choice
          conditionalPanel(
            condition = "input.has_pilot !== undefined && input.has_pilot === 'Yes'",
            fileInput("pilot_data", "Pilot Data (CSV)",
                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
          ),

          # Conditional panel for selecting a cell type, only once a choice is made
          conditionalPanel(
            condition = "input.has_pilot !== undefined && input.has_pilot === 'No'",
            selectInput("cell_type", "Cell Type",
                        choices = c("Immune cells", "Neurons", "Stem cells"))
          )
        ),

        box(
          id = "box2", title = "gRNA Design",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,  # Start collapsed
          width = 12,

          # Multiplicity of Infection input
          numericInput("moi", "Multiplicity of Infection", value = 1, min = 1, max = 100, step = 1),
          checkboxInput("help_moi", "Help me choose", value = FALSE),

          # Radio buttons: Have you already designed gRNAs?
          radioButtons("has_grna_design", "Have you already designed gRNAs?", choices = c("Yes", "No"), inline = TRUE, selected = character(0)),

          # Conditional: If "Yes" to gRNA design
          conditionalPanel(
            condition = "input.has_grna_design === 'Yes'",

            # File upload for designed gRNAs
            fileInput("grna_file", "Upload gRNA design file (CSV)",
                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),

            # Radio buttons: Do you have estimates of gRNA efficiencies?
            radioButtons("has_grna_efficiency", "Does this file include estimates of gRNA efficiencies?", choices = c("Yes", "No"), inline = TRUE, selected = character(0)),

            # Conditional: If "No" to gRNA efficiencies
            conditionalPanel(
              condition = "input.has_grna_efficiency === 'No'",
              numericInput("grna_eff", "gRNA efficiency (%)", value = 50, min = 0, max = 100, step = 1)
            )
          ),

          # Conditional: If "No" to gRNA design, show numeric inputs for gRNA design parameters
          conditionalPanel(
            condition = "input.has_grna_design === 'No'",

            numericInput("num_elements", "Number of Elements Targeted", value = 50, min = 1, max = 1e3, step = 1),
            checkboxInput("help_num_elements", "Help me choose", value = FALSE),

            numericInput("gRNAs_per_element", "Number of gRNAs per Element", value = 100, min = 10, max = 1e4, step = 10),
            checkboxInput("help_gRNAs_per_element", "Help me choose", value = FALSE),

            numericInput("nt_gRNAs", "Number of Non-Targeting gRNAs", value = 10, min = 1, max = 100, step = 1),
            checkboxInput("help_nt_gRNAs", "Help me choose", value = FALSE),

            numericInput("grna_eff", "Predicted average gRNA efficiency (%)", value = 50, min = 0, max = 100, step = 1)
          )
        ),

        # Box 3: Library prep
        box(
          id = "box3", title = "Library preparation",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,  # Start collapsed
          width = 12,
          numericInput("num_cells", "Number of Cells", value = 1000, min = 100, max = 1e6, step = 100),
          checkboxInput("help_num_cells", "Help me choose", value = FALSE),
          numericInput("num_cells", "Cost per cell ($)", value = 0.01, min = 0, max = 1, step = 0.01),
          numericInput("doublet_rate", "Doublet Rate (%)", value = 5, min = 0, max = 100, step = 0.1),
          # Radio button: Do you plan to directly capture transcripts?
          radioButtons("capture_transcripts", "Do you plan to directly capture transcripts?", choices = c("Yes", "No"), inline = TRUE, selected = character(0)),

          # Conditional panel: Only show file input if the answer is "Yes"
          conditionalPanel(
            condition = "input.capture_transcripts === 'Yes'",

            # Label for gene list upload with "i" icon immediately after the label
            div(
              strong("Upload Gene List (CSV)"),
              tags$i(icon("info-circle"), id = "gene_info"),
              bsPopover("gene_info", "Gene List Info", "Upload a list of genes for direct transcript capture.", placement = "right", trigger = "hover"),
              style = "display: inline-block;"
            ),

            # File upload input
            fileInput("gene_list", NULL, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
          )
        ),

        # Box 4: Sequencing
        box(
          id = "box4", title = "Sequencing",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,  # Start collapsed
          width = 12,
          numericInput("num_reads", "Number of Reads", value = 1e6, min = 1e4, max = 1e9, step = 1e4),
          checkboxInput("help_num_reads", "Help me choose", value = FALSE),
          numericInput("cost_per_read", "Cost per Read ($)", value = 0.002, min = 0.001, max = 1, step = 0.001),
          numericInput("mapping_efficiency", "Mapping Efficiency (%)", value = 90, min = 0, max = 100, step = 1)
        ),

        box(
          id = "box5", title = "Perturbation-Gene Testing",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,  # Start collapsed
          width = 12,

          # First, ask if they know which perturbation-gene pairs to test
          radioButtons("know_pairs", "Do you know which perturbation-gene pairs to test?", choices = c("Yes", "No"), inline = TRUE, selected = character(0)),

          # Conditional panel: If "Yes" to knowing perturbation-gene pairs
          conditionalPanel(
            condition = "input.know_pairs === 'Yes'",

            # Ask for file upload
            fileInput("pairs_file", "Upload Perturbation-Gene Pairs File (CSV)",
                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),

            # Ask if the file contains estimated effect sizes, unselected by default
            radioButtons("has_effect_size", "Does the above file contain estimated effect sizes?", choices = c("Yes", "No"), inline = TRUE, selected = character(0)),

            # Conditional: If "No" to effect sizes in the file, ask for average effect size
            conditionalPanel(
              condition = "input.has_effect_size === 'No'",
              numericInput("avg_effect_size", "Average Effect Size", value = 0.5, min = 0.01, max = 5, step = 0.01)
            )
          ),

          # Conditional panel: If "No" to knowing perturbation-gene pairs, ask for number of pairs to test, effect size, and minimum gene expression level in TPM
          conditionalPanel(
            condition = "input.know_pairs === 'No'",
            numericInput("pairs_to_test", "Number of Pairs to Test", value = 100, min = 1, max = 1000, step = 1),
            checkboxInput("help_pairs_to_test", "Help me choose", value = FALSE),

            numericInput("effect_size", "Effect Size", value = 0.5, min = 0.01, max = 5, step = 0.01),

            numericInput("min_gene_expr", "Minimum Gene Expression Level (TPM)", value = 10, min = 0, max = 1000, step = 1),
            checkboxInput("help_min_gene_expr", "Help me choose", value = FALSE)
          ),

          # Show Multiple Testing Correction Level only after both radio buttons have been selected
          conditionalPanel(
            condition = "input.know_pairs !== undefined && input.know_pairs !== '' &&
                 (input.know_pairs === 'No' || (input.has_effect_size !== undefined && input.has_effect_size !== ''))",
            numericInput("multiple_testing", "Multiple Testing Correction Level", value = 0.05, min = 0, max = 1, step = 0.01)
          )
        ),




        # Plan button at the bottom of all boxes
        actionButton("submit_button", "Plan", class = "btn-primary")
      ),

      # Right column: Outputs
      box(
        width = 8,  # Half of the page width for outputs
        title = "Cost and Power", status = "info", solidHeader = TRUE
      )
    )
  )
)

server <- function(input, output, session) {

  # Helper function to collapse only open boxes
  collapse_if_open <- function(box_id) {
    if (!input[[box_id]]$collapsed) {
      updateBox(box_id, action = "toggle", session = session)
    }
  }

  # Observe when box 1 is opened and close others if they are open
  observeEvent(input$box1$collapsed, {
    if (!input$box1$collapsed) {
      collapse_if_open("box2")
      collapse_if_open("box3")
      collapse_if_open("box4")
      collapse_if_open("box5")
    }
  })

  # Observe when box 2 is opened and close others if they are open
  observeEvent(input$box2$collapsed, {
    if (!input$box2$collapsed) {
      collapse_if_open("box1")
      collapse_if_open("box3")
      collapse_if_open("box4")
      collapse_if_open("box5")
    }
  })

  # Similarly for the rest of the boxes
  observeEvent(input$box3$collapsed, {
    if (!input$box3$collapsed) {
      collapse_if_open("box1")
      collapse_if_open("box2")
      collapse_if_open("box4")
      collapse_if_open("box5")
    }
  })

  observeEvent(input$box4$collapsed, {
    if (!input$box4$collapsed) {
      collapse_if_open("box1")
      collapse_if_open("box2")
      collapse_if_open("box3")
      collapse_if_open("box5")
    }
  })

  observeEvent(input$box5$collapsed, {
    if (!input$box5$collapsed) {
      collapse_if_open("box1")
      collapse_if_open("box2")
      collapse_if_open("box3")
      collapse_if_open("box4")
    }
  })
}

shinyApp(ui, server, options = list(height = 1080))
