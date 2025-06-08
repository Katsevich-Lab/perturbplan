# Server module for power calculations and reactive values

create_power_server <- function(input, output, session) {
  # Reactive values
  planned <- reactiveVal(FALSE)
  output$need_plan <- reactive(!planned())
  outputOptions(output, "need_plan", suspendWhenHidden = FALSE)
  observeEvent(input$plan_btn, planned(TRUE))

  # Real power calculation using perturbplan package
  power_results <- reactive({
    req(planned())
    # Call the package function with current input values
    perturbplan::calculate_power_grid(
      num_targets = input$num_targets,
      gRNAs_per_target = input$gRNAs_per_target, 
      non_targeting_gRNAs = input$non_targeting_gRNAs,
      tmp_threshold = input$tpm_threshold,
      fdr_target = input$fdr_target,
      fc_mean = input$fc_mean,
      fc_sd = input$fc_sd,
      prop_non_null = input$prop_non_null,
      MOI = input$MOI,
      biological_system = input$biological_system,
      experimental_platform = input$experimental_platform,
      side = input$side,
      control_group = input$control_group
    )
  })

  # Grid setup - get from power calculation results
  cells_seq <- reactive({
    req(planned())
    power_results()$cells_seq
  })
  
  reads_seq <- reactive({
    req(planned())
    power_results()$reads_seq
  })
  
  dcells <- reactive({
    seq_vals <- cells_seq()
    diff(seq_vals)[1]
  })
  
  dreads <- reactive({
    seq_vals <- reads_seq()
    diff(seq_vals)[1]
  })

  gridDF <- reactive({
    req(planned())
    power_results()$power_grid
  })

  # Return all reactive values and functions for use by other modules
  return(list(
    planned = planned,
    power_results = power_results,
    cells_seq = cells_seq,
    reads_seq = reads_seq,
    dcells = dcells,
    dreads = dreads,
    gridDF = gridDF
  ))
}