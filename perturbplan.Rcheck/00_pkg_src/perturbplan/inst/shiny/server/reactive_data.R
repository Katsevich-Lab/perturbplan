# Reactive data and power calculations
# Handles the interface between UI inputs and power calculation functions

# Main power calculation - calls the package function
power_results <- reactive({
  req(planned())
  
  # Call the package function with current input values
  perturbplan::calculate_power_grid(
    num_targets = input$num_targets,
    gRNAs_per_target = input$gRNAs_per_target, 
    non_targeting_gRNAs = input$non_targeting_gRNAs,
    num_pairs = input$num_pairs,
    tpm_threshold = input$tpm_threshold,
    fdr_target = input$fdr_target,
    fc_mean = input$fc_mean,
    fc_sd = input$fc_sd,
    prop_non_null = input$prop_non_null,
    MOI = input$MOI,
    biological_system = input$biological_system,
    experimental_platform = input$experimental_platform
  )
})

# Grid data for plots
gridDF <- reactive({
  req(power_results())
  power_results()$power_grid
})

# Sequences for grid
cells_seq <- reactive({
  req(power_results())
  power_results()$cells_seq
})

reads_seq <- reactive({
  req(power_results())
  power_results()$reads_seq
})

# Grid spacing for plot rectangles
dcells <- reactive({
  seq_vals <- cells_seq()
  diff(seq_vals)[1]
})

dreads <- reactive({
  seq_vals <- reads_seq()
  diff(seq_vals)[1]
})