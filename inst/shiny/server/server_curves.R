# Server module for power curves and download functionality

create_curves_server <- function(input, output, session, power_data, selection_data) {
  
  # Compute detailed power curves only when needed (for selected tiles)
  selected_power_curves <- reactive({
    req(power_data$planned(), power_data$fc_expression_info(), selection_data$is_sel("tile"), nrow(selection_data$sel$tiles) > 0)
    
    # Create selected tiles data frame
    selected_tiles <- data.frame(
      cells = power_data$cells_seq()[selection_data$sel$tiles$row],
      reads = power_data$reads_seq()[selection_data$sel$tiles$col]
    )
    
    # Compute power curves only for selected tiles using the new workflow
    perturbplan::calculate_power_curves(
      selected_tiles = selected_tiles,
      fc_expression_info = power_data$fc_expression_info(),  # Use extracted fc_expression_info
      num_targets = input$num_targets,
      gRNAs_per_target = input$gRNAs_per_target,
      non_targeting_gRNAs = input$non_targeting_gRNAs,
      tpm_threshold = input$tpm_threshold,
      fdr_target = input$fdr_target,
      prop_non_null = input$prop_non_null,
      MOI = input$MOI,
      experimental_platform = input$experimental_platform,
      side = input$side,
      control_group = input$control_group
    )
  })

  # Per-pair plots with real power curves
  output$pp_combined <- renderPlot({
    req(power_data$planned(), selection_data$is_sel("tile"))

    # Get power curves for selected tiles
    curves_results <- selected_power_curves()
    tiles_info <- curves_results$tiles_info
    
    # Extract power curve data for selected tiles  
    create_real_plot_data <- function(curve_type) {
      power_curves <- curves_results$power_curves
      
      dfs <- Map(function(i) {
        # Get curve data for this tile
        if (curve_type == "fc") {
          curve_data <- power_curves$fc_curves[[i]]
        } else {
          curve_data <- power_curves$expr_curves[[i]]
        }
        
        if (!is.null(curve_data) && nrow(curve_data) > 0) {
          if (curve_type == "fc") {
            data.frame(
              x = curve_data$fold_change,
              Power = curve_data$power,
              label = sprintf("%d × %d", tiles_info$cells[i], tiles_info$reads[i])
            )
          } else {
            # Convert relative_expression to TPM scale
            tpm_values <- curve_data$relative_expression * 1e6
            data.frame(
              x = tpm_values,
              Power = curve_data$power,
              label = sprintf("%d × %d", tiles_info$cells[i], tiles_info$reads[i])
            )
          }
        } else {
          data.frame(x = numeric(0), Power = numeric(0), label = character(0))
        }
      }, seq_len(nrow(tiles_info)))
      
      # Combine all valid data frames
      valid_dfs <- dfs[sapply(dfs, nrow) > 0]
      if (length(valid_dfs) > 0) {
        do.call(rbind, valid_dfs)
      } else {
        data.frame(x = numeric(0), Power = numeric(0), label = character(0))
      }
    }

    # First plot: Power vs TPM (real data)
    dfs1 <- create_real_plot_data("expr")
    p1 <- ggplot(dfs1, aes(x, Power, colour = label)) +
      geom_line() +
      geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
      scale_x_log10() +  # Log scale for TPM values
      theme_bw(base_size = 16) +
      theme(aspect.ratio = 1) +
      labs(x = "Expression Level (TPM)", y = "Power", colour = "Design (cells × reads/cell)")

    # Second plot: Power vs Fold‐change (real data)
    dfs2 <- create_real_plot_data("fc")
    p2 <- ggplot(dfs2, aes(x, Power, colour = label)) +
      geom_line() +
      geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
      theme_bw(base_size = 16) +
      theme(aspect.ratio = 1) +
      labs(x = "Fold Change",
           y = "Power",
           colour = "Design (cells × reads/cell)")

    # Combine with shared legend below
    (p1 + p2) +
      plot_layout(ncol = 2, guides = "collect") &
      theme(legend.position = "bottom")
  })
  
  # Download handler for results
  output$download_results <- downloadHandler(
    filename = function() {
      paste("perturbplan_results_", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      req(power_data$planned())
      
      # Get power results
      power_data_results <- power_data$power_results()
      
      # Create a workbook
      wb <- openxlsx::createWorkbook()
      
      # Add power grid sheet
      openxlsx::addWorksheet(wb, "Power_Grid")
      openxlsx::writeData(wb, "Power_Grid", power_data_results$power_grid)
      
      # Add parameters sheet
      params_df <- data.frame(
        Parameter = c(
          "Number of targets", "gRNAs per target", "Non-targeting gRNAs",
          "TPM threshold", "FDR target", "Fold-change mean", "Fold-change SD",
          "Proportion non-null", "MOI", "Biological system", "Experimental platform",
          "Test side", "Control group"
        ),
        Value = c(
          input$num_targets, input$gRNAs_per_target, input$non_targeting_gRNAs,
          input$tpm_threshold, input$fdr_target, input$fc_mean, input$fc_sd,
          input$prop_non_null, input$MOI, input$biological_system, 
          input$experimental_platform, input$side, input$control_group
        ),
        stringsAsFactors = FALSE
      )
      
      openxlsx::addWorksheet(wb, "Parameters")
      openxlsx::writeData(wb, "Parameters", params_df)
      
      # Add gene list if uploaded
      if (!is.null(power_data$gene_list()) && length(power_data$gene_list()) > 0) {
        gene_list_df <- data.frame(
          Gene_Name = power_data$gene_list(),
          stringsAsFactors = FALSE
        )
        openxlsx::addWorksheet(wb, "Gene_List")
        openxlsx::writeData(wb, "Gene_List", gene_list_df)
      }
      
      # Add power curves if available
      if (selection_data$is_sel("tile") && nrow(selection_data$sel$tiles) > 0) {
        curves_data <- selected_power_curves()
        
        # Add tiles info
        openxlsx::addWorksheet(wb, "Selected_Tiles")
        openxlsx::writeData(wb, "Selected_Tiles", curves_data$tiles_info)
        
        # Add fold-change curves
        if (!is.null(curves_data$power_curves$fc_curves)) {
          fc_combined <- do.call(rbind, lapply(seq_along(curves_data$power_curves$fc_curves), function(i) {
            fc_data <- curves_data$power_curves$fc_curves[[i]]
            if (!is.null(fc_data) && nrow(fc_data) > 0) {
              fc_data$tile_index <- i
              fc_data$cells <- curves_data$tiles_info$cells[i]
              fc_data$reads <- curves_data$tiles_info$reads[i]
              return(fc_data)
            }
            return(NULL)
          }))
          
          if (!is.null(fc_combined) && nrow(fc_combined) > 0) {
            openxlsx::addWorksheet(wb, "FC_Power_Curves")
            openxlsx::writeData(wb, "FC_Power_Curves", fc_combined)
          }
        }
        
        # Add expression curves
        if (!is.null(curves_data$power_curves$expr_curves)) {
          expr_combined <- do.call(rbind, lapply(seq_along(curves_data$power_curves$expr_curves), function(i) {
            expr_data <- curves_data$power_curves$expr_curves[[i]]
            if (!is.null(expr_data) && nrow(expr_data) > 0) {
              expr_data$tile_index <- i
              expr_data$cells <- curves_data$tiles_info$cells[i]
              expr_data$reads <- curves_data$tiles_info$reads[i]
              # Convert to TPM scale
              expr_data$TPM <- expr_data$relative_expression * 1e6
              return(expr_data)
            }
            return(NULL)
          }))
          
          if (!is.null(expr_combined) && nrow(expr_combined) > 0) {
            openxlsx::addWorksheet(wb, "TPM_Power_Curves")
            openxlsx::writeData(wb, "TPM_Power_Curves", expr_combined)
          }
        }
      }
      
      # Save workbook
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}