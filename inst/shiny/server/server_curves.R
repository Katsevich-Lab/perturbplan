# Server module for power curves and download functionality

create_curves_server <- function(input, output, session, power_data, selection_data) {
  
  # Show curves control box when curves are available
  output$curves_available <- reactive({
    power_data$planned() && selection_data$is_sel("tile") && nrow(selection_data$sel$tiles) > 0
  })
  outputOptions(output, "curves_available", suspendWhenHidden = FALSE)
  
  # Compute detailed power curves only when needed (for selected tiles)
  selected_power_curves <- reactive({
    req(power_data$planned(), power_data$fc_expression_info(), power_data$library_info(), selection_data$is_sel("tile"), nrow(selection_data$sel$tiles) > 0)
    
    # Create selected tiles data frame with pre-computed control cells
    # Get selected cell and read values
    selected_cells <- power_data$cells_seq()[selection_data$sel$tiles$row]
    selected_reads <- power_data$reads_seq()[selection_data$sel$tiles$col]
    
    # Debug: check if we have valid selections
    if (length(selected_cells) == 0 || length(selected_reads) == 0) {
      stop("No valid tile selections found")
    }
    
    # Find matching rows in power grid
    grid_df <- power_data$gridDF()
    
    # Use a more robust matching approach that handles potential rounding differences
    match_indices <- integer(length(selected_cells))
    for (i in seq_along(selected_cells)) {
      # Find exact matches for both cells and reads
      cell_matches <- which(abs(grid_df$cells - selected_cells[i]) < 1e-6)
      read_matches <- which(abs(grid_df$reads - selected_reads[i]) < 1e-6)
      common_matches <- intersect(cell_matches, read_matches)
      
      if (length(common_matches) == 0) {
        # If no exact match, find closest
        distances <- sqrt((grid_df$cells - selected_cells[i])^2 + (grid_df$reads - selected_reads[i])^2)
        match_indices[i] <- which.min(distances)
      } else {
        match_indices[i] <- common_matches[1]
      }
    }
    
    # Verify we found valid matches
    if (any(match_indices <= 0 | match_indices > nrow(grid_df))) {
      stop("Invalid match indices found")
    }
    
    selected_tiles <- data.frame(
      cells = selected_cells,
      reads = selected_reads,
      num_cntrl_cells = grid_df$num_cntrl_cells[match_indices]
    )
    
    # Debug: Print selected_tiles structure
    cat("DEBUG: selected_tiles structure:\n")
    print(str(selected_tiles))
    cat("DEBUG: selected_tiles data:\n")
    print(selected_tiles)
    cat("DEBUG: grid_df columns:\n")
    print(colnames(grid_df))
    cat("DEBUG: grid_df first few rows:\n")
    print(head(grid_df))
    
    # Compute power curves only for selected tiles using the new workflow
    perturbplan::calculate_power_curves(
      selected_tiles = selected_tiles,
      fc_expression_info = power_data$fc_expression_info(),  # Use extracted fc_expression_info
      library_info = power_data$library_info(),  # Use extracted library_info
      num_targets = input$num_targets,
      gRNAs_per_target = input$gRNAs_per_target,
      non_targeting_gRNAs = input$non_targeting_gRNAs,
      tpm_threshold = input$tpm_threshold,
      fdr_target = input$fdr_target,
      prop_non_null = input$prop_non_null,
      MOI = input$MOI,
      side = input$side,
      control_group = input$control_group
    )
  })

  # Helper function to extract power curve data for plots
  create_real_plot_data <- function(curve_type, curves_results) {
    power_curves <- curves_results$power_curves
    tiles_info <- curves_results$tiles_info
    
    # Ensure tiles_info columns are numeric to prevent factor/numeric mismatch
    tiles_info$cells <- as.numeric(as.character(tiles_info$cells))
    tiles_info$reads <- as.numeric(as.character(tiles_info$reads))
    
    dfs <- Map(function(i) {
      # Get curve data for this tile
      if (curve_type == "fc") {
        curve_data <- power_curves$fc_curves[[i]]
      } else {
        curve_data <- power_curves$expr_curves[[i]]
      }
      
      if (!is.null(curve_data) && nrow(curve_data) > 0) {
        # Extract values and ensure they're numeric
        cells_val <- round(as.numeric(tiles_info$cells[i]))
        reads_val <- as.numeric(tiles_info$reads[i])
        
        if (curve_type == "fc") {
          # Create data frame with explicit numeric columns
          df <- data.frame(
            x = as.numeric(curve_data$fold_change),
            Power = as.numeric(curve_data$power),
            label = as.character(sprintf("%d × %d", cells_val, reads_val)),
            cells = rep(cells_val, length(curve_data$fold_change)),
            reads = rep(reads_val, length(curve_data$fold_change)),
            stringsAsFactors = FALSE
          )
          # Ensure columns are the right type
          df$cells <- as.numeric(df$cells)
          df$reads <- as.numeric(df$reads)
          df
        } else {
          # Convert relative_expression to TPM scale
          tpm_values <- curve_data$relative_expression * 1e6
          df <- data.frame(
            x = as.numeric(tpm_values),
            Power = as.numeric(curve_data$power),
            label = as.character(sprintf("%d × %d", cells_val, reads_val)),
            cells = rep(cells_val, length(tpm_values)),
            reads = rep(reads_val, length(tpm_values)),
            stringsAsFactors = FALSE
          )
          # Ensure columns are the right type
          df$cells <- as.numeric(df$cells)
          df$reads <- as.numeric(df$reads)
          df
        }
      } else {
        data.frame(x = numeric(0), Power = numeric(0), label = character(0), 
                  cells = numeric(0), reads = numeric(0), stringsAsFactors = FALSE)
      }
    }, seq_len(nrow(tiles_info)))
    
    # Combine all valid data frames
    valid_dfs <- dfs[sapply(dfs, nrow) > 0]
    if (length(valid_dfs) > 0) {
      do.call(rbind, valid_dfs)
    } else {
      data.frame(x = numeric(0), Power = numeric(0), label = character(0), 
                cells = numeric(0), reads = numeric(0), stringsAsFactors = FALSE)
    }
  }

  # TPM (Expression) plot
  output$pp_tpm <- renderPlot({
    req(power_data$planned(), selection_data$is_sel("tile"))

    # Get power curves for selected tiles
    curves_results <- selected_power_curves()
    
    # Get baseline expression data for marginal distribution
    fc_expression_info <- power_data$fc_expression_info()
    baseline_expression <- fc_expression_info$fc_expression_df
    
    # Get power curve data for expression plot
    dfs1 <- create_real_plot_data("expr", curves_results)
    
    # Get baseline expression data for marginal distribution
    expr_data <- data.frame(
      tpm = baseline_expression$relative_expression * 1e6
    )
    
    # Create power curves with ggside marginal distributions
    library(ggside)
    
    # Expression power curve with marginal histogram
    if (nrow(dfs1) > 0) {
      display_mode <- if(is.null(input$curves_display_mode)) "all_together" else input$curves_display_mode
      
      if (display_mode == "all_together") {
        # Create factor variables for cells (color) and reads per cell (linetype and shape)
        dfs1$cells_factor <- factor(dfs1$cells,
                                   levels = sort(unique(dfs1$cells)),
                                   labels = paste0(sort(unique(dfs1$cells)), " treatment cells"))
        dfs1$reads_factor <- factor(dfs1$reads,
                                   levels = sort(unique(dfs1$reads)),
                                   labels = sort(unique(dfs1$reads)))
        
        # Create shape mapping (cycle through available shapes)
        unique_reads <- sort(unique(dfs1$reads))
        shape_values <- c(16, 17, 15, 18, 3, 4, 8, 7, 9, 10, 11, 12, 13, 14)[1:length(unique_reads)]
        names(shape_values) <- unique_reads
        
        # Create linetype mapping
        linetype_values <- rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), length.out = length(unique_reads))
        names(linetype_values) <- unique_reads
        
        ggplot(dfs1, aes(x, Power, colour = cells_factor, linetype = reads_factor, shape = reads_factor)) +
          geom_line(linewidth = 1) +
          geom_point(size = 2) +
          geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
          geom_xsidehistogram(data = expr_data, aes(x = tpm), 
                             bins = 60, fill = "lightblue", alpha = 0.7, 
                             inherit.aes = FALSE) +
          scale_x_log10() +
          scale_xsidey_continuous() +
          scale_shape_manual(values = shape_values) +
          scale_linetype_manual(values = linetype_values) +
          theme_bw(base_size = 16) +
          theme(aspect.ratio = 1,
                ggside.panel.scale = 0.3,
                legend.position = "right") +
          labs(x = "Expression Level (TPM)", y = "Power", 
               colour = "Number of treatment cells", linetype = "Reads per cell", shape = "Reads per cell")
        
      } else if (display_mode == "facet_cells") {
        # Facet by number of cells (horizontal panels)
        # Convert cells to factor for faceting to avoid ggplot2 internal issues
        dfs1$cells_factor <- factor(dfs1$cells, 
                                   levels = sort(unique(dfs1$cells)),
                                   labels = paste0("treatment cells: ", sort(unique(dfs1$cells))))
        # Create reads factor for coloring
        dfs1$reads_factor <- factor(dfs1$reads,
                                   levels = sort(unique(dfs1$reads)),
                                   labels = sort(unique(dfs1$reads)))
        
        ggplot(dfs1, aes(x, Power, colour = reads_factor)) +
          geom_line(linewidth = 1) +
          geom_point() +
          geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
          geom_xsidehistogram(data = expr_data, aes(x = tpm), 
                             bins = 60, fill = "lightblue", alpha = 0.7, 
                             inherit.aes = FALSE) +
          facet_grid(. ~ cells_factor) +
          scale_x_log10() +
          scale_xsidey_continuous() +
          theme_bw(base_size = 14) +
          theme(ggside.panel.scale = 0.2,
                legend.position = "bottom",
                strip.text = element_text(size = 12),
                aspect.ratio = 1) +
          labs(x = "Expression Level (TPM)", y = "Power", colour = "Reads per cell")
        
      } else if (display_mode == "facet_reads") {
        # Facet by reads per cell (horizontal panels)
        # Convert reads to factor for faceting to avoid ggplot2 internal issues
        dfs1$reads_factor <- factor(dfs1$reads,
                                   levels = sort(unique(dfs1$reads)),
                                   labels = paste0("reads: ", sort(unique(dfs1$reads))))
        # Create cells factor for coloring
        dfs1$cells_factor <- factor(dfs1$cells,
                                   levels = sort(unique(dfs1$cells)),
                                   labels = sort(unique(dfs1$cells)))
        
        ggplot(dfs1, aes(x, Power, colour = cells_factor)) +
          geom_line(linewidth = 1) +
          geom_point() +
          geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
          geom_xsidehistogram(data = expr_data, aes(x = tpm), 
                             bins = 60, fill = "lightblue", alpha = 0.7, 
                             inherit.aes = FALSE) +
          facet_grid(. ~ reads_factor) +
          scale_x_log10() +
          scale_xsidey_continuous() +
          theme_bw(base_size = 14) +
          theme(ggside.panel.scale = 0.2,
                legend.position = "bottom",
                strip.text = element_text(size = 12),
                aspect.ratio = 1) +
          labs(x = "Expression Level (TPM)", y = "Power", colour = "Number of treatment cells")
      }
    } else {
      ggplot() + theme_void()
    }
  }, height = 480)

  # Fold Change plot  
  output$pp_fc <- renderPlot({
    req(power_data$planned(), selection_data$is_sel("tile"))

    # Get power curves for selected tiles
    curves_results <- selected_power_curves()
    
    # Get power curve data for fold change plot
    dfs2 <- create_real_plot_data("fc", curves_results)
    
    # Create power curves with ggside marginal distributions
    library(ggside)

    # Fold-change power curve with marginal histogram
    if (nrow(dfs2) > 0) {
      # Use actual sampled fold changes from fc_expression_info
      fc_expression_df <- power_data$fc_expression_info()$fc_expression_df
      
      # Filter fold changes based on test side
      if (input$side == "left") {
        # Left-sided test: only show fold changes <= 1 (downregulation)
        filtered_fc <- fc_expression_df$fold_change[fc_expression_df$fold_change <= 1]
        x_limits <- c(min(dfs2$x, na.rm = TRUE), 1)
      } else if (input$side == "right") {
        # Right-sided test: only show fold changes >= 1 (upregulation)
        filtered_fc <- fc_expression_df$fold_change[fc_expression_df$fold_change >= 1]
        x_limits <- c(1, max(dfs2$x, na.rm = TRUE))
      } else {
        # Both-sided test: use same limits as either left or right based on fold_change_mean
        fold_change_mean <- power_data$fc_expression_info()$fold_change_mean
        if (fold_change_mean < 1) {
          # Use left-side limits (knockdown range)
          filtered_fc <- fc_expression_df$fold_change[fc_expression_df$fold_change <= 1]
          x_limits <- c(min(dfs2$x, na.rm = TRUE), 1)
        } else {
          # Use right-side limits (overexpression range)
          filtered_fc <- fc_expression_df$fold_change[fc_expression_df$fold_change >= 1]
          x_limits <- c(1, max(dfs2$x, na.rm = TRUE))
        }
      }
      
      fc_sample_data <- data.frame(fc = filtered_fc)
      display_mode <- if(is.null(input$curves_display_mode)) "all_together" else input$curves_display_mode
      
      base_plot <- ggplot(dfs2, aes(x, Power, colour = label)) +
        geom_line(linewidth = 1) +
        geom_point() +
        geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
        geom_xsidehistogram(data = fc_sample_data, aes(x = fc), 
                           bins = 60, fill = "pink", alpha = 0.7, inherit.aes = FALSE) +
        scale_x_continuous(limits = x_limits) +
        scale_xsidey_continuous(name = "Pair count") +
        labs(x = "Fold Change", y = "Power", colour = "Design (cells × reads/cell)")
      
      # Add vertical line at fold change = 1 for reference
      if (input$side %in% c("left", "right")) {
        base_plot <- base_plot + geom_vline(xintercept = 1, linetype = "dotted", colour = "darkgrey", alpha = 0.8)
      }
      
      if (display_mode == "all_together") {
        # Create factor variables for cells (color) and reads per cell (linetype and shape)
        dfs2$cells_factor <- factor(dfs2$cells,
                                   levels = sort(unique(dfs2$cells)),
                                   labels = paste0(sort(unique(dfs2$cells)), " treatment cells"))
        dfs2$reads_factor <- factor(dfs2$reads,
                                   levels = sort(unique(dfs2$reads)),
                                   labels = sort(unique(dfs2$reads)))
        
        # Create shape mapping (cycle through available shapes)
        unique_reads <- sort(unique(dfs2$reads))
        shape_values <- c(16, 17, 15, 18, 3, 4, 8, 7, 9, 10, 11, 12, 13, 14)[1:length(unique_reads)]
        names(shape_values) <- unique_reads
        
        # Create linetype mapping
        linetype_values <- rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), length.out = length(unique_reads))
        names(linetype_values) <- unique_reads
        
        ggplot(dfs2, aes(x, Power, colour = cells_factor, linetype = reads_factor, shape = reads_factor)) +
          geom_line(linewidth = 1) +
          geom_point(size = 2) +
          geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
          geom_xsidehistogram(data = fc_sample_data, aes(x = fc), 
                             bins = 60, fill = "pink", alpha = 0.7, inherit.aes = FALSE) +
          scale_x_continuous(limits = x_limits) +
          scale_xsidey_continuous(name = "Pair count") +
          scale_shape_manual(values = shape_values) +
          scale_linetype_manual(values = linetype_values) +
          labs(x = "Fold Change", y = "Power", 
               colour = "Number of treatment cells", linetype = "Reads per cell", shape = "Reads per cell") +
          {if (input$side %in% c("left", "right")) geom_vline(xintercept = 1, linetype = "dotted", colour = "darkgrey", alpha = 0.8)} +
          theme_bw(base_size = 16) +
          theme(aspect.ratio = 1,
                ggside.panel.scale = 0.3,
                legend.position = "right")
        
      } else if (display_mode == "facet_cells") {
        # Facet by number of cells (horizontal panels)
        # Convert cells to factor for faceting to avoid ggplot2 internal issues
        dfs2$cells_factor <- factor(dfs2$cells,
                                   levels = sort(unique(dfs2$cells)),
                                   labels = paste0("treatment cells: ", sort(unique(dfs2$cells))))
        # Create reads factor for coloring
        dfs2$reads_factor <- factor(dfs2$reads,
                                   levels = sort(unique(dfs2$reads)),
                                   labels = sort(unique(dfs2$reads)))
        
        base_plot <- ggplot(dfs2, aes(x, Power, colour = reads_factor)) +
          geom_line(linewidth = 1) +
          geom_point() +
          geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
          geom_xsidehistogram(data = fc_sample_data, aes(x = fc), 
                             bins = 60, fill = "pink", alpha = 0.7, inherit.aes = FALSE) +
          scale_x_continuous(limits = x_limits) +
          scale_xsidey_continuous() +
          labs(x = "Fold Change", y = "Power", colour = "Reads per cell")
        
        # Add vertical line at fold change = 1 for reference
        if (input$side %in% c("left", "right")) {
          base_plot <- base_plot + geom_vline(xintercept = 1, linetype = "dotted", colour = "darkgrey", alpha = 0.8)
        }
        
        base_plot +
          facet_grid(. ~ cells_factor) +
          theme_bw(base_size = 14) +
          theme(ggside.panel.scale = 0.2,
                legend.position = "bottom",
                strip.text = element_text(size = 12),
                aspect.ratio = 1)
        
      } else if (display_mode == "facet_reads") {
        # Facet by reads per cell (horizontal panels)
        # Convert reads to factor for faceting to avoid ggplot2 internal issues
        dfs2$reads_factor <- factor(dfs2$reads,
                                   levels = sort(unique(dfs2$reads)),
                                   labels = paste0("reads: ", sort(unique(dfs2$reads))))
        # Create cells factor for coloring
        dfs2$cells_factor <- factor(dfs2$cells,
                                   levels = sort(unique(dfs2$cells)),
                                   labels = sort(unique(dfs2$cells)))
                                   
        base_plot <- ggplot(dfs2, aes(x, Power, colour = cells_factor)) +
          geom_line(linewidth = 1) +
          geom_point() +
          geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
          geom_xsidehistogram(data = fc_sample_data, aes(x = fc), 
                             bins = 60, fill = "pink", alpha = 0.7, inherit.aes = FALSE) +
          scale_x_continuous(limits = x_limits) +
          scale_xsidey_continuous() +
          labs(x = "Fold Change", y = "Power", colour = "Number of treatment cells")
        
        # Add vertical line at fold change = 1 for reference
        if (input$side %in% c("left", "right")) {
          base_plot <- base_plot + geom_vline(xintercept = 1, linetype = "dotted", colour = "darkgrey", alpha = 0.8)
        }
        
        base_plot +
          facet_grid(. ~ reads_factor) +
          theme_bw(base_size = 14) +
          theme(ggside.panel.scale = 0.2,
                legend.position = "bottom",
                strip.text = element_text(size = 12),
                aspect.ratio = 1)
      }
    } else {
      ggplot() + theme_void()
    }
  }, height = 480)
  
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
      
      # 1. Add parameters sheet first (most important for understanding the analysis)
      # Include baseline expression information
      baseline_info <- if (!is.null(input$baseline_choice) && input$baseline_choice == "custom") {
        "Custom (uploaded RDS file)"
      } else {
        paste("Default (", input$biological_system, ")", sep = "")
      }
      
      library_info <- if (!is.null(input$library_choice) && input$library_choice == "custom") {
        "Custom (uploaded RDS file)"
      } else {
        paste("Default (", input$biological_system, ")", sep = "")
      }
      
      gene_list_info <- if (!is.null(input$gene_list_mode) && input$gene_list_mode == "custom") {
        "Custom (uploaded CSV file)"
      } else {
        "Random sampling"
      }
      
      params_df <- data.frame(
        Parameter = c(
          "Number of targets", "gRNAs per target", "Non-targeting gRNAs",
          "TPM threshold", "FDR target", "Fold-change mean", "Fold-change SD",
          "Proportion non-null", "MOI", "Biological system", "Experimental platform",
          "Test side", "Control group", "Baseline expression", "Library parameters", "Gene list mode"
        ),
        Value = c(
          input$num_targets, input$gRNAs_per_target, input$non_targeting_gRNAs,
          input$tpm_threshold, input$fdr_target, input$fc_mean, input$fc_sd,
          input$prop_non_null, input$MOI, input$biological_system, 
          input$experimental_platform, input$side, input$control_group,
          baseline_info, library_info, gene_list_info
        ),
        stringsAsFactors = FALSE
      )
      
      openxlsx::addWorksheet(wb, "1_Parameters")
      openxlsx::writeData(wb, "1_Parameters", params_df)
      
      # 2. Add power grid sheet (main results)
      # Rename columns to be more descriptive
      power_grid_clean <- power_data_results$power_grid
      colnames(power_grid_clean) <- c("Treatment_Cells", "Reads_per_Cell", "Power")
      
      openxlsx::addWorksheet(wb, "2_Power_Grid")
      openxlsx::writeData(wb, "2_Power_Grid", power_grid_clean)
      
      # 3. Add gene list if uploaded
      if (!is.null(power_data$gene_list()) && length(power_data$gene_list()) > 0) {
        gene_list_df <- data.frame(
          Gene_ID = power_data$gene_list(),
          stringsAsFactors = FALSE
        )
        openxlsx::addWorksheet(wb, "3_Gene_List")
        openxlsx::writeData(wb, "3_Gene_List", gene_list_df)
      }
      
      # 4. Add power curves if available (detailed drill-down results)
      if (selection_data$is_sel("tile") && nrow(selection_data$sel$tiles) > 0) {
        curves_data <- selected_power_curves()
        
        # Add tiles info first with logical column order
        tiles_info_clean <- data.frame(
          Design = paste0(round(curves_data$tiles_info$cells), " × ", curves_data$tiles_info$reads),
          Treatment_Cells = round(curves_data$tiles_info$cells),
          Reads_per_Cell = curves_data$tiles_info$reads,
          stringsAsFactors = FALSE
        )
        openxlsx::addWorksheet(wb, "4_Selected_Designs")
        openxlsx::writeData(wb, "4_Selected_Designs", tiles_info_clean)
        
        # 5. Add fold-change power curves with logical column order
        if (!is.null(curves_data$power_curves$fc_curves)) {
          fc_combined <- do.call(rbind, lapply(seq_along(curves_data$power_curves$fc_curves), function(i) {
            fc_data <- curves_data$power_curves$fc_curves[[i]]
            if (!is.null(fc_data) && nrow(fc_data) > 0) {
              # Reorder columns logically
              fc_clean <- data.frame(
                Design = paste0(round(curves_data$tiles_info$cells[i]), " × ", curves_data$tiles_info$reads[i]),
                Treatment_Cells = round(curves_data$tiles_info$cells[i]),
                Reads_per_Cell = curves_data$tiles_info$reads[i],
                Fold_Change = fc_data$fold_change,
                Power = fc_data$power,
                stringsAsFactors = FALSE
              )
              return(fc_clean)
            }
            return(NULL)
          }))
          
          if (!is.null(fc_combined) && nrow(fc_combined) > 0) {
            openxlsx::addWorksheet(wb, "5_Fold_Change_Power")
            openxlsx::writeData(wb, "5_Fold_Change_Power", fc_combined)
          }
        }
        
        # 6. Add expression power curves with logical column order
        if (!is.null(curves_data$power_curves$expr_curves)) {
          expr_combined <- do.call(rbind, lapply(seq_along(curves_data$power_curves$expr_curves), function(i) {
            expr_data <- curves_data$power_curves$expr_curves[[i]]
            if (!is.null(expr_data) && nrow(expr_data) > 0) {
              # Reorder columns logically
              expr_clean <- data.frame(
                Design = paste0(round(curves_data$tiles_info$cells[i]), " × ", curves_data$tiles_info$reads[i]),
                Treatment_Cells = round(curves_data$tiles_info$cells[i]),
                Reads_per_Cell = curves_data$tiles_info$reads[i],
                Expression_TPM = expr_data$relative_expression * 1e6,
                Power = expr_data$power,
                stringsAsFactors = FALSE
              )
              return(expr_clean)
            }
            return(NULL)
          }))
          
          if (!is.null(expr_combined) && nrow(expr_combined) > 0) {
            openxlsx::addWorksheet(wb, "6_Expression_Power")
            openxlsx::writeData(wb, "6_Expression_Power", expr_combined)
          }
        }
      }
      
      # Save workbook
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}