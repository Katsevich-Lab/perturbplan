# Plot rendering functions
# Handles all plot generation and visualization

# Helper function for drawing selection rectangles
draw_rects <- function(df,p)
  p + geom_rect(data=df,
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                inherit.aes=FALSE,colour="red",linewidth=1,fill=NA)

# Main heatmap
output$heat <- renderPlot({
  req(planned())
  p <- ggplot(gridDF(),aes(reads,cells,fill=power))+
    geom_tile(colour=NA)+
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw(base_size = 16)+
    theme(panel.grid=element_blank(),
          aspect.ratio = 1)+
    labs(x="Reads per cell",y="Number of cells",fill="Power")
  
  if (is_sel("row") && length(sel$idx)) {
    rect <- data.frame(
      xmin=min(reads_seq())-dreads()/2,
      xmax=max(reads_seq())+dreads()/2,
      ymin=cells_seq()[sel$idx]-dcells()/2,
      ymax=cells_seq()[sel$idx]+dcells()/2
    )
    p <- draw_rects(rect,p)
  } else if (is_sel("col") && length(sel$idx)) {
    rect <- data.frame(
      xmin=reads_seq()[sel$idx]-dreads()/2,
      xmax=reads_seq()[sel$idx]+dreads()/2,
      ymin=min(cells_seq())-dcells()/2,
      ymax=max(cells_seq())+dcells()/2
    )
    p <- draw_rects(rect,p)
  } else if (is_sel("tile") && nrow(sel$tiles)) {
    rect <- data.frame(
      xmin=reads_seq()[sel$tiles$col]-dreads()/2,
      xmax=reads_seq()[sel$tiles$col]+dreads()/2,
      ymin=cells_seq()[sel$tiles$row]-dcells()/2,
      ymax=cells_seq()[sel$tiles$row]+dcells()/2
    )
    p <- draw_rects(rect,p)
  }
  p
})

# Slice title & plot
output$slice_title <- renderText({
  if (identical(slice_mode(),"row"))
    sprintf("Power versus reads per cell")
  else if (identical(slice_mode(),"col"))
    sprintf("Power vs number of cells")
  else "Slice view"
})

output$slice_plot <- renderPlot({
  req(planned(), !is.null(slice_mode()), length(sel$idx))
  df <- gridDF()
  if (slice_mode()=="row") {
    sub <- subset(df, cells %in% cells_seq()[sel$idx])
    p <- ggplot(sub,aes(reads,power,colour=factor(cells)))+
      geom_line()+
      geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
      theme_bw(base_size = 16)+
      theme(aspect.ratio = 1) +
      labs(x="Reads per cell",y="Power",colour="Cells")
    
    # Add vertical lines for all selected points
    if (length(slice_x()) > 0) {
      p <- p + geom_vline(xintercept=slice_x(),colour="red")
    }
    p
  } else {
    sub <- subset(df, reads %in% reads_seq()[sel$idx])
    p <- ggplot(sub,aes(cells,power,colour=factor(reads)))+
      geom_line()+
      geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
      theme_bw(base_size = 16)+
      theme(aspect.ratio = 1) +
      labs(x="Number of cells",y="Power",colour="Reads")
    
    # Add vertical lines for all selected points
    if (length(slice_x()) > 0) {
      p <- p + geom_vline(xintercept=slice_x(),colour="red")
    }
    p
  }
})

# Per-pair title & plots
output$pair_title <- renderText({
  "Per-pair power"
})

output$pp_combined <- renderPlot({
  req(planned(), is_sel("tile"))
  
  # Get real power curves from calculate_power_grid results
  power_curves <- power_results()$power_curves
  
  # Extract power curve data for selected tiles  
  create_real_plot_data <- function(curve_type) {
    # Get the full power results from calculate_power_grid
    full_power_results <- power_results()
    
    # Map selected tiles to corresponding rows in the power computation results
    dfs <- Map(function(r, c) {
      # Get the actual cell and read values for this tile
      cell_val <- cells_seq()[r]
      read_val <- reads_seq()[c]
      
      # Find matching row in the grid - need to match against the original computation
      # The power_results() comes from calculate_power_grid which calls compute_power_grid_efficient
      power_grid <- full_power_results$power_grid
      matching_row <- which(power_grid$cells == cell_val & power_grid$reads == read_val)
      
      if (length(matching_row) > 0) {
        # Get curve data from the power_curves section
        if (curve_type == "fc") {
          curve_data <- full_power_results$power_curves$fc_curves[[matching_row[1]]]
        } else {
          curve_data <- full_power_results$power_curves$expr_curves[[matching_row[1]]]
        }
        
        if (!is.null(curve_data) && nrow(curve_data) > 0) {
          # The curve_data is already a data frame with proper column names
          if (curve_type == "fc") {
            # For fold change: scale to percentage and use meaningful labels
            data.frame(
              x = curve_data$fold_change,  # Already in fold change units
              Power = curve_data$power,     # Already in power units (0-1)
              design = factor(sprintf("%d × %d", cell_val, read_val))  # Combined cells-reads label
            )
          } else {
            # For expression: convert relative expression to TPM scale
            # The relative_expression is in relative units, convert to TPM
            tpm_values <- curve_data$relative_expression * 1e6  # Convert to TPM scale
            data.frame(
              x = tpm_values,              # Expression in TPM units  
              Power = curve_data$power,     # Power values (0-1)
              design = factor(sprintf("%d × %d", cell_val, read_val))  # Combined cells-reads label
            )
          }
        } else {
          # Return empty data frame if no curve data
          data.frame(
            x = numeric(0),
            Power = numeric(0),
            design = factor(character(0))
          )
        }
      } else {
        # Return empty data frame if no matching grid point
        data.frame(
          x = numeric(0), 
          Power = numeric(0),
          design = factor(character(0))
        )
      }
    }, sel$tiles$row, sel$tiles$col)
    
    # Combine all valid data frames
    valid_dfs <- dfs[sapply(dfs, nrow) > 0]
    if (length(valid_dfs) > 0) {
      do.call(rbind, valid_dfs)
    } else {
      # Return empty data frame with correct structure
      data.frame(
        x = numeric(0),
        Power = numeric(0), 
        design = factor(character(0))
      )
    }
  }

  # First plot: Power vs Expression Level (TPM scale)
  df1 <- create_real_plot_data("expr")
  
  p1 <- ggplot(df1, aes(x, Power, colour = design)) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
    scale_x_log10() +  # Log scale for TPM values
    theme_bw(base_size = 16) +
    theme(aspect.ratio = 1) +
    labs(x = "Expression Level (TPM)", y = "Power", colour = "Design (cells × reads/cell)")

  # Second plot: Power vs Fold-change
  df2 <- create_real_plot_data("fc")
  
  p2 <- ggplot(df2, aes(x, Power, colour = design)) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
    theme_bw(base_size = 16) +
    theme(aspect.ratio = 1) +
    labs(x = "Fold Change", y = "Power", colour = "Design (cells × reads/cell)")

  # Combine with shared legend below
  (p1 + p2) +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
})