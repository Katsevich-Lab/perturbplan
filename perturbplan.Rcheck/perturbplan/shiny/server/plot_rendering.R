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

make_curve <- function(x, base, fun, label)
  data.frame(x=x, Power=fun(x, base), label=label)

output$pp_combined <- renderPlot({
  req(planned(), is_sel("tile"))
  
  # Get max power for scaling (from placeholder function)
  max_pow <- max(log10(cells_seq()) * log10(reads_seq()))

  # Create data with consistent color/linetype mapping
  # Color = cell number, Linetype = reads per cell
  create_plot_data <- function(x_vals, power_fun) {
    dfs <- Map(function(r, c){
      base <- log10(cells_seq()[r] * reads_seq()[c])
      data.frame(
        x = x_vals,
        Power = power_fun(x_vals, base),
        cells = factor(cells_seq()[r]),
        reads = factor(reads_seq()[c])
      )
    }, sel$tiles$row, sel$tiles$col)
    do.call(rbind, dfs)
  }

  # First plot: Power vs TPM
  tpm  <- seq(0, 10, length.out = 100)
  df1 <- create_plot_data(tpm, function(v, b) 1 - exp(-((v/5)^1.8) * (b/max_pow)^10 * 5000))
  
  p1 <- ggplot(df1, aes(x, Power, colour = cells, linetype = reads)) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
    theme_bw(base_size = 16) +
    theme(aspect.ratio = 1) +
    labs(x = "TPM", y = "Power", colour = "Cells", linetype = "Reads/cell")

  # Second plot: Power vs Fold-change
  fc   <- seq(0, 50, length.out = 100)
  df2 <- create_plot_data(fc, function(v, b) 1 - exp(-(((v*4/50)/2)^1.8) * (b/max_pow)^8 * 2000))
  
  p2 <- ggplot(df2, aes(x, Power, colour = cells, linetype = reads)) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
    theme_bw(base_size = 16) +
    theme(aspect.ratio = 1) +
    labs(x = "Fold-change (percent)", y = "Power", colour = "Cells", linetype = "Reads/cell")

  # Combine with shared legend below
  (p1 + p2) +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
})