# Server module for generating plots and visualizations

create_plots_server <- function(input, output, session, power_data, selection_data) {
  
  # Render heatmap
  draw_rects <- function(df,p)
    p + geom_rect(data=df,
                  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                  inherit.aes=FALSE,colour="red",linewidth=1,fill=NA)
  
  # Base heatmap (no selection dependency to prevent flashing)
  base_heatmap <- reactive({
    req(power_data$planned())
    ggplot(power_data$gridDF(),aes(reads,cells,fill=power))+
      geom_tile(colour=NA)+
      scale_x_log10(expand = c(0.02, 0.02), labels = scales::comma_format()) +
      scale_y_log10(expand = c(0.02, 0.02), labels = scales::comma_format()) +
      theme_bw(base_size = 16)+
      theme(panel.grid=element_blank(),
            aspect.ratio = 1)+
      labs(x="Reads per cell (log scale)",y="Number of cells per target (log scale)",fill="Power")
  })
  
  # Helper function to calculate logarithmic tile boundaries
  get_log_tile_bounds <- function(values, idx) {
    if (length(values) == 1) {
      # Single value - create small bounds around it
      factor <- 1.2  # 20% extension
      return(list(min = values / factor, max = values * factor))
    }
    
    # Multiple values - use geometric spacing
    if (idx == 1) {
      # First tile
      ratio <- values[2] / values[1]
      min_val <- values[1] / sqrt(ratio)
      max_val <- (values[1] + values[2]) / 2
    } else if (idx == length(values)) {
      # Last tile
      ratio <- values[idx] / values[idx-1]
      min_val <- (values[idx-1] + values[idx]) / 2
      max_val <- values[idx] * sqrt(ratio)
    } else {
      # Middle tiles
      min_val <- (values[idx-1] + values[idx]) / 2
      max_val <- (values[idx] + values[idx+1]) / 2
    }
    
    return(list(min = min_val, max = max_val))
  }

  # Selection overlay (isolated from base heatmap)
  selection_overlay <- reactive({
    req(power_data$planned())
    
    reads_seq <- power_data$reads_seq()
    cells_seq <- power_data$cells_seq()
    
    if (selection_data$is_sel("row") && length(selection_data$sel$idx)) {
      # Row selection - highlight entire rows
      y_bounds <- lapply(selection_data$sel$idx, function(i) {
        get_log_tile_bounds(cells_seq, i)
      })
      
      data.frame(
        xmin = min(reads_seq) * 0.8,  # Extend slightly beyond data range
        xmax = max(reads_seq) * 1.2,
        ymin = sapply(y_bounds, function(b) b$min),
        ymax = sapply(y_bounds, function(b) b$max)
      )
    } else if (selection_data$is_sel("col") && length(selection_data$sel$idx)) {
      # Column selection - highlight entire columns
      x_bounds <- lapply(selection_data$sel$idx, function(i) {
        get_log_tile_bounds(reads_seq, i)
      })
      
      data.frame(
        xmin = sapply(x_bounds, function(b) b$min),
        xmax = sapply(x_bounds, function(b) b$max),
        ymin = min(cells_seq) * 0.8,  # Extend slightly beyond data range
        ymax = max(cells_seq) * 1.2
      )
    } else if (selection_data$is_sel("tile") && nrow(selection_data$sel$tiles)) {
      # Tile selection - highlight individual tiles
      x_bounds <- lapply(selection_data$sel$tiles$col, function(i) {
        get_log_tile_bounds(reads_seq, i)
      })
      y_bounds <- lapply(selection_data$sel$tiles$row, function(i) {
        get_log_tile_bounds(cells_seq, i)
      })
      
      data.frame(
        xmin = sapply(x_bounds, function(b) b$min),
        xmax = sapply(x_bounds, function(b) b$max),
        ymin = sapply(y_bounds, function(b) b$min),
        ymax = sapply(y_bounds, function(b) b$max)
      )
    } else {
      NULL
    }
  })
  
  # Main heatmap output (combines base + overlay efficiently)
  output$heat <- renderPlot({
    p <- base_heatmap()
    overlay <- selection_overlay()
    if (!is.null(overlay)) {
      p <- draw_rects(overlay, p)
    }
    p
  })

  # Check if slice is available (has selections)
  output$slice_available <- reactive({
    power_data$planned() && !is.null(selection_data$slice_mode()) && length(selection_data$sel$idx) > 0
  })
  outputOptions(output, "slice_available", suspendWhenHidden = FALSE)

  # Slice tab UI
  output$slice_box_ui <- renderUI({
    req(power_data$planned(), !is.null(selection_data$slice_mode()))
    lab <- if (identical(selection_data$slice_mode(),"row"))
      "Drill down by number of reads per cell:<br/>(click to select multiple points)"
    else "Drill down by number of cells per target:<br/>(click to select multiple points)"
    textInput("slice_points", HTML(lab),
              if (length(selection_data$slice_x())) paste(selection_data$slice_x(), collapse=", ") else "")
  })

  # Slice title & plot
  output$slice_title <- renderText({
    if (identical(selection_data$slice_mode(),"row"))
      sprintf("Power versus reads per cell")
    else if (identical(selection_data$slice_mode(),"col"))
      sprintf("Power vs number of cells per target")
    else "Slice view"
  })
  
  output$slice_plot <- renderPlot({
    req(power_data$planned(), !is.null(selection_data$slice_mode()), length(selection_data$sel$idx))
    df <- power_data$gridDF()
    if (selection_data$slice_mode()=="row") {
      sub <- subset(df, cells %in% power_data$cells_seq()[selection_data$sel$idx])
      ggplot(sub,aes(reads,power,colour=factor(cells)))+
        geom_line()+
        geom_point()+
        geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
        geom_vline(xintercept=selection_data$slice_x(),colour="red")+
        scale_x_log10(labels = scales::comma_format()) +
        theme_bw(base_size = 16)+
        theme(aspect.ratio = 1) +
        labs(x="Reads per cell (log scale)",y="Power",colour="Cells per target")
    } else {
      sub <- subset(df, reads %in% power_data$reads_seq()[selection_data$sel$idx])
      ggplot(sub,aes(cells,power,colour=factor(reads)))+
        geom_line()+
        geom_point()+
        geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
        geom_vline(xintercept=selection_data$slice_x(),colour="red")+
        scale_x_log10(labels = scales::comma_format()) +
        theme_bw(base_size = 16)+
        theme(aspect.ratio = 1) +
        labs(x="Number of cells per target (log scale)",y="Power",colour="Reads per cell")
    }
  })

  # Slice interactions - support multiple selections
  observeEvent(input$slice_click,{
    req(selection_data$slice_mode()%in%c("row","col"))
    x <- if (selection_data$slice_mode()=="row")
      power_data$reads_seq()[which.min(abs(power_data$reads_seq()-input$slice_click$x))]
    else
      power_data$cells_seq()[which.min(abs(power_data$cells_seq()-input$slice_click$x))]
    
    # Toggle selection
    current_vals <- selection_data$slice_x()
    if (x %in% current_vals) {
      new_vals <- setdiff(current_vals, x)
    } else {
      new_vals <- c(current_vals, x)
    }
    
    selection_data$slice_x(new_vals)
    updateTextInput(session,"slice_points",value=paste(new_vals, collapse=", "))
  })
  
  observeEvent(input$slice_points,ignoreInit=TRUE,{
    txt <- gsub("\\s","", input$slice_points)
    if (nzchar(txt)) {
      v <- as.numeric(strsplit(txt,",")[[1]])
      v <- v[!is.na(v)]
      selection_data$slice_x(v)
    } else {
      selection_data$slice_x(numeric(0))
    }
  })
  
  observe({
    toggleState("go_slice", length(selection_data$slice_x())>=1)
  })
  
  observeEvent(input$slice_clear,{
    selection_data$slice_x(numeric(0)); updateTextInput(session,"slice_points",value="")
  })

  # Slice -> Per-pair (support multiple selections)
  observeEvent(input$go_slice,{
    req(length(selection_data$slice_x())>=1)
    if (selection_data$slice_mode()=="row")
      selection_data$sel$tiles <- expand.grid(row=selection_data$sel$idx,
                               col=match(selection_data$slice_x(),power_data$reads_seq()))
    else
      selection_data$sel$tiles <- expand.grid(row=match(selection_data$slice_x(),power_data$cells_seq()),
                               col=selection_data$sel$idx)
    selection_data$sel$type <- "tile"
    updateTabsetPanel(session, "main_tabs", selected = "per_pair")
  })

}