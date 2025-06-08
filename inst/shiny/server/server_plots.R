# Server module for generating plots and visualizations

create_plots_server <- function(input, output, session, power_data, selection_data) {
  
  # Render heatmap
  draw_rects <- function(df,p)
    p + geom_rect(data=df,
                  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                  inherit.aes=FALSE,colour="red",linewidth=1,fill=NA)
  
  output$heat <- renderPlot({
    req(power_data$planned())
    p <- ggplot(power_data$gridDF(),aes(reads,cells,fill=power))+
      geom_tile(colour=NA)+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_bw(base_size = 16)+
      theme(panel.grid=element_blank(),
            aspect.ratio = 1)+
      labs(x="Reads per cell",y="Number of cells",fill="Power")
    if (power_data$is_sel("row") && length(selection_data$sel$idx)) {
      rect <- data.frame(
        xmin=min(power_data$reads_seq())-power_data$dreads()/2,
        xmax=max(power_data$reads_seq())+power_data$dreads()/2,
        ymin=power_data$cells_seq()[selection_data$sel$idx]-power_data$dcells()/2,
        ymax=power_data$cells_seq()[selection_data$sel$idx]+power_data$dcells()/2
      )
      p <- draw_rects(rect,p)
    } else if (power_data$is_sel("col") && length(selection_data$sel$idx)) {
      rect <- data.frame(
        xmin=power_data$reads_seq()[selection_data$sel$idx]-power_data$dreads()/2,
        xmax=power_data$reads_seq()[selection_data$sel$idx]+power_data$dreads()/2,
        ymin=min(power_data$cells_seq())-power_data$dcells()/2,
        ymax=max(power_data$cells_seq())+power_data$dcells()/2
      )
      p <- draw_rects(rect,p)
    } else if (power_data$is_sel("tile") && nrow(selection_data$sel$tiles)) {
      rect <- data.frame(
        xmin=power_data$reads_seq()[selection_data$sel$tiles$col]-power_data$dreads()/2,
        xmax=power_data$reads_seq()[selection_data$sel$tiles$col]+power_data$dreads()/2,
        ymin=power_data$cells_seq()[selection_data$sel$tiles$row]-power_data$dcells()/2,
        ymax=power_data$cells_seq()[selection_data$sel$tiles$row]+power_data$dcells()/2
      )
      p <- draw_rects(rect,p)
    }
    p
  })

  # Slice tab UI
  output$slice_box_ui <- renderUI({
    req(power_data$planned(), !is.null(selection_data$slice_mode()))
    lab <- if (identical(selection_data$slice_mode(),"row"))
      "Drill down by number of reads / cell:<br/>(click to select multiple points)"
    else "Drill down by number of cells:<br/>(click to select multiple points)"
    textInput("slice_points", HTML(lab),
              if (length(selection_data$slice_x())) paste(selection_data$slice_x(), collapse=", ") else "")
  })

  # Slice title & plot
  output$slice_title <- renderText({
    if (identical(selection_data$slice_mode(),"row"))
      sprintf("Power versus reads per cell")
    else if (identical(selection_data$slice_mode(),"col"))
      sprintf("Power vs number of cells")
    else "Slice view"
  })
  
  output$slice_plot <- renderPlot({
    req(power_data$planned(), !is.null(selection_data$slice_mode()), length(selection_data$sel$idx))
    df <- power_data$gridDF()
    if (selection_data$slice_mode()=="row") {
      sub <- subset(df, cells %in% power_data$cells_seq()[selection_data$sel$idx])
      ggplot(sub,aes(reads,power,colour=factor(cells)))+
        geom_line()+
        geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
        geom_vline(xintercept=selection_data$slice_x(),colour="red")+
        theme_bw(base_size = 16)+
        theme(aspect.ratio = 1) +
        labs(x="Reads per cell",y="Power",colour="Cells")
    } else {
      sub <- subset(df, reads %in% power_data$reads_seq()[selection_data$sel$idx])
      ggplot(sub,aes(cells,power,colour=factor(reads)))+
        geom_line()+
        geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
        geom_vline(xintercept=selection_data$slice_x(),colour="red")+
        theme_bw(base_size = 16)+
        theme(aspect.ratio = 1) +
        labs(x="Number of cells",y="Power",colour="Reads")
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