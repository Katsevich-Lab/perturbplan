# Main server logic with real power calculations
# Updated to use perturbplan package for statistical analysis

server <- function(input, output, session) {
  # Helpers
  is_sel <- function(x) identical(sel$type, x)
  toggle  <- function(v, x) if (x %in% v) setdiff(v, x) else c(v, x)
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
      tpm_threshold = input$tpm_threshold,
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

  # Selection state
  sel <- reactiveValues(type=NULL, idx=integer(0),
                        tiles=data.frame(row=integer(0), col=integer(0)))
  slice_mode <- reactiveVal(NULL)
  slice_x    <- reactiveVal(numeric(0))

  # Clear selections
  observeEvent(input$clear, {
    sel$type <- NULL; sel$idx <- integer(0); sel$tiles <- sel$tiles[0,]
    slice_mode(NULL); slice_x(numeric(0))
    updateTextInput(session,"overall_points","")
    updateTextInput(session,"slice_points","")
  })

  # Overall UI box
  output$overall_box_ui <- renderUI({
    req(planned())
    if (input$mode=="cells") {
      vals <- cells_seq()[sel$idx]
      textInput("overall_points","Number of cells (comma-separated):",paste(vals,collapse=", "))
    } else if (input$mode=="reads") {
      vals <- reads_seq()[sel$idx]
      textInput("overall_points","Number of reads per cell (comma-separated):",paste(vals,collapse=", "))
    } else {
      if (nrow(sel$tiles)) {
        vals <- sprintf("%d×%d",
                        cells_seq()[sel$tiles$row],
                        reads_seq()[sel$tiles$col])
        txt  <- paste(vals, collapse=", ")
      } else txt <- ""
      textInput("overall_points","Tiles (cells × reads):", txt)
    }
  })

  # Parse overall_points edits
  observeEvent(input$overall_points, ignoreInit=TRUE, {
    req(planned())
    txt <- gsub("\\s","", input$overall_points)
    if (input$mode=="cells") {
      nums <- as.numeric(strsplit(txt,",")[[1]])
      sel$type <- "row"
      sel$idx  <- which(cells_seq() %in% nums)
      slice_mode("row")
    } else if (input$mode=="reads") {
      nums <- as.numeric(strsplit(txt,",")[[1]])
      sel$type <- "col"
      sel$idx  <- which(reads_seq() %in% nums)
      slice_mode("col")
    } else {
      if (nzchar(txt)) {
        pairs <- strsplit(txt,",")[[1]]
        new   <- do.call(rbind, lapply(pairs, function(p) {
          nums <- as.numeric(strsplit(p,"[×x]")[[1]])
          if (length(nums)==2)
            data.frame(row=match(nums[1],cells_seq()),
                       col=match(nums[2],reads_seq()))
        }))
        sel$type  <- "tile"
        sel$tiles <- na.omit(unique(new))
      } else {
        sel$type  <- NULL
        sel$tiles <- sel$tiles[0,]
      }
    }
  })

  # Heatmap clicks
  observeEvent(input$heat_click, {
    req(planned())
    r <- which.min(abs(cells_seq() - input$heat_click$y))
    c <- which.min(abs(reads_seq() - input$heat_click$x))
    if (input$mode=="cells") {
      sel$type <- "row"; sel$idx <- toggle(sel$idx,r); slice_mode("row")
      updateTextInput(session,"overall_points",
                      value=paste(cells_seq()[sel$idx],collapse=", "))
    } else if (input$mode=="reads") {
      sel$type <- "col"; sel$idx <- toggle(sel$idx,c); slice_mode("col")
      updateTextInput(session,"overall_points",
                      value=paste(reads_seq()[sel$idx],collapse=", "))
    } else {
      sel$type <- "tile"
      hit <- with(sel$tiles, which(row==r & col==c))
      if (length(hit)) sel$tiles <- sel$tiles[-hit,]
      else sel$tiles <- rbind(sel$tiles, data.frame(row=r,col=c))
      txt <- if (nrow(sel$tiles)) {
        paste(sprintf("%d×%d",
                      cells_seq()[sel$tiles$row],
                      reads_seq()[sel$tiles$col]),
              collapse=", ")
      } else ""
      updateTextInput(session,"overall_points",value=txt)
    }
  })

  # Enable Go button
  observe({
    ok <- planned() && (
      (is_sel("row")  && length(sel$idx)>=1) ||
        (is_sel("col")  && length(sel$idx)>=1) ||
        (is_sel("tile") && nrow(sel$tiles)>=1)
    )
    toggleState("go_overall", ok)
  })
  observeEvent(input$go_overall, {
    dest <- if (is_sel("tile")) "per_pair" else "overall_slice"
    updateTabItems(session, "sidebar", selected = dest)
  })

  # Render heatmap
  draw_rects <- function(df,p)
    p + geom_rect(data=df,
                  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                  inherit.aes=FALSE,colour="red",linewidth=1,fill=NA)
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

  # Slice tab UI
  output$slice_box_ui <- renderUI({
    req(planned(), !is.null(slice_mode()))
    lab <- if (identical(slice_mode(),"row"))
      "Drill down by number of reads / cell:<br/>(click to select multiple points)"
    else "Drill down by number of cells:<br/>(click to select multiple points)"
    textInput("slice_points", HTML(lab),
              if (length(slice_x())) paste(slice_x(), collapse=", ") else "")
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
      ggplot(sub,aes(reads,power,colour=factor(cells)))+
        geom_line()+
        geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
        geom_vline(xintercept=slice_x(),colour="red")+
        theme_bw(base_size = 16)+
        theme(aspect.ratio = 1) +
        labs(x="Reads per cell",y="Power",colour="Cells")
    } else {
      sub <- subset(df, reads %in% reads_seq()[sel$idx])
      ggplot(sub,aes(cells,power,colour=factor(reads)))+
        geom_line()+
        geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
        geom_vline(xintercept=slice_x(),colour="red")+
        theme_bw(base_size = 16)+
        theme(aspect.ratio = 1) +
        labs(x="Number of cells",y="Power",colour="Reads")
    }
  })

  # Slice interactions - support multiple selections
  observeEvent(input$slice_click,{
    req(slice_mode()%in%c("row","col"))
    x <- if (slice_mode()=="row")
      reads_seq()[which.min(abs(reads_seq()-input$slice_click$x))]
    else
      cells_seq()[which.min(abs(cells_seq()-input$slice_click$x))]
    
    # Toggle selection
    current_vals <- slice_x()
    if (x %in% current_vals) {
      new_vals <- setdiff(current_vals, x)
    } else {
      new_vals <- c(current_vals, x)
    }
    
    slice_x(new_vals)
    updateTextInput(session,"slice_points",value=paste(new_vals, collapse=", "))
  })
  observeEvent(input$slice_points,ignoreInit=TRUE,{
    txt <- gsub("\\s","", input$slice_points)
    if (nzchar(txt)) {
      v <- as.numeric(strsplit(txt,",")[[1]])
      v <- v[!is.na(v)]
      slice_x(v)
    } else {
      slice_x(numeric(0))
    }
  })
  observe({
    toggleState("go_slice", length(slice_x())>=1)
  })
  observeEvent(input$slice_clear,{
    slice_x(numeric(0)); updateTextInput(session,"slice_points",value="")
  })

  # Slice -> Per-pair (support multiple selections)
  observeEvent(input$go_slice,{
    req(length(slice_x())>=1)
    if (slice_mode()=="row")
      sel$tiles <- expand.grid(row=sel$idx,
                               col=match(slice_x(),reads_seq()))
    else
      sel$tiles <- expand.grid(row=match(slice_x(),cells_seq()),
                               col=sel$idx)
    sel$type <- "tile"
    updateTabItems(session, "sidebar", selected = "per_pair")
  })

  # Per-pair plots with real power curves
  output$pp_combined <- renderPlot({
    req(planned(), is_sel("tile"))

    # Get real power curves from calculate_power_grid results
    full_power_results <- power_results()
    power_grid <- full_power_results$power_grid
    
    # Extract power curve data for selected tiles  
    create_real_plot_data <- function(curve_type) {
      dfs <- Map(function(r, c) {
        # Get the actual cell and read values for this tile
        cell_val <- cells_seq()[r]
        read_val <- reads_seq()[c]
        
        # Find matching row in the grid
        matching_row <- which(power_grid$cells == cell_val & power_grid$reads == read_val)
        
        if (length(matching_row) > 0) {
          # Get curve data from the power_curves section
          if (curve_type == "fc") {
            curve_data <- full_power_results$power_curves$fc_curves[[matching_row[1]]]
          } else {
            curve_data <- full_power_results$power_curves$expr_curves[[matching_row[1]]]
          }
          
          if (!is.null(curve_data) && nrow(curve_data) > 0) {
            if (curve_type == "fc") {
              data.frame(
                x = curve_data$fold_change,
                Power = curve_data$power,
                label = sprintf("%d × %d", cell_val, read_val)
              )
            } else {
              # Convert relative_expression to TPM scale
              tpm_values <- curve_data$relative_expression * 1e6
              data.frame(
                x = tpm_values,
                Power = curve_data$power,
                label = sprintf("%d × %d", cell_val, read_val)
              )
            }
          } else {
            data.frame(x = numeric(0), Power = numeric(0), label = character(0))
          }
        } else {
          data.frame(x = numeric(0), Power = numeric(0), label = character(0))
        }
      }, sel$tiles$row, sel$tiles$col)
      
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
}