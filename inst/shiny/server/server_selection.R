# Server module for handling selections and interactions

create_selection_server <- function(input, output, session, power_data) {
  # Selection state
  sel <- reactiveValues(type=NULL, idx=integer(0),
                        tiles=data.frame(row=integer(0), col=integer(0)))
  slice_mode <- reactiveVal(NULL)
  slice_x    <- reactiveVal(numeric(0))
  
  # Helper functions that use sel
  is_sel <- function(x) identical(sel$type, x)
  toggle  <- function(v, x) if (x %in% v) setdiff(v, x) else c(v, x)

  # Clear selections
  observeEvent(input$clear, {
    sel$type <- NULL; sel$idx <- integer(0); sel$tiles <- sel$tiles[0,]
    slice_mode(NULL); slice_x(numeric(0))
    updateTextInput(session,"overall_points","")
    updateTextInput(session,"slice_points","")
  })

  # Overall UI box
  output$overall_box_ui <- renderUI({
    req(power_data$planned())
    if (input$mode=="cells") {
      vals <- power_data$cells_seq()[sel$idx]
      textInput("overall_points","Number of cells (comma-separated):",paste(vals,collapse=", "))
    } else if (input$mode=="reads") {
      vals <- power_data$reads_seq()[sel$idx]
      textInput("overall_points","Number of reads per cell (comma-separated):",paste(vals,collapse=", "))
    } else {
      if (nrow(sel$tiles)) {
        vals <- sprintf("%d×%d",
                        power_data$cells_seq()[sel$tiles$row],
                        power_data$reads_seq()[sel$tiles$col])
        txt  <- paste(vals, collapse=", ")
      } else txt <- ""
      textInput("overall_points","Tiles (cells × reads):", txt)
    }
  })

  # Parse overall_points edits
  observeEvent(input$overall_points, ignoreInit=TRUE, {
    req(power_data$planned())
    txt <- gsub("\\s","", input$overall_points)
    if (input$mode=="cells") {
      nums <- as.numeric(strsplit(txt,",")[[1]])
      sel$type <- "row"
      sel$idx  <- which(power_data$cells_seq() %in% nums)
      slice_mode("row")
    } else if (input$mode=="reads") {
      nums <- as.numeric(strsplit(txt,",")[[1]])
      sel$type <- "col"
      sel$idx  <- which(power_data$reads_seq() %in% nums)
      slice_mode("col")
    } else {
      if (nzchar(txt)) {
        pairs <- strsplit(txt,",")[[1]]
        new   <- do.call(rbind, lapply(pairs, function(p) {
          nums <- as.numeric(strsplit(p,"[×x]")[[1]])
          if (length(nums)==2)
            data.frame(row=match(nums[1],power_data$cells_seq()),
                       col=match(nums[2],power_data$reads_seq()))
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
    req(power_data$planned())
    r <- which.min(abs(power_data$cells_seq() - input$heat_click$y))
    c <- which.min(abs(power_data$reads_seq() - input$heat_click$x))
    if (input$mode=="cells") {
      sel$type <- "row"; sel$idx <- toggle(sel$idx,r); slice_mode("row")
      updateTextInput(session,"overall_points",
                      value=paste(power_data$cells_seq()[sel$idx],collapse=", "))
    } else if (input$mode=="reads") {
      sel$type <- "col"; sel$idx <- toggle(sel$idx,c); slice_mode("col")
      updateTextInput(session,"overall_points",
                      value=paste(power_data$reads_seq()[sel$idx],collapse=", "))
    } else {
      sel$type <- "tile"
      hit <- with(sel$tiles, which(row==r & col==c))
      if (length(hit)) sel$tiles <- sel$tiles[-hit,]
      else sel$tiles <- rbind(sel$tiles, data.frame(row=r,col=c))
      txt <- if (nrow(sel$tiles)) {
        paste(sprintf("%d×%d",
                      power_data$cells_seq()[sel$tiles$row],
                      power_data$reads_seq()[sel$tiles$col]),
              collapse=", ")
      } else ""
      updateTextInput(session,"overall_points",value=txt)
    }
  })

  # Enable Go button
  observe({
    ok <- power_data$planned() && (
      (is_sel("row")  && length(sel$idx)>=1) ||
        (is_sel("col")  && length(sel$idx)>=1) ||
        (is_sel("tile") && nrow(sel$tiles)>=1)
    )
    toggleState("go_overall", ok)
  })
  
  observeEvent(input$go_overall, {
    dest <- if (is_sel("tile")) "per_pair" else "overall_slice"
    updateTabsetPanel(session, "main_tabs", selected = dest)
  })

  # Return selection state and helper functions for use by other modules
  return(list(
    sel = sel,
    slice_mode = slice_mode,
    slice_x = slice_x,
    is_sel = is_sel,
    toggle = toggle
  ))
}