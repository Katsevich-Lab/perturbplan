# Selection and interaction logic
# Handles user clicks and selections on plots

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
updating_from_click <- reactiveVal(FALSE)
observeEvent(input$overall_points, ignoreInit=TRUE, {
  if (updating_from_click()) return()
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
  updating_from_click(TRUE)
  
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
  
  # Reset flag after a short delay  
  Sys.sleep(0.01)
  updating_from_click(FALSE)
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
  dest <- if (is_sel("tile")) "Per-pair power" else "Overall power (slice)"
  updateTabsetPanel(session,"tabs",selected=dest)
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

# Slice interactions - now supports multiple selections
updating_slice_from_click <- reactiveVal(FALSE)
observeEvent(input$slice_click,{
  req(slice_mode()%in%c("row","col"))
  updating_slice_from_click(TRUE)
  
  x <- if (slice_mode()=="row")
    reads_seq()[which.min(abs(reads_seq()-input$slice_click$x))]
  else
    cells_seq()[which.min(abs(cells_seq()-input$slice_click$x))]
  
  # Toggle selection - add if not present, remove if already selected
  current_vals <- slice_x()
  if (x %in% current_vals) {
    new_vals <- setdiff(current_vals, x)
  } else {
    new_vals <- c(current_vals, x)
  }
  
  slice_x(new_vals)
  updateTextInput(session,"slice_points",value=paste(new_vals, collapse=", "))
  
  # Reset flag after a short delay
  Sys.sleep(0.01)
  updating_slice_from_click(FALSE)
})

observeEvent(input$slice_points,ignoreInit=TRUE,{
  if (updating_slice_from_click()) return()
  txt <- gsub("\\s","", input$slice_points)
  if (nzchar(txt)) {
    v <- as.numeric(strsplit(txt,",")[[1]])
    v <- v[!is.na(v)]  # Remove any NAs
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

# Slice -> Per-pair
observeEvent(input$go_slice,{
  req(length(slice_x())>=1)
  if (slice_mode()=="row") {
    # slice_x contains reads_per_cell values, sel$idx contains row indices (cell counts)
    sel$tiles <- expand.grid(row=sel$idx,
                             col=match(slice_x(),reads_seq()))
  } else {
    # slice_x contains cell values, sel$idx contains col indices (reads_per_cell)  
    sel$tiles <- expand.grid(row=match(slice_x(),cells_seq()),
                             col=sel$idx)
  }
  sel$type <- "tile"
  updateTabsetPanel(session,"tabs",selected="Per-pair power")
})