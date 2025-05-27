# Main server logic
# Extracted from app-prototype.R and modularized

server <- function(input, output, session) {
  
  # Source server modules
  source("server/reactive_data.R", local = TRUE)
  source("server/selection_logic.R", local = TRUE)  
  source("server/plot_rendering.R", local = TRUE)
  
  # Core reactive values and helpers
  is_sel <- function(x) identical(sel$type, x)
  toggle  <- function(v, x) if (x %in% v) setdiff(v, x) else c(v, x)
  planned <- reactiveVal(FALSE)
  
  # UI state management
  output$need_plan <- reactive(!planned())
  outputOptions(output,  "need_plan", suspendWhenHidden = FALSE)
  observeEvent(input$plan_btn, planned(TRUE))
  
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
  
  # Make reactive values available to other modules
  values <- reactiveValues(
    sel = sel,
    slice_mode = slice_mode,
    slice_x = slice_x,
    planned = planned,
    is_sel = is_sel,
    toggle = toggle
  )
}