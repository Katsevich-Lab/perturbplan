# DEPRECATED: Use app.R instead
# This file contains old placeholder power calculations
# The current app with real statistical power analysis is in app.R
# 
# app.R  –  Perturb-seq power planner (single interactive Overall box)

library(shiny)
library(shinyBS)
library(shinyjs)
library(ggplot2)
library(patchwork)

# UI
ui <- fluidPage(
  useShinyjs(),

  sidebarLayout(
    # ------------------------------------------------------------------
    sidebarPanel(
      width = 3,
      tags$div(
        style = "overflow-y:auto; height:92vh; position:sticky; top:0;",
        bsCollapse(
          id   = "panels",
          open = "Experimental choices",

          # --- Experimental choices -----------------------------------
          bsCollapsePanel(
            "Experimental choices",
            selectInput("biological_system",   "Biological system:",   c("K562")),
            selectInput("experimental_platform","Experimental platform:", c("10x Chromium v3")),
            numericInput("MOI",               "MOI:",                 10, 1, 30, 0.5),
            numericInput("num_targets",       "Number of targets:",   100, 1, 1000),
            numericInput("gRNAs_per_target",  "Number of gRNAs per target:",     4, 1, 10),
            numericInput("non_targeting_gRNAs","Number of non-targeting gRNAs:",10, 0, 100)
          ),

          # --- Analysis choices ---------------------------------------
          bsCollapsePanel(
            "Analysis choices",
            numericInput("num_pairs",     "Number of pairs analyzed:",  1000, 10, 10000),
            numericInput("tpm_threshold", "Minimum TPM threshold:",      10, 0, 10, 0.5),
            numericInput("fdr_target",    "FDR target level:", 0.05, 0.001, 0.1, 0.001)
          ),

          # --- Assumed effect sizes -----------------------------------
          bsCollapsePanel(
            "Assumed effect sizes",
            numericInput("fc_mean", "Fold-change mean:", 0.85, 1.1, 10, 0.05),
            numericInput("fc_sd",   "Fold-change SD:",   0.15, 0.1,  5, 0.05),
            numericInput("prop_non_null", "Proportion of non-null pairs:", 0.1, 0, 1, 0.01)
          )
        ),

        hr(),
        actionButton("plan_btn", "Plan", class = "btn-success", width = "100%")
      )
    ),

    # ------------------------------------------------------------------
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        # === Overall tab ==============================================
        tabPanel(
          "Overall power (heatmap)",

          # Message before planning
          conditionalPanel(
            condition = "output.need_plan",
            h4("Select your parameters and click \"Plan\"")
          ),

          # Heat-map and controls (after Plan)
          conditionalPanel(
            condition = "!output.need_plan",
            h4("Power versus number of cells and reads per cell"),
            fluidRow(
              column(
                8,
                plotOutput("heat", height = 460, click = "heat_click")
              ),
              column(
                4,
                wellPanel(
                  radioButtons(
                    "mode", "Drill down by:",
                    c("Number of cells (click one or more rows)"   = "cells",
                      "Reads per cell (click one or more columns)" = "reads",
                      "Both (click one or more tiles)"             = "tile"),
                    selected = "cells"
                  ),
                  uiOutput("overall_box_ui"),
                  actionButton("go_overall", "Go",   class = "btn-primary", width = "100%"),
                  actionButton("clear",      "Clear",                       width = "100%")
                )
              )
            )
          )
        ),

        # === Slice tab ===============================================
        tabPanel(
          "Overall power (slice)",
          h4(textOutput("slice_title")),
          fluidRow(
            column(
              8,
              plotOutput("slice_plot", height = 460, click = "slice_click")
            ),
            column(
              4,
              wellPanel(
                uiOutput("slice_box_ui"),
                actionButton("go_slice",    "Go",   class = "btn-primary", width = "100%"),
                actionButton("slice_clear", "Clear",                       width = "100%")
              )
            )
          )
        ),

        # === Per-pair tab ============================================
        tabPanel(
          "Per-pair power",
          h4(textOutput("pair_title")),
          fluidRow(
            column(
              8,
              plotOutput("pp_combined", height = "400px", width = "100%")
            ),
            column(
              4
              # empty spacer
            )
          )
        )
      )
    )
  )
)



# Server
server <- function(input, output, session) {
  # Helpers
  is_sel <- function(x) identical(sel$type, x)
  toggle  <- function(v, x) if (x %in% v) setdiff(v, x) else c(v, x)
  planned <- reactiveVal(FALSE)
  output$need_plan <- reactive(!planned())
  outputOptions(output,  "need_plan", suspendWhenHidden = FALSE)
  observeEvent(input$plan_btn, planned(TRUE))

  # Grid setup
  cells_seq <- round(seq(200, 10000, length.out=20))
  reads_seq <- round(seq(2000, 50000, length.out=20))
  dcells <- diff(cells_seq)[1]
  dreads <- diff(reads_seq)[1]
  max_pow <- max(log10(cells_seq)*log10(reads_seq))
  gridDF <- reactive({
    df <- expand.grid(cells=cells_seq, reads=reads_seq)
    df$power <- pmin(1,(log10(df$cells)*log10(df$reads))/max_pow)
    df
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
    req(planned())      # <- hides the whole UI box until planned() is TRUE
    if (input$mode=="cells") {
      vals <- cells_seq[sel$idx]
      textInput("overall_points","Number of cells (comma-separated):",paste(vals,collapse=", "))
    } else if (input$mode=="reads") {
      vals <- reads_seq[sel$idx]
      textInput("overall_points","Number of reads per cell (comma-separated):",paste(vals,collapse=", "))
    } else {
      if (nrow(sel$tiles)) {
        vals <- sprintf("%d×%d",
                        cells_seq[sel$tiles$row],
                        reads_seq[sel$tiles$col])
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
      sel$idx  <- which(cells_seq %in% nums)
      slice_mode("row")
    } else if (input$mode=="reads") {
      nums <- as.numeric(strsplit(txt,",")[[1]])
      sel$type <- "col"
      sel$idx  <- which(reads_seq %in% nums)
      slice_mode("col")
    } else {
      if (nzchar(txt)) {
        pairs <- strsplit(txt,",")[[1]]
        new   <- do.call(rbind, lapply(pairs, function(p) {
          nums <- as.numeric(strsplit(p,"[×x]")[[1]])
          if (length(nums)==2)
            data.frame(row=match(nums[1],cells_seq),
                       col=match(nums[2],reads_seq))
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
    r <- which.min(abs(cells_seq - input$heat_click$y))
    c <- which.min(abs(reads_seq - input$heat_click$x))
    if (input$mode=="cells") {
      sel$type <- "row"; sel$idx <- toggle(sel$idx,r); slice_mode("row")
      updateTextInput(session,"overall_points",
                      value=paste(cells_seq[sel$idx],collapse=", "))
    } else if (input$mode=="reads") {
      sel$type <- "col"; sel$idx <- toggle(sel$idx,c); slice_mode("col")
      updateTextInput(session,"overall_points",
                      value=paste(reads_seq[sel$idx],collapse=", "))
    } else {
      sel$type <- "tile"
      hit <- with(sel$tiles, which(row==r & col==c))
      if (length(hit)) sel$tiles <- sel$tiles[-hit,]
      else sel$tiles <- rbind(sel$tiles, data.frame(row=r,col=c))
      txt <- if (nrow(sel$tiles)) {
        paste(sprintf("%d×%d",
                      cells_seq[sel$tiles$row],
                      reads_seq[sel$tiles$col]),
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
    dest <- if (is_sel("tile")) "Per-pair power" else "Overall power (slice)"
    updateTabsetPanel(session,"tabs",selected=dest)
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
        xmin=min(reads_seq)-dreads/2,
        xmax=max(reads_seq)+dreads/2,
        ymin=cells_seq[sel$idx]-dcells/2,
        ymax=cells_seq[sel$idx]+dcells/2
      )
      p <- draw_rects(rect,p)
    } else if (is_sel("col") && length(sel$idx)) {
      rect <- data.frame(
        xmin=reads_seq[sel$idx]-dreads/2,
        xmax=reads_seq[sel$idx]+dreads/2,
        ymin=min(cells_seq)-dcells/2,
        ymax=max(cells_seq)+dcells/2
      )
      p <- draw_rects(rect,p)
    } else if (is_sel("tile") && nrow(sel$tiles)) {
      rect <- data.frame(
        xmin=reads_seq[sel$tiles$col]-dreads/2,
        xmax=reads_seq[sel$tiles$col]+dreads/2,
        ymin=cells_seq[sel$tiles$row]-dcells/2,
        ymax=cells_seq[sel$tiles$row]+dcells/2
      )
      p <- draw_rects(rect,p)
    }
    p
  })

  # Slice tab UI
  output$slice_box_ui <- renderUI({
    req(planned(), !is.null(slice_mode()))   # only after overall plot + slice mode chosen
    lab <- if (identical(slice_mode(),"row"))
      "Drill down by number of reads / cell:<br/>(click the plot)"
    else "Drill down by number of cells:<br/>(click the plot)"
    textInput("slice_points", HTML(lab),
              if (length(slice_x())) slice_x() else "")
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
      sub <- subset(df, cells %in% cells_seq[sel$idx])
      ggplot(sub,aes(reads,power,colour=factor(cells)))+
        geom_line()+
        geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
        geom_vline(xintercept=slice_x(),colour="red")+
        theme_bw(base_size = 16)+
        theme(aspect.ratio = 1) +
        labs(x="Reads per cell",y="Power",colour="Cells")
    } else {
      sub <- subset(df, reads %in% reads_seq[sel$idx])
      ggplot(sub,aes(cells,power,colour=factor(reads)))+
        geom_line()+
        geom_hline(yintercept=0.8,linetype="dashed",colour="grey") +
        geom_vline(xintercept=slice_x(),colour="red")+
        theme_bw(base_size = 16)+
        theme(aspect.ratio = 1) +
        labs(x="Number of cells",y="Power",colour="Reads")
    }
  })

  # Slice interactions
  observeEvent(input$slice_click,{
    req(slice_mode()%in%c("row","col"))
    x <- if (slice_mode()=="row")
      reads_seq[which.min(abs(reads_seq-input$slice_click$x))]
    else
      cells_seq[which.min(abs(cells_seq-input$slice_click$x))]
    slice_x(x); updateTextInput(session,"slice_points",value=x)
  })
  observeEvent(input$slice_points,ignoreInit=TRUE,{
    v <- as.numeric(strsplit(input$slice_points,",")[[1]])
    if (!is.na(v[1])) slice_x(v[1])
  })
  observe({
    toggleState("go_slice", length(slice_x())==1)
  })
  observeEvent(input$slice_clear,{
    slice_x(numeric(0)); updateTextInput(session,"slice_points",value="")
  })

  # Slice -> Per-pair
  observeEvent(input$go_slice,{
    req(length(slice_x())==1)
    if (slice_mode()=="row")
      sel$tiles <- expand.grid(row=sel$idx,
                               col=match(slice_x(),reads_seq))
    else
      sel$tiles <- expand.grid(row=match(slice_x(),cells_seq),
                               col=sel$idx)
    sel$type <- "tile"
    updateTabsetPanel(session,"tabs",selected="Per-pair power")
  })

  # Per-pair title & plots
  output$pair_title <- renderText({
    "Per-pair power"
  })
  make_curve <- function(x, base, fun, label)
    data.frame(x=x, Power=fun(x, base), label=label)

  output$pp_combined <- renderPlot({
    req(planned(), is_sel("tile"))

    # First plot: Power vs TPM
    tpm  <- seq(0, 10, length.out = 100)
    dfs1 <- Map(function(r, c){
      base <- log10(cells_seq[r] * reads_seq[c])
      make_curve(tpm, base,
                 function(v, b) 1 - exp(-((v/5)^1.8) * (b/max_pow)^10 * 5000),
                 sprintf("%d × %d", cells_seq[r], reads_seq[c]))
    }, sel$tiles$row, sel$tiles$col)
    p1 <- ggplot(do.call(rbind, dfs1), aes(x, Power, colour = label)) +
      geom_line() +
      geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
      theme_bw(base_size = 16) +
      theme(aspect.ratio = 1) +
      labs(x = "TPM", y = "Power", colour = "Design (cells × reads/cell)")

    # Second plot: Power vs Fold‐change
    fc   <- seq(1, 5, length.out = 100)
    fc   <- seq(0, 50, length.out = 100)
    dfs2 <- Map(function(r, c){
      base <- log10(cells_seq[r] * reads_seq[c])
      make_curve(fc, base,
                 function(v, b) 1 - exp(-(((v*4/50)/2)^1.8) * (b/max_pow)^8 * 2000),
                 sprintf("%d × %d", cells_seq[r], reads_seq[c]))
    }, sel$tiles$row, sel$tiles$col)
    p2 <- ggplot(do.call(rbind, dfs2), aes(x, Power, colour = label)) +
      geom_line() +
      geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
      theme_bw(base_size = 16) +
      theme(aspect.ratio = 1) +
      labs(x = "Fold-change (percent)",
           y = "Power",
           colour = "Design (cells × reads/cell)")

    # Combine with shared legend below
    (p1 + p2) +
      plot_layout(ncol = 2, guides = "collect") &
      theme(legend.position = "bottom")
  })
}

shinyApp(ui, server)
