# Header UI module for PerturbPlan dashboard

create_header <- function() {
  dashboardHeader(
    title = "PerturbPlan",
    tags$li(
      class = "dropdown",
      style = "margin: 8px 10px 0 0;",
      conditionalPanel(
        condition = "!output.need_plan",
        downloadButton(
          "download_results",
          "Download Results",
          class = "btn-primary",
          style = "background-color: white; color: #3c8dbc; border: 1px solid #3c8dbc; margin-top: 7px;"
        )
      )
    )
  )
}