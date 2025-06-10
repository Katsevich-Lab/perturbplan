# Header UI module for PerturbPlan dashboard

create_header <- function() {
  dashboardHeader(
    title = "PerturbPlan",
    tags$li(
      class = "dropdown",
      style = "float: right; margin-right: 20px;",
      conditionalPanel(
        condition = "!output.need_plan",
        downloadButton(
          "download_results",
          "Download Results",
          class = "btn-primary"
        )
      )
    )
  )
}