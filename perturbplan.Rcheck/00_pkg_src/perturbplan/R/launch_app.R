#' Launch the PerturbPlan Shiny App
#'
#' @return None
#' @export
launch_app <- function() {
  app_file <- system.file("shiny", "app.R", package = "perturbplan")
  if (app_file == "") {
    stop("Could not find app.R. Try reinstalling the `perturbplan` package.", call. = FALSE)
  }
  shiny::runApp(app_file)
}

