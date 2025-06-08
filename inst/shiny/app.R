# Main Shiny app entry point
# Loads modular components and launches the app

library(shiny)
library(shinydashboard)
library(shinyBS)
library(shinyjs)
library(ggplot2)
library(patchwork)

# Source modular components
source("ui/ui_main.R")
source("server/server_main.R")

# Launch the app
shinyApp(ui = ui, server = server)