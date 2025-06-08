# Modular server logic for PerturbPlan dashboard
# Load all server modules

source("server/server_power.R")
source("server/server_selection.R")
source("server/server_plots.R")
source("server/server_curves.R")

# Main server composition
server <- function(input, output, session) {
  # Initialize core power calculation module
  power_data <- create_power_server(input, output, session)
  
  # Initialize selection handling module  
  selection_data <- create_selection_server(input, output, session, power_data)
  
  # Initialize plotting module
  create_plots_server(input, output, session, power_data, selection_data)
  
  # Initialize curves and download module
  create_curves_server(input, output, session, power_data, selection_data)
}