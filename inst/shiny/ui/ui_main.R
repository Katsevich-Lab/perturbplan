# Modular UI definition for PerturbPlan dashboard
# Load all UI modules

source("ui/ui_header.R")
source("ui/ui_sidebar.R") 
source("ui/ui_styles.R")
source("ui/ui_tabs.R")

# Main UI composition
ui <- dashboardPage(
  header = create_header(),
  sidebar = create_sidebar(),
  body = dashboardBody(
    create_styles(),
    create_main_tabs()
  )
)