
source("tab_description.R")
source("tab_tissue_specific.R")
source("tab_combined_analysis.R")
#source("tab_LRdb.R")

ui_scAgeCom <- fluidPage(
  titlePanel(
    title = "A Murine ATLAS of Age-related Intercellular Communication Changes."
    ),
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  ),
  navbarPage(
    title = "",
    tab_description,
    tab_tissue_specific,
    tab_combined_analysis#,
    #tab_LRdb
  )
)

