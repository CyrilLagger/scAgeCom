
source("tab_tissue_specific.R")
#source("tab_global.R")
source("tab_LR6db.R")

ui_scAgeCom <- fluidPage(
  navbarPage(
    title = "A Murine ATLAS of Age-related Intercellular Communication Changes.",
    tab_tissue_specific,
    #tab_global,
    tab_LR6db
  )
)

