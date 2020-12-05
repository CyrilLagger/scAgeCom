####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - November 2020
##
####################################################
##

## Source code ####
source("global_scAgeCom.R", local = TRUE)
source("ui_scAgeCom.R", local = TRUE)
source("server_scAgeCom.R", local = TRUE)

## Options ####
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))
## Main shinyApp call ####

shinyApp(
  ui = ui_scAgeCom,
  server = server_scAgeCom
)