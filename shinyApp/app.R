####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - December 2020
##
####################################################
##

## Source code ####
source("global_scAgeCom.R", local = TRUE)
source("ui_scAgeCom.R", local = TRUE)
source("server_scAgeCom.R", local = TRUE)

## Options ####
#options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))
#options(htmlwidgets.TOJSON_ARGS = NULL)
#options("DT.TOJSON_ARGS" = NULL)

## Main shinyApp call ####

shinyApp(
  ui = ui_scAgeCom,
  server = server_scAgeCom
)