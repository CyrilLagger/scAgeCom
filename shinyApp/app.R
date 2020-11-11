####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - November 2020
##
####################################################
##

## Libraries ####

library(shiny)
library(shinyWidgets)
library(DT)
library(data.table)
library(ggplot2)
library(scDiffCom)
library(ComplexUpset)

## Source code ####

source("ui_scAgeCom.R", local = TRUE)
source("server_scAgeCom.R", local = TRUE)

## Global Data ####

LR6db_curated <- readRDS("data/LR6db_curated.rds")
cols_to_show_LR6db <- c(
  "LIGAND_1", "LIGAND_2",
  "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
  "DATABASE", "SOURCE"
)

DATASETS_light <- readRDS("data/scdiffcom_objects_shiny.rds")
cols_to_show_DATA <- c("Ligand-Receptor Genes", "Emitter Cell Type", "Receiver Cell Type",
                            "LOG2FC", "Adj. P-Value", "REGULATION")
cols_numeric_DATA <- c("LOG2FC", "Adj. P-Value")

## Global functions ####

show_DT <- function(
  data,
  cols_to_show,
  cols_numeric
) {
  res <- DT::datatable(
    data = data[, cols_to_show, with = FALSE],
    options = list(
      pageLength = 10
    )
  ) 
  if(!is.null(cols_numeric)) {
    res <- DT::formatSignif(
      table = res,
      columns = cols_numeric,
      digits = 3
    )
  }
  return(res)
}

show_volcano <- function(
  data
) {
  ggplot(data, aes(
    x = LOG2FC,
    y = -log10(`Adj. P-Value` + 1E-4
    ),
    color = ifelse(
      `Adj. P-Value` <= 0.05 & LOG2FC >= log2(1.2),
      "red",
      ifelse(
        `Adj. P-Value` <= 0.05 & LOG2FC <= -log2(1.2),
        "blue",
        "green"
      )
      )
  )) +
    geom_point() +
    scale_color_identity() +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = log2(1.2)) +
    geom_vline(xintercept = -log2(1.2)) +
    xlab(expression(paste(Log[2], "FC"))) +
    ylab(expression(paste(-Log[10], " ", p[BH]))) +
    theme(text=element_text(size=20))
}

## Options ####
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

## Main shinyApp call ####

shinyApp(
  ui = ui_scAgeCom,
  server = server_scAgeCom
)