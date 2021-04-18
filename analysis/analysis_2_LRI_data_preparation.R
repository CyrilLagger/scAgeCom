####################################################
##
## Project: scAgeCom
##
## Last update - April 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## LRI data preparation
##
####################################################
##

## Dataset Overview ####

# The collections of ligand-receptor interactions (for both human and
# mouse) are built and stored in the package scDiffCom. See e.g. the
# scDiffCom vignette or the scDiffCom file R/utils_LRI.R

## Libraries ####

library(scDiffCom)
library(data.table)
library(ggplot2)
library(ComplexUpset)

## Prepare mouse LRI database for scAgeComShiny ####

LRI_mouse_curated <- copy(scDiffCom::LRI_mouse$LRI_curated)
LRI_mouse_curated[
  ,
  SOURCE := gsub("SCT:SingleCellSignalR|SCT:CellPhoneDB|SCT:DLRP", "", SOURCE)
  ]
LRI_mouse_curated[, SOURCE := gsub(";;", ";", SOURCE)]
LRI_mouse_curated[, SOURCE := sub(";$", "", SOURCE)]
LRI_mouse_curated[, SOURCE := sub("^;", "", SOURCE)]
LRI_cols_to_keep <- c(
  "LIGAND_1", "LIGAND_2",
  "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
  "DATABASE", "SOURCE"
)
LRI_mouse_curated <- LRI_mouse_curated[, LRI_cols_to_keep, with = FALSE]
setnames(
  LRI_mouse_curated,
  old = LRI_cols_to_keep,
  new = c(
    "Ligand (1)", "Ligand (2)",
    "Receptor (1)", "Receptor (2)", "Receptor (3)",
    "Database(s) of Origin", "Retrieved Sources"
  )
)

LRI_DATABASES <- sort(
  unique(
    unlist(
      strsplit(
        LRI_mouse_curated$`Database(s) of Origin`,
        ";")
    )
  )
)

LRI_mouse_curated[, COMPLEX := !is.na(`Ligand (2)`) | !is.na(`Receptor (2)`)]
LRI_mouse_curated[, c(LRI_DATABASES) := lapply(LRI_DATABASES, function(i) {
  ifelse(grepl(i, `Database(s) of Origin`), TRUE, FALSE)
})]

## Prepare LRI utility functions for scAgeComShiny ####

plot_lri_upset <- function(
  LRI_table,
  groups,
  min_size
) {
  p <- ComplexUpset::upset(
    as.data.frame(LRI_table),
    groups,
    base_annotations = list(
      'Intersection size' = intersection_size(
        mapping = aes(fill = COMPLEX),
        counts = TRUE,
        bar_number_threshold = 100
      )
    ),
    themes = upset_default_themes(text = element_text(size = 16)),
    min_size = min_size
  ) +
    ggtitle("Number of Ligand-Receptor Interactions")
  return(p)
}

build_LRI_display <- function(
  LRI_table,
  LRI_database
) {
  dt <- LRI_table[
    apply(
      sapply(
        LRI_database,
        function(i) {
          grepl(i, `Database(s) of Origin`)
        }
      ),
      MARGIN = 1,
      any
    )
  ]
  DT::datatable(
    data = dt[, 1:7],
    class = "display compact",
    options = list(
      pageLength = 10,
      columnDefs = list(
        list(
          targets = c(6,7),
          render = htmlwidgets::JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 20 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
            "}")
        )
      )
    ),
    caption = tags$caption(
      style = 'caption-side: top; text-align: center; color:black; font-size:130% ;',
      "Table of Mouse Ligand-Receptor Interactions"
    )
  )
}

## Prepare LRI figures for manuscript ####

LRI_mouse_upset <- plot_lri_upset(
  LRI_table = LRI_mouse_curated,
  groups = LRI_DATABASES,
  min_size = 40
)
LRI_mouse_upset

## save all results ####

data_2_LRI_data_preparation <- list(
  LRI_mouse_curated = LRI_mouse_curated,
  LRI_DATABASES = LRI_DATABASES,
  build_LRI_display = build_LRI_display,
  plot_lri_upset = plot_lri_upset
)

saveRDS(
  data_2_LRI_data_preparation,
  "../data_scAgeCom/analysis/outputs_data/data_2_LRI_data_preparation.rds"
)

