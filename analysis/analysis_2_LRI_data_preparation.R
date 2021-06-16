####################################################
##
## Project: scAgeCom
##
## Last update - June 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## LRI data preparation
##
####################################################
##

## Libraries ####

library(scDiffCom)
library(data.table)
library(openxlsx)
library(ggplot2)
library(kableExtra)

## Dataset Overview ####

# The collections of ligand-receptor interactions (for both human and
# mouse) are built and stored in the package scDiffCom. See e.g. the
# scDiffCom vignette or the scDiffCom file R/utils_LRI.R

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

LRI_mouse_curated[
  ,
  Type := ifelse(
    !is.na(`Ligand (2)`) | !is.na(`Receptor (2)`),
    "Complex",
    "Simple"
  )
]
LRI_mouse_curated[, c(LRI_DATABASES) := lapply(LRI_DATABASES, function(i) {
  ifelse(grepl(i, `Database(s) of Origin`), TRUE, FALSE)
})]

## Save results for scAgeComShiny ####

data_2_LRI_data_preparation <- list(
  LRI_mouse_curated = LRI_mouse_curated,
  LRI_DATABASES = LRI_DATABASES
)

#saveRDS(
#  data_2_LRI_data_preparation,
#  "../data_scAgeCom/analysis/outputs_data/data_2_LRI_data_preparation.rds"
#)

## Prepare Supplemental Data 1 for the manuscript ####

LRI_human_to_xlsx <- scDiffCom::LRI_human[c(1,2,3)]
names(LRI_human_to_xlsx) <- paste0(names(LRI_human_to_xlsx), "_human")
LRI_mouse_to_xlsx <- scDiffCom::LRI_mouse[c(1,2,3)]
names(LRI_mouse_to_xlsx) <- paste0(names(LRI_mouse_to_xlsx), "_mouse")

#openxlsx::write.xlsx(
#  c(LRI_mouse_to_xlsx, LRI_human_to_xlsx),
#  file = "../data_scAgeCom/analysis/outputs_data/Supplemental_Data_1.xlsx"
#)

## Prepare Figure 1.a for the manuscript ####

LRI_upset_mouse <- ComplexUpset::upset(
  as.data.frame(LRI_mouse_curated),
  LRI_DATABASES,
  name = "Database",
  base_annotations = list(
    'Intersection size' = ComplexUpset::intersection_size(
      mapping = ggplot2::aes(fill = Type),
      counts = TRUE,
      bar_number_threshold = 100,
      text = list(size = 8)
    #) + scale_fill_grey(
    #  start = 0.2, end = 0.6
    ) + scale_fill_manual(
      values = c("purple", "coral")
    )
  ),
  themes = ComplexUpset::upset_default_themes(
    text = ggplot2::element_text(size = 26)
  ),
  min_size = 40
) + ggplot2::ggtitle(
  "Number of curated mouse ligand-receptor interactions"
)
LRI_upset_mouse
#manul save: 2000x1200
