####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - November 2020
##
## Perfom ORA analysis and result interpretation
## per tissue.
##
####################################################
##

## Libraries ####
library(scDiffCom)
library(data.table)
library(ggplot2)

## Specify data directory ####
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Load scDiffCom (light) results ####
DATASETS_light <- readRDS(paste0(dir_data_analysis, "a4_data_results_all_cases.rds"))

names(DATASETS_light)
DATASETS_light <- DATASETS_light[c(1,4,8)]

## Add specific terms to cell-types and genes ####

#cell-type families
cell_types_dt <- setDT(read.csv(paste0(dir_data_analysis, "scDiffCom_cell_types.csv"), stringsAsFactors = FALSE))

#genage longevity genes
genage_mouse <- setDT(read.csv(paste0(dir_data_analysis, "genage_mouse.tsv"), sep = "\t", header = TRUE))

#add the terms to the data.table
DATASETS_light <- lapply(
  DATASETS_light,
  function(dataset) {
    lapply(
      dataset,
      function(tiss) {
        dt <- scDiffCom:::get_cci_table_filtered(tiss)
        if(identical(dt, list())) return(tiss)
        dt[cell_types_dt, on = "L_CELLTYPE==scDiffCom.cell.type", L_CELL_FAMILY := i.Family...broad]
        dt[cell_types_dt, on = "R_CELLTYPE==scDiffCom.cell.type", R_CELL_FAMILY := i.Family...broad]
        dt[, LR_CELL_FAMILY := paste(L_CELL_FAMILY, R_CELL_FAMILY, sep = "_")]
        dt[,GENAGE := ifelse(
          LIGAND_1 %in% genage_mouse$Gene.Symbol | LIGAND_2 %in% genage_mouse$Gene.Symbol |
            RECEPTOR_1 %in% genage_mouse$Gene.Symbol | RECEPTOR_2 %in% genage_mouse$Gene.Symbol,
          "longevity_associated",
          "not_longevity_associated"
        )
        ]
        tiss <- scDiffCom:::set_cci_table_filtered(tiss, dt)
      }
    )
  }
)

## Redo ORA with new categories ####

DATASETS_light <- lapply(
  DATASETS_light,
  function(dataset) {
    lapply(
      dataset,
      function(tiss) {
        scDiffCom::run_ORA(
          object = tiss,
          categories = c("LR_CELLTYPE", "LR_CELL_FAMILY", "LR_NAME", "GO", "GENAGE"),
          overwrite = TRUE,
          logfc_threshold = log(1.2)
        )
      }
    )
  }
)

## Save results with new ORA ####

saveRDS(DATASETS_light, paste0(dir_data_analysis, "a4_data_results_all_cases.rds"))


