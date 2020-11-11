####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - November 2020
##
####################################################
##

## Load libraries ####

library(data.table)
library(scDiffCom)

## Specify data directory ####
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Load scDiffCom (light) results ####
DATASETS_light <- readRDS(paste0(dir_data_analysis, "a4_data_results_all_cases.rds"))

names(DATASETS_light)
names(DATASETS_light) <- c("Calico Data", "TMS Droplet Data", "TMS FACS Data")

## Process data for tissue-specific analysis ####

# Add columns and change colnames

DATASETS_light <- lapply(
  DATASETS_light,
  function(dataset) {
    lapply(
      dataset,
      function(tiss) {
        cci_table_filtered <- scDiffCom:::get_cci_table_filtered(tiss)
        cci_table_filtered[, c("LOG2FC") := list(LOGFC*log2(exp(1)))]
        setnames(
          cci_table_filtered,
          old = c("L_CELLTYPE", "R_CELLTYPE", "BH_PVAL_DIFF", "LR_NAME"),
          new = c("Emitter Cell Type", "Receiver Cell Type", "Adj. P-Value", "Ligand-Receptor Genes")
        )
        new_object <- scDiffCom:::set_cci_table_filtered(tiss, cci_table_filtered)
        ora_tables <- scDiffCom:::get_ora_tables(new_object)
        ora_tables <- ora_tables[c("GO", "LR_NAME", "LR_CELLTYPE", "LR_CELL_FAMILY")]
        names(ora_tables) <- c("GO Terms", "Ligand-Receptor Genes", "Cell Types", "Cell Families")
        ora_tables <- lapply(
          ora_tables,
          function(ORA_dt) {
            setnames(
              ORA_dt,
              old = c("OR_UP", "pval_adjusted_UP", "OR_DOWN", "pval_adjusted_DOWN",
                      "OR_FLAT", "pval_adjusted_FLAT"),
              new = c("Odds Ratio Up", "Adj. P-Value Up", "Odds Ratio Down", "Adj. P-Value Down",
                      "Odds Ratio Stable", "Adj. P-Value Stable")
            )
            if("Value_NAME" %in% colnames(ORA_dt)) {
              setnames(
                ORA_dt,
                old = c("Value", "Value_NAME"),
                new = c("Value_ID", "Value")
              )
            }
            return(ORA_dt)
          }
        )
        new_object <- scDiffCom:::set_ora_tables(new_object, ora_tables)
        return(new_object)
      }
    )
  }
)

saveRDS(DATASETS_light, "shinyApp/data/scdiffcom_objects_shiny.rds")

