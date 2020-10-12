####################################################
##
## Project: scAgeCom
## 
## September 2020
##
## anais.equey@etu.univ-amu.fr 
## cyril.lagger@liverpool.ac.uk 
##
## Script to explore some results of scAgeCom
##
####################################################
##

## 1.  Load libraries ####

library(scDiffCom)
library(data.table)
library(ggplot2)

## 2. Load the results with all cell-cell interactions (CCI) ####

CCI_all <- readRDS("../data_scAgeCom/analysis/analysis_4_data_diffcom_filter.rds") 
#note change the path to where the rds files is stored on your computer

#results are stored as a list (one item per dataset)
#for each dataset the results are in a data.table with one CCI per row
CCI_tmsFACS <- CCI_all$tms_facs
CCI_tmsDroplet <- CCI_all$tms_droplet
CCI_calico <- CCI_all$calico

#we don't need to keep all the columns here
cols_to_keep <- c("TISSUE", "L_CELLTYPE", "R_CELLTYPE", "LR_NAME", paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3),
                  "BH_PVAL_DIFF", "LOGFC", "LOGFC_ABS", "CASE_TYPE"#,
                  #"LR_SCORE_young", "LR_SCORE_old", "LR_DETECTED_young", "LR_DETECTED_old"
)
CCI_tmsFACS <- CCI_tmsFACS[, cols_to_keep, with = FALSE]
CCI_tmsDroplet <- CCI_tmsDroplet[, cols_to_keep, with = FALSE]
CCI_calico <- CCI_calico[, cols_to_keep, with = FALSE]

## 3. Look at how each CCI are distributed among the various "cases"
table(CCI_tmsFACS$CASE_TYPE)
table(CCI_tmsDroplet$CASE_TYPE)
table(CCI_calico$CASE_TYPE)

#volcano plot of pvalues and logfc, optional, uncomment if you want to see it
#ggplot(CCI_tmsFACS, aes(x = LOGFC, y = - log10(BH_PVAL_DIFF))) + geom_point()
#ggplot(CCI_tmsDroplet, aes(x = LOGFC, y = - log10(BH_PVAL_DIFF))) + geom_point()
#ggplot(CCI_calico, aes(x = LOGFC, y = - log10(BH_PVAL_DIFF))) + geom_point()

## 4. Filtering: we only keep the CCIs that change with age
CCI_tmsFACS_age <- CCI_tmsFACS[CASE_TYPE %in% c("FTTU", "TFTD", "TTTD", "TTTU")]
CCI_tmsDroplet_age <- CCI_tmsDroplet[CASE_TYPE %in% c("FTTU", "TFTD", "TTTD", "TTTU")]
CCI_calico_age <- CCI_calico[CASE_TYPE %in% c("FTTU", "TFTD", "TTTD", "TTTU")]

#test <- CCI_tmsFACS_age[L_NCELLS_old >= 20 & R_NCELLS_old >= 20 & L_NCELLS_young >= 20 & L_NCELLS_young >= 20]

## 5. Look at the top up and down regulated CCI

#note this top-X is global, namely is not done in functions of tissues or cell-types

n_top <- 1000

top_LR_up_FACS <- unique(CCI_tmsFACS_age[order(-LOGFC)][1:n_top]$LR_NAME)
top_LR_up_FACS

top_genes_up_FACS <- sort(unique(unlist(CCI_tmsFACS_age[order(-LOGFC)][1:n_top, c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3) )])))
top_genes_up_FACS

top_LR_down_FACS <- unique(CCI_tmsFACS_age[order(LOGFC)][1:n_top]$LR_NAME)
top_LR_down_FACS

top_genes_down_FACS <- sort(unique(unlist(CCI_tmsFACS_age[order(LOGFC)][1:n_top, c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3) )])))
top_genes_down_FACS


## 6. Look at the results of the global ORA (overrepresentation analysis)

ORA_all <- readRDS("../data_scAgeCom/analysis/analysis_4_data_ora.rds")
ORA_tmsFACS <- ORA_all$tms_facs

top_ORA_up_FACS <- ORA_tmsFACS[Tissue == "All" & Category == "LR_NAME" & pval_UP <= 0.05 & OR_UP >= 1][order(pval_UP)]$Value
top_ORA_up_FACS

top_ORA_down_FACS <- ORA_tmsFACS[Tissue == "All" & Category == "LR_NAME" & pval_DOWN <= 0.05 & OR_DOWN >= 1][order(pval_DOWN)]$Value
top_ORA_down_FACS


