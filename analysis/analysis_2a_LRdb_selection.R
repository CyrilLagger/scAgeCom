####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - September 2020
##
## Investigate the 6 LR databases and do a final
## choice regading the curated interations.
##
####################################################
##

## Libraries ####

library(scDiffCom)
library(UpSetR)
library(ComplexUpset)
library(ggplot2)


## Analysis data path ####
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Load the LR database ####
LR6db_full <- scDiffCom::LR6db$LR6db_all
LR6db_curated <- scDiffCom::LR6db$LR6db_curated
LR6db_source <- scDiffCom::LR6db$LR6db_source

category_sources <- c("PPI", "IUPHAR", "KEGG", "PMID", "CPDB", "ramilowski",
                      "cellsignal.com", "HPMR", "HPRD", "reactome")

category_DBs <- names(LR6db_source)

## Save the curated database as a csv file ####
write.csv(LR6db_curated, paste0(dir_data_analysis, "a2_LR6db_curated.csv"), row.names = FALSE)

## Clean up the SOURCE of each interaction ####
LR_rm_sctensor <- c("SWISSPROT_STRING", "TREMBL_STRING")
LR_rm_nichenet <- c("ppi_bidir_bidir", "ppi_bidir_bidir_go", "ppi_bidir_r",
                    "ppi_bidir_r_go", "ppi_l_bidir", "ppi_l_bidir_go",
                    "ppi_lr", "ppi_lr_go")

LR6db_full[, SOURCE_CLEAN := gsub(paste0(c(LR_rm_nichenet, LR_rm_sctensor, "uniprot"), collapse = "|"), "PPI", SOURCE)]
LR6db_full[, SOURCE_CLEAN := gsub("pharmacology", "IUPHAR", SOURCE_CLEAN)]
LR6db_full[, SOURCE_CLEAN := gsub("kegg", "KEGG", SOURCE_CLEAN)]
LR6db_full[, SOURCE_CLEAN := gsub("fantom5", "ramilowski", SOURCE_CLEAN)]
LR6db_full[, SOURCE_no_digit := gsub(" ", "", gsub('[[:digit:]]+', '', SOURCE))]
LR6db_full[SOURCE_no_digit %in% c("", "; ", " ;", ";") | nchar(SOURCE_no_digit) <= 2, SOURCE_CLEAN := paste0("PMID:", SOURCE_CLEAN)]
sum(!grepl(paste0(category_sources, collapse = "|"), LR6db_full$SOURCE_CAT))

LR6db_curated[, SOURCE_CLEAN := gsub(paste0(c(LR_rm_nichenet, LR_rm_sctensor, "uniprot"), collapse = "|"), "PPI", SOURCE)]
LR6db_curated[, SOURCE_CLEAN := gsub("pharmacology", "IUPHAR", SOURCE_CLEAN)]
LR6db_curated[, SOURCE_CLEAN := gsub("kegg", "KEGG", SOURCE_CLEAN)]
LR6db_curated[, SOURCE_CLEAN := gsub("fantom5", "ramilowski", SOURCE_CLEAN)]

LR6db_curated[, SOURCE_no_digit := gsub(" ", "", gsub('[[:digit:]]+', '', SOURCE))]
LR6db_curated[SOURCE_no_digit %in% c("", "; ", " ;", ";") | nchar(SOURCE_no_digit) <= 2, SOURCE_CLEAN := paste0("PMID:", SOURCE_CLEAN)]


sum(!grepl(paste0(category_sources, collapse = "|"), LR6db_curated$SOURCE_CAT))

## Produce Upsetplot of full and curated LR interactions based on their source and database of origin ####

#produce correct format for upset function
LR_full_upset <- LR6db_full
LR_full_upset[, c(category_sources) := lapply(category_sources, function(i) {
  ifelse(grepl(i, SOURCE_CLEAN), 1, 0)
})]
LR_full_upset[, c(category_DBs) := lapply(category_DBs, function(i) {
  ifelse(grepl(i, DATABASE), 1, 0)
})]
LR_full_upset[, complex := !is.na(LIGAND_2) | !is.na(RECEPTOR_2)]

LR_curated_upset <- LR6db_curated
LR_curated_upset[, c(category_sources) := lapply(category_sources, function(i) {
  ifelse(grepl(i, SOURCE_CLEAN), 1, 0)
})]
LR_curated_upset[, c(category_DBs) := lapply(category_DBs, function(i) {
  ifelse(grepl(i, DATABASE), 1, 0)
})]
LR_curated_upset[, complex := !is.na(LIGAND_2) | !is.na(RECEPTOR_2)]

#upset plot of all interactions
plot_upset_full_dbs <- ComplexUpset::upset(
  LR_full_upset,
  category_DBs,
  min_size = 50,
  base_annotations=list(
    'Intersection size'=intersection_size(
      bar_number_threshold = 100
    )
  )
)
plot_upset_full_dbs
ggsave(filename = paste0(dir_data_analysis, "a2_plot_upset_full_dbs.png"),
       plot = plot_upset_full_dbs, scale = 1.5)

plot_upset_full_sources <- ComplexUpset::upset(
  LR_full_upset,
  category_sources,
  min_size = 50,
  base_annotations=list(
    'Intersection size'=intersection_size(
      bar_number_threshold = 100
    )
  )
)
plot_upset_full_sources
ggsave(filename = paste0(dir_data_analysis, "a2_plot_upset_full_sources.png"),
       plot = plot_upset_full_sources, scale = 1.5)

plot_upset_curated_dbs <- ComplexUpset::upset(
  LR_curated_upset,
  category_DBs,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      aes=aes(fill=complex),
      text = list(position = position_stack(vjust = 0.0)),
      bar_number_threshold = 100
    )
  ),
  min_size = 20
)
plot_upset_curated_dbs
ggsave(filename = paste0(dir_data_analysis, "a2_plot_upset_curated_dbs.png"),
       plot = plot_upset_curated_dbs, scale = 1.5)

plot_upset_curated_sources <- ComplexUpset::upset(
  LR_curated_upset,
  category_sources,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      aes=aes(fill=complex),
      text = list(position = position_stack(vjust = 0.0)),
      bar_number_threshold = 100
    )
  ),
  min_size = 30
)
plot_upset_curated_sources
ggsave(filename = paste0(dir_data_analysis, "a2_plot_upset_curated_sources.png"),
       plot = plot_upset_curated_sources, scale = 1.5)

#UpSetR::upset(LR_full_upset[, ..category_DBs], nsets = 6, order.by = "freq", nintersects = 35)
#UpSetR::upset(LR_full_upset[, ..category_sources], nsets = 10, order.by = "freq", nintersects = 35)
#UpSetR::upset(LR_curated_upset[, ..category_DBs], nsets = 6, order.by = "freq", nintersects = 35)
#UpSetR::upset(LR_curated_upset[, ..category_sources], nsets = 10, order.by = "freq", nintersects = 35)

## Check orthology confidence of all interactions 
ftable(LR6db_curated$RECEPTOR_1_CONF, LR6db_curated$LIGAND_1_CONF)






