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
  min_size = 50
)
ggsave(filename = paste0(dir_data_analysis, "analysis_2_plot_upset_full_dbs.png"),
       plot = plot_upset_full_dbs, scale = 2)

plot_upset_full_sources <- ComplexUpset::upset(
  LR_full_upset,
  category_sources,
  min_size = 50
)
ggsave(filename = paste0(dir_data_analysis, "analysis_2_plot_upset_full_sources.png"),
       plot = plot_upset_full_sources, scale = 2)

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
ggsave(filename = paste0(dir_data_analysis, "analysis_2_plot_upset_curated_dbs.png"),
       plot = plot_upset_curated_dbs, scale = 2)

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
ggsave(filename = paste0(dir_data_analysis, "analysis_2_plot_upset_curated_sources.png"),
       plot = plot_upset_curated_sources, scale = 2)

#UpSetR::upset(LR_full_upset[, ..category_DBs], nsets = 6, order.by = "freq", nintersects = 35)
#UpSetR::upset(LR_full_upset[, ..category_sources], nsets = 10, order.by = "freq", nintersects = 35)
#UpSetR::upset(LR_curated_upset[, ..category_DBs], nsets = 6, order.by = "freq", nintersects = 35)
#UpSetR::upset(LR_curated_upset[, ..category_sources], nsets = 10, order.by = "freq", nintersects = 35)





#look at orthology
#OK:
table(LRdb_filtered$LIGAND_2_CONF)
table(LRdb_filtered$RECEPTOR_3_CONF)
# a few 0
table(LRdb_filtered$LIGAND_1_CONF)
table(LRdb_filtered$RECEPTOR_1_CONF)
table(LRdb_filtered$RECEPTOR_2_CONF)
ftable(LRdb_filtered$LIGAND_1_CONF, LRdb_filtered$RECEPTOR_1_CONF)
ftable(LRdb_filtered$LIGAND_1_CONF, LRdb_filtered$RECEPTOR_1_CONF, LRdb_filtered$RECEPTOR_2_CONF)



library(data.table)
fwrite(LRdb_filtered, file = "../../../../../test.csv")

#annotation with GO
library(clusterProfiler)
library(org.Mm.eg.db)

LR_genes <- unique(unlist(LR6db$LR6db_curated[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
LR_genes <- LR_genes[!is.na(LR_genes)]


LR_ggo_bc_l3 <- groupGO(
  gene = LR_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  level = 3,
  readable = FALSE
)@result

LR_ggo_bc_l2 <- groupGO(
  gene = LR_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  level = 2,
  readable = FALSE
)@result

LR_ggo_bc_l4 <- groupGO(
  gene = LR_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  level = 4,
  readable = FALSE
)@result

####

library(biomaRt)

mart <- biomaRt::useMart(
  "ensembl",
  dataset = "mmusculus_gene_ensembl"
)
LR_genes_go <- biomaRt::getBM(
  attributes = c(
    "mgi_symbol",
    "go_id",
    "name_1006", 
    "namespace_1003"
  ),
  filters = "mgi_symbol",
  mart = mart,
  values = LR_genes
)
setDT(LR_genes_go)

LR_genes_go_comb <- LR_genes_go[, list(text = paste(name_1006, collapse=",")), by = mgi_symbol]

LR_genes_go_comb2 <- LR_genes_go[, list(text = paste(mgi_symbol, collapse=",")), by = name_1006]

library(GO.db)
Term("GO:0016021")
t <- c(GOBPOFFSPRING[["GO:0042110"]], "GO:0042110")
Term("GO:0042110")
Term(t[2])

library(data.table)
setDT(LR_genes_go)

LR_genes_go[, count := .N, by = go_id]


LR_genes_go[ go_id == "GO:0050789"]

test <- LRall
test2 <- test[ cpdb == TRUE | scsr == TRUE]

#####
upset_list_full_sources <- sapply(category_sources, function(i) {
  LR6db_full[grepl(i, SOURCE_CLEAN)]$LR_SORTED
})
upset_list_full_DBs <- sapply(category_DBs, function(i) {
  LR6db_full[grepl(i, DATABASE)]$LR_SORTED
})

upset_list_curated_sources <- sapply(category_sources, function(i) {
  LR6db_curated[grepl(i, SOURCE_CLEAN)]$LR_SORTED
})
upset_list_curated_DBs <- sapply(category_DBs, function(i) {
  LR6db_curated[grepl(i, DATABASE)]$LR_SORTED
})
