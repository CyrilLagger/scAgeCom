####################################################
##
## Project: scAgeCom
##
## Last update - December 2020
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## Investigate the 8 LR databases and do a final
## choice regading the curated interations.
##
####################################################
##

## Libraries ####

library(scDiffCom)
library(UpSetR)
library(ComplexUpset)
library(ggplot2)
library(data.table)



test <- LRI_mouse$LRI_curated

test$

## Analysis data path ####
path_directory_analysis <- "../data_scAgeCom/analysis/"

## Load the LR database ####
LRdb <- list(
  human_notcurated = copy(scDiffCom::LRdb_human$LRdb_not_curated),
  human_curated = copy(scDiffCom::LRdb_human$LRdb_curated),
  mouse_notcurated = copy(scDiffCom::LRdb_mouse$LRdb_not_curated),
  mouse_curated = copy(scDiffCom::LRdb_mouse$LRdb_curated)
)

## Add some annotations ####

LRdb_DBS <- sort(unique(unlist(strsplit(LRdb$human_notcurated$DATABASE, ";"))))
LRdb_DBS

LRdb_sources <- c(
  "FANTOM5", "HPMR", "HPRD", "PMID",
  "CellPhoneDB", "KEGG", "IUPHAR", "reactome",
  "cellsignal.com", "PPI"
)

LRdb <- lapply(
  LRdb,
  function(i) {
    i[, COMPLEX := !is.na(LIGAND_2) | !is.na(RECEPTOR_2)]
    i[, c(LRdb_DBS) := lapply(LRdb_DBS, function(i) {
      ifelse(grepl(i, DATABASE), TRUE, FALSE)
    })]
    i[, c(LRdb_sources) := lapply(LRdb_sources, function(i) {
      ifelse(grepl(i, SOURCE), TRUE, FALSE)
    })]
    return(i)
  }
)

## pmid or pmc of the databases (used to determine cross-reference below) ####
LRdb_DBS_pmid <- c("32103204", "30429548", "33147626", "31819264",
                   "26198319", "4525178", "33024107", "7538930",
                   "32196115", "7261168")

## Below we investigate each DB in more datail (see also summary table on google drive) ####

## Study CellChat ####

# retrieve from internal function
cellchat_human <- scDiffCom:::prepare_LR_CellChat(species = "human") 
cellchat_mouse <- scDiffCom:::prepare_LR_CellChat(species = "mouse")

# all evidences are from KEGG or PMID or PMC
cellchat_human[!grepl("KEGG|PMID|PMC", SOURCE)]
cellchat_mouse[!grepl("KEGG|PMID|PMC", SOURCE)]

# check if cross-dependencies with other DBS
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, cellchat_human$SOURCE))
)
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, cellchat_mouse$SOURCE))
)

# number of interactions by types
LRdb$human_notcurated[grepl("CELLCHAT", DATABASE)][,.N, by = "COMPLEX"]
LRdb$mouse_notcurated[grepl("CELLCHAT", DATABASE)][,.N, by = "COMPLEX"]

## Study CellPhoneDB ####

# retrieve from internal function
cpdb_human <- scDiffCom:::prepare_LR_cpdb(species = "human") 
cpdb_mouse <- scDiffCom:::prepare_LR_cpdb(species = "mouse") 

# we have not retrieve the evidence, so we just refer as coming from cpdb
unique(cpdb_human$SOURCE)
unique(cpdb_mouse$SOURCE)

# we don't check for cross-dependency, maybe with FANTOM5?

# number of interactions by types
LRdb$human_notcurated[grepl("CELLPHONEDB", DATABASE)][,.N, by = "COMPLEX"]
LRdb$mouse_notcurated[grepl("CELLPHONEDB", DATABASE)][,.N, by = "COMPLEX"]

## Study CellTalkDB ####

# retrieve from internal function
celltalk_human <- scDiffCom:::prepare_LR_CellTalkDB(species = "human") 
celltalk_mouse <- scDiffCom:::prepare_LR_CellTalkDB(species = "mouse")

# all evidences have a pmid attached
celltalk_human[!grepl("PMID", SOURCE)]
celltalk_mouse[!grepl("PMID", SOURCE)]

# but there are some cross-dependencies
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, celltalk_human$SOURCE))
)
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, celltalk_mouse$SOURCE))
)

# number of interactions by types
LRdb$human_notcurated[grepl("CELLTALK", DATABASE)][,.N, by = "COMPLEX"]
LRdb$mouse_notcurated[grepl("CELLTALK", DATABASE)][,.N, by = "COMPLEX"]

## Study NATMI connectomeDB ####

# retrieve from internal function
natmi_human <- scDiffCom:::prepare_LR_CONNECTOMEDB(species = "human") 
natmi_mouse <- scDiffCom:::prepare_LR_CONNECTOMEDB(species = "mouse")

# all evidences have a pmid attached
natmi_human[!grepl("PMID", SOURCE)]
natmi_mouse[!grepl("PMID", SOURCE)]

# no clear cross-dependencies
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, natmi_human$SOURCE))
)
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, natmi_mouse$SOURCE))
)

# number of interactions by types
LRdb$human_notcurated[grepl("CONNECTOMEDB", DATABASE)][,.N, by = "COMPLEX"]
LRdb$mouse_notcurated[grepl("CONNECTOMEDB", DATABASE)][,.N, by = "COMPLEX"]

## Study ICELLNET ####

# retrieve from internal function
icellnet_human <- scDiffCom:::prepare_LR_ICELLNET(species = "human") 
icellnet_mouse <- scDiffCom:::prepare_LR_ICELLNET(species = "mouse")

# all evidences have a pmid attached
icellnet_human[!grepl("PMID", SOURCE)]
icellnet_mouse[!grepl("PMID", SOURCE)]

# no clear cross-dependencies
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, icellnet_human$SOURCE))
)
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, icellnet_mouse$SOURCE))
)

# number of interactions by types
LRdb$human_notcurated[grepl("ICELLNET", DATABASE)][,.N, by = "COMPLEX"]
LRdb$mouse_notcurated[grepl("ICELLNET", DATABASE)][,.N, by = "COMPLEX"]

## Study NicheNet ####

# retrieve from internal function
nichenet_human <- scDiffCom:::prepare_LR_nichenet(species = "human") 
nichenet_mouse <- scDiffCom:::prepare_LR_nichenet(species = "mouse")

# evidences come from  KEGG, FANTOM5 and PPI
unique(nichenet_human$SOURCE)
table(nichenet_human$SOURCE)
unique(nichenet_mouse$SOURCE)
table(nichenet_mouse$SOURCE)

# number of interactions by types
LRdb$human_notcurated[grepl("NICHENET", DATABASE)][,.N, by = "COMPLEX"]
LRdb$human_curated[grepl("NICHENET", DATABASE)][,.N, by = "COMPLEX"]

LRdb$mouse_notcurated[grepl("ICELLNET", DATABASE)][,.N, by = "COMPLEX"]
LRdb$mouse_curated[grepl("NICHENET", DATABASE)][,.N, by = "COMPLEX"]

## Study SingleCellSignalR ####

# retrieve from internal function
scsr_human <- scDiffCom:::prepare_LR_scsr(species = "human") 
scsr_mouse <- scDiffCom:::prepare_LR_scsr(species = "mouse")

# evidences come from PMID, reactome, PPI and cellsignal
scsr_human[!grepl("PMID|reactome|Ramilowski2015|PPI|cellsignal", SOURCE)]
scsr_mouse[!grepl("PMID|reactome|Ramilowski2015|PPI|cellsignal", SOURCE)]


# no clear cross-dependencies
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, scsr_human$SOURCE))
)
sapply(
  LRdb_DBS_pmid,
  function(i) sum(grepl(i, scsr_mouse$SOURCE))
)

#
scsr_human[!grepl("PMID|reactome|Ramilowski2015|cellsignal", SOURCE)]
scsr_mouse[!grepl("PMID|reactome|Ramilowski2015|cellsignal", SOURCE)]

# number of interactions by types
LRdb$human_notcurated[grepl("SCSR", DATABASE)][,.N, by = "COMPLEX"]
LRdb$human_curated[grepl("SCSR", DATABASE)][,.N, by = "COMPLEX"]

LRdb$mouse_notcurated[grepl("SCSR", DATABASE)][,.N, by = "COMPLEX"]
LRdb$mouse_curated[grepl("SCSR", DATABASE)][,.N, by = "COMPLEX"]

## Study scTensor ####

# retrieve from internal function
sct_human <- scDiffCom:::prepare_LR_scTensor(species = "human") 
sct_mouse <- scDiffCom:::prepare_LR_scTensor(species = "mouse")

# evidences
unique(sct_human$SOURCE)
unique(sct_mouse$SOURCE)

# number of interactions by types
LRdb$human_notcurated[grepl("SCTENSOR", DATABASE)][,.N, by = "COMPLEX"]
LRdb$human_curated[grepl("SCTENSOR", DATABASE)][,.N, by = "COMPLEX"]

LRdb$mouse_notcurated[grepl("SCTENSOR", DATABASE)][,.N, by = "COMPLEX"]
LRdb$mouse_curated[grepl("SCTENSOR", DATABASE)][,.N, by = "COMPLEX"]

## Plot by Database of origin ####

analysis_2_figure_mouse_upset_db <- ComplexUpset::upset(
  LRdb$mouse_curated,
  LRdb_DBS,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      aes=aes(fill=COMPLEX),
      text = list(position = position_stack(vjust = 0.0)),
      bar_number_threshold = 100
    )
  ),
  themes=upset_default_themes(text=element_text(size=20)),
  min_size = 39
)
analysis_2_figure_mouse_upset_db

## Plot by Evidence of origin ####

analysis_2_figure_mouse_upset_evid <- ComplexUpset::upset(
  LRdb$mouse_curated,
  LRdb_sources,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      aes=aes(fill=COMPLEX),
      text = list(position = position_stack(vjust = 0.0)),
      bar_number_threshold = 100
    )
  ),
  themes=upset_default_themes(text=element_text(size=20)),
  min_size = 30
)
analysis_2_figure_mouse_upset_evid

## old code ####
LRdb_mouse_curated[, c(LRdb_SOURCES) := lapply(LRdb_SOURCES, function(i) {
  ifelse(grepl(i, SOURCE), 1, 0)
})]

ComplexUpset::upset(
  LRdb_mouse_curated,
  LRdb_SOURCES,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      aes=aes(fill=COMPLEX),
      text = list(position = position_stack(vjust = 0.0)),
      bar_number_threshold = 100
    )
  ),
  min_size = 20
)


category_sources <- c("PPI", "IUPHAR", "KEGG", "PMID", "CPDB", "ramilowski",
                      "cellsignal.com", "HPMR", "HPRD", "reactome")

category_DBs <- names(LR6db_source)

## Save the curated database as a csv file ####
write.csv(LR6db_curated, paste0(path_directory_analysis, "a2_LR6db_curated.csv"), row.names = FALSE)

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
ggsave(filename = paste0(path_directory_analysis, "a2_plot_upset_full_dbs.png"),
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
ggsave(filename = paste0(path_directory_analysis, "a2_plot_upset_full_sources.png"),
       plot = plot_upset_full_sources, scale = 1.5)


plot_upset_curated_dbs
ggsave(filename = paste0(path_directory_analysis, "a2_plot_upset_curated_dbs.png"),
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
ggsave(filename = paste0(path_directory_analysis, "a2_plot_upset_curated_sources.png"),
       plot = plot_upset_curated_sources, scale = 1.5)

#UpSetR::upset(LR_full_upset[, ..category_DBs], nsets = 6, order.by = "freq", nintersects = 35)
#UpSetR::upset(LR_full_upset[, ..category_sources], nsets = 10, order.by = "freq", nintersects = 35)
#UpSetR::upset(LR_curated_upset[, ..category_DBs], nsets = 6, order.by = "freq", nintersects = 35)
#UpSetR::upset(LR_curated_upset[, ..category_sources], nsets = 10, order.by = "freq", nintersects = 35)

## Check orthology confidence of all interactions 
ftable(LR6db_curated$RECEPTOR_1_CONF, LR6db_curated$LIGAND_1_CONF)

####################################################
##
## Project: scAgeCom
## 
## October 2020
##
## cyril.lagger@liverpool.ac.uk
## anais.equey@etu.univ-amu.fr
##
## Script to annotate the (curated) ligand-receptor
## interactions.
##
####################################################
##
## SUMMARY
## 1 to 3 - we start by listing all the genes we have in our database
## 4 - then we find the go terms associated with each of our genes
##     (based on ensembl)
## 4 - For each go term, we also add its parents (based on Gene Ontology)
## 5 - Each interaction involves several genes. We take both the union and the
##     intersection of their go terms and add the result to separate columns
##
####################################################

#--------- 1. Load libraries --------#

library(scDiffCom)
library(data.table)
library(ontologyIndex)
library(ontoProc)
library(biomaRt)

#--------- 2. Load the curated ligand-receptor interactions --------#

LR_interactions <- scDiffCom::LR6db$LR6db_curated

#--------- 3. Get all unique genes in the data.table --------#

# Ligand genes
L_genes <- unique(unlist(LR_interactions[, c("LIGAND_1", "LIGAND_2")]))
L_genes <- L_genes[!is.na(L_genes)] # remove NA (empty) values
# Receptor genes
R_genes <- unique(unlist(LR_interactions[, c("RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
R_genes <- R_genes[!is.na(R_genes)] # remove NA (empty) values
# All genes (some genes may appear as both ligand and receptors)
LR_genes <- unique(c(L_genes, R_genes))

#--------- 4. Get the Ensembl terms (GO, family or else) --------#

# create an access to ensembl database via biomart (might take a minute or so)
mart <- biomaRt::useMart(
  "ensembl",
  dataset = "mmusculus_gene_ensembl"
)
# retrieve the information we want for our genes of interest
LR_genes_info <- biomaRt::getBM(
  attributes = c(
    "mgi_symbol",
    "go_id",
    "name_1006",
    #"definition_1006",
    "namespace_1003",
    "family",
    "family_description"
    #"goslim_goa_accession",
    #"goslim_goa_description"
  ),
  filters = "mgi_symbol",
  mart = mart,
  values = LR_genes
)
setDT(LR_genes_info)

#--------- 5. adding ancestors from Gene Ontology --------#

# Load go terms from gene ontology
onto_go_terms = getGeneOnto()
go_names <- onto_go_terms$name
LR_genes_go <- sapply(
  LR_genes,
  function(gene) {
    temp_go <- unique(LR_genes_info[mgi_symbol == gene]$name_1006)
    get_ancestors(onto_go_terms, names(go_names[go_names %in% temp_go]))
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

#--------- 6. Create table of go terms associated to each LR-interaction --------#

LR_interactions_go_union <- rbindlist(
  apply(
    LR_interactions,
    MARGIN = 1,
    function(row) {
      LIGAND_GO <- unique(c(
        LR_genes_go[[row[["LIGAND_1"]]]],
        LR_genes_go[[row[["LIGAND_2"]]]]
      ))
      RECEPTOR_GO <- unique(c(
        LR_genes_go[[row[["RECEPTOR_1"]]]],
        LR_genes_go[[row[["RECEPTOR_2"]]]],
        LR_genes_go[[row[["RECEPTOR_3"]]]]
      ))
      res_union <- unique(c(LIGAND_GO, RECEPTOR_GO))
      res_union <- data.table(
        LR_SORTED = rep(row[["LR_SORTED"]], length(res_union)),
        GO_union = res_union
      )
    }
  )
)

LR_interactions_go_intersection <- rbindlist(
  apply(
    LR_interactions,
    MARGIN = 1,
    function(row) {
      LIGAND_GO <- unique(c(
        LR_genes_go[[row[["LIGAND_1"]]]],
        LR_genes_go[[row[["LIGAND_2"]]]]
      ))
      RECEPTOR_GO <- unique(c(
        LR_genes_go[[row[["RECEPTOR_1"]]]],
        LR_genes_go[[row[["RECEPTOR_2"]]]],
        LR_genes_go[[row[["RECEPTOR_3"]]]]
      ))
      res_inter <- intersect(LIGAND_GO, RECEPTOR_GO)
      if(length(res_inter) > 0) {
        res_inter <- data.table(
          LR_SORTED = rep(row[["LR_SORTED"]], length(res_inter)),
          GO_intersection = res_inter
        )
      } else {
        res_inter <- NULL
      }
      return(res_inter)
    }
  )
)






