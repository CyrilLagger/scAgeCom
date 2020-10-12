####################################################
##
## Project: scAgeCom
## 
## September 2020
##
## anais.equey@etu.univ-amu.fr 
## cyril.lagger@liverpool.ac.uk 
##
## Script to annotate the (curated) ligand-receptor
## interactions.
##
####################################################
##

## 1. Load libraries ####

library(scDiffCom)
library(data.table)
library(biomaRt)

# Analysis data path
dir_data_analysis <- "../data_scAgeCom/analysis/"

## 2. Load the curated ligand-receptor interactions ####

#the database of LR interactions is included in scDiffCom and can be loaded as follow
LR_interactions <- scDiffCom::LR6db$LR6db_curated
#this gives a data.table of 4507 interactions (rows) with 22 columns
head(LR_interactions, 10) #display the 10 first rows in the console

## 3. Get all unique genes in the data.table ####

#Ligand genes
L_genes <- unique(unlist(LR_interactions[, c("LIGAND_1", "LIGAND_2")]))
L_genes <- L_genes[!is.na(L_genes)] # needed to remove NA (empty) values
#Receptor genes
R_genes <- unique(unlist(LR_interactions[, c("RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
R_genes <- R_genes[!is.na(R_genes)] # needed to remove NA (empty) values
#All genes (note that some genes may appear as both ligand and receptors)
LR_genes <- unique(c(L_genes, R_genes))

## 4. Get all previous annotation and family categories ####

#we just need to split the characters that have been pasted by a comma (e.g. "A,B" becomes "A" and "B")
#and remove some blanck spaces.
LR_annotation <- sort(unique(gsub("^\\s+|\\s+$", "", unlist(strsplit(unique(LR_interactions$ANNOTATION), split =  "\\,|\\/")))))
LR_family <- sort(unique(gsub("^\\s+|\\s+$", "", unlist(strsplit(unique(LR_interactions$FAMILY), split =  "\\,|\\/")))))
LR_subfamily <- sort(unique(gsub("^\\s+|\\s+$", "", unlist(strsplit(unique(LR_interactions$SUBFAMILY), split =  "\\,|\\/")))))
LR_terms <- sort(unique(c(LR_annotation, LR_family, LR_subfamily)))

#we can display them in the console
LR_terms

## 5. Get the Ensembl terms (GO, family or else) ####

#create an access to ensembl database via biomart (might take a minute or so)
mart <- biomaRt::useMart(
  "ensembl",
  dataset = "mmusculus_gene_ensembl"
)
#retrieve the information we want for our genes of interest
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

saveRDS(LR_genes_info, paste0(dir_data_analysis, "analysis_2_LR_genes_info.rds"))

#display the first 10 rows
head(LR_genes_info, 10)

## 6. Create a table of the GO/family term frequency ####

LR_go_frequency <- LR_genes_info[, .N, by = c("name_1006", "namespace_1003", "go_id")]
LR_go_frequency <- LR_go_frequency[order(-N)]
LR_go_frequency_top <- LR_go_frequency[N >= 30]
LR_go_frequency_top <- LR_go_frequency_top[-c(1:11)]

LR_family_frequency <- LR_genes_info[, .N, by = c("family_description")]
LR_family_frequency <- LR_family_frequency[order(-N)]

LR_genes_info_filtered <- LR_genes_info[name_1006 %in% LR_go_frequency_top$name_1006]

## 7. Combine the GO and family terms (one row per gene) ####

LR_genes_go <- as.data.table(
  sapply(unique(LR_genes_info$mgi_symbol), function(i) {
    temp <- LR_genes_info[mgi_symbol == i]
    res <- paste0(sort(unique(temp$name_1006)), collapse = ",")
    return(res)
  },
  USE.NAMES = TRUE),
  keep.rownames = TRUE
)
colnames(LR_genes_go) <- c("gene", "GO_terms")

LR_genes_go_filtered <- as.data.table(
  sapply(unique(LR_genes_info_filtered$mgi_symbol), function(i) {
    temp <- LR_genes_info_filtered[mgi_symbol == i]
    res <- paste0(sort(unique(temp$name_1006)), collapse = ",")
    return(res)
  },
  USE.NAMES = TRUE),
  keep.rownames = TRUE
)
colnames(LR_genes_go_filtered) <- c("gene", "GO_terms")

LR_genes_family <- as.data.table(
  sapply(unique(LR_genes_info$mgi_symbol), function(i) {
    temp <- LR_genes_info[mgi_symbol == i]
    res <- paste0(sort(unique(temp$family_description)), collapse = ",")
    return(res)
  },
  USE.NAMES = TRUE),
  keep.rownames = TRUE
)
colnames(LR_genes_family) <- c("gene", "family_terms")

LR_genes_go_family <- merge.data.table(
  LR_genes_go,
  LR_genes_family,
  by = "gene",
  all = TRUE
)
saveRDS(LR_genes_go_family, paste0(dir_data_analysis, "analysis_2_LR_GO-family_terms.rds"))



## 8. #####
existing_colnames <- c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3))
new_colnames_2 <- paste0(existing_colnames, "_GO_terms_2")

LR_interactions[, (new_colnames) := lapply(existing_colnames, function(col) {
    LR_genes_go[.SD, on = paste0("gene==", col), x.GO_terms]
})]

LR_interactions[, (new_colnames_2) := lapply(existing_colnames, function(col) {
  LR_genes_go[.SD, on = paste0("gene==", col), x.GO_terms]
})]

LR_interactions[, GO_intersection := sapply(1:nrow(.SD), function(i) {
  temp_L <- sapply(paste0("LIGAND_", 1:2, "_GO_terms"), function(col) {
    unlist(strsplit(get(col)[[i]], ",", fixed = TRUE))
  })
  temp_L <- temp_L[!is.na(temp_L)]
  temp_L <- unlist(temp_L)
  temp_R <- sapply(paste0("RECEPTOR_", 1:3, "_GO_terms"), function(col) {
    unlist(strsplit(get(col)[[i]], ",", fixed = TRUE))
  })
  temp_R <- temp_R[!is.na(temp_R)]
  temp_R <- unlist(temp_R)
  temp_inter <- intersect(temp_L, temp_R)
  temp_inter <- paste0(temp_inter, collapse = ",")
}) ]

LR_interactions[, GO_union := sapply(1:nrow(.SD), function(i) {
  temp_L <- sapply(paste0("LIGAND_", 1:2, "_GO_terms"), function(col) {
    unlist(strsplit(get(col)[[i]], ",", fixed = TRUE))
  })
  temp_L <- temp_L[!is.na(temp_L)]
  temp_L <- unlist(temp_L)
  temp_R <- sapply(paste0("RECEPTOR_", 1:3, "_GO_terms"), function(col) {
    unlist(strsplit(get(col)[[i]], ",", fixed = TRUE))
  })
  temp_R <- temp_R[!is.na(temp_R)]
  temp_R <- unlist(temp_R)
  temp_union <- union(temp_L, temp_R)
  temp_union <- paste0(temp_union, collapse = ",")
  
}) ]

LR_interactions[GO_intersection == ""]$LR_SORTED
LR_interactions[GO_union == ""]

test <- LR_interactions[, c("LR_SORTED", "LIGAND_1_GO_terms", "LIGAND_1_GO_terms_2", "GO_intersection")]

unique(LR_genes[!(LR_genes %in% LR_genes_go$gene)])

######################

## 5. Add new columns with the intersection and union of the GO terms (takes one minute or so)



##### Below: code for testing ####

## 8. On how to extract group of genes with specific terms ####

#Example: let's retrieve all the go terms with "immune" and all the genes with "immune"
go_immune <- LR_go_frequency[grepl("immune", name_1006)]
LR_genes_immune_1 <- LR_genes_go[grepl("immune", GO_terms)]
LR_genes_immune_2 <- LR_genes_go[grepl("immune response", GO_terms)]

#inflammation
go_inflam <- LR_go_frequency[grepl("inflam", name_1006)]
LR_genes_inflam_1 <- LR_genes_go[grepl("inflam", GO_terms)]

#collagen
go_colla <- LR_go_frequency[grepl("collagen", name_1006)]
LR_genes_colla_1 <- LR_genes_go[grepl("collagen", GO_terms)]


## Some test
LR_go_frequency_short <- LR_go_frequency[N >= 10]
LR_genes_go_short <- LR_genes_go[name_1006 %in% LR_go_frequency_short$name_1006]

LR_genes_go_combined_short <- as.data.table(
  sapply(unique(LR_genes_go_short$mgi_symbol), function(i) {
    temp <- LR_genes_go_short[mgi_symbol == i]
    res <- paste0(sort(temp$name_1006), collapse = ",")
    return(res)
  },
  USE.NAMES = TRUE),
  keep.rownames = TRUE
)
colnames(LR_genes_go_combined_short) <- c("gene", "GO_terms")

## Some test with omnipathR
library(OmnipathR)
get_annotation_resources()
test <- import_omnipath_annotations(proteins = LR_genes)

LR_genes[!(LR_genes %in% human_LR$mouse_symbol)]

human_LR <- get_orthologs(LR_genes, input_species = "mouse")

annotations <- import_omnipath_annotations(
  proteins = human_LR$human_symbol[1:500],
  resources = c("NetPath")
  #resources = c("SignaLink_pathway")
  #resources = c("GO_Intercell", "SignaLink_pathway")
)

annotations2 <- import_omnipath_annotations(
  proteins = human_LR$human_symbol[501:1000],
  resources = c("NetPath", "SignaLink_pathway")
  #resources = c("SignaLink_pathway")
  #resources = c("GO_Intercell", "SignaLink_pathway")
)

annotations3 <- import_omnipath_annotations(
  proteins = human_LR$human_symbol[1001:1500],
  resources = c("NetPath", "SignaLink_pathway")
  #resources = c("SignaLink_pathway")
  #resources = c("GO_Intercell", "SignaLink_pathway")
)


