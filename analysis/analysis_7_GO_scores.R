####################################################
##
## Project: scAgeCom
## 
## September 2020
##
## anais.equey@etu.univ-amu.fr 
## cyril.lagger@liverpool.ac.uk 
##
## Script 
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


## 3. Load the GO terms and family associated to each GENE #####

#use LR_genes_go_family for the script analysis_2_c (for me it is saved on the disk)
#LR_genes_go_family <- readRDS("../data_scAgeCom/analysis/analysis_2_LR_GO-family_terms.rds")

#use LR_genes_info from th other script
#LR_genes_info <- readRDS("../data_scAgeCom/analysis/analysis_2_LR_genes_info.rds")

## 4. Add GO term columsn to the table of results 
existing_colnames <- c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3))
new_colnames <- paste0(existing_colnames, "_GO_terms")

CCI_all <- lapply(CCI_all, function(dataset) {
  dataset[, (new_colnames) := lapply(existing_colnames, function(col) {
    LR_genes_go_family[.SD, on = paste0("gene==", col), x.GO_terms]
  })]
})

## 5. Add new columns with the intersection and union of the GO terms (takes one minute or so)
CCI_all <- lapply(CCI_all, function(dataset) {
  dataset[, GO_intersection := sapply(1:nrow(.SD), function(i) {
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
})

CCI_all <- lapply(CCI_all, function(dataset) {
  dataset[, GO_union := sapply(1:nrow(.SD), function(i) {
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
})

## 6. Score for each GO term

#get all GO terms one by one
GO_terms <- unique(LR_genes_info$name_1006)

#test on TMS facs that change with age
CCI_tmsFACS <- CCI_all$tms_facs
CCI_tmsFACS_age <- CCI_tmsFACS[CASE_TYPE %in% c("FTTU", "TFTD", "TTTD", "TTTU")]

#the part below is very slow and not working properly, I have to think about it

GO_number <- pbsapply(GO_terms[1:500], function(i) {
  sum(grepl(i, CCI_tmsFACS_age$GO_union, fixed = TRUE))
})

hist(GO_number)

head(sort(GO_number), 100)

GO_logfc <- pbsapply(GO_terms, function(i) {
  mean(CCI_tmsFACS_age[grepl(i, CCI_tmsFACS_age$GO_intersection, fixed = TRUE)]$LOGFC)
})

hist(GO_logfc)

head(sort(GO_logfc, decreasing = TRUE))
head(sort(GO_logfc, decreasing = FALSE))


hist(CCI_tmsFACS_age[grepl("zinc ion binding", CCI_tmsFACS_age$GO_union)]$LOGFC, breaks = 100)

test <- CCI_all$calico
test$GO_intersection[[12345]]
test$GO_union[[1234]]

identical(test$GO_intersection[[1]], test$GO_intersection[[25]])
