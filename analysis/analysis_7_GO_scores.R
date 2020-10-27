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

#CCI_all <- readRDS("../data_scAgeCom/analysis/a4_data_diffcom_all_filtered.rds") 
#note change the path to where the rds files is stored on your computer

CCI_all <- readRDS("../data_scAgeCom/scdiffcom_all_results/scdiffcom_tms_facs_mixed_size_factor_logTRUE_10000iter/scdiffcom_Lung.rds")
#CCI_all <- CCI_all$scdiffcom_dt_filtered
CCI_all <- list(CCI_all)

## 3. Load the GO terms and family associated to each GENE #####

#use LR_genes_go_family for the script analysis_2_c (for me it is saved on the disk)
LR_genes_go_family <- readRDS("../data_scAgeCom/analysis/analysis_2_LR_GO-family_terms.rds")

#use LR_genes_info from th other script
LR_genes_info <- readRDS("../data_scAgeCom/analysis/analysis_2_LR_genes_info.rds")

## 4. Add GO term columsn to the table of results 
existing_colnames <- c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3))
new_colnames <- paste0(existing_colnames, "_GO_terms")

CCI_all <- lapply(CCI_all, function(dataset) {
  dataset$scdiffcom_dt_filtered[, (new_colnames) := lapply(existing_colnames, function(col) {
    LR_genes_go_family[.SD, on = paste0("gene==", col), x.GO_terms]
  })]
  return(dataset)
})

## 5. Add new columns with the intersection and union of the GO terms (takes one minute or so)
CCI_all <- lapply(CCI_all, function(dataset) {
  dataset$scdiffcom_dt_filtered[, GO_intersection := sapply(1:nrow(.SD), function(i) {
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
  return(dataset)
})

CCI_all <- lapply(CCI_all, function(dataset) {
  dataset$scdiffcom_dt_filtered[, GO_union := sapply(1:nrow(.SD), function(i) {
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
  return(dataset)
})


## 6. Score for each GO term

#get all GO terms one by one
GO_terms <- unique(LR_genes_info$name_1006)

GO_terms_keep <- GO_terms[GO_terms %in% unique(unlist(strsplit(CCI_all[[1]]$scdiffcom_dt_filtered$GO_intersection, ",")))]

go_ora <- lapply(
  GO_terms_keep,
  function(go_t) {
    CCI_all <- lapply(CCI_all, function(dataset) {
      dataset$scdiffcom_dt_filtered[, temp_col := ifelse(grepl(go_t, GO_intersection, fixed = TRUE), TRUE, FALSE)]
      return(dataset)
    })
    test <- scDiffCom::run_ORA(
      CCI_all[[1]],
      verbose = TRUE,
      categories = "temp_col"
    )
    test <- test$ORA
    test <- test[Value == TRUE]
    test[, name := go_t]
  }
)

go_scores <- rbindlist(lapply(
  GO_terms_keep,
  function(go_t) {
    temp <- CCI_all[[1]]$scdiffcom_dt_filtered[grepl(go_t, GO_intersection)]
    res <- data.table(
      GO = go_t,
      sum_up = sum(temp[REGULATION_SIMPLE == "UP"][["LOGFC"]]),
      sum_down = sum(temp[REGULATION_SIMPLE == "DOWN"][["LOGFC"]]),
      sum_flat = sum(temp[REGULATION_SIMPLE == "FLAT"][["LOGFC"]])
      )
    return(res)
  }
))

go_scores[, diff := sum_up + sum_down]

test_ora <- rbindlist(go_ora)

test_ora_up <- test_ora[pval_adjusted_UP <= 0.05 & OR_UP >= 1, c("pval_adjusted_UP", "OR_UP", "name")][order(pval_adjusted_UP)]
test_ora_down <- test_ora[pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1, c("pval_adjusted_DOWN", "OR_DOWN", "name")][order(pval_adjusted_DOWN)]


ggplot(test_ora, aes(x = log2(OR_UP), y = -log10(pval_adjusted_UP))) + geom_point()

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

######
library(data.table)
library(scDiffCom)
LR <- LR6db$LR6db_curated

LR_genes <- unique(unlist(LR[,c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
LR_genes <- LR_genes[!is.na(LR_genes)]

genage <- read.csv("../../../../../genage_models_export.tsv", header = TRUE, sep = "\t")

LR_long_genes <- genage$Gene.Symbol[genage$Gene.Symbol %in% LR_genes]


LR_genage <- LR[LIGAND_1 %in% LR_long_genes | LIGAND_2 %in% LR_long_genes |
                  RECEPTOR_1 %in% LR_long_genes | RECEPTOR_2 %in% LR_long_genes | RECEPTOR_3 %in% LR_long_genes]

cellage <- readRDS("../../P2_scAgingNetworks/scSenescence/senescence_2.rds")
setDT(cellage)
cellage <- cellage[database == "cellage"]
cellage_ortho <- get_orthologs(cellage$gene, input_species = "human")

LR_cellage_genes <- cellage_ortho$mouse_symbol[cellage_ortho$mouse_symbol %in% LR_genes]

#####

res_facs <- DATASETS$`TMS FACS Male`$results


test <- res_facs[REGULATION %in% c("UP"), .N, by = c("LR_NAME", "TISSUE", "REGULATION")]

test2 <- test[, .N, by = c("LR_NAME")]

test3 <- dcast.data.table(
  test,
  formula = LR_NAME ~ TISSUE,
  value.var = "N"
)


testd <- res_facs[REGULATION %in% c("DOWN"), .N, by = c("LR_NAME", "TISSUE", "REGULATION")]

testd2 <- testd[, .N, by = c("LR_NAME")]

testd3 <- dcast.data.table(
  testd,
  formula = LR_NAME ~ TISSUE,
  value.var = "N"
)

library(parallel)
library(foreach)

data <- 1:1e9
data_list <- list("1" = data,
                  "2" = data,
                  "3" = data,
                  "4" = data)
detectCores()
cl <- parallel::makeCluster(detectCores())
# Activate cluster for foreach library
doParallel::registerDoParallel(cl)
time_foreach <- system.time({
  r <- foreach::foreach(i = 1:length(data_list),
                        .combine = rbind) %dopar% {
                          mean(data_list[[i]])
                        }
})
time_foreach[3]
# Stop cluster to free up resources
parallel::stopCluster(cl)

