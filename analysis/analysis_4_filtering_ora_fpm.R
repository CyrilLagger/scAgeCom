####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - September 2020
##
## Filtering, ORA and FPM on scDiffCom results
##
####################################################
##

## Libraries ####

library(data.table)
library(glue)
library(scDiffCom)
library(arules)

source("src/src_1_filtering.R")

#library(igraph)
#library(purrr)
#library(pheatmap)
#library(RColorBrewer)
#library(ggplot2)
#library(gridExtra)
#library(grid)
#library(gtable)


## Specify the directory with scDiffCom results ####
dir_results <- "../data_scAgeCom/scDiffCom_all_results"

dir_data_analysis <- "../data_scAgeCom/analysis/"

## Get all results directories and all tissues in each directories ####

RESULT_PATHS <- list.dirs(dir_results, recursive = FALSE)

DATASETS <- lapply(
  RESULT_PATHS, 
  function(i) {
    if(grepl("sex", i)) {
      conds <- c("female", "male")
    } else {
      conds <- c("YOUNG", "OLD")
    }
    tissues <- gsub(".*scdiffcom_(.+)\\.rds.*", "\\1", list.files(i))
    results <- bind_tissues(
      path = i,
      list_of_tissues =  tissues, 
      conds = conds,
      pre_filtering = TRUE,
      min_cells = 5
    )
    results <- analyze_CCI_per_tissue(
      results,
      conds = conds,
      cutoff_quantile = 0.25,
      recompute_BH = TRUE
    )
    results <- results[!(CASE_TYPE %in% c("FFF"))]
    list(
      id = gsub(".*scdiffcom_", "", i),
      tissues = tissues,
      conds = conds,
      results = results
    )
  }
)

## Rename datasets ####
lapply(DATASETS, function(i) i$id)

names(DATASETS) <- c(
  "Calico",
  "Calico subtype",
  "TMS Droplet Female",
  "TMS Droplet Male",
  "TMS Droplet Mixed",
  "TMS Droplet bySex",
  "TMS FACS Female",
  "TMS FACS Male",
  "TMS FACS Mixed",
  "TMS FACS bySex"
)

## Save filtered results ####

saveRDS(DATASETS, "../data_scAgeCom/analysis/a4_data_diffcom_all_filtered.rds")

## ORA analysis ####

source("src/src_2_ora.R")

ORA_RESULTS <- lapply(
  DATASETS,
  function(i) {
    analyze_ORA(i$results)
  }
)

saveRDS(ORA_RESULTS, "../data_scAgeCom/analysis/a4_data_ora.rds")

## FPM analysis ####

#source("../src/src_3_fpm.R")

#fpm_results <- lapply(
#  DATASETS_FILTERED,
#  analyze_FreqItemSets,
#  target = "closed frequent itemsets",
#  support = 0.00001,
#  confidence = 0.01
#)


## Add regulation column ####

lapply(
  DATASETS, 
  function(i) {
    i$results[, REGULATION := ifelse(
      CASE_TYPE %in% c("TTFD", "TTFU"),
      "FLAT",
      ifelse(
        CASE_TYPE %in% c("TTTD", "TFTD"),
        "DOWN",
        "UP"
      )
    )]
  }
)

## Compare CASE_TYPE to ligand logFC and receptor logFC ####

plots_LR_log2fc <- lapply(
  c(1,2,3,4,5,7,8,9),
  function(i) {
    temp <- DATASETS[[i]]$results
    temp[, L_MIN_EXPRESSION_OLD := pmin(L1_EXPRESSION_OLD, L2_EXPRESSION_OLD, na.rm = TRUE)]
    temp[, L_MIN_EXPRESSION_YOUNG := pmin(L1_EXPRESSION_YOUNG, L2_EXPRESSION_YOUNG, na.rm = TRUE)]
    temp[, R_MIN_EXPRESSION_OLD := pmin(R1_EXPRESSION_OLD, R2_EXPRESSION_OLD, R3_EXPRESSION_OLD, na.rm = TRUE)]
    temp[, R_MIN_EXPRESSION_YOUNG := pmin(R1_EXPRESSION_YOUNG, R2_EXPRESSION_YOUNG, R3_EXPRESSION_YOUNG, na.rm = TRUE)]
    temp[, L_DIFF := L_MIN_EXPRESSION_OLD - L_MIN_EXPRESSION_YOUNG]
    temp[, R_DIFF := R_MIN_EXPRESSION_OLD - R_MIN_EXPRESSION_YOUNG]
    ggplot(temp, aes(x = L_DIFF*log2(exp(1)), y = R_DIFF*log2(exp(1)), color = REGULATION)) +
      geom_point(size = 0.5) +
      xlab("LOG2FC - Ligand") +
      ylab("LOG2FC - Receptor") +
      ggtitle(names(DATASETS)[[i]]) + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0)
  }
)

plot_full_LR_log2fc <- cowplot::plot_grid(
  plotlist = plots_LR_log2fc,
  ncol = 3,
  align = "v"
)
ggsave(filename = paste0(dir_data_analysis, "a4_plot_full_LR_log2fc.png"),
       plot = plot_full_LR_log2fc, scale = 2)


## Naive GO analysis
library(clusterProfiler)
library(org.Mm.eg.db)

universe_LR <- unique(unlist(LR6db$LR6db_curated[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
universe_LR <- universe_LR[!is.na(universe_LR)]

ego_MF_strong <- lapply(
  c(1,2,3,4,5,7,8,9),
  function(i) {
    temp <- DATASETS[[i]]$results
    up_genes <- unique(unlist(temp[REGULATION == "UP" & LOGFC > log(2), c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
    up_genes <- up_genes[!is.na(up_genes)]
    down_genes <- unique(unlist(temp[REGULATION == "DOWN" & LOGFC < -log(2), c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
    down_genes <- down_genes[!is.na(down_genes)]
    ego_up <- enrichGO(gene       = up_genes,
                       universe      = universe_LR,
                       keyType       = "SYMBOL",
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = FALSE)
    ego_down <- enrichGO(gene       = down_genes,
                         universe      = universe_LR,
                         keyType       = "SYMBOL",
                         OrgDb         = org.Mm.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
    return(list(ego_up, ego_down))
  }
  
)

plot_go_MF_droplet <- cowplot::plot_grid(
  plotlist = list(
    dotplot(ego_MF_strong[[3]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[3]][[2]], showCategory = 10),
    dotplot(ego_MF_strong[[4]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[4]][[2]], showCategory = 10),
    dotplot(ego_MF_strong[[5]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[5]][[2]], showCategory = 10)
  ),
  ncol = 2,
  align = "v",
  labels = c(paste(names(DATASETS)[3], c("UP", "DOWN")),
             paste(names(DATASETS)[4], c("UP", "DOWN")),
             paste(names(DATASETS)[5], c("UP", "DOWN")))
)
plot_go_MF_droplet
ggsave(filename = paste0(dir_data_analysis, "a4_plot_go_MF_droplet.png"),
       plot = plot_go_MF_droplet, scale = 2)


cowplot::plot_grid(
  plotlist = list(
    dotplot(ego_MF_strong[[6]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[6]][[2]], showCategory = 10),
    dotplot(ego_MF_strong[[7]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[7]][[2]], showCategory = 10),
    dotplot(ego_MF_strong[[8]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[8]][[2]], showCategory = 10)
  ),
  ncol = 2,
  align = "v"
)



#################################################################################



names(DATASETS_FILTERED) <- sapply(DATASETS, function(i) {i$id})

FULL_RESULTS_AGE <- rbindlist(
  lapply(DATASETS[c(1,2,3,4,5,7,8,9)], function(i) {i$results}),
  use.names = TRUE,
  idcol = "DATASET"
)

FULL_RESULTS_AGE_2 <- rbindlist(
  lapply(DATASETS_2[c(1,2,3,4,5,7,8,9)], function(i) {i$results}),
  use.names = TRUE,
  idcol = "DATASET"
)



## General statistics per (dataset-) tissue ####

CCI_counts <- FULL_RESULTS_AGE[, .N, by = c("DATASET", "TISSUE", "CASE_TYPE")]
CCI_counts_full <- FULL_RESULTS_AGE[, .N, by = c("DATASET", "TISSUE")]
CCI_counts_full2 <- FULL_RESULTS_AGE_2[, .N, by = c("DATASET", "TISSUE")]

CCI_counts_dc <- dcast.data.table(
  CCI_counts,
  formula = DATASET + TISSUE ~ REGULATION,
  value.var = "N"
)

test <- dcast.data.table(
  CCI_counts,
  formula = TISSUE  ~ DATASET + REGULATION,
  value.var = "N"
)

FULL_RESULTS_AGE[, CCIT := paste(TISSUE, LR_CELLTYPE, LR_NAME, sep = "_")]

head(sort(table(FULL_RESULTS_AGE$CCI), decreasing = TRUE))
head(sort(table(FULL_RESULTS_AGE[ REGULATION == "UP"]$CCIT), decreasing = TRUE), 50)
head(sort(table(FULL_RESULTS_AGE[ REGULATION == "DOWN"]$CCIT), decreasing = TRUE), 50)


test <- ORA_RESULTS$tms_facs
test2 <- test[Category == "LR_NAME", c("Tissue", "Category", "Value", "pval_UP", "OR_UP")]
test3<- test[Category == "LR_NAME", c("Tissue", "Category", "Value", "pval_DOWN", "OR_DOWN")]

LRdbcur <- LR6db$LR6db_curated
LRdbcur[LIGAND_1 == "Crlf2"]

