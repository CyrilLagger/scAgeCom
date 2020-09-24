library(data.table)
#library(igraph)
#library(purrr)
library(glue)
library(scDiffCom)
library(arules)
#library(pheatmap)
#library(RColorBrewer)
#library(ggplot2)
#library(gridExtra)
#library(grid)
#library(gtable)

source("src/src_1_filtering.R")

DIFFCOM_RESULT_PATH <- list(
  calico = "../data_scAgeCom/scDiffCom_Complex_results/diffcom_calico_size_factor_log_10000iter_curated/",
  #calico_sub = "../../data_scAgeCom/scDiffCom_results/diffcom_calico_subtype_size_factor_log_10000iter_mixed",
  tms_facs = "../data_scAgeCom/scDiffCom_Complex_results/diffcom_tms_facs_size_factor_log_10000iter_curated/",
  #tms_facs_female = "../../data_scAgeCom/scDiffCom_results_ageBySex/diffcom_tms_facs_size_factor_log_10000iter_mixed_female/",
  #tms_facs_male = "../../data_scAgeCom/scDiffCom_results_ageBySex/diffcom_tms_facs_size_factor_log_10000iter_mixed_male/",
  tms_droplet = "../data_scAgeCom/scDiffCom_Complex_results/diffcom_tms_droplet_size_factor_log_10000iter_curated/"#,
  #tms_droplet_female = "../../data_scAgeCom/scDiffCom_results_ageBySex/diffcom_tms_droplet_size_factor_log_10000iter_mixed_female/",
  #tms_droplet_male = "../../data_scAgeCom/scDiffCom_results_ageBySex/diffcom_tms_droplet_size_factor_log_10000iter_mixed_male/"
)

calico_tiss <- c("Kidney", "Lung", "Spleen")
tms_facs_base_tiss <- c("Aorta", "BAT", "Bladder", "Brain_Myeloid",
                        "Brain_Non-Myeloid", "GAT",
                        "Heart", "Kidney", "Large_Intestine",
                        "Limb_Muscle", "Liver", "Lung", 
                        "Marrow", "MAT", "Pancreas", "SCAT", "Skin",
                        "Spleen", "Thymus", "Tongue", "Trachea")
tms_droplet_base_tiss <- c("Heart_and_Aorta", "Kidney",
                           "Limb_Muscle", "Liver", "Lung",
                           "Marrow", "Spleen")
TISSUE_DATASET <- list(
  calico = calico_tiss,
  #calico_sub = calico_tiss,
  tms_facs = sort(c(tms_facs_base_tiss, "Diaphragm", "Mammary_Gland")),
  #tms_facs_female = sort(c(tms_facs_base_tiss, "Mammary_Gland")),
  #tms_facs_male = sort(c(tms_facs_base_tiss, "Diaphragm")),
  tms_droplet = sort(c(tms_droplet_base_tiss, "Bladder", "Mammary_Gland", "Thymus", "Tongue"))#,
  #tms_droplet_female = sort(c(tms_droplet_base_tiss, "Mammary_Gland", "Thymus")),
  #tms_droplet_male = sort(c(tms_droplet_base_tiss, "Bladder", "Tongue"))
)

message("Combining all tissues together")
DATASETS <- mapply(
  bind_tissues,
  DIFFCOM_RESULT_PATH,
  TISSUE_DATASET,
  MoreArgs = list(pre_filtering = TRUE),
  SIMPLIFY = FALSE
)

hist(c(DATASETS[[1]]$LR_SCORE_old, DATASETS[[1]]$LR_SCORE_young), breaks = 100)
quantile(c(DATASETS[[1]]$LR_SCORE_old, DATASETS[[1]]$LR_SCORE_young), 0.25)

DATASETS_wTYPES_new <- lapply(
  DATASETS,
  function(data) {
    data <- analyze_CCI_per_tissue(
      data,
      cutoff_quantile = 0.25
    )
  }
)

DATASETS_wTYPES_old <- lapply(
  DATASETS,
  function(data) {
    cutoff <- quantile(c(data$LR_SCORE_old, data$LR_SCORE_young), 0.25)
    data <- analyze_CCI(
      data,
      cutoff_score = cutoff
    )
  }
)

DATASETS_FILTERED_old <- lapply(
  DATASETS_wTYPES_old,
  function(data) {
    data[!(CASE_TYPE %in% c("FFF"))]
  }
)

DATASETS_FILTERED_new<- lapply(
  DATASETS_wTYPES_new,
  function(data) {
    data[!(CASE_TYPE %in% c("FFF"))]
  }
)

saveRDS(DATASETS_FILTERED, "../data_scAgeCom/analysis/analysis_4_data_diffcom_filter.rds")


test1 <- DATASETS_FILTERED_new$tms_droplet
test2 <- DATASETS_FILTERED_old$tms_droplet
test1[, id := paste(LR_CELLTYPE, LR_NAME, sep = "_")]
test2[, id := paste(LR_CELLTYPE, LR_NAME, sep = "_")]

length(intersect(test1$id, test2$id))

######

source("src/src_2_ora_new.R")

ORA_RESULTS <- lapply(
  DATASETS_FILTERED,
  analyze_ORA
)

saveRDS(ORA_RESULTS, "../data_scAgeCom/analysis/analysis_4_data_ora.rds")

test <- ORA_RESULTS$tms_facs
test2 <- test[Category == "LR_NAME", c("Tissue", "Category", "Value", "pval_UP", "OR_UP")]
test3<- test[Category == "LR_NAME", c("Tissue", "Category", "Value", "pval_DOWN", "OR_DOWN")]

LRdbcur <- LR6db$LR6db_curated
LRdbcur[LIGAND_1 == "Crlf2"]

#######

source("../src/src_3_fpm.R")

fpm_results <- lapply(
  DATASETS_FILTERED,
  analyze_FreqItemSets,
  target = "closed frequent itemsets",
  support = 0.00001,
  confidence = 0.01
)
