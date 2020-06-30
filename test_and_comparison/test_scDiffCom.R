####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - June 2020
##
## Test that the package scDiffCom returns
## correct result when applied to aging research.
##
####################################################
##

library(Seurat)
library(scDiffCom)
#library(profvis)
#library(microbenchmark)

#load a Seurat objects
#seurat_tms_test_file.rds corresponds to the Liver tissue from TMS FACS data.

seurat_test <- readRDS("../data_seurat_example.rds")
seurat_test <- NormalizeData(seurat_test, assay = "RNA")
seurat_test$age_group <- ifelse(seurat_test$age %in% c('1m', '3m'), 'young', 'old' )
seurat_test$cell_ontology_class <- as.character(seurat_test$cell_ontology_class)

#load LR database from scDiffCom
LR_test <- LRall$LRall_one2one
LR_test <- LR_test[LR_test$sCsR, c(2,3,1)]

#run scDiffCom (in non-log space with 1000 permutations)
start_time <- Sys.time()
diffcom_test <- run_diffcom(
  seurat_obj = seurat_test,
  LR_data = LR_test,
  seurat_cell_type_id = "cell_ontology_class",
  condition_id = "age_group",
  assay = "RNA",
  slot = "data",
  log_scale = TRUE,
  min_cells = 5,
  threshold = 0.1,
  permutation_analysis = TRUE,
  iterations = 10000,
  return_distr = FALSE,
  one_sided = FALSE
)
end_time <- Sys.time()
end_time - start_time

#save the file for future reference
#saveRDS(object = diffcom_test, file = "test_and_comparison/data_results_diffcom.rds")
#saveRDS(object = diffcom_test, file = "test_and_comparison/data_results_diffcom_10000iter_log.rds")
diffcom_test <- readRDS("test_and_comparison/data_results_diffcom.rds")

#Pick 5 random CCI for comparison below
cci_test_list <- list(diffcom_test[LR_DETECTED_young & LR_DETECTED_old][28],
                      diffcom_test[!LR_DETECTED_young & LR_DETECTED_old][120],
                      diffcom_test[LR_DETECTED_young & !LR_DETECTED_old][345],
                      diffcom_test[!LR_DETECTED_young & !LR_DETECTED_old][90],
                      diffcom_test[LR_DETECTED_young & LR_DETECTED_old][60])

#Manual check that the LR scores are correct
compare_cci_score_test <- function(cci_test, condition) {
  seurat_sub_l <- subset(seurat_test, features = cci_test[["L_GENE"]],
                         subset = (cell_ontology_class == cci_test[["L_CELLTYPE"]] &
                                     age_group == condition))
  seurat_sub_r <- subset(seurat_test, features = cci_test[["R_GENE"]],
                         subset = (cell_ontology_class == cci_test[["R_CELLTYPE"]] &
                                     age_group == condition))
  LR_avg_expr <- (mean(expm1(seurat_sub_l[["RNA"]]@data[1,])) + mean(expm1(seurat_sub_r[["RNA"]]@data[1,])))/2
  return(LR_avg_expr - cci_test[[paste0("LR_SCORE_", condition)]])
}

lapply(cci_test_list, compare_cci_score_test, "young")
lapply(cci_test_list, compare_cci_score_test, "old")

#Manual check that the LR differential p-values are correct
compare_cci_pvalue_diff_test <- function(seurat_obj, cci_test, iterations) {
  seurat_temp <- subset(seurat_obj, features = c(cci_test[["L_GENE"]], cci_test[["R_GENE"]]),
                        subset = cell_ontology_class %in% c(cci_test[["L_CELLTYPE"]], cci_test[["R_CELLTYPE"]]))
  get_shuffled_LR <- function(permut) {
    if(permut) {
      seurat_temp$age_group[seurat_temp$cell_ontology_class == cci_test[["L_CELLTYPE"]] ] <- 
        sample(seurat_temp$age_group[seurat_temp$cell_ontology_class == cci_test[["L_CELLTYPE"]] ])
      seurat_temp$age_group[seurat_temp$cell_ontology_class == cci_test[["R_CELLTYPE"]] ] <- 
        sample(seurat_temp$age_group[seurat_temp$cell_ontology_class == cci_test[["R_CELLTYPE"]] ])
    }
    seurat_sub_l_young <- subset(seurat_temp, features = cci_test[["L_GENE"]],
                           subset = (cell_ontology_class == cci_test[["L_CELLTYPE"]] &
                                       age_group == "young"))
    seurat_sub_r_young <- subset(seurat_temp, features = cci_test[["R_GENE"]],
                           subset = (cell_ontology_class == cci_test[["R_CELLTYPE"]] &
                                       age_group == "young"))
    seurat_sub_l_old <- subset(seurat_temp, features = cci_test[["L_GENE"]],
                                 subset = (cell_ontology_class == cci_test[["L_CELLTYPE"]] &
                                             age_group == "old"))
    seurat_sub_r_old <- subset(seurat_temp, features = cci_test[["R_GENE"]],
                                 subset = (cell_ontology_class == cci_test[["R_CELLTYPE"]] &
                                             age_group == "old"))
    LR_diff <- (mean(expm1(seurat_sub_l_old[["RNA"]]@data[1,])) + mean(expm1(seurat_sub_r_old[["RNA"]]@data[1,])))/2 - 
      (mean(expm1(seurat_sub_l_young[["RNA"]]@data[1,])) + mean(expm1(seurat_sub_r_young[["RNA"]]@data[1,])))/2
    return(LR_diff)
  }
  all_LR <- replicate(iterations, get_shuffled_LR(TRUE))
  all_LR <- c(all_LR, get_shuffled_LR(FALSE))
  pval <- sum(abs(all_LR[1:iterations]) >= abs(all_LR[(iterations + 1)])) / iterations
  return(list(pval = pval,
              distr = all_LR))
}

start_time <- Sys.time()
test_pval_diff <- compare_cci_pvalue_diff_test(seurat_test, cci_test_list[[4]], iterations = 1000)
end_time <- Sys.time()
end_time - start_time

cci_test_list[[4]]$pvals_diff
#test$distr
test_pval_diff$pval
hist(test_pval_diff$distr)
mean(test_pval_diff$distr)
mean(test_pval_diff$distr)/sd(test_pval_diff$distr)

#Manual check that the LR differential p-values are correct
compare_cci_pvalue_spec_test <- function(seurat_obj, cci_test, iterations, condition) {
  seurat_temp <- subset(seurat_obj, features = c(cci_test[["L_GENE"]], cci_test[["R_GENE"]]),
                        subset =  age_group == condition)
  get_shuffled_LR_ct <- function(permut) {
    if(permut) {
      seurat_temp$cell_ontology_class <- sample(as.character(seurat_temp$cell_ontology_class))
    }
    
    seurat_sub_l <- subset(seurat_temp, features = cci_test[["L_GENE"]],
                                 subset = (cell_ontology_class == cci_test[["L_CELLTYPE"]]))
    seurat_sub_r <- subset(seurat_temp, features = cci_test[["R_GENE"]],
                                 subset = (cell_ontology_class == cci_test[["R_CELLTYPE"]]))
    (mean(expm1(seurat_sub_l[["RNA"]]@data[1,])) + mean(expm1(seurat_sub_r[["RNA"]]@data[1,])))/2
    return((mean(expm1(seurat_sub_l[["RNA"]]@data[1,])) + mean(expm1(seurat_sub_r[["RNA"]]@data[1,])))/2)
  }
  all_LR <- replicate(iterations, get_shuffled_LR_ct(TRUE))
  all_LR <- c(all_LR, get_shuffled_LR_ct(FALSE))
  pval <- sum(all_LR[1:iterations] >= all_LR[(iterations + 1)]) / iterations
  return(list(pval = pval,
              distr = all_LR))
}

start_time <- Sys.time()
test_pval_spec <- compare_cci_pvalue_spec_test(seurat_test, cci_test_list[[5]], iterations = 1000, "old")
end_time <- Sys.time()
end_time - start_time

cci_test_list[[5]]$PVAL_old
#test_pval_spec$distr
test_pval_spec$pval
hist(test_pval_spec$distr)

#volcano plot
library(ggrepel)
diffcom_1000iter <- readRDS("test_and_comparison/data_results_diffcom_10000iter_log.rds")
diffcom_1000iter[, LR_DIFF := LR_SCORE_old - LR_SCORE_young]
ggplot(diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old), ], 
       aes(x = LR_DIFF, y = -log10(BH_PVAL_DIFF + 1E-4))) +
  geom_point() +
  geom_vline(xintercept = log(1.1)) +
  geom_vline(xintercept = -log(1.1)) +
  geom_hline(yintercept = -log10(0.05))

+
  geom_text(
    data = diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old) & BH_PVAL_DIFF < 0.01 & abs(LR_DIFF) > 0.5, ],
    aes(label = LR_GENES),
    size = 5
    #box.padding = unit(0.35, "lines"),
    #point.padding = unit(0.3, "lines")
  )

hist(log10(diffcom_1000iter$LR_SCORE_old), breaks = 100)
hist(log10(diffcom_1000iter$LR_SCORE_young), breaks = 100)

ggplot(diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old) & (BH_PVAL_young <= 0.05 | BH_PVAL_old <= 0.05) & 
                          (LR_SCORE_young > 0.1 | LR_SCORE_old > 0.1), ], 
       aes(x = LR_DIFF, y = -log10(BH_PVAL_DIFF + 1E-4))) +
  geom_point() +
  geom_vline(xintercept = log(1.1)) +
  geom_vline(xintercept = -log(1.1)) +
  geom_hline(yintercept = -log10(0.05))


+
  geom_text(
    data = diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old) & BH_PVAL_DIFF < 0.01 & abs(LR_DIFF) > 0.3, ],
    aes(label = L_CELLTYPE),
    size = 5
    #box.padding = unit(0.35, "lines"),
    #point.padding = unit(0.3, "lines")
  )


