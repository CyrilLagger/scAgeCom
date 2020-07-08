####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - June 2020
##
## Test how the results of scDiffCom depends
## on various preprocessing and internal parameters
##
####################################################
##

library(Seurat)
library(scDiffCom)
library(ggplot2)
library(VennDiagram)
library(data.table)

#load a Seurat objects
#seurat_tms_test_file.rds corresponds to the Liver tissue from TMS FACS data.
seurat_test <- readRDS("../data_seurat_example.rds")
seurat_test$age_group <- ifelse(seurat_test$age %in% c('1m', '3m'), 'young', 'old' )
seurat_test$cell_ontology_class <- as.character(seurat_test$cell_ontology_class)

#load LR database from scDiffCom
LR_test <- LRall$LRall_one2one
LR_test <- LR_test[LR_test$sCsR, c(2,3,1)]

#choice of normalization
seurat_test <- NormalizeData(seurat_test, assay = "RNA")
seurat_test <- SCTransform(
  object = seurat_test,
  return.only.var.genes = FALSE
)

#list of parameters to run over
param_list <- list(
  list(assay = "RNA", log_scale = TRUE, one_sided = TRUE),
  list(assay = "RNA", log_scale = TRUE, one_sided = FALSE),
  list(assay = "RNA", log_scale = FALSE, one_sided = TRUE),
  list(assay = "RNA", log_scale = FALSE, one_sided = FALSE),
  list(assay = "SCT", log_scale = TRUE, one_sided = TRUE),
  list(assay = "SCT", log_scale = TRUE, one_sided = FALSE),
  list(assay = "SCT", log_scale = FALSE, one_sided = TRUE),
  list(assay = "SCT", log_scale = FALSE, one_sided = FALSE)
)

#run scDiffCom on all listed parameters
diffcom_bind <- rbindlist(
  lapply(
    param_list,
    function(
      param
    ) {
      run_diffcom(
        seurat_obj = seurat_test,
        LR_data = LR_test,
        seurat_cell_type_id = "cell_ontology_class",
        condition_id = "age_group",
        assay = param[["assay"]],
        log_scale = param[["log_scale"]],
        iterations = 1000,
        one_sided = param[["one_sided"]],
        permutation_analysis = TRUE
      )
    }
  ),
  idcol = "param_group"
)
#save the result for later use
saveRDS(diffcom_bind, file = "test_and_comparison/data_results_diffcom_param.rds")

#be careful that sct removes genes
diffcom_param <- diffcom_bind[LR_GENES %in% diffcom_bind[param_group == 5, LR_GENES]]
nrow(diffcom_param)/8
#diffcom_param <- diffcom_param[, c(1,2,5,6,7,8,9,10,11,12,13,14,15,16)]
diffcom_param[, LR_DIFF := LR_SCORE_old - LR_SCORE_young]
diffcom_param[, LR_LOGFC := log(LR_SCORE_old/LR_SCORE_young)]

#compare detection rate
#it should be identical for the different log and one_sided choices
identical(diffcom_param[param_group == 1, LR_DETECTED_young], diffcom_param[param_group == 4, LR_DETECTED_young])
identical(diffcom_param[param_group == 5, LR_DETECTED_young], diffcom_param[param_group == 8, LR_DETECTED_young])
#there is only a difference depending on the normalization
test_detec <- dcast(
  diffcom_param[, c("param_group", "LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "LR_DETECTED_old")],
  formula = LR_GENES + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_DETECTED_old"
)
table(test_detec$`1`, test_detec$`5`)
#we see that scTransform detects more signals than size-factor
sum(test_detec$`1`)
sum(test_detec$`5`)

#compare LR scores
test_LR_score <- dcast(
  diffcom_param[, c("param_group", "LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "LR_SCORE_old")],
  formula = LR_GENES + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_SCORE_old"
)
ggplot(test_LR_score, aes(x = `1`, y = `3`)) + geom_point()
ggplot(test_LR_score, aes(x = `1`, y = `5`)) + geom_point()


#Compare LR diff, without filtering by p-values (do that later on below)
test_LR_diff <- dcast(
  diffcom_param[!(LR_DETECTED_young == FALSE & LR_DETECTED_old == FALSE),
                c("param_group", "LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "LR_DIFF")],
  formula = LR_GENES + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_DIFF"
)

table(test_LR_diff$`1` > 0, test_LR_diff$`3` > 0)
table(test_LR_diff$`5` > 0, test_LR_diff$`7` > 0)
table(test_LR_diff$`1` > 0, test_LR_diff$`5` > 0)
table(test_LR_diff$`1` > 0, test_LR_diff$`7` > 0)
table(test_LR_diff$`3` > 0, test_LR_diff$`5` > 0)
table(test_LR_diff$`3` > 0, test_LR_diff$`7` > 0)


ggplot(diffcom_param[(LR_DETECTED_young | LR_DETECTED_old) & param_group == 2, ], aes(x = LR_DIFF, y = -log10(BH_PVAL_DIFF))) +
  geom_point() +
  geom_vline(xintercept = log(1.1)) +
  geom_vline(xintercept = -log(1.1)) +
  geom_hline(yintercept = -log10(0.05))

ggplot(diffcom_param[(LR_DETECTED_young | LR_DETECTED_old) & param_group == 1, ], aes(x = LR_DIFF, y = -log10(BH_PVAL_DIFF))) + geom_point()
ggplot(diffcom_param[(LR_DETECTED_young | LR_DETECTED_old) & param_group == 3, ], aes(x = LR_LOGFC, y = -log10(BH_PVAL_DIFF))) + geom_point()
ggplot(diffcom_param[(LR_DETECTED_young | LR_DETECTED_old) & param_group == 4, ], aes(x = LR_LOGFC, y = -log10(BH_PVAL_DIFF))) + geom_point()

ggplot(diffcom_param[(LR_DETECTED_young | LR_DETECTED_old) & param_group == 5, ], aes(x = LR_DIFF, y = -log10(BH_PVAL_DIFF))) + geom_point()
ggplot(diffcom_param[(LR_DETECTED_young | LR_DETECTED_old) & param_group == 6, ], aes(x = LR_LOGFC, y = -log10(BH_PVAL_DIFF))) + geom_point()
ggplot(diffcom_param[(LR_DETECTED_young | LR_DETECTED_old) & param_group == 7, ], aes(x = LR_LOGFC, y = -log10(BH_PVAL_DIFF))) + geom_point()
ggplot(diffcom_param[(LR_DETECTED_young | LR_DETECTED_old) & param_group == 8, ], aes(x = LR_LOGFC, y = -log10(BH_PVAL_DIFF))) + geom_point()


#####################old stuff below##############################


#Compare p-values diff
test_BH_diff_sf <- data.table::merge.data.table(
  x = diffcom_sf_log[, c(1,4,5,7)],
  y = diffcom_sf_nlog[, c(1,4,5,7)],
  by = c("LR_pair", "L_CELLTYPE", "R_CELLTYPE"),
  suffixes = c("_sf_log", "_sf_nlog"), 
  all = FALSE
)
test_BH_diff_sct <- data.table::merge.data.table(
  x = diffcom_sct_log[, c(1,4,5,7)],
  y = diffcom_sct_nlog[, c(1,4,5,7)],
  by = c("LR_pair", "L_CELLTYPE", "R_CELLTYPE"),
  suffixes = c("_sct_log", "_sct_nlog"), 
  all = FALSE
)
test_BH_diff <- data.table::merge.data.table(
  x = test_BH_diff_sf,
  y = test_BH_diff_sct,
  by = c("LR_pair", "L_CELLTYPE", "R_CELLTYPE"),
  all = FALSE
)
ggplot(test_BH_diff_sf, aes(x = BH_pvals_diff_sf_log, y = BH_pvals_diff_sf_nlog)) + geom_point()
ggplot(test_BH_diff_sct, aes(x = BH_pvals_diff_sct_log, y = BH_pvals_diff_sct_nlog)) + geom_point()
ggplot(test_BH_diff, aes(x = BH_pvals_diff_sct_log, y = BH_pvals_diff_sf_log)) + geom_point()
ggplot(test_BH_diff, aes(x = BH_pvals_diff_sct_nlog, y = BH_pvals_diff_sf_nlog)) + geom_point()


cols <- names(test_BH_diff)[sapply(test_BH_diff, is.numeric)]
test_BH_sig <- test_BH_diff[, lapply(.SD, function(x) { ifelse(x<=0.05, TRUE, FALSE)}), .SDcols = cols]
table(test_BH_sig$BH_pvals_diff_sf_log, test_BH_sig$BH_pvals_diff_sf_nlog)
table(test_BH_sig$BH_pvals_diff_sf_log, test_BH_sig$BH_pvals_diff_sct_log)

#venn diagram here
venn.diagram(
  x = list(
    sf_log = which(test_BH_sig$BH_pvals_diff_sf_log),
    sf_nlog = which(test_BH_sig$BH_pvals_diff_sf_nlog),
    sct_log = which(test_BH_sig$BH_pvals_diff_sct_log),
    sct_nlog = which(test_BH_sig$BH_pvals_diff_sct_nlog)
  ),
  filename = "../../../../../diff_BH.tiff"#,
  #print.mode = "percent"
)

sum(test_BH_sig$BH_pvals_diff_sf_log)
sum(test_BH_sig$BH_pvals_diff_sf_nlog)
sum(test_BH_sig$BH_pvals_diff_sct_log)
sum(test_BH_sig$BH_pvals_diff_sct_nlog)

#only keep CCI which are detected in pass the specificity p-value
diff_thres <- log(1.)
spec_thres <- 0.05
diffcom_sf_log$LR_diff <- diffcom_sf_log$LR_score_old - diffcom_sf_log$LR_score_young
diffcom_sf_log[, sig := ifelse(BH_pvals_diff > 0.05, FALSE, TRUE)]
diffcom_sf_log[, keep := ifelse(abs(LR_diff) < diff_thres |
                                  BH_pvals_diff > 0.05 |
                                  (BH_pvals_young > spec_thres  & BH_pvals_old > spec_thres ),
                                "flat", ifelse(LR_diff < 0, "down", "up"))]

diffcom_sf_log_1s$LR_diff <- diffcom_sf_log_1s$LR_score_old - diffcom_sf_log_1s$LR_score_young
diffcom_sf_log_1s[, sig := ifelse(BH_pvals_diff > 0.05, FALSE, TRUE)]
diffcom_sf_log_1s[, keep := ifelse(abs(LR_diff) < diff_thres |
                                     BH_pvals_diff > 0.05 |
                                     (BH_pvals_young > spec_thres  & BH_pvals_old > spec_thres ),
                                   "flat", ifelse(LR_diff < 0, "down", "up"))]

diffcom_sct_log$LR_diff <- diffcom_sct_log$LR_score_old - diffcom_sct_log$LR_score_young
diffcom_sct_log[, keep := ifelse(abs(LR_diff) < diff_thres |
                                   BH_pvals_diff > 0.05 | 
                                   (BH_pvals_young > spec_thres  & BH_pvals_old > spec_thres ),
                                 "flat", ifelse(LR_diff < 0, "down", "up"))]
diffcom_sf_nlog$LR_fc <- log(diffcom_sf_nlog$LR_score_old/diffcom_sf_nlog$LR_score_young)
diffcom_sf_nlog[, keep := ifelse(abs(LR_fc) < diff_thres |
                                   BH_pvals_diff > 0.05 |
                                   (BH_pvals_young > spec_thres  & BH_pvals_old > spec_thres ),
                                 "flat", ifelse(LR_fc < 0, "down", "up"))]
diffcom_sct_nlog$LR_fc <- log(diffcom_sct_nlog$LR_score_old/diffcom_sct_nlog$LR_score_young)
diffcom_sct_nlog[, keep := ifelse(abs(LR_fc) < diff_thres |
                                    BH_pvals_diff > 0.05 | 
                                    (BH_pvals_young > spec_thres  & BH_pvals_old > spec_thres ),
                                  "flat", ifelse(LR_fc < 0, "down", "up"))]


table(diffcom_sf_log$keep)
table(diffcom_sf_nlog$keep)
table(diffcom_sct_log$keep)
table(diffcom_sct_nlog$keep)

diffcom_sf_log[, cci := paste(LR_pair, L_CELLTYPE, R_CELLTYPE, sep = "_")]
diffcom_sf_nlog[, cci := paste(LR_pair, L_CELLTYPE, R_CELLTYPE, sep = "_")]
diffcom_sct_log[, cci := paste(LR_pair, L_CELLTYPE, R_CELLTYPE, sep = "_")]
diffcom_sct_nlog[, cci := paste(LR_pair, L_CELLTYPE, R_CELLTYPE, sep = "_")]

venn.diagram(
  x = list(
    sf_log = diffcom_sf_log[keep == "up", cci],
    sf_nlog = diffcom_sf_nlog[keep == "up", cci],
    sct_log = diffcom_sct_log[keep == "up", cci],
    sct_nlog = diffcom_sct_nlog[keep == "up", cci]
  ),
  filename = "../../../../../cci_up.tiff"#,
  #print.mode = "percent"
)

venn.diagram(
  x = list(
    sf_log = diffcom_sf_log[keep == "down", cci],
    sf_nlog = diffcom_sf_nlog[keep == "down", cci],
    sct_log = diffcom_sct_log[keep == "down", cci],
    sct_nlog = diffcom_sct_nlog[keep == "down", cci]
  ),
  filename = "../../../../../cci_down.tiff"#,
  #print.mode = "percent"
)





head(sort(table(diffcom_sf_log[keep == "down", LR_pair]), decreasing = TRUE), 10)
head(sort(table(diffcom_sf_nlog[keep == "down", LR_pair]), decreasing = TRUE), 10)
head(sort(table(diffcom_sct_log[keep == "down", LR_pair]), decreasing = TRUE), 10)
head(sort(table(diffcom_sct_nlog[keep == "down", LR_pair]), decreasing = TRUE), 10)

head(sort(table(diffcom_sf_log[keep == "down", L_CELLTYPE]), decreasing = TRUE), 10)
head(sort(table(diffcom_sf_nlog[keep == "down", L_CELLTYPE]), decreasing = TRUE), 10)
head(sort(table(diffcom_sct_log[keep == "down", L_CELLTYPE]), decreasing = TRUE), 10)
head(sort(table(diffcom_sct_nlog[keep == "down", L_CELLTYPE]), decreasing = TRUE), 10)

head(sort(table(diffcom_sf_log[keep == "up", LR_pair]), decreasing = TRUE), 10)
head(sort(table(diffcom_sf_nlog[keep == "up", LR_pair]), decreasing = TRUE), 10)
head(sort(table(diffcom_sct_log[keep == "up", LR_pair]), decreasing = TRUE), 10)
head(sort(table(diffcom_sct_nlog[keep == "up", LR_pair]), decreasing = TRUE), 10)


#compare 2-sided vs 1-sided BH-diff p-values
ggplot(data.frame(x = diffcom_sf_log$BH_pvals_diff, y = diffcom_sf_log_1s$BH_pvals_diff ), aes(x=x, y=y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
ggplot(data.frame(x = diffcom_sf_log$pvals_diff, y = diffcom_sf_log_1s$pvals_diff ), aes(x=x, y=y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

table(diffcom_sf_log$sig)
table(diffcom_sf_log_1s$sig)
