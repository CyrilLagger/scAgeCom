####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - June 2020
##
## Compare the results from CellPhoneDB to 
## scDiffCom for the example Seurat object.
##
####################################################
##

library(data.table)
library(ggplot2)

#read results of both packages
#both results have been obtained from the liver tissue of TMS FACS data, size-factor normalization, not log-normalized
#both results obtained from 1000 iterations for the shake of comparison
diffcom_test <- readRDS("test_and_comparison/data_results_diffcom.rds")
cpdb_test <- read.table("test_and_comparison/data_results_cpdb.txt",
                        header = TRUE,
                        sep = "\t")

#convert the two tables to comparable tables
orthologs <- read.table(file = "test_and_comparison/data_ortholgos_cpdb.txt",
                    header = TRUE,
                    sep = "\t")
#only keep monodimer ligand-receptor interactions for this particular comparison
cpdb_test <- cpdb_test[is.na(cpdb_test$name_a_2) & is.na(cpdb_test$name_b_2),]
cpdb_test <- merge.data.table(cpdb_test, orthologs, by.x = "name_a_1", by.y = "human_symbol", all.x = TRUE)
cpdb_test <- merge.data.table(cpdb_test, orthologs, by.x = "name_b_1", by.y = "human_symbol", all.x = TRUE)
cpdb_test$LR_GENES <- paste(cpdb_test$mouse_symbol.x, cpdb_test$mouse_symbol.y, sep = "_")
#common LR pairs
common_pairs <- intersect(unique(cpdb_test$LR_GENES), unique(diffcom_test$LR_GENES))
#subset the results tables and reorder the columns
diffcom_test <- diffcom_test[LR_GENES %in% common_pairs]
diffcom_test <- diffcom_test[, c(1, 4,5,6,7,8,9, 10, 11, 12,13,14,15)]
diffcom_test$L_CELLTYPE <- gsub(" ", ".", diffcom_test$L_CELLTYPE)
diffcom_test$R_CELLTYPE <- gsub(" ", ".", diffcom_test$R_CELLTYPE)
cpdb_test <- cpdb_test[cpdb_test$LR_GENES %in% common_pairs, ]
cpdb_test <- cpdb_test[, c("LR_GENES", "cell_type_a", "cell_type_b", "pvalue_old", "pvalue_young",
                                       "score_old", "score_young")]
setDT(cpdb_test)
#create comparison dt
comp_dt <- merge.data.table(
  x = diffcom_test,
  y = cpdb_test,
  by.x = c("LR_GENES", "L_CELLTYPE", "R_CELLTYPE"),
  by.y = c("LR_GENES", "cell_type_a", "cell_type_b"),
  all.x = TRUE
)
comp_dt <- na.omit(comp_dt)
comp_dt$sig_young_cpdb <- ifelse(comp_dt$pvalue_young <= 0.05, TRUE, FALSE)
comp_dt$sig_young_diffcom <- ifelse(comp_dt$PVAL_young <= 0.05, TRUE, FALSE)
comp_dt$sig_old_cpdb <- ifelse(comp_dt$pvalue_old <= 0.05, TRUE, FALSE)
comp_dt$sig_old_diffcom <- ifelse(comp_dt$PVAL_old <= 0.05, TRUE, FALSE)
#compare LR scores for old cells (only when the CCI is detected)
ggplot(comp_dt[LR_DETECTED_old == TRUE], aes(x = LR_SCORE_old, y = score_old)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10()
#compare LR scores for young cells (only when the CCI is detected)
ggplot(comp_dt[LR_DETECTED_young == TRUE], aes(x = LR_SCORE_young, y = score_young)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10()
#compare LR pvalues for old cells (only when the CCI is detected)
ggplot(comp_dt[LR_DETECTED_old == TRUE], aes(x = PVAL_old, y = pvalue_old)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) 
#compare LR pvalues for young cells (only when the CCI is detected)
ggplot(comp_dt[LR_DETECTED_young == TRUE], aes(x = PVAL_young, y = pvalue_young)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 0.05) +
  ylim(0, 0.05)
#Show how many significant p-values for young/old  are returned by each method (only considering detected)
table(comp_dt[LR_DETECTED_young == TRUE]$sig_young_cpdb, comp_dt[LR_DETECTED_young == TRUE]$sig_young_diffcom)
table(comp_dt[LR_DETECTED_old == TRUE]$sig_old_cpdb, comp_dt[LR_DETECTED_old == TRUE]$sig_old_diffcom)
#considering all CCI, so does make less sense (because methods are not exactly similar regarding non-detected CCI) but still quite good
table(comp_dt$sig_young_cpdb, comp_dt$sig_young_diffcom)
table(comp_dt$sig_old_cpdb, comp_dt$sig_old_diffcom)



