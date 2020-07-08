####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - June 2020
##
## Test how the permutation tests of scDiffCom
## behaves in various terms. Compare one-tailed
## to two-tailed p-values, results from log-normalized
## data vs nonlog-normalized data, etc.
##
####################################################
##

library(Seurat)
library(scDiffCom)
library(data.table)
library(ggplot2)
library(e1071)
#library(profvis)
#library(microbenchmark)

#load a Seurat objects
#seurat_tms_test_file.rds corresponds to the Liver tissue from TMS FACS data.
seurat_test <- readRDS("../data_seurat_example.rds")
#seurat_test <- readRDS("../../../../../seurat_droplet_example_kidney.rds")
seurat_test <- NormalizeData(seurat_test, assay = "RNA")
seurat_test$age_group <- ifelse(seurat_test$age %in% c('1m', '3m'), 'young', 'old' )
seurat_test$cell_ontology_class <- as.character(seurat_test$cell_ontology_class)

#load LR database from scDiffCom
LR_test <- LRall$LRall_one2one
LR_test <- LR_test[LR_test$sCsR, c(2,3,1)]

#run scDiffCom (in non-log space) and return the distributions from the permutations
start_time <- Sys.time()
diffcom_distr_test <- run_diffcom(
  seurat_obj = seurat_test,
  LR_data = LR_test,
  seurat_cell_type_id = "cell_ontology_class",
  condition_id = "age_group",
  assay = "RNA",
  slot = "data",
  log_scale = FALSE,
  min_cells = 5,
  threshold = 0.1,
  permutation_analysis = TRUE,
  iterations = 1000,
  return_distr = TRUE
)
end_time <- Sys.time()
end_time - start_time
#in log-space
diffcom_distr_test_log <- run_diffcom(
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
  iterations = 1000,
  return_distr = TRUE
)


#qqplots of some CCI
qqnorm(diffcom_distr_test$distr_diff[1234,], pch = 1, frame = FALSE)
qqline(diffcom_distr_test$distr_diff[1234,], col = "steelblue", lwd = 2)
hist(diffcom_distr_test$distr_diff[1234,], breaks = 50)
shapiro.test(diffcom_distr_test$distr_diff[1234,])$p.value
skewness(diffcom_distr_test$distr_diff[1234,])
#
qqnorm(diffcom_distr_test_log$distr_diff[1234,], pch = 1, frame = FALSE)
qqline(diffcom_distr_test_log$distr_diff[1234,], col = "steelblue", lwd = 2)
hist(diffcom_distr_test_log$distr_diff[1234,], breaks = 50)
shapiro.test(diffcom_distr_test_log$distr_diff[1234,1:5000])$p.value
skewness(diffcom_distr_test_log$distr_diff[1234,])

qqnorm(diffcom_distr_test$distr_diff[18,], pch = 1, frame = FALSE)
qqline(diffcom_distr_test$distr_diff[18,], col = "steelblue", lwd = 2)
hist(diffcom_distr_test$distr_diff[18,], breaks = 50)
shapiro.test(diffcom_distr_test$distr_diff[18,])$p.value
skewness(diffcom_distr_test$distr_diff[18,])
#
qqnorm(diffcom_distr_test_log$distr_diff[18,], pch = 1, frame = FALSE)
qqline(diffcom_distr_test_log$distr_diff[18,], col = "steelblue", lwd = 2)
hist(diffcom_distr_test_log$distr_diff[18,], breaks = 50)
shapiro.test(diffcom_distr_test_log$distr_diff[18,1:5000])$p.value
skewness(diffcom_distr_test_log$distr_diff[18,])

qqnorm(diffcom_distr_test$distr_cond1[124,], pch = 1, frame = FALSE)
qqline(diffcom_distr_test$distr_cond1[124,], col = "steelblue", lwd = 2)
hist(diffcom_distr_test$distr_cond1[124,], breaks = 50)
shapiro.test(diffcom_distr_test$distr_cond1[124,])$p.value
skewness(diffcom_distr_test$distr_cond1[124,])
#
qqnorm(diffcom_distr_test_log$distr_cond1[124,], pch = 1, frame = FALSE)
qqline(diffcom_distr_test_log$distr_cond1[124,], col = "steelblue", lwd = 2)
hist(diffcom_distr_test_log$distr_cond1[124,], breaks = 50)
shapiro.test(diffcom_distr_test_log$distr_diff[124,1:5000])$p.value
skewness(diffcom_distr_test_log$distr_diff[124,])
skewness(diffcom_distr_test_log$distr_cond1[124,])


#shapiro test of normality, most distributions are not normal according to the test, but this is not a problem
# see https://stackoverflow.com/questions/7781798/seeing-if-data-is-normally-distributed-in-r/7788452#7788452
#most distributions look OK to me!
#log-normalized data returns more normally-distributed distributions
shapiro_diff <- apply(diffcom_distr_test$distr_diff, MARGIN = 1, function(x) {
  shapiro.test(x)$p.value
})
hist(log10(shapiro_diff), breaks = 100)
sum(shapiro_diff > 0.05)/length(shapiro_diff)*100
#
shapiro_diff_log <- apply(diffcom_distr_test_log$distr_diff, MARGIN = 1, function(x) {
  shapiro.test(x)$p.value
})
hist(log10(shapiro_diff_log), breaks = 100)
sum(shapiro_diff_log > 0.05)/length(shapiro_diff_log)*100

shapiro_cond1 <- apply(diffcom_distr_test$distr_cond1, MARGIN = 1, function(x) {
  shapiro.test(x)$p.value
})
hist(log10(shapiro_cond1), breaks = 100)
sum(shapiro_cond1 > 0.05)/length(shapiro_cond1)*100
#
shapiro_cond1_log <- apply(diffcom_distr_test_log$distr_cond1, MARGIN = 1, function(x) {
  shapiro.test(x)$p.value
})
hist(log10(shapiro_cond1_log), breaks = 100)
sum(shapiro_cond1_log > 0.05)/length(shapiro_cond1_log)*100

shapiro_cond2 <- apply(diffcom_distr_test$distr_cond2, MARGIN = 1, function(x) {
  shapiro.test(x)$p.value
})
hist(log10(shapiro_cond2), breaks = 100)
sum(shapiro_cond2 > 0.05)/length(shapiro_cond2)*100
#
shapiro_cond2_log <- apply(diffcom_distr_test_log$distr_cond2, MARGIN = 1, function(x) {
  shapiro.test(x)$p.value
})
hist(log10(shapiro_cond2_log), breaks = 100)
sum(shapiro_cond2_log > 0.05)/length(shapiro_cond2_log)*100

#ratio mean variance to see if the distribution is centered
mean_sd_diff <- apply(diffcom_distr_test$distr_diff, MARGIN = 1, function(x) {
  abs(mean(x))/sd(x)
})
hist(mean_sd_diff, breaks = 50)
hist(log10(mean_sd_diff), breaks = 50)
median_sd_diff <- apply(diffcom_distr_test$distr_diff, MARGIN = 1, function(x) {
  abs(median(x))/sd(x)
})
hist(median_sd_diff, breaks = 50)
hist(log10(median_sd_diff), breaks = 50)
#
mean_sd_diff_log <- apply(diffcom_distr_test_log$distr_diff, MARGIN = 1, function(x) {
  abs(mean(x))/sd(x)
})
hist(mean_sd_diff_log, breaks = 50)
hist(log10(mean_sd_diff_log), breaks = 50)
median_sd_diff_log <- apply(diffcom_distr_test_log$distr_diff, MARGIN = 1, function(x) {
  abs(median(x))/sd(x)
})
hist(median_sd_diff_log, breaks = 50)
hist(log10(median_sd_diff_log), breaks = 50)

hist(log10(abs(apply(diffcom_distr_test$distr_diff, MARGIN = 1, mean))), breaks = 40 )
hist(log10(apply(diffcom_distr_test$distr_diff, MARGIN = 1, sd)), breaks = 40 )

#skewness of distribution
skew_diff <- apply(diffcom_distr_test$distr_diff, MARGIN = 1, function(x) {
  skewness(x)
})
hist(skew_diff, breaks = 100)
#
skew_diff_log <- apply(diffcom_distr_test_log$distr_diff, MARGIN = 1, function(x) {
  skewness(x)
})
hist(skew_diff_log, breaks = 50)
#
skew_cond1 <- apply(diffcom_distr_test$distr_cond1, MARGIN = 1, function(x) {
  skewness(x)
})
hist(skew_cond1, breaks = 100)

#compute one-tailed vs two-tailed pvalues without FDR correction
distr_diff <- diffcom_distr_test$distr_diff
p_value_diff_twoT <- rowSums(abs(distr_diff[, 1:ncol(distr_diff)]) >= abs(distr_diff[, ncol(distr_diff)])) / ncol(distr_diff)
p_value_diff_posT <- rowSums(distr_diff[, 1:ncol(distr_diff)] >= distr_diff[, ncol(distr_diff)]) / ncol(distr_diff)
p_value_diff_negT <- rowSums(distr_diff[, 1:ncol(distr_diff)] <= distr_diff[, ncol(distr_diff)]) / ncol(distr_diff)
p_value_diff_combT <- pmin(p_value_diff_posT, p_value_diff_negT)
hist(p_value_diff_combT, breaks = 50)
#
distr_diff_log <- diffcom_distr_test_log$distr_diff
p_value_diff_log_twoT <- rowSums(abs(distr_diff_log[, 1:ncol(distr_diff_log)]) >= abs(distr_diff_log[, ncol(distr_diff_log)])) / ncol(distr_diff_log)
p_value_diff_log_posT <- rowSums(distr_diff_log[, 1:ncol(distr_diff_log)] >= distr_diff_log[, ncol(distr_diff_log)]) / ncol(distr_diff_log)
p_value_diff_log_negT <- rowSums(distr_diff_log[, 1:ncol(distr_diff_log)] <= distr_diff_log[, ncol(distr_diff_log)]) / ncol(distr_diff_log)
p_value_diff_log_combT <- pmin(p_value_diff_log_posT, p_value_diff_log_negT)
hist(p_value_diff_log_combT, breaks = 50)

p_value_comp <- data.table(
  twoT = p_value_diff_twoT,
  combT = p_value_diff_combT,
  log_twoT = p_value_diff_log_twoT,
  log_combT = p_value_diff_log_combT
)

ggplot(p_value_comp, aes(x = twoT, y = log_twoT)) + 
  geom_point() +
  geom_density_2d() +
  scale_x_log10() +
  scale_y_log10()

p_value_comp[, c("twoT_sig", "combT_sig", "log_twoT_sig", "log_combT_sig") := lapply(.SD, function(x) {
  ifelse(x <= 0.05, TRUE, FALSE)
})]

table(p_value_comp$twoT_sig, p_value_comp$combT_sig)
sum(p_value_comp$twoT_sig)
sum(p_value_comp$combT_sig)

table(p_value_comp$twoT_sig, p_value_comp$log_twoT_sig)*100/nrow(p_value_comp)
sum(p_value_comp$log_twoT_sig)
table(p_value_comp$combT_sig, p_value_comp$log_combT_sig)
sum(p_value_comp$log_combT_sig)

table(p_value_comp$combT_sig, p_value_comp$log_combT_sig)*100/nrow(p_value_comp)


#FDR test and comparison
#first get the total number of CCI, namely before removing the non-detected CCI
diffcom_test <- readRDS("test_and_comparison/data_results_diffcom.rds")
total_cci <- nrow(diffcom_test)
#total_cci <- nrow(p_value_comp)
#see how including all CCI vs only the detected ones changes the FDR distrib
bh_comp_twoT <- data.table(
  p_value_diff_twoT = p_value_diff_twoT,
  bh_diff_twoT = p.adjust(p_value_diff_twoT, method = "BH"),
  bh_diff_twoT_f = p.adjust(p_value_diff_twoT, method = "BH", n = total_cci),
  p_value_diff_log_twoT = p_value_diff_log_twoT,
  bh_diff_log_twoT = p.adjust(p_value_diff_log_twoT, method = "BH"),
  bh_diff_log_twoT_f = p.adjust(p_value_diff_log_twoT, method = "BH", n = total_cci)
)
bh_comp_twoT[, c("pval_sig", "bh_sig", "bh_f_sig", "log_pval_sig", "log_bh_sig", "log_bh_f_sig") := lapply(.SD, function(x) {
  ifelse(x <= 0.05, TRUE, FALSE)
})]

bh_comp_posT <- data.table(
  p_value_diff_posT = p_value_diff_posT,
  bh_diff_posT = p.adjust(p_value_diff_posT, method = "BH"),
  bh_diff_posT_f= p.adjust(p_value_diff_posT, method = "BH", n = total_cci),
  p_value_diff_log_posT = p_value_diff_log_posT,
  bh_diff_log_posT = p.adjust(p_value_diff_log_posT, method = "BH"),
  bh_diff_log_posT_f = p.adjust(p_value_diff_log_posT, method = "BH", n = total_cci)
)
bh_comp_posT[, c("pval_sig", "bh_sig", "bh_f_sig", "log_pval_sig", "log_bh_sig", "log_bh_f_sig") := lapply(.SD, function(x) {
  ifelse(x <= 0.05, TRUE, FALSE)
})]

bh_comp_negT <- data.table(
  p_value_diff_negT = p_value_diff_negT,
  bh_diff_negT = p.adjust(p_value_diff_negT, method = "BH"),
  bh_diff_negT_f = p.adjust(p_value_diff_negT, method = "BH", n = total_cci),
  p_value_diff_log_negT = p_value_diff_log_negT,
  bh_diff_log_negT = p.adjust(p_value_diff_log_negT, method = "BH"),
  bh_diff_log_negT_f = p.adjust(p_value_diff_log_negT, method = "BH", n = total_cci)
)
bh_comp_negT[, c("pval_sig", "bh_sig", "bh_f_sig", "log_pval_sig", "log_bh_sig", "log_bh_f_sig") := lapply(.SD, function(x) {
  ifelse(x <= 0.05, TRUE, FALSE)
})]

ggplot(bh_comp_twoT) +
  geom_point(aes(x = p_value_diff_twoT, y = bh_diff_twoT)) +
  geom_point(aes(x = p_value_diff_twoT, y = bh_diff_twoT_f)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(bh_comp) +
  geom_point(aes(x = p_value_diff_log_twoT, y = bh_diff_log_twoT)) +
  geom_point(aes(x = p_value_diff_log_twoT, y = bh_diff_log_twoT_f)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0, 0.01) +
  ylim(0, 0.05)

ggplot(bh_comp_posT) +
  geom_point(aes(x = p_value_diff_posT, y = bh_diff_posT)) +
  geom_point(aes(x = p_value_diff_posT, y = bh_diff_posT_f)) +
  geom_abline(intercept = 0, slope = 1)


table(bh_comp_twoT$pval_sig, bh_comp_twoT$bh_sig)
table(bh_comp_twoT$bh_sig, bh_comp_twoT$bh_f_sig)

table(bh_comp_posT$pval_sig, bh_comp_posT$bh_sig)
table(bh_comp_posT$bh_sig, bh_comp_posT$bh_f_sig)

table(bh_comp_negT$pval_sig, bh_comp_negT$bh_sig)
table(bh_comp_negT$bh_sig, bh_comp_negT$bh_f_sig)


table(bh_comp_twoT$log_bh_sig, bh_comp_twoT$log_bh_f_sig)
table(bh_comp_posT$log_bh_sig, bh_comp_posT$log_bh_f_sig)
table(bh_comp_negT$log_bh_sig, bh_comp_negT$log_bh_f_sig)

ggplot(bh_comp_negT) +
  geom_point(aes(x = p_value_diff_log_negT, y = bh_diff_log_negT)) +
  geom_point(aes(x = p_value_diff_log_negT, y = bh_diff_log_negT_f)) +
  geom_abline(intercept = 0, slope = 1)
sum(bh_comp_negT$log_bh_sig)
head(sort(bh_comp_negT$bh_diff_log_negT_f))
