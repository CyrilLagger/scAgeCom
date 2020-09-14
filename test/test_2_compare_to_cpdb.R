####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - September 2020
##
## Check that scDiffCom returns the same values as
## CellPhoneDB when used in the same conditions.
##
####################################################
##

## Libraries ####

library(Seurat)
library(scDiffCom)
library(data.table)
library(ggplot2)

## Data path ####
dir_data <- "../data_scAgeCom/"

## Load a Seurat objects ####

#here the file corresponds to the Liver tissue from TMS FACS data.
seurat_test_2 <- readRDS(paste0(dir_data, "data_seurat_example.rds"))
seurat_test_2 <- NormalizeData(seurat_test_2, assay = "RNA")
seurat_test_2$age_group <- ifelse(seurat_test_2$age %in% c('1m', '3m'), 'young', 'old' )
seurat_test_2$cell_ontology_class <- as.character(seurat_test_2$cell_ontology_class)

## Load the LR database from scDiffCom ####
LR_test_2 <- scDiffCom::LR6db$LR6db_curated
LR_test_2 <- LR_test_2[grepl("CELLPHONEDB", DATABASE)]

## Run CellPhoneDB ####

#Note: the function scDiffCom::run_cpdb_from_seurat is a R wrapper
#      around the original CPDB python package. It only works on
#      Linux and if CPDB is installed in your current python environment.
#      This function is probably not runnning anymore and should be modified
#      consequently. However we have stored the data from one example case
#      (Liver-TMS-FACS in non-log scale). We can use it for comparison.

#typical usage (for the parameters see the CPDB manual: https://github.com/Teichlab/cellphonedb)
#the output is the typical CPDB output with various files
# run_cpdb_from_seurat(
#   seurat_obj = seurat_test_2,
#   assay = "RNA",
#   slot = "data",
#   log_scale = FALSE,
#   seurat_cell_type_id = "cell_ontology_class",
#   min_cells = 5,
#   condition_id = "age_group",
#   input_dir = paste0(getwd(),"/cpdb_example"),
#   create_plots = FALSE,
#   method = 'statistical_analysis',
#   iterations = 1000,
#   threshold = NULL,
#   result_precision = NULL,
#   counts_data = "gene_name",
#   output_format = NULL,
#   verbose = TRUE,
#   debug_seed = NULL,
#   threads = NULL,
#   subsampling = FALSE,
#   subsampling_log = FALSE,
#   subsampling_num_pc = NULL,
#   subsampling_num_cells = NULL,
#   return_full_dt = TRUE
# )


## Create (or Read) the equivalent scDiffCom result file ####

#usage to create the file with correct parameters.
#no need to run if data already stored
# start_time <- Sys.time()
# diffcom_test_2 <- run_diffcom(
#   seurat_object = seurat_test_2,
#   LR_data = LR_test_2,
#   seurat_cell_type_id = "cell_ontology_class",
#   condition_id = "age_group",
#   assay = "RNA",
#   slot = "data",
#   log_scale = FALSE,
#   min_cells = 5,
#   threshold = 0.1,
#   permutation_analysis = TRUE,
#   one_sided = FALSE,
#   iterations = 1000,
#   return_distr = FALSE
# )
# end_time <- Sys.time()
# end_time - start_time

#saveRDS(object = diffcom_test_2, file = paste0(dir_data, "test/test_2_data_scDiffcom_for_cpdb.rds"))
diffcom_test_2 <- readRDS(paste0(dir_data, "test/test_2_data_scDiffcom_for_cpdb.rds"))

## Read CPDB results and process them for the comparison ####

#Note: the output of CPDB requires some processing to be compared
#      to the output of scDiffCom. Here we use a file that has 
#      already been partially preprocessed.

cpdb_test_2 <- read.table(
  paste0(dir_data, "test/test_2_data_cpdb_results.txt"),
  header = TRUE,
  sep = "\t"
)
orthologs_test_2 <- read.table(
  file = paste0(dir_data, "test/test_2_data_cpdb_orthologs.txt"),
  header = TRUE,
  sep = "\t"
)
setDT(cpdb_test_2)
setDT(orthologs_test_2)
cpdb_test_2[, c("mouse_a_1", "mouse_a_2", "mouse_b_1", "mouse_b_2") := list(
              orthologs_test_2[.SD, on = "human_symbol==name_a_1", x.mouse_symbol],
              orthologs_test_2[.SD, on = "human_symbol==name_a_2", x.mouse_symbol],
              orthologs_test_2[.SD, on = "human_symbol==name_b_1", x.mouse_symbol],
              orthologs_test_2[.SD, on = "human_symbol==name_b_2", x.mouse_symbol])]
cpdb_test_2 <- cpdb_test_2[!(!is.na(name_a_2) & is.na(mouse_a_2)) | !(!is.na(name_b_2) & is.na(mouse_b_2))]
cpdb_test_2[, c("LR_ID") := list(sapply(1:nrow(.SD), function(i) {
  temp <- c(mouse_a_1[[i]], mouse_a_2[[i]], mouse_b_1[[i]], mouse_b_2[[i]])
  temp <- temp[!is.na(temp)]
  temp <- sort(temp)
  temp <- paste0(temp, collapse = "_")
}))]
cpdb_test_2[, CCI := paste(cell_type_a, cell_type_b, LR_ID, sep = "_")]
cpdb_test_2[, CCI2 := paste(cell_type_b, cell_type_a, LR_ID, sep = "_")]
cpdb_test_2 <- cpdb_test_2[, c("CCI", "CCI2", "LR_ID", "cell_type_a", "cell_type_b",
                               "mouse_a_1", "mouse_a_2", "mouse_b_1", "mouse_b_2",
                               "score_old", "pvalue_old", "score_young", "pvalue_young",
                               "mean_a_1_old", "mean_a_2_old", "mean_b_1_old", "mean_b_2_old",
                               "mean_a_1_young", "mean_a_2_young", "mean_b_1_young", "mean_b_2_young")]


## Process scDiffCom results for comparison ####
diffcom_test_2$L_CELLTYPE <- gsub(" ", ".", diffcom_test_2$L_CELLTYPE)
diffcom_test_2$R_CELLTYPE <- gsub(" ", ".", diffcom_test_2$R_CELLTYPE)
diffcom_test_2[, CCI := paste(L_CELLTYPE, R_CELLTYPE, LR_SORTED, sep = "_")]
diffcom_test_2_old <- diffcom_test_2[LR_DETECTED_old == TRUE,]
diffcom_test_2_young <- diffcom_test_2[LR_DETECTED_young == TRUE,]
diffcom_test_2_mix <- diffcom_test_2[LR_DETECTED_young == TRUE | LR_DETECTED_old == TRUE,]

## Consider the common CCI between the two packages ####
common_CCI_old <- intersect(unique(cpdb_test_2$CCI), unique(diffcom_test_2_old$CCI))
common_CCI_young <- intersect(unique(cpdb_test_2$CCI), unique(diffcom_test_2_young$CCI))
common_CCI_mix <- intersect(unique(cpdb_test_2$CCI), unique(diffcom_test_2_mix$CCI))

## Merge the results together

test <- merge.data.table(
  cpdb_test_2[ , c("CCI", "mouse_a_1", "mouse_a_2", "mean_a_1_old", "mean_a_2_old", "mean_b_1_old", "mean_b_2_old")],
  diffcom_test_2[, c("CCI", "L1_EXPRESSION_old", "L2_EXPRESSION_old", "R1_EXPRESSION_old", "R2_EXPRESSION_old")],
  by = "CCI",
  all.x = TRUE
)


cpdb_test_2_com <- cpdb_test_2[CCI %in% common_CCI]
diffcom_test_2_com <- diffcom_test_2[CCI %in% common_CCI]

cpdb_test_2_com[duplicated(cpdb_test_2_com$CCI)]$CCI

comp_dt <- cpdb_test_2_com[diffcom_test_2_com, on = "CCI"]
comp_dt2 <- comp_dt[cpdb_test_2, on = "CCI==CCI2"]

testxk <- comp_dt2[, c("CCI", "CCI2", "i.score_old", "score_old", "LR_SCORE_old")]
testxk[, diff1 := abs(LR_SCORE_old - score_old)]
testxk[, diff2 := abs(LR_SCORE_old - i.score_old)]
testxk[, mindif := pmin(diff1, diff2)]
hist(log10(testxk$mindif), breaks = 100)

testxk2 <- comp_dt2[, c("CCI", "CCI2", "i.pvalue_old", "pvalue_old", "PVAL_old")]
testxk2[, diff1 := abs(PVAL_old - pvalue_old)]
testxk2[, diff2 := abs(PVAL_old - i.pvalue_old)]
testxk2[, mindif := pmin(diff1, diff2)]
hist(log10(testxk2$mindif), breaks = 100)


plot(comp_dt$LR_SCORE_old, comp_dt$score_old)
plot(comp_dt$pvalue_old, comp_dt$PVAL_old)
plot(comp_dt$pvalue_young, comp_dt$PVAL_young)

test <- comp_dt[, c("CCI", "score_old", "LR_SCORE_old")]
test[, diff:= abs(score_old - LR_SCORE_old)]

test_w <- subset(seurat_test_2, subset = cell_ontology_class == "hepatocyte" & age_group == "old")
test_w2 <- subset(seurat_test_2, subset = cell_ontology_class == "endothelial cell of hepatic sinusoid" & age_group == "old")
(mean(expm1(test_w$RNA@data["Fn1",]))  + min(mean(expm1(test_w2$RNA@data["Itga3",])), mean(expm1(test_w2$RNA@data["Itgb1",]))))/2
(mean(expm1(test_w2$RNA@data["Cd74",]))  + mean(expm1(test_w$RNA@data["Ccr3",])))/2

test_y <- subset(seurat_test_2, subset = cell_ontology_class == "NK cell" & age_group == "old")
test_y2 <- subset(seurat_test_2, subset = cell_ontology_class == "B cell" & age_group == "old")
(mean(expm1(test_y$RNA@data["Ccl5",]))  + mean(expm1(test_y2$RNA@data["Gpr75",])))/2
(mean(expm1(test_y$RNA@data["Ccl5",]))  + min(mean(expm1(test_w2$RNA@data["Itga3",])), mean(expm1(test_w2$RNA@data["Itgb1",]))))/2
(mean(expm1(test_w2$RNA@data["Ccl5",]))  + mean(expm1(test_w$RNA@data["Ccr3",])))/2

comparison_test_2 <- na.omit(comparison_test_2)
comparison_test_2$sig_young_cpdb <- ifelse(comparison_test_2$pvalue_young <= 0.05, TRUE, FALSE)
comparison_test_2$sig_young_diffcom <- ifelse(comparison_test_2$PVAL_young <= 0.05, TRUE, FALSE)
comparison_test_2$sig_old_cpdb <- ifelse(comparison_test_2$pvalue_old <= 0.05, TRUE, FALSE)
comparison_test_2$sig_old_diffcom <- ifelse(comparison_test_2$PVAL_old <= 0.05, TRUE, FALSE)

## Compare LR scores for detected CCI ####

g_score_old <- ggplot(comparison_test_2[LR_DETECTED_old == TRUE], aes(x = LR_SCORE_old, y = score_old)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10() +
  xlab(expression(xi["old,scdiffcom"])) +
  ylab(expression(xi["old,cpdb"])) +
  theme(text=element_text(size=20))
g_score_young <- ggplot(comparison_test_2[LR_DETECTED_young == TRUE], aes(x = LR_SCORE_young, y = score_young)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10() +
  xlab(expression(xi["young,scdiffcom"])) +
  ylab(expression(xi["young,cpdb"])) +
  theme(text=element_text(size=20)) 

g_score <- cowplot::plot_grid(
  g_score_young, g_score_old,
  ncol = 2
)
#ggsave(filename = paste0(dir_data, "test/test_2_comp_cpdb_LR_score.png"), plot = g_score, scale = 1.5)

## Compare specificity p-values for detected CCI ####

g_spec_p_old <- ggplot(comparison_test_2[LR_DETECTED_old == TRUE], aes(x = PVAL_old, y = pvalue_old)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlab(expression(p["old,scdiffcom"])) +
  ylab(expression(p["old,cpdb"])) +
  theme(text=element_text(size=20)) 
g_spec_p_young <- ggplot(comparison_test_2[LR_DETECTED_young == TRUE], aes(x = PVAL_young, y = pvalue_young)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlab(expression(p["young,scdiffcom"])) +
  ylab(expression(p["young,cpdb"])) +
  theme(text=element_text(size=20))
g_spec_p <- cowplot::plot_grid(
  g_spec_p_young, g_spec_p_old,
  ncol = 2
)
#ggsave(filename = paste0(dir_data, "test/test_2_comp_cpdb_pvalue.png"), plot = g_spec_p, scale = 1.5)

## Compare significant p-values ####

#for detected CCI only
table(comparison_test_2[LR_DETECTED_young == TRUE]$sig_young_cpdb, comparison_test_2[LR_DETECTED_young == TRUE]$sig_young_diffcom)
table(comparison_test_2[LR_DETECTED_old == TRUE]$sig_old_cpdb, comparison_test_2[LR_DETECTED_old == TRUE]$sig_old_diffcom)
#for all CCI, less relevant
table(comparison_test_2$sig_young_cpdb, comparison_test_2$sig_young_diffcom)
table(comparison_test_2$sig_old_cpdb, comparison_test_2$sig_old_diffcom)


(mean(expm1(test_w2$RNA@data["B2m",])))
(mean((test_w2$RNA@data["Cd74",])))

## Conclusion #### 

# We can observe a very good matching between the two methods.
# The scores are (as intended) exactly the sames. 
# The p-values are similar. The differences come from 
# the randomness of the permutation test (and we have only used 
# 1000 permutations) and also because we consider more cell-types
# in scDiffCom.
# All in all, the classification of significant/non-significant CCIs is
# similar at 98-99% between the two methods.