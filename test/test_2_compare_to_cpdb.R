####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - July 2020
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
#   iterations = 10,
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
#   LR_data = LRall[cpdb == TRUE, c("GENESYMB_L", "GENESYMB_R", "SYMB_LR")],
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

## Read CPDB results and preprocess ####

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
#only keep monodimer ligand-receptor interactions for this particular comparison
cpdb_test_2 <- cpdb_test_2[is.na(cpdb_test_2$name_a_2) & is.na(cpdb_test_2$name_b_2),]
cpdb_test_2 <- merge.data.table(cpdb_test_2, orthologs_test_2, by.x = "name_a_1", by.y = "human_symbol", all.x = TRUE)
cpdb_test_2 <- merge.data.table(cpdb_test_2, orthologs_test_2, by.x = "name_b_1", by.y = "human_symbol", all.x = TRUE)
cpdb_test_2$LR_GENES <- paste(cpdb_test_2$mouse_symbol.x, cpdb_test_2$mouse_symbol.y, sep = "_")
#common LR pairs
common_pairs <- intersect(unique(cpdb_test_2$LR_GENES), unique(diffcom_test_2$LR_GENES))
#subset the results
diffcom_test_2 <- diffcom_test_2[LR_GENES %in% common_pairs]
diffcom_test_2$L_CELLTYPE <- gsub(" ", ".", diffcom_test_2$L_CELLTYPE)
diffcom_test_2$R_CELLTYPE <- gsub(" ", ".", diffcom_test_2$R_CELLTYPE)
cpdb_test_2 <- cpdb_test_2[cpdb_test_2$LR_GENES %in% common_pairs, ]
cpdb_test_2 <- cpdb_test_2[, c("LR_GENES", "cell_type_a", "cell_type_b", "pvalue_old", "pvalue_young",
                           "score_old", "score_young")]
#create comparison dt
setDT(cpdb_test_2)
comparison_test_2 <- merge.data.table(
  x = diffcom_test_2,
  y = cpdb_test_2,
  by.x = c("LR_GENES", "L_CELLTYPE", "R_CELLTYPE"),
  by.y = c("LR_GENES", "cell_type_a", "cell_type_b"),
  all.x = TRUE
)
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





## Conclusion #### 

# We can observe a very good matching between the two methods.
# The scores are (as intended) exactly the sames. 
# The p-values are similar. The differences come from 
# the randomness of the permutation test (and we have only used 
# 1000 permutations) and also because we consider more cell-types
# in scDiffCom.
# All in all, the classification of significant/non-significant CCIs is
# similar at 98-99% between the two methods.