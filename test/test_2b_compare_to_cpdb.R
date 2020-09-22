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
dir_data_test <- "../data_scAgeCom/test/"

## Load a Seurat objects ####

#here the file corresponds to the Liver tissue from TMS FACS data.
seurat_t2 <- readRDS(paste0(dir_data_test, "data_seurat_example_facs_liver.rds"))
seurat_t2 <- NormalizeData(seurat_t2, assay = "RNA")
seurat_t2$age_group <- ifelse(seurat_t2$age %in% c('1m', '3m'), 'young', 'old' )
seurat_t2$cell_ontology_class <- as.character(seurat_t1$cell_ontology_class)

## Load interactions from CELLPHONEDB ####
LR_t2 <- scDiffCom::LR6db$LR6db_curated
LR_t2 <- LR_t2[grepl("CELLPHONEDB", DATABASE)]

## Create (or Read) the scDiffCom result file ####

#usage to create the file with correct parameters.
#no need to run if data already stored
# diffcom_t2 <- run_diffcom(
#   seurat_object = seurat_t2,
#   LR_data = LR_t2,
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

#saveRDS(object = diffcom_t2, file = paste0(dir_data_test, "t2_data_scDiffcom_for_cpdb.rds"))
diffcom_t2 <- readRDS(paste0(dir_data_test, "t2_data_scDiffcom_for_cpdb.rds"))

## Read CPDB results and process them for the comparison ####
cpdb_t2 <- read.table(
  paste0(dir_data_test, "cpdb_test_results/cpdb_full_table_withCond.txt"),
  header = TRUE,
  sep = "\t"
)
data.table::setDT(cpdb_t2)

## process the data.tables for comparison and create a single data.table #####

#add interaction detail to the results
cpdb_ids <- scDiffCom:::prepare_LR_cpdb(keep_id = TRUE)
cpdb_t2 <- cpdb_t2[, c("LR_SORTED") := cpdb_ids[.SD, on = "id_cp_interaction", x.LR_SORTED]]
rm(cpdb_ids)
cpdb_t2[, CCI := paste(cell_type_pair, LR_SORTED, sep = ".")]

diffcom_t2[, L_CELLTYPE2 := gsub("\\-|\\_|\\,| ", "\\.", L_CELLTYPE)]
diffcom_t2[, R_CELLTYPE2 := gsub("\\-|\\_|\\,| ", "\\.", R_CELLTYPE)]
diffcom_t2[, cell_type_pair := paste(L_CELLTYPE2, R_CELLTYPE2, sep = ".")]
diffcom_t2[, cell_type_pair2 := paste(R_CELLTYPE2, L_CELLTYPE2, sep = ".")]
diffcom_t2[, CCI := paste(cell_type_pair, LR_SORTED, sep = ".")]
diffcom_t2[, CCI2 := paste(cell_type_pair2, LR_SORTED, sep = ".")]

#check for duplicated cci
sum(duplicated(diffcom_t2$CCI)) #should be zero, no duplicate returned by scDiffCom

sum(duplicated(cpdb_t2$CCI)) #duplicates because of NA (interaction not in the LR db after orthology conversion)
cpdb_t2 <- cpdb_t2[!is.na(LR_SORTED)]
sum(duplicated(cpdb_t2$CCI))

#create comparison table
common_cci <- intersect(unique(diffcom_t2$CCI), unique(cpdb_t2$CCI))
comparison_t2 <- data.table::merge.data.table(
  diffcom_t2[CCI %in% common_cci],
  cpdb_t2[CCI %in% common_cci],
  by = "CCI"
)

## Compare new scDiffCom results and CPDB results ####

#we can see that there are some inversions in the interactions
plot(comparison_t2[LR_DETECTED_old==TRUE]$LR_SCORE_old, comparison_t2[LR_DETECTED_old == TRUE]$means_old)

#we have already checked in test_1_scDiffCom that the results of scDiffCom are consistent
# the two options are that some LR interactions are not ordered in the same way or a problem with the scores in CPDB

#we try to reorder the pairs
comparison_t2[, means_old2 := comparison_t2[.SD, on = "CCI==CCI2", x.means_old]]
comparison_t2[, means_young2 := comparison_t2[.SD, on = "CCI==CCI2", x.means_young]]
comparison_t2[, pvalues_old2 := comparison_t2[.SD, on = "CCI==CCI2", x.pvalues_old]]
comparison_t2[, pvalues_young2 := comparison_t2[.SD, on = "CCI==CCI2", x.pvalues_young]]

comparison_t2[, means_mindiff_old := pmin(abs(LR_SCORE_old - means_old), abs(LR_SCORE_old - means_old2))]
comparison_t2[, means_mindiff_young := pmin(abs(LR_SCORE_young - means_young), abs(LR_SCORE_young - means_young2))]

comparison_t2[, pvalues_mindiff_old := pmin(abs(PVAL_old - pvalues_old), abs(PVAL_old - pvalues_old2))]
comparison_t2[, pvalues_mindiff_young := pmin(abs(PVAL_young - pvalues_young), abs(PVAL_young - pvalues_young2))]

#we should get very small differences for the scors and pvalues
#we still observe some CCI with quite high differences
hist(log10(comparison_t2[LR_DETECTED_old == TRUE, means_mindiff_old]), breaks = 100)
hist(log10(comparison_t2[LR_DETECTED_old == TRUE, pvalues_mindiff_old]), breaks = 100)

hist(log10(comparison_t2[LR_DETECTED_young == TRUE, means_mindiff_young]), breaks = 100)
hist(log10(comparison_t2[LR_DETECTED_young == TRUE, pvalues_mindiff_young]), breaks = 100)

#let us look at only the non-subunit interactions
#it looks much better actually
rows_to_keep <- lengths(regmatches(comparison_t2$LR_SORTED.x, gregexpr("_", comparison_t2$LR_SORTED.x))) == 1

hist(log10(comparison_t2[LR_DETECTED_old == TRUE & rows_to_keep == TRUE, means_mindiff_old]), breaks = 100)
hist(log10(comparison_t2[LR_DETECTED_old == TRUE & rows_to_keep == TRUE, pvalues_mindiff_old]), breaks = 100)

hist(log10(comparison_t2[LR_DETECTED_young == TRUE & rows_to_keep == TRUE, means_mindiff_young]), breaks = 100)
hist(log10(comparison_t2[LR_DETECTED_young == TRUE & rows_to_keep == TRUE, pvalues_mindiff_young]), breaks = 100)

#so let us look at these strange with-subunit cases and recompute them from scratch
head(comparison_t2[LR_DETECTED_old == TRUE, c("CCI", "means_mindiff_old", "pvalues_mindiff_old")][order(-means_mindiff_old)])
head(comparison_t2[LR_DETECTED_young == TRUE, c("CCI", "means_mindiff_young", "pvalues_mindiff_young")][order(-means_mindiff_young)])

scratch_cci <- "T.cell.T.cell.Cd8a_Cd8b1_Lck"
scratch_seurat1 <- subset(seurat_t2, subset = cell_ontology_class == "T cell" & age_group == "old")
scratch_seurat2 <- subset(seurat_t2, subset = cell_ontology_class == "T cell" & age_group == "old")
scratch_LR_score <- (mean(expm1(scratch_seurat1$RNA@data["Lck",]))  +
                       pmin(mean(expm1(scratch_seurat2$RNA@data["Cd8a",])), mean(expm1(scratch_seurat2$RNA@data["Cd8b1",]))))/2
#we can see that the scratch result corresponds to scDiffCom but not CELLPHONEDB
scratch_LR_score 
comparison_t2[CCI == scratch_cci, "LR_SCORE_old"]
comparison_t2[CCI == scratch_cci, "means_old"]
#strangely, the CELLPHONEDB result corresponds to the case when we take the max of the two receptors, or only one of them
#it is like one the gene is not taken into account. 
scratch_wrong_score <- (mean(expm1(scratch_seurat1$RNA@data["Lck",]))  +
                       pmax(mean(expm1(scratch_seurat2$RNA@data["Cd8a",])), mean(expm1(scratch_seurat2$RNA@data["Cd8b1",]))))/2
scratch_wrong_score

#We can conclude that there is a problem in the CELLPHONEDB results, but that the scDiffCom results are consistent
# We can observe a good matching between the two methods for simple interactions.
