####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - October 2020
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
library(future)

## Data path ####
dir_data_test <- "../data_scAgeCom/test/"

## Load a Seurat objects ####

#here the file corresponds to the Liver tissue from TMS FACS data.
seurat_t2 <- readRDS(paste0(dir_data_test, "data_seurat_example_facs_liver.rds"))
seurat_t2 <- NormalizeData(seurat_t2, assay = "RNA")
seurat_t2$age_group <- ifelse(seurat_t2$age %in% c('1m', '3m'), 'YOUNG', 'OLD' )
seurat_t2$cell_ontology_class <- as.character(seurat_t2$cell_ontology_class)

## Load interactions from CELLPHONEDB ####
LR_t2 <- scDiffCom::LR6db$LR6db_curated
LR_t2 <- LR_t2[grepl("CELLPHONEDB", DATABASE)]

## Create scDiffCom results with same parameters as CPDB

# future::plan(sequential)
# scdiffcom_t2 <- run_scdiffcom(
#   seurat_object = seurat_t2,
#   LR_object = LR_t2,
#   celltype_col_id = "cell_ontology_class",
#   condition_col_id = "age_group",
#   cond1_name = "YOUNG",
#   cond2_name = "OLD",
#   assay = "RNA",
#   slot = "data",
#   log_scale = FALSE,
#   min_cells = 5,
#   pct_threshold = 0.1,
#   permutation_analysis = TRUE,
#   iterations = 1000,
#   cutoff_quantile_score = 0.25,
#   cutoff_pval_specificity = 0.05,
#   cutoff_pval_de = 0.05,
#   cutoff_logfc = log(1.1),
#   return_distr = TRUE,
#   seed = 42,
#   verbose = TRUE,
#   sparse = TRUE
# )

## Save or read the scDiffCom result file ####
#saveRDS(object = scdiffcom_t2, file = paste0(dir_data_test, "t2_data_scdiffcom.rds"))
scdiffcom_t2 <- readRDS(paste0(dir_data_test, "t2_data_scDiffcom.rds"))
scdiffcom_t2_filtered <- scdiffcom_t2$scdiffcom_dt_filtered

## Read cpdb data ####

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

scdiffcom_t2_filtered[, L_CELLTYPE2 := gsub("\\-|\\_|\\,| ", "\\.", L_CELLTYPE)]
scdiffcom_t2_filtered[, R_CELLTYPE2 := gsub("\\-|\\_|\\,| ", "\\.", R_CELLTYPE)]
scdiffcom_t2_filtered[, cell_type_pair := paste(L_CELLTYPE2, R_CELLTYPE2, sep = ".")]
scdiffcom_t2_filtered[, cell_type_pair2 := paste(R_CELLTYPE2, L_CELLTYPE2, sep = ".")]
scdiffcom_t2_filtered[, CCI := paste(cell_type_pair, LR_SORTED, sep = ".")]
scdiffcom_t2_filtered[, CCI2 := paste(cell_type_pair2, LR_SORTED, sep = ".")]

#check for duplicated cci
sum(duplicated(scdiffcom_t2_filtered$CCI)) #should be zero, no duplicate returned by scDiffCom

sum(duplicated(cpdb_t2$CCI)) #duplicates because of NA (interaction not in the LR db after orthology conversion)
cpdb_t2 <- cpdb_t2[!is.na(LR_SORTED)]
sum(duplicated(cpdb_t2$CCI))

#create comparison table
common_cci <- intersect(unique(scdiffcom_t2_filtered$CCI), unique(cpdb_t2$CCI))
comparison_t2 <- data.table::merge.data.table(
  scdiffcom_t2_filtered[CCI %in% common_cci],
  cpdb_t2[CCI %in% common_cci],
  by = "CCI"
)

## Compare new scDiffCom results and CPDB results ####

#we can see that there are some inversions in the interactions
plot(comparison_t2[LR_DETECTED_OLD==TRUE]$LR_SCORE_OLD, comparison_t2[LR_DETECTED_OLD == TRUE]$means_old)

#we have already checked in test_1_scDiffCom that the results of scDiffCom are consistent
# the two options are that some LR interactions are not ordered in the same way or a problem with the scores in CPDB

#we try to reorder the pairs
comparison_t2[, means_old2 := comparison_t2[.SD, on = "CCI==CCI2", x.means_old]]
comparison_t2[, means_young2 := comparison_t2[.SD, on = "CCI==CCI2", x.means_young]]
comparison_t2[, pvalues_old2 := comparison_t2[.SD, on = "CCI==CCI2", x.pvalues_old]]
comparison_t2[, pvalues_young2 := comparison_t2[.SD, on = "CCI==CCI2", x.pvalues_young]]

comparison_t2[, means_mindiff_old := pmin(abs(LR_SCORE_OLD - means_old), abs(LR_SCORE_OLD - means_old2))]
comparison_t2[, means_mindiff_young := pmin(abs(LR_SCORE_YOUNG - means_young), abs(LR_SCORE_YOUNG - means_young2))]

comparison_t2[, pvalues_mindiff_old := pmin(abs(PVAL_OLD - pvalues_old), abs(PVAL_OLD - pvalues_old2))]
comparison_t2[, pvalues_mindiff_young := pmin(abs(PVAL_YOUNG - pvalues_young), abs(PVAL_YOUNG - pvalues_young2))]

#we should get very small differences for the scors and pvalues
#we still observe some CCI with quite high differences
hist(log10(comparison_t2[LR_DETECTED_OLD== TRUE, means_mindiff_old]), breaks = 100)
hist(log10(comparison_t2[LR_DETECTED_OLD == TRUE, pvalues_mindiff_old]), breaks = 100)

hist(log10(comparison_t2[LR_DETECTED_YOUNG == TRUE, means_mindiff_young]), breaks = 100)
hist(log10(comparison_t2[LR_DETECTED_YOUNG == TRUE, pvalues_mindiff_young]), breaks = 100)

#let us look at only the non-subunit interactions
#it looks much better actually
rows_to_keep <- lengths(regmatches(comparison_t2$LR_SORTED.x, gregexpr("_", comparison_t2$LR_SORTED.x))) == 1

hist(log10(comparison_t2[LR_DETECTED_OLD == TRUE & rows_to_keep == TRUE, means_mindiff_old]), breaks = 100)
hist(log10(comparison_t2[LR_DETECTED_OLD == TRUE & rows_to_keep == TRUE, pvalues_mindiff_old]), breaks = 100)

hist(log10(comparison_t2[LR_DETECTED_YOUNG == TRUE & rows_to_keep == TRUE, means_mindiff_young]), breaks = 100)
hist(log10(comparison_t2[LR_DETECTED_YOUNG == TRUE & rows_to_keep == TRUE, pvalues_mindiff_young]), breaks = 100)

#so let us look at these strange with-subunit cases and recompute them from scratch
head(comparison_t2[LR_DETECTED_OLD == TRUE, c("CCI", "means_mindiff_old", "pvalues_mindiff_old")][order(-means_mindiff_old)])
head(comparison_t2[LR_DETECTED_YOUNG == TRUE, c("CCI", "means_mindiff_young", "pvalues_mindiff_young")][order(-means_mindiff_young)])

scratch_cci <- "myeloid.leukocyte.myeloid.leukocyte.Itga4_Itgb1_Plaur"

scratch_seurat1 <- subset(seurat_t2, subset = cell_ontology_class == "myeloid leukocyte" & age_group == "YOUNG")
scratch_seurat2 <- subset(seurat_t2, subset = cell_ontology_class == "myeloid leukocyte" & age_group == "YOUNG")
scratch_LR_score <- (mean(expm1(scratch_seurat2$RNA@data["Plaur",]))  +
                       pmin(mean(expm1(scratch_seurat1$RNA@data["Itga4",])), mean(expm1(scratch_seurat1$RNA@data["Itgb1",]))))/2
#we can see that the scratch result corresponds to scDiffCom but not CELLPHONEDB
scratch_LR_score 
comparison_t2[CCI == scratch_cci, "LR_SCORE_YOUNG"]
comparison_t2[CCI == scratch_cci, "means_young"]
#strangely, the CELLPHONEDB result corresponds to the case when we take the max of the two receptors, or only one of them
#it is like one the gene is not taken into account. 
scratch_wrong_score <- (mean(expm1(scratch_seurat2$RNA@data["Plaur",]))  +
                       pmax(mean(expm1(scratch_seurat1$RNA@data["Itga4",])), mean(expm1(scratch_seurat1$RNA@data["Itgb1",]))))/2
scratch_wrong_score

#We can conclude that there is a problem in the CELLPHONEDB results, but that the scDiffCom results are consistent
# We can observe a good matching between the two methods for simple interactions.
