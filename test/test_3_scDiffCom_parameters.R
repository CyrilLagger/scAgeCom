####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - July 2020
##
## Check the behaviour of scDiffCom in functions of
## its parameters.
##
####################################################
##

## Libraries ####
library(Seurat)
library(scDiffCom)
library(data.table)
library(ggplot2)
library(e1071)
library(VennDiagram)

## Data path ####
dir_data <- "../data_scAgeCom/"

## Load a Seurat objects ####

#here the file corresponds to the Liver tissue from TMS FACS data.
seurat_test_3 <- readRDS(paste0(dir_data, "data_seurat_example.rds"))
seurat_test_3$age_group <- ifelse(seurat_test_3$age %in% c('1m', '3m'), 'young', 'old' )
seurat_test_3$cell_ontology_class <- as.character(seurat_test_3$cell_ontology_class)

## Load the LR database from scDiffCom ####

#we only take a subset of the pairs for testing
LR_test_3 <- LRall[ scsr == TRUE | cpdb == TRUE, c("GENESYMB_L", "GENESYMB_R", "SYMB_LR")]

## Normalization choice ####
seurat_test_3 <- NormalizeData(seurat_test_3, assay = "RNA")
seurat_test_3 <- SCTransform(
  object = seurat_test_3,
  return.only.var.genes = FALSE
)

## Create a list of parameters to test ####

#Note: We want to look at normalization and log-scale.
#      We also want to check the usage of one-sided vs two-sided test for differential p-values.

param_list <- list(
  R_log_one = list(id = "R_log_one", assay = "RNA", log_scale = TRUE, one_sided = TRUE),
  R_log_two = list(id = "R_log_two",assay = "RNA", log_scale = TRUE, one_sided = FALSE),
  R_nlog_one = list(id = "R_nlog_one", assay = "RNA", log_scale = FALSE, one_sided = TRUE),
  R_nlog_two = list(id = "R_nlog_two", assay = "RNA", log_scale = FALSE, one_sided = FALSE),
  S_log_one = list(id = "S_log_one", assay = "SCT", log_scale = TRUE, one_sided = TRUE),
  S_log_two = list(id = "S_log_two", assay = "SCT", log_scale = TRUE, one_sided = FALSE),
  S_nlog_one = list(id = "S_nlog_one", assay = "SCT", log_scale = FALSE, one_sided = TRUE),
  S_nlog_two = list(id = "S_nlog_two", assay = "SCT", log_scale = FALSE, one_sided = FALSE)
)

## Run scDiffCom on all parameters ####

#create a data.table with all conditions
#it will take some times, so we use 1000 iterations as a start
# diffcom_bind_test_3 <- rbindlist(
#   lapply(
#     param_list,
#     function(
#       param
#     ) {
#       run_diffcom(
#         seurat_obj = seurat_test_3,
#         LR_data = LR_test_3,
#         seurat_cell_type_id = "cell_ontology_class",
#         condition_id = "age_group",
#         assay = param[["assay"]],
#         slot = "data",
#         log_scale = param[["log_scale"]],
#         min_cells = 5,
#         threshold = 0.1,
#         iterations = 1000,
#         one_sided = param[["one_sided"]],
#         permutation_analysis = TRUE,
#         return_distr = FALSE
#       )
#     }
#   ),
#   idcol = "param_group"
# )

## Save or load scDiffcom results ####

#saveRDS(diffcom_bind_test_3, file = paste0(dir_data, "test/test_3_data_scDiffCom_allparameters_1000iter.rds"))
diffcom_bind_test_3 <- readRDS(paste0(dir_data, "test/test_3_data_scDiffCom_allparameters_1000iter.rds"))

## Do some processing ####

#scTransform do some internal filtering and remove some genes
#we only keep comparable CCIs
diffcom_param <- diffcom_bind_test_3[LR_GENES %in% diffcom_bind_test_3[param_group == "S_log_one", LR_GENES]]
table(diffcom_param$param_group)

#as for the main downstream analysis we add some useful columns
diffcom_param[, TISSUE := "Liver"]
diffcom_param[, LR_LOGFC := ifelse(
  param_group %in% c("R_log_one", "R_log_two", "S_log_one", "S_log_two"),
  LR_SCORE_old - LR_SCORE_young,
  log(LR_SCORE_old/LR_SCORE_young)
)]
diffcom_param[, LR_LOGFC_ABS := abs(LR_LOGFC)]
diffcom_param[, LR_CELLTYPES := paste(L_CELLTYPE, R_CELLTYPE, sep = "_")]

## Compare detection rate in functions of normalization ####

#Note: we look at detection of CCIs before filtering, so only based on 10% expression and min of 5 cells.

#build the drate data.tables
comp_drate_old <- dcast(
  diffcom_param[, c("param_group", "LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "LR_DETECTED_old")],
  formula = LR_GENES + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_DETECTED_old"
)
comp_drate_young <- dcast(
  diffcom_param[, c("param_group", "LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "LR_DETECTED_young")],
  formula = LR_GENES + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_DETECTED_young"
)

#we observe that scTransform detects more CCIs. It seems like size-factor is conservative here.
table(comp_drate_old$R_log_one, comp_drate_old$S_log_one)
table(comp_drate_young$R_log_one, comp_drate_young$S_log_one)

## Compare LR scores ####

#build score data.tables
comp_LR_score_old <- dcast(
  diffcom_param[, c("param_group", "LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "LR_SCORE_old")],
  formula = LR_GENES + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_SCORE_old"
)
comp_LR_score_young <- dcast(
  diffcom_param[, c("param_group", "LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "LR_SCORE_young")],
  formula = LR_GENES + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_SCORE_young"
)

g_comp_LR_score_old <- cowplot::plot_grid(
  plotlist = list(
    ggplot(comp_LR_score_old, aes(x = R_log_one, y = R_nlog_one)) +
      geom_point() +
      xlab(expression(xi[SFlog])) +
      ylab(expression(xi[SFnlog])),
    ggplot(comp_LR_score_old, aes(x = R_log_one, y = S_log_one)) +
      geom_point() +
      xlab(expression(xi[SFlog])) +
      ylab(expression(xi[SCTlog])),
    ggplot(comp_LR_score_old, aes(x = R_log_one, y = S_nlog_one)) +
      geom_point() +
      xlab(expression(xi[SFlog])) +
      ylab(expression(xi[SCTnlog])),
    ggplot(comp_LR_score_old, aes(x = R_nlog_one, y = S_log_one)) +
      geom_point() +
      xlab(expression(xi[SFnlog])) +
      ylab(expression(xi[SCTlog])),
    ggplot(comp_LR_score_old, aes(x = R_nlog_one, y = S_nlog_one)) +
      geom_point() +
      xlab(expression(xi[SFnlog])) +
      ylab(expression(xi[SCTnlog])),
    ggplot(comp_LR_score_old, aes(x = S_log_one, y = S_nlog_one)) +
      geom_point() +
      xlab(expression(xi[SCTlog])) +
      ylab(expression(xi[SCTnlog]))
  ),
  ncol = 2,
  align = "v"
)

#take some times to save
#ggsave(paste0(dir_data, "test/test_3_plot_compParam_scores.png"), g_comp_LR_score_old, scale = 2)

#histogram for score distributions only for detected CCI

create_histo_score_detected_old <- function(group) {
  ggplot(diffcom_param[param_group == group & LR_DETECTED_old == TRUE], aes(x = LR_SCORE_old)) +
    geom_histogram(bins = 50) +
    xlab(expression(xi[old])) +
    ylab("Counts") +
    theme(text=element_text(size=20)) #+
    #scale_y_log10()
}
create_histo_score_detected_young <- function(group) {
  ggplot(diffcom_param[param_group == group & LR_DETECTED_young == TRUE], aes(x = LR_SCORE_young)) +
    geom_histogram(bins = 50) +
    xlab(expression(xi[young])) +
    ylab("Counts") +
    theme(text=element_text(size=20)) #+
  #scale_y_log10()
}

g_distr_score_old <- cowplot::plot_grid(
  plotlist = lapply(list("R_log_one", "R_nlog_one", "S_log_one", "S_nlog_one"), create_histo_score_detected_old),
  ncol = 2,
  align = "v",
  labels = c("SF-log", "SF-nlog", "SCT-log", "SCT-nlog")
)
#ggsave(paste0(dir_data, "test/test_3_plot_comp_distr_scores_old.png"), g_distr_score_old, scale = 2)
g_distr_score_young <- cowplot::plot_grid(
  plotlist = lapply(list("R_log_one", "R_nlog_one", "S_log_one", "S_nlog_one"), create_histo_score_detected_young),
  ncol = 2,
  align = "v",
  labels = c("SF-log", "SF-nlog", "SCT-log", "SCT-nlog")
)
#ggsave(paste0(dir_data, "test/test_3_plot_comp_distr_scores_young.png"), g_distr_score_young, scale = 2)

## Load filtering functions and filter CCI ####

#load the filtering functions
source("src/src_1_filtering.R")

#functions to get the 20% high scores
get_quantile_score_detected_old <- function(group, q) {
  quantile(diffcom_param[param_group == group & LR_DETECTED_old == TRUE, LR_SCORE_old], q)
}
get_quantile_score_detected_young <- function(group, q) {
  quantile(diffcom_param[param_group == group & LR_DETECTED_young == TRUE, LR_SCORE_young], q)
}

diffcom_param_filter <- rbindlist(
  lapply(
    param_list,
    function(
      param
    ) {
      dt = diffcom_param[param_group == param[["id"]]]
      dt[, param_group := NULL]
      thr_score <- get_quantile_score_detected_young(param[["id"]], 0.8)
      filter_CCI(
        dt = dt,
        CUTOFF_SCORE_YOUNG = thr_score,
        CUTOFF_SCORE_OLD = thr_score,
        CUTOFF_LOGFC = log(1.1)
      )
      reassign_CCI(
        dt = dt,
        CUTOFF_SCORE_YOUNG = thr_score,
        CUTOFF_SCORE_OLD = thr_score,
        CUTOFF_LOGFC = log(1.1)
      )
      return(dt)
    }
  ),
  idcol = "param_group"
)

diffcom_param_filter[, SIMPLE_TYPE := ifelse(
  SIG_TYPE %in% c("FFF"),
  "ND",
  ifelse(
    SIG_TYPE %in% c("TTF"),
    "FLAT",
    ifelse(
      SIG_TYPE %in% c("FTT", "TTTU"),
      "UP",
      "DOWN"
    )
  )
)]

## Compare the CCI classification ####

#build a comparison data.table
comp_class <- dcast(
  diffcom_param_filter[, c("param_group", "LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "SIG_TYPE")],
  formula = LR_GENES + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "SIG_TYPE"
)
comp_topclass <- dcast(
  diffcom_param_filter[, c("param_group", "LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "SIMPLE_TYPE")],
  formula = LR_GENES + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "SIMPLE_TYPE"
)

#actual comparisons
comp_class_summary <- comp_class[, .N, by = names(param_list)]
comp_topclass_summary <- comp_topclass[, .N, by = names(param_list)]

comp_class_normalization <- comp_class[, .N, by = c("R_log_two", "S_log_two")]
comp_topclass_normalization <- comp_topclass[, .N, by = c("R_log_two", "S_log_two")]

comp_class_log2 <- comp_class[, .N, by = c("R_log_two", "R_nlog_two")]
comp_topclass_log2 <- comp_topclass[, .N, by = c("R_log_two", "R_nlog_two")]

#Note: we can observe some differences in the classification depending on the processing methods.
#      There is a priori no clear way to define a better approach based only on this data.
#      There is also probably an effect from choosing the cutoffs in each case. 

## Run scDiffCom to return distributions on all parameters ####

#Note: returning all the distributions requires quite a lot of memory
#      If we do 1000 iterations for the 8 parameter choices, we get 5.4 GB.
#      We don't save the full results, as it takes only 16minutes to get it.
diffcom_distr_bind_test_3 <- lapply(
  param_list,
  function(
    param
  ) {
    run_diffcom(
      seurat_obj = seurat_test_3,
      LR_data = LR_test_3,
      seurat_cell_type_id = "cell_ontology_class",
      condition_id = "age_group",
      assay = param[["assay"]],
      slot = "data",
      log_scale = param[["log_scale"]],
      min_cells = 5,
      threshold = 0.1,
      iterations = 1000,
      one_sided = param[["one_sided"]],
      permutation_analysis = TRUE,
      return_distr = TRUE
    )
  }
)
diffcom_distr_bind_test_3 <- list(
  R_log = diffcom_distr_bind_test_3$R_log_two,
  R_nlog = diffcom_distr_bind_test_3$R_nlog_two,
  S_log = diffcom_distr_bind_test_3$S_log_two,
  S_nlog = diffcom_distr_bind_test_3$S_nlog_two
)

## Compute distribution properties
distr_properties <- lapply(diffcom_distr_bind_test_3, function(distr_param) {
  lapply(distr_param, function(mat) {
    res <- apply(mat, MARGIN = 1, function(mat_row) {
      m <- mean(mat_row)
      md <- median(mat_row)
      sd <- sd(mat_row)
      #shapiro <- shapiro.test(mat_row)$p.value
      skewn <- skewness(mat_row)
      return(c("mean" = m, "median" = md, "sd" = sd, "skewness" = skewn))
    })
    as.data.table(t(res))
  })
})

#create a data.table with all properties
distr_full_table <- rbindlist(
  l = lapply(
    distr_properties,
    function(x) {
      rbindlist(
        l = x,
        idcol = "distribution_type"
      )
    }
  ),
  idcol = "param_group"
)
#saveRDS(distr_full_table, paste0(dir_data, "test/test_3_data_scDiffCom_distribution_table.rds"))

## Look at mean/sd ratio ####

#Note 1: we expect it to be close to zero for distr_diff
#Note 2: there is a problem somewhere probably because of the low cell-number cell-types
param_list_short <- list(
  R_log = param_list$R_log_two,
  R_nlog = param_list$R_nlog_two,
  S_log = param_list$S_log_two,
  S_nlog = param_list$S_nlog_two
)

distr_full_table[, mean_sd := mean/sd]

g_mean_sd <- cowplot::plot_grid(
  plotlist = lapply(names(param_list_short), function(param) {
    ggplot(distr_full_table[param_group == param & distribution_type == "distr_diff" & mean_sd > -0.1], aes(x = mean_sd)) +
      geom_histogram(bins = 100) +
      xlab(expression(mu/sigma)) +
      ylab("Counts") +
      theme(text=element_text(size=20))
  }),
  ncol = 2,
  align = "v",
  labels = c("R_log", "R_nlog", "S_log", "S_nlog")
)
#ggsave(paste0(dir_data, "test/test_3_plot_comp_mean_sd.png"), g_mean_sd, scale = 2)

# g_mean_sd_cond1 <- cowplot::plot_grid(
#   plotlist = lapply(names(param_list_short), function(param) {
#     ggplot(distr_full_table[param_group == param & distribution_type == "distr_cond1"], aes(x = mean_sd)) +
#       geom_histogram(bins = 100) +
#       xlab(expression(mu/sigma)) +
#       ylab("Counts") +
#       theme(text=element_text(size=20))
#   }),
#   ncol = 2,
#   align = "v",
#   labels = c("R_log", "R_nlog", "S_log", "S_nlog")
# )

## Look at skewness ####

g_skewness <- cowplot::plot_grid(
  plotlist = lapply(names(param_list_short), function(param) {
    ggplot(distr_full_table[param_group == param & distribution_type == "distr_diff"], aes(x = skewness)) +
      geom_histogram(bins = 100) +
      xlab("skewness") +
      ylab("Counts") +
      theme(text=element_text(size=20))
  }),
  ncol = 2,
  align = "v",
  labels = c("R_log", "R_nlog", "S_log", "S_nlog")
)
ggsave(paste0(dir_data, "test/test_3_plot_comp_skewness.png"), g_skewness, scale = 2)

#Other ideas: qqplot, shapiro
# see https://stackoverflow.com/questions/7781798/seeing-if-data-is-normally-distributed-in-r/7788452#7788452



