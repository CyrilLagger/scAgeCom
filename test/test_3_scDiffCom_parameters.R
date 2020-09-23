####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - September 2020
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
library(UpSetR)

## Data path ####
dir_data_test <- "../data_scAgeCom/test/"

## Load a Seurat objects ####

#here the file corresponds to the Liver tissue from TMS FACS data.
seurat_t3 <- readRDS(paste0(dir_data_test, "data_seurat_example_facs_liver.rds"))
seurat_t3$age_group <- ifelse(seurat_t3$age %in% c('1m', '3m'), 'young', 'old' )
seurat_t3$cell_ontology_class <- as.character(seurat_t3$cell_ontology_class)

## Load interactions ####
LR_t3 <- scDiffCom::LR6db$LR6db_curated

## Normalization choice ####
seurat_t3 <- NormalizeData(seurat_t3, assay = "RNA")
seurat_t3 <- SCTransform(
  object = seurat_t3,
  return.only.var.genes = FALSE
)

## Create a list of parameters to test ####

#Note: We want to look at normalization and log-scale.
#      We also want to check the usage of one-sided vs two-sided test for differential p-values.

param_t3 <- list(
  R_log_one = list(id = "R_log_one", assay = "RNA", log_scale = TRUE, one_sided = TRUE),
  R_log_two = list(id = "R_log_two",assay = "RNA", log_scale = TRUE, one_sided = FALSE),
  R_nlog_one = list(id = "R_nlog_one", assay = "RNA", log_scale = FALSE, one_sided = TRUE),
  R_nlog_two = list(id = "R_nlog_two", assay = "RNA", log_scale = FALSE, one_sided = FALSE),
  S_log_one = list(id = "S_log_one", assay = "SCT", log_scale = TRUE, one_sided = TRUE),
  S_log_two = list(id = "S_log_two", assay = "SCT", log_scale = TRUE, one_sided = FALSE),
  S_nlog_one = list(id = "S_nlog_one", assay = "SCT", log_scale = FALSE, one_sided = TRUE),
  S_nlog_two = list(id = "S_nlog_two", assay = "SCT", log_scale = FALSE, one_sided = FALSE)
)

param_t3_short <- list(
  R_log = param_t3$R_log_two,
  R_nlog = param_t3$R_nlog_two,
  S_log = param_t3$S_log_two,
  S_nlog = param_t3$S_nlog_two
)

## Run scDiffCom on all parameters ####

#create a data.table with all conditions
#it will take some times, so we use 1000 iterations as a start
# diffcom_bind_t3 <- rbindlist(
#   lapply(
#     param_t3,
#     function(
#       param
#     ) {
#       run_diffcom(
#         seurat_obj = seurat_t3,
#         LR_data = LR_t3,
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


## Preprocessing, filtering and saving results ####
#not needed if already saved

#only keep detected interactions in either young or old samples
#diffcom_bind_t3 <- diffcom_bind_t3[LR_DETECTED_young == TRUE | LR_DETECTED_old == TRUE]
#diffcom_bind_t3[, TISSUE := "Liver"]

#apply filtering analysis to get the different regulation cases
# source("src/src_1_filtering.R")
# diffcom_bind_t3 <- rbindlist(
#   l = lapply(
#     unique(diffcom_bind_t3$param_group),
#     function(param) {
#       temp <- diffcom_bind_t3[param_group == param]
#       temp <- analyze_CCI_per_tissue(
#         data = temp,
#         cutoff_quantile = 0.25,
#         cutoff_logFC_abs = log(1.1),
#         is_log = param_t3[[param]]$log_scale
#       )
#     }
#   ),
#   use.names = TRUE
# )

#save results
#saveRDS(diffcom_bind_t3, file = paste0(dir_data_test, "t3_data_scDiffCom_allparameters_1000iter.rds"))

## Load previously saved results ####

#read files
diffcom_bind_t3 <- readRDS(paste0(dir_data_test, "t3_data_scDiffCom_allparameters_1000iter.rds"))

## Compare classification for each parameter ####

diffcom_bind_t3[ , CASE_TYPE_2 := ifelse(CASE_TYPE %in% c("FTTU", "TTTU"),
                                         "UP",
                                         ifelse(CASE_TYPE %in% c("TFTD", "TTTD"),
                                                "DOWN",
                                                ifelse(CASE_TYPE == "FFF",
                                                       "FFF",
                                                       "FLAT")))]
diffcom_bind_t3[, CCI := paste(LR_CELLTYPE, LR_NAME, sep = "_")]

#look how the CCIs are distributed
ftable(diffcom_bind_t3$param_group, diffcom_bind_t3$CASE_TYPE)
ftable(diffcom_bind_t3[grepl("two", param_group)]$param_group, diffcom_bind_t3[grepl("two", param_group)]$CASE_TYPE)
ftable(diffcom_bind_t3[grepl("two", param_group)]$param_group, diffcom_bind_t3[grepl("two", param_group)]$CASE_TYPE_2)

diffcom_bind_t3[, param_case := paste(param_group, CASE_TYPE_2, sep = "_")]
param_case_list <- unique(diffcom_bind_t3$param_case)
CCI_per_param_list <- sapply(param_case_list, function(i) {
  diffcom_bind_t3[param_case == i]$CCI
})
CCI_per_param_list_2 <- sapply(param_case_list[grepl("two", param_case_list)], function(i) {
  diffcom_bind_t3[param_case == i]$CCI
})
UpSetR::upset(fromList(CCI_per_param_list), nsets = 32, order.by = "freq", nintersects = 35)
UpSetR::upset(fromList(CCI_per_param_list_2), nsets = 16, order.by = "freq", nintersects = 35)

#look and remove the CCI that are consistent over all parameter cases
dcast_param_t3 <- dcast(
  diffcom_bind_t3[, c("CCI", "param_group", "CASE_TYPE_2")],
  formula = CCI ~ param_group,
  value.var = "CASE_TYPE_2"
)
dcast_param_t3[, is_eq := R_log_two == R_nlog_two & R_nlog_two == S_log_two & S_log_two == S_nlog_two]

CCI_conserved <- dcast_param_t3[is_eq == TRUE]$CCI
CCI_non_conserved <- dcast_param_t3[is_eq == FALSE | is.na(is_eq)]$CCI
CCI_non_conserved_noNA <- dcast_param_t3[is_eq == FALSE]$CCI
CCI_per_param_list_2_conserved <- sapply(param_case_list[grepl("two", param_case_list)], function(i) {
  diffcom_bind_t3[param_case == i & CCI %in% CCI_conserved]$CCI
})
CCI_per_param_list_2_non_conserved <- sapply(param_case_list[grepl("two", param_case_list)], function(i) {
  diffcom_bind_t3[param_case == i & CCI %in% CCI_non_conserved]$CCI
})
CCI_per_param_list_2_non_conserved_noNA <- sapply(param_case_list[grepl("two", param_case_list)], function(i) {
  diffcom_bind_t3[param_case == i & CCI %in% CCI_non_conserved_noNA]$CCI
})
CCI_per_param_list_3_conserved <- sapply(param_case_list[grepl("_log_two", param_case_list)], function(i) {
  diffcom_bind_t3[param_case == i & CCI %in% CCI_conserved]$CCI
})


UpSetR::upset(fromList(CCI_per_param_list_2_conserved), nsets = 16, order.by = "freq", nintersects = 35)
UpSetR::upset(fromList(CCI_per_param_list_2_non_conserved), nsets = 16, order.by = "freq", nintersects = 35)
UpSetR::upset(fromList(CCI_per_param_list_2_non_conserved_noNA), nsets = 16, order.by = "freq", nintersects = 35)
UpSetR::upset(fromList(CCI_per_param_list_3_conserved), nsets = 8, order.by = "freq")

#Note: we can observe some differences in the classification depending on the processing methods.
#      There is a priori no clear way to define a better approach based only on this data.

## Compare ORA results for each parameter ####
source("src/src_2_ora.R")

ora_list_t3 <- sapply(
  unique(diffcom_bind_t3$param_group),
  function(param) {
    res <- analyze_ORA(diffcom_bind_t3[param_group == param])
  },
  simplify = FALSE
)

ora_up_LR_NAME <- rbindlist(
  l = lapply(
    ora_list_t3,
    function(i) {
      res <- i[Tissue == "Liver" & Category == "LR_NAME" & pval_UP <= 0.05 & OR_UP >= 1]
      res <- res[, c("Tissue", "Category", "Value", "pval_UP", "OR_UP")]
    }
  ),
  use.names = TRUE,
  idcol = "param"
)
ora_down_LR_NAME <- rbindlist(
  l = lapply(
    ora_list_t3,
    function(i) {
      res <- i[Tissue == "Liver" & Category == "LR_NAME" & pval_DOWN <= 0.05 & OR_DOWN >= 1]
      res <- res[, c("Tissue", "Category", "Value", "pval_DOWN", "OR_DOWN")]
    }
  ),
  use.names = TRUE,
  idcol = "param"
)

ora_up_LR_CELLTYPE <- rbindlist(
  l = lapply(
    ora_list_t3,
    function(i) {
      res <- i[Tissue == "Liver" & Category == "LR_CELLTYPE" & pval_UP <= 0.05 & OR_UP >= 1]
      res <- res[, c("Tissue", "Category", "Value", "pval_UP", "OR_UP")]
    }
  ),
  use.names = TRUE,
  idcol = "param"
)
ora_down_LR_CELLTYPE <- rbindlist(
  l = lapply(
    ora_list_t3,
    function(i) {
      res <- i[Tissue == "Liver" & Category == "LR_CELLTYPE" & pval_DOWN <= 0.05 & OR_DOWN >= 1]
      res <- res[, c("Tissue", "Category", "Value", "pval_DOWN", "OR_DOWN")]
    }
  ),
  use.names = TRUE,
  idcol = "param"
)

param_list <- unique(ora_up_LR_NAME$param)
param_list <- param_list[grepl("two", param_list)]

ora_up_LR_NAME_per_param <- sapply(param_list, function(i) {
  ora_up_LR_NAME[param == i]$Value
})
ora_down_LR_NAME_per_param <- sapply(param_list, function(i) {
  ora_down_LR_NAME[param == i]$Value
})
ora_up_LR_CELLTYPE_per_param <- sapply(param_list, function(i) {
  ora_up_LR_CELLTYPE[param == i]$Value
})
ora_down_LR_CELLTYPE_per_param <- sapply(param_list, function(i) {
  ora_down_LR_CELLTYPE[param == i]$Value
})
UpSetR::upset(fromList(ora_up_LR_NAME_per_param), nsets = 4, order.by = "freq")
UpSetR::upset(fromList(ora_down_LR_NAME_per_param), nsets = 4, order.by = "freq")
UpSetR::upset(fromList(ora_up_LR_CELLTYPE_per_param), nsets = 4, order.by = "freq")
UpSetR::upset(fromList(ora_down_LR_CELLTYPE_per_param), nsets = 4, order.by = "freq")

## Compare LR scores ####

#build score data.tables
comp_LR_score_old <- dcast(
  diffcom_bind_t3[, c("param_group", "LR_SORTED", "L_CELLTYPE", "R_CELLTYPE", "LR_SCORE_old")],
  formula = LR_SORTED + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_SCORE_old"
)
comp_LR_score_young <- dcast(
  diffcom_bind_t3[, c("param_group", "LR_SORTED", "L_CELLTYPE", "R_CELLTYPE", "LR_SCORE_young")],
  formula = LR_SORTED + L_CELLTYPE + R_CELLTYPE ~ param_group,
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
g_comp_LR_score_old
ggsave(paste0(dir_data_test, "t3_plot_compParam_scores.png"), g_comp_LR_score_old, scale = 2)

#histogram for score distributions only for detected CCI 

create_histo_score_detected_old <- function(group) {
  ggplot(diffcom_bind_t3[param_group == group & LR_DETECTED_old == TRUE], aes(x = LR_SCORE_old)) +
    geom_histogram(bins = 50) +
    xlab(expression(xi[old])) +
    ylab("Counts") +
    theme(text=element_text(size=20)) #+
    #scale_y_log10()
}
create_histo_score_detected_young <- function(group) {
  ggplot(diffcom_bind_t3[param_group == group & LR_DETECTED_young == TRUE], aes(x = LR_SCORE_young)) +
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
g_distr_score_old 
ggsave(paste0(dir_data_test, "t3_plot_comp_distr_scores_old.png"), g_distr_score_old, scale = 2)
g_distr_score_young <- cowplot::plot_grid(
  plotlist = lapply(list("R_log_one", "R_nlog_one", "S_log_one", "S_nlog_one"), create_histo_score_detected_young),
  ncol = 2,
  align = "v",
  labels = c("SF-log", "SF-nlog", "SCT-log", "SCT-nlog")
)
g_distr_score_young
ggsave(paste0(dir_data_test, "t3_plot_comp_distr_scores_young.png"), g_distr_score_young, scale = 2)

## Run scDiffCom to return distributions on all parameters ####

#Note: returning all the distributions requires quite a lot of memory
#      If we do 1000 iterations for the 8 parameter choices, we get 5.4 GB.
#      We don't save the full results, as it takes only 16 minutes to get it. 
#      But we will save other downstream results later on.
diffcom_distr_bind_t3 <- lapply(
  param_t3_short,
  function(
    param
  ) {
    run_diffcom(
      seurat_obj = seurat_t3,
      LR_data = LR_t3,
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

## Compute distribution properties
distr_properties <- lapply(diffcom_distr_bind_t3, function(distr_param) {
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
distr_t3_full_table <- rbindlist(
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
saveRDS(distr_t3_full_table, paste0(dir_data_test, "t3_data_scDiffCom_distribution_table.rds"))

distr_t3_full_table <- readRDS(paste0(dir_data_test, "t3_data_scDiffCom_distribution_table.rds"))

## Look at mean/sd ratio ####

#Note 1: we expect it to be close to zero for distr_diff
#Note 2: We need to correct for the cases when there are zero cell-types (it is not a problem in itsel as
#        these cases are treated separately or filtered out in our main analysis.)

distr_t3_full_table[, mean_sd := mean/sd]

g_mean_sd <- cowplot::plot_grid(
  plotlist = lapply(names(param_t3_short), function(param) {
    ggplot(distr_t3_full_table[param_group == param & distribution_type == "distr_diff" & mean_sd > -0.1], aes(x = mean_sd)) +
      geom_histogram(bins = 100) +
      xlab(expression(mu/sigma)) +
      ylab("Counts") +
      theme(text=element_text(size=20))
  }),
  ncol = 2,
  align = "v",
  labels = c("R_log", "R_nlog", "S_log", "S_nlog")
)
g_mean_sd 
ggsave(paste0(dir_data_test, "t3_plot_comp_mean_sd.png"), g_mean_sd, scale = 2)


## Look at the skewness ####

g_skewness <- cowplot::plot_grid(
  plotlist = lapply(names(param_t3_short), function(param) {
    ggplot(distr_t3_full_table[param_group == param & distribution_type == "distr_diff"], aes(x = skewness)) +
      geom_histogram(bins = 100) +
      xlab("skewness") +
      ylab("Counts") +
      theme(text=element_text(size=20))
  }),
  ncol = 2,
  align = "v",
  labels = c("R_log", "R_nlog", "S_log", "S_nlog")
)
g_skewness
ggsave(paste0(dir_data_test, "t3_plot_comp_skewness.png"), g_skewness, scale = 2)

ggplot(distr_t3_full_table[distribution_type == "distr_diff" & param_group %in% c("S_log", "S_nlog")], aes(x = skewness, color = param_group)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.4) +
  xlab("skewness") +
  ylab("Counts") +
  theme(text=element_text(size=20))
ggplot(distr_t3_full_table[distribution_type == "distr_diff" & param_group %in% c("S_log", "R_log")], aes(x = skewness, color = param_group)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.4) +
  xlab("skewness") +
  ylab("Counts") +
  theme(text=element_text(size=20))

#note: we cleary see that the log-normalized cases are less skewed. Therefore we will work only with log-normalized data
# in our main analysis. 
# We also decide to work with size-factor normalization because there is no clear argument to choose between this one and scTransform
# and in this way the data are normalized as in the original TMS paper. Moreover, it is not clear is scTransform
# can be applied to non-UMI data (so to the FACS dataset)



