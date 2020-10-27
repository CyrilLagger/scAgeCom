####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - October 2020
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
library(future)

## Data path ####
dir_data_test <- "../data_scAgeCom/test/"

## Load Seurat objects ####

seurat_objects_t3 <- list(
  facs_liver = readRDS(paste0(dir_data_test, "seurat_testing_tms_facs_liver.rds")),
  facs_spleen = readRDS(paste0(dir_data_test, "seurat_testing_tms_facs_spleen.rds")),
  droplet_spleen = readRDS(paste0(dir_data_test, "seurat_testing_tms_droplet_spleen.rds")),
  calico_spleen = readRDS(paste0(dir_data_test, "seurat_testing_calico_spleen.rds"))
)
seurat_objects_t3 <- lapply(
  seurat_objects_t3,
  function(obj) {
    obj$age_group <- ifelse(obj$age %in% c('1m', '3m', "young"), 'YOUNG', 'OLD' )
    return(obj)
  }
)
seurat_objects_t3 <- lapply(
  seurat_objects_t3,
  NormalizeData,
  normalization.method = "LogNormalize",
  assay = "RNA",
  scale.factor = 10000
)

## Load interactions ####
LR_t3 <- scDiffCom::LR6db$LR6db_curated

## Create a list of parameters to test ####

#Note: We want to look at normalization and log-scale.

param_t3 <- list(
  SF_log = list(id = "SF_log", assay = "RNA", log_scale = TRUE),
  SF_nlog = list(id = "SF_nlog", assay = "RNA", log_scale = FALSE),
  SCT_log = list(id = "SCT_log", assay = "SCT", log_scale = TRUE),
  SCT_nlog = list(id = "SCT_nlog", assay = "SCT", log_scale = FALSE)
)

## Run scDiffCom on all parameters ####

#Note: We also return the distribution, but this requires quite a lot of memory
#      So we perform it one seurat object at a time

seurat_object_temp <- seurat_objects_t3$droplet_spleen
seurat_object_temp <- SCTransform(seurat_object_temp, return.only.var.genes = TRUE)

plan(sequential)
# scdiffcom_t3_temp <- lapply(
#   param_t3,
#   function(
#     param
#   ) {
#     run_scdiffcom(
#       seurat_object = seurat_object_temp,
#       LR_object = LR_t3,
#       celltype_col_id = "cell_ontology_scdiffcom",
#       condition_col_id = "age_group",
#       cond1_name = "YOUNG",
#       cond2_name = "OLD",
#       assay = param[["assay"]],
#       slot = "data",
#       log_scale = param[["log_scale"]],
#       min_cells = 5,
#       pct_threshold = 0.1,
#       permutation_analysis = TRUE,
#       iterations = 1000,
#       cutoff_quantile_score = 0.25,
#       cutoff_pval_specificity = 0.05,
#       cutoff_pval_de = 0.05,
#       cutoff_logfc = log(1.1),
#       return_distr = TRUE,
#       seed = 42,
#       verbose = TRUE,
#       sparse = TRUE
#     )
#   }
# )

saveRDS(scdiffcom_t3_temp, file = paste0(dir_data_test, "t3_data_scDiffCom_allparameters_droplet_spleen.rds"))

## Reading the results ####

scdiffcom_t3 <-  list(
  facs_liver = readRDS(paste0(dir_data_test, "t3_data_scDiffCom_allparameters_facs_liver.rds")),
  facs_spleen = readRDS(paste0(dir_data_test, "t3_data_scDiffCom_allparameters_facs_spleen.rds")),
  droplet_spleen = readRDS(paste0(dir_data_test, "t3_data_scDiffCom_allparameters_droplet_spleen.rds")),
  calico_spleen = readRDS(paste0(dir_data_test, "t3_data_scDiffCom_allparameters_calico_spleen.rds"))
)

## Compare LR scores ####

scdiffcom_t3_comp <- lapply(
  scdiffcom_t3,
  function(dataset) {
    common_cci <- Reduce(
      intersect,
      lapply(
        dataset,
        function (param) {
          param$scdiffcom_dt_filtered[, CCI := paste(LR_CELLTYPE, LR_NAME)]
          param$scdiffcom_dt_filtered[["CCI"]]
        }
      )
    )
    dt <- rbindlist(
      lapply(
        dataset,
        function(param) {
          param$scdiffcom_dt_filtered[CCI %in% common_cci]
        }
      ),
      use.names = TRUE,
      idcol = "param_group"
    )
    dcast(
      dt[, c("param_group", "CCI", "LR_SCORE_OLD", "LR_SCORE_YOUNG", "LOGFC", "REGULATION")],
      formula = CCI ~ param_group,
      value.var = c("LR_SCORE_OLD", "LR_SCORE_YOUNG", "LOGFC", "REGULATION")
    )
  }
)

scdiffcom_t3_logfc_table <- dcast(
  rbindlist(
    lapply(
      scdiffcom_t3,
      function(dataset) {
        rbindlist(
          lapply(
            dataset,
            function(param) {
              param$scdiffcom_dt_filtered[, CCI := paste(LR_CELLTYPE, LR_NAME)]
              param$scdiffcom_dt_filtered[,c("CCI", "LOGFC")]
            }
          ),
          use.names = TRUE,
          idcol = "param_group"
        )
      }
    ),
    use.names = TRUE,
    idcol = "dataset"
  ),
  formula = CCI + dataset ~ param_group,
  value.var = "LOGFC"
)

ggplot(scdiffcom_t3_logfc_table, aes(x = SCT_log, y = SCT_nlog)) + geom_point() + geom_smooth()
ggplot(scdiffcom_t3_logfc_table, aes(x = SF_nlog, y = SCT_nlog)) + geom_point() + geom_smooth()


scdiffcom_t3_distr_params <- lapply(
  scdiffcom_t3,
  function(dataset) {
    dt <- rbindlist(
      lapply(
        dataset,
        function(param) {
          distr <- param$scdiffcom_distributions$distr_diff
          mat <- apply(
            distr,
            MARGIN = 1,
            function(mat_row) {
              m <- mean(mat_row)
              md <- median(mat_row)
              sd <- sd(mat_row)
              ratio <- m/sd
              skewn <- skewness(mat_row)
              return(c("mean" = m, "median" = md, "sd" = sd, "ratio" = ratio, "skewness" = skewn))
            }
          )
          as.data.table(t(mat))
        }
      ),
      use.names = TRUE,
      idcol = "param_group"
    )
  }
)

scdiffcom_t3_distr_params <- rbindlist(
  scdiffcom_t3_distr_params,
  use.names = TRUE,
  idcol = "dataset"
)


ggplot(scdiffcom_t3_distr_params, aes(x = ratio)) +
  geom_boxplot(aes(color = param_group)) +
  facet_wrap(~dataset)
ggplot(scdiffcom_t3_distr_params, aes(x = skewness)) +
  geom_boxplot(aes(color = param_group)) +
  facet_wrap(~dataset)


scdiffcom_t3_pval_comp <- rbindlist(
  lapply(
    scdiffcom_t3,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(param) {
            pval_2sided <- rowSums(
              abs(param$scdiffcom_distributions$distr_diff[,1:1000]) >= 
                abs(param$scdiffcom_distributions$distr_diff[,1001])
            ) /1000
            temp1 <- rowSums(
              param$scdiffcom_distributions$distr_diff[,1:1000] >= 
                param$scdiffcom_distributions$distr_diff[,1001]) /1000
            temp2 <- rowSums(
              param$scdiffcom_distributions$distr_diff[,1:1000] <= 
                param$scdiffcom_distributions$distr_diff[,1001]) /1000
            pval_1sided <- pmin(temp1, temp2)
            data.table(
              pval_2sided = pval_2sided,
              pval_1sided = pval_1sided
            )
          }
        ),
        use.names = TRUE,
        idcol = "param_group"
      )
    }
  ),
  use.names = TRUE,
  idcol = "dataset"
)
  

ggplot(scdiffcom_t3_pval_comp, aes(x = pval_2sided, y = 2*pval_1sided)) + 
  geom_point() + 
  geom_vline(xintercept = 0.05) +
  geom_hline(yintercept = 0.05) +
  facet_grid(param_group ~ dataset)

scdiffcom_t3_pval_comp[, sig_2sided := ifelse(pval_2sided <= 0.05, "sig2", "non-sig2")]
scdiffcom_t3_pval_comp[, sig_1sided := ifelse(2*pval_1sided <= 0.05, "sig1", "non-sig1")]

scdiffcom_t3_signif_table <- dcast(
  scdiffcom_t3_pval_comp[,.N, by = c("sig_1sided", "sig_2sided", "dataset", "param_group")],
  formula = dataset + param_group ~ sig_1sided + sig_2sided,
  value.var = "N"
)





##build score data.tables ####
comp_LR_score_OLD <- dcast(
  scdiffcom_t3_dt_filtered[, c("param_group", "LR_SORTED", "L_CELLTYPE", "R_CELLTYPE", "LR_SCORE_OLD")],
  formula = LR_SORTED + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_SCORE_OLD"
)
comp_LR_score_YOUNG <- dcast(
  scdiffcom_t3_dt_filtered[, c("param_group", "LR_SORTED", "L_CELLTYPE", "R_CELLTYPE", "LR_SCORE_YOUNG")],
  formula = LR_SORTED + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LR_SCORE_YOUNG"
)

comp_LR_logfc_abs <- dcast(
  scdiffcom_t3_dt_filtered[, c("param_group", "LR_SORTED", "L_CELLTYPE", "R_CELLTYPE", "LOGFC_ABS")],
  formula = LR_SORTED + L_CELLTYPE + R_CELLTYPE ~ param_group,
  value.var = "LOGFC_ABS"
)

ggplot(comp_LR_logfc_abs, aes(x = SF_log, y = SCT_log)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(comp_LR_logfc_abs, aes(x = SF_nlog, y = SCT_nlog)) +
  geom_point() +
  geom_smooth(method = "lm")

g_comp_LR_logfc_abs <- cowplot::plot_grid(
  plotlist = list(
    ggplot(comp_LR_logfc_abs, aes(x = SF_log, y = SF_nlog)) +
      geom_point() +
      xlab(expression(xi[SFlog])) +
      ylab(expression(xi[SFnlog])) +
      geom_smooth(method = "lm"),
    ggplot(comp_LR_logfc_abs, aes(x = SF_log, y = SCT_log)) +
      geom_point() +
      xlab(expression(xi[SFlog])) +
      ylab(expression(xi[SCTlog]))+
      geom_smooth(method = "lm"),
    ggplot(comp_LR_logfc_abs, aes(x = SF_log, y = SCT_nlog)) +
      geom_point() +
      xlab(expression(xi[SFlog])) +
      ylab(expression(xi[SCTnlog]))+
      geom_smooth(method = "lm"),
    ggplot(comp_LR_logfc_abs, aes(x = SF_nlog, y = SCT_log)) +
      geom_point() +
      xlab(expression(xi[SFnlog])) +
      ylab(expression(xi[SCTlog]))+
      geom_smooth(method = "lm"),
    ggplot(comp_LR_logfc_abs, aes(x = SF_nlog, y = SCT_nlog)) +
      geom_point() +
      xlab(expression(xi[SFnlog])) +
      ylab(expression(xi[SCTnlog]))+
      geom_smooth(method = "lm"),
    ggplot(comp_LR_logfc_abs, aes(x = SCT_log, y = SCT_nlog)) +
      geom_point() +
      xlab(expression(xi[SCTlog])) +
      ylab(expression(xi[SCTnlog]))+
      geom_smooth(method = "lm")
  ),
  ncol = 2,
  align = "v"
)

g_comp_LR_score_OLD <- cowplot::plot_grid(
  plotlist = list(
    ggplot(comp_LR_score_OLD, aes(x = SF_log, y = SF_nlog)) +
      geom_point() +
      xlab(expression(xi[SFlog])) +
      ylab(expression(xi[SFnlog])),
    ggplot(comp_LR_score_OLD, aes(x = SF_log, y = SCT_log)) +
      geom_point() +
      xlab(expression(xi[SFlog])) +
      ylab(expression(xi[SCTlog])),
    ggplot(comp_LR_score_OLD, aes(x = SF_log, y = SCT_nlog)) +
      geom_point() +
      xlab(expression(xi[SFlog])) +
      ylab(expression(xi[SCTnlog])),
    ggplot(comp_LR_score_OLD, aes(x = SF_nlog, y = SCT_log)) +
      geom_point() +
      xlab(expression(xi[SFnlog])) +
      ylab(expression(xi[SCTlog])),
    ggplot(comp_LR_score_OLD, aes(x = SF_nlog, y = SCT_nlog)) +
      geom_point() +
      xlab(expression(xi[SFnlog])) +
      ylab(expression(xi[SCTnlog])),
    ggplot(comp_LR_score_OLD, aes(x = SCT_log, y = SCT_nlog)) +
      geom_point() +
      xlab(expression(xi[SCTlog])) +
      ylab(expression(xi[SCTnlog]))
  ),
  ncol = 2,
  align = "v"
)
g_comp_LR_score_OLD
ggsave(paste0(dir_data_test, "t3_plot_compParam_scores.png"), g_comp_LR_score_OLD, scale = 2)

#histogram for score distributions only for detected CCI 

create_histo_score_detected_OLD <- function(group) {
  ggplot(diffcom_bind_t3[param_group == group & LR_DETECTED_OLD == TRUE], aes(x = LR_SCORE_OLD)) +
    geom_histogram(bins = 50) +
    xlab(expression(xi[OLD])) +
    ylab("Counts") +
    theme(text=element_text(size=20)) #+
  #scale_y_log10()
}
create_histo_score_detected_YOUNG <- function(group) {
  ggplot(diffcom_bind_t3[param_group == group & LR_DETECTED_YOUNG == TRUE], aes(x = LR_SCORE_YOUNG)) +
    geom_histogram(bins = 50) +
    xlab(expression(xi[YOUNG])) +
    ylab("Counts") +
    theme(text=element_text(size=20)) #+
  #scale_y_log10()
}

g_distr_score_OLD <- cowplot::plot_grid(
  plotlist = lapply(list("SF_log", "SF_nlog", "SCT_log", "SCT_nlog"), create_histo_score_detected_OLD),
  ncol = 2,
  align = "v",
  labels = c("SF-log", "SF-nlog", "SCT-log", "SCT-nlog")
)
g_distr_score_OLD 
ggsave(paste0(dir_data_test, "t3_plot_comp_distr_scores_OLD.png"), g_distr_score_OLD, scale = 2)
g_distr_score_YOUNG <- cowplot::plot_grid(
  plotlist = lapply(list("SF_log", "SF_nlog", "SCT_log", "SCT_nlog"), create_histo_score_detected_YOUNG),
  ncol = 2,
  align = "v",
  labels = c("SF-log", "SF-nlog", "SCT-log", "SCT-nlog")
)
g_distr_score_YOUNG
ggsave(paste0(dir_data_test, "t3_plot_comp_distr_scores_YOUNG.png"), g_distr_score_YOUNG, scale = 2)


## Create a copy of the results but with larger logfc threshold ####

scdiffcom_t3_bis <- lapply(
  scdiffcom_t3,
  scDiffCom::run_filtering_and_ORA,
  new_cutoff_logfc = log(1.5)
)

## Create data.tables of raw CCIs, filtered CCIs and ORA for each parameter ####

scdiffcom_t3_dt_raw <- rbindlist(
  lapply(scdiffcom_t3, function(i) i$scdiffcom_dt_raw),
  idcol = "param_group"
)

scdiffcom_t3_dt_filtered <- rbindlist(
  lapply(scdiffcom_t3, function(i) i$scdiffcom_dt_filtered),
  idcol = "param_group"
)

scdiffcom_t3_dt_ORA <- rbindlist(
  lapply(scdiffcom_t3, function(i) i$ORA),
  idcol = "param_group"
)

scdiffcom_t3_dt_raw_bis <- rbindlist(
  lapply(scdiffcom_t3_bis, function(i) i$scdiffcom_dt_raw),
  idcol = "param_group"
)

scdiffcom_t3_dt_filtered_bis <- rbindlist(
  lapply(scdiffcom_t3_bis, function(i) i$scdiffcom_dt_filtered),
  idcol = "param_group"
)

scdiffcom_t3_dt_ORA_bis <- rbindlist(
  lapply(scdiffcom_t3_bis, function(i) i$ORA),
  idcol = "param_group"
)

## Compare classification for each parameter ####

#look how the CCIs are distributed
ftable(scdiffcom_t3_dt_filtered$param_group, scdiffcom_t3_dt_filtered$REGULATION)
ftable(scdiffcom_t3_dt_filtered_bis$param_group, scdiffcom_t3_dt_filtered_bis$REGULATION)
ftable(scdiffcom_t3_dt_filtered$param_group, scdiffcom_t3_dt_filtered$REGULATION_SIMPLE)
ftable(scdiffcom_t3_dt_filtered_bis$param_group, scdiffcom_t3_dt_filtered_bis$REGULATION_SIMPLE)


hist(scdiffcom_t3_dt_filtered[param_group == "SCT_log"]$LOGFC_ABS, breaks = 50)
hist(scdiffcom_t3_dt_filtered[param_group == "SCT_nlog"]$LOGFC_ABS, breaks = 50)

hist(scdiffcom_t3_dt_filtered[param_group == "SF_log"]$LOGFC_ABS, breaks = 50)
hist(scdiffcom_t3_dt_filtered[param_group == "SF_nlog"]$LOGFC_ABS, breaks = 50)


log(2)





scdiffcom_t3_dt_filtered[, CCI := paste(LR_CELLTYPE, LR_NAME, sep = "_")]
scdiffcom_t3_dt_filtered_bis[, CCI := paste(LR_CELLTYPE, LR_NAME, sep = "_")]

scdiffcom_t3_dt_filtered[, param_case := paste(param_group, REGULATION_SIMPLE, sep = "_")]
scdiffcom_t3_dt_filtered_bis[, param_case := paste(param_group, REGULATION_SIMPLE, sep = "_")]

param_case_list <- unique(scdiffcom_t3_dt_filtered$param_case)

CCI_per_param_list <- sapply(param_case_list, function(i) {
  scdiffcom_t3_dt_filtered[param_case == i]$CCI
})
UpSetR::upset(fromList(CCI_per_param_list), nsets = 12, order.by = "freq", nintersects = 35)
CCI_per_param_list_bis <- sapply(param_case_list, function(i) {
  scdiffcom_t3_dt_filtered_bis[param_case == i]$CCI
})
UpSetR::upset(fromList(CCI_per_param_list_bis), nsets = 12, order.by = "freq", nintersects = 35)


#look and remove the CCI that are consistent over all parameter cases
dcast_param_t3 <- dcast(
  scdiffcom_t3_dt_filtered[, c("CCI", "param_group", "REGULATION_SIMPLE")],
  formula = CCI ~ param_group,
  value.var = "REGULATION_SIMPLE"
)
dcast_param_t3[, is_eq := SF_log == SF_nlog & SF_nlog == SCT_log & SCT_log == SCT_nlog]

CCI_conserved <- dcast_param_t3[is_eq == TRUE]$CCI
CCI_non_conserved <- dcast_param_t3[is_eq == FALSE | is.na(is_eq)]$CCI
CCI_non_conserved_noNA <- dcast_param_t3[is_eq == FALSE]$CCI
CCI_per_param_list_2_conserved <- sapply(param_case_list[grepl("two", param_case_list)], function(i) {
  scdiffcom_t3_dt_filtered[param_case == i & CCI %in% CCI_conserved]$CCI
})
CCI_per_param_list_2_non_conserved <- sapply(param_case_list, function(i) {
  scdiffcom_t3_dt_filtered[param_case == i & CCI %in% CCI_non_conserved]$CCI
})
CCI_per_param_list_2_non_conserved_noNA <- sapply(param_case_list, function(i) {
  scdiffcom_t3_dt_filtered[param_case == i & CCI %in% CCI_non_conserved_noNA]$CCI
})
CCI_per_param_list_2_conserved <- sapply(param_case_list, function(i) {
  scdiffcom_t3_dt_filtered[param_case == i & CCI %in% CCI_conserved]$CCI
})


UpSetR::upset(fromList(CCI_per_param_list_2_conserved), nsets = 12, order.by = "freq", nintersects = 35)
UpSetR::upset(fromList(CCI_per_param_list_2_non_conserved), nsets = 12, order.by = "freq", nintersects = 35)
UpSetR::upset(fromList(CCI_per_param_list_2_non_conserved_noNA), nsets = 12, order.by = "freq", nintersects = 35)


#Note: we can observe some differences in the classification depending on the processing methods.
#      It is not clear which is the best method, but most of the time 3 of the 4 groups match together. 
#      Usually it is up (or down) that goes to flat in one of the 4 cases.

## Compare ORA results for each parameter ####


ora_up_LR_NAME <- scdiffcom_t3_dt_ORA[Category == "LR_NAME" & pval_adjusted_UP <= 0.05 & OR_UP >= 1,
                                      c("Category", "Value", "pval_adjusted_UP", "OR_UP", "param_group")]
ora_down_LR_NAME <- scdiffcom_t3_dt_ORA[Category == "LR_NAME" & pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1,
                                      c("Category", "Value", "pval_adjusted_UP", "OR_UP", "param_group")]
ora_up_LR_CELLTYPE <- scdiffcom_t3_dt_ORA[Category == "LR_CELLTYPE" & pval_adjusted_UP <= 0.05 & OR_UP >= 1,
                                      c("Category", "Value", "pval_adjusted_UP", "OR_UP", "param_group")]
ora_down_LR_CELLTYPE <- scdiffcom_t3_dt_ORA[Category == "LR_CELLTYPE" & pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1,
                                        c("Category", "Value", "pval_adjusted_UP", "OR_UP", "param_group")]


ora_up_LR_NAME_per_param <- sapply(param_t3, function(i) {
  ora_up_LR_NAME[param_group %in% i]$Value
})
ora_down_LR_NAME_per_param <- sapply(param_t3, function(i) {
  ora_down_LR_NAME[param_group %in% i]$Value
})
ora_up_LR_CELLTYPE_per_param <- sapply(param_t3, function(i) {
  ora_up_LR_CELLTYPE[param_group %in% i]$Value
})
ora_down_LR_CELLTYPE_per_param <- sapply(param_t3, function(i) {
  ora_down_LR_CELLTYPE[param_group %in% i]$Value
})
UpSetR::upset(fromList(ora_up_LR_NAME_per_param), nsets = 4, order.by = "freq")
UpSetR::upset(fromList(ora_down_LR_NAME_per_param), nsets = 4, order.by = "freq")
UpSetR::upset(fromList(ora_up_LR_CELLTYPE_per_param), nsets = 4, order.by = "freq")
UpSetR::upset(fromList(ora_down_LR_CELLTYPE_per_param), nsets = 4, order.by = "freq")
