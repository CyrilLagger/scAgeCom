####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - October 2020
##
## Double-check that scDiffCom returns correct
## results and distributions.
##
####################################################
##

## Libraries ####
library(Seurat)
library(data.table)
library(scDiffCom)
library(ggplot2)
library(pbapply)
library(future)

## Data path and parameters####
dir_data_test <- "../data_scAgeCom/test/"

n_iter <- 1000 #number of permutations for the double-check

## Load and process Seurat objects for testing ####

seurat_objects_t1 <- list(
  facs_liver = readRDS(paste0(dir_data_test, "seurat_testing_tms_facs_liver.rds")),
  facs_spleen = readRDS(paste0(dir_data_test, "seurat_testing_tms_facs_spleen.rds")),
  droplet_spleen = readRDS(paste0(dir_data_test, "seurat_testing_tms_droplet_spleen.rds")),
  calico_spleen = readRDS(paste0(dir_data_test, "seurat_testing_calico_spleen.rds"))
)
?NormalizeData
seurat_objects_t1 <- lapply(
  seurat_objects_t1,
  NormalizeData,
  normalization.method = "LogNormalize",
  assay = "RNA",
  scale.factor = 10000
)
seurat_objects_t1 <- lapply(
  seurat_objects_t1,
  function(obj) {
    obj$age_group <- ifelse(obj$age %in% c('1m', '3m', "young"), 'YOUNG', 'OLD' )
    return(obj)
  }
)
lapply(
  seurat_objects_t1,
  function(i) {
    table(i$age, i$age_group)
  }
)


## Load the LR database from scDiffCom ####
LR_t1 <- scDiffCom::LR6db$LR6db_curated

## Typical usage of scDiffCom ####

#NOTES: --running scDiffcom takes a few minutes for 1000 iterations (e.g: around 3 min on a laptop for the parameters below)
#       --the command does not need to be run each time if the results have been saved already
#       --uncomment the lines if you want to test it
#       --scdiffcom can be run in parallel through the future package

?run_scdiffcom
future::plan(sequential)
scdiffcom_t1 <- lapply(
  seurat_objects_t1,
  function(i) {
    run_scdiffcom(
      seurat_object = i,
      LR_object = LR_t1,
      celltype_col_id = "cell_ontology_scdiffcom",
      condition_col_id = "age_group",
      cond1_name = "YOUNG",
      cond2_name = "OLD",
      assay = "RNA",
      slot = "data",
      log_scale = TRUE,
      min_cells = 5,
      pct_threshold = 0.1,
      permutation_analysis = TRUE,
      iterations = n_iter,
      cutoff_quantile_score = 0.25,
      cutoff_pval_specificity = 0.05,
      cutoff_pval_de = 0.05,
      cutoff_logfc = log(1.1),
      return_distr = TRUE,
      seed = 42,
      verbose = TRUE,
      sparse = TRUE
    )
  }
) 

## Save or load results ####

#either save
#saveRDS(object = scdiffcom_t1, file = paste0(dir_data_test, "t1_data_scdiffcom.rds"))
#or read
scdiffcom_t1 <- readRDS(paste0(dir_data_test, "t1_data_scDiffcom.rds"))

## Select random CCIs for further comparison ####

#4 random CCIs per dataset (read them first if already saved)

cci_random_filtered_t1 <- readRDS(paste0(dir_data_test, "t1_data_cci_random_filtered.rds"))

# cci_random_filtered_t1 <- lapply(
#   scdiffcom_t1,
#   function(i) {
#     rbindlist(
#       list(
#         i$scdiffcom_dt_filtered[sample(.N, 2)],
#         i$scdiffcom_dt_filtered[grepl("_", LR_NAME)][sample(.N, 2)]
#       )
#     )
#   }
# )
# lapply(
#   cci_random_filtered_t1,
#   function(i) i$REGULATION
# )
# 
# saveRDS(cci_random_filtered_t1, file = paste0(dir_data_test, "t1_data_cci_random_filtered.rds"))

cci_random_filtered_distr_t1 <- lapply(
  1:length(cci_random_filtered_t1),
  function(i) {
    lapply(
      scdiffcom_t1[[i]]$scdiffcom_distributions,
      function(j) {
        j[paste(cci_random_filtered_t1[[i]]$LR_SORTED, cci_random_filtered_t1[[i]]$LR_CELLTYPE, sep = "_"),]
      }
    )
  }
)
names(cci_random_filtered_distr_t1) <- names(cci_random_filtered_t1)

## Select specific CCIs for further comparison ####

cci_specific_filtered_t1 <- readRDS(paste(dir_data_test, "t1_data_cci_specific_filtered.rds"))

# cci_specific_filtered_t1 <- lapply(
#   scdiffcom_t1,
#   function(i) {
#     rbindlist(
#       list(
#         i$scdiffcom_dt_filtered[REGULATION == "UP_APPEARS"][3],
#         i$scdiffcom_dt_filtered[REGULATION == "UP"][10],
#         i$scdiffcom_dt_filtered[REGULATION == "FLAT"][7],
#         i$scdiffcom_dt_filtered[REGULATION == "DOWN_DISAPPEARS"][3],
#         i$scdiffcom_dt_filtered[REGULATION == "DOWN"][14]
#       )
#     )
#   }
# )
# lapply(
#   cci_specific_filtered_t1,
#   function(i) i$REGULATION
# )
# saveRDS(cci_specific_filtered_t1, file = paste(dir_data_test, "t1_data_cci_specific_filtered.rds"))

cci_specific_filtered_distr_t1 <- lapply(
  1:length(cci_specific_filtered_t1),
  function(i) {
    lapply(
      scdiffcom_t1[[i]]$scdiffcom_distributions,
      function(j) {
        j[paste(cci_specific_filtered_t1[[i]]$LR_SORTED, cci_specific_filtered_t1[[i]]$LR_CELLTYPE, sep = "_"),]
      }
    )
  }
)
names(cci_specific_filtered_distr_t1) <- names(cci_specific_filtered_t1)

## Double-check that the LR scores are correct ####

# We create a (slow) function that compute the LR scores in a different way than scDiffCom, such that
# we can double check the results.
compare_cci_score_test <- function(
  cci_test_dt,
  seurat_object,
  condition,
  log_scale
) {
  res <- lapply(1:nrow(cci_test_dt), function(i) {
    cci_test <- cci_test_dt[i,]
    seurat_sub_l <- subset(seurat_object, features = na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])),
                           subset = (cell_ontology_scdiffcom == cci_test[["L_CELLTYPE"]] &
                                       age_group == condition))
    seurat_sub_r <- subset(seurat_object, features = na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")])),
                           subset = (cell_ontology_scdiffcom == cci_test[["R_CELLTYPE"]] &
                                       age_group == condition))
    if(log_scale) {
      LR_avg_expr <- (min(sapply(1:nrow(seurat_sub_l), function(i) {
        mean((seurat_sub_l[["RNA"]]@data[i,]))
      })) +
        min(sapply(1:nrow(seurat_sub_r), function(i) {
          mean((seurat_sub_r[["RNA"]]@data[i,]))
        }))
      )/2
    } else {
      LR_avg_expr <- (min(sapply(1:nrow(seurat_sub_l), function(i) {
        mean(expm1(seurat_sub_l[["RNA"]]@data[i,]))
      })) +
        min(sapply(1:nrow(seurat_sub_r), function(i) {
          mean(expm1(seurat_sub_r[["RNA"]]@data[i,]))
        }))
      )/2
    }
    return(LR_avg_expr - cci_test[[paste0("LR_SCORE_", condition)]])
  }
  )
  return(res)
}


#should return zero or numerically close to zero
mapply(compare_cci_score_test, cci_random_filtered_t1, seurat_objects_t1,
       MoreArgs = list(condition = "YOUNG", log_scale = TRUE))
mapply(compare_cci_score_test, cci_random_filtered_t1, seurat_objects_t1,
       MoreArgs = list(condition = "OLD", log_scale = TRUE))
mapply(compare_cci_score_test, cci_specific_filtered_t1, seurat_objects_t1,
       MoreArgs = list(condition = "YOUNG", log_scale = TRUE))
mapply(compare_cci_score_test, cci_specific_filtered_t1, seurat_objects_t1,
       MoreArgs = list(condition = "OLD", log_scale = TRUE))


## Functions to double check the p-values from the permutation test ####

#We create a (very slow) function that compute the LR differential p-value (and also returns the distribution)
#in a different way than scDiffCom.
compare_cci_pvalue_diff_test <- function(
  seurat_obj,
  cci_test_dt,
  iterations,
  log_scale
) {
  res <- lapply(1:nrow(cci_test_dt), function(i) {
    cci_test <- cci_test_dt[i,]
    cells_use <- colnames(seurat_obj)[seurat_obj$cell_ontology_scdiffcom  %in% c(cci_test[["L_CELLTYPE"]], cci_test[["R_CELLTYPE"]])]
    seurat_temp <- subset(seurat_obj, features = c(na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])), na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")]))),
                          cells = cells_use)
    get_shuffled_LR <- function(permut) {
      if(permut) {
        seurat_temp$age_group[seurat_temp$cell_ontology_scdiffcom == cci_test[["L_CELLTYPE"]] ] <- 
          sample(seurat_temp$age_group[seurat_temp$cell_ontology_scdiffcom == cci_test[["L_CELLTYPE"]] ])
        seurat_temp$age_group[seurat_temp$cell_ontology_scdiffcom == cci_test[["R_CELLTYPE"]] ] <- 
          sample(seurat_temp$age_group[seurat_temp$cell_ontology_scdiffcom == cci_test[["R_CELLTYPE"]] ])
      }
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_scdiffcom %in% cci_test[["L_CELLTYPE"]] &
                                            seurat_temp$age_group == "YOUNG"]
      seurat_sub_l_young <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])),
                                   cells = cells_use1)
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_scdiffcom %in% cci_test[["R_CELLTYPE"]] &
                                            seurat_temp$age_group == "YOUNG"]
      seurat_sub_r_young <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")])),
                                   cells = cells_use1)
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_scdiffcom %in% cci_test[["L_CELLTYPE"]] &
                                            seurat_temp$age_group == "OLD"]
      seurat_sub_l_old <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])),
                                 cells = cells_use1)
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_scdiffcom %in% cci_test[["R_CELLTYPE"]] &
                                            seurat_temp$age_group == "OLD"]
      seurat_sub_r_old <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")])),
                                 cells = cells_use1)
      if(log_scale) {
        LR_diff <- (min(sapply(1:nrow(seurat_sub_l_old), function(i) {
          mean((seurat_sub_l_old[["RNA"]]@data[i,]))
        })) +
          min(sapply(1:nrow(seurat_sub_r_old), function(i) {
            mean((seurat_sub_r_old[["RNA"]]@data[i,]))
          }))
        )/2 -
          (min(sapply(1:nrow(seurat_sub_l_young), function(i) {
            mean((seurat_sub_l_young[["RNA"]]@data[i,]))
          })) +
            min(sapply(1:nrow(seurat_sub_r_young), function(i) {
              mean((seurat_sub_r_young[["RNA"]]@data[i,]))
            }))
          )/2
      } else {
        LR_diff <- (min(sapply(1:nrow(seurat_sub_l_old), function(i) {
          mean(expm1(seurat_sub_l_old[["RNA"]]@data[i,]))
        })) +
          min(sapply(1:nrow(seurat_sub_r_old), function(i) {
            mean(expm1(seurat_sub_r_old[["RNA"]]@data[i,]))
          }))
        )/2 -
          (min(sapply(1:nrow(seurat_sub_l_young), function(i) {
            mean(expm1(seurat_sub_l_young[["RNA"]]@data[i,]))
          })) +
            min(sapply(1:nrow(seurat_sub_r_young), function(i) {
              mean(expm1(seurat_sub_r_young[["RNA"]]@data[i,]))
            }))
          )/2
      }
      return(LR_diff)
    }
    all_LR <- pbreplicate(iterations, get_shuffled_LR(TRUE))
    all_LR <- c(all_LR, get_shuffled_LR(FALSE))
    pval <- sum(abs(all_LR[1:iterations]) >= abs(all_LR[(iterations + 1)])) / iterations
    return(list(pval = pval,
                distr = all_LR))
  }
  )
}

#same for specificity p-value
compare_cci_pvalue_specific_test <- function(
  seurat_obj,
  cci_test_dt,
  iterations,
  condition,
  log_scale
) {
  res <- lapply(1:nrow(cci_test_dt), function(i) {
    cci_test <- cci_test_dt[i,]
    cells_use <- colnames(seurat_obj)[seurat_obj$age_group == condition]
    seurat_temp <- subset(seurat_obj, features = c(na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])), na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")]))),
                          cells = cells_use)
    get_shuffled_LR_ct <- function(permut) {
      if(permut) {
        seurat_temp$cell_ontology_scdiffcom <- sample(as.character(seurat_temp$cell_ontology_scdiffcom))
      }
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_scdiffcom == cci_test[["L_CELLTYPE"]]]
      seurat_sub_l <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])),
                             cells = cells_use1)
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_scdiffcom == cci_test[["R_CELLTYPE"]]]
      seurat_sub_r <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")])),
                             cells = cells_use1)
      if(log_scale) {
        score <- (mean((seurat_sub_l[["RNA"]]@data[1,])) + mean((seurat_sub_r[["RNA"]]@data[1,])))/2
      } else {
        score <- (mean(expm1(seurat_sub_l[["RNA"]]@data[1,])) + mean(expm1(seurat_sub_r[["RNA"]]@data[1,])))/2
      }
      return(score)
    }
    all_LR <- pbreplicate(iterations, get_shuffled_LR_ct(TRUE))
    all_LR <- c(all_LR, get_shuffled_LR_ct(FALSE))
    pval <- sum(all_LR[1:iterations] >= all_LR[(iterations + 1)]) / iterations
    return(list(pval = pval,
                distr = all_LR))
  }
  )
}

#We create a function that can display two distributions on the same plot
comp_distr <- function(
  distr1,
  distr2,
  case,
  name
) {
  if(case == "young") {
    expr <- expression(xi[young])
  } else if(case == "old") {
    expr <- expression(xi[old])
  } else {
    expr <- expression(xi[old]-xi[young])
  }
  ggplot(data.frame(x = c(distr1, distr2),
                    Case = c(rep("scDiffCom", length(distr1)), rep("double-check", length(distr2)))),
         aes(x=x, fill = Case)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
    geom_vline(xintercept = distr1[length(distr1)] ) +
    xlab(expr) +
    ylab("Counts") +
    theme(text=element_text(size=20)) +
    ggtitle(name) +
    theme(plot.title = element_text(size = 10))
}

## Need to be redone below ##Computation of the p-values ####

# note: the function for the double check are very slow (as used only here as a double check), so we run it only for a few cases

pval_diff_1 <- compare_cci_pvalue_diff_test(
  seurat_objects_t1$facs_liver,
  cci_random_filtered_t1$facs_liver[3],
  iterations = 10,
  log_scale = TRUE
)
pval_diff_2 <- compare_cci_pvalue_diff_test(
  seurat_objects_t1$facs_spleen,
  cci_random_filtered_t1$facs_spleen[2],
  iterations = 10,
  log_scale = TRUE
)
pval_diff_3 <- compare_cci_pvalue_diff_test(
  seurat_objects_t1$droplet_spleen,
  cci_random_filtered_t1$droplet_spleen[1],
  iterations = 10,
  log_scale = TRUE
)
pval_diff_4 <- compare_cci_pvalue_diff_test(
  seurat_objects_t1$calico_spleen,
  cci_random_filtered_t1$calico_spleen[4],
  iterations = 10,
  log_scale = TRUE
)

#compute the p-values and distributions
#pv_diff_random_t1 <- compare_cci_pvalue_diff_test(seurat_t1, cci_random_filtered_t1, iterations = n_iter, log_scale = TRUE)
#pv_diff_specific_t1 <- compare_cci_pvalue_diff_test(seurat_t1, cci_specific_filtered_t1, iterations = n_iter, log_scale = TRUE)

#save them or read them
#saveRDS(object = pv_diff_random_t1, file = paste0(dir_data_test, "t1_data_pv_diff_random.rds"))
#saveRDS(object = pv_diff_specific_t1, file = paste0(dir_data_test, "t1_data_pv_diff_specific.rds"))

pv_diff_random_t1 <- readRDS(paste0(dir_data_test, "t1_data_pv_diff_random.rds"))
pv_diff_specific_t1 <- readRDS(paste0(dir_data_test, "t1_data_pv_diff_specific.rds"))

#compare the p-values 
lapply(seq_along(pv_diff_random_t1), function(i) {
  c(cci_random_filtered_t1[i,]$PVAL_DIFF, pv_diff_random_t1[[i]]$pval)
})
lapply(seq_along(pv_diff_specific_t1), function(i) {
  c(cci_specific_filtered_t1[i,]$PVAL_DIFF, pv_diff_specific_t1[[i]]$pval)
})

#compare the distributions
t1_plot_random_pval_diff_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_diff_random_t1), function(i) {
    comp_distr(cci_random_filtered_distr_t1$distr_diff[i,], pv_diff_random_t1[[i]]$distr, case = "diff",
               name = paste(cci_random_filtered_t1$LR_CELLTYPE, cci_random_filtered_t1$LR_NAME, sep = "_")[[i]])
  }),
  ncol = 2,
  align = "v"
)
t1_plot_random_pval_diff_comp
ggsave(filename = paste0(dir_data_test, "t1_plot_random_distr_permutations_comp.png"), plot = t1_plot_random_pval_diff_comp, scale = 2)
t1_plot_specific_pval_diff_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_diff_specific_t1), function(i) {
    comp_distr(cci_specific_filtered_distr_t1$distr_diff[i,], pv_diff_specific_t1[[i]]$distr, case = "diff",
               name = paste(cci_specific_filtered_t1$LR_CELLTYPE, cci_specific_filtered_t1$LR_NAME, sep = "_")[[i]])
  }),
  ncol = 2,
  align = "v"
)
t1_plot_specific_pval_diff_comp
ggsave(filename = paste0(dir_data_test, "t1_plot_specific_distr_permutations_comp.png"), plot = t1_plot_specific_pval_diff_comp, scale = 2)

## Double-check that the LR specificity p-values are correct ####

#We create a (very slow) function that compute the LR specificity p-value (and also returns the distribution)
#in a different way than scDiffCom.

#compute the p-values and distributions
pv_spec_young_random_t1 <- compare_cci_pvalue_specific_test(seurat_t1, cci_random_filtered_t1, iterations = n_iter, condition = "YOUNG", log_scale = TRUE)
pv_spec_old_random_t1 <- compare_cci_pvalue_specific_test(seurat_t1, cci_random_filtered_t1, iterations = n_iter, condition = "OLD", log_scale = TRUE)
pv_spec_young_specific_t1 <- compare_cci_pvalue_specific_test(seurat_t1, cci_specific_filtered_t1, iterations = n_iter, condition = "YOUNG", log_scale = TRUE)
pv_spec_old_specific_t1 <- compare_cci_pvalue_specific_test(seurat_t1, cci_specific_filtered_t1, iterations = n_iter, condition = "OLD", log_scale = TRUE)

#save them!!!

#compare the p-values 
lapply(seq_along(pv_spec_young_random_t1), function(i) {
  c(cci_random_filtered_t1[i,]$PVAL_YOUNG, pv_spec_young_random_t1[[i]]$pval)
})
lapply(seq_along(pv_spec_old_random_t1), function(i) {
  c(cci_random_t1[i,]$PVAL_OLD, pv_spec_old_random_t1[[i]]$pval)
})
lapply(seq_along(pv_spec_young_specific_t1), function(i) {
  c(cci_specific_t1[i,]$PVAL_YOUNG, pv_spec_young_specific_t1[[i]]$pval)
})
lapply(seq_along(pv_spec_old_specific_t1), function(i) {
  c(cci_specific_t1[i,]$PVAL_OLD, pv_spec_old_specific_t1[[i]]$pval)
})

#compare the distributions
t1_plot_random_pval_spec_young_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_spec_young_random_t1), function(i) {
    comp_distr(cci_random_filtered_distr_t1$distr_YOUNG[i,], pv_spec_young_random_t1[[i]]$distr, case = "young")
  }),
  ncol = 2,
  align = "v"
)
ggsave(filename = paste0(dir_data_test, "t1_plot_random_pval_spec_young_comp.png"), plot = t1_plot_random_pval_spec_young_comp, scale = 2)
t1_plot_random_pval_spec_old_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_spec_old_random_t1), function(i) {
    comp_distr(cci_random_distr_t1$distr_cond1[i,], pv_spec_old_random_t1[[i]]$distr, case = "old")
  }),
  ncol = 2,
  align = "v"
)
ggsave(filename = paste0(dir_data_test, "t1_plot_random_pval_spec_old_comp.png"), plot = t1_plot_random_pval_spec_old_comp, scale = 2)
t1_plot_specific_pval_spec_young_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_spec_young_specific_t1), function(i) {
    comp_distr(cci_specific_distr_t1$distr_cond2[i,], pv_spec_young_specific_t1[[i]]$distr, case = "young")
  }),
  ncol = 2,
  align = "v"
)
ggsave(filename = paste0(dir_data_test, "t1_plot_specific_pval_spec_young_comp.png"), plot = t1_plot_specific_pval_spec_young_comp, scale = 2)
t1_plot_specific_pval_spec_old_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_spec_old_specific_t1), function(i) {
    comp_distr(cci_specific_distr_t1$distr_cond1[i,], pv_spec_old_specific_t1[[i]]$distr, case = "old")
  }),
  ncol = 2,
  align = "v"
)
ggsave(filename = paste0(dir_data_test, "t1_plot_specific_pval_spec_young_comp.png"), plot = t1_plot_specific_pval_spec_young_comp, scale = 2)
