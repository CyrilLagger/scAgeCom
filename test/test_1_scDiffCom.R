####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - September 2020
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

## Data path and parameters####
dir_data_test <- "../data_scAgeCom/test/"

n_iter <- 1000 #number of permutations for the double-check

## Load a Seurat objects ####

#here the test file corresponds to the Liver from TMS FACS data.
seurat_t1 <- readRDS(paste0(dir_data_test, "data_seurat_example_facs_liver.rds"))
seurat_t1 <- NormalizeData(seurat_t1, assay = "RNA")
seurat_t1$age_group <- ifelse(seurat_t1$age %in% c('1m', '3m'), 'young', 'old' )
seurat_t1$cell_ontology_class <- as.character(seurat_t1$cell_ontology_class)

## Load the LR database from scDiffCom ####
LR_t1 <- scDiffCom::LR6db$LR6db_curated

## Typical usage of scDiffCom (1000 iterations) ####

#NOTES: --running scDiffcom takes a few minutes for 1000 iterations (e.g: around 4 min on a laptop for the parameters below)
#       --the command does not need to be run each time if the results have been saved already
#       --uncomment the lines if you want to test it

# ?run_diffcom
# start_time <- Sys.time()
# diffcom_t1 <- run_diffcom(
#   seurat_object = seurat_t1,
#   LR_data = LR_t1,
#   seurat_cell_type_id = "cell_ontology_class",
#   condition_id = "age_group",
#   assay = "RNA",
#   slot = "data",
#   log_scale = TRUE,
#   min_cells = 5,
#   threshold = 0.1,
#   permutation_analysis = TRUE,
#   one_sided = FALSE,
#   iterations = 1000,
#   return_distr = FALSE
# )
# end_time <- Sys.time()
# end_time - start_time

#run_diffcom but returning the distributions of the permutation test (only the detected ones)
# diffcom_t1_distr <- run_diffcom(
#   seurat_object = seurat_t1,
#   LR_data = LR_t1,
#   seurat_cell_type_id = "cell_ontology_class",
#   condition_id = "age_group",
#   assay = "RNA",
#   slot = "data",
#   log_scale = TRUE,
#   min_cells = 5,
#   threshold = 0.1,
#   permutation_analysis = TRUE,
#   one_sided = FALSE,
#   iterations = 1000,
#   return_distr = TRUE
# )

## Preprocessing, filtering and saving results ####
#not needed if already saved

#only keep detected interactions in either young or old samples
#diffcom_t1 <- diffcom_t1[LR_DETECTED_young == TRUE | LR_DETECTED_old == TRUE]
#diffcom_t1[, TISSUE := "Liver"]

#apply filtering analysis to get the different regulation cases
#source("src/src_1_filtering.R")
#quantile(c(diffcom_t1$LR_SCORE_old, diffcom_t1$LR_SCORE_young), 0.25)
#hist(c(diffcom_t1$LR_SCORE_old, diffcom_t1$LR_SCORE_young), breaks = 100)
#abline(v = quantile(c(diffcom_t1$LR_SCORE_old, diffcom_t1$LR_SCORE_young), 0.25) )

#diffcom_t1 <- analyze_CCI(
#  data = diffcom_t1,
#  cutoff_score = quantile(c(diffcom_t1$LR_SCORE_old, diffcom_t1$LR_SCORE_young), 0.25)
#)

#save results
#saveRDS(object = diffcom_t1, file = paste0(dir_data_test, "t1_data_scDiffcom_1000iter.rds"))
#saveRDS(object = diffcom_t1_distr, file = paste0(dir_data_test, "t1_data_scDiffcom_distr_1000iter.rds"))

## Load previously saved results ####

#read files
diffcom_t1 <- readRDS(paste0(dir_data_test, "t1_data_scDiffcom_1000iter.rds"))
diffcom_t1_distr <- readRDS(paste0(dir_data_test, "t1_data_scDiffcom_distr_1000iter.rds"))

#check tables are in the same order
identical(diffcom_t1$LR_SCORE_young - diffcom_t1$LR_SCORE_old,
          diffcom_t1_distr$distr_diff[, 1001])
identical(diffcom_t1$LR_SCORE_old,
          diffcom_t1_distr$distr_cond1[, 1001])
identical(diffcom_t1$LR_SCORE_young,
          diffcom_t1_distr$distr_cond2[, 1001])

diffcom_t1$rn <- as.numeric(rownames(diffcom_t1))

## Select (random) CCIs for further comparison ####

table(diffcom_t1$CASE_TYPE)

#5 random CCIs
cci_random_t1 <- rbindlist(
  list(
  diffcom_t1[
  L_NCELLS_old >=5 & R_NCELLS_old >= 5 & L_NCELLS_young >=5 & R_NCELLS_young >=5
  ][sample(.N, 2)],
  diffcom_t1[grepl("_", LR_NAME) & CASE_TYPE != "FFF" &
    L_NCELLS_old >=5 & R_NCELLS_old >= 5 & L_NCELLS_young >=5 & R_NCELLS_young >=5
    ][sample(.N, 2)]
  )
)

cci_random_distr_t1 <- lapply(diffcom_t1_distr, function(i) {
  i[cci_random_t1$rn,]
})
  
#"special cases" CCI
cci_choosen_t1 <- rbindlist(
  list(
    diffcom_t1[CASE_TYPE == "FFF" & L_NCELLS_young < 5 & L_NCELLS_young > 0 & R_NCELLS_young >0][3],
    diffcom_t1[CASE_TYPE != "FFF" & L_NCELLS_young < 5 & L_NCELLS_young > 0 & R_NCELLS_young >0][1],
    diffcom_t1[CASE_TYPE == "FTTU"][100],
    diffcom_t1[CASE_TYPE == "TFTD"][2],
    diffcom_t1[CASE_TYPE == "TTTD"][2],
    diffcom_t1[CASE_TYPE == "TTTU"][2],
    diffcom_t1[CASE_TYPE == "TTFD"][2],
    diffcom_t1[CASE_TYPE == "TTFU"][2]
  )
)

cci_choosen_distr_t1 <- lapply(diffcom_t1_distr, function(i) {
  i[cci_choosen_t1$rn,]
})

## Double-check that the LR scores are correct ####

# We create a (slow) function that compute the LR scores in a different way than scDiffCom, such that
# we can double check the results.
compare_cci_score_test <- function(
  cci_test_dt,
  condition,
  log_scale
) {
  res <- lapply(1:nrow(cci_test_dt), function(i) {
    cci_test <- cci_test_dt[i,]
    seurat_sub_l <- subset(seurat_t1, features = na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])),
                           subset = (cell_ontology_class == cci_test[["L_CELLTYPE"]] &
                                       age_group == condition))
    seurat_sub_r <- subset(seurat_t1, features = na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")])),
                           subset = (cell_ontology_class == cci_test[["R_CELLTYPE"]] &
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
compare_cci_score_test(cci_random_t1, "young", TRUE)
compare_cci_score_test(cci_random_t1, "old", TRUE)
compare_cci_score_test(cci_choosen_t1, "young", TRUE)
compare_cci_score_test(cci_choosen_t1, "old", TRUE)

## Double-check that the LR differential p-values are correct ####

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
    cells_use <- colnames(seurat_obj)[seurat_obj$cell_ontology_class  %in% c(cci_test[["L_CELLTYPE"]], cci_test[["R_CELLTYPE"]])]
    seurat_temp <- subset(seurat_obj, features = c(na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])), na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")]))),
                          cells = cells_use)
    get_shuffled_LR <- function(permut) {
      if(permut) {
        seurat_temp$age_group[seurat_temp$cell_ontology_class == cci_test[["L_CELLTYPE"]] ] <- 
          sample(seurat_temp$age_group[seurat_temp$cell_ontology_class == cci_test[["L_CELLTYPE"]] ])
        seurat_temp$age_group[seurat_temp$cell_ontology_class == cci_test[["R_CELLTYPE"]] ] <- 
          sample(seurat_temp$age_group[seurat_temp$cell_ontology_class == cci_test[["R_CELLTYPE"]] ])
      }
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_class %in% cci_test[["L_CELLTYPE"]] &
                                            seurat_temp$age_group == "young"]
      seurat_sub_l_young <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])),
                                   cells = cells_use1)
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_class %in% cci_test[["R_CELLTYPE"]] &
                                            seurat_temp$age_group == "young"]
      seurat_sub_r_young <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")])),
                                   cells = cells_use1)
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_class %in% cci_test[["L_CELLTYPE"]] &
                                            seurat_temp$age_group == "old"]
      seurat_sub_l_old <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])),
                                 cells = cells_use1)
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_class %in% cci_test[["R_CELLTYPE"]] &
                                            seurat_temp$age_group == "old"]
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

#We create a function that can display two distributions on the same plot
comp_distr <- function(
  distr1,
  distr2,
  case
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
    theme(text=element_text(size=20)) 
}

#compute the p-values and distributions
pv_diff_random_t1 <- compare_cci_pvalue_diff_test(seurat_t1, cci_random_t1, iterations = n_iter, log_scale = TRUE)
pv_diff_choosen_t1 <- compare_cci_pvalue_diff_test(seurat_t1, cci_choosen_t1, iterations = n_iter, log_scale = TRUE)

#compare the p-values 
lapply(seq_along(pv_diff_random_t1), function(i) {
  c(cci_random_t1[i,]$PVAL_DIFF, pv_diff_random_t1[[i]]$pval)
})
lapply(seq_along(pv_diff_choosen_t1), function(i) {
  c(cci_choosen_t1[i,]$PVAL_DIFF, pv_diff_choosen_t1[[i]]$pval)
})

#compare the distributions
t1_plot_random_pval_diff_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_diff_random_t1), function(i) {
    comp_distr(cci_random_distr_t1$distr_diff[i,], pv_diff_random_t1[[i]]$distr, case = "diff")
  }),
  ncol = 2,
  align = "v"
)
ggsave(filename = paste0(dir_data_test, "t1_plot_random_distr_permutations_comp.png"), plot = t1_plot_random_pval_diff_comp, scale = 2)
t1_plot_choosen_pval_diff_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_diff_choosen_t1), function(i) {
    comp_distr(cci_choosen_distr_t1$distr_diff[i,], pv_diff_choosen_t1[[i]]$distr, case = "diff")
  }),
  ncol = 2,
  align = "v"
)
ggsave(filename = paste0(dir_data_test, "t1_plot_choosen_distr_permutations_comp.png"), plot = t1_plot_choosen_pval_diff_comp, scale = 2)

## Double-check that the LR specificity p-values are correct ####

#We create a (very slow) function that compute the LR specificity p-value (and also returns the distribution)
#in a different way than scDiffCom.
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
        seurat_temp$cell_ontology_class <- sample(as.character(seurat_temp$cell_ontology_class))
      }
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_class == cci_test[["L_CELLTYPE"]]]
      seurat_sub_l <- subset(seurat_temp, features = na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])),
                             cells = cells_use1)
      cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_class == cci_test[["R_CELLTYPE"]]]
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

#compute the p-values and distributions
pv_spec_young_random_t1 <- compare_cci_pvalue_specific_test(seurat_t1, cci_random_t1, iterations = n_iter, condition = "young", log_scale = TRUE)
pv_spec_old_random_t1 <- compare_cci_pvalue_specific_test(seurat_t1, cci_random_t1, iterations = n_iter, condition = "old", log_scale = TRUE)
pv_spec_young_choosen_t1 <- compare_cci_pvalue_specific_test(seurat_t1, cci_choosen_t1, iterations = n_iter, condition = "young", log_scale = TRUE)
pv_spec_old_choosen_t1 <- compare_cci_pvalue_specific_test(seurat_t1, cci_choosen_t1, iterations = n_iter, condition = "old", log_scale = TRUE)

#compare the p-values 
lapply(seq_along(pv_spec_young_random_t1), function(i) {
  c(cci_random_t1[i,]$PVAL_young, pv_spec_young_random_t1[[i]]$pval)
})
lapply(seq_along(pv_spec_old_random_t1), function(i) {
  c(cci_random_t1[i,]$PVAL_old, pv_spec_old_random_t1[[i]]$pval)
})
lapply(seq_along(pv_spec_young_choosen_t1), function(i) {
  c(cci_choosen_t1[i,]$PVAL_young, pv_spec_young_choosen_t1[[i]]$pval)
})
lapply(seq_along(pv_spec_old_choosen_t1), function(i) {
  c(cci_choosen_t1[i,]$PVAL_old, pv_spec_old_choosen_t1[[i]]$pval)
})

#compare the distributions
t1_plot_random_pval_spec_young_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_spec_young_random_t1), function(i) {
    comp_distr(cci_random_distr_t1$distr_cond2[i,], pv_spec_young_random_t1[[i]]$distr, case = "young")
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
t1_plot_choosen_pval_spec_young_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_spec_young_choosen_t1), function(i) {
    comp_distr(cci_choosen_distr_t1$distr_cond2[i,], pv_spec_young_choosen_t1[[i]]$distr, case = "young")
  }),
  ncol = 2,
  align = "v"
)
ggsave(filename = paste0(dir_data_test, "t1_plot_choosen_pval_spec_young_comp.png"), plot = t1_plot_choosen_pval_spec_young_comp, scale = 2)
t1_plot_choosen_pval_spec_old_comp <- cowplot::plot_grid(
  plotlist = lapply(seq_along(pv_spec_old_choosen_t1), function(i) {
    comp_distr(cci_choosen_distr_t1$distr_cond1[i,], pv_spec_old_choosen_t1[[i]]$distr, case = "old")
  }),
  ncol = 2,
  align = "v"
)
ggsave(filename = paste0(dir_data_test, "t1_plot_choosen_pval_spec_young_comp.png"), plot = t1_plot_choosen_pval_spec_young_comp, scale = 2)
