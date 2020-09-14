####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - September 2020
##
## Double-check that scDiffCom returns correct results.
##
####################################################
##

## Libraries ####
library(Seurat)
library(data.table)
library(scDiffCom)
library(ggplot2)
library(pbapply)

## Data path ####
dir_data <- "../data_scAgeCom/"

## Load a Seurat objects ####

#here the test file corresponds to the Liver tissue from TMS FACS data.
seurat_test_1 <- readRDS(paste0(dir_data, "data_seurat_example.rds"))
seurat_test_1 <- NormalizeData(seurat_test_1, assay = "RNA")
seurat_test_1$age_group <- ifelse(seurat_test_1$age %in% c('1m', '3m'), 'young', 'old' )
seurat_test_1$cell_ontology_class <- as.character(seurat_test_1$cell_ontology_class)

## Load the LR database from scDiffCom ####
LR_test_1 <- scDiffCom::LR6db$LR6db_curated

## Typical usage of scDiffCom (1000 iterations) ####

#NOTES: --running scDiffcom takes a few minutes for 1000 iterations (e.g: around 4 min on a laptop for the parameters below)
#       --the command does not need to be ran each time if the results have been saved already
#       --uncomment the lines if you want to test it

# ?run_diffcom
# start_time <- Sys.time()
# diffcom_test_1 <- run_diffcom(
#   seurat_object = seurat_test_1,
#   LR_data = LR_test_1,
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
# diffcom_test_1_distr <- run_diffcom(
#   seurat_object = seurat_test_1,
#   LR_data = LR_test_1,
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

## Save the results or load previously saved results ####

#saveRDS(object = diffcom_test_1, file = paste0(dir_data, "test/test_1_data_scDiffcom_1000iter.rds"))
#saveRDS(object = diffcom_test_1_distr, file = paste0(dir_data, "test/test_1_data_scDiffcom_distr_1000iter.rds"))

diffcom_test_1 <- readRDS(paste0(dir_data, "test/test_1_data_scDiffcom_1000iter.rds"))
diffcom_test_1_distr <- readRDS(paste0(dir_data, "test/test_1_data_scDiffcom_distr_1000iter.rds"))

#consider only detected interactions and check that they correspond to the distributions
diffcom_test_1_detected <- diffcom_test_1[LR_DETECTED_young == TRUE | LR_DETECTED_old == TRUE]
identical(diffcom_test_1_detected$LR_SCORE_young - diffcom_test_1_detected$LR_SCORE_old,
          diffcom_test_1_distr$distr_diff[, 1001])
identical(diffcom_test_1_detected$LR_SCORE_old,
          diffcom_test_1_distr$distr_cond1[, 1001])
identical(diffcom_test_1_detected$LR_SCORE_young,
          diffcom_test_1_distr$distr_cond2[, 1001])

## Select groups of CCIs for further comparison ####

cci_test_choosen_list <- list(
  diffcom_test_1[LR_DETECTED_young & LR_DETECTED_old][28],
  diffcom_test_1[!LR_DETECTED_young & LR_DETECTED_old][16103],
  diffcom_test_1[LR_DETECTED_young & !LR_DETECTED_old][345],
  diffcom_test_1[LR_DETECTED_young & !LR_DETECTED_old][1234],
  diffcom_test_1[!LR_DETECTED_young & !LR_DETECTED_old][90]
)

cci_test_detected <- diffcom_test_1_detected[sample(.N, 5)]


## Double-check that the LR scores are correct ####

# We create a (slow) function that compute the LR scores in a different way than scDiffCom, such that
# we can double check the results.
compare_cci_score_test <- function(
  cci_test,
  condition,
  log_scale
) {
  seurat_sub_l <- subset(seurat_test_1, features = na.omit(unlist(cci_test[, c("LIGAND_1", "LIGAND_2")])),
                         subset = (cell_ontology_class == cci_test[["L_CELLTYPE"]] &
                                     age_group == condition))
  seurat_sub_r <- subset(seurat_test_1, features = na.omit(unlist(cci_test[, c("RECEPTOR_1", "RECEPTOR_2")])),
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

#should return zero or numerically close to zero
lapply(cci_test_choosen_list, compare_cci_score_test, "young", TRUE)
lapply(cci_test_list, compare_cci_score_test, "old", TRUE)
lapply(1:nrow(cci_test_detected), function(i) {
  compare_cci_score_test(cci_test_detected[i], "young", TRUE)
})
lapply(1:nrow(cci_test_detected), function(i) {
  compare_cci_score_test(cci_test_detected[i], "old", TRUE)
})


## Double-check that the LR differential p-values are correct ####

#We create a (very slow) function that compute the LR differential p-value (and also returns the distribution)
#in a different way than scDiffCom.
compare_cci_pvalue_diff_test <- function(
  seurat_obj,
  cci_test,
  iterations,
  log_scale
) {
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

test_pval_diff <- lapply(1:nrow(cci_test_detected), function(i) {
  compare_cci_pvalue_diff_test(seurat_test_1, cci_test_detected[i], iterations = 100, TRUE)
})

#compare the p-values 
mapply(function(x,y) {x$pval-cci_test_detected[y]$PVAL_DIFF}, test_pval_diff, 1:nrow(cci_test_detected))
lapply(test_pval_diff, function(x) {x$pval})
cci_test_detected$PVAL_DIFF

#compare the distributions
distr_plots <- lapply(
  1:length(test_pval_diff),
  function(i) {
    ggplot(data.frame(x=test_pval_diff[[i]]$distr), aes(x=x)) + geom_histogram(bins=50) +
      geom_vline(xintercept = cci_test_detected[i]$LR_SCORE_old - cci_test_detected[i]$LR_SCORE_young ) +
      xlab(expression(xi[old]-xi[young])) +
      ylab("Counts")+
      theme(text=element_text(size=20)) 
  }
)
g_distr_plots <- cowplot::plot_grid(
  plotlist = distr_plots,
  ncol = 2,
  align = "v"
)
#ggsave(filename = paste0(dir_data, "test/test_1_plot_distr_permutations.png"), plot = g_distr_plots)

#compare the distributions when existing
distr_id <- lapply(1:nrow(cci_test_detected), function(x) {
  id <- which(diffcom_test_1_distr$distr_cond1[,1001] == cci_test_detected[x]$LR_SCORE_old)
  if(length(id) != 1){
    return(NULL)
  } else
    return(id)
})

comp_distr <- function(distr1, distr2, case) {
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
    xlab(expr) +
    ylab("Counts") +
    theme(text=element_text(size=20)) 
}

plot_comp_distr_diff <- lapply(1:nrow(cci_test_detected), function(i) {
  comp_distr(-diffcom_test_1_distr$distr_diff[distr_id[[i]],], test_pval_diff[[i]]$distr, "diff")
})
g_plot_comp_distr_diff <- cowplot::plot_grid(plotlist = plot_comp_distr_diff, ncol = 2, align = "v")
#ggsave(filename = paste0(dir_data, "test/test_1_plot_comp_distr_permutations.png"), plot = g_plot_comp_distr_diff, scale = 2)


## Double-check that the LR specificity p-values are correct ####

#for the comparison to be equivalent we need the same cell-types
seurat_test_1_sub <- subset(seurat_test_1, subset = cell_ontology_class %in% unique(diffcom_test$L_CELLTYPE))

#create a function that compute the LR specificity p-values (and also returns the distribution) in a different way than scDiffCom
compare_cci_pvalue_spec_test <- function(seurat_obj, cci_test, iterations, condition, log_scale) {
  cells_use <- colnames(seurat_obj)[seurat_obj$age_group == condition]
  seurat_temp <- subset(seurat_obj, features = c(cci_test[["L_GENE"]], cci_test[["R_GENE"]]),
                        cells = cells_use)
  get_shuffled_LR_ct <- function(permut) {
    if(permut) {
      seurat_temp$cell_ontology_class <- sample(as.character(seurat_temp$cell_ontology_class))
    }
    cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_class == cci_test[["L_CELLTYPE"]]]
    seurat_sub_l <- subset(seurat_temp, features = cci_test[["L_GENE"]],
                           cells = cells_use1)
    cells_use1 <- colnames(seurat_temp)[seurat_temp$cell_ontology_class == cci_test[["R_CELLTYPE"]]]
    seurat_sub_r <- subset(seurat_temp, features = cci_test[["R_GENE"]],
                           cells = cells_use1)
    if(log_scale) {
      score <- (mean((seurat_sub_l[["RNA"]]@data[1,])) + mean((seurat_sub_r[["RNA"]]@data[1,])))/2
    } else {
      score <- (mean(expm1(seurat_sub_l[["RNA"]]@data[1,])) + mean(expm1(seurat_sub_r[["RNA"]]@data[1,])))/2
    }
    return(score)
  }
  all_LR <- replicate(iterations, get_shuffled_LR_ct(TRUE))
  all_LR <- c(all_LR, get_shuffled_LR_ct(FALSE))
  pval <- sum(all_LR[1:iterations] >= all_LR[(iterations + 1)]) / iterations
  return(list(pval = pval,
              distr = all_LR))
}

#takes a bit of time, do 1000 iterations to be comparable
test_pval_spec_old <- lapply(cci_test_list, function(x) {
  compare_cci_pvalue_spec_test(seurat_test_1_sub, x, iterations = 1000, "old", TRUE)
})
test_pval_spec_young <- lapply(cci_test_list, function(x) {
  compare_cci_pvalue_spec_test(seurat_test_1_sub, x, iterations = 1000, "young", TRUE)
})

#compare the p-values 
mapply(function(x,y) {x$pval-y$PVAL_old}, test_pval_spec_old, cci_test_list)
lapply(test_pval_spec_old, function(x) {x$pval})
lapply(cci_test_list, function(x) {x$PVAL_old})

mapply(function(x,y) {x$pval-y$PVAL_young}, test_pval_spec_young, cci_test_list)
lapply(test_pval_spec_young, function(x) {x$pval})
lapply(cci_test_list, function(x) {x$PVAL_young})

#compare the distributions when existing

plot_comp_distr_old <- lapply(1:6, function(i) {
  comp_distr(diffcom_test_distr$distr_cond1[distr_id[[i]],], test_pval_spec_old[[i]]$distr, "old")
})
g_plot_comp_distr_old <- cowplot::plot_grid(plotlist = plot_comp_distr_old, ncol = 2, align = "v")
#ggsave(filename = paste0(dir_data, "test/test_1_plot_comp_distr_old_permutations.png"), plot = g_plot_comp_distr_old, scale = 2)

plot_comp_distr_young <- lapply(1:6, function(i) {
  comp_distr(diffcom_test_distr$distr_cond2[distr_id[[i]],], test_pval_spec_young[[i]]$distr, "young")
})
g_plot_comp_distr_young <- cowplot::plot_grid(plotlist = plot_comp_distr_young, ncol = 2, align = "v")
#ggsave(filename = paste0(dir_data, "test/test_1_plot_comp_distr_young_permutations.png"), plot = g_plot_comp_distr_young, scale = 2)
