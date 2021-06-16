####################################################
##
## Project: scAgeCom
##
## Last update - June 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## scDiffCom objects preparation
##
####################################################
##

## Libraries ####

library(Seurat)
library(scDiffCom)
library(future)

options(future.globals.maxSize = 10000 * 1024 ^ 2)

## Dataset path on Server ####

server_path <- "/home/nis/lcyril/work/lcyril_data/scRNA_seq/seurat_processed/"

dataset_path <- c(
  tms_facs = paste0(server_path, "seurat_final_tms_facs.rds"),
  tms_droplet = paste0(server_path, "seurat_final_tms_droplet.rds"),
  calico_kidney = paste0(server_path, "seurat_final_calico_kidney.rds"),
  calico_lung = paste0(server_path, "seurat_final_calico_lung.rds"),
  calico_spleen = paste0(server_path, "seurat_final_calico_spleen.rds")
)

## Output path ####

dir_data_output <- getwd()

## List of analysis to do over datasets and sex ####

analysis_list <- list(
  "calico_male",
  "facs_female",
  "facs_male",
  "droplet_female",
  "droplet_male"
)

## Do the analysis ####

for (analysis in analysis_list) {
  message(paste0("Starting the analysis of ", analysis))
  output_dir <- paste0(
    dir_data_output,
    "/scdiffcom_",
    analysis
  )
  if (!dir.exists(output_dir)) {
    message("Create output directory.")
    dir.create(output_dir)
  }
  message("Reading seurat object.")
  if (analysis == "calico_male") {
    seurat_kidney <- readRDS(dataset_path[["calico_kidney"]])
    seurat_kidney$age_group <- ifelse(
      seurat_kidney$age == "young",
      'YOUNG',
      'OLD'
    )
    seurat_lung <- readRDS(dataset_path[["calico_lung"]])
    seurat_lung$age_group <- ifelse(
      seurat_lung$age == "young",
      'YOUNG',
      'OLD'
    )
    seurat_spleen <- readRDS(dataset_path[["calico_spleen"]])
    seurat_spleen$age_group <- ifelse(
      seurat_spleen$age == "young",
      'YOUNG',
      'OLD'
    )
    tissue_list <- c("Kidney", "Lung", "Spleen")
    seurat_list <- list(
      Kidney = seurat_kidney,
      Lung = seurat_lung,
      Spleen = seurat_spleen
    )
    n_tissue <- 3
  } else if (analysis == "facs_female") {
    seurat_obj <- readRDS(dataset_path[["tms_facs"]])
    seurat_obj <- subset(seurat_obj, subset = sex == "female")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c('3m'),
      'YOUNG',
      'OLD'
    )
    tissue_list <- c(
      "BAT", "Bladder", "Brain", "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle",
      "Lung", "MAT", "Mammary_Gland", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "facs_male") {
    seurat_obj <- readRDS(dataset_path[["tms_facs"]])
    seurat_obj <- subset(seurat_obj, subset = sex == "male")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c('3m'),
      'YOUNG',
      'OLD'
    )
    tissue_list <- c(
      "Aorta", "BAT", "Bladder", "Brain", "Diaphragm",
      "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle", "Liver",
      "Lung", "MAT", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_female") {
    seurat_obj <- readRDS(dataset_path[["tms_droplet"]])
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj <- subset(seurat_obj, subset = sex == "female")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c('3m'),
      'YOUNG',
      'OLD'
    )
    tissue_list <- c(
      "Heart_and_Aorta", "Kidney", "Limb_Muscle",
      "Liver", "Lung", "Mammary_Gland", "Marrow",
      "Spleen", "Thymus"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_male") {
    seurat_obj <- readRDS(dataset_path[["tms_droplet"]])
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj <- subset(seurat_obj, subset = sex == "male")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c('3m'),
      'YOUNG',
      'OLD'
    )
    tissue_list <- c(
      "Bladder", "Kidney", "Liver", "Lung",
      "Spleen", "Tongue"
    )
    n_tissue <- length(tissue_list)
  }
  for (i in 1:n_tissue) {
    tiss <- tissue_list[i]
    message(paste0(
      "Analysis of the ",
      tiss,
      ". Tissue ",
      i,
      " out of ",
      n_tissue,
      "."
    ))
    if(analysis == "calico_male") {
      seurat_tiss <- seurat_list[[tiss]]
    } else {
      if (tiss == "Brain") {
        seurat_tiss <- subset(
          seurat_obj, subset = tissue %in% c(
            "Brain_Myeloid",
            "Brain_Non-Myeloid"
          )
        )
      } else {
        cells_tiss <- colnames(seurat_obj)[which(seurat_obj$tissue == tiss)]
        seurat_tiss <- subset(seurat_obj, cells = cells_tiss)
      }
    }
    cell_type_id <- "cell_ontology_final"
    condition_id <- list(
      column_name = "age_group",
      cond1_name = "YOUNG",
      cond2_name = "OLD"
    )
    DefaultAssay(seurat_tiss) <- "RNA"
    n_ct <- length(unique(seurat_tiss[[cell_type_id]][[1]]))
    if (n_ct > 1) {
      message("Size-factor normalization:")
      seurat_tiss <- NormalizeData(seurat_tiss, assay = "RNA")
      saveRDS(seurat_tiss[[]], file = paste0(output_dir, "/md_", tiss, ".rds"))
      future::plan(multicore, workers = 30)
      scdf_res <- scDiffCom::run_interaction_analysis(
        seurat_object = seurat_tiss,
        LRI_species = "mouse",
        seurat_celltype_id = cell_type_id,
        seurat_condition_id = condition_id,
        iterations = 10000, 
        scdiffcom_object_name = tiss,      
        seurat_assay = "RNA",
        seurat_slot = "data",
        log_scale = FALSE,
        score_type = "geometric_mean",
        threshold_min_cells = 5,
        threshold_pct = 0.1,
        threshold_quantile_score = 0.2,
        threshold_p_value_specificity = 0.05,
        threshold_p_value_de = 0.05,
        threshold_logfc = log(1.5),
        return_distributions = FALSE,
        seed = 42,
        verbose = TRUE
      )
      message(paste0("Saving results for the ", tiss, "."))
      saveRDS(scdf_res, file = paste0(output_dir, "/scdiffcom_", tiss, ".rds"))
    } else {
      message("Not enough cell types, not performing the analysis.")
    }
    future::plan(sequential)
  }
}
