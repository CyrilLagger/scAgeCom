####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## scDiffCom objects preparation
##
####################################################
##

## Add libraries ####

library(future)

## Set options ####

options(future.globals.maxSize = 15 * 1024^3)

## Check Seurat analysis objects ####

seurats_analysis

## scDiffCom output path ####

path_scagecom_output_scdiffcom <- paste0(
  path_scagecom_output,
  "scDiffCom_14_05_2022"
)

## List of analysis to do over datasets and sex ####

analysis_list <- list(
  "calico_diffage_male",
  "facs_diffage_female",
  "facs_diffage_male",
  "facs_diffage_mixed",
  "droplet_diffage_female",
  "droplet_diffage_male",
  "droplet_diffage_mixed"
)

## Do the analysis ####

for (analysis in analysis_list) {
  message(paste0("Starting the analysis of ", analysis))
  output_dir <- paste0(
    path_scagecom_output_scdiffcom,
    "/scdiffcom_",
    analysis
  )
  if (!dir.exists(output_dir)) {
    message("Create output directory.")
    dir.create(output_dir)
  }
  message("Reading seurat object.")
  if (analysis == "calico_diffage_male") {
    seurat_kidney <- seurats_analysis$calico_kidney
    seurat_kidney$age_group <- ifelse(
      seurat_kidney$age == "young",
      "YOUNG",
      "OLD"
    )
    seurat_lung <- seurats_analysis$calico_lung
    seurat_lung$age_group <- ifelse(
      seurat_lung$age == "young",
      "YOUNG",
      "OLD"
    )
    seurat_spleen <- seurats_analysis$calico_spleen
    seurat_spleen$age_group <- ifelse(
      seurat_spleen$age == "young",
      "YOUNG",
      "OLD"
    )
    tissue_list <- c("Kidney", "Lung", "Spleen")
    seurat_list <- list(
      Kidney = seurat_kidney,
      Lung = seurat_lung,
      Spleen = seurat_spleen
    )
    n_tissue <- 3
  } else if (analysis == "facs_diffage_female") {
    seurat_obj <- seurats_analysis$tms_facs
    seurat_obj <- subset(seurat_obj, subset = sex == "female")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "BAT", "Bladder", "Brain", "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle",
      "Lung", "MAT", "Mammary_Gland", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "facs_diffage_male") {
    seurat_obj <- seurats_analysis$tms_facs
    seurat_obj <- subset(seurat_obj, subset = sex == "male")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
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
  } else if (analysis == "facs_diffage_mixed") {
    seurat_obj <- seurats_analysis$tms_facs
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Aorta", "BAT", "Bladder", "Brain", "Diaphragm",
      "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle", "Liver",
      "Lung", "Mammary_Gland", "MAT", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_diffage_female") {
    seurat_obj <- seurats_analysis$tms_droplet
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj <- subset(seurat_obj, subset = sex == "female")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Heart_and_Aorta", "Kidney", "Limb_Muscle",
      "Liver", "Lung", "Mammary_Gland", "Marrow",
      "Spleen", "Thymus"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_diffage_male") {
    seurat_obj <- seurats_analysis$tms_droplet
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj <- subset(seurat_obj, subset = sex == "male")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Bladder", "Kidney", "Liver", "Lung",
      "Spleen", "Tongue"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_diffage_mixed") {
    seurat_obj <- seurats_analysis$tms_droplet
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Bladder", "Heart_and_Aorta", "Kidney", "Liver", "Limb_Muscle",
      "Liver", "Lung", "Mammary_Gland", "Marrow",
      "Spleen", "Tongue", "Thymus"
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
    if (analysis == "calico_diffage_male") {
      seurat_tiss <- seurat_list[[tiss]]
    } else {
      if (tiss == "Brain") {
        seurat_tiss <- subset(
          seurat_obj,
          subset = tissue %in% c(
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
      future::plan(multicore, workers = 24)
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

## Perform ICC differential analysis between male and female #####

analysis_list_sex <- list(
  "facs_diffsex_young",
  "facs_diffsex_old",
  "facs_diffsex_combined",
  "droplet_diffsex_young",
  "droplet_diffsex_old",
  "droplet_diffsex_all"
)

dt_metadata[
  ,
  age_group := ifelse(
    age %in% c("3m"), "YOUNG",
    ifelse(
      age %in% c("18m", "21m", "24m"), "OLD",
      "NOT"
    )
  )
]

dt_md_tms_age_sex <- dt_metadata[
  dataset %in% c("tms_facs", "tms_droplet") & age_group != "NOT",
  .N,
  by = c("dataset", "tissue", "cell_ontology_final", "sex", "age_group")
][order(dataset, tissue, sex, age_group)]

fwrite(
  dt_md_tms_age_sex,
  paste0(
    path_scagecom_output,
    "dt_md_tms_age_sex.csv"
  )
)

for (analysis in analysis_list_sex) {
  message(paste0("Starting the analysis of ", analysis))
  output_dir <- paste0(
    path_scagecom_output_scdiffcom,
    "/scdiffcom_",
    analysis
  )
  if (!dir.exists(output_dir)) {
    message("Create output directory.")
    dir.create(output_dir)
  }
  message("Reading seurat object.")
  if (analysis == "facs_diffsex_young") {
    seurat_obj <- seurats_analysis$tms_facs
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    seurat_obj <- subset(seurat_obj, subset = age_group == "YOUNG")
     tissue_list <- c(
      "Aorta", "BAT", "Bladder", "Brain", "Diaphragm",
      "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle",
      "Lung", "MAT", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "facs_diffsex_old") {
    seurat_obj <- seurats_analysis$tms_facs
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    seurat_obj <- subset(seurat_obj, subset = age_group == "OLD")
    tissue_list <- c(
      "Aorta", "BAT", "Bladder", "Brain",
      "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle", "Liver",
      "Lung", "MAT", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "facs_diffsex_combined") {
    seurat_obj <- seurats_analysis$tms_facs
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
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
  } else if (analysis == "droplet_diffsex_young") {
    seurat_obj <- seurats_analysis$tms_droplet
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    seurat_obj <- subset(seurat_obj, subset = age_group == "YOUNG")
    tissue_list <- c(
      "Bladder", "Kidney",
      "Liver", "Lung",
      "Spleen", "Tongue"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_diffsex_old") {
    seurat_obj <- seurats_analysis$tms_droplet
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    seurat_obj <- subset(seurat_obj, subset = age_group == "OLD")
    tissue_list <- c(
      "Heart_and_Aorta", "Kidney", "Limb_Muscle",
      "Liver", "Lung", "Marrow", "Pancreas", "Skin",
      "Spleen", "Thymus"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_diffsex_all") {
    seurat_obj <- seurats_analysis$tms_droplet
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Bladder", "Heart_and_Aorta", "Kidney", "Limb_Muscle",
      "Liver", "Lung", "Marrow", "Pancreas", "Skin",
      "Spleen", "Thymus", "Tongue"
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
    if (tiss == "Brain") {
      seurat_tiss <- subset(
        seurat_obj,
        subset = tissue %in% c(
          "Brain_Myeloid",
          "Brain_Non-Myeloid"
        )
      )
    } else {
      cells_tiss <- colnames(seurat_obj)[which(seurat_obj$tissue == tiss)]
      seurat_tiss <- subset(seurat_obj, cells = cells_tiss)
    }
    cell_type_id <- "cell_ontology_final"
    condition_id <- list(
      column_name = "sex",
      cond1_name = "female",
      cond2_name = "male"
    )
    DefaultAssay(seurat_tiss) <- "RNA"
    n_ct <- length(unique(seurat_tiss[[cell_type_id]][[1]]))
    if (n_ct > 1) {
      message("Size-factor normalization:")
      seurat_tiss <- NormalizeData(seurat_tiss, assay = "RNA")
      saveRDS(seurat_tiss[[]], file = paste0(output_dir, "/md_", tiss, ".rds"))
      future::plan(multicore, workers = 24)
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
