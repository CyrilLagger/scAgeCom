##
## Project: scAgeCom
##
## Last update - March 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## Apply scDiffCom to all three datasets.
##
####################################################
##

## Libraries ####

library(Seurat)
library(scDiffCom)
library(future)

options(future.globals.maxSize = 10000 * 1024 ^ 2)

## Common parameters ####

normalization <- "size_factor"
min_cells <- 5
is_log <- FALSE
n_iter <- 10000
dir_data_output <- getwd()
run_test <- FALSE

## Path where the Seurat files are stored ####

dataset_paths <- c(
  tms_facs = "/home/nis/lcyril/work/lcyril_data/scRNA_seq/seurat_processed/seurat_final_tms_facs.rds",
  tms_droplet = "/home/nis/lcyril/work/lcyril_data/scRNA_seq/seurat_processed/seurat_final_tms_droplet.rds",
  calico_kidney = "/home/nis/lcyril/work/lcyril_data/scRNA_seq/seurat_processed/seurat_final_calico_kidney.rds",
  calico_lung = "/home/nis/lcyril/work/lcyril_data/scRNA_seq/seurat_processed/seurat_final_calico_lung.rds",
  calico_spleen = "/home/nis/lcyril/work/lcyril_data/scRNA_seq/seurat_processed/seurat_final_calico_spleen.rds"
)

## Load the curated LR interactions ####

LRI_curated <- scDiffCom::LRI_mouse$LRI_curated

## List of analysis to do over datasets and sex ####

analysis_list <- list(
  list(dataset = "calico", type = "default"),
  list(dataset = "tms_facs", type = "female"),
  list(dataset = "tms_facs", type = "male"),
  list(dataset = "tms_droplet", type = "female"),
  list(dataset = "tms_droplet", type = "male")#,
  #list(dataset = "tms_facs", type = "mixed"),
  #list(dataset = "tms_facs", type = "female"),
  #list(dataset = "tms_facs", type = "male"),
  #list(dataset = "tms_facs", type = "sex"),
  #list(dataset = "tms_droplet", type = "mixed"),
  #list(dataset = "tms_droplet", type = "sex")
  #list(dataset = "tms_facs", type = "overall")#,
  #list(dataset = "tms_droplet", type = "overall")
)

## Do the analysis ####

for(analysis in analysis_list) {
  message(paste0("Starting the analysis of ", analysis$dataset, "_", analysis$type))
  if(!run_test) {
    output_dir <- paste0(
      dir_data_output, "/scdiffcom_",
      analysis$dataset, "_", analysis$type, "_",
      normalization, "_log", is_log, "_", n_iter, "iter"
    )
  } else {
    output_dir <- paste0(dir_data_output, "/scdiffcom_test_run")
  }
  if(!dir.exists(output_dir)) {
    message("Create output directory.")
    dir.create(output_dir)
  }
  message("Read seurat object.")
  if(analysis$dataset %in% c("tms_facs", "tms_droplet")) {
    seurat_obj <- readRDS(dataset_paths[[analysis$dataset]])
    if (analysis$dataset == "tms_droplet") {
      seurat_obj <- subset(seurat_obj, subset = age %in% c("3m", "18m", "21m", "24m"))
    }
    if (any(!unique(seurat_obj$age) %in% c("3m", "18m", "21m", "24m"))) {
      stop("Error somewhere with the ages")
    }
    seurat_obj$age_group <- ifelse(seurat_obj$age %in% c('3m'), 'YOUNG', 'OLD')
    if(analysis$type == "female") {
      seurat_obj <- subset(seurat_obj, subset = sex == "female")
    }
    if(analysis$type == "male") {
      seurat_obj <- subset(seurat_obj, subset = sex == "male")
    }
    if(analysis$type == "sex") {
      group_filter <- "sex"
    } else {
      group_filter <- "age_group"
    }
    if(analysis$type == "overall") {
      seurat_obj$tissue_cell_ontology_scdiffcom <- paste(
        seurat_obj$tissue,
        seurat_obj$cell_ontology_scdiffcom,
        sep = "_")
      tissue_list <- "overall"
      n_tissue <- 1
    } else {
      md_temp <- seurat_obj@meta.data
      tokeep <- apply(
        table(
          as.character(md_temp$tissue),
          as.character(md_temp[[group_filter]])
          ) >= min_cells,
        MARGIN = 1,
        FUN = all
      )
      tissue_list <- names(tokeep[tokeep])
      n_tissue <- length(tissue_list)
    }
  }
  if(analysis$dataset == "calico") {
    seurat_kidney <- readRDS(dataset_paths[["calico_kidney"]])
    seurat_kidney$age_group <- ifelse(seurat_kidney$age == "young", 'YOUNG', 'OLD')
    seurat_lung <- readRDS(dataset_paths[["calico_lung"]])
    seurat_lung$age_group <- ifelse(seurat_lung$age == "young", 'YOUNG', 'OLD')
    seurat_spleen <- readRDS(dataset_paths[["calico_spleen"]])
    seurat_spleen$age_group <- ifelse(seurat_spleen$age == "young", 'YOUNG', 'OLD')
    tissue_list <- c("Kidney", "Lung", "Spleen")
    seurat_list <- list(
      Kidney = seurat_kidney,
      Lung = seurat_lung,
      Spleen = seurat_spleen
    )
    n_tissue <- 3
  }
  if(run_test){
    #tissue_list <- tissue_list[1:3]
    tissue_list <- tissue_list[tissue_list == "Lung"]
    n_tissue <- 1
  }
  for(i in 1:n_tissue) {
    tiss <- tissue_list[i]
    message(paste0("Analysis of the ", tiss, ". Tissue ", i, " out of ", n_tissue, "."))
    if(analysis$dataset %in% c("tms_facs", "tms_droplet")) {
      if(analysis$type == "overall") {
        seurat_tiss <- seurat_obj
        cell_type_id <- "tissue_cell_ontology_scdiffcom"
      } else {
        message("Subset Seurat object")
        cells_tiss <- colnames(seurat_obj)[which(seurat_obj$tissue == tiss)]
        seurat_tiss <- subset(seurat_obj, cells = cells_tiss)
        cell_type_id <- "cell_ontology_final"
      }
      if(analysis$type == "sex") {
        condition_id <- list(
          column_name = "sex",
          cond1_name = "female",
          cond2_name = "male"
        )
      } else {
        condition_id <- list(
          column_name = "age_group",
          cond1_name = "YOUNG",
          cond2_name = "OLD"
        )
      }
    } else if(analysis$dataset == "calico") {
      seurat_tiss <- seurat_list[[tiss]]
      condition_id <- list(
        column_name = "age_group",
        cond1_name = "YOUNG",
        cond2_name = "OLD"
      )
      if(analysis$type == "subtype") {
        cell_type_id <- "subtype"
      } else {
        cell_type_id <- "cell_ontology_final"
      }
    }
    DefaultAssay(seurat_tiss) <- "RNA"
    n_ct <- length(unique(seurat_tiss[[cell_type_id]]))
    
    if (n_ct > 1) {
      message("Size-factor normalization:")
      seurat_tiss <- NormalizeData(seurat_tiss, assay = "RNA")
      future::plan(multicore, workers = 27)
      dt_res <- scDiffCom::run_interaction_analysis(
        seurat_object = seurat_tiss,
        LRI_species = "mouse",
        seurat_celltype_id = cell_type_id,
        seurat_condition_id = condition_id,
        iterations = n_iter, 
        scdiffcom_object_name = tiss,      
        seurat_assay = "RNA",
        seurat_slot = "data",
        log_scale = is_log,
        score_type = "geometric_mean",
        threshold_min_cells = min_cells,
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
      saveRDS(dt_res, file = paste0(output_dir, "/scdiffcom_", tiss, ".rds"))
    } else {
      message("Not enough cell types, not performing the analysis.")
    }
    future::plan(sequential)
  }
}




