####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - November 2020
##
## Apply scDiffCom to all three datasets.
##
####################################################
##

## Libraries ####

library(Seurat)
library(scDiffCom)
library(future)
library(progressr)

future::plan(multiprocess)
options(future.globals.maxSize = 64 * 1000 ^ 3)
options(progressr.enable = TRUE)

## Set up fixed parameters that should not change ####

normalization <- "size_factor"
min_cells <- 5
is_log <- FALSE
n_iter <- 10000
dir_data_output <- getwd()
run_test <- FALSE

## Path where the Seurat files are stored ####

dataset_paths <- c(
  tms_facs = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_tms_facs.rds",
  tms_droplet = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_tms_droplet.rds",
  calico_kidney = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_calico_kidney.rds",
  calico_lung = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_calico_lung.rds",
  calico_spleen = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_calico_spleen.rds"
)

## Load the curated LR interactions ####

LR6db_curated <- scDiffCom::LR6db$LR6db_curated

## List of analysis to do over datasets and sex

analysis_list <- list(
  #list(dataset = "calico", type = "default"),
  #list(dataset = "calico", type = "subtype")#,
  #list(dataset = "tms_facs", type = "mixed"),
  #list(dataset = "tms_facs", type = "female"),
  #list(dataset = "tms_facs", type = "male"),
  #list(dataset = "tms_facs", type = "sex"),
  #list(dataset = "tms_droplet", type = "mixed")#,
  #list(dataset = "tms_droplet", type = "female"),
  #list(dataset = "tms_droplet", type = "male"),
  #list(dataset = "tms_droplet", type = "sex"),
  list(dataset = "tms_facs", type = "overall"),
  list(dataset = "tms_droplet", type = "overall")
)

## Do the analysis ####

for(analysis in analysis_list) {
  message(paste0("Start analysis of ", analysis$dataset, "_", analysis$type))
  if(!run_test) {
    output_dir <- paste0(dir_data_output, "/scdiffcom_",
                         analysis$dataset, "_", analysis$type, "_",
                         normalization, "_log", is_log, "_", n_iter, "iter")
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
    seurat_obj$age_group <- ifelse(seurat_obj$age %in% c('1m', '3m'), 'YOUNG', 'OLD')
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
      seurat_obj$tissue_cell_ontology_scdiffcom <- paste(seurat_obj$tissue, seurat_obj$cell_ontology_scdiffcom, sep = "_")
      tissue_list <- "overall"
      n_tissue <- 1
    } else {
      md_temp <- seurat_obj@meta.data
      tokeep <- apply(
        table(as.character(md_temp$tissue), as.character(md_temp[[group_filter]])) >= min_cells,
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
        cell_type_id <- "cell_ontology_scdiffcom"
      }
      if(analysis$type == "sex") {
        condition_id <- "sex"
        cond1_name <- "female"
        cond2_name <- "male"
      } else {
        condition_id <- "age_group"
        cond1_name <- "YOUNG"
        cond2_name <- "OLD"
      }
    } else if(analysis$dataset == "calico") {
      seurat_tiss <- seurat_list[[tiss]]
      condition_id <- "age_group"
      cond1_name <- "YOUNG"
      cond2_name <- "OLD"
      if(analysis$type == "subtype") {
        cell_type_id <- "subtype"
      } else {
        cell_type_id <- "cell_ontology_scdiffcom"
      }
    }
    DefaultAssay(seurat_tiss) <- "RNA"
    message("Size-factor normalization:")
    seurat_tiss <- NormalizeData(seurat_tiss, assay = "RNA")
    dt_res <- scDiffCom::run_interaction_analysis(
      seurat_obj = seurat_tiss,
      LR_object = LR6db_curated,
      celltype_column_id = cell_type_id,
      condition_column_id = condition_id,
      cond1_name = cond1_name,
      cond2_name = cond2_name,
      object_name = tiss,
      assay = "RNA",
      slot = "data",
      log_scale = is_log,
      min_cells = min_cells,
      pct_threshold = 0.1,
      permutation_analysis = TRUE,
      iterations = n_iter,
      cutoff_quantile_score = 0.25,
      cutoff_pval_specificity = 0.05,
      cutoff_pval_de = 0.05,
      cutoff_logfc = log(1.2),
      return_distr = FALSE,
      seed = 42,
      verbose = TRUE,
      sparse = TRUE
    )
    message(paste0("Saving results for the ", tiss, "."))
    saveRDS(dt_res, file = paste0(output_dir, "/scdiffcom_", tiss, ".rds"))
  }
}




