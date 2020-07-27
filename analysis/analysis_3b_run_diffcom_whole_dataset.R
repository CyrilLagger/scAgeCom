####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - July 2020
##
## Apply scDiffCom to all three datasets globally!
##
####################################################
##

## Libraries ####

library(Seurat)
library(scDiffCom)
#library(foreach)
#library(doParallel)

## Parameters ####

#Note: change each time before running the analysis

dataset <- "tms_facs" # c("tms_facs", "tms_droplet", "calico")
LR_db <- "mixed" # c("all", "scsr", "sctensor", "nichenet", "mixed")
normalization <- "size_factor" #c("size_factor", "sctransform")
is_log <- TRUE
n_iter <- 10
run_test <- FALSE
calico_subtype <- FALSE
dir_data_output <- getwd()

## Data path ####

paths <- c(
  tms_facs = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_tms_facs.rds",
  tms_droplet = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_tms_droplet.rds",
  calico_kidney = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_kidney.rds",
  calico_lung = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_lung.rds",
  calico_spleen = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_spleen.rds"
)

## Create output path ####

if(is_log) {
  if(calico_subtype) {
    output_dir <- paste0(dir_data_output, "/diffcom_global_", dataset, "_subtype_", normalization, "_log_", n_iter, "iter_", LR_db)
  } else {
    output_dir <- paste0(dir_data_output, "/diffcom_global_", dataset, "_", normalization, "_log_", n_iter, "iter_", LR_db)
  }
} else {
  if(calico_subtype) {
    output_dir <- paste0(dir_data_output, "/diffcom_global_", dataset, "_subtype_", normalization, "_nlog_", n_iter, "iter_", LR_db)
  } else {
    output_dir <- paste0(dir_data_output, "/diffcom_global_", dataset, "_", normalization, "_nlog_", n_iter, "iter_", LR_db)
  }
}
if(run_test) {
  output_dir <- paste0(output_dir, "_test")
}

if(!dir.exists(output_dir)) {
  message("Create output directory.")
  dir.create(output_dir)
}

## Read LR-pairs ####

message("Read ligand-receptor database.")
LR <- LRall
if(LR_db == "all") {
  LR <- LR[, c("GENESYMB_L","GENESYMB_R", "SYMB_LR")]
} else if(LR_db == "mixed") {
  LR <- LR[LR$scsr == TRUE | LR$cpdb == TRUE, c("GENESYMB_L","GENESYMB_R", "SYMB_LR")]
} else {
  LR <- LR[LR[[LR_db]], c("GENESYMB_L","GENESYMB_R", "SYMB_LR")]
}

## Read Seurat object ####

if(dataset %in% c("tms_facs", "tms_droplet")) {
  message("Read seurat object.")
  seurat_obj <- readRDS(paths[[dataset]])
  seurat_obj$age_group <- ifelse(seurat_obj$age %in% c('1m', '3m'), 'young', 'old')
  seurat_obj$tissue_cell_type <- paste(seurat_obj$tissue, seurat_obj$cell_ontology_class, sep = "_")
  cell_type_id <- "tissue_cell_type"
  #we only keep the tissues with young/old cells
  if(dataset == "tms_facs") {
    #tissue_list <- unique(as.character(seurat_obj$tissue))
    tissue_list <- c("Pancreas", "Heart", "MAT", "SCAT", "GAT")
    #c("Aorta", "BAT", "Bladder", "Brain_Myeloid",
    #  "Brain_Non-Myeloid", "Diaphragm", "GAT",
    #  "Heart", "Kidney", "Large_Intestine",
    #  "Limb_Muscle", "Liver", "Lung", "Mammary_Gland",
    #  "Marrow", "MAT", "Pancreas", "SCAT", "Skin",
    #  "Spleen", "Thymus", "Tongue", "Trachea")
  } else {
    tissue_list <- c(
      "Bladder", "Heart_and_Aorta", "Kidney",
      "Limb_Muscle", "Liver", "Lung",
      "Mammary_Gland", "Marrow", "Spleen",
      "Thymus", "Tongue"
    )
  }
  seurat_obj <- subset(seurat_obj, subset = tissue %in% tissue_list)
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  message("Size-factor normalization:")
  seurat_obj <- Seurat::NormalizeData(seurat_obj, assay = "RNA")
  message("scDiffCom analysis:")
  dt_res <- scDiffCom::run_diffcom(
    seurat_obj = seurat_obj,
    LR_data = LR,
    seurat_cell_type_id = cell_type_id,
    condition_id = "age_group",
    assay = "RNA",
    slot = "data",
    log_scale = is_log,
    min_cells = 5,
    threshold = 0.1,
    permutation_analysis = TRUE,
    one_sided = FALSE,
    iterations = n_iter,
    return_distr = FALSE
  )
  message("Saving global results.")
  saveRDS(dt_res, file = paste0(output_dir, "/diffcom_global.rds"))
}
