####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - June 2020
##
## Main scAgeCom analysis by running scDiffCom.
## Creates CCI data.frames to be used for
## downstream analysis.
##
####################################################
##

####################################################
## Source the libraries and parameters 
####################################################
##

library(Seurat)
library(scDiffCom)
library(foreach)
library(doParallel)


#parameters
dataset <- "tms_droplet" # c("tms_facs", "tms_droplet", "calico")
LR_db <- "mixed" # c("all", "sCsR", "scTensor", "nichenet", "mixed")
normalization <-  "size_factor" #c("size_factor", "sctransform")
is_log <- TRUE
n_iter <- 10000
run_test <- FALSE

#paths
paths <- c(tms_facs = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_tms_facs.rds",
           tms_droplet = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_tms_droplet.rds",
           calico_kidney = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_kidney.rds",
           calico_lung = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_lung.rds",
           calico_spleen = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_spleen.rds" )
if(is_log) {
  output_dir <- paste0(getwd(), "/diffcom_", dataset, "_", normalization, "_log_", n_iter, "iter_", LR_db)
} else {
  output_dir <- paste0(getwd(), "/diffcom_", dataset, "_", normalization, "_nlog_", n_iter, "iter_", LR_db)
}
if(run_test) {
  output_dir <- paste0(output_dir, "_test")
}

if(!dir.exists(output_dir)) {
  message("Create output directory.")
  dir.create(output_dir)
}

####################################################
## Code starts
####################################################
##

message("Read ligand-receptor database.")
LR <- LRall$LRall_one2one
if(LR_db == "all") {
  LR <- LR[, c(2,3,1)]
} else if(LR_db == "mixed") {
  LR <- LR[!(LR$scTensor & !LR$sCsR & !LR$nichenet & !LR$cpdb), c(2,3,1)]
} else {
  LR <- LR[LR[[LR_db]], c(2,3,1)]
}

message("Read seurat object.")
if(dataset %in% c("tms_facs", "tms_droplet")) {
  seurat_obj <- readRDS(paths[[dataset]])
  seurat_obj$age_group <- ifelse(seurat_obj$age %in% c('1m', '3m'), 'young', 'old')
  tissue_list <- unique(as.character(seurat_obj$tissue))
  #tissue_list <- c( "GAT","Lung", "Marrow")
  n_tissue <- length(tissue_list)
  
} else if(dataset == "calico") {
  seurat_kidney <- readRDS(paths[["calico_kidney"]])
  seurat_kidney$age_group <- as.character(seurat_kidney$age)
  seurat_lung <- readRDS(paths[["calico_lung"]])
  seurat_lung$age_group <- as.character(seurat_lung$age)
  seurat_spleen <- readRDS(paths[["calico_spleen"]])
  seurat_spleen$age_group <- as.character(seurat_spleen$age)
  tissue_list <- c("kidney", "lung", "spleen")
  seurat_list <- list(kidney = seurat_kidney,
                      lung = seurat_lung,
                      spleen = seurat_spleen)
  n_tissue <- 3
} else {
  stop("Dataset not supported.")
}

if(run_test){
  tissue_list <- tissue_list[1:2]
  n_tissue <- 2
}

message("Setup parallel environment.")
registerDoParallel(cores= min(n_tissue, 30))

message("Start tissue per tissue scDiffCom analysis")
foreach(i = 1:n_tissue, .packages = c("Seurat", "scDiffCom")) %dopar% {
  tiss <- tissue_list[i]
  message(paste0("Analysis of the ", tiss, ". Tissue ", i, " out of ", n_tissue, "."))
  if(dataset %in% c("tms_facs", "tms_droplet")) {
    message("Subset Seurat object")
    cells_tiss <- colnames(seurat_obj)[which(seurat_obj$tissue == tiss)]
    seurat_tiss <- subset(seurat_obj, cells = cells_tiss)
    cell_type_id <- "cell_ontology_class"
  } else if(dataset == "calico") {
    seurat_tiss <- seurat_list[[tiss]]
    cell_type_id <- "cell_type"
  } else {
    stop("Dataset not supported.")
  }
  if(normalization == "size_factor") {
    Seurat::DefaultAssay(seurat_tiss) <- "RNA"
    message("Size-factor normalization:")
    seurat_tiss <- Seurat::NormalizeData(seurat_tiss, assay = "RNA")
    dt_res <- scDiffCom::run_diffcom(seurat_obj = seurat_tiss,
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
                                     return_distr = FALSE)
    message(paste0("Saving results for the ", tiss, "."))
    saveRDS(dt_res, file = paste0(output_dir, "/diffcom_", tiss, ".rds"))
    
  } else if(normalization == "sct") {
    stop("To do later on.")
  }
  NULL
}
