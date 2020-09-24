####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - September 2020
##
## Check the scRNA-seq dataset preprocessing.
##
####################################################
##

# Note: works on the server only due to large file-size

## Libraries ####

library(Seurat)

## Paths of the Seurat objects ####

seurat_path <- c(
  tms_facs = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_tms_facs.rds",
  tms_droplet = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_tms_droplet.rds",
  calico_kidney = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_kidney.rds",
  calico_lung = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_lung.rds",
  calico_spleen = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_spleen.rds"
)

## Analysis data path ####
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Load the Seurat objects (several GBs) #####

seurat_objects <- list(
  tms_facs = readRDS(seurat_path[["tms_facs"]]),
  tms_droplet = readRDS(seurat_path[["tms_droplet"]]),
  calico_kidney = readRDS(seurat_path[["calico_kidney"]]),
  calico_lung = readRDS(seurat_path[["calico_lung"]]),
  calico_spleen = readRDS(seurat_path[["calico_spleen"]])
)

## Check basic contents of Seurat objects ####

# Available assays: "RNA"
lapply(seurat_objects, function(x) {
  x@assays
})

# Content of "counts" slot: integers
lapply(seurat_objects, function(x) {
  x$RNA@counts[1:5,1:5]
})

# Content of metadata
lapply(seurat_objects, function(x) {
  colnames(x@meta.data)
})

## Check QC ####

# Number of unique genes per cell
unique_genes <- lapply(seurat_objects, function(x) {
  Matrix::colSums(x$RNA@counts > 0)
})
lapply(unique_genes, function(i) {head(sort(i))})

min_genes_per_cell <- lapply(unique_genes, min)
min_genes_per_cell

# Total number of counts per cell
total_counts <- lapply(seurat_objects, function(x) {
  Matrix::colSums(x$RNA@counts)
})
lapply(total_counts, function(i) {head(sort(i))})

min_total_count_per_cell <- lapply(total_counts, min)
min_total_count_per_cell

# Compare counts to previously stored values (some difference due to early filtering)
mapply(function(x,y) {identical(x$n_genes,y)}, seurat_objects, unique_genes, SIMPLIFY = FALSE)
identical(seurat_objects$tms_droplet$n_genes, unique_genes$tms_droplet)
droplet_n_genes_diff <- seurat_objects$tms_droplet$n_genes - unique_genes$tms_droplet
droplet_n_genes_diff[droplet_n_genes_diff != 0]


## Save metadata to be used in other scripts not on the server ####

seurat_md <- lapply(seurat_objects, function(i) {
  i@meta.data
})

saveRDS(seurat_md, paste0(dir_data_analysis, "analysis_1_data_seurat_md.rds"))
