####################################################
##
## Project: scAgeCom
##
## Last update - December 2020
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## Check integrity of TMS and Calico Seurat files.
## Rename cell-types for scDiffCom analysis.
##
####################################################
##

# Note: works on the server only due to large file-size

## Libraries ####

library(Seurat)
library(data.table)

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

## Add new cell-type naming convention (only do it once, then save the new Seurat objects) ####

celltype_conversion <- read.csv(paste0(dir_data_analysis, "scDiffCom_cell_types.csv"), stringsAsFactors = FALSE)
celltype_conversion <- celltype_conversion[, c("Tissue", "TMS.Calico.cell.type", "scDiffCom.cell.type")]
colnames(celltype_conversion) <- c("tissue", "cell_ontology_class", "cell_ontology_scdiffcom")
setDT(celltype_conversion)
celltype_conversion <- unique(celltype_conversion)

celltype_check <- celltype_conversion[, c("tissue", "cell_ontology_class")]
any(duplicated(celltype_check))

#
md_tms_facs <- seurat_objects$tms_facs[[]]
md_tms_facs$cell_ontology_class <- as.character(md_tms_facs$cell_ontology_class)
setDT(md_tms_facs)
md_tms_facs[
  celltype_conversion,
  on = c("tissue", "cell_ontology_class"),
  cell_ontology_scdiffcom := i.cell_ontology_scdiffcom
  ]
sort(table(md_tms_facs$cell_ontology_scdiffcom), decreasing = TRUE)
sort(unique(md_tms_facs$cell_ontology_scdiffcom))
anyNA(md_tms_facs$cell_ontology_scdiffcom)
md_tms_facs[, cell_ontology_scdiffcom := ifelse(is.na(cell_ontology_scdiffcom),
                                                cell_ontology_class,
                                                cell_ontology_scdiffcom)]
sort(table(md_tms_facs$cell_ontology_scdiffcom), decreasing = TRUE)
sort(unique(md_tms_facs$cell_ontology_scdiffcom))
anyNA(md_tms_facs$cell_ontology_scdiffcom)
identical(as.character(md_tms_facs$cell_ontology_class), as.character(seurat_objects$tms_facs$cell_ontology_class))
seurat_objects$tms_facs$cell_ontology_scdiffcom <- md_tms_facs$cell_ontology_scdiffcom

#
md_tms_droplet <- seurat_objects$tms_droplet[[]]
md_tms_droplet$cell_ontology_class <- as.character(md_tms_droplet$cell_ontology_class)
setDT(md_tms_droplet)
md_tms_droplet[
  celltype_conversion,
  on = c("tissue", "cell_ontology_class"),
  cell_ontology_scdiffcom := i.cell_ontology_scdiffcom
  ]
sort(table(md_tms_droplet$cell_ontology_scdiffcom), decreasing = TRUE)
sort(unique(md_tms_droplet$cell_ontology_scdiffcom))
anyNA(md_tms_droplet$cell_ontology_scdiffcom)
md_tms_droplet[, cell_ontology_scdiffcom := ifelse(is.na(cell_ontology_scdiffcom),
                                                   cell_ontology_class,
                                                   cell_ontology_scdiffcom)]
sort(table(md_tms_droplet$cell_ontology_scdiffcom), decreasing = TRUE)
sort(unique(md_tms_droplet$cell_ontology_scdiffcom))
anyNA(md_tms_droplet$cell_ontology_scdiffcom)
identical(as.character(md_tms_droplet$cell_ontology_class), as.character(seurat_objects$tms_droplet$cell_ontology_class))
seurat_objects$tms_droplet$cell_ontology_scdiffcom <- md_tms_droplet$cell_ontology_scdiffcom

#
md_calico_kidney <- seurat_objects$calico_kidney[[]]
md_calico_kidney$cell_type <- as.character(md_calico_kidney$cell_type)
setDT(md_calico_kidney)
md_calico_kidney[
  celltype_conversion[tissue == "Kidney"],
  on = c("cell_type==cell_ontology_class"),
  cell_ontology_scdiffcom := i.cell_ontology_scdiffcom
  ]
sort(table(md_calico_kidney$cell_ontology_scdiffcom), decreasing = TRUE)
sort(unique(md_calico_kidney$cell_ontology_scdiffcom))
anyNA(md_calico_kidney$cell_ontology_scdiffcom)
md_calico_kidney[, cell_ontology_scdiffcom := ifelse(is.na(cell_ontology_scdiffcom),
                                                     cell_type,
                                                     cell_ontology_scdiffcom)]
anyNA(md_calico_kidney$cell_ontology_scdiffcom)
md_calico_kidney[
  celltype_conversion[tissue == "Kidney"],
  on = c("subtype==cell_ontology_class"),
  subtype_ontology_scdiffcom := i.cell_ontology_scdiffcom
  ]
md_calico_kidney[, subtype_ontology_scdiffcom := ifelse(is.na(subtype_ontology_scdiffcom),
                                                        subtype,
                                                        subtype_ontology_scdiffcom)]
identical(md_calico_kidney$cell_ontology_scdiffcom, md_calico_kidney$subtype_ontology_scdiffcom)
identical(as.character(md_calico_kidney$cell_type), as.character(seurat_objects$calico_kidney$cell_type))
seurat_objects$calico_kidney$cell_ontology_scdiffcom <- md_calico_kidney$cell_ontology_scdiffcom

#
md_calico_lung <- seurat_objects$calico_lung[[]]
md_calico_lung$cell_type <- as.character(md_calico_lung$cell_type)
setDT(md_calico_lung)
md_calico_lung[
  celltype_conversion[tissue == "Lung"],
  on = c("cell_type==cell_ontology_class"),
  cell_ontology_scdiffcom := i.cell_ontology_scdiffcom
  ]
sort(table(md_calico_lung$cell_ontology_scdiffcom), decreasing = TRUE)
sort(unique(md_calico_lung$cell_ontology_scdiffcom))
anyNA(md_calico_lung$cell_ontology_scdiffcom)
md_calico_lung[, cell_ontology_scdiffcom := ifelse(is.na(cell_ontology_scdiffcom),
                                                   cell_type,
                                                   cell_ontology_scdiffcom)]
anyNA(md_calico_lung$cell_ontology_scdiffcom)
md_calico_lung[
  celltype_conversion[tissue == "Lung"],
  on = c("subtype==cell_ontology_class"),
  subtype_ontology_scdiffcom := i.cell_ontology_scdiffcom
  ]
md_calico_lung[, subtype_ontology_scdiffcom := ifelse(is.na(subtype_ontology_scdiffcom),
                                                      subtype,
                                                      subtype_ontology_scdiffcom)]
identical(md_calico_lung$cell_ontology_scdiffcom, md_calico_lung$subtype_ontology_scdiffcom)
identical(as.character(md_calico_lung$cell_type), as.character(seurat_objects$calico_lung$cell_type))
seurat_objects$calico_lung$cell_ontology_scdiffcom <- md_calico_lung$cell_ontology_scdiffcom

#
md_calico_spleen <- seurat_objects$calico_spleen[[]]
md_calico_spleen$cell_type <- as.character(md_calico_spleen$cell_type)
setDT(md_calico_spleen)
md_calico_spleen[
  celltype_conversion[tissue == "Spleen"],
  on = c("cell_type==cell_ontology_class"),
  cell_ontology_scdiffcom := i.cell_ontology_scdiffcom
  ]
sort(table(md_calico_spleen$cell_ontology_scdiffcom), decreasing = TRUE)
sort(unique(md_calico_spleen$cell_ontology_scdiffcom))
anyNA(md_calico_spleen$cell_ontology_scdiffcom)
md_calico_spleen[, cell_ontology_scdiffcom := ifelse(is.na(cell_ontology_scdiffcom),
                                                     cell_type,
                                                     cell_ontology_scdiffcom)]
anyNA(md_calico_spleen$cell_ontology_scdiffcom)
md_calico_spleen[
  celltype_conversion[tissue == "Spleen"],
  on = c("subtype==cell_ontology_class"),
  subtype_ontology_scdiffcom := i.cell_ontology_scdiffcom
  ]
md_calico_spleen[, subtype_ontology_scdiffcom := ifelse(is.na(subtype_ontology_scdiffcom),
                                                        subtype,
                                                        subtype_ontology_scdiffcom)]
identical(md_calico_spleen$cell_ontology_scdiffcom, md_calico_spleen$subtype_ontology_scdiffcom)
identical(as.character(md_calico_spleen$cell_type), as.character(seurat_objects$calico_spleen$cell_type))
seurat_objects$calico_spleen$cell_ontology_scdiffcom <- md_calico_spleen$cell_ontology_scdiffcom

##
saveRDS(seurat_objects$tms_facs,
        "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_tms_facs.rds")
saveRDS(seurat_objects$tms_droplet,
        "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_tms_droplet.rds")
saveRDS(seurat_objects$calico_kidney,
        "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_calico_kidney.rds")
saveRDS(seurat_objects$calico_lung,
        "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_calico_lung.rds")
saveRDS(seurat_objects$calico_spleen,
        "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_scdiffcom_calico_spleen.rds")

## Save metadata to be used in other scripts not on the server ####

seurat_md <- lapply(seurat_objects, function(i) {
  i@meta.data
})

saveRDS(seurat_md, paste0(dir_data_analysis, "analysis_1_data_seurat_md.rds"))

## Save the spleen tissue for testing
tms_facs_spleen <- subset(seurat_objects$tms_facs, subset = tissue == "Spleen")
tms_droplet_spleen <- subset(seurat_objects$tms_droplet, subset = tissue == "Spleen")
tms_facs_liver <- subset(seurat_objects$tms_facs, subset = tissue == "Liver")
saveRDS(tms_facs_spleen,
        "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_testing_tms_facs_spleen.rds")
saveRDS(tms_droplet_spleen,
        "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_testing_tms_droplet_spleen.rds")
saveRDS(tms_facs_liver,
        "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_testing_tms_facs_liver.rds")


