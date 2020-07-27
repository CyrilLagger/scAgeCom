####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - July 2020
##
## Check the scRNA-seq dataset preprocessing.
## Provide some useful information about tissue,
## cell-types, genes, LR-pair etc
##
####################################################
##

# Note: works on the server only due to large file-size

## Libraries ####

library(Seurat)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)

## Paths of the Seurat objects ####

seurat_path <- c(
  tms_facs = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_tms_facs.rds",
  tms_droplet = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_tms_droplet.rds",
  calico_kidney = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_kidney.rds",
  calico_lung = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_lung.rds",
  calico_spleen = "/home/nis/zabidi/work/lcyril_data/scRNA_seq/seurat_processed/seurat_calico_spleen.rds"
)

## Other data path ####
dir_data <- "../data_scAgeCom/"

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

min_genes_per_cell <- lapply(unique_genes, min)
head(sort(unique_genes$calico_kidney))

# Total number of counts per cell
total_counts <- lapply(seurat_objects, function(x) {
  Matrix::colSums(x$RNA@counts)
})

min_total_count_per_cell <- lapply(total_counts, min)
head(sort(total_counts$calico_kidney))

# Compare counts to previously stored values (some differnce due to early filtering)
mapply(function(x,y) {identical(x$n_genes,y)}, seurat_objects, unique_genes, SIMPLIFY = FALSE)
identical(seurat_objects$tms_droplet$n_genes, unique_genes$tms_droplet)

## Prepare data for age/tissue/cell-types comparison ####

# add age groups
anyNA(seurat_objects$tms_facs$age)
unique(seurat_objects$tms_facs$age)
seurat_objects$tms_facs$age_group <- ifelse(seurat_objects$tms_facs$age %in% c('1m', '3m'), 'young', 'old')
anyNA(seurat_objects$tms_droplet$age)
unique(seurat_objects$tms_droplet$age)
seurat_objects$tms_droplet$age_group <- ifelse(seurat_objects$tms_droplet$age %in% c('1m', '3m'), 'young', 'old')
seurat_objects$calico_kidney$age_group <- seurat_objects$calico_kidney$age
seurat_objects$calico_lung$age_group <- seurat_objects$calico_lung$age
seurat_objects$calico_spleen$age_group <- seurat_objects$calico_spleen$age

# add tissue-cell types
seurat_objects[1:2] <- lapply(seurat_objects[1:2], function(x) {
  x$tissue_cell_type <- paste(x$tissue, x$cell_ontology_class, sep = "_")
  return(x)
})
seurat_objects$calico_kidney$tissue <- "Kidney"
seurat_objects$calico_lung$tissue <- "Lung"
seurat_objects$calico_spleen$tissue <- "Spleen"
seurat_objects[3:5] <- lapply(seurat_objects[3:5], function(x) {
  x$tissue_cell_type <- paste(x$tissue, x$cell_type, sep = "_")
  return(x)
})

# change some categories to character
seurat_objects <- lapply(seurat_objects, function(x) {
  x$tissue <- as.character(x$tissue)
  x$age_group <- as.character(x$age_group)
  x$tissue_cell_type <- as.character(x$tissue_cell_type)
  return(x)
})
seurat_objects$tms_facs$sex <- as.character(seurat_objects$tms_facs$sex)
seurat_objects$tms_droplet$sex <- as.character(seurat_objects$tms_droplet$sex)

# consider tissues with more than 5 cells 
tissue_toKeep <- lapply(seurat_objects, function(obj) {
  tokeep <- apply(
    table(obj$tissue, obj$age_group) >= 5,
    MARGIN = 1,
    FUN = all
  )
  names(tokeep[tokeep])
})

seurat_objects_filtered <- mapply(
  FUN = function(x,y) {
    cells_keep <- colnames(x)[x$tissue %in% y]
    return(subset(x, cells = cells_keep))
  },
  seurat_objects,
  tissue_toKeep,
  SIMPLIFY = FALSE
)

tissue_sex_toKeep <- lapply(seurat_objects[c(1,2)], function(obj) {
  tokeep <- apply(
    table(obj$tissue, obj$sex) >= 5,
    MARGIN = 1,
    FUN = all
  )
  names(tokeep[tokeep])
})

# consider cell types with more than 5 cells per age
tissue_cell_type_toKeep <- lapply(seurat_objects_filtered, function(obj) {
  tokeep <- apply(
    table(obj$tissue_cell_type, obj$age_group) >= 5,
    MARGIN = 1,
    FUN = any
  )
  names(tokeep[tokeep])
})

seurat_objects_filtered <- mapply(
  FUN = function(x,y) {
    cells_keep <- colnames(x)[x$tissue_cell_type %in% y]
    return(subset(x, cells = cells_keep))
  },
  seurat_objects_filtered,
  tissue_cell_type_toKeep,
  SIMPLIFY = FALSE
)

# number of tissue-cell type after filtering
number_tissue_ct <- sapply(seurat_objects_filtered, function(x) {
  length(unique(x$tissue_cell_type))
})
number_tissue_ct <- c(number_tissue_ct[1:2], sum(number_tissue_ct[3:5]))

## Summary data.frame ####

seurat_summary <- data.frame(
  name = c("Tabula Muris Senis - FACS", "Tabula Muris Senis - Droplet", "Calico"),
  tissue_number = c(length(unique(seurat_objects_filtered$tms_facs$tissue)),
                   length(unique(seurat_objects_filtered$tms_droplet$tissue)),
                   3),
  cell_type_number = number_tissue_ct
)
seurat_summary<- transpose(seurat_summary, keep.names = "")
colnames(seurat_summary) <- seurat_summary[1,]
seurat_summary <- seurat_summary[-c(1),]
seurat_summary$name <- c(
  "Number of tissues",
  "Number of cell types"
)
colnames(seurat_summary)[[1]] <- ""
g_seurat_summary <- tableGrob(seurat_summary, rows = NULL)
grid.newpage()
grid.draw(g_seurat_summary)
ggsave(filename = paste0(dir_data, "analysis/analysis_1_plot_seurat_summary.png"),
       plot = g_seurat_summary, scale = 1.5)

## Prepare list of meta.data ####

seurat_md <- list(
  tms_facs = seurat_objects_filtered$tms_facs@meta.data,
  tms_droplet = seurat_objects_filtered$tms_droplet@meta.data,
  calico = rbindlist(
    l = list(seurat_objects_filtered$calico_kidney@meta.data,
             seurat_objects_filtered$calico_lung@meta.data,
             seurat_objects_filtered$calico_spleen@meta.data
    ),
    use.names = TRUE
  )
)

# order tissue alphabetically
seurat_md <- lapply(seurat_md, function(md) {
  md$tissue <- factor(md$tissue, levels = sort(unique(md$tissue), decreasing = TRUE))
  return(md)
})

## Plot number of cell-types ####

g_cells_per_tissue <- cowplot::plot_grid(
  plotlist = lapply(seurat_md, function(md) {
    ggplot(md, aes(x = tissue, fill = age_group)) + geom_bar(position = "dodge") +
      theme(text=element_text(size=28)) +
      scale_y_log10() +
      xlab("Tissue") +
      ylab("Number of cells") +
      labs(fill = "Age") +
      coord_flip()
  }),
  ncol = 2,
  labels = c("TMS FACS", "TMS Droplet", "Calico"),
  align = "v",
  rel_heights = c(2, 1.2)#,
  #label_x = 0, label_y = 0,
  #hjust = -0.5, vjust = -0.5
)
#ggsave(paste0(dir_data, "analysis/analysis_1_plot_cell_per_tissue.png"),
#       plot = g_cells_per_tissue, scale = 2.5)

# number of cell-types per tissues
celltype_distr <- lapply(seurat_md, function(md) {
  setDT(md)
  unique(md[,c("tissue", "tissue_cell_type")])
})

g_celltype_per_tissue <- cowplot::plot_grid(
  plotlist = lapply(celltype_distr, function(md) {
    ggplot(md, aes(x = tissue)) + geom_bar() +
      geom_text(stat='count', aes(label=..count..), hjust = -0.1) +
      theme(text=element_text(size=28)) +
      xlab("Tissue") +
      ylab("Number of cell-types") +
      coord_flip()
  }),
  ncol = 2,
  labels = c("TMS FACS", "TMS Droplet", "Calico"),
  align = "v",
  rel_heights = c(2, 1.2)#,
  #label_x = 0, label_y = 0,
  #hjust = -0.5, vjust = -0.5
)
ggsave(paste0(dir_data, "analysis/analysis_1_plot_celltypes_per_tissue.png"),
       plot = g_celltype_per_tissue, scale = 2.4)

## Comparison of common tissues ####

seurat_common_md <- lapply(seurat_md, function(x) {
  x <- x[x$tissue %in% c("Kidney", "Lung", "Spleen"),]
  x$tissue <- factor(x$tissue, levels = sort(unique(as.character(x$tissue))))
  return(x)
})

g_cells_per_celltype <- cowplot::plot_grid(
  plotlist = lapply(seurat_common_md, function(md) {
    ggplot(md, aes(x = tissue_cell_type, fill = age_group)) + geom_bar(position = "dodge") +
      theme(text=element_text(size=16),
            axis.text.y=element_blank()) +
      xlab("Cell type (Kidney, Lung, Spleen)") +
      scale_y_log10() +
      ylab("Number of cells") +
      labs(fill = "Age") +
      coord_flip()
  })
  ,
  ncol = 2,
  labels = c("TMS FACS", "TMS Droplet", "Calico"),
  align = "v",
  rel_heights = c(1.3, 1)#,
  #label_x = 0, label_y = 0,
  #hjust = -0.5, vjust = -0.5
)
ggsave(paste0(dir_data, "analysis/analysis_1_plot_cell_per_celltypes.png"),
       plot = g_cells_per_celltype, scale = 2.4)

kidney_celltypes <- list(
  tms_facs = unique(seurat_common_md$tms_facs[tissue == "Kidney",]$tissue_cell_type),
  tms_droplet = unique(seurat_common_md$tms_droplet[tissue == "Kidney",]$tissue_cell_type),
  calico = unique(seurat_common_md$calico[tissue == "Kidney",]$tissue_cell_type)
)
kidney_celltypes

lung_celltypes <- list(
  tms_facs = unique(seurat_common_md$tms_facs[tissue == "Lung",]$tissue_cell_type),
  tms_droplet = unique(seurat_common_md$tms_droplet[tissue == "Lung",]$tissue_cell_type),
  calico = unique(seurat_common_md$calico[tissue == "Lung",]$tissue_cell_type)
)
lung_celltypes

spleen_celltypes <- list(
  tms_facs = sort(unique(as.character(seurat_common_md$tms_facs[tissue == "Spleen",]$cell_ontology_class))),
  tms_droplet = sort(unique(as.character(seurat_common_md$tms_droplet[tissue == "Spleen",]$cell_ontology_class))),
  calico = sort(unique(as.character(seurat_common_md$calico[tissue == "Spleen",]$cell_type))),
  calico_sub = sort(unique(as.character(seurat_common_md$calico[tissue == "Spleen",]$subtype)))
)
spleen_celltypes

# spleen data.frame for presentation

spleen_ct_df <- data.frame(
  tms_facs = c(spleen_celltypes$tms_facs, rep("",6)),
  tms_droplet = spleen_celltypes$tms_droplet,
  calico = c(spleen_celltypes$calico, rep("", 8)),
  calico_sub = c(spleen_celltypes$calico_sub, rep("", 4))
)
colnames(spleen_ct_df) <- c("TMS FACS", "TMS Droplet", "Calico", "Calico (subtypes)")
g_spleen <- tableGrob(spleen_ct_df, rows = NULL)
grid.newpage()
grid.draw(g_spleen)

ggsave(paste0(dir_data, "analysis/analysis_1_plot_spleen_celltypes.png"),
       plot = g_spleen, scale = 1.5)
