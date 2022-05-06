####################################################
##
## Project: scAgeCom
##
## Last update - May 2022
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## scRNA-seq data preparation
##
####################################################
##

## Libraries ####

library(Seurat)
library(data.table)

## Specific paths ####

path_project <- "/workspace/postdoc_aging/projects/all_projects/"
path_scagecom_input <- paste0(
  path_project,
  "P1_scInterComAging/data_scAgeCom/input/"
)
path_scagecom_output <- paste0(
  path_project,
  "P1_scInterComAging/data_scAgeCom/output/"
)

## Dataset Overview ####

# Tabula Muris Senis scRNA-seq datasets have been retrieved from Amazon S3
# following the link given here https://github.com/czbiohub/tabula-muris-senis

# Calico scRNA-seq datasets have been downloaded from here
# https://mca.research.calicolabs.com/data

# There were both stored as h5ad scanpy files, that we converted to
# Seurat objects and stored on our lab server here:

path_seurat <- "/workspace/lcyril_data/scRNA_seq/seurat_processed/"

# The above folder contains both original and processed Seurat objects.
# The processed objects contains cell type names and metadata to be used for our
# intercellular communication analysis.
# They have been obtained according to the code below.

## Loading unprocessed datasets ####

paths_datased_unprocessed <- c(
  tms_facs = paste0(path_seurat, "seurat_tms_facs.rds"),
  tms_droplet = paste0(path_seurat, "seurat_tms_droplet.rds"),
  calico_kidney = paste0(path_seurat, "seurat_calico_kidney.rds"),
  calico_lung = paste0(path_seurat, "seurat_calico_lung.rds"),
  calico_spleen = paste0(path_seurat, "seurat_calico_spleen.rds")
)

# requires several GB of RAM
seurats_unprocessed <- list(
  tms_facs = readRDS(paths_datased_unprocessed[["tms_facs"]]),
  tms_droplet = readRDS(paths_datased_unprocessed[["tms_droplet"]]),
  calico_kidney = readRDS(paths_datased_unprocessed[["calico_kidney"]]),
  calico_lung = readRDS(paths_datased_unprocessed[["calico_lung"]]),
  calico_spleen = readRDS(paths_datased_unprocessed[["calico_spleen"]])
)

## Double check basic contents of Seurat unprocessed objects ####

# Available assays: "RNA"
lapply(seurats_unprocessed, function(x) {
  x@assays
})

# Content of "counts" slot: integers
lapply(seurats_unprocessed, function(x) {
  x$RNA@counts[1:5,1:5]
})

# Content of metadata
lapply(seurats_unprocessed, function(x) {
  colnames(x@meta.data)
})

## Load file that contain new cell-type annotations ####

dt_celltype_conversion <- fread(
  paste0(
    path_scagecom_input,
    "table_tms_cell_type_conversion.csv"
  )
)
dt_celltype_conversion <- dt_celltype_conversion[
  ,
  c(
    "Tissue",
    "Original_annotation",
    "Final_annotation",
    "Family_broad",
    "Abbreviation"
  )
]
setnames(
  dt_celltype_conversion,
  old = colnames(dt_celltype_conversion),
  new = c(
    "tissue",
    "cell_ontology_class",
    "cell_ontology_final",
    "cell_family",
    "cell_abbreviation"
  )
)

any(duplicated(dt_celltype_conversion[, c("tissue", "cell_ontology_class")]))

## Change FACS cell-type annotations ####

md_tms_facs <- copy(seurats_unprocessed$tms_facs[[]])
md_tms_facs$cell_ontology_class <- as.character(
  md_tms_facs$cell_ontology_class
)
setDT(md_tms_facs)
md_tms_facs[
  dt_celltype_conversion,
  on = c("tissue", "cell_ontology_class"),
  names(dt_celltype_conversion)[3:5] := mget(
    paste0("i.", names(dt_celltype_conversion)[3:5])
  )
]
sort(table(md_tms_facs$cell_ontology_final), decreasing = TRUE)
sort(unique(md_tms_facs$cell_ontology_final))
anyNA(md_tms_facs$cell_ontology_final)
md_tms_facs[
  ,
  cell_ontology_final := ifelse(
    is.na(cell_ontology_final),
    cell_ontology_class,
    cell_ontology_final
  )
]
sort(table(md_tms_facs$cell_ontology_final), decreasing = TRUE)
sort(unique(md_tms_facs$cell_ontology_final))
anyNA(md_tms_facs$cell_ontology_final)
identical(
  as.character(md_tms_facs$cell_ontology_class),
  as.character(seurats_unprocessed$tms_facs$cell_ontology_class)
)
seurats_unprocessed$tms_facs$cell_ontology_final <- md_tms_facs$cell_ontology_final
seurats_unprocessed$tms_facs$cell_family <- md_tms_facs$cell_family
seurats_unprocessed$tms_facs$cell_abbreviation <- md_tms_facs$cell_abbreviation

## Change Droplet cell-type annotations ####

md_tms_droplet <- copy(seurats_unprocessed$tms_droplet[[]])
md_tms_droplet$cell_ontology_class <- as.character(
  md_tms_droplet$cell_ontology_class
)
setDT(md_tms_droplet)
md_tms_droplet[
  dt_celltype_conversion,
  on = c("tissue", "cell_ontology_class"),
  names(dt_celltype_conversion)[3:5] := mget(
    paste0("i.", names(dt_celltype_conversion)[3:5])
  )
]

sort(table(md_tms_droplet$cell_ontology_final), decreasing = TRUE)
sort(unique(md_tms_droplet$cell_ontology_final))
anyNA(md_tms_droplet$cell_ontology_final)
md_tms_droplet[
  , cell_ontology_final := ifelse(
    is.na(cell_ontology_final),
    cell_ontology_class,
    cell_ontology_final
  )
]
sort(table(md_tms_droplet$cell_ontology_final), decreasing = TRUE)
sort(unique(md_tms_droplet$cell_ontology_final))
anyNA(md_tms_droplet$cell_ontology_final)
identical(
  as.character(md_tms_droplet$cell_ontology_class),
  as.character(seurats_unprocessed$tms_droplet$cell_ontology_class)
)
seurats_unprocessed$tms_droplet$cell_ontology_final <- md_tms_droplet$cell_ontology_final
seurats_unprocessed$tms_droplet$cell_family <- md_tms_droplet$cell_family
seurats_unprocessed$tms_droplet$cell_abbreviation <- md_tms_droplet$cell_abbreviation

## Change Calico Kidney cell-type annotations ####

md_calico_kidney <- copy(seurats_unprocessed$calico_kidney[[]])
md_calico_kidney$cell_type <- as.character(md_calico_kidney$cell_type)
setDT(md_calico_kidney)
md_calico_kidney[
  dt_celltype_conversion[tissue == "Kidney"],
  on = c("cell_type==cell_ontology_class"),
  names(dt_celltype_conversion)[3:5] := mget(
    paste0("i.", names(dt_celltype_conversion)[3:5])
  )
]
sort(table(md_calico_kidney$cell_ontology_final), decreasing = TRUE)
sort(unique(md_calico_kidney$cell_ontology_final))
anyNA(md_calico_kidney$cell_ontology_final)
md_calico_kidney[
  , cell_ontology_final := ifelse(
    is.na(cell_ontology_final),
    cell_type,
    cell_ontology_final
  )
]
anyNA(md_calico_kidney$cell_ontology_final)
md_calico_kidney[
  dt_celltype_conversion[tissue == "Kidney"],
  on = c("subtype==cell_ontology_class"),
  subtype_ontology_final := i.cell_ontology_final
]
md_calico_kidney[
  ,
  subtype_ontology_final := ifelse(
    is.na(subtype_ontology_final),
    subtype,
    subtype_ontology_final
  )
]
identical(md_calico_kidney$cell_ontology_final, md_calico_kidney$subtype_ontology_final)
identical(as.character(md_calico_kidney$cell_type), as.character(seurats_unprocessed$calico_kidney$cell_type))
seurats_unprocessed$calico_kidney$cell_ontology_final <- md_calico_kidney$cell_ontology_final
seurats_unprocessed$calico_kidney$cell_family <- md_calico_kidney$cell_family
seurats_unprocessed$calico_kidney$cell_abbreviation <- md_calico_kidney$cell_abbreviation

## Change Calico Lung cell-type annotations ####

md_calico_lung <- copy(seurats_unprocessed$calico_lung[[]])
md_calico_lung$cell_type <- as.character(md_calico_lung$cell_type)
setDT(md_calico_lung)
md_calico_lung[
  dt_celltype_conversion[tissue == "Lung"],
  on = c("cell_type==cell_ontology_class"),
  names(dt_celltype_conversion)[3:5] := mget(paste0("i.", names(dt_celltype_conversion)[3:5]))
]
sort(table(md_calico_lung$cell_ontology_final), decreasing = TRUE)
sort(unique(md_calico_lung$cell_ontology_final))
anyNA(md_calico_lung$cell_ontology_final)
md_calico_lung[
  , cell_ontology_final := ifelse(
    is.na(cell_ontology_final),
    cell_type,
    cell_ontology_final
  )
]
anyNA(md_calico_lung$cell_ontology_final)

md_calico_lung[
  dt_celltype_conversion[tissue == "Lung"],
  on = c("subtype==cell_ontology_class"),
  subtype_ontology_final := i.cell_ontology_final
]
md_calico_lung[
  ,
  subtype_ontology_final := ifelse(
    is.na(subtype_ontology_final),
    subtype,
    subtype_ontology_final
  )
]
identical(md_calico_lung$cell_ontology_final, md_calico_lung$subtype_ontology_final)
identical(as.character(md_calico_lung$cell_type), as.character(seurats_unprocessed$calico_lung$cell_type))
seurats_unprocessed$calico_lung$cell_ontology_final <- md_calico_lung$cell_ontology_final
seurats_unprocessed$calico_lung$cell_family <- md_calico_lung$cell_family
seurats_unprocessed$calico_lung$cell_abbreviation <- md_calico_lung$cell_abbreviation

## Change Calico Spleen cell-type annotations ####

md_calico_spleen <- copy(seurats_unprocessed$calico_spleen[[]])
md_calico_spleen$cell_type <- as.character(md_calico_spleen$cell_type)
setDT(md_calico_spleen)
md_calico_spleen[
  dt_celltype_conversion[tissue == "Spleen"],
  on = c("cell_type==cell_ontology_class"),
  names(dt_celltype_conversion)[3:5] := mget(paste0("i.", names(dt_celltype_conversion)[3:5]))
]

sort(table(md_calico_spleen$cell_ontology_final), decreasing = TRUE)
sort(unique(md_calico_spleen$cell_ontology_final))
anyNA(md_calico_spleen$cell_ontology_final)
md_calico_spleen[
  , cell_ontology_final := ifelse(
    is.na(cell_ontology_final),
    cell_type,
    cell_ontology_final
  )
]
anyNA(md_calico_spleen$cell_ontology_final)
md_calico_spleen[
  dt_celltype_conversion[tissue == "Spleen"],
  on = c("subtype==cell_ontology_class"),
  subtype_ontology_final := i.cell_ontology_final
]
md_calico_spleen[
  , subtype_ontology_final := ifelse(
    is.na(subtype_ontology_final),
    subtype,
    subtype_ontology_final
  )
]
identical(md_calico_spleen$cell_ontology_final, md_calico_spleen$subtype_ontology_final)
identical(as.character(md_calico_spleen$cell_type), as.character(seurats_unprocessed$calico_spleen$cell_type))
seurats_unprocessed$calico_spleen$cell_ontology_final <- md_calico_spleen$cell_ontology_final
seurats_unprocessed$calico_spleen$cell_family <- md_calico_spleen$cell_family
seurats_unprocessed$calico_spleen$cell_abbreviation <- md_calico_spleen$cell_abbreviation


## save new Seurat objects ####

# saveRDS(
#   seurats_unprocessed$tms_facs,
#   paste0(
#     path_seurat,
#     "seurat_final_tms_facs.rds"
#   )
# )
# saveRDS(
#   seurats_unprocessed$tms_droplet,
#   paste0(
#     path_seurat,
#     "seurat_final_tms_droplet.rds"
#   )
# )
# saveRDS(
#   seurats_unprocessed$calico_kidney,
#   paste0(
#     path_seurat,
#     "seurat_final_calico_kidney.rds"
#   )
# )
# saveRDS(
#   seurats_unprocessed$calico_lung,
#   paste0(
#     path_seurat,
#     "seurat_final_calico_lung.rds"
#   )
# )
# saveRDS(
#   seurats_unprocessed$calico_spleen,
#   paste0(
#     path_seurat,
#     "seurat_final_calico_spleen.rds"
#   )
# )

## Load processed datasets ####

paths_dataset_processed <- c(
  tms_facs = paste0(path_seurat, "seurat_final_tms_facs.rds"),
  tms_droplet = paste0(path_seurat, "seurat_final_tms_droplet.rds"),
  calico_kidney = paste0(path_seurat, "seurat_final_calico_kidney.rds"),
  calico_lung = paste0(path_seurat, "seurat_final_calico_lung.rds"),
  calico_spleen = paste0(path_seurat, "seurat_final_calico_spleen.rds")
)

seurats_processed <- list(
  tms_facs = readRDS(paths_dataset_processed[["tms_facs"]]),
  tms_droplet = readRDS(paths_dataset_processed[["tms_droplet"]]),
  calico_kidney = readRDS(paths_dataset_processed[["calico_kidney"]]),
  calico_lung = readRDS(paths_dataset_processed[["calico_lung"]]),
  calico_spleen = readRDS(paths_dataset_processed[["calico_spleen"]])
)

## Create metadata data.table ####

dt_metadata <- rbindlist(
  lapply(
    seurats_processed,
    function(i) {
      as.data.table(i[[]])
    }
  )
  ,
  idcol = "dataset",
  fill = TRUE
)

saveRDS(
  dt_metadata,
  paste0(
    path_scagecom_output,
    "a1_dt_metadata.rds"
  )
)

## Check which cell-type annotations are all in cell-ontology ####

onto_cell <- ontoProc::getCellOnto()

table(
  unique(dt_metadata$cell_ontology_final) %in% onto_cell$name
)
unique(dt_metadata$cell_ontology_final)[
  !unique(dt_metadata$cell_ontology_final) %in% onto_cell$name
]
