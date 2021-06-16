####################################################
##
## Project: scAgeCom
##
## Last update - June 2021
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## scRNA-seq data preparation
##
## Note: file intended to be used from our server
##
####################################################
##

## Libraries ####

library(Seurat)
library(data.table)

## Dataset Overview ####

# Tabula Muris Senis scRNA-seq datasets have been retrieved from Amazon S3 
# following the link given here https://github.com/czbiohub/tabula-muris-senis

# Calico scRNA-seq datasets have been downloaded from here
# https://mca.research.calicolabs.com/data

# There were both stored as h5ad scanpy files, that we converted to 
# Seurat objects and stored on our lab server here:

server_path <- "/home/nis/lcyril/work/lcyril_data/scRNA_seq/seurat_processed/"

dataset_unprocessed_path <- c(
  tms_facs = paste0(server_path, "seurat_tms_facs.rds"),
  tms_droplet = paste0(server_path, "seurat_tms_droplet.rds"),
  calico_kidney = paste0(server_path, "seurat_calico_kidney.rds"),
  calico_lung = paste0(server_path, "seurat_calico_lung.rds"),
  calico_spleen = paste0(server_path, "seurat_calico_spleen.rds")
)

## Loading unprocessed datasets ####

# requires several GB of RAM
seurat_unprocessed_objects <- list(
  tms_facs = readRDS(dataset_unprocessed_path[["tms_facs"]]),
  tms_droplet = readRDS(dataset_unprocessed_path[["tms_droplet"]]),
  calico_kidney = readRDS(dataset_unprocessed_path[["calico_kidney"]]),
  calico_lung = readRDS(dataset_unprocessed_path[["calico_lung"]]),
  calico_spleen = readRDS(dataset_unprocessed_path[["calico_spleen"]])
)

## Double check basic contents of Seurat unprocessed objects ####

# Available assays: "RNA"
lapply(seurat_unprocessed_objects, function(x) {
  x@assays
})

# Content of "counts" slot: integers
lapply(seurat_unprocessed_objects, function(x) {
  x$RNA@counts[1:5,1:5]
})

# Content of metadata
lapply(seurat_unprocessed_objects, function(x) {
  colnames(x@meta.data)
})

## Load file that contain new cell-type annotations ####

data_analysis_path <- "/home/nis/lcyril/postdoc_aging/projects/all_projects/P1_scInterComAging/data_scAgeCom/analysis/"

celltype_conversion <- read.csv(
  paste0(
    data_analysis_path,
    "scDiffCom_cell_types_clean.csv"
  ),
  stringsAsFactors = FALSE
)

celltype_conversion <- celltype_conversion[
  ,
  c(
    "Tissue",
    "Original_annotation",
    "Final_annotation",
    "Family_broad",
    "Abbreviation"
  )
]

colnames(celltype_conversion) <- c(
  "tissue",
  "cell_ontology_class",
  "cell_ontology_final",
  "cell_family",
  "cell_abbreviation"
)

setDT(celltype_conversion)

celltype_conversion <- unique(celltype_conversion)

celltype_check <- celltype_conversion[, c("tissue", "cell_ontology_class")]
any(duplicated(celltype_check))

## Change FACS cell-type annotations ####

md_tms_facs <- copy(seurat_unprocessed_objects$tms_facs[[]])
md_tms_facs$cell_ontology_class <- as.character(
  md_tms_facs$cell_ontology_class
)
setDT(md_tms_facs)
md_tms_facs[
  celltype_conversion,
  on = c("tissue", "cell_ontology_class"),
  names(celltype_conversion)[3:5] := mget(
    paste0("i.", names(celltype_conversion)[3:5])
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
  as.character(seurat_unprocessed_objects$tms_facs$cell_ontology_class)
)
seurat_unprocessed_objects$tms_facs$cell_ontology_final <- md_tms_facs$cell_ontology_final
seurat_unprocessed_objects$tms_facs$cell_family <- md_tms_facs$cell_family
seurat_unprocessed_objects$tms_facs$cell_abbreviation <- md_tms_facs$cell_abbreviation

## Change Droplet cell-type annotations ####

md_tms_droplet <- copy(seurat_unprocessed_objects$tms_droplet[[]])
md_tms_droplet$cell_ontology_class <- as.character(
  md_tms_droplet$cell_ontology_class
)
setDT(md_tms_droplet)
md_tms_droplet[
  celltype_conversion,
  on = c("tissue", "cell_ontology_class"),
  names(celltype_conversion)[3:5] := mget(
    paste0("i.", names(celltype_conversion)[3:5])
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
  as.character(seurat_unprocessed_objects$tms_droplet$cell_ontology_class)
)
seurat_unprocessed_objects$tms_droplet$cell_ontology_final <- md_tms_droplet$cell_ontology_final
seurat_unprocessed_objects$tms_droplet$cell_family <- md_tms_droplet$cell_family
seurat_unprocessed_objects$tms_droplet$cell_abbreviation <- md_tms_droplet$cell_abbreviation

## Change Calico Kidney cell-type annotations ####

md_calico_kidney <- copy(seurat_unprocessed_objects$calico_kidney[[]])
md_calico_kidney$cell_type <- as.character(md_calico_kidney$cell_type)
setDT(md_calico_kidney)
md_calico_kidney[
  celltype_conversion[tissue == "Kidney"],
  on = c("cell_type==cell_ontology_class"),
  names(celltype_conversion)[3:5] := mget(
    paste0("i.", names(celltype_conversion)[3:5])
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
  celltype_conversion[tissue == "Kidney"],
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
identical(as.character(md_calico_kidney$cell_type), as.character(seurat_unprocessed_objects$calico_kidney$cell_type))
seurat_unprocessed_objects$calico_kidney$cell_ontology_final <- md_calico_kidney$cell_ontology_final
seurat_unprocessed_objects$calico_kidney$cell_family <- md_calico_kidney$cell_family
seurat_unprocessed_objects$calico_kidney$cell_abbreviation <- md_calico_kidney$cell_abbreviation

## Change Calico Lung cell-type annotations ####

md_calico_lung <- copy(seurat_unprocessed_objects$calico_lung[[]])
md_calico_lung$cell_type <- as.character(md_calico_lung$cell_type)
setDT(md_calico_lung)
md_calico_lung[
  celltype_conversion[tissue == "Lung"],
  on = c("cell_type==cell_ontology_class"),
  names(celltype_conversion)[3:5] := mget(paste0("i.", names(celltype_conversion)[3:5]))
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
  celltype_conversion[tissue == "Lung"],
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
identical(as.character(md_calico_lung$cell_type), as.character(seurat_unprocessed_objects$calico_lung$cell_type))
seurat_unprocessed_objects$calico_lung$cell_ontology_final <- md_calico_lung$cell_ontology_final
seurat_unprocessed_objects$calico_lung$cell_family <- md_calico_lung$cell_family
seurat_unprocessed_objects$calico_lung$cell_abbreviation <- md_calico_lung$cell_abbreviation

## Change Calico Spleen cell-type annotations ####

md_calico_spleen <- copy(seurat_unprocessed_objects$calico_spleen[[]])
md_calico_spleen$cell_type <- as.character(md_calico_spleen$cell_type)
setDT(md_calico_spleen)
md_calico_spleen[
  celltype_conversion[tissue == "Spleen"],
  on = c("cell_type==cell_ontology_class"),
  names(celltype_conversion)[3:5] := mget(paste0("i.", names(celltype_conversion)[3:5]))
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
  celltype_conversion[tissue == "Spleen"],
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
identical(as.character(md_calico_spleen$cell_type), as.character(seurat_unprocessed_objects$calico_spleen$cell_type))
seurat_unprocessed_objects$calico_spleen$cell_ontology_final <- md_calico_spleen$cell_ontology_final
seurat_unprocessed_objects$calico_spleen$cell_family <- md_calico_spleen$cell_family
seurat_unprocessed_objects$calico_spleen$cell_abbreviation <- md_calico_spleen$cell_abbreviation


## save new Seurat objects ####

# saveRDS(
#   seurat_unprocessed_objects$tms_facs,
#   paste0(
#     server_path,
#     "seurat_final_tms_facs.rds"
#   )
# )
# saveRDS(
#   seurat_unprocessed_objects$tms_droplet,
#   paste0(
#     server_path,
#     "seurat_final_tms_droplet.rds"
#   )
# )
# saveRDS(
#   seurat_unprocessed_objects$calico_kidney,
#   paste0(
#     server_path,
#     "seurat_final_calico_kidney.rds"
#   )
# )
# saveRDS(
#   seurat_unprocessed_objects$calico_lung,
#   paste0(
#     server_path,
#     "seurat_final_calico_lung.rds"
#   )
# )
# saveRDS(
#   seurat_unprocessed_objects$calico_spleen,
#   paste0(
#     server_path,
#     "seurat_final_calico_spleen.rds"
#   )
# )

## Loading processed datasets

dataset_processed_path <- c(
  tms_facs = paste0(server_path, "seurat_final_tms_facs.rds"),
  tms_droplet = paste0(server_path, "seurat_final_tms_droplet.rds"),
  calico_kidney = paste0(server_path, "seurat_final_calico_kidney.rds"),
  calico_lung = paste0(server_path, "seurat_final_calico_lung.rds"),
  calico_spleen = paste0(server_path, "seurat_final_calico_spleen.rds")
)

seurat_processed_objects <- list(
  tms_facs = readRDS(dataset_processed_path[["tms_facs"]]),
  tms_droplet = readRDS(dataset_processed_path[["tms_droplet"]]),
  calico_kidney = readRDS(dataset_processed_path[["calico_kidney"]]),
  calico_lung = readRDS(dataset_processed_path[["calico_lung"]]),
  calico_spleen = readRDS(dataset_processed_path[["calico_spleen"]])
)

## Create meta.data object ####

seurat_md <- lapply(seurat_processed_objects, function(i) {
  i@meta.data
})
seurat_md <- lapply(
  seurat_md,
  setDT
)

#save it for possible later use
#saveRDS(seurat_md, paste0(data_analysis_path, "analysis_1_data_seurat_md.rds"))

## Check that scDiffCom cell-type annotations are all in cell-ontology ####

cl <- ontoProc::getCellOnto()
cl_names <- cl$name

all_scdiffcom_ontology <- unique(
  unlist(
    lapply(
      seurat_md,
      function(i) i$cell_ontology_final
    )
  )
)

all_scdiffcom_ontology %in% cl_names
all_scdiffcom_ontology[!all_scdiffcom_ontology %in% cl_names]
# there is some difference because of undetermined cell-types





