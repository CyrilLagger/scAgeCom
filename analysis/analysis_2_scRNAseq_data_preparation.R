####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## scRNA-seq data preparation
##
####################################################
##

## Add libraries ####

library(Seurat)
library(org.Mm.eg.db)

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

paths_dataset_unprocessed <- c(
  tms_facs = paste0(path_seurat, "seurat_tms_facs.rds"),
  tms_droplet = paste0(path_seurat, "seurat_tms_droplet.rds"),
  calico_kidney = paste0(path_seurat, "seurat_calico_kidney_unfiltered.rds"),
  calico_lung = paste0(path_seurat, "seurat_calico_lung_unfiltered.rds"),
  calico_spleen = paste0(path_seurat, "seurat_calico_spleen_unfiltered.rds")
)

# requires several GB of RAM
seurats_unprocessed <- list(
  tms_facs = readRDS(paths_dataset_unprocessed[["tms_facs"]]),
  tms_droplet = readRDS(paths_dataset_unprocessed[["tms_droplet"]]),
  calico_kidney = readRDS(paths_dataset_unprocessed[["calico_kidney"]]),
  calico_lung = readRDS(paths_dataset_unprocessed[["calico_lung"]]),
  calico_spleen = readRDS(paths_dataset_unprocessed[["calico_spleen"]])
)

## Double check basic contents of Seurat unprocessed objects ####

# Available assays: "RNA"
lapply(seurats_unprocessed, function(x) {
  x@assays
})

# Content of "counts" slot: integers
lapply(seurats_unprocessed, function(x) {
  x$RNA@counts[1:5, 1:5]
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
seurats_unprocessed$tms_facs$cell_ontology_final <-
  md_tms_facs$cell_ontology_final
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
seurats_unprocessed$tms_droplet$cell_ontology_final <-
  md_tms_droplet$cell_ontology_final
seurats_unprocessed$tms_droplet$cell_family <- md_tms_droplet$cell_family
seurats_unprocessed$tms_droplet$cell_abbreviation <-
  md_tms_droplet$cell_abbreviation

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
identical(
  md_calico_kidney$cell_ontology_final,
  md_calico_kidney$subtype_ontology_final
)
identical(
  as.character(md_calico_kidney$cell_type),
  as.character(seurats_unprocessed$calico_kidney$cell_type)
)
seurats_unprocessed$calico_kidney$cell_ontology_final <-
  md_calico_kidney$cell_ontology_final
seurats_unprocessed$calico_kidney$cell_family <-
  md_calico_kidney$cell_family
seurats_unprocessed$calico_kidney$cell_abbreviation <-
  md_calico_kidney$cell_abbreviation

## Change Calico Lung cell-type annotations ####

md_calico_lung <- copy(seurats_unprocessed$calico_lung[[]])
md_calico_lung$cell_type <- as.character(md_calico_lung$cell_type)
setDT(md_calico_lung)
md_calico_lung[
  dt_celltype_conversion[tissue == "Lung"],
  on = c("cell_type==cell_ontology_class"),
  names(dt_celltype_conversion)[3:5] := mget(
    paste0("i.", names(dt_celltype_conversion)[3:5])
  )
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
identical(
  md_calico_lung$cell_ontology_final,
  md_calico_lung$subtype_ontology_final
)
identical(
  as.character(md_calico_lung$cell_type),
  as.character(seurats_unprocessed$calico_lung$cell_type)
)
seurats_unprocessed$calico_lung$cell_ontology_final <-
  md_calico_lung$cell_ontology_final
seurats_unprocessed$calico_lung$cell_family <-
  md_calico_lung$cell_family
seurats_unprocessed$calico_lung$cell_abbreviation <-
  md_calico_lung$cell_abbreviation

## Change Calico Spleen cell-type annotations ####

md_calico_spleen <- copy(seurats_unprocessed$calico_spleen[[]])
md_calico_spleen$cell_type <- as.character(md_calico_spleen$cell_type)
setDT(md_calico_spleen)
md_calico_spleen[
  dt_celltype_conversion[tissue == "Spleen"],
  on = c("cell_type==cell_ontology_class"),
  names(dt_celltype_conversion)[3:5] := mget(
    paste0("i.", names(dt_celltype_conversion)[3:5])
  )
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
identical(
  md_calico_spleen$cell_ontology_final,
  md_calico_spleen$subtype_ontology_final
)
identical(
  as.character(md_calico_spleen$cell_type),
  as.character(seurats_unprocessed$calico_spleen$cell_type)
)
seurats_unprocessed$calico_spleen$cell_ontology_final <-
  md_calico_spleen$cell_ontology_final
seurats_unprocessed$calico_spleen$cell_family <-
  md_calico_spleen$cell_family
seurats_unprocessed$calico_spleen$cell_abbreviation <-
  md_calico_spleen$cell_abbreviation

## save Processed Seurat objects ####

saveRDS(
  seurats_unprocessed$tms_facs,
  paste0(
    path_seurat,
    "seurat_final_tms_facs.rds"
  )
)
saveRDS(
  seurats_unprocessed$tms_droplet,
  paste0(
    path_seurat,
    "seurat_final_tms_droplet.rds"
  )
)
saveRDS(
  seurats_unprocessed$calico_kidney,
  paste0(
    path_seurat,
    "seurat_final_calico_kidney_unfiltered.rds"
  )
)
saveRDS(
  seurats_unprocessed$calico_lung,
  paste0(
    path_seurat,
    "seurat_final_calico_lung_unfiltered.rds"
  )
)
saveRDS(
  seurats_unprocessed$calico_spleen,
  paste0(
    path_seurat,
    "seurat_final_calico_spleen_unfiltered.rds"
  )
)

## Load processed datasets ####

paths_dataset_processed <- c(
  tms_facs = paste0(path_seurat, "seurat_final_tms_facs.rds"),
  tms_droplet = paste0(path_seurat, "seurat_final_tms_droplet.rds"),
  calico_kidney = paste0(
    path_seurat,
    "seurat_final_calico_kidney_unfiltered.rds"
  ),
  calico_lung = paste0(
    path_seurat,
    "seurat_final_calico_lung_unfiltered.rds"
  ),
  calico_spleen = paste0(
    path_seurat,
    "seurat_final_calico_spleen_unfiltered.rds"
  )
)

seurats_processed <- list(
  tms_facs = readRDS(paths_dataset_processed[["tms_facs"]]),
  tms_droplet = readRDS(paths_dataset_processed[["tms_droplet"]]),
  calico_kidney = readRDS(paths_dataset_processed[["calico_kidney"]]),
  calico_lung = readRDS(paths_dataset_processed[["calico_lung"]]),
  calico_spleen = readRDS(paths_dataset_processed[["calico_spleen"]])
)

## Mapping gene names to MGI ####

genes_seurat_unmapped <- lapply(
  seurats_processed,
  rownames
)

mart_seurats <- biomaRt::getBM(
  attributes = c(
    "mgi_symbol",
    "entrezgene_id",
    "ensembl_gene_id"
  ),
  filters = "mgi_symbol",
  mart = mart_db_mouse,
  values = unique(c(unlist(genes_seurat_unmapped), genes_lri_mouse))
)
setDT(mart_seurats)

columns(org.Mm.eg.db)
dt_org_mm <- select(
  org.Mm.eg.db,
  keys = keys(org.Mm.eg.db),
  columns = c("SYMBOL", "ENSEMBL", "ENTREZID", "MGI", "ALIAS")
)
setDT(dt_org_mm)

table(genes_lri_mouse %in% dt_org_mm$SYMBOL)
table(genes_lri_mouse %in% dt_org_mm$ALIAS)

table(genes_seurat_unmapped$tms_facs %in% dt_org_mm$ALIAS)

dt_genes_seurat <- rbindlist(
  l = lapply(
    genes_seurat_unmapped,
    function(i) {
      data.table(gene_orig = i)
    }
  ),
  idcol = "dataset"
)

dt_seurat_org <- merge.data.table(
  dt_genes_seurat,
  dt_org_mm[, c("SYMBOL", "ALIAS")],
  by.x = "gene_orig",
  by.y = "ALIAS",
  all.x = TRUE,
  all.y = TRUE
)
dt_seurat_org <- unique(dt_seurat_org[!is.na(dataset)])
dt_seurat_org[
  ,
  symbol_in_lri := SYMBOL %in% genes_lri_mouse
]
dt_seurat_org[
  ,
  n_symbol_by_orig := .N,
  by = c("gene_orig", "dataset")
]
dt_seurat_org[
  ,
  n_orig_by_symbol := .N,
  by = c("SYMBOL", "dataset")
]

fun_map_genes <- function(dt) {
  n_symbol_by_orig <- n_orig_by_symbol <- gene_orig <-
  SYMBOL <- gene_res <- NULL
  dt_one2one <- dt[n_symbol_by_orig == 1 & n_orig_by_symbol == 1]
  dt_one2one[, gene_res := SYMBOL]
  dt_na <- dt[is.na(SYMBOL)]
  dt_na[, gene_res := gene_orig]
  dt_res <- rbindlist(
    l = list(
      dt_one2one[, c("gene_orig", "gene_res")],
      dt_na[, c("gene_orig", "gene_res")]
    )
  )
  dt_temp <- dt[!gene_orig %in% dt_res$gene_orig]
  dt_sim <- dt_temp[gene_orig == SYMBOL]
  dt_sim[, gene_res := gene_orig]
  dt_res <- rbindlist(
    l = list(
      dt_res,
      dt_sim[, c("gene_orig", "gene_res")]
    )
  )
  dt_temp <- dt_temp[gene_orig != SYMBOL]
  dt_temp <- dt_temp[!gene_orig %in% dt_res$gene_orig]
  for (orig in unique(dt_temp$gene_orig)) {
    if (orig %in% dt_res$gene_orig | orig %in% dt_res$gene_res) {
        print(orig)
        stop("internal error")
      } else {
        symbs <- dt_temp[gene_orig == orig]$SYMBOL
        symbs <- symbs[!symbs %in% c(dt_res$gene_orig, dt_res$gene_res)]
        if (length(symbs) == 0) {
          dt_res <- rbindlist(
            l = list(
              dt_res,
              data.table(gene_orig = orig, gene_res = orig)
            )
          )
        } else if (length(symbs) == 1) {
          dt_res <- rbindlist(
            l = list(
              dt_res,
              data.table(gene_orig = orig, gene_res = symbs)
            )
          )
        } else {
          dt_temp_2 <- dt_temp[gene_orig == orig & SYMBOL %in% symbs][
            symbol_in_lri == TRUE
          ]
          if (nrow(dt_temp_2) == 0) {
            dt_res <- rbindlist(
              l = list(
                dt_res,
                data.table(gene_orig = orig, gene_res = symbs[[1]])
              )
            )
          } else {
            dt_res <- rbindlist(
              l = list(
                dt_res,
                dt_temp_2[1][, gene_res := SYMBOL][, c("gene_orig", "gene_res")]
              )
            )
          }
        }
      }
  }
  return(dt_res)
}

genes_seurat_mapped <- lapply(
  names(genes_seurat_unmapped),
  function(dat) {
    fun_map_genes(dt_seurat_org[dataset == dat])
  }
)
names(genes_seurat_mapped) <- names(genes_seurat_unmapped)

## Function inspired from Seurat.utils to rename genes ####

rename_seurat_genes <- function(
  obj,
  newnames) {
  RNA <- obj@assays$RNA
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
  } else {
    stop("Unequal gene sets: nrow(RNA) != nrow(newnames)")
  }
  obj@assays$RNA <- RNA
  return(obj)
}

seurats_processed_renamed <- list()

## Check and change gene names for TMS FACS ####

table(
  genes_seurat_unmapped$tms_facs %in% genes_seurat_mapped$tms_facs$gene_orig
)
table(
  genes_seurat_mapped$tms_facs$gene_orig %in% genes_seurat_unmapped$tms_facs
)
any(duplicated(genes_seurat_mapped$tms_facs$gene_orig))
any(duplicated(genes_seurat_mapped$tms_facs$gene_res))
anyNA(genes_seurat_mapped$tms_facs$gene_res)

table(genes_seurat_unmapped$tms_facs %in% genes_lri_mouse)
table(genes_seurat_mapped$tms_facs$gene_orig %in% genes_lri_mouse)
table(genes_seurat_mapped$tms_facs$gene_res %in% genes_lri_mouse)

genes_seurat_mapped$tms_facs <- genes_seurat_mapped$tms_facs[
  order(match(gene_orig, rownames(seurats_processed$tms_facs)))
]
identical(
  genes_seurat_mapped$tms_facs$gene_orig,
  rownames(seurats_processed$tms_facs)
)

seurats_processed_renamed$tms_facs <- rename_seurat_genes(
  seurats_processed$tms_facs,
  genes_seurat_mapped$tms_facs$gene_res
)

identical(
  rownames(seurats_processed_renamed$tms_facs),
  genes_seurat_mapped$tms_facs$gene_res
)

## Check and change gene names for TMS Droplet ####

table(
  genes_seurat_unmapped$tms_droplet %in%
   genes_seurat_mapped$tms_droplet$gene_orig
)
table(
  genes_seurat_mapped$tms_droplet$gene_orig %in%
   genes_seurat_unmapped$tms_droplet
)
any(duplicated(genes_seurat_mapped$tms_droplet$gene_orig))
any(duplicated(genes_seurat_mapped$tms_droplet$gene_res))
anyNA(genes_seurat_mapped$tms_droplet$gene_res)

table(genes_seurat_unmapped$tms_droplet %in% genes_lri_mouse)
table(genes_seurat_mapped$tms_droplet$gene_orig %in% genes_lri_mouse)
table(genes_seurat_mapped$tms_droplet$gene_res %in% genes_lri_mouse)

genes_seurat_mapped$tms_droplet <- genes_seurat_mapped$tms_droplet[
  order(match(gene_orig, rownames(seurats_processed$tms_droplet)))
]
identical(
  genes_seurat_mapped$tms_droplet$gene_orig,
  rownames(seurats_processed$tms_droplet)
)

seurats_processed_renamed$tms_droplet <- rename_seurat_genes(
  seurats_processed$tms_droplet,
  genes_seurat_mapped$tms_droplet$gene_res
)

identical(
  rownames(seurats_processed_renamed$tms_droplet),
  genes_seurat_mapped$tms_droplet$gene_res
)

## Check and change gene names for Calico kidney ####

table(
  genes_seurat_unmapped$calico_kidney %in%
  genes_seurat_mapped$calico_kidney$gene_orig
)
table(
  genes_seurat_mapped$calico_kidney$gene_orig %in%
  genes_seurat_unmapped$calico_kidney
)
any(duplicated(genes_seurat_mapped$calico_kidney$gene_orig))
any(duplicated(genes_seurat_mapped$calico_kidney$gene_res))
anyNA(genes_seurat_mapped$calico_kidney$gene_res)

table(genes_seurat_unmapped$calico_kidney %in% genes_lri_mouse)
table(genes_seurat_mapped$calico_kidney$gene_orig %in% genes_lri_mouse)
table(genes_seurat_mapped$calico_kidney$gene_res %in% genes_lri_mouse)

genes_seurat_mapped$calico_kidney <- genes_seurat_mapped$calico_kidney[
  order(match(gene_orig, rownames(seurats_processed$calico_kidney)))
]
identical(
  genes_seurat_mapped$calico_kidney$gene_orig,
  rownames(seurats_processed$calico_kidney)
)

seurats_processed_renamed$calico_kidney <- rename_seurat_genes(
  seurats_processed$calico_kidney,
  genes_seurat_mapped$calico_kidney$gene_res
)

identical(
  rownames(seurats_processed_renamed$calico_kidney),
  genes_seurat_mapped$calico_kidney$gene_res
)

## Check and change gene names for Calico lung ####

table(
  genes_seurat_unmapped$calico_lung %in%
  genes_seurat_mapped$calico_lung$gene_orig
)
table(
  genes_seurat_mapped$calico_lung$gene_orig %in%
  genes_seurat_unmapped$calico_lung
)
any(duplicated(genes_seurat_mapped$calico_lung$gene_orig))
any(duplicated(genes_seurat_mapped$calico_lung$gene_res))
anyNA(genes_seurat_mapped$calico_lung$gene_res)

table(genes_seurat_unmapped$calico_lung %in% genes_lri_mouse)
table(genes_seurat_mapped$calico_lung$gene_orig %in% genes_lri_mouse)
table(genes_seurat_mapped$calico_lung$gene_res %in% genes_lri_mouse)

genes_seurat_mapped$calico_lung <- genes_seurat_mapped$calico_lung[
  order(match(gene_orig, rownames(seurats_processed$calico_lung)))
]
identical(
  genes_seurat_mapped$calico_lung$gene_orig,
  rownames(seurats_processed$calico_lung)
)

seurats_processed_renamed$calico_lung <- rename_seurat_genes(
  seurats_processed$calico_lung,
  genes_seurat_mapped$calico_lung$gene_res
)

identical(
  rownames(seurats_processed_renamed$calico_lung),
  genes_seurat_mapped$calico_lung$gene_res
)

## Check and change gene names for Calico spleen ####

table(
  genes_seurat_unmapped$calico_spleen %in%
  genes_seurat_mapped$calico_spleen$gene_orig
)
table(
  genes_seurat_mapped$calico_spleen$gene_orig %in%
  genes_seurat_unmapped$calico_spleen
)
any(duplicated(genes_seurat_mapped$calico_spleen$gene_orig))
any(duplicated(genes_seurat_mapped$calico_spleen$gene_res))
anyNA(genes_seurat_mapped$calico_spleen$gene_res)

table(genes_seurat_unmapped$calico_spleen %in% genes_lri_mouse)
table(genes_seurat_mapped$calico_spleen$gene_orig %in% genes_lri_mouse)
table(genes_seurat_mapped$calico_spleen$gene_res %in% genes_lri_mouse)

genes_seurat_mapped$calico_spleen <- genes_seurat_mapped$calico_spleen[
  order(match(gene_orig, rownames(seurats_processed$calico_spleen)))
]
identical(
  genes_seurat_mapped$calico_spleen$gene_orig,
  rownames(seurats_processed$calico_spleen)
)

seurats_processed_renamed$calico_spleen <- rename_seurat_genes(
  seurats_processed$calico_spleen,
  genes_seurat_mapped$calico_spleen$gene_res
)

identical(
  rownames(seurats_processed_renamed$calico_spleen),
  genes_seurat_mapped$calico_spleen$gene_res
)

## save renamed Seurat objects ####

saveRDS(
  seurats_processed_renamed$tms_facs,
  paste0(
    path_seurat,
    "seurat_final_and_renamed_tms_facs.rds"
  )
)
saveRDS(
  seurats_processed_renamed$tms_droplet,
  paste0(
    path_seurat,
    "seurat_final_and_renamed_tms_droplet.rds"
  )
)
saveRDS(
  seurats_processed_renamed$calico_kidney,
  paste0(
    path_seurat,
    "seurat_final_and_renamed_calico_kidney_unfiltered.rds"
  )
)
saveRDS(
  seurats_processed_renamed$calico_lung,
  paste0(
    path_seurat,
    "seurat_final_and_renamed_calico_lung_unfiltered.rds"
  )
)
saveRDS(
  seurats_processed_renamed$calico_spleen,
  paste0(
    path_seurat,
    "seurat_final_and_renamed_calico_spleen_unfiltered.rds"
  )
)

## Load datasets to analyse ####

paths_dataset_analysis <- c(
  tms_facs = paste0(path_seurat, "seurat_final_and_renamed_tms_facs.rds"),
  tms_droplet = paste0(path_seurat, "seurat_final_and_renamed_tms_droplet.rds"),
  calico_kidney = paste0(
    path_seurat,
    "seurat_final_and_renamed_calico_kidney_unfiltered.rds"
  ),
  calico_lung = paste0(
    path_seurat,
    "seurat_final_and_renamed_calico_lung_unfiltered.rds"
  ),
  calico_spleen = paste0(
    path_seurat,
    "seurat_final_and_renamed_calico_spleen_unfiltered.rds"
  )
)

seurats_analysis <- list(
  tms_facs = readRDS(paths_dataset_analysis[["tms_facs"]]),
  tms_droplet = readRDS(paths_dataset_analysis[["tms_droplet"]]),
  calico_kidney = readRDS(paths_dataset_analysis[["calico_kidney"]]),
  calico_lung = readRDS(paths_dataset_analysis[["calico_lung"]]),
  calico_spleen = readRDS(paths_dataset_analysis[["calico_spleen"]])
)

## Create metadata data.table ####

dt_metadata <- rbindlist(
  lapply(
    seurats_analysis,
    function(i) {
      as.data.table(i[[]])
    }
  ),
  idcol = "dataset",
  fill = TRUE
)

## Check which cell-type annotations are all in cell-ontology ####

table(
  unique(dt_metadata$cell_ontology_final) %in% ontoProc::getCellOnto()$name
)
unique(dt_metadata$cell_ontology_final)[
  !unique(dt_metadata$cell_ontology_final) %in% ontoProc::getCellOnto()$name
]
