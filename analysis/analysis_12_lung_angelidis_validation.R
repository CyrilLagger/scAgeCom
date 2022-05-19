
## Load libraries ####

library(Seurat)
library(data.table)
library(scDiffCom)
library(future)

## Load scRNA-seq data and metadata ####

load(
    "data_scAgeCom_11_04_2022_processed/GSE124872_raw_counts_single_cell.RData"
)

angel_md <- fread(
    "data_scAgeCom_11_04_2022_processed/GSE124872_Angelidis_2018_metadata.csv"
)
setDF(angel_md, rownames = angel_md$cell_name)

## Clean cell names #####

identical(
    rownames(angel_md),
    colnames(raw_counts)
)
identical(
    sub(".*:", "", rownames(angel_md)),
    sub(".*:", "", colnames(raw_counts))
)
unique(
    data.table(
        data = sub(":[^:]*$", "", colnames(raw_counts)),
        md = sub(":[^:]*$", "", rownames(angel_md))
    )
)

colnames(raw_counts) <- paste(
    sub(":[^:]*$", "", rownames(angel_md)),
    sub(".*:", "", colnames(raw_counts)),
    sep = ":"
)
identical(
    rownames(angel_md),
    colnames(raw_counts)
)

## Create Seurat object ####

angel_seurat <- CreateSeuratObject(
    counts = raw_counts,
    project = "LungAngelidis",
    meta.data = angel_md
)

angel_seurat
rownames(angel_seurat)
colnames(angel_seurat)
str(angel_seurat[[]])
angel_seurat@assays$RNA@counts

## Normalize Seurat object ####

angel_seurat <- NormalizeData(angel_seurat)
angel_seurat@assays$RNA@data

## Compare gene names to scDiffCom TODO #####

angel_genes <- data.table(
    gene = rownames(angel_seurat)
)

scd_genes <- unique(
    unlist(
        LRI_mouse$LRI_curated[, 2:6]
    )
)
scd_genes <- scd_genes[!is.na(scd_genes)]

table(scd_genes %in% rownames(angel_seurat))
sort(scd_genes[!scd_genes %in% rownames(angel_seurat)])

"R75078" %in% angel_genes$gene
"Adam2" %in% angel_genes

## Rename cell types for scDiffcom ####

angel_seurat$cell_type_scd <- ifelse(
    angel_seurat$celltype %in% c(
        "Ccl17+/Cd103-/Cd11b-_dendritic_cells",
        "Cd103+/Cd11b-_dendritic_cells",
        "CD209+/Cd11b+_dendritic_cells"
    ),
    "dendritic_cells",
    angel_seurat$celltype
)
angel_seurat$cell_type_scd <- ifelse(
    angel_seurat$celltype %in% c(
        "Cd4+_T_cells",
        "CD8+_T_cells",
        "Gamma-Delta_T_cells"
    ),
    "T_cells",
    angel_seurat$cell_type_scd
)
angel_seurat$cell_type_scd <- ifelse(
    angel_seurat$celltype %in% c(
        "classical_monocyte_(Ly6c2+)"
    ),
    "classical_monocytes",
    angel_seurat$cell_type_scd
)
angel_seurat$cell_type_scd <- ifelse(
    angel_seurat$celltype %in% c(
        "non-classical_monocyte_(Ly6c2-)"
    ),
    "non-classical_monocyte",
    angel_seurat$cell_type_scd
)
angel_seurat$cell_type_scd <- ifelse(
    angel_seurat$celltype %in% c(
        "Capillary_endothelial_cells",
        "lymphatic_endothelial_cells",
        "vascular_endothelial_cells",
        "Vcam1+_endothelial_cells"
    ),
    "endothelial_cells",
    angel_seurat$cell_type_scd
)
angel_seurat$cell_type_scd <- ifelse(
    angel_seurat$celltype %in% c(
        "Alveolar_macrophage",
        "Fn1+_macrophage",
        "Interstitial_macrophages"
    ),
    "macrophage",
    angel_seurat$cell_type_scd
)
angel_seurat$cell_type_scd <- ifelse(
    angel_seurat$celltype %in% c(
        "Interstitial_Fibroblast",
        "Lipofibroblast"
    ),
    "fibroblast",
    angel_seurat$cell_type_scd
)

## Remove cell types ####

ct_remove <- c(
    "low_quality_cells",
    "Mki67+_proliferating_cells",
    "red_blood_cells",
    "Megakaryocytes"
)
ct_keep <- setdiff(unique(angel_seurat$cell_type_scd), ct_remove)

angel_seurat_scd <- subset(
    angel_seurat,
    subset = cell_type_scd %in% ct_keep
)
angel_seurat_scd

## Check important metadata distribution #####

ftable(
    angel_seurat_scd$cell_type_scd,
    angel_seurat_scd$grouping
)


## Run scDiffCom ####

plan(multicore, workers = 24)
options(future.globals.maxSize = 15 * 1024^3)

angel_scd <- run_interaction_analysis(
    seurat_object = angel_seurat_scd,
    LRI_species = "mouse",
    seurat_celltype_id = "cell_type_scd",
    seurat_condition_id = list(
        column_name = "grouping",
        cond1_name = "3m",
        cond2_name = "24m"
    ),
    iterations = 1000
)
angel_scd

BuildNetwork(angel_scd)
PlotORA(
    angel_scd,
    category = "LRI",
    regulation = "UP",
    max_terms_show = 30
)
PlotORA(
    angel_scd,
    category = "LRI",
    regulation = "DOWN",
    max_terms_show = 30
)

PlotORA(
    angel_scd,
    category = "GO_TERMS",
    regulation = "UP",
    max_terms_show = 30
)

BuildShiny(angel_scd)
