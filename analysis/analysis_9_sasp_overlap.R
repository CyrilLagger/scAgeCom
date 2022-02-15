####################################################
##
## Project: scAgeCom
##
## Last update - February 2022
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## overlap scAgeCom and SASP ATLAS
##
####################################################
##

## Libraries ####

library(data.table)
library(scDiffCom)

## Load SASP atlas data ####

sasp_atlas <- fread("../data_scAgeCom/sasp_atlas.csv")

## Build reference gene annotations ####

genes_mart_nov2020 <- biomaRt::useMart(
    "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "https://nov2020.archive.ensembl.org",
    verbose = TRUE
)
genes_ensembl_nov2020 <- biomaRt::getBM(
    attributes = c(
        "ensembl_gene_id",
        "hgnc_symbol",
        "hgnc_id"
    ),
    mart = genes_mart_nov2020,
)
setDT(genes_ensembl_nov2020)
genes_ensembl_nov2020[
    ,
    hgnc_symbol := ifelse(
        hgnc_symbol == "",
        NA,
        hgnc_symbol
    )
]
genes_ensembl_nov2020[
    ,
    hgnc_id := ifelse(
        hgnc_id == "",
        NA,
        hgnc_id
    )
]

genes_ortho_nov2020 <- biomaRt::getBM(
    attributes = c(
        "ensembl_gene_id",
        "mmusculus_homolog_associated_gene_name",
        "mmusculus_homolog_orthology_confidence",
        "mmusculus_homolog_orthology_type"
    ),
    filters = "ensembl_gene_id",
    mart = genes_mart_nov2020,
    value = genes_ensembl_nov2020$ensembl_gene_id
)
setDT(genes_ortho_nov2020)

genes_ortho_nov2020[
    ,
    mmusculus_homolog_associated_gene_name := ifelse(
        mmusculus_homolog_associated_gene_name == "",
        NA,
        mmusculus_homolog_associated_gene_name
    )
]
genes_ortho_nov2020[
    ,
    mmusculus_homolog_orthology_type := ifelse(
        mmusculus_homolog_orthology_type == "",
        NA,
        mmusculus_homolog_orthology_type
    )
]

## Check gene names ####

table(
    sasp_atlas$gene %in% genes_ensembl_nov2020$hgnc_symbol
)

## Add orthologs ####

sasp_atlas <- merge.data.table(
    sasp_atlas,
    genes_ortho_nov2020,
    by.x = "ensembl_id",
    by.y = "ensembl_gene_id",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
)
table(is.na(sasp_atlas$mmusculus_homolog_associated_gene_name))
table(sasp_atlas$mmusculus_homolog_orthology_confidence)
setnames(
    sasp_atlas,
    old = "mmusculus_homolog_associated_gene_name",
    new = "mgi_symbol"
)

##  ####

sasp_atlas[, c("gene", "mgi_symbol")]
sasp_atlas[, 1:4]

table(unique(sasp_atlas$mgi_symbol) %in% LRI_mouse$LRI_curated$LIGAND_1)

unique(sasp_atlas$mgi_symbol)[
    unique(sasp_atlas$mgi_symbol) %in% LRI_mouse$LRI_curated$RECEPTOR_3
]