####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## LRI data preparation and aging annotation
##
####################################################
##

## Intro ####
# This script loads the ligand-receptor database from scDiffCom and
# does some preprocessing such as adding aging-annotations to each gene

## Load libraries ####

library(scDiffCom)
library(data.table)
library(biomaRt)
library(rentrez)

## Define paths ####

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

# The collections of ligand-receptor interactions (for both human and
# mouse) are stored in the package scDiffCom. See e.g. the
# scDiffCom vignette or the scDiffCom file R/utils_LRI.R

## Load human/mouse LRI database and genes ####

dt_lri_human <- copy(scDiffCom::LRI_human$LRI_curated)
genes_lri_human <- unique(
  unlist(
    dt_lri_human[, 2:6]
  )
)
genes_lri_human <- genes_lri_human[!is.na(genes_lri_human)]

dt_lri_mouse <- copy(scDiffCom::LRI_mouse$LRI_curated)
genes_lri_mouse <- unique(
  unlist(
    dt_lri_mouse[, 2:6]
  )
)
genes_lri_mouse <- genes_lri_mouse[!is.na(genes_lri_mouse)]

lri_dbs <- sort(
  unique(
    unlist(
      strsplit(
        dt_lri_mouse$DATABASE,
        ";")
    )
  )
)

## Aging annotation strategy ####

# The mouse LRI database will be annotated by the number of
# pmid related to aging and each gene (including human orthologs),
# as well as if (or its ortholog) it belongs to HAGR.

## Orthology relationship ####

mart_db_mouse <- biomaRt::useMart(
  "ensembl",
  host = "https://nov2020.archive.ensembl.org",
  dataset = "mmusculus_gene_ensembl",
  verbose = TRUE
)
mart_db_human <- biomaRt::useMart(
  "ensembl",
  host = "https://nov2020.archive.ensembl.org",
  dataset = "hsapiens_gene_ensembl",
  verbose = TRUE
)

mart_lri_mouse <- biomaRt::getBM(
  attributes = c(
    "mgi_symbol",
    "entrezgene_id",
    "ensembl_gene_id"
  ),
  filters = "mgi_symbol",
  mart = mart_db_mouse,
  values = genes_lri_mouse
)
mart_ortho_mouse <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_confidence",
    "hsapiens_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  mart = mart_db_mouse,
  values = mart_lri_mouse$ensembl_gene_id
)
mart_lri_human <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "entrezgene_id",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  mart = mart_db_human,
  values = unique(mart_ortho_mouse$hsapiens_homolog_associated_gene_name)
)
colnames(mart_lri_human) <- paste0(colnames(mart_lri_human), "_human")

setDT(mart_lri_mouse)
setDT(mart_ortho_mouse)
setDT(mart_lri_human)

mart_lri_mouse_anno <- merge.data.table(
  mart_lri_mouse,
  mart_ortho_mouse,
  by = "ensembl_gene_id",
  all.x = TRUE,
  all.y = TRUE
)
mart_lri_mouse_anno <- merge.data.table(
  mart_lri_mouse_anno,
  mart_lri_human[, c(1, 2)],
  by.x = "hsapiens_homolog_associated_gene_name",
  by.y = "hgnc_symbol_human",
  all.x = TRUE,
  all.y = TRUE
)
mart_lri_mouse_anno <- unique(mart_lri_mouse_anno)

anyNA(mart_lri_mouse_anno$mgi_symbol)
table(mart_lri_mouse_anno$mgi_symbol %in% genes_lri_mouse)
table(genes_lri_mouse %in% mart_lri_mouse_anno$mgi_symbol)
genes_lri_mouse[!genes_lri_mouse %in% mart_lri_mouse_anno$mgi_symbol]

## Load gene2pubmed data ####

gene2pubmed <- fread(
  paste0(
    path_scagecom_input,
    "gene2pubmed"
  )
)
setnames(
  gene2pubmed,
  old = colnames(gene2pubmed)[[1]],
  new = c("tax_id")
)
gene2pubmed_mouse <- gene2pubmed[tax_id == "10090"]
gene2pubmed_human <- gene2pubmed[tax_id == "9606"]
gene2pubmed_mouse[
  ,
  `:=`(count_GeneID = .N),
  by = "GeneID"
]
gene2pubmed_mouse[
  ,
  `:=`(count_PubMed_ID = .N),
  by = "PubMed_ID"
]
gene2pubmed_human[
  ,
  `:=`(count_GeneID = .N),
  by = "GeneID"
]
gene2pubmed_human[
  ,
  `:=`(count_PubMed_ID = .N),
  by = "PubMed_ID"
]

## Get aging pubmed ids #####

pmid_aging <- rentrez::entrez_search(
  "pubmed",
  paste0(
    "aging[TIAB] OR longevity[TIAB] OR senescence[TIAB]",
    " OR age-related[TIAB] OR dementia[TIAB] OR alzheimer[TIAB]",
    " OR parkinson[TIAB] OR atherosclerosis[TIAB] OR stroke[TIAB]",
    " OR arthritis[TIAB] OR osteoporosis[TIAB] OR cataract[TIAB]"
  ),
  retmax = 10000001
)
any(duplicated(pmid_aging$ids))

#if not working because of new NCBI retmax
pmid_aging_ids <- fread(
  paste0(
    path_scagecom_input,
    "aging_pmid_ncbi.txt"
  )
)

gene2pubmed_mouse_aging <- gene2pubmed_mouse[
  PubMed_ID %in% pmid_aging_ids$V1 &
    count_PubMed_ID <= 50
]
gene2pubmed_human_aging <- gene2pubmed_human[
  PubMed_ID %in% pmid_aging_ids$V1 &
    count_PubMed_ID <= 50
]

## Add number of pmids to LRIs ####

pmid_aging_lri <- lapply(
  genes_lri_mouse,
  function(i) {
    temp1 <- gene2pubmed_mouse_aging[
      GeneID %in%
        mart_lri_mouse_anno[
          mgi_symbol == i
        ]$entrezgene_id
    ]$PubMed_ID
    temp2 <- gene2pubmed_human_aging[
      GeneID %in%
        mart_lri_mouse_anno[
          mgi_symbol == i
        ]$entrezgene_id_human
    ]$PubMed_ID
    unique(c(temp1, temp2))
  }
)
names(pmid_aging_lri) <- genes_lri_mouse

table(sapply(pmid_aging_lri, length))

pmid_aging_lri_n <- data.table(
  gene = names(sapply(pmid_aging_lri, length)),
  count = sapply(pmid_aging_lri, length)
)

dt_lri_mouse[
  pmid_aging_lri_n,
  on = "LIGAND_1==gene",
  L1_N_agepmid := i.count
]
dt_lri_mouse[
  pmid_aging_lri_n,
  on = "LIGAND_2==gene",
  L2_N_agepmid := i.count
]
dt_lri_mouse[
  pmid_aging_lri_n,
  on = "RECEPTOR_1==gene",
  R1_N_agepmid := i.count
]
dt_lri_mouse[
  pmid_aging_lri_n,
  on = "RECEPTOR_2==gene",
  R2_N_agepmid := i.count
]
dt_lri_mouse[
  pmid_aging_lri_n,
  on = "RECEPTOR_3==gene",
  R3_N_agepmid := i.count
]

pmid_colnames <- c(
  "L1_N_agepmid", "L2_N_agepmid",
  "R1_N_agepmid", "R2_N_agepmid", "R3_N_agepmid"
)

## Load HAGR data ####

hagr_genage_human <- fread(
  paste0(
    path_scagecom_input,
    "hagr/genage_human.csv"
  )
)
hagr_genage_mouse <- fread(
  paste0(
    path_scagecom_input,
    "hagr/genage_models.csv"
  )
)
hagr_genage_mouse <- hagr_genage_mouse[
  organism == "Mus musculus"
]
hagr_genage_genes <- unique(
  c(
    mart_lri_mouse_anno[
      hsapiens_homolog_associated_gene_name %in% hagr_genage_human$symbol
    ]$mgi_symbol,
    mart_lri_mouse_anno[
      mgi_symbol %in% hagr_genage_mouse$symbol
    ]$mgi_symbol
  )
)

hagr_microarray <- fread(
  paste0(
    path_scagecom_input,
    "hagr/genage_microarray.csv"
  )
)
hagr_microarray_genes <- unique(
  c(
    mart_lri_mouse_anno[
      hsapiens_homolog_associated_gene_name %in%
        hagr_microarray[Species == "human"]$Symbol
    ]$mgi_symbol,
    mart_lri_mouse_anno[
      mgi_symbol %in% hagr_microarray[Species == "mouse"]$Symbol
    ]$mgi_symbol
  )
)

hagr_cellage <- fread(
  paste0(
    path_scagecom_input,
    "hagr/cellage_v2.csv"
  )
)
hagr_cellage_genes <- unique(
  mart_lri_mouse_anno[
    hsapiens_homolog_associated_gene_name %in%
      hagr_cellage$gene_name
  ]$mgi_symbol
)

hagr_longevity_map <- fread(
  paste0(
    path_scagecom_input,
    "hagr/longevity_map.csv"
  )
)
hagr_longevity_map <- hagr_longevity_map[Association == "significant"]
hagr_longevity_map_genes <- unique(
  mart_lri_mouse_anno[
    hsapiens_homolog_associated_gene_name %in%
    unique(unlist(strsplit(hagr_longevity_map[["Gene(s)"]], ",")))
  ]$mgi_symbol
)

hagr_genes <- unique(
  c(
    hagr_genage_genes,
    hagr_microarray_genes,
    hagr_cellage_genes,
    hagr_longevity_map_genes
  )
)

## Annotate LRI with HAGR ####

hagr_colnames <- c(
  "HAGR_L1", "HAGR_L2", "HAGR_R1", "HAGR_R2", "HAGR_R3"
)
dt_lri_mouse[
  ,
  (hagr_colnames) := lapply(
    .SD,
    function(i) {
      ifelse(
        is.na(i),
        NA,
        ifelse(
          i %in% hagr_genes, TRUE, FALSE
        )
      )
    }
  ),
  .SDcols = c(
    "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
  )
]
