####################################################
##
## Project: scAgeCom
##
## Last update - April 2022
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## Overlap scAgeCom results
## with PubMed, GenAge, CellAge, etc
##
####################################################
##

## Libraries ####

library(data.table)
library(scDiffCom)
library(biomaRt)
library(rentrez)
library(HAGRR)

## Load scDiffCom LRI and genes ####

LRI_mouse_dt <- scDiffCom::LRI_mouse$LRI_curated
LRI_mouse_genes <- unique(unlist(LRI_mouse_dt[, 2:6]))
LRI_mouse_genes <- LRI_mouse_genes[!is.na(LRI_mouse_genes)]

LRI_human_dt <- scDiffCom::LRI_human$LRI_curated
LRI_human_genes <- unique(unlist(LRI_human_dt[, 2:6]))
LRI_human_genes <- LRI_human_genes[!is.na(LRI_human_genes)]

## Orthology relationship ####

mart_mouse <- biomaRt::useMart(
  "ensembl",
  host = "https://nov2020.archive.ensembl.org",
  dataset = "mmusculus_gene_ensembl",
  verbose = TRUE
)
mart_human <- biomaRt::useMart(
  "ensembl",
  host = "https://nov2020.archive.ensembl.org",
  dataset = "hsapiens_gene_ensembl",
  verbose = TRUE
)

mouse_genes_info <- biomaRt::getBM(
  attributes = c(
    "mgi_symbol",
    "entrezgene_id",
    "ensembl_gene_id"
  ),
  filters = "mgi_symbol",
  mart = mart_mouse,
  values = LRI_mouse_genes
)
mouse_genes_ortho <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_confidence",
    "hsapiens_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  mart = mart_mouse,
  values = mouse_genes_info$ensembl_gene_id
)
human_genes_info <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "entrezgene_id",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  mart = mart_human,
  values = unique(mouse_genes_ortho$hsapiens_homolog_associated_gene_name)
)
colnames(human_genes_info) <- paste0(colnames(human_genes_info), "_human")

setDT(mouse_genes_info)
setDT(mouse_genes_ortho)
setDT(human_genes_info)

mouse_genes <- merge.data.table(
    mouse_genes_info,
    mouse_genes_ortho,
    by = "ensembl_gene_id",
    all.x = TRUE,
    all.y = TRUE
)
mouse_genes <- merge.data.table(
    mouse_genes,
    human_genes_info[, c(1,2)],
    by.x = "hsapiens_homolog_associated_gene_name",
    by.y = "hgnc_symbol_human",
    all.x = TRUE,
    all.y = TRUE
)
mouse_genes <- unique(mouse_genes)

anyNA(mouse_genes$mgi_symbol)
table(mouse_genes$mgi_symbol %in% LRI_mouse_genes)
table(LRI_mouse_genes %in% mouse_genes_info$mgi_symbol)

## Load gene2pubmed data ####

gene2pubmed <- fread("/workspace/lcyril_data/gene2pubmed/gene2pubmed")
setnames(
    gene2pubmed,
    old = colnames(gene2pubmed)[[1]],
    new = c("tax_id")
)
gene2pubmed_mouse <- gene2pubmed[tax_id == "10090"]
gene2pubmed_human <- gene2pubmed[tax_id == "9606"]

gene2pubmed_mouse[
    ,
    `:=` (count_GeneID = .N),
    by = "GeneID"
]
gene2pubmed_mouse[
    ,
    `:=` (count_PubMed_ID = .N),
    by = "PubMed_ID"
]

gene2pubmed_human[
    ,
    `:=` (count_GeneID = .N),
    by = "GeneID"
]
gene2pubmed_human[
    ,
    `:=` (count_PubMed_ID = .N),
    by = "PubMed_ID"
]

## Get aging pubmed ids #####

pmid_aging <- rentrez::entrez_search(
    "pubmed",
    paste0(
        "aging[TIAB] OR longevity[TIAB] OR senescence[TIAB]",
        " OR age-related[TIAB] OR dementia[TIAB]"
    ),
    retmax = 1000001
)
any(duplicated(pmid_aging$ids))

gene2pubmed_mouse_aging <- gene2pubmed_mouse[
    PubMed_ID %in% pmid_aging$ids &
    count_PubMed_ID <= 50
]
gene2pubmed_human_aging <- gene2pubmed_human[
    PubMed_ID %in% pmid_aging$ids &
    count_PubMed_ID <= 50
]

## Add number of pmid to LRI ####

LRI_mouse_pmid_aging <- lapply(
    LRI_mouse_genes,
    function(i) {
        gene2pubmed_mouse_aging[
            GeneID %in%
            mouse_genes_info[
                mgi_symbol == i
            ]$entrezgene_id
        ]$PubMed_ID
    }
)
names(LRI_mouse_pmid_aging) <- LRI_mouse_genes

LRI_combined_pmid_aging <- lapply(
    LRI_mouse_genes,
    function(i) {
        temp1 <- gene2pubmed_mouse_aging[
            GeneID %in%
            mouse_genes[
                mgi_symbol == i
            ]$entrezgene_id
        ]$PubMed_ID
        temp2 <- gene2pubmed_human_aging[
            GeneID %in%
            mouse_genes[
                mgi_symbol == i
            ]$entrezgene_id_human
        ]$PubMed_ID
        unique(c(temp1, temp2))
    }
)
names(LRI_combined_pmid_aging) <- LRI_mouse_genes

table(sapply(LRI_mouse_pmid_aging, length))
table(sapply(LRI_combined_pmid_aging, length))

LRI_combined_pmid_aging_n <- data.table(
    gene = names(sapply(LRI_combined_pmid_aging, length)),
    count = sapply(LRI_combined_pmid_aging, length)
)

LRI_mouse_dt[
    LRI_combined_pmid_aging_n,
    on = "LIGAND_1==gene",
    LIGAND_1_agepmid := i.count
]
LRI_mouse_dt[
    LRI_combined_pmid_aging_n,
    on = "LIGAND_2==gene",
    LIGAND_2_agepmid := i.count
]
LRI_mouse_dt[
    LRI_combined_pmid_aging_n,
    on = "RECEPTOR_1==gene",
    RECEPTOR_1_agepmid := i.count
]
LRI_mouse_dt[
    LRI_combined_pmid_aging_n,
    on = "RECEPTOR_2==gene",
    RECEPTOR_2_agepmid := i.count
]
LRI_mouse_dt[
    LRI_combined_pmid_aging_n,
    on = "RECEPTOR_3==gene",
    RECEPTOR_3_agepmid := i.count
]

pmid_colnames <- c(
    "LIGAND_1_agepmid", "LIGAND_2_agepmid",
    "RECEPTOR_1_agepmid", "RECEPTOR_2_agepmid", "RECEPTOR_3_agepmid"
)

## Load HAGR data ####

genage_human <- fread("/workspace/lcyril_data/hagr/genage_human.csv")
genage_mouse <- fread("/workspace/lcyril_data/hagr/genage_models.csv")
genage_mouse <- genage_mouse[
    organism == "Mus musculus"
]
genage_genes <- unique(
    c(
        mouse_genes[
          hsapiens_homolog_associated_gene_name %in% genage_human$symbol
        ]$mgi_symbol,
        mouse_genes[
          mgi_symbol %in% genage_mouse$symbol
        ]$mgi_symbol
    )
)

genage_microarray <- fread(
    "/workspace/lcyril_data/hagr/genage_microarray.csv"
)
microarray_genes <- unique(
    c(
        mouse_genes[
          hsapiens_homolog_associated_gene_name %in% 
          genage_microarray[Species == "human"]$Symbol
        ]$mgi_symbol,
        mouse_genes[
          mgi_symbol %in% genage_microarray[Species == "mouse"]$Symbol
        ]$mgi_symbol
    )
)

cellage <- fread("/workspace/lcyril_data/hagr/cellage_v2.csv")
cellage_genes <- unique(
    mouse_genes[
          hsapiens_homolog_associated_gene_name %in% 
          cellage$gene_name
        ]$mgi_symbol
)

hagr_genes <- unique(
    c(
        genage_genes,
        microarray_genes,
        cellage_genes
    )
)

## Annotate LRI with HAGR ####

hagr_colnames <- c(
    "HAGR_L1", "HAGR_L2", "HAGR_R1", "HAGR_R2", "HAGR_R3"
)
LRI_mouse_dt[
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

LRI_mouse_dt[
    ,
    HAGR_full_LRI := HAGR_L1 & HAGR_R1 &
    (is.na(HAGR_L2) | HAGR_L2) & (is.na(HAGR_R2) | HAGR_R2) &
    (is.na(HAGR_R3) | HAGR_R3)
]

LRI_mouse_dt[, .N, by = hagr_colnames][order(-N)][1:10]
sort(LRI_mouse_dt[HAGR_full_LRI == TRUE]$LRI)

LRI_mouse_dt[HAGR_full_LRI == TRUE][
    ,
     c("LIGAND_1_agepmid", "RECEPTOR_1_agepmid")
]

## Annotate ORA LRI with HAGR and PMID ####

data_4_tissue_specific_results <- readRDS(
  "data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds"
)

ORA_table <- data_4_tissue_specific_results$ORA_table

table(ORA_table$ORA_REGULATION)
table(ORA_table$ORA_CATEGORY)

ORA_LRI <- ORA_table[
    ORA_CATEGORY == "LRI" & ORA_REGULATION != "Not Over-represented",
    c("Dataset", "Tissue", "VALUE", "ORA_REGULATION")
]
ORA_LRI[
    LRI_mouse_dt,
    on = "VALUE==LRI",
    (pmid_colnames) := mget(paste0("i.", pmid_colnames))
]
ORA_LRI[
    LRI_mouse_dt,
    on = "VALUE==LRI",
    (hagr_colnames) := mget(paste0("i.", hagr_colnames))
]
ORA_LRI_N <- ORA_LRI[
    ,
    .N,
    by = eval(colnames(ORA_LRI)[3:14])
][order(-N)]

ORA_LIGAND <- ORA_table[
    ORA_CATEGORY == "LIGAND_COMPLEX" &
     ORA_REGULATION != "Not Over-represented",
    c("Dataset", "Tissue", "VALUE", "ORA_REGULATION")
]
ORA_LIGAND[
    LRI_mouse_dt,
    on = "VALUE==LIGAND_1",
    LIGAND_1_agepmid := i.LIGAND_1_agepmid
]
ORA_LIGAND[
    LRI_mouse_dt,
    on = "VALUE==LIGAND_1",
    HAGR_L1 := i.HAGR_L1
]
ORA_LIGAND_N <- ORA_LIGAND[
    ,
    .N,
    by = eval(colnames(ORA_LIGAND)[3:6])
][order(-N)]

ORA_RECEPTOR <- ORA_table[
    ORA_CATEGORY == "RECEPTOR_COMPLEX" &
     ORA_REGULATION != "Not Over-represented",
    c("Dataset", "Tissue", "VALUE", "ORA_REGULATION")
]
ORA_RECEPTOR[
    LRI_mouse_dt,
    on = "VALUE==RECEPTOR_1",
    RECEPTOR_1_agepmid := i.RECEPTOR_1_agepmid
]
ORA_RECEPTOR[
    LRI_mouse_dt,
    on = "VALUE==RECEPTOR_1",
    HAGR_R1 := i.HAGR_R1
]
ORA_RECEPTOR_N <- ORA_RECEPTOR[
    ,
    .N,
    by = eval(colnames(ORA_RECEPTOR)[3:6])
][order(-N)]

## Statistics of LRI ####

LRI_combined_pmid_aging_n[order(-count)][1:10]

table(
    LRI_combined_pmid_aging_n[order(-count)][1:10]$gene %in%
    hagr_genes
)


## New LRI of interest ####

ORA_LIGAND_N[
    ORA_REGULATION != "FLAT" &
    LIGAND_1_agepmid <= 5 &
    HAGR_L1 == FALSE &
    N >= 5
][1:20]
ORA_RECEPTOR_N[
    ORA_REGULATION != "FLAT" &
    RECEPTOR_1_agepmid <= 5 &
    HAGR_R1 == FALSE &
    N >= 5
][1:20]

## Egfl7 #####

ORA_LRI_N[
    grepl("Egfl", VALUE)
]

# Egfl7:Notch1/2/3 upregulated.
# Egfl7 not aging characterized
# Notch is linked to pubmed aging (thymocyte development, angiogenesis
# down implies less regeneration in muscle, senescence, etc, LRI_mouse_pmid_aging$Notch1). Notch
# is related to many processes.
# Here: angiogenesis (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3286203/) and
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3031397/)
# neural stem cell renewal (https://www.nature.com/articles/ncb1896)
# Egfl7 acts as a Notch antagonist (reduce target gene expression)/opposite as well
#/depent on the micorenviornment
# Overexpression of Egfl7 can lead to abnormal vessel patterning
# In NSC, Egfl7 reduces proliferation and trigger differentiation,
# here maybe a similar effect on MSC in adipose tissues

ccis_egfl7_up <- data_4_tissue_specific_results$CCI_table[
    LIGAND_1 == "Egfl7" &
    `Age Regulation` == "UP"
]

#appears up-regulaed in adipose tissues, brain, heart and mammary gland a lot
ccis_egfl7_up[, .N, by = c("Dataset", "Tissue")][order(-N)]
#secreted by endothelial cells
ccis_egfl7_up[, .N, by = c("Tissue", "Emitter Cell Type")][order(-N)]
#towards MSC in adipose, endocardial cells and fibroblast of cardiac tissues in heart, T cells
ccis_egfl7_up[, .N, by = c("Tissue", "Receiver Cell Type")][order(-N)][1:10]
#towards pericyte and oligodendrocyte precursor cells in Brain
ccis_egfl7_up[, .N, by = c("Tissue", "Receiver Cell Type")][order(-N)][11:20]

ccis_egfl7_up[, .N,
 by = c("Dataset", "Tissue", "Emitter Cell Type", "Receiver Cell Type")
][order(-N)][1:10]

## F13a1 ####

ORA_LRI_N[
    grepl("F13a1", VALUE)
]

#downregulated with integrins (https://pubmed.ncbi.nlm.nih.gov/10816592/)
#The interaction of α4β1 and α9β1 with FXIII could be biologically significant.
# Both α4β1 and α9β1 are highly expressed on specific populations of leukocytes,
# where they play a prominent role in transendothelial leukocyte migration.
# As a member of the coagulation cascade, FXIII is likely to be enriched at sites
# of vascular injury and inflammation, where its interaction with α4β1 and α9β1 
#could promote leukocyte extravasation.
LRI_mouse_dt[LIGAND_1 == "F13a1"]

#not well documented with aging (AMD papers)
LRI_combined_pmid_aging$F13a1

ccis_f13a1_down <- data_4_tissue_specific_results$CCI_table[
    LIGAND_1 == "F13a1" &
    `Age Regulation` == "DOWN"
]

#appears down-regulaed in Marrow, adipose tissue and other
ccis_f13a1_down[, .N, by = c("Dataset", "Tissue")][order(-N)]
#secreted by (pro-)monocyte/macrophages/myeloid principally
ccis_f13a1_down[, .N, by = c("Tissue", "Emitter Cell Type")][order(-N)]
#towards immune cell types
ccis_f13a1_down[, .N, by = c("Tissue", "Receiver Cell Type")][order(-N)]

#not well documented interaction with integrins

## Gpi1 ####

ORA_LRI_N[
    grepl("Gpi1", VALUE)
]

LRI_combined_pmid_aging$Gpi1
LRI_mouse_dt[LIGAND_1 == "Gpi1"]

## Tln1 (probably false positive as intracellular with Itgb1) ####

## Col6a2, H2-Q6 # to do ####

## Slpi Plscr1 ####
ORA_LRI_N[
    grepl("Slpi", VALUE)
]

#Slpi with Plscr1 not linked to aging

LRI_mouse_dt[LRI == "Slpi:Pls"]

ccis_slpi_up <- data_4_tissue_specific_results$CCI_table[
    LIGAND_1 == "Slpi" &
    `Age Regulation` == "UP"
]

## Rpsa
ORA_LRI_N[
    grepl("Rpsa", VALUE)
]

ccis_rpsa_up <- data_4_tissue_specific_results$CCI_table[
    RECEPTOR_1 == "Rpsa" &
    `Age Regulation` == "UP"
]

## Gpc1

## Tnfrsf21