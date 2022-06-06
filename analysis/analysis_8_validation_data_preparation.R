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

## Load Tushaus secretomics brain data ####

tushaus_data <- fread(
  "data_scAgeCom_11_04_2022_processed/tushaus_2020_secretome.csv"
)
tushaus_data <- tushaus_data[PG.Genes != ""]

## Cell types to consider in TMS/scAgeCom ####

tushaus_cell_type <- data.table(
    tushaus = c("Astro", "Microglia", "Neurons", "Oligodendrocytes"),
    tms = c("astrocyte", "microglial cell", "neuron", "oligodendrocyte")
)

## Clean Tushaus gene/protein names ####

table(tushaus_data$PG.Genes %in% unique(unlist(seurat_genes)))
tushaus_data$PG.Genes[
  !tushaus_data$PG.Genes %in% unique(unlist(seurat_genes))
]

tushaus_data[, PG.Genes := ifelse(PG.Genes == "Acta2;Actg2", "Acta2", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Adgrb3", "Bai3", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Akr1b1", "Akr1b3", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Aoc1", "Abp1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Arf3;Arf1", "Arf1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "B3GNT2", "B3gnt2", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Cntnap5b;Cntnap5a", "Cntnap5a", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "C5", "Hc", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Ca3", "Car3", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Ces1;Ces1f", "Ces1f", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Dnase2", "Dnase2a", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Dsg1a;Dsg1b", "Dsg1b", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Ero1a", "Ero1l", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Galnt17", "Wbscr17", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Gpi", "Gpi1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Gstp1;Gstp2", "Gstp1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "H3f3c;Hist1h3a;Hist1h3b;H3f3a", "H3f3a", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Hist1h1c;Hist1h1e;Hist1h1d", "Hist1h1c", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Hist1h2bf;Hist1h2bm;Hist1h2bb;Hist1h2bh;Hist2h2bb;Hist1h2bc;Hist1h2bk;Hist1h2bp", "Hist1h2bf", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Hnrnph1;Hnrnph2", "Hnrnph1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Ica", "1300017J02Rik", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Large1", "Large", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Man1a1", "Man1a", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Pea15", "Pea15a", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Rala;Ralb", "Rala", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Rnaset2a;Rnaset2b", "Rnaset2a", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Selenbp1;Selenbp2", "Selenbp1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Tenm1", "Odz1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Tenm2", "Odz2", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Tf", "Trf", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "D1Pas1;Ddx3y;Ddx3x", "D1Pas1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Eef1b", "Eef1b2", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Gsto1;Gsto2", "Gsto1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Hist1h2ab;Hist1h2ac;Hist1h2ad;Hist1h2ae;Hist1h2ag;Hist1h2ai;Hist1h2an;Hist1h2ao;Hist1h2ap;Hist2h2ac;Hist2h2aa1;Hist3h2a;Hist1h2af;Hist1h2ah;Hist1h2ak;H2afj", "Hist1h2ab", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Hprt1", "Hprt", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Hspa1b;Hspa1a", "Hspa1b", PG.Genes)]
tushaus_data <- rbindlist(
  list(
    tushaus_data,
    tushaus_data[PG.Genes == "Hspa1b"]
  )
)
tushaus_data[PG.Genes == "Hspa1b", PG.Genes := c("Hspa1a", "Hspa1b")]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Krt72", "Krt72-ps", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Nme1;Nme2", "Nme1", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Ppp1ca;Ppp1cc", "Ppp1ca", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Rab6a;Rab6b;Rab39a", "Rab6a", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Rab7a", "Rab7", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Sprr2d;Sprr2h", "Sprr2d", PG.Genes)]
tushaus_data[, PG.Genes := ifelse(PG.Genes == "Ubb;Ubc;Rps27a;Uba52", "Rps27a", PG.Genes)]

table(tushaus_data$PG.Genes %in% unique(unlist(seurat_genes)))
tushaus_data$PG.Genes[
  !tushaus_data$PG.Genes %in% unique(unlist(seurat_genes))
]

#"Rps27a" %in% unique(unlist(seurat_genes))
#"Rps27a" %in% LRI_mouse_genes

## Annotate Tushaus data by number of detected protein/cell type ####

tushaus_data[
  ,
  paste0("is_detected_Astro_", 1:6) := lapply(
    .SD,
    function(i) {!is.nan(i)}
  ),
  .SDcols = paste0("log2_LFQ_Astro_", 1:6)
]
tushaus_data[
  ,
  nrep_Astro := rowSums(.SD),
  .SDcols = paste0("is_detected_Astro_", 1:6)
]

tushaus_data[
  ,
  paste0("is_detected_Microglia_", 1:6) := lapply(
    .SD,
    function(i) {!is.nan(i)}
  ),
  .SDcols = paste0("log2_LFQ_Microglia_", 1:6)
]
tushaus_data[
  ,
  nrep_Microglia := rowSums(.SD),
  .SDcols = paste0("is_detected_Microglia_", 1:6)
]

tushaus_data[
  ,
  paste0("is_detected_Neurons_", 1:6) := lapply(
    .SD,
    function(i) {!is.nan(i)}
  ),
  .SDcols = paste0("log2_LFQ_Neurons_", 1:6)
]
tushaus_data[
  ,
  nrep_Neurons := rowSums(.SD),
  .SDcols = paste0("is_detected_Neurons_", 1:6)
]

tushaus_data[
  ,
  paste0("is_detected_Oligodendrocytes_", 1:6) := lapply(
    .SD,
    function(i) {!is.nan(i)}
  ),
  .SDcols = paste0("log2_LFQ_Oligodendrocytes_", 1:6)
]
tushaus_data[
  ,
  nrep_Oligodendrocytes := rowSums(.SD),
  .SDcols = paste0("is_detected_Oligodendrocytes_", 1:6)
]

ftable(
  tushaus_data$nrep_Astro> 4,
  tushaus_data$nrep_Microglia > 4,
  tushaus_data$nrep_Neurons > 4,
  tushaus_data$nrep_Oligodendrocytes > 4
)

## Look at RNA expression in TMS facs ####

tms_facs <- readRDS(dataset_path[[1]])
tms_facs_brain <- subset(
    tms_facs,
    subset = tissue %in% c("Brain_Non-Myeloid", "Brain_Myeloid")
)
lapply(
    seq_len(nrow(tushaus_cell_type)),
    function(i) {
        temp_seurat <- subset(
            tms_facs_brain,
            subset = cell_ontology_final == tushaus_cell_type[i]$tms
        )
        expr <- rowMeans(
            expm1(
                GetAssayData(
                    temp_seurat
                )[
                    tushaus_data[
                        get(paste0("nrep_", tushaus_cell_type[i]$tushaus)) > 4 &
                        PG.Genes %in% rownames(temp_seurat)
                    ]$PG.Genes,
                ]
            )
        )
        dt <- data.table(
            gene = names(expr),
            avgExpr = expr
        )
        tushaus_data[
            dt,
            on = "PG.Genes==gene",
            paste0("tms_avgExpr_", tushaus_cell_type[i]$tushaus) := i.avgExpr
        ]
    }
)

ggplot(
    tushaus_data,
    aes(
        x = log2(tms_avgExpr_Microglia),
        y = log2Avg_Microglia
    )
) + geom_point() + geom_smooth(method="lm")

ggplot(
    tushaus_data,
    aes(
        x = log2(tms_avgExpr_Astro),
        y = log2Avg_Astro
    )
) + geom_point() + geom_smooth(method="lm")

ggplot(
    tushaus_data,
    aes(
        x = log2(tms_avgExpr_Neurons),
        y = log2Avg_Neurons
    )
) + geom_point() + geom_smooth(method="lm")

ggplot(
    tushaus_data,
    aes(
        x = log2(tms_avgExpr_Oligodendrocytes),
        y = log2Avg_Oligodendrocytes
    )
) + geom_point() + geom_smooth(method="lm")


## Annotate CCIs with Tushaus genes ####

cci_tushaus <- copy(CCI_table)

cci_tushaus[
    ,
    L1_th := ifelse(
        LIGAND_1 %in% tushaus_data$PG.Genes, TRUE, FALSE
    )
]
cci_tushaus[
    ,
    L2_th := ifelse(
        LIGAND_2 %in% tushaus_data$PG.Genes, TRUE, FALSE
    )
]
cci_tushaus[
    ,
    R1_th := ifelse(
        RECEPTOR_1 %in% tushaus_data$PG.Genes, TRUE, FALSE
    )
]
cci_tushaus[
    ,
    R2_th := ifelse(
        RECEPTOR_2 %in% tushaus_data$PG.Genes, TRUE, FALSE
    )
]
cci_tushaus[
    ,
    R3_th := ifelse(
        RECEPTOR_3 %in% tushaus_data$PG.Genes, TRUE, FALSE
    )
]

sort(table(cci_tushaus$LIGAND_1))
sort(table(cci_tushaus[L1_th == TRUE]$LIGAND_1))
sort(table(cci_tushaus[R1_th == TRUE]$RECEPTOR_1))

## ####

tushaus_spec_neuron <- tushaus_data[
    nrep_Neurons > 4 &
    nrep_Astro < 5 &
    nrep_Microglia < 5 & 
    nrep_Oligodendrocytes < 5
]$PG.Genes

cci_tushaus[L1_th == TRUE & LIGAND_1 %in% tushaus_spec_neuron]

test <- cci_tushaus[L1_th == TRUE & LIGAND_1 %in% tushaus_spec_neuron][
    ,
    sum(CCI_SCORE_YOUNG),
    by = c("Dataset", "Tissue", "EMITTER_CELLTYPE", "LIGAND_1")
][order(-V1)]

test <- test[Tissue == "Brain"]

test <- cci_tushaus[L1_th == TRUE][
    ,
    sum(CCI_SCORE_YOUNG),
    by = c("Dataset", "Tissue", "EMITTER_CELLTYPE", "LIGAND_1")
][order(-V1)][LIGAND_1 == "App"]

ggplot(
    test,
    aes(
        x = V1,
        y = Tissue
    )
) + geom_boxplot(
) + geom_jitter()

test[Tissue == "Lung"]

cci_tushaus[L1_th == TRUE]

cci_tushaus[, .N, by = c("Dataset", "R3_th")]

scd_brain_facs <- CCI_table[
  Tissue == "Brain"
]
scd_brain_facs_male <- CCI_table[
  Dataset == "TMS FACS (male)" & Tissue == "Brain"
]
scd_brain_facs_female <- CCI_table[
  Dataset == "TMS FACS (female)" & Tissue == "Brain"
]

##

## CCIs originating from Tushaus cells ####

cci_tushaus <- scd_brain_facs[
    EMITTER_CELLTYPE %in% tushaus_cell_type$tms
]

table(cci_tushaus$EMITTER_CELLTYPE)

scd_brain_facs_male_from_astro <- scd_brain_facs_male[
  EMITTER_CELLTYPE == "astrocyte"
]

scd_brain_facs_male_fromto_astro <- scd_brain_facs_male[
  EMITTER_CELLTYPE == "astrocyte" &
  RECEIVER_CELLTYPE == "astrocyte"
]


scd_liver_facs_male_from_nk <- CCI_table[
  Dataset == "TMS FACS (male)" & Tissue == "Liver" &
    EMITTER_CELLTYPE == "natural killer cell"
]

ligand_background <- unique(c(
  unique(LRI_mouse_dt$LIGAND_1),
  tushaus_data$PG.Genes
))

test <- sapply(
  unique(CCI_table$Dataset),
  function(dataset) {
    sapply(
      unique(CCI_table[Dataset == dataset]$Tissue),
      function(tissue) {
        sapply(
          unique(CCI_table[
            Dataset == dataset & Tissue == tissue
          ]$EMITTER_CELLTYPE),
          function(emitter) {
            ligs <- unique(
              CCI_table[
                Dataset == dataset &
                  Tissue == tissue &
                  EMITTER_CELLTYPE == emitter
              ]$LIGAND_1
            )
            t1 <- intersect(
              ligs,
              tushaus_data[nrep_Microglia > 4]$PG.Genes
            )
            t2 <- setdiff(
              ligs,
              tushaus_data[nrep_Microglia > 4]$PG.Genes
            )
            t3 <- setdiff(
              tushaus_data[nrep_Microglia > 4]$PG.Genes,
              ligs
            )
            fisher.test(
              matrix(
                c(
                  length(t1),
                  length(t2),
                  length(t3),
                  length(ligand_background) - length(t1) - length(t2) - length(t3)
                ),
                nrow = 2
              )
            )$estimate
            #sum(ligs %in% tushaus_data[nrep_Astrocytes > 4]$PG.Genes) /
            #  length(ligs)
            #table(ligs %in% tushaus_data[nrep_Astrocytes > 4]$PG.Genes)
          }
        )
      }
    )
  }
)

test$`TMS Droplet (male)`$Liver

test$`TMS FACS (male)`$Brain

hist(unlist(test), breaks = 50)

ggplot(
  scd_brain_facs_male_from_astro,
  aes(
    x = log10(CCI_SCORE_YOUNG),
    y = log10(CCI_SCORE_OLD)
  )
) + geom_point() + geom_abline()

scd_brain_facs_male_from_astro[, .N, by = c("LIGAND_1")][order(-N)]

out_scores_astro <- unique(
    scd_brain_facs_male_from_astro[, c("LIGAND_1", "L1_EXPRESSION_YOUNG")]
)

out_scores_astro[
  tushaus_data,
  on = "LIGAND_1==PG.Genes",
  secretome_score := i.log2Avg_Asstrocytes
]

table(
    out_scores_astro$LIGAND_1 %in% tushaus_data[nrep_Astrocytes > 4]$PG.Genes
)

sum(out_scores_astro$secretome_score > 0, na.rm = TRUE)

table(!is.na(out_scores_astro$secretome_score) & !is.nan(out_scores_astro$secretome_score))

64/(64+111)

table(!is.na(out_scores_liver_nk$secretome_score) & !is.nan(out_scores_liver_nk$secretome_score))

15/75


out_scores_astro$ligand

ligand_background <- unique(c(
  unique(LRI_mouse_dt$LIGAND_1),
  tushaus_data$PG.Genes
))

out_scores_astro$ligand %in% ligand_background
tushaus_data[nrep_Astrocytes > 4]$PG.Genes %in% ligand_background

t1 <- intersect(
  out_scores_astro$ligand,
  tushaus_data[nrep_Astrocytes > 4]$PG.Genes
)
t2 <- setdiff(
  out_scores_astro$ligand,
  tushaus_data[nrep_Astrocytes > 4]$PG.Genes
)
t3 <- setdiff(
  tushaus_data[nrep_Astrocytes > 4]$PG.Genes,
  out_scores_astro$ligand
)
ft <- fisher.test(
  matrix(
    c(
      length(t1),
      length(t2),
      length(t3),
      length(ligand_background) - length(t1) - length(t2) - length(t3)
    ),
    nrow = 2
  )
)
ft$estimate
fisher.test(
  matrix(
    c(
      100,
      10,
      13,
      1244
    ),
    nrow = 2
  )
)

ggplot(
  out_scores_astro,
  aes(
    x = log2(L1_EXPRESSION_YOUNG),
    y = secretome_score
  )
) + geom_point() + geom_smooth(method = "lm")

table(is.na(out_scores_astro$secretome_score))

sort(test, decreasing = TRUE)[1:10]

tushaus_data[nrep_Astrocytes > 4][, ][, c("PG.Genes", "log2Avg_Asstrocytes")]

table(scd_brain_facs_male$REGULATION)

LRI_mouse_dt[
  ,
  L1_tushaus := ifelse(
    LIGAND_1 %in% tushaus_data$PG.Gene,
    TRUE, FALSE
  )
]

table(LRI_mouse_dt$L1_tushaus)

CCI_table[Tissue == "Brain" & EMITTER_CELLTYPE == "astrocyte"][, .N, by = "LRI"][order(-N)][
  LRI_mouse_dt,
  on = "LRI",
  L1_tushaus := i.L1_tushaus
][1:20]



tushaus_data$PG.Genes[
  !tushaus_data$PG.Genes %in% LRI_mouse_genes
]

## Load SASP atlas data ####

sasp_atlas <- fread("data_scAgeCom_11_04_2022_processed/sasp_atlas.csv")

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

## SASP genes of interest ####

sasp_atlas_up <- sasp_atlas[
  log2_ratio >= log2(1.5) & q_value <= 0.05 & !is.na(mgi_symbol),
  c("cell_type", "treatment", "mgi_symbol")
]
sasp_atlas_up[
  ,
  condition :=do.call(paste, c(.SD, sep = "_")),
  .SDcols = c("cell_type", "treatment")
]
sasp_atlas_up[, .N, by = c("cell_type", "treatment")]


sasp_atlas_down <- sasp_atlas[
  log2_ratio <= -log2(1.5) & q_value <= 0.05 & !is.na(mgi_symbol),
  c("cell_type", "treatment", "mgi_symbol")
]
sasp_atlas_down[
  ,
  condition :=do.call(paste, c(.SD, sep = "_")),
  .SDcols = c("cell_type", "treatment")
]
sasp_atlas_down[, .N, by = c("cell_type", "treatment")]


## Annotate CCI with SASP  ####

CCI_table[
  sasp_atlas_up[cell_type == "fibroblast" & treatment == "xir_d10"],
  on = "LIGAND_1==mgi_symbol",
  sasp_up_fib_xir := i.condition
]
CCI_table$sasp_up_fib_xir
CCI_table[
  grepl("fibroblast", EMITTER_CELLTYPE),
  .N,
  by = c("sasp_up_fib_xir", "REGULATION", "Dataset")
]


unique(CCI_table[Tissue == "Skin"]$EMITTER_CELLTYPE)

CCI_table[Tissue == "Skin", .N, by = c("sasp_up_fib_xir", "REGULATION", "Dataset")]

test <- CCI_table[
  grepl("fibroblast", EMITTER_CELLTYPE) &
    sasp_up_fib_xir == "fibroblast_xir_d10" &
    REGULATION == "UP"
][, c("LRI", "ER_CELLTYPES", "Tissue")]

unique(CCI_table$EMITTER_CELLTYPE)[
  grepl("fibroblast", unique(CCI_table$EMITTER_CELLTYPE))
]

sasp_genes_up <- sasp_genes_up[!is.na(sasp_genes_up)]
anyNA(sasp_genes_up)

sasp_genes_epi_up <- unique(
  sasp_atlas[log2_ratio >= log2(1.5) & q_value <= 0.05 & cell_type == "epithelial"]$mgi_symbol
)
sasp_genes_epi_up <- sasp_genes_epi_up[!is.na(sasp_genes_epi_up)]

CCI_table

table(CCI_table[
  (LIGAND_1 %in% sasp_genes_up |
     LIGAND_2 %in% sasp_genes_up |
     RECEPTOR_1 %in% sasp_genes_up |
     RECEPTOR_2 %in% sasp_genes_up |
     RECEPTOR_3 %in% sasp_genes_up)
]$REGULATION)

10584 / (33574 + 132592 + 55751 + 10584)
18256 / (50492 + 214024 + 99750 + 18256)

fisher.test(matrix(c(10584, 221917, 7672, 142349), ncol = 2, byrow = TRUE))
(10584 / 221917) / (7672 / 142349)

table(CCI_table$REGULATION)

sort(
  table(CCI_table[
    (LIGAND_1 %in% sasp_genes_up |
       LIGAND_2 %in% sasp_genes_up |
       RECEPTOR_1 %in% sasp_genes_up |
       RECEPTOR_2 %in% sasp_genes_up |
       RECEPTOR_3 %in% sasp_genes_up) &
      REGULATION == "UP"
  ]$Tissue),
  decreasing = TRUE
)[1:8]


sort(
  table(CCI_table[
    (LIGAND_1 %in% sasp_genes_epi_up |
       LIGAND_2 %in% sasp_genes_epi_up |
       RECEPTOR_1 %in% sasp_genes_epi_up |
       RECEPTOR_2 %in% sasp_genes_epi_up |
       RECEPTOR_3 %in% sasp_genes_epi_up) &
      REGULATION == "UP" &
      Tissue == "Kidney"
  ]$LRI),
  decreasing = TRUE
)[1:8]

###########


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
