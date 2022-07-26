####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## Process secretomics datasets used for validation
##
####################################################
##

## Pancreas (human) secretomics (Li et al. 2022) ####

# read spreadsheet
val_pancreas <- openxlsx::read.xlsx(
  paste0(
    path_scagecom_input,
    "Li_pancreas_secretome.xlsx"
  ),
  sheet = "S4",
  startRow = 5
)
setDT(val_pancreas)

val_pancreas_diff <- openxlsx::read.xlsx(
  paste0(
    path_scagecom_input,
    "Li_pancreas_secretome.xlsx"
  ),
  sheet = "S6",
  startRow = 6
)
setDT(val_pancreas_diff)

anyNA(val_pancreas$Accession)
any(duplicated(val_pancreas$Accession))
anyNA(val_pancreas_diff$Protein.accession)
any(duplicated(val_pancreas_diff$Protein.accession))

# extract genes of proteins detected
val_pancreas_genes <- val_pancreas[
  ,
  c("Accession", "Gene.Symbol", "Gene.ID")
]
val_pancreas_genes <- merge.data.table(
  x = val_pancreas_genes[
    ,
    lapply(.SD, function(x) unlist(tstrsplit(x, "; ", fixed = TRUE))),
    by = "Accession",
    .SDcols = "Gene.Symbol"
  ],
  y = val_pancreas_genes[
    ,
    lapply(.SD, function(x) unlist(tstrsplit(x, "; ", fixed = TRUE))),
    by = "Accession",
    .SDcols = "Gene.ID"
  ],
  by = "Accession",
  all = TRUE
)

# look for homologs
mart_pancreas_human_from_ens <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "ensembl_gene_id",
  mart = mart_db_human,
  values = sort(unique(val_pancreas_genes$Gene.Symbol))
)
setDT(mart_pancreas_human_from_ens)
mart_pancreas_human_from_symb <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  mart = mart_db_human,
  values = sort(unique(val_pancreas_genes$Gene.ID))
)
setDT(mart_pancreas_human_from_symb)

val_pancreas_genes <- merge.data.table(
  x = val_pancreas_genes,
  y = mart_pancreas_human_from_ens,
  by.x = "Gene.Symbol",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)
val_pancreas_genes <- merge.data.table(
  x = val_pancreas_genes,
  y = mart_pancreas_human_from_symb,
  by.x = "Gene.ID",
  by.y = "hgnc_symbol",
  all.x = TRUE
)

mart_pancreas_ortho <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "mmusculus_homolog_associated_gene_name",
    "mmusculus_homolog_orthology_confidence",
    "mmusculus_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  mart = mart_db_human,
  values = sort(
    unique(
      c(
        val_pancreas_genes$Gene.Symbol,
        val_pancreas_genes$ensembl_gene_id
      )
    )
  )
)
setDT(mart_pancreas_ortho)
anyNA(mart_pancreas_ortho$mmusculus_homolog_associated_gene_name)
mart_pancreas_ortho <- mart_pancreas_ortho[
  mmusculus_homolog_associated_gene_name != ""
]
val_pancreas_genes <- merge.data.table(
  val_pancreas_genes,
  mart_pancreas_ortho[, 1:2],
  by.x = "Gene.Symbol",
  by.y = "ensembl_gene_id",
  all = TRUE
)
setnames(
  val_pancreas_genes,
  old = "mmusculus_homolog_associated_gene_name",
  new = "mmusculus_homolog_1"
)
val_pancreas_genes <- merge.data.table(
  val_pancreas_genes,
  mart_pancreas_ortho[, 1:2],
  by.x = "ensembl_gene_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)
setnames(
  val_pancreas_genes,
  old = "mmusculus_homolog_associated_gene_name",
  new = "mmusculus_homolog_2"
)

# select genes that are in LRI
val_pancreas_lr <- sort(
  unique(
    c(
      val_pancreas_genes$mmusculus_homolog_1,
      val_pancreas_genes$mmusculus_homolog_2
    )
  )
)
val_pancreas_lr <- val_pancreas_lr[val_pancreas_lr %in% genes_lri_mouse]

# select genes that are in LRI only in hpde
val_pancreas_lr_hpde <- sort(
  unique(
    c(
      val_pancreas_genes[
        Accession %in% val_pancreas[
          Found.in.HPDE.replicate.1 != "Not Found" |
          Found.in.HPDE.replicate.2 != "Not Found"
        ]$Accession
      ]$mmusculus_homolog_1,
      val_pancreas_genes[
        Accession %in% val_pancreas[
          Found.in.HPDE.replicate.1 != "Not Found" |
          Found.in.HPDE.replicate.2 != "Not Found"
        ]$Accession
      ]$mmusculus_homolog_2
    )
  )
)
val_pancreas_lr_hpde <- val_pancreas_lr_hpde[
  val_pancreas_lr_hpde %in% genes_lri_mouse
]

# select genes that are differentially expressed with cancer
val_pancreas_lr_can_up <- sort(
  unique(
    c(
      val_pancreas_genes[
        Accession %in%
          val_pancreas_diff[log2.change > 0]$Protein.accession
      ]$mmusculus_homolog_1,
      val_pancreas_genes[
        Accession %in%
          val_pancreas_diff[log2.change > 0]$Protein.accession
      ]$mmusculus_homolog_2
    )
  )
)
val_pancreas_lr_can_up <- val_pancreas_lr_can_up[
  val_pancreas_lr_can_up %in% genes_lri_mouse
]

val_pancreas_lr_can_down <- sort(
  unique(
    c(
      val_pancreas_genes[
        Accession %in%
          val_pancreas_diff[log2.change < 0]$Protein.accession
      ]$mmusculus_homolog_1,
      val_pancreas_genes[
        Accession %in%
          val_pancreas_diff[log2.change < 0]$Protein.accession
      ]$mmusculus_homolog_2
    )
  )
)
val_pancreas_lr_can_down <- val_pancreas_lr_can_down[
  val_pancreas_lr_can_down %in% genes_lri_mouse
]

## Cardiomyocytes (rat) proteomics (Kuhn et al. 2020) ####

# read spreadsheet
val_cardio <- list(
  openxlsx::read.xlsx(
    paste0(
      path_scagecom_input,
      "kuhn_cardio_secretome.xlsx"
    ),
    sheet = "hypoxia_12",
    startRow = 1
  ),
  openxlsx::read.xlsx(
    paste0(
      path_scagecom_input,
      "kuhn_cardio_secretome.xlsx"
    ),
    sheet = "hypoxia_24",
    startRow = 1
  ),
  openxlsx::read.xlsx(
    paste0(
      path_scagecom_input,
      "kuhn_cardio_secretome.xlsx"
    ),
    sheet = "hypoxia_24_30",
    startRow = 1
  )
)
val_cardio <- lapply(
  val_cardio,
  function(i) {
    setDT(i)
    setnames(
      i,
      old = colnames(i),
      new = c("gene", "protein", "log2", "sdlog2",
      "pval", "bh", "fasta", "secreted")
    )
    i <- i[, c("gene", "log2", "pval", "bh")]
  }
)
val_cardio <- merge.data.table(
  merge.data.table(
    val_cardio[[1]],
    val_cardio[[2]],
    by = "gene",
    all = TRUE,
    suffixes = c("_h12", "_h24")
  ),
  val_cardio[[3]],
  by = "gene",
  all = TRUE
)
setnames(
  val_cardio,
  old = c("log2", "pval", "bh"),
  new = c("log2_h24_30", "pval_h24_30", "bh_h24_30")
)
val_cardio <- val_cardio[!is.na(gene)]

val_cardio[
  ,
  diff_h12 := ifelse(
    bh_h12 <= 0.05 & log2_h12 > 0,
    "UP",
    ifelse(
      bh_h12 <= 0.05 & log2_h12 < 0,
      "DOWN",
      "NS"
    )
  )
]
val_cardio[
  ,
  diff_h24 := ifelse(
    bh_h24 <= 0.05 & log2_h24 > 0,
    "UP",
    ifelse(
      bh_h24 <= 0.05 & log2_h24 < 0,
      "DOWN",
      "NS"
    )
  )
]
val_cardio[
  ,
  diff_h24_30 := ifelse(
    bh_h24_30 <= 0.05 & log2_h24_30 > 0,
    "UP",
    ifelse(
      bh_h24_30 <= 0.05 & log2_h24_30 < 0,
      "DOWN",
      "NS"
    )
  )
]
val_cardio[
  ,
  diff_h12 := ifelse(is.na(diff_h12), "NS", diff_h12)
]
val_cardio[
  ,
  diff_h24 := ifelse(is.na(diff_h24), "NS", diff_h24)
]
val_cardio[
  ,
  diff_h24_30 := ifelse(is.na(diff_h24_30), "NS", diff_h24_30)
]

val_cardio[
  ,
  diff_total := ifelse(
    diff_h12 == "UP" | diff_h24 == "UP" |
    diff_h24_30 == "UP",
    "UP",
    ifelse(
      diff_h12 == "DOWN" | diff_h24 == "DOWN" |
      diff_h24_30 == "DOWN",
      "DOWN",
      "NS"
    )
  )
]

# extract genes of proteins detected
val_cardio_genes <- val_cardio[
  ,
  c("gene")
]
anyNA(val_cardio_genes$gene)
any(duplicated(val_cardio_genes$gene))

# look for mouse-rat homologs

mart_db_rat <- biomaRt::useMart(
  "ensembl",
  dataset = "rnorvegicus_gene_ensembl",
  verbose = TRUE
)
mart_cardio_rat_from_ens <- biomaRt::getBM(
  attributes = c(
    "external_gene_name",
    "ensembl_gene_id"
  ),
  filters = "external_gene_name",
  mart = mart_db_rat,
  values = sort(unique(val_cardio_genes$gene))
)
table(
  sort(unique(val_cardio_genes$gene)) %in%
  mart_cardio_rat_from_ens$external_gene_name
)
sort(unique(val_cardio_genes$gene))[
  !sort(unique(val_cardio_genes$gene)) %in%
  mart_cardio_rat_from_ens$external_gene_name
]
setDT(mart_cardio_rat_from_ens)

val_cardio_genes <- merge.data.table(
  x = val_cardio_genes,
  y = mart_cardio_rat_from_ens,
  by.x = "gene",
  by.y = "external_gene_name",
  all.x = TRUE
)

mart_cardio_ortho <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "mmusculus_homolog_associated_gene_name",
    "mmusculus_homolog_orthology_confidence",
    "mmusculus_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  mart = mart_db_rat,
  values = unique(val_cardio_genes$ensembl_gene_id)
)
setDT(mart_cardio_ortho)
anyNA(mart_cardio_ortho$mmusculus_homolog_associated_gene_name)
mart_cardio_ortho <- mart_cardio_ortho[
  mmusculus_homolog_associated_gene_name != ""
]
val_cardio_genes <- merge.data.table(
  val_cardio_genes,
  mart_cardio_ortho[, 1:2],
  by = "ensembl_gene_id",
  all = TRUE
)
setnames(
  val_cardio_genes,
  old = "mmusculus_homolog_associated_gene_name",
  new = "mmusculus_homolog"
)
val_cardio_genes[
  ,
  mgi_symbol := ifelse(
    is.na(mgi_symbol),
    gene,
    mgi_symbol
  )
]
val_cardio_genes[mgi_symbol != gene]

# select genes that are in LRI
val_cardio_lr <- sort(val_cardio_genes[
  mgi_symbol %in% genes_lri_mouse
]$mgi_symbol)

# select genes that are in LRI and diff with hypoxia
val_cardio_lr_hyp_up <- unique(val_cardio_genes[
  gene %in% val_cardio[diff_total == "UP"]$gene &
  mgi_symbol %in% genes_lri_mouse
]$mgi_symbol)
val_cardio_lr_hyp_down <- unique(val_cardio_genes[
  gene %in% val_cardio[diff_total == "DOWN"]$gene &
  mgi_symbol %in% genes_lri_mouse
]$mgi_symbol)

#do brain, sasp similar as above
#do cci annotations


## code below in progress #####
########################################



dt_cci_rel[
  ,
  in_pancreas_lr_all := LIGAND_1 %in% pancreas_lr_all |
    LIGAND_2 %in% pancreas_lr_all |
    RECEPTOR_1 %in% pancreas_lr_all |
    RECEPTOR_2 %in% pancreas_lr_all |
    RECEPTOR_3 %in% pancreas_lr_all
]

dt_cci_rel[
  ,
  in_pancreas_lr_hpde := LIGAND_1 %in% val_pancreas_lr_hpde |
    LIGAND_2 %in% val_pancreas_lr_hpde |
    RECEPTOR_1 %in% val_pancreas_lr_hpde |
    RECEPTOR_2 %in% val_pancreas_lr_hpde |
    RECEPTOR_3 %in% val_pancreas_lr_hpde
]

dt_cci_rel[
  ,
  in_pancreas_lr_cancer_up := LIGAND_1 %in% val_pancreas_lr_can_up |
    LIGAND_2 %in% val_pancreas_lr_can_up |
    RECEPTOR_1 %in% val_pancreas_lr_can_up |
    RECEPTOR_2 %in% val_pancreas_lr_can_up |
    RECEPTOR_3 %in% val_pancreas_lr_can_up
]

dt_cci_rel[
  ,
  in_pancreas_lr_cancer_down := LIGAND_1 %in% val_pancreas_lr_can_down |
    LIGAND_2 %in% val_pancreas_lr_can_down |
    RECEPTOR_1 %in% val_pancreas_lr_can_down |
    RECEPTOR_2 %in% val_pancreas_lr_can_down |
    RECEPTOR_3 %in% val_pancreas_lr_can_down
]

test <- dt_cci_rel[
  tissue == "Pancreas" &
  in_pancreas_lr_cancer_up == TRUE
][, .N, by = c("dataset", "REGULATION", "EMITTER_CELLTYPE")][
  order(-N)
]

test2 <- dt_cci_rel[
  tissue == "Pancreas" &
  in_pancreas_lr_cancer_up == TRUE &
  EMITTER_CELLTYPE == "pancreatic ductal cell" &
  REGULATION == "UP"
]

test3 <- dt_cci_rel[
  tissue == "Pancreas" &
  in_pancreas_lr_cancer_up == TRUE &
  EMITTER_CELLTYPE == "pancreatic ductal cell" &
  REGULATION == "DOWN"
]

dt_cci_rel[, .N, by = c(
  "dataset", "tissue", "in_pancreas_lr_all"
)]

dt_cci_rel[
  , {
    totwt = .N
    .SD[, .(frac = .N / totwt), by = in_pancreas_lr_all]
  },
  by = c("dataset", "tissue")
][order(-frac)][in_pancreas_lr_all == TRUE]

dt_cci_pancreas <- dt_cci_rel[tissue == "Pancreas"]
table(dt_cci_pancreas$dataset, dt_cci_pancreas$in_pancreas_lr_all)

dt_cci_pancreas[, .N, by = c("EMITTER_CELLTYPE", "in_pancreas_lr_all")][
  order(-N)
]
dt_cci_pancreas[
  , {
    totwt <- .N
    .SD[, .(frac = .N / totwt), by = in_pancreas_lr_all]
  },
  by = c("EMITTER_CELLTYPE")
][order(-frac)]

table(dt_cci_rel$in_pancreas_lr_all)

dt_cci_rel[
  LIGAND_1 %in% pancreas_lr_all
][, .N, by = "tissue"][order(-N)][1:10]
dt_cci_rel[
  LIGAND_1 %in% pancreas_lr_all
][, .N, by = c("EMITTER_CELLTYPE", "tissue")][order(-N)][1:30]

test4 <- dt_cci_full_diffsex[
  dataset == "TMS FACS (young)" &
  tissue == "Pancreas" &
  LRI == "Tgm2:Adgrg1"
]

## Brain mouse secretome (bms) data preparation ####

bms_data <- fread(
  paste0(
    path_scagecom_input,
    "tushaus_2020_secretome.csv"
  )
)
bms_data <- bms_data[PG.Genes != ""]
table(
  bms_data$tms_genes %in% unique(unlist(seurat_genes))
)
bms_data$tms_genes[
  !bms_data$tms_genes %in% unique(unlist(seurat_genes))
]

bms_cell_type <- data.table(
    bms = c("Astro", "Microglia", "Neurons", "Oligodendrocytes"),
    tms = c("astrocyte", "microglial cell", "neuron", "oligodendrocyte")
)

bms_data[
  ,
  paste0("is_detected_Astro_", 1:6) := lapply(
    .SD,
    function(i) {!is.nan(i)}
  ),
  .SDcols = paste0("log2_LFQ_Astro_", 1:6)
]
bms_data[
  ,
  nrep_Astro := rowSums(.SD),
  .SDcols = paste0("is_detected_Astro_", 1:6)
]
bms_data[
  ,
  paste0("is_detected_Microglia_", 1:6) := lapply(
    .SD,
    function(i) {!is.nan(i)}
  ),
  .SDcols = paste0("log2_LFQ_Microglia_", 1:6)
]
bms_data[
  ,
  nrep_Microglia := rowSums(.SD),
  .SDcols = paste0("is_detected_Microglia_", 1:6)
]
bms_data[
  ,
  paste0("is_detected_Neurons_", 1:6) := lapply(
    .SD,
    function(i) {!is.nan(i)}
  ),
  .SDcols = paste0("log2_LFQ_Neurons_", 1:6)
]
bms_data[
  ,
  nrep_Neurons := rowSums(.SD),
  .SDcols = paste0("is_detected_Neurons_", 1:6)
]
bms_data[
  ,
  paste0("is_detected_Oligodendrocytes_", 1:6) := lapply(
    .SD,
    function(i) {!is.nan(i)}
  ),
  .SDcols = paste0("log2_LFQ_Oligodendrocytes_", 1:6)
]
bms_data[
  ,
  nrep_Oligodendrocytes := rowSums(.SD),
  .SDcols = paste0("is_detected_Oligodendrocytes_", 1:6)
]

ftable(
  bms_data$nrep_Astro > 4,
  bms_data$nrep_Microglia > 4,
  bms_data$nrep_Neurons > 4,
  bms_data$nrep_Oligodendrocytes > 4
)

## Annotate scAgeCom results with bms data ####

dt_cci_rel[
    ,
    L1_bms := ifelse(
        LIGAND_1 %in% bms_data$tms_genes, TRUE, FALSE
    )
]
dt_cci_rel[
    ,
    L2_bms := ifelse(
        LIGAND_2 %in% bms_data$tms_genes, TRUE, FALSE
    )
]
dt_cci_rel[
    ,
    R1_bms := ifelse(
        RECEPTOR_1 %in% bms_data$tms_genes, TRUE, FALSE
    )
]
dt_cci_rel[
    ,
    R2_bms := ifelse(
        RECEPTOR_2 %in% bms_data$tms_genes, TRUE, FALSE
    )
]
dt_cci_rel[
    ,
    R3_bms := ifelse(
        RECEPTOR_3 %in% bms_data$tms_genes, TRUE, FALSE
    )
]

table(
  dt_ora_rel[
    tissue == "Brain" & ORA_CATEGORY == "LIGAND_COMPLEX"
  ]$VALUE %in% bms_data$tms_genes
)

table(
  dt_cci_rel[
    tissue == "Brain" & EMITTER_CELLTYPE %in% bms_cell_type$tms
  ]$L1_bms
)
table(dt_cci_rel$L1_bms)

table(unique(dt_cci_rel[
  tissue == "Brain" & EMITTER_CELLTYPE == "neuron"
]$LIGAND_1) %in% bms_data[
    nrep_Neurons > 4
]$tms_genes)




table(unique(dt_cci_rel[
  tissue == "Liver" & EMITTER_CELLTYPE == "hepatocyte"
]$LIGAND_1) %in% bms_data[
    nrep_Neurons > 4
]$tms_genes)


test <- sapply(
  unique(dt_cci_rel$dataset),
  function(i) {
    sapply(
      unique(dt_cci_rel[dataset == i]$tissue),
      function(j) {
        sapply(
          unique(dt_cci_rel[dataset == i & tissue == j]$EMITTER_CELLTYPE),
          function(k) {
            temp <- unique(
              dt_cci_rel[
                dataset == i & tissue == j & EMITTER_CELLTYPE == k
              ]$LIGAND_1
            )
            sum(temp %in% bms_data[nrep_Neurons > 4]$tms_genes)/length(temp)
          }
        )
      }
    )
  }
)

table(
  unique(dt_lri_mouse$LIGAND_1) %in% bms_data$tms_genes
)

table(
  unique(bms_data$tms_genes) %in% unique(dt_lri_mouse$RECEPTOR_1)
)

testx <- sapply(
  unique(dt_cci_rel$dataset),
  function(i) {
    sapply(
      unique(dt_cci_rel[dataset == i]$tissue),
      function(j) {
        sapply(
          unique(dt_cci_rel[dataset == i & tissue == j]$EMITTER_CELLTYPE),
          function(k) {
            temp <- unique(
              dt_cci_rel[
                dataset == i & tissue == j & EMITTER_CELLTYPE == k
              ]$LIGAND_1
            )
            sum(temp %in% bms_data[nrep_Microglia > 4]$tms_genes)/length(temp)
          }
        )
      }
    )
  }
)

test2x <- data.table(
  x = unlist(testx),
  y = names(unlist(testx))
)


max(unlist(test))

test2 <- data.table(
  x = unlist(test),
  y = names(unlist(test))
)

hist(test2$x)

summary(test2x$x)

## ####

bms_spec_neuron <- bms_data[
    nrep_Neurons > 4 &
    nrep_Astro < 5 &
    nrep_Microglia < 5 &
    nrep_Oligodendrocytes < 5
]$PG.Genes

cci_bms[L1_th == TRUE & LIGAND_1 %in% bms_spec_neuron]

test <- cci_bms[L1_th == TRUE & LIGAND_1 %in% bms_spec_neuron][
    ,
    sum(CCI_SCORE_YOUNG),
    by = c("Dataset", "Tissue", "EMITTER_CELLTYPE", "LIGAND_1")
][order(-V1)]

test <- test[Tissue == "Brain"]

test <- cci_bms[L1_th == TRUE][
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

cci_bms[L1_th == TRUE]

cci_bms[, .N, by = c("Dataset", "R3_th")]

scd_brain_facs <- CCI_table[
  Tissue == "Brain"
]
scd_brain_facs_male <- CCI_table[
  Dataset == "TMS FACS (male)" & Tissue == "Brain"
]
scd_brain_facs_female <- CCI_table[
  Dataset == "TMS FACS (female)" & Tissue == "Brain"
]

## CCIs originating from bms cells ####

cci_bms <- scd_brain_facs[
    EMITTER_CELLTYPE %in% bms_cell_type$tms
]

table(cci_bms$EMITTER_CELLTYPE)

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
  bms_data$PG.Genes
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
              bms_data[nrep_Microglia > 4]$PG.Genes
            )
            t2 <- setdiff(
              ligs,
              bms_data[nrep_Microglia > 4]$PG.Genes
            )
            t3 <- setdiff(
              bms_data[nrep_Microglia > 4]$PG.Genes,
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
            #sum(ligs %in% bms_data[nrep_Astrocytes > 4]$PG.Genes) /
            #  length(ligs)
            #table(ligs %in% bms_data[nrep_Astrocytes > 4]$PG.Genes)
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
  bms_data,
  on = "LIGAND_1==PG.Genes",
  secretome_score := i.log2Avg_Asstrocytes
]

table(
    out_scores_astro$LIGAND_1 %in% bms_data[nrep_Astrocytes > 4]$PG.Genes
)

sum(out_scores_astro$secretome_score > 0, na.rm = TRUE)

table(!is.na(out_scores_astro$secretome_score) & !is.nan(out_scores_astro$secretome_score))

64/(64+111)

table(!is.na(out_scores_liver_nk$secretome_score) & !is.nan(out_scores_liver_nk$secretome_score))

15/75


out_scores_astro$ligand

ligand_background <- unique(c(
  unique(LRI_mouse_dt$LIGAND_1),
  bms_data$PG.Genes
))

out_scores_astro$ligand %in% ligand_background
bms_data[nrep_Astrocytes > 4]$PG.Genes %in% ligand_background

t1 <- intersect(
  out_scores_astro$ligand,
  bms_data[nrep_Astrocytes > 4]$PG.Genes
)
t2 <- setdiff(
  out_scores_astro$ligand,
  bms_data[nrep_Astrocytes > 4]$PG.Genes
)
t3 <- setdiff(
  bms_data[nrep_Astrocytes > 4]$PG.Genes,
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

bms_data[nrep_Astrocytes > 4][, ][, c("PG.Genes", "log2Avg_Asstrocytes")]

table(scd_brain_facs_male$REGULATION)

LRI_mouse_dt[
  ,
  L1_bms := ifelse(
    LIGAND_1 %in% bms_data$PG.Gene,
    TRUE, FALSE
  )
]

table(LRI_mouse_dt$L1_bms)

CCI_table[Tissue == "Brain" & EMITTER_CELLTYPE == "astrocyte"][, .N, by = "LRI"][order(-N)][
  LRI_mouse_dt,
  on = "LRI",
  L1_bms := i.L1_bms
][1:20]



bms_data$PG.Genes[
  !bms_data$PG.Genes %in% LRI_mouse_genes
]

## Load SASP atlas data ####

sasp_atlas <- fread(
  paste0(
    path_scagecom_input,
    "sasp_atlas.csv"
  )
)

## SASP Atlas orthologs ####

mart_sasp <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "entrezgene_id",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  mart = mart_db_human,
  values = unique(sasp_atlas$gene)
)
table(sasp_atlas$gene %in% mart_sasp$hgnc_symbol)

mart_sasp_ortho <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "mmusculus_homolog_associated_gene_name",
    "mmusculus_homolog_orthology_confidence",
    "mmusculus_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  mart = mart_db_human,
  values = mart_sasp$ensembl_gene_id
)
setDT(mart_sasp_ortho)
setDT(mart_sasp)

mart_sasp_ortho[
  mart_sasp,
  on = "ensembl_gene_id",
  (c("hgnc_symbol", "entrezgene_id")) := mget(
    c("i.hgnc_symbol", "i.entrezgene_id")
  )
]

table(
  sasp_atlas$gene %in% mart_sasp_ortho$hgnc_symbol
)

sasp_atlas <- merge.data.table(
  sasp_atlas,
  mart_sasp_ortho,
  by.x = "ensembl_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

setnames(
  sasp_atlas,
  old = "mmusculus_homolog_associated_gene_name",
  new = "mgi_symbol"
)

anyNA(sasp_atlas$mgi_symbol)
table(sasp_atlas$mgi_symbol == "")

sort(unique(sasp_atlas[
  mgi_symbol == ""
]$hgnc_symbol))

#annotation of remaining genes
hom_jax <- fread(
  paste0(
    path_scagecom_input,
    "HOM_MouseHumanSequence.rpt"
  )
)

table(sasp_atlas[
  mgi_symbol == ""
]$entrezgene_id %in% hom_jax$`EntrezGene ID`,
sasp_atlas[
  mgi_symbol == ""
]$hgnc_symbol %in% hom_jax$Symbol)

sasp_atlas$entrezgene_id
sasp_ortho_manual <- data.table(
  c("AKR1C1")
)

dt_lri_human[LIGAND_1 %in% unique(sasp_atlas[
  mgi_symbol == ""
]$hgnc_symbol)][ , 1:2]

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
