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

## Library ####

library(openxlsx)
library(ggplot2)

## hPDE secretomics (Li et al. 2022) ####

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
  mmusculus_homolog := ifelse(
    is.na(mmusculus_homolog),
    gene,
    mmusculus_homolog
  )
]
val_cardio_genes[mmusculus_homolog != gene]

# select genes that are in LRI
val_cardio_lr <- sort(val_cardio_genes[
  mmusculus_homolog %in% genes_lri_mouse
]$mmusculus_homolog)

# select genes that are in LRI and diff with hypoxia
val_cardio_lr_hyp_up <- unique(val_cardio_genes[
  gene %in% val_cardio[diff_total == "UP"]$gene &
  mmusculus_homolog %in% genes_lri_mouse
]$mmusculus_homolog)
val_cardio_lr_hyp_down <- unique(val_cardio_genes[
  gene %in% val_cardio[diff_total == "DOWN"]$gene &
  mmusculus_homolog %in% genes_lri_mouse
]$mmusculus_homolog)

## Brain mouse secretomics (Tushaus et al. 2020) ####

val_brain <- fread(
  paste0(
    path_scagecom_input,
    "tushaus_2020_secretome.csv"
  )
)
val_brain <- val_brain[PG.Genes != ""]

val_brain_celltype <- data.table(
    bms = c("Astro", "Microglia", "Neurons", "Oligodendrocytes"),
    tms = c("astrocyte", "microglial cell", "neuron", "oligodendrocyte")
)

val_brain[
  ,
  paste0("is_detected_Astro_", 1:6) := lapply(
    .SD,
    function(i) {
      !is.nan(i)
    }
  ),
  .SDcols = paste0("log2_LFQ_Astro_", 1:6)
]
val_brain[
  ,
  nrep_Astro := rowSums(.SD),
  .SDcols = paste0("is_detected_Astro_", 1:6)
]
val_brain[
  ,
  paste0("is_detected_Microglia_", 1:6) := lapply(
    .SD,
    function(i) {
      !is.nan(i)
    }
  ),
  .SDcols = paste0("log2_LFQ_Microglia_", 1:6)
]
val_brain[
  ,
  nrep_Microglia := rowSums(.SD),
  .SDcols = paste0("is_detected_Microglia_", 1:6)
]
val_brain[
  ,
  paste0("is_detected_Neurons_", 1:6) := lapply(
    .SD,
    function(i) {
      !is.nan(i)
    }
  ),
  .SDcols = paste0("log2_LFQ_Neurons_", 1:6)
]
val_brain[
  ,
  nrep_Neurons := rowSums(.SD),
  .SDcols = paste0("is_detected_Neurons_", 1:6)
]
val_brain[
  ,
  paste0("is_detected_Oligodendrocytes_", 1:6) := lapply(
    .SD,
    function(i) {
      !is.nan(i)
    }
  ),
  .SDcols = paste0("log2_LFQ_Oligodendrocytes_", 1:6)
]
val_brain[
  ,
  nrep_Oligodendrocytes := rowSums(.SD),
  .SDcols = paste0("is_detected_Oligodendrocytes_", 1:6)
]

# extract genes of proteins detected
val_brain_genes <- list(
  astro = unique(val_brain[nrep_Astro > 0]$tms_genes),
  micro = unique(val_brain[nrep_Microglia > 0]$tms_genes),
  neurons = unique(val_brain[nrep_Neurons > 0]$tms_genes),
  oligo = unique(val_brain[nrep_Oligodendrocytes > 0]$tms_genes)
)
sapply(val_brain_genes, anyNA)

# select genes that are in LRI
val_brain_lr <- lapply(
  val_brain_genes,
  function(i) {
    sort(i[i %in% genes_lri_mouse])
  }
)

val_brain_lr$glia <- unique(
  c(
    val_brain_lr$astro,
    val_brain_lr$micro,
    val_brain_lr$oligo
  )
)
names(val_brain_lr)

## SASP atlas (human) secretomics (Basisty et al., 2020) ####

val_sasp <- fread(
  paste0(
    path_scagecom_input,
    "sasp_atlas.csv"
  )
)

# extract genes of proteins detected
val_sasp_genes <- unique(val_sasp[
  ,
  c("gene", "ensembl_id")
])
any(duplicated(val_sasp_genes$gene))
any(duplicated(val_sasp_genes$ensembl_id))
anyNA(val_sasp_genes$gene)
anyNA(val_sasp_genes$ensembl_id)

# look for homologs

mart_sasp_from_ens <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "ensembl_gene_id",
  mart = mart_db_human,
  values = val_sasp_genes$ensembl_id
)
setDT(mart_sasp_from_ens)
mart_sasp_from_symb <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  mart = mart_db_human,
  values = val_sasp_genes$gene
)
setDT(mart_sasp_from_symb)

val_sasp_genes <- merge.data.table(
  x = val_sasp_genes,
  y = mart_sasp_from_ens,
  by.x = "ensembl_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)
val_sasp_genes <- merge.data.table(
  x = val_sasp_genes,
  y = mart_sasp_from_symb,
  by.x = "gene",
  by.y = "hgnc_symbol",
  all.x = TRUE
)

mart_sasp_ortho <- biomaRt::getBM(
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
        val_sasp_genes$ensembl_id,
        val_sasp_genes$ensembl_gene_id
      )
    )
  )
)
setDT(mart_sasp_ortho)
anyNA(mart_sasp_ortho$mmusculus_homolog_associated_gene_name)
mart_sasp_ortho <- mart_sasp_ortho[
  mmusculus_homolog_associated_gene_name != ""
]
val_sasp_genes <- merge.data.table(
  val_sasp_genes,
  mart_sasp_ortho[, 1:2],
  by.x = "ensembl_id",
  by.y = "ensembl_gene_id",
  all = TRUE
)
setnames(
  val_sasp_genes,
  old = "mmusculus_homolog_associated_gene_name",
  new = "mmusculus_homolog_1"
)
val_sasp_genes <- merge.data.table(
  val_sasp_genes,
  mart_sasp_ortho[, 1:2],
  by.x = "ensembl_gene_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)
setnames(
  val_sasp_genes,
  old = "mmusculus_homolog_associated_gene_name",
  new = "mmusculus_homolog_2"
)

hom_jax <- fread(
  paste0(
    path_scagecom_input,
    "HOM_MouseHumanSequence.rpt"
  )
)
colnames(hom_jax) <- make.names(colnames(hom_jax))
hom_jax <- hom_jax[
  DB.Class.Key %in% hom_jax[
    Symbol %in% val_sasp_genes[
      is.na(mmusculus_homolog_1) |
      is.na(mmusculus_homolog_2)
    ]$gene
  ]$DB.Class.Key
][, c(1, 2, 4)]
hom_jax <- hom_jax[
  Common.Organism.Name == "mouse, laboratory" |
  (Common.Organism.Name == "human" & Symbol %in%
  val_sasp_genes$gene)
]
hom_jax <- merge.data.table(
  hom_jax[
    Common.Organism.Name == "human"
  ][, c(1, 3)],
  hom_jax[
    Common.Organism.Name == "mouse, laboratory"
  ][, c(1, 3)],
  by = "DB.Class.Key",
  all = TRUE,
  suffixes = c("_human", "_mouse")
)

val_sasp_genes <- merge.data.table(
  val_sasp_genes,
  hom_jax[, c(2, 3)],
  by.x = "gene",
  by.y = "Symbol_human",
  all.x = TRUE,
  allow.cartesian = TRUE
)

# select genes that are in LRI
val_sasp_lr <- sort(
  unique(
    c(
      val_sasp_genes$mmusculus_homolog_1,
      val_sasp_genes$mmusculus_homolog_2,
      val_sasp_genes$Symbol_mouse
    )
  )
)
val_sasp_lr <- val_sasp_lr[val_sasp_lr %in% genes_lri_mouse]

# select genes that are in different conditions
select_sasp_genes <- function(
  ct,
  cond,
  regulation
) {
  if (regulation == "UP") {
    dt <- val_sasp[
          log2_ratio >= log2(1.5) & q_value <= 0.05 &
          cell_type == ct &
          treatment == cond
        ]
  } else if (regulation == "DOWN") {
    dt <- val_sasp[
          log2_ratio <= -log2(1.5) & q_value <= 0.05 &
          cell_type == ct &
          treatment == cond
        ]
  } else {
    stop("Wrong regulation!")
  }
  res <- unique(
    c(
      sort(
        unique(
          unlist(
            val_sasp_genes[
              gene %in% dt$gene
              ][
                ,
                c(
                  "mmusculus_homolog_1",
                  "mmusculus_homolog_2",
                  "Symbol_mouse"
                )
              ]
              )
            )
          ),
      sort(
        unique(
          unlist(
            val_sasp_genes[
              ensembl_id %in% dt$ensembl_id
            ][
              ,
              c(
                "mmusculus_homolog_1",
                "mmusculus_homolog_2",
                "Symbol_mouse"
              )
            ]
          )
        )
      )
    )
  )
  res[res %in% genes_lri_mouse]
}

val_sasp_lr_cond <- list(
  fibro_ras_up = select_sasp_genes("fibroblast", "ras_d7", "UP"),
  fibro_ataz_up = select_sasp_genes("fibroblast", "atazanavir", "UP"),
  fibro_xir_up = select_sasp_genes("fibroblast", "xir_d10", "UP"),
  epi_xir_up = select_sasp_genes("epithelial", "xir_d10", "UP"),
  fibro_ras_down = select_sasp_genes("fibroblast", "ras_d7", "DOWN"),
  fibro_ataz_down = select_sasp_genes("fibroblast", "atazanavir", "DOWN"),
  fibro_xir_down = select_sasp_genes("fibroblast", "xir_d10", "DOWN"),
  epi_xir_down = select_sasp_genes("epithelial", "xir_d10", "DOWN")
)
sapply(val_sasp_lr_cond, length)
intersect(
  intersect(
    val_sasp_lr_cond$fibro_ras_up,
    val_sasp_lr_cond$fibro_xir_up
  ),
  val_sasp_lr_cond$fibro_ataz_up
)

## BM macrophages (mouse) activated (Meissner et al. 2013) ####

val_macro <- openxlsx::read.xlsx(
  paste0(
    path_scagecom_input,
    "Meissner_macrophage.xlsx"
  ),
  startRow = 1
)
setDT(val_macro)
val_macro[
  ,
  proteins.kendall.up := ifelse(
    is.na(proteins.kendall.up),
    "-",
    proteins.kendall.up
  )
]
val_macro[
  ,
  proteins.kendall.down := ifelse(
    is.na(proteins.kendall.down),
    "-",
    proteins.kendall.down
  )
]
val_macro[
  ,
  direction := ifelse(
    proteins.kendall.up == "+",
    "up",
    ifelse(
      proteins.kendall.down == "+",
      "down",
      "ns"
    )
  )
]
table(val_macro$direction)

val_macro <- val_macro[
  !is.na(Gene.names)
]

# extract genes of proteins detected
val_macro_genes <- val_macro[
  ,
  c("Protein.IDs", "Gene.names")
]
anyNA(val_macro_genes$Gene.names)
which(duplicated(val_macro_genes$Gene.names))
anyNA(val_macro_genes$Protein.IDs)
any(duplicated(val_macro_genes$Protein.IDs))

val_macro_genes <- val_macro_genes[
  ,
  lapply(
    .SD,
    function(x) unlist(tstrsplit(x, ";", fixed = TRUE))
  ),
  by = Protein.IDs
]

# check gene symbols
table(
  val_macro_genes$Gene.names %in% unique(unlist(seurat_genes)),
  val_macro_genes$Gene.names %in% dt_org_mm$SYMBOL
)
table(
  val_macro_genes$Gene.names %in% unique(unlist(seurat_genes)),
  val_macro_genes$Gene.names %in% dt_org_mm$ALIAS
)
val_macro_genes[
  ,
  in_tms := Gene.names %in% unique(unlist(seurat_genes))
]
val_macro_genes[
  ,
  in_orgdb := Gene.names %in% dt_org_mm$SYMBOL
]

val_macro_genes <- merge.data.table(
  val_macro_genes,
  dt_org_mm[
    ALIAS %in% val_macro_genes[
      in_tms == FALSE
    ]$Gene.names
  ][, c("SYMBOL", "ALIAS")],
  by.x = "Gene.names",
  by.y = "ALIAS",
  all.x = TRUE
)

val_macro_genes[
  ,
  mgi_symbol := ifelse(
    in_tms == TRUE,
    Gene.names,
    SYMBOL
  )
]
sum(is.na(val_macro_genes$mgi_symbol))
table(val_macro_genes$mgi_symbol %in% unique(unlist(seurat_genes)))
 
# select genes that are in LRI
val_macro_lr <- val_macro_genes[
  mgi_symbol %in% genes_lri_mouse
]$mgi_symbol

# select genes that are regulated by activation
val_macro_lr_up <- val_macro_genes[
  Protein.IDs %in% val_macro[
    direction == "up"
  ]$Protein.IDs &
  mgi_symbol %in% genes_lri_mouse
]$mgi_symbol

val_macro_lr_down <- val_macro_genes[
  Protein.IDs %in% val_macro[
    direction == "down"
  ]$Protein.IDs &
  mgi_symbol %in% genes_lri_mouse
]$mgi_symbol

## Endo cells (human) secretomics (Brioschi et al. 2013) ####

# read spreadsheet
val_endo <- read.xlsx(
  paste0(
    path_scagecom_input,
    "Brioschi_endo_secretome.xlsx"
  ),
  startRow = 2
)
setDT(val_endo)
val_endo[
  ,
  UniProt.Accession := sub("*", "", UniProt.Accession, fixed = TRUE)
]

# extract protein ids
val_endo_genes <- val_endo[
  ,
  c("UniProt.Accession", "Description")
]
anyNA(val_endo_genes$UniProt.Accession)
any(duplicated(val_endo_genes$UniProt.Accession))

# look for gene names and homologs

mart_endo_human <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "uniprotswissprot"
  ),
  filters = "uniprotswissprot",
  mart = mart_db_human,
  values = sort(unique(val_endo_genes$UniProt.Accession))
)
setDT(mart_endo_human)
table(val_endo_genes$UniProt.Accession %in% mart_endo_human$uniprotswissprot)

val_endo_genes <- merge.data.table(
  x = val_endo_genes,
  y = mart_endo_human,
  by.x = "UniProt.Accession",
  by.y = "uniprotswissprot",
  all = TRUE
)

mart_endo_ortho <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "mmusculus_homolog_associated_gene_name",
    "mmusculus_homolog_orthology_confidence",
    "mmusculus_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  mart = mart_db_human,
  values = unique(val_endo_genes$ensembl_gene_id)
)
setDT(mart_endo_ortho)
anyNA(mart_endo_ortho$mmusculus_homolog_associated_gene_name)
mart_endo_ortho <- mart_endo_ortho[
  mmusculus_homolog_associated_gene_name != ""
]
val_endo_genes <- merge.data.table(
  val_endo_genes,
  mart_endo_ortho[, 1:2],
  by.x = "ensembl_gene_id",
  by.y = "ensembl_gene_id",
  all = TRUE
)
setnames(
  val_endo_genes,
  old = "mmusculus_homolog_associated_gene_name",
  new = "mmusculus_homolog"
)
table(
  val_endo$UniProt.Accession %in%
  val_endo_genes[!is.na(mmusculus_homolog)]$UniProt.Accession
)

# select genes that are in LRI
val_endo_lr <- sort(val_endo_genes[
  mmusculus_homolog %in% genes_lri_mouse
]$mmusculus_homolog)

## Huvec secretomics (Zhao et al. 2020) ####
val_huvec <- read.xlsx(
  paste0(
    path_scagecom_input,
    "zhao_huvec_secretomics.xlsx"
  ),
  startRow = 3
)
setDT(val_huvec)

# extract genes of proteins detected
val_huvec_genes <- val_huvec[
  ,
  c("Protein.IDs", "Gene.names")
]

# look for homologs
mart_huvec_human <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  mart = mart_db_human,
  values = sort(unique(val_huvec_genes$Gene.names))
)
setDT(mart_huvec_human)
val_huvec_genes <- merge.data.table(
  x = val_huvec_genes,
  y = mart_huvec_human,
  by.x = "Gene.names",
  by.y = "hgnc_symbol",
  all.x = TRUE
)
mart_huvec_ortho <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "mmusculus_homolog_associated_gene_name",
    "mmusculus_homolog_orthology_confidence",
    "mmusculus_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  mart = mart_db_human,
  values = sort(unique(val_huvec_genes$ensembl_gene_id))
)
setDT(mart_huvec_ortho)
anyNA(mart_huvec_ortho$mmusculus_homolog_associated_gene_name)
mart_huvec_ortho <- mart_huvec_ortho[
  mmusculus_homolog_associated_gene_name != ""
]
val_huvec_genes <- merge.data.table(
  val_huvec_genes,
  mart_huvec_ortho[, 1:2],
  by.x = "ensembl_gene_id",
  by.y = "ensembl_gene_id",
  all = TRUE
)
setnames(
  val_huvec_genes,
  old = "mmusculus_homolog_associated_gene_name",
  new = "mmusculus_homolog_1"
)

# select genes that are in LRI
val_huvec_lr <- sort(val_huvec_genes[
  mmusculus_homolog_1 %in% genes_lri_mouse
]$mmusculus_homolog_1)

## MSC AT mouse secretomics (Acar et al. 2020) ####

val_mscat <- read.xlsx(
  paste0(
    path_scagecom_input,
    "acar_mscat_secretome.xlsx"
  ),
  startRow = 1
)
setDT(val_mscat)

# extract proteins detected
val_mscat_genes <- val_mscat[, c("Protein_id")]

# look for gene symbols
mart_db_mouse_2 <- biomaRt::useMart(
  "ensembl",
  host = "https://dec2017.archive.ensembl.org",
  dataset = "mmusculus_gene_ensembl",
  verbose = TRUE
)
mart_mscat_mouse <- biomaRt::getBM(
  attributes = c(
    "mgi_symbol",
    "ensembl_gene_id",
    "uniprotswissprot"
  ),
  filters = "uniprotswissprot",
  mart = mart_db_mouse,
  values = sort(unique(val_mscat_genes$Protein_id))
)
setDT(mart_mscat_mouse)
table(val_mscat_genes$Protein_id %in% mart_mscat_mouse$uniprotswissprot)

val_mscat_genes <- merge.data.table(
  x = val_mscat_genes,
  y = mart_mscat_mouse,
  by.x = "Protein_id",
  by.y = "uniprotswissprot",
  all = TRUE
)

table(is.na(val_mscat_genes$mgi_symbol))

# select genes that are in LRI
val_mscat_lr <- sort(val_mscat_genes[
  mgi_symbol %in% genes_lri_mouse
]$mgi_symbol)

# select genes upregulated with HFD
val_mscat_lr_hfd <- sort(val_mscat_genes[
  mgi_symbol %in% genes_lri_mouse &
  Protein_id %in% val_mscat[Condition == "HFD"]$Protein_id
]$mgi_symbol)

## MSC human secretomics (Shin et al. 2021) ####

# read spreadsheet
val_msccomb <- read.xlsx(
  paste0(
    path_scagecom_input,
    "shin_msc_secretome.xlsx"
  ),
  startRow = 3
)
setDT(val_msccomb)

# extract genes
val_msccomb_genes <- val_msccomb[
  ,
  c("Accession", "SYM")
]
anyNA(val_msccomb_genes$Accession)
any(duplicated(val_msccomb_genes$Accession))

# look for homologs
mart_msccomb_human <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  mart = mart_db_human,
  values = sort(unique(val_msccomb_genes$SYM))
)
setDT(mart_msccomb_human)

val_msccomb_genes <- merge.data.table(
  x = val_msccomb_genes,
  y = mart_msccomb_human,
  by.x = "SYM",
  by.y = "hgnc_symbol",
  all = TRUE
)
mart_msccomb_ortho <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "mmusculus_homolog_associated_gene_name",
    "mmusculus_homolog_orthology_confidence",
    "mmusculus_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  mart = mart_db_human,
  values = unique(val_msccomb_genes$ensembl_gene_id)
)
setDT(mart_msccomb_ortho)
anyNA(mart_msccomb_ortho$mmusculus_homolog_associated_gene_name)
mart_msccomb_ortho <- mart_msccomb_ortho[
  mmusculus_homolog_associated_gene_name != ""
]
val_msccomb_genes <- merge.data.table(
  val_msccomb_genes,
  mart_msccomb_ortho[, 1:2],
  by.x = "ensembl_gene_id",
  by.y = "ensembl_gene_id",
  all = TRUE
)
setnames(
  val_msccomb_genes,
  old = "mmusculus_homolog_associated_gene_name",
  new = "mmusculus_homolog"
)

# select genes that are in LRI
val_msccomb_lr <- sort(val_msccomb_genes[
  mmusculus_homolog %in% genes_lri_mouse
]$mmusculus_homolog)

## list of all conditions ####

val_all <- c(
  list(
    pancreas = val_pancreas_lr,
    pancreas_hpde = val_pancreas_lr_hpde,
    pancreas_can_up = val_pancreas_lr_can_up,
    pancreas_can_down = val_pancreas_lr_can_down,
    cardio = val_cardio_lr,
    cardio_hyp_up = val_cardio_lr_hyp_up,
    cardio_hyp_down = val_cardio_lr_hyp_down,
    macro = val_macro_lr,
    macro_up = val_macro_lr_up,
    macro_down = val_macro_lr_down,
    endo = val_endo_lr,
    huvec = val_huvec_lr,
    mscat = val_mscat_lr,
    mscat_hfd = val_mscat_lr_hfd,
    msccomb = val_msccomb_lr
  ),
  val_brain_lr,
  val_sasp_lr_cond
)

val_all <- val_all[lapply(val_all, length) > 0]
names(val_all)

## Upset plot of secreted lris ####

val_dt <- data.table(
  gene = sort(unique(unlist(val_all)))
)
val_dt[
  ,
  (names(val_all)) := lapply(
    val_all,
    function(i) {
      gene %in% i
    }
  )
]

val_upset_keep <- c(
  "pancreas", "huvec", "macro", "mscat",
  "neurons", "glia", "cardio"
)

val_dt_selection <- val_dt[
  ,
  c("gene", val_upset_keep),
  with = FALSE
][,
  tot := Reduce(`|`, lapply(.SD, function(x) {
    x
  })),
  .SDcols = val_upset_keep
][tot == TRUE]
val_dt_selection[
  ,
  N := rowSums(.SD),
  .SDcols = val_upset_keep
]

sort(val_dt_selection[N >= 6]$gene)

ComplexUpset::upset(
  as.data.frame(
    val_dt[
      ,
      c("gene", val_upset_keep),
      with = FALSE
    ][,
     tot := Reduce(`|`, lapply(.SD, function(x) {x})),
      .SDcols = val_upset_keep
    ][tot == TRUE]
  ),
  val_upset_keep,
  name = "xx",
  set_sizes = ComplexUpset::upset_set_size(),
  min_size = 6
) + ggtitle(
  "xx"
)

## Annotate LRIs and CCIs with secretomics genes ####

# annoated LRI table first
dt_lri_mouse_val <- copy(dt_lri_mouse)
dt_lri_mouse_val[
  ,
  (names(val_all)) := lapply(
    val_all,
    function(i) {
      ifelse(
        LIGAND_1 %in% i,
        "L1",
        ifelse(
          LIGAND_2 %in% i,
          "L2",
          ifelse(
            RECEPTOR_1 %in% i,
            "R1",
            ifelse(
              RECEPTOR_2 %in% i,
              "R2",
              ifelse(
                RECEPTOR_3 %in% i,
                "R3",
                "NO"
              )
            )
          )
        )
      )
    }
  )
]
table(dt_lri_mouse_val$pancreas)
table(dt_lri_mouse_val$macro)

# annotate CCIs
dt_cci_rel <- merge.data.table(
  dt_cci_rel,
  dt_lri_mouse_val[, c("LRI", names(val_all)), with = FALSE],
  by = "LRI",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

dt_cci_rel[
  ,
  dataset_tissue := paste(
    dataset,
    tissue,
    sep = "_"
  )
]

dt_cci_sex <- merge.data.table(
  dt_cci_sex,
  dt_lri_mouse_val[, c("LRI", names(val_all)), with = FALSE],
  by = "LRI",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

dt_cci_sex[
  ,
  dataset_tissue := paste(
    dataset,
    tissue,
    sep = "_"
  )
]

## ORA function on validation #####

run_ora_validation <- function(
  data,
  category,
  val_set,
  by_group = "none"
) {
  if (by_group == "none") {
    bg <- NULL
  } else {
    bg <- by_group
  }
  dt1 <- data[
    get(val_set) != "NO"
  ][, .(CAT_LR = .N), by = c(category, bg)]
  dt1[, notCAT_LR := sum(CAT_LR) - CAT_LR, by = bg]
  dt2 <- data[
    get(val_set) == "NO"
  ][, .(CAT_notLR = .N), by = c(category, bg)]
  dt2[, notCAT_notLR := sum(CAT_notLR) - CAT_notLR, by = bg]
  dt <- merge.data.table(
    dt1,
    dt2,
    by = c(category, bg)
  )
  dt[
    ,
    c("OR", "pval") :=
    scDiffCom:::vfisher_2sided(
      CAT_LR, CAT_notLR, notCAT_LR, notCAT_notLR
    )
  ]
  dt[
    ,
    BH := p.adjust(pval, method = "BH")
  ]
  dt[, category_type := category]
  dt[, group_type := by_group]
  if (by_group == "none") {
    dt[, group := "none"]
    setnames(
      dt,
      old = c(category),
      new = c("category")
    )
  } else {
    setnames(
      dt,
      old = c(category, by_group),
      new = c("category", "group")
    )
  }
  dt[, validation_set := val_set]
  return(dt)
}

## Compute ORA validation ####

dt_ora_validation_groups <- data.table(
  category = c(
    "dataset", "tissue", "tissue",
    "EMITTER_CELLTYPE", "EMITTER_CELLTYPE", "EMITTER_CELLTYPE",
    "RECEIVER_CELLTYPE", "RECEIVER_CELLTYPE", "RECEIVER_CELLTYPE",
    "EMITTER_CELLFAMILY", "EMITTER_CELLFAMILY", "EMITTER_CELLFAMILY",
    "RECEIVER_CELLFAMILY", "RECEIVER_CELLFAMILY", "RECEIVER_CELLFAMILY",
    "EMITTER_CELLFAMILY_2", "EMITTER_CELLFAMILY_2", "EMITTER_CELLFAMILY_2",
    "RECEIVER_CELLFAMILY_2", "RECEIVER_CELLFAMILY_2", "RECEIVER_CELLFAMILY_2",
    "ER_CELLTYPES", "ER_CELLTYPES", "ER_CELLTYPES",
    "ER_CELLFAMILY", "ER_CELLFAMILY", "ER_CELLFAMILY",
    "ER_CELLFAMILY_2", "ER_CELLFAMILY_2", "ER_CELLFAMILY_2"
  ),
  by_group = c(
    "none", "none", "dataset",
    "none", "dataset", "dataset_tissue",
    "none", "dataset", "dataset_tissue",
    "none", "dataset", "dataset_tissue",
    "none", "dataset", "dataset_tissue",
    "none", "dataset", "dataset_tissue",
    "none", "dataset", "dataset_tissue",
    "none", "dataset", "dataset_tissue",
    "none", "dataset", "dataset_tissue",
    "none", "dataset", "dataset_tissue"
  )
)

dt_ora_validation_groups <- setkey(
  data.table(val_set = names(val_all))[, c(k = 1, .SD)],
  k
)[dt_ora_validation_groups[, c(k = 1, .SD)],
allow.cartesian = TRUE][, k := NULL]

options(future.globals.maxSize = 15 * 1024^3)
future::plan(multicore, workers = 20)

dt_ora_validation_groups_p <- dt_ora_validation_groups[
  category %in% c(
    "EMITTER_CELLTYPE",
    "RECEIVER_CELLTYPE",
    #"ER_CELLTYPES",
    "ER_CELLFAMILY",
    "ER_CELLFAMILY_2"
  ) &
  by_group %in% c(
    "dataset_tissue",
    "none",
    "dataset"
  )
]

dt_ora_validation_groups_p2 <- dt_ora_validation_groups[
  category %in% c(
    "EMITTER_CELLFAMILY",
    "EMITTER_CELLFAMILY_2",
    "RECEIVER_CELLFAMILY",
    "RECEIVER_CELLFAMILY_2"
  ) &
  by_group %in% c(
    "dataset_tissue",
    "none",
    "dataset"
  )
]

dt_cci_secr_validation <- rbindlist(
  future.apply::future_lapply(
    seq_len(nrow(dt_ora_validation_groups_p)),
    function(i) {
      run_ora_validation(
        data = dt_cci_rel,
        category = dt_ora_validation_groups_p$category[[i]],
        val_set = dt_ora_validation_groups_p$val_set[[i]],
        by_group = dt_ora_validation_groups_p$by_group[[i]]
      )
    }
  ),
  use.names = TRUE
)

dt_cci_secr_validation2 <- rbindlist(
  future.apply::future_lapply(
    seq_len(nrow(dt_ora_validation_groups_p2)),
    function(i) {
      run_ora_validation(
        data = dt_cci_rel,
        category = dt_ora_validation_groups_p2$category[[i]],
        val_set = dt_ora_validation_groups_p2$val_set[[i]],
        by_group = dt_ora_validation_groups_p2$by_group[[i]]
      )
    }
  ),
  use.names = TRUE
)

dt_cci_secr_validation_full <- rbindlist(
  list(
    dt_cci_secr_validation,
    dt_cci_secr_validation2
  )
)

dt_cci_secr_validation_full[
   ,
  recall := CAT_LR / (CAT_LR + CAT_notLR)
]

fwrite(
  dt_cci_secr_validation_full,
  paste0(
    path_scagecom_output,
    "cci_validation.csv"
  )
)

## Specific results of interest to select ####

# huvec
dt_val_huvec_fam <- dt_cci_secr_validation_full[
  validation_set == "huvec" &
  group_type == "none" &
  category_type == "EMITTER_CELLFAMILY_2"
][category != "leukocyte"]

#macro
dt_val_macro_fam <- dt_cci_secr_validation_full[
  validation_set == "macro" &
  group_type == "none" &
  category_type == "EMITTER_CELLFAMILY_2"
][category != "leukocyte"]

#mscat
dt_val_mscat_fam <- dt_cci_secr_validation_full[
  validation_set == "mscat" &
  group_type == "none" &
  category_type == "EMITTER_CELLFAMILY_2"
][category != "leukocyte"]

dt_val_mscat_at <-   dt_cci_secr_validation_full[
  validation_set == "mscat" &
  group_type == "dataset_tissue" &
  category_type == "EMITTER_CELLTYPE"
][grepl("Adipose", group)]

#neurons
dt_val_neuron_fam <-   dt_cci_secr_validation_full[
  validation_set == "neurons" &
  group_type == "none" &
  category_type == "EMITTER_CELLFAMILY_2"
][category != "leukocyte"]

#glial
dt_val_glial_fam <-   dt_cci_secr_validation_full[
  validation_set == "glia" &
  group_type == "none" &
  category_type == "EMITTER_CELLFAMILY_2"
][category != "leukocyte"]

#cardio
dt_val_cardio_fam <-   dt_cci_secr_validation_full[
  validation_set == "cardio" &
  group_type == "none" &
  category_type == "EMITTER_CELLFAMILY_2"
][category != "leukocyte"]

dt_val_cardio_heart <- dt_cci_secr_validation_full[
  validation_set == "cardio" &
  group_type == "dataset_tissue" &
  category_type == "EMITTER_CELLTYPE"
][grepl("Heart", group)]

#pancreas
dt_val_pancreas_fam <-   dt_cci_secr_validation_full[
  validation_set == "pancreas" &
  group_type == "none" &
  category_type == "EMITTER_CELLFAMILY_2"
][category != "leukocyte"]

dt_val_pancreas_pct <-   dt_cci_secr_validation_full[
  validation_set == "pancreas" &
  group_type == "dataset_tissue" &
  category_type == "EMITTER_CELLTYPE"
][grepl("Pancreas", group)]

## Rename validation in cci and related datasets ####

dt_cci_age <- copy(dt_cci_rel)
dt_cci_age[
  ,
  c(
    "pancreas_hpde", "pancreas_can_up", "pancreas_can_down",
    "cardio_hyp_up", "cardio_hyp_down", "macro_up", "macro_down",
    "endo", "mscat_hfd", "msccomb", "astro", "micro", "oligo",
    "fibro_ras_up", "fibro_ataz_up", "fibro_xir_up", "fibro_ras_down",
    "fibro_xir_down", "epi_xir_up", "epi_xir_down"
  ) := NULL
]
setnames(
  dt_cci_age,
  old = c("pancreas", "cardio", "macro", "huvec", "mscat", "neurons", "glia"),
  new = c("hPDE", "rCM", "mBBM", "hUVEC", "mMSC-AT", "mNeuron", "mGlial")
)

dt_cci_sex[
  ,
  c(
    "pancreas_hpde", "pancreas_can_up", "pancreas_can_down",
    "cardio_hyp_up", "cardio_hyp_down", "macro_up", "macro_down",
    "endo", "mscat_hfd", "msccomb", "astro", "micro", "oligo",
    "fibro_ras_up", "fibro_ataz_up", "fibro_xir_up", "fibro_ras_down",
    "fibro_xir_down", "epi_xir_up", "epi_xir_down"
  ) := NULL
]
setnames(
  dt_cci_sex,
  old = c("pancreas", "cardio", "macro", "huvec", "mscat", "neurons", "glia"),
  new = c("hPDE", "rCM", "mBBM", "hUVEC", "mMSC-AT", "mNeuron", "mGlial")
)

dt_secr_validation_final <- dt_cci_secr_validation_full[
  validation_set %in% c(
    "pancreas", "cardio", "macro", "huvec", "mscat", "neurons", "glia"
  )
]
dt_secr_validation_final[
  ,
  validation_set := ifelse(
    validation_set == "pancreas",
    "hPDE",
    ifelse(
      validation_set == "cardio",
      "rCM",
      ifelse(
        validation_set == "macro",
        "mBBM",
        ifelse(
          validation_set == "huvec",
          "hUVEC",
          ifelse(
            validation_set == "mscat",
            "mMSC-AT",
            ifelse(
              validation_set == "neurons",
              "mNeuron",
              "mGlial"
            )
          )
        )
      )
    )
  )
]
table(dt_secr_validation_final$validation_set)

saveRDS(
  dt_secr_validation_final,
  paste0(
    path_scagecom_output,
    "dt_secr_validation_final.rds"
  )
)
fwrite(
  dt_secr_validation_final,
  paste0(
    path_scagecom_output,
    "dt_secr_validation_final.csv"
  )
)

saveRDS(
  dt_cci_age,
  paste0(
    path_scagecom_output,
    "dt_cci_age.rds"
  )
)
fwrite(
  dt_cci_age,
  paste0(
    path_scagecom_output,
    "dt_cci_age.csv"
  )
)

saveRDS(
  dt_cci_sex,
  paste0(
    path_scagecom_output,
    "dt_cci_sex.rds"
  )
)
fwrite(
  dt_cci_sex,
  paste0(
    path_scagecom_output,
    "dt_cci_sex.csv"
  )
)

dt_lri_mouse_val_clean <- copy(dt_lri_mouse_val)
dt_lri_mouse_val_clean[
  ,
  c(
    "pancreas_hpde", "pancreas_can_up", "pancreas_can_down",
    "cardio_hyp_up", "cardio_hyp_down", "macro_up", "macro_down",
    "endo", "mscat_hfd", "msccomb", "astro", "micro", "oligo",
    "fibro_ras_up", "fibro_ataz_up", "fibro_xir_up", "fibro_ras_down",
    "fibro_xir_down", "epi_xir_up", "epi_xir_down"
  ) := NULL
]
setnames(
  dt_lri_mouse_val_clean,
  old = c("pancreas", "cardio", "macro", "huvec", "mscat", "neurons", "glia"),
  new = c("hPDE", "rCM", "mBBM", "hUVEC", "mMSC-AT", "mNeuron", "mGlial")
)

dt_lri_mouse_val_clean[
  dcast.data.table(
    melt.data.table(
      dt_lri_mouse_val_clean[
        ,
        c(
          "LRI", "hPDE", "rCM", "mBBM", "hUVEC",
          "mMSC-AT", "mNeuron", "mGlial"
        )
      ],
      id.vars = "LRI"
    )[
      ,
      value := value != "NO"
    ][
      value == TRUE
    ],
    LRI ~ .,
    value.var = "variable",
    fun.aggregate = function(i) {
      paste(sort(i), collapse = "/")
    }
  ),
  on = "LRI",
  summary_val := i..
]

dt_lri_mouse_val_clean[
  dcast.data.table(
    melt.data.table(
      dt_lri_mouse_val_clean[
        ,
        c(
          "LRI", "L1_N_agepmid", "L2_N_agepmid",
          "R1_N_agepmid", "R2_N_agepmid", "R3_N_agepmid"
        )
      ],
      id.vars = "LRI"
    )[
      !is.na(value)
    ],
    LRI ~ .,
    value.var = "value",
    fun.aggregate = function(i) {
      paste(i, collapse = "_")
    }
  ),
  on = "LRI",
  pubmed := i..
]
