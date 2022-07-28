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
    function(i) {!is.nan(i)}
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
    function(i) {!is.nan(i)}
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
    function(i) {!is.nan(i)}
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
    function(i) {!is.nan(i)}
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

ComplexUpset::upset(
  as.data.frame(val_dt),
  c(
    "pancreas", "cardio", "macro", "astro",
    "micro", "neurons", "oligo", "fibro_xir_up",
    "epi_xir_up"
  ),
  name = "xx",
  set_sizes = ComplexUpset::upset_set_size(),
  min_size = 5
) + ggtitle(
  "xx"
)

## Annotate LRIs and CCIs with secretomics genes ####

# list of all conditions
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
    macro_down = val_macro_lr_down
  ),
  val_brain_lr,
  val_sasp_lr_cond
)

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

## Frequency table of CCIs by validation group ####

dt_cci_val_tissue <- rbindlist(
  lapply(
    setNames(
      names(val_all),
      names(val_all)
    ),
    function(val_name) {
      dt <- dt_cci_rel[
      , {
        totwt <- .N
        .SD[, .(frac = .N / totwt, N = .N), by = val_name]
      },
      by = c("dataset", "tissue")
    ]
    setnames(
      dt,
      old = val_name,
      new = "lr_type"
    )
    }
  ),
  idcol = "validation_set"
)

dt_cci_val_tissue_reg <- rbindlist(
  lapply(
    setNames(
      names(val_all),
      names(val_all)
    ),
    function(val_name) {
      dt <- dt_cci_rel[
      , {
        totwt <- .N
        .SD[, .(frac = .N / totwt, N = .N), by = val_name]
      },
      by = c("dataset", "tissue", "REGULATION")
    ]
    setnames(
      dt,
      old = val_name,
      new = "lr_type"
    )
    }
  ),
  idcol = "validation_set"
)

dt_cci_val_emmitter <- rbindlist(
  lapply(
    setNames(
      names(val_all),
      names(val_all)
    ),
    function(val_name) {
      dt <- dt_cci_rel[
      , {
        totwt <- .N
        .SD[, .(frac = .N / totwt, N = .N), by = val_name]
      },
      by = c("dataset", "tissue", "EMITTER_CELLTYPE")
    ]
    setnames(
      dt,
      old = val_name,
      new = "lr_type"
    )
    }
  ),
  idcol = "validation_set"
)

dt_cci_val_emmitter_reg <- rbindlist(
  lapply(
    setNames(
      names(val_all),
      names(val_all)
    ),
    function(val_name) {
      dt <- dt_cci_rel[
      , {
        totwt <- .N
        .SD[, .(frac = .N / totwt, N = .N), by = val_name]
      },
      by = c("dataset", "tissue", "EMITTER_CELLTYPE", "REGULATION")
    ]
    setnames(
      dt,
      old = val_name,
      new = "lr_type"
    )
    }
  ),
  idcol = "validation_set"
)

dt_cci_val_receiver <- rbindlist(
  lapply(
    setNames(
      names(val_all),
      names(val_all)
    ),
    function(val_name) {
      dt <- dt_cci_rel[
      , {
        totwt <- .N
        .SD[, .(frac = .N / totwt, N = .N), by = val_name]
      },
      by = c("dataset", "tissue", "RECEIVER_CELLTYPE")
    ]
    setnames(
      dt,
      old = val_name,
      new = "lr_type"
    )
    }
  ),
  idcol = "validation_set"
)

dt_cci_val_receiver_reg <- rbindlist(
  lapply(
    setNames(
      names(val_all),
      names(val_all)
    ),
    function(val_name) {
      dt <- dt_cci_rel[
      , {
        totwt <- .N
        .SD[, .(frac = .N / totwt, N = .N), by = val_name]
      },
      by = c("dataset", "tissue", "RECEIVER_CELLTYPE", "REGULATION")
    ]
    setnames(
      dt,
      old = val_name,
      new = "lr_type"
    )
    }
  ),
  idcol = "validation_set"
)

write.xlsx(
  x = list(
    dt_cci_val_tissue,
    dt_cci_val_tissue_reg,
    dt_cci_val_emmitter,
    dt_cci_val_emmitter_reg,
    dt_cci_val_receiver,
    dt_cci_val_receiver_reg
  ),
  paste0(
    path_scagecom_output,
    "cci_validation_tables.xlsx"
  )
)


dt_cci_val_tissue[
  lr_type == "NO" &
  validation_set == "pancreas"
]

ggplot(
  dt_cci_val_tissue[
    lr_type == "NO" &
    validation_set == "macro"
  ],
  aes(
    x = reorder(paste(dataset, tissue), frac),
    y = (1 - frac) * 100
  )
) + geom_point() + coord_flip()

ggplot(
  dt_cci_val_tissue_reg[
    lr_type == "NO" &
    validation_set == "cardio_hyp_down"
  ][frac < 0.95],
  aes(
    x = reorder(paste(dataset, tissue, REGULATION), frac),
    y = (1 - frac) * 100
  )
) + geom_point(
) + coord_flip()

ggplot(
  dt_cci_val_emmitter[
    lr_type == "NO" &
    validation_set == "macro"
  ][frac < 0.2],
  aes(
    x = reorder(paste(dataset, tissue, EMITTER_CELLTYPE), frac),
    y = (1 - frac) * 100
  )
) + geom_point() + coord_flip()

ggplot(
  dt_cci_val_emmitter_reg[
    lr_type == "NO" &
    validation_set == "macro_up"
  ][frac < 0.1],
  aes(
    x = reorder(paste(dataset, tissue, EMITTER_CELLTYPE, REGULATION), frac),
    y = (1 - frac) * 100
  )
) + geom_point(
) + coord_flip(
)

ggplot(
  dt_cci_val_emmitter_reg[
    lr_type == "NO" &
    validation_set == "cardio" &
    tissue == "Heart" &
    #EMITTER_CELLTYPE == "pancreatic ductal cell" &
    dataset == "TMS FACS (female)"
  ],
  aes(
    x = reorder(paste(dataset, tissue, EMITTER_CELLTYPE, REGULATION), frac),
    y = (1 - frac) * 100
  )
) + geom_point() + coord_flip()


ggplot(
  dt_cci_val_receiver[
    lr_type == "NO" &
    validation_set == "pancreas_can_down"
  ][frac < 0.97],
  aes(
    x = reorder(paste(dataset, tissue, RECEIVER_CELLTYPE), frac),
    y = 1 - frac
  )
) + geom_point() + coord_flip()

names(val_all)
val_pancreas_lr_can_down


dt_cci_val_tissue[lr_type == "NO", .N, by = "validation_set"]

###########
dt_cci_rel[
  , {
    totwt = .N
    .SD[, .(frac = .N / totwt), by = cardio_hyp_up]
  },
  by = c("dataset", "tissue", "REGULATION")
][order(-frac)][cardio_hyp_up != "NO"][
  tissue == "Heart_and_Aorta"
]

dt_cci_rel[
  , {
    totwt = .N
    .SD[, .(frac = .N / totwt), by = macro]
  },
  by = c("dataset", "tissue", "EMITTER_CELLTYPE")
][order(-frac)][macro != "NO"][1:20]

