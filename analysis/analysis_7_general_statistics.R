####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## General statistics about all results
##
####################################################
##

## Keep only relevant CCIs for analysis ####

dt_cci_rel <- dt_cci_full[!grepl("mixed", dataset)]

dt_cci_rel[
  ,
  scDiffCom_regulation := ifelse(
    !REGULATION %in% c("UP", "DOWN"),
    "NO",
    REGULATION
  )
]

dt_cci_rel[
  dt_celltype_conversion[
    ,
    c("new_tissue", "cell_ontology_final", "cell_family")
  ],
  on = c("tissue==new_tissue", "EMITTER_CELLTYPE==cell_ontology_final"),
  EMITTER_CELLFAMILY := i.cell_family
]
dt_cci_rel[
  dt_celltype_conversion[
    ,
    c("new_tissue", "cell_ontology_final", "cell_family")
  ],
  on = c("tissue==new_tissue", "RECEIVER_CELLTYPE==cell_ontology_final"),
  RECEIVER_CELLFAMILY := i.cell_family
]
dt_cci_rel[
  ,
  ER_CELLFAMILY := paste(
    EMITTER_CELLFAMILY,
    RECEIVER_CELLFAMILY,
    sep = "_"
  )
]

dt_cci_rel[
  dt_celltype_conversion[
    ,
    c("new_tissue", "cell_ontology_final", "cell_family_mid_2")
  ],
  on = c("tissue==new_tissue", "EMITTER_CELLTYPE==cell_ontology_final"),
  EMITTER_CELLFAMILY_2 := i.cell_family_mid_2
]
dt_cci_rel[
  dt_celltype_conversion[
    ,
    c("new_tissue", "cell_ontology_final", "cell_family_mid_2")
  ],
  on = c("tissue==new_tissue", "RECEIVER_CELLTYPE==cell_ontology_final"),
  RECEIVER_CELLFAMILY_2 := i.cell_family_mid_2
]
dt_cci_rel[
  ,
  ER_CELLFAMILY_2 := paste(
    EMITTER_CELLFAMILY_2,
    RECEIVER_CELLFAMILY_2,
    sep = "_"
  )
]

dt_ora_rel <- dt_ora_full[!grepl("mixed", dataset)]

## Keep only relevant sex CCIs for analysis ####

dt_cci_sex <- dt_cci_full_diffsex[!grepl("combined", dataset)]

dt_cci_sex[
  ,
  scDiffCom_regulation := ifelse(
    !REGULATION %in% c("UP", "DOWN"),
    "NO",
    REGULATION
  )
]

dt_cci_sex[
  dt_celltype_conversion[
    ,
    c("new_tissue", "cell_ontology_final", "cell_family")
  ],
  on = c("tissue==new_tissue", "EMITTER_CELLTYPE==cell_ontology_final"),
  EMITTER_CELLFAMILY := i.cell_family
]
dt_cci_sex[
  dt_celltype_conversion[
    ,
    c("new_tissue", "cell_ontology_final", "cell_family")
  ],
  on = c("tissue==new_tissue", "RECEIVER_CELLTYPE==cell_ontology_final"),
  RECEIVER_CELLFAMILY := i.cell_family
]
dt_cci_sex[
  ,
  ER_CELLFAMILY := paste(
    EMITTER_CELLFAMILY,
    RECEIVER_CELLFAMILY,
    sep = "_"
  )
]

dt_cci_sex[
  dt_celltype_conversion[
    ,
    c("new_tissue", "cell_ontology_final", "cell_family_mid_2")
  ],
  on = c("tissue==new_tissue", "EMITTER_CELLTYPE==cell_ontology_final"),
  EMITTER_CELLFAMILY_2 := i.cell_family_mid_2
]
dt_cci_sex[
  dt_celltype_conversion[
    ,
    c("new_tissue", "cell_ontology_final", "cell_family_mid_2")
  ],
  on = c("tissue==new_tissue", "RECEIVER_CELLTYPE==cell_ontology_final"),
  RECEIVER_CELLFAMILY_2 := i.cell_family_mid_2
]
dt_cci_sex[
  ,
  ER_CELLFAMILY_2 := paste(
    EMITTER_CELLFAMILY_2,
    RECEIVER_CELLFAMILY_2,
    sep = "_"
  )
]

dt_ora_sex <- dt_ora_full_diffsex[!grepl("combined", dataset)]

## Full list of Seurat genes ####

seurat_genes <- lapply(
  seurats_analysis,
  rownames
)
seurat_genes$calico <- unique(
  c(seurat_genes$calico_kidney,
  seurat_genes$calico_lung,
  seurat_genes$calico_spleen
  )
)

## General LRI statistics ####

dt_lri_human
dt_lri_human[
  ,
  Type := ifelse(
    !is.na(LIGAND_2) | !is.na(RECEPTOR_2),
    "Complex",
    "Simple"
  )
]
table(dt_lri_human$Type)

dt_lri_mouse
dt_lri_mouse[
  ,
  Type := ifelse(
    !is.na(LIGAND_2) | !is.na(RECEPTOR_2),
    "Complex",
    "Simple"
  )
]
table(dt_lri_mouse$Type)

## Aging LRI statistics ####

top10_lr_aging <- pmid_aging_lri_n[order(-count)][1:10]$gene
top10_lr_aging
top10_lr_aging[top10_lr_aging %in% hagr_genes]

## Number of detected CCIs ####

# total
dt_cci_full[, .N]
dt_cci_rel[, .N]

# by datasets
dt_cci_full[, .N, by = c("dataset")]

# by regulation
dt_cci_full[, .N, by = c("REGULATION")]
dt_cci_full[, .N, by = c("REGULATION")][, N / sum(N) * 100]
dt_cci_rel[, .N, by = c("REGULATION")]
dt_cci_rel[, .N, by = c("REGULATION")][
  , N / sum(N) * 100
]

dt_cci_rel[ REGULATION == "NSC", .N, by = c("IS_CCI_DETECTED_YOUNG", "IS_CCI_DETECTED_OLD")]

dt_cci_rel[ REGULATION == "NSC", .N, by = c("IS_CCI_DETECTED_YOUNG", "IS_CCI_DETECTED_OLD")][
  , N / sum(N) * 100
]

dt_cci_rel[dataset == "TMS FACS (male)" & tissue == "Lung", .N]
dt_cci_rel[dataset == "TMS Droplet (male)" & tissue == "Lung", .N]

## Average LRI per ER_CELLTYPES ####

dt_cci_rel[
  ,
   DTER_CELLTYPES := paste(
    dataset,
    tissue,
    EMITTER_CELLTYPE,
    RECEIVER_CELLTYPE,
    sep = "_"
    )
]

NLRI_template <- rbindlist(
  lapply(
    unique(dt_cci_rel$dataset),
    function(data) {
      rbindlist(
        lapply(
          unique(dt_cci_rel[dataset == data]$tissue),
          function(tiss) {
            dt_temp <- CJ(
              unique(
                dt_cci_rel[
                  dataset == data & tissue == tiss
                ]$EMITTER_CELLTYPE
              ),
              unique(
                dt_cci_rel[
                  dataset == data & tissue == tiss
                ]$RECEIVER_CELLTYPE
              )
            )
            dt_temp[, dataset := data]
            dt_temp[, tissue := tiss]
            dt_temp[
              ,
              DTER_CELLTYPES := paste(
                dataset,
                tissue,
                V1,
                V2,
                sep = "_"
                )
            ]
          }
        )
      )
    }
  )
)
NLRI_template[, V1 := NULL]
NLRI_template[, V2 := NULL]

NLRI_template[
  dt_cci_rel[
    ,
    .N,
    by = c("dataset", "tissue", "DTER_CELLTYPES")
  ],
  on = c("dataset", "tissue", "DTER_CELLTYPES"),
  N := i.N
]
NLRI_template[is.na(NLRI_template)] <- 0
NLRI_template[
  ,
  dataset := factor(
    dataset,
    c(
      "TMS FACS (male)",
      "TMS FACS (female)",
      "TMS Droplet (male)",
      "TMS Droplet (female)",
      "Calico Droplet (male)")
  )
]

NLRI_template[, mean(N)]
NLRI_template[, sd(N)]

## LRI not detected at all and counts of detection for each LRI per dataset ####

LRI_template <- LRI_mouse$LRI_curated[, 1:6]

LRI_in_seurat <- lapply(
  seurat_genes,
  function(i) {
    LRI_mouse$LRI_curated[
      LIGAND_1 %in% i &
        RECEPTOR_1 %in% i &
        LIGAND_2 %in% c(i, NA) &
        RECEPTOR_2 %in% c(i, NA) &
        RECEPTOR_3 %in% c(i, NA)
    ][]
  }
)
LRI_in_seurat <- unique(
  rbindlist(LRI_in_seurat)
)

LRI_not_detected <- LRI_in_seurat[
  !LRI %in% unique(dt_cci_rel$LRI)
]

LRI_not_detected_genes <- sort(unique(
  unlist(
    LRI_not_detected[
      ,
      c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")
    ]
  )
))
LRI_not_detected_genes <- LRI_not_detected_genes[!is.na(LRI_not_detected_genes)]

## Compare SDEA with scDiffCom ####

dt_cci_rel[
  ,
  l1_deg :=  ifelse(
    L1_BH_P_VALUE_DE <= 0.05 & L1_LOGFC >= log(1.5),
    "UP",
    ifelse(
      L1_BH_P_VALUE_DE <= 0.05 & L1_LOGFC <= -log(1.5),
      "DOWN",
      "NOT DE"
    )
  )
]
dt_cci_rel[
  ,
  l2_deg :=  ifelse(
    is.na(L2_BH_P_VALUE_DE),
    NA,
    ifelse(
      L2_BH_P_VALUE_DE <= 0.05 & L2_LOGFC >= log(1.5),
      "UP",
      ifelse(
        L2_BH_P_VALUE_DE <= 0.05 & L2_LOGFC <= -log(1.5),
        "DOWN",
        "NOT DE"
      )
    )
  )
]
dt_cci_rel[
  ,
  r1_deg :=  ifelse(
    R1_BH_P_VALUE_DE <= 0.05 & R1_LOGFC >= log(1.5),
    "UP",
    ifelse(
      R1_BH_P_VALUE_DE <= 0.05 & R1_LOGFC <= -log(1.5),
      "DOWN",
      "NOT DE"
    )
  )
]
dt_cci_rel[
  ,
  r2_deg :=  ifelse(
    is.na(R2_BH_P_VALUE_DE),
    NA,
    ifelse(
      R2_BH_P_VALUE_DE <= 0.05 & R2_LOGFC >= log(1.5),
      "UP",
      ifelse(
        R2_BH_P_VALUE_DE <= 0.05 & R2_LOGFC <= -log(1.5),
        "DOWN",
        "NOT DE"
      )
    )
  )
]
dt_cci_rel[
  ,
  r3_deg :=  ifelse(
    is.na(R3_BH_P_VALUE_DE),
    NA,
    ifelse(
      R3_BH_P_VALUE_DE <= 0.05 & R3_LOGFC >= log(1.5),
      "UP",
      ifelse(
        R3_BH_P_VALUE_DE <= 0.05 & R3_LOGFC <= -log(1.5),
        "DOWN",
        "NOT DE"
      )
    )
  )
]

table(dt_cci_rel$l1_deg)
table(dt_cci_rel$l2_deg)
table(dt_cci_rel$r1_deg)
table(dt_cci_rel$r2_deg)
table(dt_cci_rel$r3_deg)

dt_cci_sdea_comp <- dt_cci_rel[
  ,
  c(
  "dataset",
  "scDiffCom_regulation",
  "l1_deg",
  "l2_deg",
  "r1_deg",
  "r2_deg",
  "r3_deg"
  )
]
dt_cci_sdea_comp[
  ,
  ligand_deg := ifelse(
    is.na(l2_deg),
    l1_deg,
    ifelse(
      l1_deg == l2_deg,
      l1_deg,
      ifelse(
        l1_deg == "NOT DE",
        l2_deg,
        ifelse(
          l2_deg == "NOT DE",
          l1_deg,
          "NOT DE"
        )
      )
    )
  )
]
dt_cci_sdea_comp[
  ,
  receptor_deg := ifelse(
    is.na(r2_deg),
    r1_deg,
    ifelse(
      r1_deg == r2_deg,
      r1_deg,
      ifelse(
        r1_deg == "NOT DE",
        r2_deg,
        ifelse(
          r2_deg == "NOT DE",
          r1_deg,
          "NOT DE"
        )
      )
    )
  )
]

table(dt_cci_sdea_comp$ligand_deg)
table(dt_cci_sdea_comp$receptor_deg)
anyNA(dt_cci_sdea_comp$receptor_deg)

dt_sdea_comp_count <- dt_cci_sdea_comp[
  ,
  .N,
  by = c(
    "scDiffCom_regulation",
    "ligand_deg",
    "receptor_deg"
  )
]
dt_sdea_comp_count[
  ,
  pct := N / sum(N) * 100#,
  #by = c("dataset")
]
setnames(
  dt_sdea_comp_count,
  old = colnames(dt_sdea_comp_count),
  new = c("scDiffCom_regulation", "ligand_SDEA_regulation",
  "receptor_SDEA_regulation", "Number_of_LRIs", "Pct_of_LRIs")
)

fwrite(
  dt_sdea_comp_count,
  paste0(
    path_scagecom_output,
    "Supplementary_Data_SDEA_comp.csv"
  )
)

# dt_sdea_comp_count_list <- lapply(
#   unique(dt_sdea_comp_count$dataset),
#   function(i) {
#     dt_sdea_comp_count[dataset == i]
#   }
# )
# names(dt_sdea_comp_count_list) <- unique(dt_sdea_comp_count$dataset)

## Distribution of OR GO Terms by GO levels ####

hist(
dt_ora_rel[
  ORA_CATEGORY == "GO_TERMS" & ORA_REGULATION %in% c("UP", "DOWN") #&
  #dataset == "TMS FACS (male)" & tissue == "Liver"
]$LEVEL, breaks = 100)
