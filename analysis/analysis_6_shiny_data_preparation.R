####################################################
##
## Project: scAgeCom
##
## Last update - May 2022
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## Process all results for the shiny app
##
####################################################
##

## Add libraries ####

library(rrvgo)
library(GOSemSim)

## Shortcut (read saved data) ####

shiny_list_full <- readRDS(
  paste0(
    path_scagecom_output,
    "scAgeComShiny_data.rds"
  )
)

shiny_dt_go_reduced <- readRDS(
  paste0(
    path_scagecom_output,
    "shiny_dt_go_reduced.rds"
  )
)

## Prepare LRI table for scAgeComShiny ####

shiny_dt_lri_mouse <- copy(dt_lri_mouse)
shiny_dt_lri_mouse[, SOURCE := gsub(";;", ";", SOURCE)]
shiny_dt_lri_mouse[, SOURCE := sub(";$", "", SOURCE)]
shiny_dt_lri_mouse[, SOURCE := sub("^;", "", SOURCE)]

shiny_cols_lri <- c(
  "LIGAND_1", "LIGAND_2",
  "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
  "DATABASE", "SOURCE"
)
shiny_dt_lri_mouse <- shiny_dt_lri_mouse[, shiny_cols_lri, with = FALSE]
setnames(
  shiny_dt_lri_mouse,
  old = shiny_cols_lri,
  new = c(
    "Ligand (1)", "Ligand (2)",
    "Receptor (1)", "Receptor (2)", "Receptor (3)",
    "Database(s) of Origin", "Retrieved Sources"
  )
)

## Process the full CCI table for shiny ####

shiny_dt_cci_full <- copy(dt_cci_full)

# round numeric values

shiny_dt_cci_full[, CCI_SCORE_YOUNG := signif(CCI_SCORE_YOUNG, 4)]
shiny_dt_cci_full[, CCI_SCORE_OLD := signif(CCI_SCORE_OLD, 4)]
shiny_dt_cci_full[, LOG2FC := signif(LOG2FC, 3)]
shiny_dt_cci_full[, BH_P_VALUE_DE := signif(BH_P_VALUE_DE, 3)]
shiny_dt_cci_full[, LOG2FC_L := signif(LOG2FC_L, 3)]
shiny_dt_cci_full[, LOG2FC_R := signif(LOG2FC_R, 3)]

# only keep relevant columns for shiny

colnames(shiny_dt_cci_full)
shiny_cci_cols_informative <- c(
  "dataset",
  "tissue",
  "LRI",
  "EMITTER_CELLTYPE",
  "RECEIVER_CELLTYPE",
  "LOG2FC",
  "BH_P_VALUE_DE",
  "REGULATION",
  "CCI_SCORE_YOUNG",
  "CCI_SCORE_OLD",
  "LOG2FC_L",
  "LOG2FC_R",
  "NCELLS_EMITTER_YOUNG",
  "NCELLS_EMITTER_OLD",
  "NCELLS_RECEIVER_YOUNG",
  "NCELLS_RECEIVER_OLD",
  "IS_CCI_EXPRESSED_YOUNG",
  "IS_CCI_EXPRESSED_OLD",
  "LIGAND_1",
  "LIGAND_2",
  "RECEPTOR_1",
  "RECEPTOR_2",
  "RECEPTOR_3",
  "GO_NAMES",
  "KEGG_NAMES"
)
shiny_new_cci_cols_informative <- c(
  "Dataset",
  "Tissue",
  "Ligand-Receptor Interaction",
  "Emitter Cell Type",
  "Receiver Cell Type",
  "Log2 FC",
  "Adj. p-value",
  "Age Regulation",
  "Young CCI Score",
  "Old CCI Score",
  "Ligand Log2 FC",
  "Receptor Log2 FC",
  "NCELLS_EMITTER_YOUNG",
  "NCELLS_EMITTER_OLD",
  "NCELLS_RECEIVER_YOUNG",
  "NCELLS_RECEIVER_OLD",
  "IS_CCI_EXPRESSED_YOUNG",
  "IS_CCI_EXPRESSED_OLD",
  "LIGAND_1",
  "LIGAND_2",
  "RECEPTOR_1",
  "RECEPTOR_2",
  "RECEPTOR_3",
  "GO_NAMES",
  "KEGG_NAMES"
)
shiny_dt_cci_full <- shiny_dt_cci_full[
  ,
  shiny_cci_cols_informative,
  with = FALSE
]
setnames(
  shiny_dt_cci_full,
  old = shiny_cci_cols_informative,
  new = shiny_new_cci_cols_informative
)

# set keys and order

setkey(shiny_dt_cci_full)
setorder(
  shiny_dt_cci_full,
  -`Age Regulation`,
  -`Log2 FC`,
  `Adj. p-value`,
  -`Old CCI Score`,
  -`Young CCI Score`
)

## Select datasets for shiny CCI results ####

shiny_dt_cci_full <- shiny_dt_cci_full[
  !grepl("mixed", Dataset)
]

## Create a table of counts summary for shiny #####

shiny_tissue_counts_summary <- dcast.data.table(
  shiny_dt_cci_full[
    ,
    .N,
    by = c("Dataset", "Tissue", "Age Regulation")
  ],
  Dataset + Tissue ~ `Age Regulation`,
  value.var = "N",
  fill = 0
)[
  ,
  "Total CCI" := DOWN + FLAT + NSC + UP
][
  shiny_dt_cci_full[
    ,
    list("NCT" = uniqueN(`Emitter Cell Type`)),
    by = c("Dataset", "Tissue")
  ],
  on = c("Dataset", "Tissue"),
  "Total Cell Types" := i.NCT
]

setnames(
  shiny_tissue_counts_summary,
  colnames(shiny_tissue_counts_summary),
  c(
    "Dataset", "Tissue",
    "Down CCIs",
    "Flat CCIs",
    "NSC CCIs",
    "UP CCIs",
    "Total Cell-Cell Interactions",
    "Total Cell Types"
  )
)

setcolorder(
  shiny_tissue_counts_summary,
  c(
    "Tissue",
    "Dataset",
    "Total Cell Types",
    "Total Cell-Cell Interactions",
    "Flat CCIs",
    "Down CCIs",
    "UP CCIs",
    "NSC CCIs"
  )
)

## Process full ORA table for shiny ####

shiny_dt_ora_full <- copy(dt_ora_full)

# only keep relevant columns

colnames(shiny_dt_ora_full)
shiny_ora_cols_informative <- c(
  "dataset",
  "tissue",
  "ORA_CATEGORY",
  "VALUE",
  "VALUE_BIS",
  "OR_UP",
  "BH_P_VALUE_UP",
  "ORA_SCORE_UP",
  "OR_DOWN",
  "BH_P_VALUE_DOWN",
  "ORA_SCORE_DOWN",
  "OR_FLAT",
  "BH_P_VALUE_FLAT",
  "ORA_SCORE_FLAT",
  "OR_DIFF",
  "BH_P_VALUE_DIFF",
  "ORA_SCORE_DIFF",
  "ASPECT",
  "LEVEL",
  "EMITTER_CELLTYPE",
  "RECEIVER_CELLTYPE",
  "IS_UP",
  "IS_DOWN",
  "IS_FLAT",
  "ORA_REGULATION"
)
shiny_new_ora_cols_informative <- c(
  "Dataset",
  "Tissue",
  "ORA_CATEGORY",
  "VALUE",
  "VALUE_BIS",
  "OR_UP",
  "BH_P_VALUE_UP",
  "ORA_SCORE_UP",
  "OR_DOWN",
  "BH_P_VALUE_DOWN",
  "ORA_SCORE_DOWN",
  "OR_FLAT",
  "BH_P_VALUE_FLAT",
  "ORA_SCORE_FLAT",
  "OR_DIFF",
  "BH_P_VALUE_DIFF",
  "ORA_SCORE_DIFF",
  "ASPECT",
  "GO Level",
  "EMITTER_CELLTYPE",
  "RECEIVER_CELLTYPE",
  "IS_UP",
  "IS_DOWN",
  "IS_FLAT",
  "ORA_REGULATION"
)
shiny_dt_ora_full <- shiny_dt_ora_full[
  ,
  shiny_ora_cols_informative,
  with = FALSE
]
setnames(
  shiny_dt_ora_full,
  old = shiny_ora_cols_informative,
  new = shiny_new_ora_cols_informative
)
setkey(shiny_dt_ora_full)

# round numeric values

shiny_dt_ora_full[, ORA_SCORE_UP := signif(ORA_SCORE_UP, 3)]
shiny_dt_ora_full[, ORA_SCORE_DOWN := signif(ORA_SCORE_DOWN, 3)]
shiny_dt_ora_full[, ORA_SCORE_FLAT := signif(ORA_SCORE_FLAT, 3)]
shiny_dt_ora_full[, OR_UP := signif(OR_UP, 3)]
shiny_dt_ora_full[, OR_DOWN := signif(OR_DOWN, 3)]
shiny_dt_ora_full[, OR_FLAT := signif(OR_FLAT, 3)]
shiny_dt_ora_full[, BH_P_VALUE_UP := signif(BH_P_VALUE_UP, 3)]
shiny_dt_ora_full[, BH_P_VALUE_DOWN := signif(BH_P_VALUE_DOWN, 3)]
shiny_dt_ora_full[, BH_P_VALUE_FLAT := signif(BH_P_VALUE_FLAT, 3)]

## Create a ORA GO reduced table for TreeMap ####

semd_BP <- GOSemSim::godata(
  OrgDb = "org.Mm.eg.db",
  ont = "BP"
)
semd_MF <- GOSemSim::godata(
  OrgDb = "org.Mm.eg.db",
  ont = "MF"
)
semd_CC <- GOSemSim::godata(
  OrgDb = "org.Mm.eg.db",
  ont = "CC"
)

semd_all_go <- unique(
  c(
    names(semd_BP@IC),
    names(semd_CC@IC),
    names(semd_MF@IC)
  )
)

# Warning: the computation goes over xxx*9 cases and takes more than one hour

shiny_dt_go_reduced <- rbindlist(
  l = lapply(
    setNames(
      sort(unique(shiny_dt_ora_full$Dataset)),
      sort(unique(shiny_dt_ora_full$Dataset))
    ),
    function(dataset) {
      print(dataset)
      dt <- copy(shiny_dt_ora_full)[
        Dataset == dataset &
          ORA_CATEGORY == "GO_TERMS" &
          VALUE_BIS %in% semd_all_go
      ]
      if (nrow(dt) == 0) return(NULL)
      rbindlist(
        l = lapply(
          setNames(
            sort(unique(dt$Tissue)),
            sort(unique(dt$Tissue))
          ),
          function(tissue) {
            print(tissue)
            rbindlist(
              l = lapply(
                setNames(
                  sort(unique(dt$ASPECT)),
                  sort(unique(dt$ASPECT))
                ),
                function(aspect) {
                  print(aspect)
                  rbindlist(
                    l = lapply(
                      list(UP = "UP", DOWN = "DOWN", FLAT = "FLAT"),
                      function(regulation) {
                        print(regulation)
                        dt_intern <- dt[Tissue == tissue & ASPECT == aspect]
                        if (regulation == "UP") {
                          dt_intern <- dt_intern[IS_UP == TRUE]
                          OR_intern <- dt_intern$OR_UP
                          BH_intern <- dt_intern$BH_P_VALUE_UP
                        } else if (regulation == "DOWN") {
                          dt_intern <- dt_intern[IS_DOWN == TRUE]
                          OR_intern <- dt_intern$OR_DOWN
                          BH_intern <- dt_intern$BH_P_VALUE_DOWN
                        } else if (regulation == "FLAT") {
                          dt_intern <- dt_intern[IS_FLAT == TRUE]
                          OR_intern <- dt_intern$OR_FLAT
                          BH_intern <- dt_intern$BH_P_VALUE_FLAT
                        } else {
                          stop("Error")
                        }
                        if (aspect == "biological_process") {
                          ont_intern <- "BP"
                          semd_intern <- semd_BP
                        } else if (aspect == "molecular_function") {
                          ont_intern <- "MF"
                          semd_intern <- semd_MF
                        } else if (aspect == "cellular_component") {
                          ont_intern <- "CC"
                          semd_intern <- semd_CC
                        }
                        if (nrow(dt_intern) == 0) return(NULL)
                        if (any(is.infinite(OR_intern))) {
                          if (all(is.infinite(OR_intern))) {
                            OR_intern <- rep(2, length(OR_intern))
                          } else {
                            max_finite <- max(OR_intern[is.finite(OR_intern)])
                            OR_intern[is.infinite(OR_intern)] <- max_finite
                          }
                        }
                        scores_intern <- -log10(BH_intern) * log2(OR_intern)
                        scores_intern <- setNames(
                          scores_intern,
                          dt_intern$VALUE_BIS
                        )
                        simMatrix <- calculateSimMatrix(
                          x = dt_intern$VALUE_BIS,
                          orgdb = "org.Mm.eg.db",
                          semdata = semd_intern,
                          ont = ont_intern,
                          method = "Rel"
                        )
                        if (is.null(dim(simMatrix))) return(NULL)
                        reducedTerms <- reduceSimMatrix(
                          simMatrix = simMatrix,
                          scores = scores_intern,
                          threshold = 0.7,
                          orgdb = "org.Mm.eg.db"
                        )
                        if (nrow(reducedTerms) == 0) return(NULL)
                        setDT(reducedTerms)
                        reducedTerms
                      }
                    ),
                    idcol = "REGULATION"
                  )
                }
              ),
              idcol = "ASPECT"
            )
          }
        ),
        idcol = "Tissue"
      )
    }
  ),
  idcol = "Dataset"
)

## Save ORA GO reduced table as it is time consuming to reproduce ####

saveRDS(
  shiny_dt_go_reduced,
  paste0(
    path_scagecom_output,
    "shiny_dt_go_reduced.rds"
  )
)

## Select datasets for ORA results ####

shiny_dt_ora_full <- shiny_dt_ora_full[
  !grepl("mixed", Dataset)
]
shiny_dt_go_reduced <- shiny_dt_go_reduced[
  !grepl("mixed", Dataset)
]

## Process the cross-tissue keyword counts table ####

shiny_dt_ora_key_counts <- copy(dt_ora_key_counts)

shiny_dt_ora_key_counts[
  data.table(
    old_category = c(
      "LRI",
      "LIGAND_COMPLEX",
      "RECEPTOR_COMPLEX",
      "ER_CELLTYPES",
      "EMITTER_CELLTYPE",
      "RECEIVER_CELLTYPE",
      "GO_TERMS",
      "KEGG_PWS",
      "ER_CELLFAMILIES",
      "EMITTER_CELLFAMILY",
      "RECEIVER_CELLFAMILY"
    ),
    new_category = c(
      "Ligand-Receptor Interaction",
      "Ligand",
      "Receptor",
      "Emitter-Receiver Cell Type",
      "Emitter Cell Type",
      "Receiver Cell Type",
      "GO Term",
      "KEGG Pathway",
      "Emitter-Receiver Cell Type Family",
      "Emitter Cell Type Family",
      "Receiver Cell Type Family"
    )
  ),
  on = c("ORA_CATEGORY==old_category"),
  ORA_CATEGORY := i.new_category
]

setcolorder(
  shiny_dt_ora_key_counts,
  c(
    "ORA_CATEGORY",
    "ORA_REGULATION",
    "VALUE",
    "Overall (Union)",
    "TMS FACS (male)",
    "TMS FACS (female)",
    "TMS Droplet (male)",
    "TMS Droplet (female)",
    "Calico Droplet (male)"
  )
)

shiny_dt_ora_key_counts[
  unique(
    dt_ora_full[
      ORA_CATEGORY == "GO_TERMS",
      c("VALUE", "ASPECT", "LEVEL")
    ]
  ),
  on = "VALUE",
  c("ASPECT", "GO Level") := mget(paste0("i.", c("ASPECT", "LEVEL")))
]

shiny_dt_ora_key_counts[, `GO Level` := factor(`GO Level`)]

setkey(shiny_dt_ora_key_counts)
setorder(
  shiny_dt_ora_key_counts,
  -`Overall (Union)`
)

shiny_dt_ora_key_counts[
  ,
  `Overall (Union)` := factor(
    ifelse(
      `Overall (Union)` < 10 & `Overall (Union)` > 0,
      paste0("0", `Overall (Union)`, "/23"),
      paste0(`Overall (Union)`, "/23")
    )
  ),
]
shiny_dt_ora_key_counts[
  ,
  `TMS FACS (male)` := factor(
    ifelse(
      `TMS FACS (male)` < 10 & `TMS FACS (male)` > 0,
      paste0("0", `TMS FACS (male)`, "/21"),
      paste0(`TMS FACS (male)`, "/21")
    )
  ),
]
shiny_dt_ora_key_counts[
  ,
  `TMS FACS (female)` := factor(
    ifelse(
      `TMS FACS (female)` < 10 & `TMS FACS (female)` > 0,
      paste0("0", `TMS FACS (female)`, "/19"),
      paste0(`TMS FACS (female)`, "/19")
    )
  ),
]
shiny_dt_ora_key_counts[
  ,
  `TMS Droplet (male)` := factor(
    paste0(`TMS Droplet (male)`, "/6")
  ),
]
shiny_dt_ora_key_counts[
  ,
  `TMS Droplet (female)` :=  factor(
    paste0(`TMS Droplet (female)`, "/9")
  ),
]
shiny_dt_ora_key_counts[
  ,
  `Calico Droplet (male)` := factor(
    paste0(`Calico Droplet (male)`, "/3")
  ),
]

shiny_dt_ora_key_counts[
  ,
  `TMS FACS (mixed)` :=  NULL
]
shiny_dt_ora_key_counts[
  ,
  `TMS Droplet (mixed)` :=  NULL
]

## Process the cross-tissue keyword summary table ####

shiny_dt_ora_key_summary <- copy(dt_ora_key_summary)
setnames(
  shiny_dt_ora_key_summary,
  old = c("tissue", "dataset", "dataset_tissue"),
  new = c("Tissue", "Dataset", "Dataset_Tissue")
)

shiny_dt_ora_key_summary[
  data.table(
    old_category = c(
      "LRI",
      "LIGAND_COMPLEX",
      "RECEPTOR_COMPLEX",
      "ER_CELLTYPES",
      "EMITTER_CELLTYPE",
      "RECEIVER_CELLTYPE",
      "GO_TERMS",
      "KEGG_PWS",
      "ER_CELLFAMILIES",
      "EMITTER_CELLFAMILY",
      "RECEIVER_CELLFAMILY"
    ),
    new_category = c(
      "Ligand-Receptor Interaction",
      "Ligand",
      "Receptor",
      "Emitter-Receiver Cell Type",
      "Emitter Cell Type",
      "Receiver Cell Type",
      "GO Term",
      "KEGG Pathway",
      "Emitter-Receiver Cell Type Family",
      "Emitter Cell Type Family",
      "Receiver Cell Type Family"
    )
  ),
  on = c("ORA_CATEGORY==old_category"),
  ORA_CATEGORY := i.new_category
]

shiny_dt_ora_key_summary[
  unique(
    dt_ora_full[
      ORA_CATEGORY == "GO_TERMS",
      c("VALUE", "ASPECT", "LEVEL")
    ]
  ),
  on = "VALUE",
  c("ASPECT", "GO Level") := mget(paste0("i.", c("ASPECT", "LEVEL")))
]
shiny_dt_ora_key_summary[, `GO Level` := factor(`GO Level`)]
shiny_dt_ora_key_summary <- shiny_dt_ora_key_summary[
  !grepl("mixed", Dataset)
]

shiny_dt_ora_key_template <- unique(
  shiny_dt_ora_key_summary[
    ORA_REGULATION != "No Data"
    ,
    c("Tissue", "Dataset")
  ]
)
shiny_dt_ora_key_template <- shiny_dt_ora_key_template[
  !grepl("mixed", Dataset)
]

## Create vectors to access categories in shiny ####

shiny_all_tissues <- sort(
  unique(
    shiny_dt_cci_full$Tissue
  )
)
names(shiny_all_tissues) <- shiny_all_tissues

shiny_all_celltypes <- shiny_dt_cci_full[
  ,
  list(
    CELLTYPE = unique(`Emitter Cell Type`)
  ),
  by = c("Dataset", "Tissue")
]

shiny_all_lris <- shiny_dt_cci_full[
  ,
  list(
    LRI = unique(`Ligand-Receptor Interaction`)
  ),
  by = c("Dataset", "Tissue")
]

shiny_all_genes <- rbindlist(
  lapply(
    scds_datasets,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            temp <- tissue@cci_table_detected
            cols_to_keep <- colnames(temp)[
              grepl("LIGAND|RECEPTOR", colnames(temp))
            ]
            data.table(
              GENE = sort(
                unique(
                  unlist(
                    temp[, cols_to_keep, with = FALSE]
                  )
                )
              )
            )
          }
        ),
        idcol = "Tissue"
      )
    }
  ),
  idcol = "Dataset"
)
shiny_all_genes <- shiny_all_genes[
  !grepl("mixed", Dataset)
]

shiny_all_go_terms <- rbindlist(
  lapply(
    shiny_all_tissues,
    function(tiss) {
      dt <- unique(shiny_dt_cci_full[Tissue == tiss, c("Dataset", "GO_NAMES")])
      datasets <- sort(unique(dt$Dataset))
      names(datasets) <- datasets
      rbindlist(
        lapply(
          datasets,
          function(dataset) {
            dt2 <- dt[Dataset == dataset]
            data.table(
              GO_NAMES = sort(
                unique(
                  unlist(
                    strsplit(
                      dt2$GO_NAMES,
                      ";"
                    )
                  )
                )
              )
            )
          }
        ),
        idcol = "Dataset"
      )
    }
  ),
  idcol = "Tissue"
)
shiny_all_go_terms <- shiny_all_go_terms[GO_NAMES != ""]

shiny_all_kegg_pws <- rbindlist(
  lapply(
    shiny_all_tissues,
    function(tiss) {
      dt <- unique(
        shiny_dt_cci_full[Tissue == tiss, c("Dataset", "KEGG_NAMES")]
      )
      datasets <- sort(unique(dt$Dataset))
      names(datasets) <- datasets
      rbindlist(
        lapply(
          datasets,
          function(dataset) {
            dt2 <- dt[Dataset == dataset]
            data.table(
              KEGG_NAMES = sort(
                unique(
                  unlist(
                    strsplit(
                      dt2$KEGG_NAMES,
                      ";"
                    )
                  )
                )
              )
            )
          }
        ),
        idcol = "Dataset"
      )
    }
  ),
  idcol = "Tissue"
)
shiny_all_kegg_pws <- shiny_all_kegg_pws[KEGG_NAMES != ""]
shiny_all_kegg_pws <- shiny_all_kegg_pws[KEGG_NAMES != "NA"]

shiny_abbr_celltype <- lapply(
  scd_dataset_names[
    !grepl("mixed", scd_dataset_names)
  ],
  function(dataset) {
    cell_types <- unique(
      shiny_dt_cci_full[Dataset == dataset]$`Emitter Cell Type`
    )
    dt <- dt_celltype_conversion[
      cell_ontology_final %in% cell_types,
      c("cell_ontology_final", "cell_abbreviation")
    ]
    setnames(
      dt,
      old = colnames(dt),
      new = c("ORIGINAL_CELLTYPE", "ABBR_CELLTYPE")
    )
    unique(dt)
  }
)

shiny_all_ora_categories_specific <- c(
  "By Cell Types",
  "By GO/KEGG",
  "By Genes"
)

shiny_all_ora_types <- c(
  "UP",
  "DOWN",
  "FLAT"
)

shiny_all_ora_go_aspects <- c(
  "Biological Process",
  "Molecular Function",
  "Cellular Component"
)

shiny_ora_categories_global <- c(
  "By GO/KEGG",
  "By Genes",
  "By Cell Type Families"
)

shiny_ora_categories_keyword <- c(
  "Ligand-Receptor Interaction",
  "Ligand",
  "Receptor",
  "GO Term",
  "KEGG Pathway",
  "Emitter-Receiver Cell Type Family",
  "Emitter Cell Type Family",
  "Receiver Cell Type Family"
)

## Remove heavy GO/KEGG column ####

shiny_dt_cci_full[, GO_NAMES := NULL]
shiny_dt_cci_full[, KEGG_NAMES := NULL]

## Combine all shiny data and save them ####

shiny_list_full <- list(
  LRI_mouse_curated = shiny_dt_lri_mouse,
  LRI_DATABASES = lri_dbs,
  CCI_table = shiny_dt_cci_full,
  ORA_table = shiny_dt_ora_full,
  TISSUE_COUNTS_SUMMARY = shiny_tissue_counts_summary,
  GO_REDUCED_table = shiny_dt_go_reduced,
  ALL_TISSUES = shiny_all_tissues,
  ALL_CELLTYPES = shiny_all_celltypes,
  ALL_LRIs = shiny_all_lris,
  ALL_GENES = shiny_all_genes,
  ALL_GO_TERMS = shiny_all_go_terms,
  ALL_KEGG_PWS = shiny_all_kegg_pws,
  ALL_ORA_CATEGORIES_SPECIFIC = shiny_all_ora_categories_specific,
  ALL_ORA_GO_ASPECTS = shiny_all_ora_go_aspects,
  ALL_ORA_TYPES = shiny_all_ora_types,
  ABBR_CELLTYPE = shiny_abbr_celltype,
  REFERENCE_GO_TERMS = copy(scDiffCom::LRI_mouse$LRI_curated_GO)[
    unique(scDiffCom::gene_ontology_level[, c("ID", "NAME")]),
    on = "GO_ID==ID",
    GO_NAME := i.NAME
  ][, c(1, 3)],
  REFERENCE_KEGG_PWS = scDiffCom::LRI_mouse$LRI_curated_KEGG[, c(1, 3)],
  ORA_KEYWORD_COUNTS = shiny_dt_ora_key_counts,
  ORA_KEYWORD_SUMMARY = shiny_dt_ora_key_summary,
  ORA_KEYWORD_TEMPLATE = shiny_dt_ora_key_template,
  ALL_ORA_CATEGORIES_GLOBAL = shiny_ora_categories_global,
  ALL_ORA_CATEGORIES_KEYWORD = shiny_ora_categories_keyword
)

saveRDS(
  shiny_list_full,
  paste0(
    path_scagecom_output,
    "scAgeComShiny_data.rds"
  )
)
