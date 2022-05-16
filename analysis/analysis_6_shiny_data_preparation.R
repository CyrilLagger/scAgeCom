####################################################
##
## Project: scAgeCom
##
## Last update - June 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## collect all results for shiny
##
####################################################
##

## Libraries ####

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

shiny_lri_dbs <- sort(
  unique(
    unlist(
      strsplit(
        shiny_dt_lri_mouse$`Database(s) of Origin`,
        ";")
    )
  )
)

shiny_dt_lri_mouse[
  ,
  Type := ifelse(
    !is.na(`Ligand (2)`) | !is.na(`Receptor (2)`),
    "Complex",
    "Simple"
  )
]
shiny_dt_lri_mouse[
  ,
  c(shiny_lri_dbs) := lapply(
    shiny_lri_dbs,
    function(i) {
      ifelse(grepl(i, `Database(s) of Origin`), TRUE, FALSE)
    }
  )
]

## Save LRI table for scAgeComShiny ####

shiny_ls_a2 <- list(
  shiny_dt_lri_mouse = shiny_dt_lri_mouse,
  shiny_lri_dbs = shiny_lri_dbs
)

saveRDS(
  shiny_ls_a2,
  paste0(
    path_scagecom_output,
    "shiny_ls_a2"
  )
)

## Create a full CCI table for shiny ####

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

## Create a full ORA table for shiny ####

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
  "RECEIVER_CELLTYPE"
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
  "RECEIVER_CELLTYPE"
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

# add ORA regulation annotations
shiny_dt_ora_full[
  ,
  ':='(
    IS_UP = ifelse(
      OR_UP >= 1 & BH_P_VALUE_UP <= 0.05,
      TRUE,
      FALSE
    ),
    IS_DOWN = ifelse(
      OR_DOWN >= 1 & BH_P_VALUE_DOWN <= 0.05,
      TRUE,
      FALSE
    ),
    IS_FLAT = ifelse(
      OR_FLAT >= 1 & BH_P_VALUE_FLAT <= 0.05,
      TRUE,
      FALSE
    )
  )
]

shiny_dt_ora_full[
  ,
  ORA_REGULATION := ifelse(
    !IS_UP & !IS_DOWN & !IS_FLAT,
    "Not Over-represented",
    ifelse(
      !IS_UP & !IS_DOWN & IS_FLAT,
      "FLAT",
      ifelse(
        !IS_UP & IS_DOWN & !IS_FLAT,
        "DOWN",
        ifelse(
          IS_UP & !IS_DOWN & !IS_FLAT,
          "UP",
          ifelse(
            IS_UP & !IS_DOWN & IS_FLAT,
            "UP",
            ifelse(
              !IS_UP & IS_DOWN & IS_FLAT,
              "DOWN",
              ifelse(
                IS_UP & IS_DOWN & !IS_FLAT,
                "UP:DOWN",
                "UP:DOWN"
              )
            )
          )
        )
      )
    )
  )
]

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

saveRDS(
  shiny_dt_go_reduced,
  "../data_scAgeCom/analysis/outputs_data/data_4_shiny_dt_go_reduced_newINF.rds"
)

# treemapPlot(
#   shiny_dt_go_reduced[
#     Dataset == "TMS FACS (male)" &
#       Tissue == "Kidney" &
#       ASPECT == "biological_process" &
#       REGULATION == "UP"
#   ][, -c(1,2,3,4)]
# )


## create vectors to access categories in shiny ####

ALL_TISSUES <- sort(
  unique(
    CCI_table$Tissue
  )
)

ALL_CELLTYPES <- CCI_table[
  ,
  list(
    CELLTYPE = unique(`Emitter Cell Type`)
  ),
  by = c("Dataset", "Tissue")
]

ALL_LRIs <- CCI_table[
  ,
  list(
    LRI = unique(`Ligand-Receptor Interaction`)
  ),
  by = c("Dataset", "Tissue")
]

ALL_GENES <- rbindlist(
  lapply(
    scds_datasets,
    function(dataset){
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            temp <- tissue@cci_table_detected
            cols_to_keep <- colnames(temp)[
              grepl("LIGAND|RECEPTOR", colnames(temp))
            ]
            data.table(
              GENE =sort(
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

names(ALL_TISSUES) <- ALL_TISSUES

ALL_GO_TERMS <- rbindlist(
  lapply(
    ALL_TISSUES,
    function(tiss) {
      dt <- unique(CCI_table[Tissue == tiss, c("Dataset", "GO_NAMES")])
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
ALL_GO_TERMS <- ALL_GO_TERMS[GO_NAMES != ""]

ALL_KEGG_PWS <- rbindlist(
  lapply(
    ALL_TISSUES,
    function(tiss) {
      dt <- unique(CCI_table[Tissue == tiss, c("Dataset", "KEGG_NAMES")])
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
ALL_KEGG_PWS <- ALL_KEGG_PWS[KEGG_NAMES != ""]
ALL_KEGG_PWS <- ALL_KEGG_PWS[KEGG_NAMES != "NA"]

ABBR_CELLTYPE <- lapply(
  scd_dataset_names,
  function(dataset){
    cell_types <- unique(
      CCI_table[Dataset == dataset]$`Emitter Cell Type`
    )
    dt <- meta_data_cell_types[
      Final_annotation %in% cell_types,
      c("Final_annotation", "Abbreviation")
    ]
    setnames(
      dt,
      old = colnames(dt),
      new = c("ORIGINAL_CELLTYPE", "ABBR_CELLTYPE")
    )
    unique(dt)
  }
)

ALL_ORA_CATEGORIES_SPECIFIC <- c(
  "By Cell Types",
  "By GO/KEGG",
  "By Genes"
)

ALL_ORA_TYPES <- c(
  "UP", 
  "DOWN",
  "FLAT"
)

ALL_ORA_GO_ASPECTS <- c(
  "Biological Process",
  "Molecular Function",
  "Cellular Component"
)

## removed heavly GO/KEGG column

CCI_table[, GO_NAMES := NULL]
CCI_table[, KEGG_NAMES := NULL]
colnames(CCI_table)

## save all results for shiny ####

data_4_tissue_specific_results <- list(
  CCI_table = CCI_table,
  ORA_table = ORA_table,
  shiny_tissue_counts_summary = shiny_tissue_counts_summary,
  shiny_dt_go_reduced = shiny_dt_go_reduced,
  ALL_TISSUES = ALL_TISSUES,
  ALL_CELLTYPES = ALL_CELLTYPES,
  ALL_LRIs = ALL_LRIs,
  ALL_GENES = ALL_GENES,
  ALL_GO_TERMS = ALL_GO_TERMS,
  ALL_KEGG_PWS = ALL_KEGG_PWS,
  ALL_ORA_CATEGORIES_SPECIFIC = ALL_ORA_CATEGORIES_SPECIFIC,
  ALL_ORA_GO_ASPECTS = ALL_ORA_GO_ASPECTS,
  ALL_ORA_TYPES = ALL_ORA_TYPES,
  ABBR_CELLTYPE = ABBR_CELLTYPE,
  REFERENCE_GO_TERMS = scDiffCom::LRI_mouse$LRI_curated_GO[, c(1, 3)],
  REFERENCE_KEGG_PWS = scDiffCom::LRI_mouse$LRI_curated_KEGG[, c(1, 3)]
)

saveRDS(
  data_4_tissue_specific_results,
  "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds"
)




## save all previous results ####

scAgeCom_shiny_data <- c(
  readRDS(
    "../data_scAgeCom/analysis/outputs_data/data_2_LRI_data_preparation.rds"
  ),
  readRDS(
    "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds"
  ),
  readRDS(
    "../data_scAgeCom/analysis/outputs_data/data_5_tissue_shared_results.rds"
  )
)

saveRDS(
  scAgeCom_shiny_data,
  "../data_scAgeCom/analysis/outputs_data/scAgeCom_shiny_data.rds"
)

