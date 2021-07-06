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
## preparation of tissue-specific results
##
####################################################
##

## Libraries ####

library(scDiffCom)
library(data.table)
library(rrvgo)
library(openxlsx)

## relative path of scDiffCom results ####

scAgeCom_path <- "../data_scAgeCom/scDiffCom_results_29_04_2021"

# it contains 5 datasets

dataset_paths <- list.dirs(scAgeCom_path, recursive = FALSE)
dataset_paths

dataset_names <- c(
  "Calico Droplet (male)",
  "TMS Droplet (female)",
  "TMS Droplet (male)",
  "TMS FACS (female)",
  "TMS FACS (male)"
)
names(dataset_names) <- dataset_names

## load meta.data cell type families annotation #####

meta_data_cell_types <- setDT(
  read.csv(
    "../data_scAgeCom/analysis/inputs_data/scDiffCom_cell_types_clean.csv",
    stringsAsFactors = FALSE
  )
)

meta_data_cell_types <- meta_data_cell_types[
  ,
  c("Tissue", "Final_annotation", "Family_broad", "Abbreviation")
]

meta_data_cell_types[Tissue == "BAT"]$Tissue <- "Adipose_Brown"
meta_data_cell_types[Tissue == "GAT"]$Tissue <- "Adipose_Gonadal"
meta_data_cell_types[Tissue == "MAT"]$Tissue <- "Adipose_Mesenteric"
meta_data_cell_types[Tissue == "SCAT"]$Tissue <- "Adipose_Subcutaneous"
meta_data_cell_types[Tissue == "Brain_Myeloid"]$Tissue <- "Brain"
meta_data_cell_types[Tissue == "Brain_Non-Myeloid"]$Tissue <- "Brain"

## clean and process results ####

process_dataset <- function(
  dataset_path
) {
  # retrieve and load each object in the dataset
  tissues <- gsub(".*scdiffcom_(.+)\\.rds.*", "\\1", list.files(dataset_path))
  dataset <- lapply(
    X = tissues,
    FUN = function(
      tiss
    ) {
      res <- readRDS(paste0(dataset_path, "/scdiffcom_", tiss, ".rds"))
    }
  )
  # change some tissue names
  tissues[tissues == "BAT"] <- "Adipose_Brown"
  tissues[tissues == "GAT"] <- "Adipose_Gonadal"
  tissues[tissues == "MAT"] <- "Adipose_Mesenteric"
  tissues[tissues == "SCAT"] <- "Adipose_Subcutaneous"
  dataset <- lapply(
    seq_along(dataset), 
    function (i) {
      dataset[[i]]@parameters$object_name <- tissues[[i]]
      dataset[[i]]
    }
  )
  names(dataset) <- tissues
  # erase raw CCI to make the object lighter
  dataset <- lapply(
    dataset,
    function(i) {
      new_object <- EraseRawCCI(i)
    }
  )
  # add cell type families ORA
  dataset <- lapply(
    dataset,
    function(i) {
      print(GetParameters(i)$object_name)
      temp_cell_types <- unique(
        GetTableCCI(
          object = i,
          type = "detected",
          simplified = FALSE
        )$EMITTER_CELLTYPE
      )
      temp_meta_family <- unique(
        meta_data_cell_types[
          Tissue == GetParameters(i)$object_name &
            Final_annotation %in% temp_cell_types
        ][, c(2,3)]
      )
      EMITTER_dt <- copy(temp_meta_family)
      setnames(
        EMITTER_dt,
        old = colnames(EMITTER_dt),
        new = c("EMITTER_CELLTYPE", "EMITTER_CELLFAMILY")
      )
      RECEIVER_dt <- copy(temp_meta_family)
      setnames(
        RECEIVER_dt,
        old = colnames(RECEIVER_dt),
        new = c("RECEIVER_CELLTYPE", "RECEIVER_CELLFAMILY")
      )
      ER_dt <- CJ(
        EMITTER_TYPE = unique(EMITTER_dt$EMITTER_CELLTYPE),
        RECEIVER_TYPE = unique(RECEIVER_dt$RECEIVER_CELLTYPE)
      )
      ER_dt[
        EMITTER_dt,
        on = "EMITTER_TYPE==EMITTER_CELLTYPE",
        EMITTER_FAMILY := i.EMITTER_CELLFAMILY
      ]
      ER_dt[
        RECEIVER_dt,
        on = "RECEIVER_TYPE==RECEIVER_CELLTYPE",
        RECEIVER_FAMILY := i.RECEIVER_CELLFAMILY
      ]
      ER_dt[, ER_CELLTYPES := paste(
        EMITTER_TYPE,
        RECEIVER_TYPE,
        sep = "_"
      )]
      ER_dt[, ER_CELLFAMILIES := paste(
        EMITTER_FAMILY,
        RECEIVER_FAMILY,
        sep = "_"
      )]
      ER_dt <- ER_dt[, c("ER_CELLTYPES", "ER_CELLFAMILIES")]
      print(EMITTER_dt)
      new_object <- RunORA(
        object = i,
        categories = c(
          "LRI", "LIGAND_COMPLEX", "RECEPTOR_COMPLEX", "ER_CELLTYPES",
          "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "GO_TERMS", "KEGG_PWS"
        ),
        extra_annotations = list(
          EMITTER_dt,
          RECEIVER_dt,
          ER_dt
        ),
        overwrite = FALSE,
        verbose = TRUE
      )
    }
  )
  return(dataset)
}

dataset_processed <- lapply(
  dataset_paths,
  function(path) {
    process_dataset(path)
  }
)
names(dataset_processed) <- dataset_names

#saveRDS(dataset_processed, "../data_scAgeCom/analysis/outputs_data/scAgeCom_results_processed_29_04_2021.rds")
dataset_processed <- readRDS("../data_scAgeCom/analysis/outputs_data/scAgeCom_results_processed_29_04_2021.rds")

## create a single CCI table for shiny ####

# create full table
CCI_table <- rbindlist(
  lapply(
    dataset_processed,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            GetTableCCI(
              object = tissue,
              type = "detected",
              simplified = FALSE
            )
          }
        ),
        fill = TRUE,
        idcol = "Tissue"
      )
    }
  ),
  idcol = "Dataset"
)

# add log2fc on top of logfc and deal with infinite values
CCI_table[, LOG2FC_BASE := LOGFC*log2(exp(1))]
CCI_table[
  ,
  LOG2FC := {
    temp <- LOG2FC_BASE
    temp_max <- ceiling(max(temp[is.finite(temp)]))
    temp_min <- floor(min(temp[is.finite(temp)]))
    temp_max <- max(temp_max, -temp_min)
    temp_min <- min(-temp_max, temp_min)
    ifelse(
      is.infinite(LOG2FC_BASE) & LOG2FC_BASE > 0,
      temp_max,
      ifelse(
        is.infinite(LOG2FC_BASE) & LOG2FC_BASE < 0,
        temp_min,
        LOG2FC_BASE
      )
    )},
  by = c(
    "Dataset",
    "Tissue"
  )
]

# add ligand (resp. receptor) log2fc
CCI_table[
  ,
  LOG2FC_L := log2(
    pmin(L1_EXPRESSION_OLD, L2_EXPRESSION_OLD, na.rm = TRUE)
    /
      pmin(L1_EXPRESSION_YOUNG, L2_EXPRESSION_YOUNG, na.rm = TRUE)
  )
]
CCI_table[
  ,
  LOG2FC_R := log2(
    pmin(R1_EXPRESSION_OLD, R2_EXPRESSION_OLD, R3_EXPRESSION_OLD, na.rm = TRUE)
    /
      pmin(R1_EXPRESSION_YOUNG, R2_EXPRESSION_YOUNG, R3_EXPRESSION_YOUNG, na.rm = TRUE)
  )
]
CCI_table[
  ,
  LOG2FC_L := {
    max_L <- ceiling(max(.SD[is.finite(LOG2FC_L)][["LOG2FC_L"]]))
    min_L <- floor(min(.SD[is.finite(LOG2FC_L)][["LOG2FC_L"]]))
    max_L <- max(max_L, -min_L)
    min_L <- min(-max_L, min_L)
    ifelse(
      is.infinite(LOG2FC_L) & LOG2FC_L > 0,
      max_L,
      ifelse(
        is.infinite(LOG2FC_L) & LOG2FC_L < 0,
        min_L,
        LOG2FC_L
      )
    )
  },
  by = c(
    "Dataset",
    "Tissue"
  )
]
CCI_table[
  ,
  LOG2FC_R := {
    max_R <- ceiling(max(.SD[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
    min_R <- floor(min(.SD[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
    max_R <- max(max_R, -min_R)
    min_R <- min(-max_R, min_R)
    ifelse(
      is.infinite(LOG2FC_R) & LOG2FC_R > 0,
      max_R,
      ifelse(
        is.infinite(LOG2FC_R) & LOG2FC_R < 0,
        min_R,
        LOG2FC_R
      )
    )
  },
  by = c(
    "Dataset",
    "Tissue"
  )
]

# add all GO names attached to each LRI 

CCI_table[
  dcast.data.table(
    scDiffCom::LRI_mouse$LRI_curated_GO[, c(1, 3)],
    LRI ~ .,
    value.var = "GO_NAME",
    fun.aggregate = paste0,
    collapse = ";"
  ),
  on = "LRI",
  GO_NAMES := i..
][
  ,
  GO_NAMES := paste0(";", GO_NAMES, ";")
]

# add all KEGG names attached to each LRI 

CCI_table[
  dcast.data.table(
    scDiffCom::LRI_mouse$LRI_curated_KEGG[, c(1, 3)],
    LRI ~ .,
    value.var = "KEGG_NAME",
    fun.aggregate = paste0,
    collapse = ";"
  ),
  on = "LRI",
  KEGG_NAMES := i..
][
  ,
  KEGG_NAMES := paste0(";", KEGG_NAMES, ";")
]

# round numeric values
CCI_table[, CCI_SCORE_YOUNG := signif(CCI_SCORE_YOUNG, 4)]
CCI_table[, CCI_SCORE_OLD := signif(CCI_SCORE_OLD, 4)]
CCI_table[, LOG2FC := signif(LOG2FC, 3)]
CCI_table[, BH_P_VALUE_DE := signif(BH_P_VALUE_DE, 3)]
CCI_table[, LOG2FC_L := signif(LOG2FC_L, 3)]
CCI_table[, LOG2FC_R := signif(LOG2FC_R, 3)]

# only keep relevant columns 
colnames(CCI_table)
CCI_cols_informative <- c(
  "Dataset",
  "Tissue",
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
CCI_cols_informative_shiny <- c(
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
CCI_table <- CCI_table[, CCI_cols_informative, with = FALSE]
setnames(
  CCI_table,
  old = CCI_cols_informative,
  new = CCI_cols_informative_shiny
)

# set keys and order
setkey(CCI_table)
setorder(
  CCI_table,
  -`Age Regulation`,
  -`Log2 FC`,
  `Adj. p-value`,
  -`Old CCI Score`,
  -`Young CCI Score`
)

## create a table of counts summary for shiny #####

TISSUE_COUNTS_SUMMARY <- dcast.data.table(
  CCI_table[
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
  CCI_table[
    ,
    list( "NCT" = uniqueN(`Emitter Cell Type`)),
    by = c("Dataset", "Tissue")
  ],
  on = c("Dataset", "Tissue"),
  "Total Cell Types" := i.NCT
]

setnames(
  TISSUE_COUNTS_SUMMARY,
  colnames(TISSUE_COUNTS_SUMMARY),
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
  TISSUE_COUNTS_SUMMARY,
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

## create a single ORA table for shiny ####

# create full table
ORA_table <- rbindlist(
  lapply(
    dataset_processed,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            rbindlist(
              GetTableORA(
                object = tissue,
                categories = "all",
                simplified = FALSE
              ),
              idcol = "ORA_CATEGORY",
              fill = TRUE
            )
          }
        ),
        fill = TRUE,
        idcol = "Tissue"
      )
    }
  ),
  idcol = "Dataset"
)

# add cell types for ERI
ORA_table[
  ORA_CATEGORY == "ER_CELLTYPES",
  c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE") := list(
    sub("_.*", "", VALUE),
    sub(".*_", "", VALUE)
  )]

# only keep relevant columns 
colnames(ORA_table)
ORA_cols_informative <- c(
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
  "LEVEL",
  "EMITTER_CELLTYPE",
  "RECEIVER_CELLTYPE"
)
ORA_cols_informative_shiny <- c(
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
ORA_table <- ORA_table[, ORA_cols_informative, with = FALSE]
setnames(
  ORA_table,
  old = ORA_cols_informative,
  new = ORA_cols_informative_shiny
)
setkey(ORA_table)

# round numeric values
ORA_table[, ORA_SCORE_UP := signif(ORA_SCORE_UP, 3)]
ORA_table[, ORA_SCORE_DOWN := signif(ORA_SCORE_DOWN, 3)]
ORA_table[, ORA_SCORE_FLAT := signif(ORA_SCORE_FLAT, 3)]
ORA_table[, OR_UP := signif(OR_UP, 3)]
ORA_table[, OR_DOWN := signif(OR_DOWN, 3)]
ORA_table[, OR_FLAT := signif(OR_FLAT, 3)]
ORA_table[, BH_P_VALUE_UP := signif(BH_P_VALUE_UP, 3)]
ORA_table[, BH_P_VALUE_DOWN := signif(BH_P_VALUE_DOWN, 3)]
ORA_table[, BH_P_VALUE_FLAT := signif(BH_P_VALUE_FLAT, 3)]

# add ORA regulation annotations
ORA_table[
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

ORA_table[
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

## create a ORA GO reduced table for TreeMap ####

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

all_org_go <- unique(c(names(semd_BP@IC), names(semd_CC@IC), names(semd_MF@IC)))

# Warning: the computation goes over 58*9 cases and takes more than one hour
GO_REDUCED_table <- rbindlist(
  l = lapply(
    setNames(
      sort(unique(ORA_table$Dataset)),
      sort(unique(ORA_table$Dataset))
    ),
    function(dataset) {
      print(dataset)
      dt <- copy(ORA_table)[
        Dataset == dataset &
          ORA_CATEGORY == "GO_TERMS" & 
          VALUE_BIS %in% all_org_go
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
  GO_REDUCED_table,
  "../data_scAgeCom/analysis/outputs_data/data_4_GO_REDUCED_table_newINF.rds"
)

# treemapPlot(
#   GO_REDUCED_table[
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
    dataset_processed,
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
  dataset_names,
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
  TISSUE_COUNTS_SUMMARY = TISSUE_COUNTS_SUMMARY,
  GO_REDUCED_table = GO_REDUCED_table,
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

## Prepare meta.data for summary picture

# scAgeCom_md_path <- "../data_scAgeCom/scDiffCom_results_md_25_04_2021/"
# 
# md_paths <- list.dirs(scAgeCom_md_path, recursive = FALSE)
# md_paths
# 
# md_names <- c(
#   "Calico Droplet (male)",
#   "TMS Droplet (female)",
#   "TMS Droplet (male)",
#   "TMS FACS (female)",
#   "TMS FACS (male)"
# )
# names(md_names) <- md_names
# 
# 
# process_md <- function(
#   md_path
# ) {
#   tissues <- gsub(".*md_(.+)\\.rds.*", "\\1", list.files(md_path))
#   md <- lapply(
#     X = tissues,
#     FUN = function(
#       tiss
#     ) {
#       res <- readRDS(paste0(md_path, "/md_", tiss, ".rds"))
#     }
#   )
#   # change some tissue names
#   tissues[tissues == "BAT"] <- "Adipose_Brown"
#   tissues[tissues == "GAT"] <- "Adipose_Gonadal"
#   tissues[tissues == "MAT"] <- "Adipose_Mesenteric"
#   tissues[tissues == "SCAT"] <- "Adipose_Subcutaneous"
#   names(md) <- tissues
#   return(md)
# }
# 
# test <- process_md(md_paths[[1]])

## Prepare Figure 3 for the manuscript ####

#retrieve the age information
md_paths <- list.dirs("../data_scAgeCom/scDiffCom_md_29_04_2021", recursive = FALSE)
md_paths

process_md <- function(
  md_path
) {
  # retrieve and load each object in the dataset
  tissues <- gsub(".*md_(.+)\\.rds.*", "\\1", list.files(md_path))
  mds <- lapply(
    X = tissues,
    FUN = function(
      tiss
    ) {
      res <- readRDS(paste0(md_path, "/md_", tiss, ".rds"))
      sort(unique(res$age))
    }
  )
  # change some tissue names
  tissues[tissues == "BAT"] <- "Adipose_Brown"
  tissues[tissues == "GAT"] <- "Adipose_Gonadal"
  tissues[tissues == "MAT"] <- "Adipose_Mesenteric"
  tissues[tissues == "SCAT"] <- "Adipose_Subcutaneous"
  names(mds) <- tissues
  return(mds)
}

md_processed <- lapply(
  md_paths,
  function(path) {
    process_md(path)
  }
)
names(md_processed) <- dataset_names

md_age <- rbindlist(
  lapply(
    md_processed,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            data.table(
              age_group = tissue
            )
          }
        ),
        idcol = "Tissue"
      )
    }
  ),
  idcol = "Dataset"
)

md_age_dc <- dcast.data.table(
  md_age,
  Dataset + Tissue ~ age_group,
  value.var = "age_group"
)

md_age_dc[
  ,
  age_group := ifelse(
    !is.na(young),
    "7-8m vs 22-23m",
    ifelse(
      is.na(`18m`),
      "3m vs 24m",
      ifelse(
        is.na(`21m`) & is.na(`24m`),
        "3m vs 18m",
        ifelse(
          is.na(`21m`),
          "3m vs 18/24m",
          "3m vs 18/21m"
        )
      )
    )
  )
]

md_age_dc[
  ,
  age_group2 := ifelse(
    Dataset == "TMS FACS (female)" & Tissue == "Mammary_Gland",
    "(3m vs 18/21m)",
    ifelse(
      Dataset == "TMS Droplet (male)" & Tissue %in% c("Liver", "Spleen"),
      "(3m vs 24m)",
      ifelse(
        Dataset == "TMS Droplet (male)" & Tissue %in% c("Lung"),
        "(3m vs 18m)",
        ""
      )
    )
  )
]


fig3_data <- readRDS(
  "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds"
)

dataset_summary <- fig3_data$CCI_table[
  ,
  c("Dataset", "Tissue", "Emitter Cell Type")
][
  ,
  uniqueN(`Emitter Cell Type`),
  by = c("Dataset", "Tissue")
]

dataset_summary[
  md_age_dc,
  on = c("Dataset", "Tissue"),
  age_group2 := i.age_group2
]


ggplot2::ggplot(
  dataset_summary,
  ggplot2::aes(
    Dataset,
    Tissue
  )
) + ggplot2::geom_tile(
  ggplot2::aes(
    width = 0.9,
    height = 0.9
  ),
  colour = "black",
  fill = "coral",
  alpha = 1
) + ggplot2::ggtitle(
  "Number of cell types per dataset (additional age information)"
) + ggplot2::geom_text(
  #ggplot2::aes(label = paste(V1, "Cell Types", age_group2)),
  ggplot2::aes(label = paste(V1, age_group2)),
  size = 8,
  fontface = "bold"
) + ggplot2::scale_x_discrete(
  limits = c(
    "TMS FACS (male)",
    "TMS FACS (female)" ,
    "TMS Droplet (male)",
    "TMS Droplet (female)",
    "Calico Droplet (male)"
  ),
  labels = c(
    "TMS\nFACS\n(male)\n3m vs 18/24m",
    "TMS\nFACS\n(female)\n3m vs 18m",
    "TMS\nDroplet\n(male)\n3m vs 18/24m",
    "TMS\nDroplet\n(female)\n3m vs 18/21m",
    "Calico\nDroplet\n(male)\n7-8m vs 22-23m"
  )
) + ggplot2::scale_y_discrete(
  limits = sort(
    unique(dataset_summary$Tissue),
    decreasing = TRUE
  )
) + ggplot2::xlab(
  ""
) + ggplot2::ylab(
  ""
) + ggplot2::theme(
  text = ggplot2::element_text(size = 28),
  axis.text = ggplot2::element_text(size = 28, face = "bold")
)
#manual save 2000x1400

## Prepare supplemental data 2 ####

supp2_table <- lapply(
  dataset_paths,
  function(path) {
    tissues <- gsub(".*scdiffcom_(.+)\\.rds.*", "\\1", list.files(path))
    dataset <- lapply(
      X = tissues,
      FUN = function(
        tiss
      ) {
        res <- readRDS(paste0(path, "/scdiffcom_", tiss, ".rds"))
      }
    )
    tissues[tissues == "BAT"] <- "Adipose_Brown"
    tissues[tissues == "GAT"] <- "Adipose_Gonadal"
    tissues[tissues == "MAT"] <- "Adipose_Mesenteric"
    tissues[tissues == "SCAT"] <- "Adipose_Subcutaneous"
    dataset <- lapply(
      seq_along(dataset), 
      function (i) {
        dataset[[i]]@parameters$object_name <- tissues[[i]]
        dataset[[i]]
      }
    )
    names(dataset) <- tissues
    regulation_dt <- rbindlist(
      lapply(
        dataset,
        function(i) {
          cci_dt <- copy(i@cci_table_raw)
          dt <- cci_dt[
            IS_CCI_EXPRESSED_YOUNG == TRUE | IS_CCI_EXPRESSED_OLD == TRUE
          ]
          dt[
            ,
            c(
              "BH_P_VALUE_YOUNG",
              "BH_P_VALUE_OLD",
              "BH_P_VALUE_DE"
            ) := list(
              stats::p.adjust(P_VALUE_YOUNG, method = "BH"),
              stats::p.adjust(P_VALUE_OLD, method = "BH"),
              stats::p.adjust(P_VALUE_DE, method = "BH")
            )
          ]
          dt[
            ,
            c(
              "IS_CCI_SCORE_YOUNG",
              "IS_CCI_SCORE_OLD",
              "IS_CCI_SPECIFIC_YOUNG",
              "IS_CCI_SPECIFIC_OLD",
              "IS_DE_LOGFC",
              "IS_DE_SIGNIFICANT",
              "DE_DIRECTION",
              "IS_CCI_DETECTED_YOUNG",
              "IS_CCI_DETECTED_OLD",
              "IS_CCI_DE"
            ) := {
              threshold_score_temp <- stats::quantile(
                x = c(
                  .SD[IS_CCI_EXPRESSED_YOUNG == TRUE][["CCI_SCORE_YOUNG"]],
                  .SD[IS_CCI_EXPRESSED_OLD == TRUE][["CCI_SCORE_OLD"]]
                ),
                probs = 0.2
              )
              is_cci_score_1 <- (CCI_SCORE_YOUNG >= threshold_score_temp)
              is_cci_score_2 <- (CCI_SCORE_OLD >= threshold_score_temp)
              is_cci_specific_1 <- BH_P_VALUE_YOUNG <= 0.05
              is_cci_specific_2 <- BH_P_VALUE_OLD <= 0.05
              is_de_logfc <- LOGFC_ABS >= log(1.5)
              is_de_significant <- BH_P_VALUE_DE <= 0.05
              de_direction <- fifelse(LOGFC > 0, "UP", "DOWN")
              is_cci_detected_1 <- (IS_CCI_EXPRESSED_YOUNG == TRUE) &
                is_cci_score_1 & is_cci_specific_1
              is_cci_detected_2 <- (IS_CCI_EXPRESSED_OLD == TRUE) &
                is_cci_score_2 & is_cci_specific_2
              is_cci_de <- is_de_logfc & is_de_significant
              list(
                is_cci_score_1, is_cci_score_2,
                is_cci_specific_1, is_cci_specific_2,
                is_de_logfc, is_de_significant, de_direction,
                is_cci_detected_1, is_cci_detected_2,
                is_cci_de
              )
            }
          ]
          dt[
            ,
            REGULATION :=
              ifelse(
                !IS_CCI_DETECTED_YOUNG & !IS_CCI_DETECTED_OLD,
                "NOT_DETECTED",
                ifelse(
                  IS_DE_LOGFC & IS_DE_SIGNIFICANT & DE_DIRECTION == "UP",
                  "UP",
                  ifelse(
                    IS_DE_LOGFC & IS_DE_SIGNIFICANT & DE_DIRECTION == "DOWN",
                    "DOWN",
                    ifelse(
                      !IS_DE_LOGFC,
                      "FLAT",
                      ifelse(
                        IS_DE_LOGFC & !IS_DE_SIGNIFICANT,
                        "NSC",
                        "There is a problem here!"
                      )
                    )
                  )
                )
              )
          ]
          dt <- dt[
            ,
            c(
              "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LRI",
              "IS_CCI_EXPRESSED_YOUNG",
              "IS_CCI_EXPRESSED_OLD",
              "IS_CCI_SCORE_YOUNG",
              "IS_CCI_SCORE_OLD",
              "IS_CCI_SPECIFIC_YOUNG",
              "IS_CCI_SPECIFIC_OLD",
              "IS_DE_LOGFC",
              "IS_DE_SIGNIFICANT",
              "DE_DIRECTION",
              "IS_CCI_DETECTED_YOUNG",
              "IS_CCI_DETECTED_OLD",
              "IS_CCI_DE",
              "REGULATION"
            )
          ]
          cci_dt <- cci_dt[
            ,
            c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LRI")
          ]
          cci_dt <- merge.data.table(
            cci_dt,
            dt,
            by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LRI"),
            all.x = TRUE,
            sort = FALSE
          )
          cci_dt[
            ,
            ER_CELLTYPES := paste(
              EMITTER_CELLTYPE,
              RECEIVER_CELLTYPE,
              sep = "_")]
          cci_dt[, CCI := paste(ER_CELLTYPES, LRI, sep = "_")]
          cci_dt[, DE_DIRECTION := ifelse(is.na(DE_DIRECTION), "NONE", DE_DIRECTION)]
          cci_dt[, REGULATION := ifelse(is.na(REGULATION), "NOT_DETECTED", REGULATION)]
          cci_dt[is.na(cci_dt)] <- FALSE
          cci_dt[
            ,
            .N,
            by =  c(
              "IS_CCI_EXPRESSED_YOUNG",
              "IS_CCI_SCORE_YOUNG",
              "IS_CCI_SPECIFIC_YOUNG",
              "IS_CCI_EXPRESSED_OLD",
              "IS_CCI_SCORE_OLD",
              "IS_CCI_SPECIFIC_OLD",
              "IS_DE_LOGFC",
              "IS_DE_SIGNIFICANT",
              "DE_DIRECTION",
              "IS_CCI_DETECTED_YOUNG",
              "IS_CCI_DETECTED_OLD",
              "IS_CCI_DE",
              "REGULATION"
            )
          ]
        }
      ),
      idcol = "Tissue"
    )
  }
)
names(supp2_table) <- dataset_names
supp2_table <- rbindlist(
  supp2_table,
  idcol = "Dataset"
)

supp2_table_non_detected <- supp2_table[REGULATION == "NOT_DETECTED"]
supp2_table_detected <- supp2_table[REGULATION != "NOT_DETECTED"]
supp2_table_detected[
  ,
  PCT := N/sum(N)*100,
  by = c("Dataset", "Tissue")
]

fwrite(
  supp2_table_detected,
  "../data_scAgeCom/analysis/outputs_data/data_4_supp2_table_detected.csv"
)


# supp2_table_final <- dcast.data.table(
#   supp2_table_clean[, -c("N")],
#   ... ~ Dataset + Tissue,
#   value.var = "PCT"
# )
# supp2_table_final[is.na(supp2_table_final)] <- 0
# supp2_table_final[
#   ,
#   AVG_PCT := rowMeans(.SD),
#   .SDcols = 14:71
# ]
# supp2_table_final <- supp2_table_final[,c(1:13,72)]
# 
# supp2_template <- CJ(
#   IS_CCI_EXPRESSED_YOUNG = c(TRUE, FALSE),
#   IS_CCI_SCORE_YOUNG = c(TRUE, FALSE),
#   IS_CCI_SPECIFIC_YOUNG = c(TRUE, FALSE),
#   IS_CCI_EXPRESSED_OLD = c(TRUE, FALSE),
#   IS_CCI_SCORE_OLD = c(TRUE, FALSE),
#   IS_CCI_SPECIFIC_OLD = c(TRUE, FALSE),
#   IS_DE_LOGFC = c(TRUE, FALSE),
#   IS_DE_SIGNIFICANT = c(TRUE, FALSE),
#   DE_DIRECTION = c("NONE", "UP", "DOWN"),
#   IS_CCI_DETECTED_YOUNG = c(TRUE, FALSE),
#   IS_CCI_DETECTED_OLD = c(TRUE, FALSE),
#   IS_CCI_DE = c(TRUE, FALSE)
# )
