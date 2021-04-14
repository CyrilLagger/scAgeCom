####################################################
##
## Project: scAgeCom
##
## Last update - April 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## preparation of tissue specific results
##
####################################################
##

## Libraries ####

library(scDiffCom)
library(data.table)
library(ggplot2)
library(htmltools)
library(plotly)

## relative path of scDiffCom results ####

scAgeCom_path <- "../data_scAgeCom/scDiffCom_results_12_04_2021"

# it contains 5 datasets

dataset_paths <- list.dirs(scAgeCom_path, recursive = FALSE)
dataset_paths

dataset_names <- c(
  "Calico2019",
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
        categories = NULL,
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
    max_L <- max(.SD[is.finite(LOG2FC_L)][["LOG2FC_L"]])
    min_L <- min(.SD[is.finite(LOG2FC_L)][["LOG2FC_L"]])
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
  "RECEPTOR_3"
)
CCI_cols_informative_shiny <- c(
  "Dataset",
  "Tissue",
  "Ligand-Receptor Interaction",
  "Emitter Cell Type",
  "Receiver Cell Type",
  "Log2 FC",
  "Adj. P-value",
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
  "RECEPTOR_3"
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
  `Adj. P-value`,
  -`Old CCI Score`,
  -`Young CCI Score`
)

## Create a table of counts summary for shiny #####

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
    "Down CCI", "Flat CCI", "NSC CCI", "Up CCI",
    "Total CCI", "Total Cell Types"
  )
)

setcolorder(
  TISSUE_COUNTS_SUMMARY,
  c(
    "Tissue", "Dataset", "Total Cell Types",
    "Total CCI", "Flat CCI", "Down CCI", "Up CCI", 
    "NSC CCI"
  )
)

## create utility function to diplay tissue counts summary on shiny #####

build_tissue_counts_display <- function(
  tissue_counts_summary
) {
  dt <- tissue_counts_summary[, -c(8)]
  DT::datatable(
    data = dt,
    options =list(
      pageLength = 10,
      dom = "t"
    ),
    rownames = FALSE,
    callback = htmlwidgets::JS(
      "var tips = [
        'Selected Tissue',
        'Dataset on which the interaction analysis has been performed',
        'Number of cell-types in the tissue of interest',
        'Number of cell-cell interactions detected in the tissue of interest',
        'Number of stable CCIs with age',
        'Number of down-regulated CCIs with age',
        'Number of up-regulated CCIs with age'
        ],
        header = table.columns().header();
        for (var i = 0; i < tips.length; i++) {
        $(header[i]).attr('title', tips[i]);
        }"
    )
  )
}


## create utility functions for CCI table on shiny ####

subset_CCI_table <- function(
  CCI_table,
  dataset_choice,
  tissue_choice,
  emitter_choice = NULL,
  receiver_choice = NULL,
  LRI_choice = NULL,
  GENE_choice = NULL,
  filter
) {
  if (filter) {
    if ('All LRIs' %in% LRI_choice) {
      if ("All Genes" %in% GENE_choice) {
        dt <- CCI_table[
          Dataset == dataset_choice &
            Tissue == tissue_choice &
            `Emitter Cell Type` %in% emitter_choice &
            `Receiver Cell Type` %in% receiver_choice
        ]
      } else {
        dt <- CCI_table[
          Dataset == dataset_choice &
            Tissue == tissue_choice &
            `Emitter Cell Type` %in% emitter_choice &
            `Receiver Cell Type` %in% receiver_choice &
            (LIGAND_1 %in% GENE_choice | LIGAND_2 %in% GENE_choice |
               RECEPTOR_1 %in% GENE_choice | RECEPTOR_2 %in% GENE_choice |
               RECEPTOR_3 %in% GENE_choice)
        ]
      }
    } else {
      if ("All Genes" %in% GENE_choice) {
        dt <- CCI_table[
          Dataset == dataset_choice &
            Tissue == tissue_choice &
            `Emitter Cell Type` %in% emitter_choice &
            `Receiver Cell Type` %in% receiver_choice &
            `Ligand-Receptor Interaction` %in% LRI_choice
        ]
      } else {
        dt <- CCI_table[
          Dataset == dataset_choice &
            Tissue == tissue_choice &
            `Emitter Cell Type` %in% emitter_choice &
            `Receiver Cell Type` %in% receiver_choice &
            `Ligand-Receptor Interaction` %in% LRI_choice &
            (LIGAND_1 %in% GENE_choice | LIGAND_2 %in% GENE_choice |
               RECEPTOR_1 %in% GENE_choice | RECEPTOR_2 %in% GENE_choice |
               RECEPTOR_3 %in% GENE_choice)
        ]
      }
    }
  } else {
    dt <- CCI_table[
      Dataset == dataset_choice &
        Tissue == tissue_choice
    ]
  }
  dt
}

build_CCI_display <- function(
  CCI_table
) {
  dt <- CCI_table[
    ,
    -c(1,2,13,14,15,16,17,18,19,20,21,22,23)
  ]
  CCI_DT <- DT::datatable(
    data = dt[, -c(9, 10)],
    options =list(
      pageLength = 10
    ),
    caption = tags$caption(
      style = paste0(
        'caption-side: top; text-align: center; ',
        'color:black; font-size:150% ;'
      ),
      "Table of Cell-Cell Interactions"
    ),
    rownames = rownames,
    extensions = c("Buttons")
  )
  vline <- function(x = 0, color = "black") {
    list(
      type = "line", 
      y0 = 0, 
      y1 = 1, 
      yref = "paper",
      x0 = x, 
      x1 = x, 
      line = list(color = color)
    )
  }
  hline <- function(y = 0, color = "black") {
    list(
      type = "line", 
      x0 = 0, 
      x1 = 1, 
      xref = "paper",
      y0 = y, 
      y1 = y, 
      line = list(color = color)
    )
  }
  dt$`Age Regulation` <- factor(
    dt$`Age Regulation`,
    levels = c("UP", "DOWN", "FLAT", "NSC")
  )
  m <- list(
    l = 100,
    r = 100,
    b = 100,
    t = 100,
    pad = 8
  )
  CCI_VOLCANO_PLOT <- plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~`Log2 FC`,
    y = ~-log10(`Adj. P-value`),
    text = ~paste(
      "LRI: ",
      `Ligand-Receptor Interaction`, 
      '$<br>Emitter:',
      `Emitter Cell Type`,
      '$<br>Receiver:',
      `Receiver Cell Type`
    ),
    color = ~`Age Regulation`,
    colors = setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    )#,
    #width = 900,
    #height = 500
  ) %>% plotly::layout(
    title = "Interactive Aging Volcano Plot",
    font = list(size = 20),
    xaxis = list(
      title = "Log2(FC)",
      titlefont = list(size = 18)
    ),
    yaxis = list(
      title = "-Log10(Adj. p-value)",
      titlefont = list(size = 18)
    ),
    shapes = list(
      vline(log2(1.5)),
      vline(-log2(1.5)),
      hline(-log10(0.05))
    ),
    legend = list(
      x = 100,
      y = 0.5
    ),
    margin = m
  )  %>% toWebGL()
  CCI_SCORE_PLOT <- plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~log10(`Young CCI Score`),
    y = ~log10(`Old CCI Score`),
    text = ~paste(
      "LRI: ",
      `Ligand-Receptor Interaction`, 
      '$<br>Emitter:',
      `Emitter Cell Type`,
      '$<br>Receiver:',
      `Receiver Cell Type`
    ),
    color = ~`Age Regulation`,
    colors = setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    )
  ) %>% plotly::layout(
    title = "Interactive Score Plot",
    font = list(size = 20),
    xaxis = list(
      title = "Log10(Young CCI Score)",
      titlefont = list(size = 18)
    ),
    yaxis = list(
      title = "Log10(Old CCI Score)",
      titlefont = list(size = 18)
    ),
    legend = list(
      x = 100,
      y = 0.5
    ),
    margin = m
  ) %>% toWebGL()
  CCI_LRFC_PLOT <- plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~`Ligand Log2 FC`,
    y = ~`Receptor Log2 FC`,
    text = ~paste(
      "LRI: ",
      `Ligand-Receptor Interaction`, 
      '$<br>Emitter:',
      `Emitter Cell Type`,
      '$<br>Receiver:',
      `Receiver Cell Type`
    ),
    color = ~`Age Regulation`,
    colors = setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    )
  )  %>% plotly::layout(
    title = "Interactive 'Ligand-FC vs Receptor-FC' Plot",
    font = list(size = 20),
    xaxis = list(
      title = "Ligand Log2(FC)",
      titlefont = list(size = 18)
    ),
    yaxis = list(
      title = "Receptor Log2(FC)",
      titlefont = list(size = 18)
    ),
    legend = list(
      x = 100,
      y = 0.5
    ),
    margin = m
  )  %>% toWebGL()
  CCI_PLOTS <- subplot(
    CCI_VOLCANO_PLOT,
    CCI_SCORE_PLOT,
    CCI_LRFC_PLOT,
    nrows = 3,
    titleX = TRUE,
    titleY = TRUE
  )
  list(
    CCI_DT = CCI_DT,
    CCI_VOLCANO_PLOT = CCI_VOLCANO_PLOT,
    CCI_SCORE_PLOT = CCI_SCORE_PLOT,
    CCI_LRFC_PLOT = CCI_LRFC_PLOT
  )
}

# build_CCI_display(
#   CCI_table = subset_CCI_table(
#     CCI_table,
#     dataset_choice = dataset_names[[5]],
#     tissue_choice = "Liver",
#     filter = FALSE
#   )
# )$CCI_DT



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

## create utility functions for CCI table on shiny ####

subset_ORA_table <- function(
  ORA_table,
  dataset_choice,
  tissue_choice
) {
  dt <- ORA_table[
    Dataset == dataset_choice &
      Tissue == tissue_choice
  ]
  category_keep <- c(
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
  )
  dt <- dt[ORA_CATEGORY %in% category_keep]
  dt[data.table(
    category_old = category_keep,
    category_new = c(
      "LRIs",
      "Ligand Gene(s)",
      "Receptor Gene(s)",
      "ER Cell Types",
      "Emitter Cell Types",
      "Receiver Cell Types",
      "GO Terms",
      "KEGG Pathways",
      "ER Cell Families",
      "Emitter Cell Families",
      "Receiver Cell Families"
    )
  ),
  on = c("ORA_CATEGORY==category_old"),
  ORA_CATEGORY := i.category_new
  ]
  dt[
    data.table(
      old_aspect = c(
        "biological_process",
        "cellular_component",
        "molecular_function"
      ),
      new_aspect = c(
        "Biological Process",
        "Cellular Component",
        "Molecular Function"
      )
    ),
    on = "ASPECT==old_aspect",
    ASPECT := i.new_aspect
  ]
  dt
}

build_ORA_display <- function(
  ORA_table,
  category_choice,
  go_aspect_choice,
  type_choice
) {
  dt <- ORA_table[ORA_CATEGORY == category_choice]
  if (category_choice == "GO Terms") {
    dt <- dt[ASPECT == go_aspect_choice]
    level_str <- "GO Level"
    options = list(
      #columnDefs = list(list(targets = c(1,2,), searchable = FALSE)),
      pageLength = 10
    )
    filter <- "top"
  } else {
    level_str <- NULL
    options <-list(pageLength = 10)
    filter <- "none"
  }
  if(type_choice == "UP") {
    cols_to_keep <- c(
      "VALUE",
      "ORA_SCORE_UP",
      "OR_UP",
      "BH_P_VALUE_UP",
      level_str
    )
    dt <- dt[
      OR_UP >= 1 & BH_P_VALUE_UP <= 0.05,
      cols_to_keep,
      with = FALSE
    ]
  } else if(type_choice == "DOWN") {
    cols_to_keep <- c(
      "VALUE",
      "ORA_SCORE_DOWN",
      "OR_DOWN",
      "BH_P_VALUE_DOWN",
      level_str
    )
    dt <- dt[
      `OR_DOWN` >= 1 & BH_P_VALUE_DOWN <= 0.05,
      cols_to_keep,
      with = FALSE
    ]
  } else if(type_choice == "FLAT") {
    cols_to_keep <- c(
      "VALUE",
      "ORA_SCORE_FLAT",
      "OR_FLAT",
      "BH_P_VALUE_FLAT",
      level_str
    )
    dt <- dt[
      `OR_FLAT` >= 1 & BH_P_VALUE_FLAT <= 0.05,
      cols_to_keep,
      with = FALSE]
  }
  if (category_choice == "GO Terms") {
    dt[, `GO Level` := as.factor(`GO Level`)]
  }
  setnames(
    dt,
    old = colnames(dt),
    new = c(
      category_choice,
      "ORA Score",
      "Odds Ratio",
      "Adj. p-value",
      level_str
    )
  )
  setorder(dt, -`ORA Score`)
  DT::datatable(
    data = dt,
    options = options,
    caption = tags$caption(
      style = 'caption-side: top; text-align: center; color:black; font-size:150% ;',
      paste0(
        category_choice,
        " over-represented among ",
        type_choice,
        "-regulated CCIs"
      )
    ),
    rownames = FALSE,
    filter = filter
  )
}

build_ORA_visnetwork <- function(
  CCI_table,
  ORA_table,
  tissue_choice,
  dataset_choice,
  abbr_celltype
) {
  CCI_dt <- copy(CCI_table)
  setnames(
    CCI_dt,
    old = c("Emitter Cell Type", "Receiver Cell Type", "Age Regulation"),
    new = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION")
  )
  ora_table_ER <- ORA_table[
    Dataset == dataset_choice &
      Tissue == tissue_choice &
      ORA_CATEGORY == "ER_CELLTYPES"
  ]
  cci_table_detected <- CCI_dt[
    Dataset == dataset_choice &
      Tissue == tissue_choice
  ]
  actual_celltypes <- union(
    cci_table_detected[["EMITTER_CELLTYPE"]],
    cci_table_detected[["RECEIVER_CELLTYPE"]]
  )
  abbreviation_table <- abbr_celltype[[dataset_choice]][
    ORIGINAL_CELLTYPE %in% actual_celltypes
  ]
  if (!identical(
    sort(actual_celltypes),
    sort(abbreviation_table[["ORIGINAL_CELLTYPE"]])
  )) {
    stop(
      paste0(
        "No abbreviation will be used:",
        " `abbreviation table` must contain",
        " a column with the original cell-types")
    )
  } else if (sum(duplicated(abbreviation_table)) > 0) {
    stop(
      paste0(
        "No abbreviation will be used:",
        " `abbreviation table` must not contain duplicated rows"))
  } else {
    cci_table_detected[
      abbreviation_table,
      on = "EMITTER_CELLTYPE==ORIGINAL_CELLTYPE",
      "EMITTER_CELLTYPE" := i.ABBR_CELLTYPE]
    cci_table_detected[
      abbreviation_table,
      on = "RECEIVER_CELLTYPE==ORIGINAL_CELLTYPE",
      "RECEIVER_CELLTYPE" := i.ABBR_CELLTYPE]
    ora_table_ER[
      abbreviation_table,
      on = "EMITTER_CELLTYPE==ORIGINAL_CELLTYPE",
      "EMITTER_CELLTYPE" := i.ABBR_CELLTYPE]
    ora_table_ER[
      abbreviation_table,
      on = "RECEIVER_CELLTYPE==ORIGINAL_CELLTYPE",
      "RECEIVER_CELLTYPE" := i.ABBR_CELLTYPE]
  }
  scDiffCom:::interactive_from_igraph(
    cci_table_detected = cci_table_detected,
    conds = c("YOUNG", "OLD"),
    ora_table_ER = ora_table_ER,
    ora_table_LR = ORA_table[
      Dataset == dataset_choice &
        Tissue == tissue_choice &
        ORA_CATEGORY == "LRI"
    ],
    network_type = "ORA_network",
    layout_type = "bipartite",
    object_name = tissue_choice
  )
}

build_ORA_plot <- function(
  ORA_table,
  tissue_choice,
  dataset_choice,
  category_choice,
  type_choice,
  go_aspect_choice = "biological_process"
) {
  dt <- ORA_table[
    Dataset == dataset_choice &
      Tissue == tissue_choice &
      ORA_CATEGORY == category_choice
  ]
  setnames(
    dt,
    "GO Level",
    "LEVEL"
  )
  scDiffCom:::plot_ora(
    ora_dt = dt,
    category = category_choice,
    regulation = type_choice,
    max_terms_show = 20,
    GO_aspect = go_aspect_choice,
    OR_threshold = 1,
    bh_p_value_threshold = 0.05
  )
}

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

ALL_ORA_CATEGORIES <- c(
  "LRIs",
  "Ligand Gene(s)",
  "Receptor Gene(s)",
  "ER Cell Types",
  "Emitter Cell Types",
  "Receiver Cell Types",
  "GO Terms",
  "KEGG Pathways",
  "ER Cell Families",
  "Emitter Cell Families",
  "Receiver Cell Families"
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

## save all results ####

data_4_tissue_specific_results <- list(
  #dataset_processed = dataset_processed,
  CCI_table = CCI_table,
  ORA_table = ORA_table,
  TISSUE_COUNTS_SUMMARY = TISSUE_COUNTS_SUMMARY,
  ALL_TISSUES = ALL_TISSUES,
  ALL_CELLTYPES = ALL_CELLTYPES,
  ALL_LRIs = ALL_LRIs,
  ALL_GENES = ALL_GENES,
  ALL_ORA_CATEGORIES = ALL_ORA_CATEGORIES,
  ALL_ORA_CATEGORIES_SPECIFIC = ALL_ORA_CATEGORIES_SPECIFIC,
  ALL_ORA_GO_ASPECTS = ALL_ORA_GO_ASPECTS,
  ALL_ORA_TYPES = ALL_ORA_TYPES,
  ABBR_CELLTYPE = ABBR_CELLTYPE,
  build_CCI_display = build_CCI_display,
  build_ORA_display = build_ORA_display,
  build_ORA_plot = build_ORA_plot,
  build_ORA_visnetwork = build_ORA_visnetwork,
  build_tissue_counts_display = build_tissue_counts_display,
  subset_CCI_table = subset_CCI_table,
  subset_ORA_table = subset_ORA_table
)

saveRDS(
  data_4_tissue_specific_results,
  "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds"
)
