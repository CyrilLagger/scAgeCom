
## Load libraries ####

library(scDiffCom)
library(htmltools)
library(data.table)
library(latex2exp)
library(ggplot2)
library(plotly)

## Loading and cleaning mouse LRI ####

LRI_mouse_curated <- copy(scDiffCom::LRI_mouse$LRI_curated)
LRI_mouse_curated[, SOURCE := gsub("SCT:SingleCellSignalR|SCT:CellPhoneDB|SCT:DLRP", "", SOURCE)]
LRI_mouse_curated[, SOURCE := gsub(";;", ";", SOURCE)]
LRI_mouse_curated[, SOURCE := sub(";$", "", SOURCE)]
LRI_mouse_curated[, SOURCE := sub("^;", "", SOURCE)]
LRI_cols_to_keep <- c(
  "LIGAND_1", "LIGAND_2",
  "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
  "DATABASE", "SOURCE"
)
LRI_mouse_curated <- LRI_mouse_curated[, LRI_cols_to_keep, with = FALSE]
setnames(
  LRI_mouse_curated,
  old = LRI_cols_to_keep,
  new = c(
    "Ligand (1)", "Ligand (2)",
    "Receptor (1)", "Receptor (2)", "Receptor (3)",
    "Database of Origin", "Retrieved Sources"
  )
)

LRI_DATABASES <- sort(
  unique(
    unlist(
      strsplit(
        LRI_mouse_curated$`Database of Origin`,
        ";")
    )
  )
)

LRI_mouse_curated[, COMPLEX := !is.na(`Ligand (2)`) | !is.na(`Receptor (2)`)]
LRI_mouse_curated[, c(LRI_DATABASES) := lapply(LRI_DATABASES, function(i) {
  ifelse(grepl(i, `Database of Origin`), TRUE, FALSE)
})]

## Loading and cleaning human LRI ####

LRI_human_curated <- copy(scDiffCom::LRI_human$LRI_curated)
LRI_human_curated[, SOURCE := gsub("SCT:SingleCellSignalR|SCT:CellPhoneDB|SCT:DLRP", "", SOURCE)]
LRI_human_curated[, SOURCE := gsub(";;", ";", SOURCE)]
LRI_human_curated[, SOURCE := sub(";$", "", SOURCE)]
LRI_human_curated[, SOURCE := sub("^;", "", SOURCE)]
LRI_cols_to_keep <- c(
  "LIGAND_1", "LIGAND_2",
  "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
  "DATABASE", "SOURCE"
)
LRI_human_curated <- LRI_human_curated[, LRI_cols_to_keep, with = FALSE]
setnames(
  LRI_human_curated,
  old = LRI_cols_to_keep,
  new = c(
    "Ligand (1)", "Ligand (2)",
    "Receptor (1)", "Receptor (2)", "Receptor (3)",
    "Database of Origin", "Retrieved Sources"
  )
)

LRI_DATABASES <- sort(
  unique(
    unlist(
      strsplit(
        LRI_human_curated$`Database of Origin`,
        ";")
    )
  )
)

LRI_human_curated[, COMPLEX := !is.na(`Ligand (2)`) | !is.na(`Receptor (2)`)]
LRI_human_curated[, c(LRI_DATABASES) := lapply(LRI_DATABASES, function(i) {
  ifelse(grepl(i, `Database of Origin`), TRUE, FALSE)
})]

## LRI functions ####

plot_lri_upset <- function(
  LRI_table,
  groups
) {
  p <- ComplexUpset::upset(
    as.data.frame(LRI_table),
    groups,
    base_annotations = list(
      'Intersection size' = intersection_size(
        mapping = aes(fill = COMPLEX),
        counts = TRUE,
        bar_number_threshold = 100
      )
    ),
    themes = upset_default_themes(text = element_text(size = 20)),
    min_size = 35
  ) +
    ggtitle("LRI Overlap by Databases of Origin")
  return(p)
}

build_LRI_display <- function(
  LRI_table,
  LRI_database
) {
  dt <- LRI_table[
    apply(
      sapply(
        LRI_database,
        function(i) {
          grepl(i, `Database of Origin`)
        }
      ),
      MARGIN = 1,
      any
    )
  ]
  DT::datatable(
    data = dt[, 1:7],
    options = list(
      pageLength = 10,
      columnDefs = list(
        list(
          targets = c(6,7),
          render = htmlwidgets::JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 20 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
            "}")
        )
      )
    ),
    caption = tags$caption(
      style = 'caption-side: top; text-align: center; color:black; font-size:150% ;',
      "Table of Ligand-Receptor Interactions"
    )
  )
}

## Loading scDiffCom objects ####

DATASETS_COMBINED <- readRDS(
  "../data_scAgeCom/analysis/outputs_data/scagecom_results_processed.rds"
)

shiny_dataset_names <- c(
  "Calico2019",
  "TMS Droplet (female)",
  "TMS Droplet (male)",
  "TMS FACS (female)",
  "TMS FACS (male)"
)

names(DATASETS_COMBINED)
names(DATASETS_COMBINED) <- shiny_dataset_names
names(DATASETS_COMBINED)

## Create a single CCI table  with only informative columns ####

CCI_table <- rbindlist(
  lapply(
    DATASETS_COMBINED,
    GetTableCCI,
    type = "detected",
    simplified = FALSE
  ),
  idcol = "DATASET"
)

# Add and round LOG2FC
CCI_table[, LOG2FC_BASE := LOGFC*log2(exp(1))]
CCI_table[, LOG2FC_ID := {
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
  )}, by = c("ID", "DATASET")]
CCI_table[
  ,
  LOG2FC_ALL := {
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
    )}
]

# Add ligand (resp. receptor) logfc

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
  by = c("ID", "DATASET")
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
  by = c("ID", "DATASET")
]

# Round numeric 
CCI_table[, CCI_SCORE_YOUNG := signif(CCI_SCORE_YOUNG, 4)]
CCI_table[, CCI_SCORE_OLD := signif(CCI_SCORE_OLD, 4)]
CCI_table[, LOG2FC_ID := signif(LOG2FC_ID, 3)]
CCI_table[, BH_P_VALUE_DE := signif(BH_P_VALUE_DE, 3)]
CCI_table[, LOG2FC_L := signif(LOG2FC_L, 3)]
CCI_table[, LOG2FC_R := signif(LOG2FC_R, 3)]

# only keep relevant columns 

colnames(CCI_table)

CCI_cols_informative <- c(
  "DATASET",
  "ID",
  "LRI",
  "EMITTER_CELLTYPE",
  "RECEIVER_CELLTYPE",
  "LOG2FC_ID",
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
  "IS_CCI_EXPRESSED_OLD"
)

CCI_cols_informative_user <- c(
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
  "IS_CCI_EXPRESSED_OLD"
)

CCI_table <- CCI_table[, CCI_cols_informative, with = FALSE]
setnames(
  CCI_table,
  old = CCI_cols_informative,
  new = CCI_cols_informative_user
)
setkey(CCI_table)

setorder(
  CCI_table,
  -`Log2 FC`,
  `Adj. P-value`,
  -`Old CCI Score`,
  -`Young CCI Score`
)

## Display CCI Details for Shiny ####

subset_cci_table <- function(
  CCI_table,
  dataset_choice,
  tissue_choice,
  emitter_choice = NULL,
  receiver_choice = NULL,
  LRI_choice = NULL,
  filter
) {
  if (filter) {
    if ('All LRIs' %in% LRI_choice) {
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
          `Ligand-Receptor Interaction` %in% LRI_choice
      ]
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
    -c(1,2,13,14,15,16,17,18)
  ]
  CCI_DT <- DT::datatable(
    data = dt[, -c(9, 10)],
    options =list(
      pageLength = 10
    ),
    caption = tags$caption(
      style = 'caption-side: top; text-align: center; color:black; font-size:150% ;',
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
    title = "Interactive Volcano Plot",
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
    title = "Interactive LRFC Plot",
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
#   CCI_table = subset_cci_table(
#     CCI_table,
#     dataset_choice = shiny_dataset_names[[5]],
#     tissue_choice = "Liver",
#     filter = FALSE
#   )
# )$CCI_LRFC_PLOT

# microbenchmark::microbenchmark(
#   a = build_CCI_display(
#     CCI_table = subset_cci_table(
#       CCI_table,
#       dataset_choice = shiny_dataset_names[[5]],
#       tissue_choice = "Lung",
#       filter = FALSE
#     )
#   ),
#   times = 10
# )

## Create a single ORA table with only informative columns ####

ORA_table <- lapply(
  DATASETS_COMBINED,
  GetTableORA,
  categories = "all",
  simplified = FALSE
)

ORA_table <- rbindlist(
  lapply(
    ORA_table,
    function(dataset) {
      rbindlist(
        dataset,
        fill = TRUE,
        idcol = "ORA_CATEGORY"
      )
    }
  ),
  idcol = "DATASET"
)

# Add cell types for ERI

ORA_table[
  ORA_CATEGORY == "ER_CELLTYPES",
  c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE") := list(
  sub("_.*", "", VALUE),
  sub(".*_", "", VALUE)
)]

# only keep relevant columns 

colnames(ORA_table)

ORA_cols_informative <- c(
  "DATASET",
  "ORA_CATEGORY",
  "ID",
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

ORA_cols_informative_user <- c(
  "Dataset",
  "ORA_CATEGORY",
  "Tissue",
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
  new = ORA_cols_informative_user
)
setkey(ORA_table)

# Round numeric 
ORA_table[, ORA_SCORE_UP := signif(ORA_SCORE_UP, 3)]
ORA_table[, ORA_SCORE_DOWN := signif(ORA_SCORE_DOWN, 3)]
ORA_table[, ORA_SCORE_FLAT := signif(ORA_SCORE_FLAT, 3)]
ORA_table[, OR_UP := signif(OR_UP, 3)]
ORA_table[, OR_DOWN := signif(OR_DOWN, 3)]
ORA_table[, OR_FLAT := signif(OR_FLAT, 3)]
ORA_table[, BH_P_VALUE_UP := signif(BH_P_VALUE_UP, 3)]
ORA_table[, BH_P_VALUE_DOWN := signif(BH_P_VALUE_DOWN, 3)]
ORA_table[, BH_P_VALUE_FLAT := signif(BH_P_VALUE_FLAT, 3)]


##

ALL_ORA_CATEGORIES <- c(
  "LRIs",
  "Cell Types",
  "GO Terms",
  "KEGG Pathways",
  "Cell Families"
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
    "KEGG_PWS",
    "GO_TERMS",
    "LRI",
    "ER_CELLTYPES",
    "ER_CELLFAMILIES"
    )
  dt <- dt[ORA_CATEGORY %in% category_keep]
  dt[data.table(
    category_old = category_keep,
    category_new = c(
      "KEGG Pathways",
      "GO Terms",
      "LRIs",
      "Cell Types",
      "Cell Families"
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
  dataset_choice
) {
  CCI_dt <- copy(CCI_table)
  setnames(
    CCI_dt,
    old = c("Emitter Cell Type", "Receiver Cell Type"),
    new = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
  )
  ora_table_ER <- ORA_table[
    Dataset == dataset_choice &
      Tissue == tissue_choice &
      ORA_CATEGORY == "ER_CELLTYPES"
  ]
  scDiffCom:::interactive_from_igraph(
    cci_table_detected = CCI_dt[
      Dataset == dataset_choice &
        Tissue == tissue_choice
    ],
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

## Loading and cleaning tissue cci counts ####

TISSUE_COUNTS_SUMMARY <- readRDS("../data_scAgeCom/analysis/outputs_data/tissue_cci_counts.rds")

TISSUE_COUNTS_SUMMARY[
  data.table(
    old_names = c(
      "calico",
      "droplet_female", "droplet_male",
      "facs_female", "facs_male"
    ),
    new_names = shiny_dataset_names
  ),
  on = "DATASET==old_names",
  DATASET := i.new_names
]

setnames(
  TISSUE_COUNTS_SUMMARY,
  colnames(TISSUE_COUNTS_SUMMARY),
  c(
    "Dataset", "Tissue",
    "Down CCI", "Flat CCI", "NSC CCI", "Up CCI",
    "Total CCI", "Total cell-types"
  )
)

setcolorder(
  TISSUE_COUNTS_SUMMARY,
  c(
    "Tissue", "Dataset", "Total cell-types",
    "Total CCI", "Flat CCI", "Down CCI", "Up CCI", 
    "NSC CCI"
  )
)

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

## All tissues of interest ####

ALL_TISSUES <- sort(unique(unlist(lapply(
  DATASETS_COMBINED,
  function(dataset) {
    dataset@cci_table_detected$ID
  }
))))

ABBR_CELLTYPE <- lapply(
  DATASETS_COMBINED,
  function(dataset) {
    dt <- unique(dataset@cci_table_detected[, c("EMITTER_CELLTYPE", "EMITTER_CELL_ABR")])
    setnames(
      dt,
      old = c("EMITTER_CELLTYPE", "EMITTER_CELL_ABR"),
      new = c("ORIGINAL_CELLTYPE", "ABBR_CELLTYPE")
    )
    dt
  }
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

## ORA dataset - tissue relationship ####

helper_dt <- rbindlist(
  lapply(
    DATASETS_COMBINED,
    function(dataset) {
      data.table(
        TISSUE = unique(dataset@cci_table_detected$ID)
      )
    }
  ),
  idcol = "DATASET"
)
helper_dt[, DATASET_TISSUE := paste(DATASET, TISSUE, sep = ":" )]

ORA_KEYWORD_TEMPLATE <- CJ(
  DATASET = c(
    "Calico2019",
    "TMS Droplet (female)", "TMS Droplet (male)",
    "TMS FACS (female)", "TMS FACS (male)"
  ),
  TISSUE = ALL_TISSUES
)
ORA_KEYWORD_TEMPLATE[, DATASET_TISSUE := paste(DATASET, TISSUE, sep = ":" )]
ORA_KEYWORD_TEMPLATE[, HAS_DATA := ifelse(DATASET_TISSUE %in% helper_dt$DATASET_TISSUE, TRUE, FALSE)]

## Loading and cleaning ORA summary ####

ORA_KEYWORD_SUMMARY <- readRDS("../data_scAgeCom/analysis/outputs_data/ora_keyword_summary.rds")
ORA_KEYWORD_SUMMARY[
  data.table(
    old_type = c("GO_TERMS", "KEGG_PWS", "LRI", "ER_CELLFAMILIES", "ER_CELLTYPES"),
    new_type = c("GO Terms", "KEGG Pathways", "LRI", "ERI Family", "ERI")
  ),
  on = c("TYPE==old_type"),
  TYPE := i.new_type
]
ORA_KEYWORD_COUNTS <- ORA_KEYWORD_SUMMARY[
  ,
  .N,
  by = c("TYPE", "REGULATION", "GENDER", "DATASET", "VALUE")
]
ORA_KEYWORD_COUNTS <- dcast.data.table(
  ORA_KEYWORD_COUNTS,
  TYPE + VALUE + REGULATION ~ DATASET + GENDER,
  value.var = "N",
  fill = 0
)[,
  c(
    "TYPE", "VALUE", "REGULATION",
    "either_either",
    "facs_male", "facs_female",
    "droplet_male", "droplet_female",
    "calico_male")
]
setnames(
  ORA_KEYWORD_COUNTS,
  old = c("either_either",
          "facs_male", "facs_female",
          "droplet_male", "droplet_female",
          "calico_male"),
  new = c("Overall (union)",
          "TMS FACS (male)", "TMS FACS (female)",
          "TMS Droplet (male)", "TMS Droplet (female)",
          "Calico2019")
)

ORA_KEYWORD_COUNTS[
  data.table(
    old_names = c("LRI", "GO_TERMS", "KEGG_PWS", "ER_CELLFAMILIES"),
    new_names = c("LRI", "GO Terms", "KEGG Pathways", "Cell-Type Families")
  ),
  on = "TYPE==old_names",
  TYPE := i.new_names
]

ORA_KEYWORD_SUMMARY_UNIQUE <- ORA_KEYWORD_SUMMARY[
  GENDER %in% c("male", "female") &
    DATASET %in% c("facs", "droplet", "calico")]
ORA_KEYWORD_SUMMARY_UNIQUE[, DATASET := paste(DATASET, GENDER, sep = "_")]
ORA_KEYWORD_SUMMARY_UNIQUE[
  data.table(
    old = c(
      "facs_male", "facs_female",
      "droplet_male", "droplet_female",
      "calico_male"),
    new = c(
      "TMS FACS (male)", "TMS FACS (female)",
      "TMS Droplet (male)", "TMS Droplet (female)",
      "Calico2019")
  ),
  on = "DATASET==old",
  DATASET := i.new
]

ALL_GLOBAL_CATEGORIES <- c(
  "LRI",
  "GO Terms",
  "KEGG Pathways",
  "ERI Family"
)

build_GLOBAL_display <- function(
  ORA_KEYWORD_COUNTS,
  global_category,
  global_type
) {
  dt <- ORA_KEYWORD_COUNTS[
    TYPE == global_category &
      REGULATION == global_type
  ][order(-`Overall (union)`)]
  setnames(dt, old = "VALUE", new = global_category)
  DT::datatable(
    data = dt[, -c(1,3)],
    options =list(
      pageLength = 10
    ),
    caption = tags$caption(
      style = 'caption-side: top; text-align: center; color:black; font-size:150% ;',
      paste0(
        "Summary over-representation for ",
        global_type,
        " by tissue and dataset."
        )
    )
  )
}

plot_keyword_tissue_vs_dataset <- function(
  ora_keyword_summary_unique,
  ora_keyword_template,
  category,
  keyword
) {
  res <- copy(ora_keyword_template)
  if (!(keyword %in% ora_keyword_summary_unique[TYPE == category]$VALUE)) {
    stop("`keyword` not found in `ora_keyword_summary_unique`")
  }
  dt <- dcast.data.table(
    ora_keyword_summary_unique[
      TYPE == category &
        VALUE == keyword],
    ID ~ DATASET,
    value.var = "REGULATION",
    fun.aggregate = paste0,
    fill = "NOT_OR",
    collapse = ":"
  )
  dt <- melt.data.table(
    dt,
    id.vars = "ID"
  )
  res[
    dt,
    on = c("DATASET==variable", "TISSUE==ID"),
    REGULATION := i.value
  ]
  res[is.na(res)] <- "NOT_OR"
  res[data.table(
    old_values = c(
      "UP", "DOWN", "FLAT",
      "UP:DOWN", "DOWN:UP",
      "UP:FLAT", "FLAT:UP",
      "DOWN:FLAT", "FLAT:DOWN",
      "NOT_OR"
    ),
    new_values = c(
      "UP", "DOWN", "FLAT",
      "UP:DOWN", "UP:DOWN",
      "UP", "UP",
      "DOWN", "DOWN",
      "NOT_OR"
    )
  ),
  on = "REGULATION==old_values",
  REGULATION := i.new_values
  ]
  res[, REGULATION := ifelse(
    HAS_DATA == FALSE,
    "NO_DATA",
    REGULATION
  )]
  setnames(
    res,
    old = c("DATASET", "TISSUE", "REGULATION"),
    new = c("Dataset", "Tissue", "Regulation")
  )
  res <- res[HAS_DATA == TRUE]
  p <- ggplot(res) +
    geom_tile(aes(
      Dataset,
      Tissue,
      fill = Regulation,
      width = 0.9,
      height = 0.9),
      colour = "black") +
    scale_fill_manual(values = c(
      "UP" = "red",
      "DOWN" = "blue",
      "FLAT" = "green",
      "UP:DOWN" = "yellow",
      "NOT_OR" = "gray", 
      "NO_DATA" = "transparent")) +
    ggtitle(paste0("Over-representation of ", keyword)) +
    scale_x_discrete(limits = c(
      "TMS FACS (male)",
      "TMS FACS (female)" ,
      "TMS Droplet (male)",
      "TMS Droplet (female)",
      "Calico2019"
    )) +
    scale_y_discrete(limits = sort(unique(res$Tissue), decreasing = TRUE)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text=element_text(size = 20)) +
    theme(axis.text=element_text(size = 18)) +
    xlab("Dataset") +
    ylab("Tissue")
  p <- plotly::ggplotly(
    p,
    source = "TCA_PLOT_KEYWORD_SUMMARY",
    tooltip = c("Dataset", "Tissue", "Regulation")
  )  %>% toWebGL()
}

## Layout config ####

shiny_layout <- list(
  color_theme = "color: rgb(20 120 206)",
  style_intro_title = paste(
    "width: 60%;",
    "margin:auto;",
    "font-size: 26px;",
    "text-align: center;"
  ),
  style_intro_text = paste(
    "width: 60%;",
    "margin:auto;",
    "font-size: 18px;",
    "text-align: justify;"
  )
)

## Shiny html ####

shiny_text <- list(
  main_title = paste(
    "A Murine ATLAS of Age-related Changes in Intercellular Communication (Alpha Version)."
  ),
  intro_title = paste(
    "Welcome to scAgeCom!"
  ),
  intro_overview_title = paste(
    "Overview of the analysis"
  )
)

shiny_html_content <- list(
  main_title = shiny_text$main_title,
  intro_title = tags$p(
    div(
      style = shiny_layout$style_intro_title,
      shiny_text$intro_title
    )
  ),
  intro_overview = tags$div(
    style = shiny_layout$style_intro_text,
    tags$h2(
      shiny_text$intro_overview_title,
      style = shiny_layout$color_theme
    ),
    tags$p(
      paste(
        "This project provides a comprehensive investigation",
        "of age-related changes in mouse intercellular communication.",
        "It combines scRNA-seq data and curated ligand-receptor",
        "interactions with a novel analysis technique that",
        "allows to statistically infere differentially expressed",
        "cell-cell interaction patterns."
      )
    ),
    tags$p(
      "Link to the future paper:..."
    ),
    tags$h3(
      "Team and Acknowledgement",
      style = shiny_layout$color_theme
    )
  ),
  intro_method = tags$div(
    style = shiny_layout$style_intro_text,
    tags$h3(
      "Code - scDiffCom",
      style = shiny_layout$color_theme
    ),
    tags$p(
      paste(
        "Inspired by several analysis techniques from",
        "the aferomentionned tools, we have build an R package,",
        "scDiffCom, that allows to statistically assess the",
        "differential expression of cell-cell interactions",
        "between two conditions of interest."
      )
    ),
    tags$p(
      "The package is currently available on GitHub",
      tags$a(
        href= "https://github.com/CyrilLagger/scDiffCom",
        "here.",
        target="_blank"
      ),
      paste(
        "It is not intended to only be used for the aging",
        "analysis presented here, but can in theory be",
        "applied to any scRNA-seq datasets."
      )
    ),
    tags$p(
      "A vignette will follow soon..."
    ),
    tags$p(
      paste(
        "The code specifically related to this aging",
        "analysis and the ShinyApp are also available on GitHub"
      ),
      tags$a(
        href= "https://github.com/CyrilLagger/scAgeCom",
        "here.",
        target="_blank"
      )
    ),
  ),
  intro_scrna_data = tags$div(
    style = shiny_layout$style_intro_text,
    tags$h3(
      "Single-cell Datasets",
      style = shiny_layout$color_theme
    ),
    tags$p(
      paste(
        "We have leveraged transcritomics single-data",
        "from two previous studies:"
      )
    ),
    tags$ol(
      tags$li(
        "Tabula Muris Senis (TMS):",
        tags$a(
          href= "https://tabula-muris-senis.ds.czbiohub.org/",
          "see their webpage",
          target="_blank"
        ),
        " and ",
        tags$a(
          href= "https://www.nature.com/articles/s41586-020-2496-1",
          "Nature article.",
          target="_blank"
        )
      ),
      tags$li(
        "Calico 2019:",
        tags$a(
          href= "https://mca.research.calicolabs.com/",
          "see their webpage",
          target="_blank"),
        " and ",
        tags$a(
          href= "https://genome.cshlp.org/content/29/12/2088",
          "Genome Research article.",
          target="_blank")
      )
    )
  ),
  intro_lri =   tags$div(
    style = "width: 60%; margin:auto; font-size: 18px; text-align: justify;",
    tags$h3(
      "Ligand-receptor Databases",
      style = shiny_layout$color_theme
    ),
    tags$p(
      "We have compared and compiled curated ligand-receptor interactions from 8 previous studies:"
    ),
    tags$ol(
      tags$li(
        "CellChat:",
        tags$a(href= "http://www.cellchat.org/", "see their webpage", target="_blank"),
        " and ",
        tags$a(href= "https://www.biorxiv.org/content/10.1101/2020.07.21.214387v1", "bioRxiv article.", target="_blank")
      ),
      tags$li(
        "CellPhoneDB:",
        tags$a(href= "https://www.cellphonedb.org/", "see their webpage", target="_blank"),
        " and ",
        tags$a(href= "https://www.nature.com/articles/s41596-020-0292-x", "Nature Protocol article.", target="_blank")
      ),
      tags$li(
        "CellTalkDB:",
        tags$a(href= "http://tcm.zju.edu.cn/celltalkdb/", "see their webpage", target="_blank"),
        " and ",
        tags$a(href= "https://academic.oup.com/bib/advance-article-abstract/doi/10.1093/bib/bbaa269/5955941", "Briefings in Bioinformatics article.", target="_blank")
      ),
      tags$li(
        "connectomeDB2020:",
        tags$a(href= "https://github.com/forrest-lab/NATMI", "see their webpage", target="_blank"),
        " and ",
        tags$a(href= "https://www.nature.com/articles/s41467-020-18873-z", "Nature Communications.", target="_blank")
      ),
      tags$li(
        "ICELLNET:",
        tags$a(href= "https://github.com/soumelis-lab/ICELLNET", "see their webpage", target="_blank"),
        " and ",
        tags$a(href= "https://www.biorxiv.org/content/10.1101/2020.03.05.976878v1", "bioRxiv article.", target="_blank")
      ),
      tags$li(
        "NicheNet:",
        tags$a(href= "https://github.com/saeyslab/nichenetr", "see their webpage", target="_blank"),
        " and ",
        tags$a(href= "https://www.nature.com/articles/s41592-019-0667-5", "Nature Methods article.", target="_blank")
      ),
      tags$li(
        "SingleCellSignalR:",
        tags$a(href= "http://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html", "see their webpage", target="_blank"),
        " and ",
        tags$a(href= "https://academic.oup.com/nar/article/48/10/e55/5810485", "Nucleic Acids Research article.", target="_blank")
      ),
      tags$li(
        "scTensor:",
        tags$a(href= "https://github.com/rikenbit/scTensor", "see their webpage", target="_blank"),
        " and ",
        tags$a(href= "https://www.biorxiv.org/content/10.1101/566182v1", "bioRxiv article.", target="_blank")
      )
    )
  )
)



## save everything for shiny ####

scAgeCom_shiny_data <- list(
  LRI_mouse_curated = LRI_mouse_curated,
  LRI_human_curated = LRI_human_curated,
  LRI_DATABASES = LRI_DATABASES,
  plot_lri_upset = plot_lri_upset,
  build_LRI_display = build_LRI_display,
  ALL_TISSUES = ALL_TISSUES,
  ALL_CELLTYPES = ALL_CELLTYPES,
  ALL_LRIs = ALL_LRIs,
  CCI_table = CCI_table,
  subset_cci_table = subset_cci_table,
  build_CCI_display = build_CCI_display,
  ALL_ORA_CATEGORIES = ALL_ORA_CATEGORIES,
  ALL_ORA_GO_ASPECTS = ALL_ORA_GO_ASPECTS,
  ALL_ORA_TYPES = ALL_ORA_TYPES,
  ORA_table = ORA_table,
  subset_ORA_table = subset_ORA_table,
  build_ORA_display = build_ORA_display,
  build_ORA_visnetwork = build_ORA_visnetwork,
  ABBR_CELLTYPE = ABBR_CELLTYPE,
  ORA_KEYWORD_COUNTS = ORA_KEYWORD_COUNTS,
  ORA_KEYWORD_SUMMARY_UNIQUE = ORA_KEYWORD_SUMMARY_UNIQUE,
  ORA_KEYWORD_TEMPLATE = ORA_KEYWORD_TEMPLATE,
  TISSUE_COUNTS_SUMMARY = TISSUE_COUNTS_SUMMARY,
  ALL_GLOBAL_CATEGORIES = ALL_GLOBAL_CATEGORIES,
  build_GLOBAL_display = build_GLOBAL_display,
  build_tissue_counts_display = build_tissue_counts_display,
  plot_keyword_tissue_vs_dataset = plot_keyword_tissue_vs_dataset,
  shiny_html_content = shiny_html_content
)

saveRDS(
  scAgeCom_shiny_data,
  "../data_scAgeCom/analysis/outputs_data/scAgeCom_shiny_data.rds"
)


