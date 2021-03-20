
## Load libraries ####

library(scDiffCom)
library(htmltools)
library(data.table)

## Loading and cleaning LRI databases ####

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

## Defining new names for Shiny ####

shiny_dataset_names <- c(
  "Calico2019",
  "TMS Droplet (female)",
  "TMS Droplet (male)",
  "TMS FACS (female)",
  "TMS FACS (male)"
)

## Loading and cleaning scDiffCom objects ####

DATASETS_COMBINED <- readRDS("../data_scAgeCom/analysis/outputs_data/scagecom_results_processed.rds")

DATASETS_COMBINED <- lapply(
  DATASETS_COMBINED,
  function(dataset) {
    dt <- dataset@cci_table_detected
    dt[, LOG2FC_BASE := LOGFC*log2(exp(1))]
    dt[, LOG2FC_ID := {
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
      )}, by = ID]
    dt[, LOG2FC_ALL := {
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
      )}]
    dataset@cci_table_detected <- dt
    return(dataset)
  }
)

names(DATASETS_COMBINED)
names(DATASETS_COMBINED) <- shiny_dataset_names

## Utilities for CCI tables ####

cols_to_show_cci_table <- c(
  "LRI",
  "Emitter Cell Type",
  "Receiver Cell Type",
  "LOG2FC",
  "Adj. p-value",
  "Age-Regulation",
  "Score Young",
  "Score Old"
)

cols_numeric_cci_table <- c(
  "LOG2FC",
  "Adj. p-value",
  "Score Young",
  "Score Old"
)

cols_old_cci_table <- c(
  "LRI",
  "EMITTER_CELLTYPE",
  "RECEIVER_CELLTYPE",
  "LOG2FC_ID",
  "BH_P_VALUE_DE",
  "REGULATION",
  "CCI_SCORE_YOUNG",
  "CCI_SCORE_OLD"
)

subset_cci_table <- function(
  cci_table,
  tissue_choice,
  emitter_choice,
  receiver_choice,
  lri_choice,
  pvalue_choice,
  log2fc_choice,
  do_order
) {
  dt <-  cci_table[ID == tissue_choice]
  if ('All LRI' %in% lri_choice) {
    dt <- dt[
      `EMITTER_CELLTYPE` %in% emitter_choice &
        `RECEIVER_CELLTYPE` %in% receiver_choice &
        `BH_P_VALUE_DE` <= pvalue_choice &
        abs(LOG2FC_ID) >= log2fc_choice
    ]
  } else {
    dt <- dt[
      `EMITTER_CELLTYPE` %in% emitter_choice &
        `RECEIVER_CELLTYPE` %in% receiver_choice &
        `LRI` %in% lri_choice &
        `BH_P_VALUE_DE` <= pvalue_choice &
        abs(LOG2FC_ID) >= log2fc_choice
    ]
  }
  if (do_order) {
    setorder(
      dt,
      -LOG2FC_ID,
      BH_P_VALUE_DE,
      -CCI_SCORE_OLD,
      -CCI_SCORE_YOUNG
    )
  }
  return(dt)
}

fix_cci_table_names <- function(
  cci_table,
  old_cols,
  new_cols
) {
  setnames(
    cci_table,
    old = old_cols,
    new = new_cols
  )
  setcolorder(cci_table, new_cols)
  cci_table <- cci_table[, new_cols, with = FALSE]
  return(cci_table)
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

## Utility functions ####

show_DT <- function(
  data,
  cols_to_show,
  cols_numeric = NULL,
  table_title = NULL,
  options = NULL,
  rownames = TRUE,
  callback = NULL
) {
  if (is.null(options)) {
    options <- list(
      pageLength = 10
    )
  }
  if (is.null(callback)) {
    callback <- htmlwidgets::JS("return table;")
  }
  res <- DT::datatable(
    data = data[, cols_to_show, with = FALSE],
    options = options,
    callback = callback,
    caption = tags$caption(style = 'caption-side: top; text-align: center; color:black; font-size:150% ;',table_title),
    rownames = rownames,
    extensions = c("Buttons")
  )
  if(!is.null(cols_numeric)) {
    res <- DT::formatSignif(
      table = res,
      columns = cols_numeric,
      digits = 3
    )
  }
  return(res)
}

show_volcano <- function(
  data,
  xlims,
  ylims
) {
  p <- ggplot(data, aes(
    x = LOG2FC,
    y = minus_log10_pval,
    color = `Age-Regulation`
  )) +
    geom_point() +
    scale_color_manual(values = c(
      "UP" = "red", "DOWN" = "blue",
      "FLAT" = "green",
      "NON_SIGNIFICANT_CHANGE" = "gray")) +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = log2(1.5)) +
    geom_vline(xintercept = -log2(1.5)) +
    xlab(expression(paste(Log[2], "FC"))) +
    ylab(expression(paste(-Log[10], " ", p[BH]))) +
    xlim(xlims) + ylim(ylims) +
    ggtitle("Volcano Plot of detected CCI (interactive)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text=element_text(size=20))
  return(p)
}

show_scores <- function(
  data,
  xlims,
  ylims
) {
  p <- ggplot(
    data,
    aes(
      x = `Score Young`,
      y = `Score Old`,
      color = `Age-Regulation`
    )
  ) +
    geom_point() +
    scale_color_manual(values = c(
      "UP" = "red", "DOWN" = "blue",
      "FLAT" = "green",
      "NON_SIGNIFICANT_CHANGE" = "gray")) +
    scale_x_log10(limits = xlims) +
    scale_y_log10(limits = ylims) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Score Young") +
    ylab("Score Old") +
    ggtitle("Old vs Young Scores of detected CCIs (interactive)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text=element_text(size=20))
  #theme(legend.title = element_blank())
  return(p)
}

show_LRIFC <- function(
  data,
  xlims,
  ylims
) {
  p <- ggplot(
    data,
    aes(
      x = LOG2FC_L,
      y = LOG2FC_R,
      color = `Age-Regulation`
    )
  ) +
    geom_point() +
    scale_color_manual(values = c(
      "UP" = "red", "DOWN" = "blue",
      "FLAT" = "green",
      "NON_SIGNIFICANT_CHANGE" = "gray")) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlim(xlims) + ylim(ylims) +
    xlab("Ligand LOG2FC") +
    ylab("Receptor LOG2FC") +
    ggtitle("Ligand vs Receptor Fold-Change of detected CCIs (interactive)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text=element_text(size=20))
  #theme(legend.title = element_blank())
  return(p)
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
  p <- ggplot(res, aes(DATASET, TISSUE)) +
    geom_tile(aes(
      fill = REGULATION,
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
    scale_y_discrete(limits = sort(unique(res$TISSUE), decreasing = TRUE)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text=element_text(size = 20)) +
    theme(axis.text=element_text(size = 18)) +
    xlab("Dataset") +
    ylab("Tissue")
  return(p)
}

plot_keyword_tissue_vs_dataset(
  ORA_KEYWORD_SUMMARY_UNIQUE,
  ORA_KEYWORD_TEMPLATE,
  "GO Terms",
  "extracellular matrix organization"
)

plot_keyword_tissue_vs_dataset(
  ORA_KEYWORD_SUMMARY_UNIQUE,
  ORA_KEYWORD_TEMPLATE,
  "LRI",
  "App:Lrp10"
)


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
    "A Murine ATLAS of Age-related Intercellular Communication Changes."
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
  LRI_DATABASES = LRI_DATABASES,
  plot_lri_upset = plot_lri_upset,
  ALL_TISSUES = ALL_TISSUES,
  DATASETS_COMBINED = DATASETS_COMBINED,
  subset_cci_table = subset_cci_table,
  fix_cci_table_names = fix_cci_table_names,
  cols_to_show_cci_table = cols_to_show_cci_table,
  cols_numeric_cci_table = cols_numeric_cci_table,
  cols_old_cci_table = cols_old_cci_table,
  ABBR_CELLTYPE = ABBR_CELLTYPE,
  ORA_KEYWORD_COUNTS = ORA_KEYWORD_COUNTS,
  ORA_KEYWORD_SUMMARY_UNIQUE = ORA_KEYWORD_SUMMARY_UNIQUE,
  ORA_KEYWORD_TEMPLATE = ORA_KEYWORD_TEMPLATE,
  TISSUE_COUNTS_SUMMARY = TISSUE_COUNTS_SUMMARY,
  plot_keyword_tissue_vs_dataset = plot_keyword_tissue_vs_dataset,
  show_DT = show_DT,
  show_LRIFC = show_LRIFC,
  show_scores = show_scores,
  show_volcano = show_volcano,
  shiny_html_content = shiny_html_content
)

saveRDS(
  scAgeCom_shiny_data,
  "../data_scAgeCom/analysis/outputs_data/scAgeCom_shiny_data.rds"
)


