## Libraries ####

library(shiny)
library(shinyWidgets)
library(DT)
library(data.table)
library(ggplot2)
library(ComplexUpset)
library(scDiffCom)
library(ggExtra)
library(visNetwork)

## Global Data ####

LRdb_curated <- readRDS("data/LRdb_mouse_curated.rds")
cols_to_show_LRdb <- c(
  "LIGAND_1", "LIGAND_2",
  "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
  "DATABASE", "SOURCE"
)

DATASETS_COMBINED <- readRDS("data/analysis4_DATASETS_COMBINED_log15_light_new_ora.rds")
DATASETS_COMBINED <- lapply(
  DATASETS_COMBINED,
  function(dataset) {
    dt <- dataset@cci_detected
    dt[, LOG2FC := LOGFC*log2(exp(1))]
    dataset@cci_detected <- dt
    return(dataset)
  }
)
names(DATASETS_COMBINED) <- c("Calico Data", "TMS Droplet Data", "TMS FACS Data")

CCI_SUMMARY <- readRDS("data/analysis5_SUMMARY_log15.rds")
names(CCI_SUMMARY) <- c("Calico Data", "TMS Droplet Data", "TMS FACS Data", "TMS FACS SUB")

## Global functions ####

show_DT <- function(
  data,
  cols_to_show,
  cols_numeric,
  table_title = NULL
) {
  res <- DT::datatable(
    data = data[, cols_to_show, with = FALSE],
    options = list(
      pageLength = 10
    ),
    caption = tags$caption( style = 'caption-side: top; text-align: center; color:black; font-size:200% ;',table_title) 
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
  data
) {
  p <- ggplot(data, aes(
    x = LOG2FC,
    y = mlog10_pval,
    color = ifelse(
      REGULATION_SIMPLE=="UP",
      "red",
      ifelse(
        REGULATION_SIMPLE=="DOWN",
        "blue",
        "green"
      )
    )
    #color = ifelse(
    #  `Adj. P-Value` <= 0.05 & LOG2FC >= log2(1.2),
    #  "red",
    #  ifelse(
    #    `Adj. P-Value` <= 0.05 & LOG2FC <= -log2(1.2),
    #    "blue",
    #    "green"
    #  )
    #)
  )) +
    geom_point() +
    scale_color_identity() +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = log2(1.5)) +
    geom_vline(xintercept = -log2(1.5)) +
    xlab(expression(paste(Log[2], "FC"))) +
    ylab(expression(paste(-Log[10], " ", p[BH]))) +
    ggtitle("Volcano Plot of all detected CCIs") +
    theme(text=element_text(size=20))
  #ggMarginal(
  #  p = p,
  #  type =  "histogram",
  #  margins = "x",
  #  xparams = list(bins = 100)
  #)
  p
}

show_scores <- function(
  data
) {
  p <- ggplot(
    data,
    aes(
      x = `CCI_SCORE_YOUNG`,
      y = `CCI_SCORE_OLD`,
      #color = ifelse(
      #  `Adj. P-Value` <= 0.05 & LOG2FC >= log2(1.2),
      #  "red",
      #  ifelse(
      #    `Adj. P-Value` <= 0.05 & LOG2FC <= -log2(1.2),
      #    "blue",
      #    "green"
      #  )
      #),
      color = ifelse(
        REGULATION_SIMPLE=="UP",
        "red",
        ifelse(
          REGULATION_SIMPLE=="DOWN",
          "blue",
          "green"
        )
      )
    )
  ) + 
    geom_point() +
    scale_color_identity() +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(slope = 1, intercept = 0) +
    xlab("CCI_SCORE_YOUNG") +
    ylab("CCI_SCORE_OLD") +
    ggtitle("Age Scores of all detected CCIs") +
    theme(text=element_text(size=20))
  
  return(p)
}