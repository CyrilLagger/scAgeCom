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
library(htmlwidgets)

## Global Data ####

LRdb_curated <- copy(scDiffCom::LRdb_mouse$LRdb_curated)
LRdb_curated[, SOURCE := gsub("SCT:SingleCellSignalR|SCT:CellPhoneDB|SCT:DLRP", "", SOURCE)]
LRdb_curated[, SOURCE := gsub(";;", ";", SOURCE)]
LRdb_curated[, SOURCE := sub(";$", "", SOURCE)]
LRdb_curated[, SOURCE := sub("^;", "", SOURCE)]

#DATASETS_COMBINED <- readRDS("data/analysis4_DATASETS_COMBINED_log15_light.rds")
DATASETS_COMBINED <- readRDS("data/TMS_scAgeCom_processed.rds")
DATASETS_COMBINED <- DATASETS_COMBINED[c(1,2,3,4,6,7,8)]
DATASETS_COMBINED <- lapply(DATASETS_COMBINED, function(i) i$dataset)
DATASETS_COMBINED <- lapply(
  DATASETS_COMBINED,
  function(dataset) {
    dt <- dataset@cci_detected
    dt[, LOG2FC := LOGFC*log2(exp(1))]
    dt[, LOG2FC := {
      temp <- LOG2FC
      temp_max <- round(max(temp[is.finite(temp)])) + 1
      temp_min <- round(min(temp[is.finite(temp)])) - 1
      ifelse(
        is.infinite(LOG2FC) & LOG2FC > 0,
        temp_max,
        ifelse(
          is.infinite(LOG2FC) & LOG2FC < 0,
          temp_min,
          LOG2FC
        )
      )}]
    dataset@cci_detected <- dt
    return(dataset)
  }
)
#names(DATASETS_COMBINED) <- c("Calico Data", "TMS Droplet Data", "TMS FACS Data")
names(DATASETS_COMBINED) <- c("Calico Data",
                              "TMS Droplet Data Female", "TMS Droplet Data Male", "TMS Droplet Data Mixed",
                              "TMS FACS Data Female", "TMS FACS Data Male", "TMS FACS Data Mixed")


CCI_SUMMARY <- readRDS("data/analysis5_SUMMARY_DATA.rds")
CCI_SUMMARY <- CCI_SUMMARY[c(1,2,3,4,6,7,8)]
#names(CCI_SUMMARY) <- c("Calico Data", "TMS Droplet Data", "TMS FACS Data")#, "TMS FACS Data (Sub)")
names(CCI_SUMMARY) <- names(DATASETS_COMBINED)

## Global functions ####

show_DT <- function(
  data,
  cols_to_show,
  cols_numeric = NULL,
  table_title = NULL,
  options = NULL,
  rownames = TRUE
) {
  if (is.null(options)) {
    options <- list(pageLength = 10)
  }
  res <- DT::datatable(
    data = data[, cols_to_show, with = FALSE],
    options = options,
    caption = tags$caption(style = 'caption-side: top; text-align: center; color:black; font-size:200% ;',table_title),
    rownames = rownames
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
    color = REGULATION
  )) +
    geom_point() +
    scale_color_manual(values = c(
      "UP" = "red", "DOWN" = "blue",
      "FLAT" = "green",
      "NON_SIGNIFICANT_CHANGE" = "black")) +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = log2(1.5)) +
    geom_vline(xintercept = -log2(1.5)) +
    xlab(expression(paste(Log[2], "FC"))) +
    ylab(expression(paste(-Log[10], " ", p[BH]))) +
    ggtitle("Interactive Volcano Plot of detected CCI") +
    theme(text=element_text(size=20)) +
    theme(legend.title = element_blank()) #+
    #theme(legend.position = c(0.8, 0.2))
  #ggMarginal(
  #  p = p,
  #  type =  "histogram",
  #  margins = "x",
  #  xparams = list(bins = 100)
  #)
  return(p)
}

show_scores <- function(
  data
) {
  p <- ggplot(
    data,
    aes(
      x = `CCI_SCORE_YOUNG` + 1E-3,
      y = `CCI_SCORE_OLD` + 1E-3,
      color = REGULATION
    )
  ) + 
    geom_point() +
    scale_color_manual(values = c(
      "UP" = "red", "DOWN" = "blue",
      "FLAT" = "green",
      "NON_SIGNIFICANT_CHANGE" = "black")) +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Young Score") +
    ylab("Old Score") +
    ggtitle("Interactive Score Plot of detected CCIs") +
    theme(text=element_text(size=20)) +
    theme(legend.title = element_blank())
  return(p)
}