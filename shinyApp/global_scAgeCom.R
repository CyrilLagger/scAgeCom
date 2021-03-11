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

#DATASETS_COMBINED <- readRDS("data/TMS_scAgeCom_processed_testing.rds")
DATASETS_COMBINED <- readRDS("data/TMS_scAgeCom_processed.rds")
DATASETS_COMBINED <- DATASETS_COMBINED[c(1,2,3,4,6,7,8)]
DATASETS_COMBINED <- lapply(DATASETS_COMBINED, function(i) i$dataset)
DATASETS_COMBINED <- lapply(
  DATASETS_COMBINED,
  function(dataset) {
    dt <- dataset@cci_detected
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
    dataset@cci_detected <- dt
    return(dataset)
  }
)

#names(DATASETS_COMBINED) <- c("Calico Data", "TMS Droplet Data", "TMS FACS Data")
names(DATASETS_COMBINED) <- c("Calico Data",
                              "TMS Droplet Data Female", "TMS Droplet Data Male", "TMS Droplet Data Mixed",
                              "TMS FACS Data Female", "TMS FACS Data Male", "TMS FACS Data Mixed")

ABBR <- lapply(
  DATASETS_COMBINED,
  function(dataset) {
    dt <- unique(dataset@cci_detected[, c("EMITTER_CELLTYPE", "EMITTER_CELL_ABR")])
    setnames(
      dt,
      old = c("EMITTER_CELLTYPE", "EMITTER_CELL_ABR"),
      new = c("ORIGINAL_CELLTYPE", "ABBR_CELLTYPE")
    )
    dt
  }
)



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
    color = `Age Regulation`
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
    theme(text=element_text(size=20)) #+
   # theme(legend.title = element_blank())
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
      color = `Age Regulation`
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
      color = `Age Regulation`
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
