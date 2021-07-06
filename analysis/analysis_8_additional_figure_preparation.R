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
## additional figure preparation
##
####################################################
##

## Libraries ####

library(Seurat)
library(data.table)
library(scDiffCom)
library(ggplot2)
library(plotly)
library(shiny)
library(kableExtra)

## Fig.2.1Data-preprocessing table ####

fig2_sc_data <- expm1(
  scDiffCom::seurat_sample_tms_liver$RNA@data[1:4, 1:4]
)
fig2_sc_data <- as.data.frame(
  as.matrix(
    fig2_sc_data
  )
)
fig2_sc_data <- round(
  fig2_sc_data,
  1
)
fig2_sc_data <- data.frame(
  lapply(
    fig2_sc_data,
    as.character
  ),
  stringsAsFactors = FALSE
)
colnames(fig2_sc_data) <- c(
  paste(
    "Cell",
    1:(ncol(fig2_sc_data)-2)
  ),
  "...", "Cell N"
)
rownames(fig2_sc_data) <- c(
  paste(
    "Gene",
    1:(nrow(fig2_sc_data)-2)
  ),
  "...", "Gene M"
)
fig2_sc_data$... <- "..."
fig2_sc_data[3, 1:4] <- "..."

group <- paste(
  seurat_sample_tms_liver$cell_type,
  seurat_sample_tms_liver$age_group,
  sep = "_"
)

fig2_sc_aggr <- t(
  DelayedArray::rowsum(
    Matrix::t(scDiffCom::seurat_sample_tms_liver$RNA@data),
    group = group
  ) / as.vector(table(group))
)
fig2_sc_aggr <- fig2_sc_aggr[1:4, c(1,2,3,7,8)]
fig2_sc_aggr <- round(
  fig2_sc_aggr,
  1
)
colnames(fig2_sc_aggr) <- c("Cond1", "Cond2", "...", "Cond1", "Cond2")
rownames(fig2_sc_aggr ) <- c(
  paste(
    "Gene",
    1:(nrow(fig2_sc_aggr)-2)
  ),
  "...",
  "Gene M"
)
fig2_sc_aggr[3,] <- "..."
fig2_sc_aggr[, 3] <- "..."


cbind(
  fig2_sc_data,
  fig2_sc_aggr
) %>% kbl(
  caption = "<span style = 'font-size: 50px;font-weight: bold;'>Averaged expression</span>",
  align = rep("c", 5)
) %>% kable_styling(
  "striped",
  full_width = FALSE
) %>% add_header_above(
  c(" " = 5, "Cell Type 1" = 2, " " = 1, "Cell Type 4" = 2)
) %>% kable_styling(
  font_size = 50
)

fig2_sc_aggr %>%
  kbl(
    caption = "<span style = 'font-size: 70px;font-weight: bold;'>Averaged expression</span>",
    align = rep("c", 5)
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(c(" " = 1, "Cell Type 1" = 2, " " = 1, "Cell Type 4" = 2)) %>%
  kable_styling(font_size = 60) %>%
  column_spec(1:6, bold = T) %>%
  save_kable(file = "../data_scAgeCom/figures/fig2_sc_aggr.png", zoom = 2, vwidth = 1200)


fig2_sc_data %>%
  kbl(
    caption = "<span style = 'font-size: 70px; font-weight: bold;'>Normalized counts</span>",
    align = rep("c", 4)
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  kable_styling(font_size = 60) %>%
  column_spec(1:5, bold = T) %>%
  save_kable(file = "../data_scAgeCom/figures/fig2_sc_data.png",  zoom = 2, vwidth = 1200)

## Fig.2.2 LRI table  ####

fig2_LRI_data <- as.data.frame(scDiffCom::LRI_mouse$LRI_curated[c(1, 643, 2, 3573), c(2, 3, 1, 4,5, 6)])
fig2_LRI_data[is.na(fig2_LRI_data)] <- ""
fig2_LRI_data[c(1,3),] <- "..."
fig2_LRI_data[, 3] <- "\u00A0 \u00A0"
colnames(fig2_LRI_data) <- c("L1", "L2", "\u00A0 \u00A0  ", "R1", "R2", "R3")


fig2_LRI_data %>% 
  kbl(
    caption = "<span style = 'font-size: 50px; font-weight: bold'>Ligand-Receptor Interactions</span>",
    align = rep("l", 5)
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(c("Ligand" = 2, "\u00A0" = 1, "Receptor" = 3)) %>%
  kable_styling(font_size = 50) %>%
  column_spec(1:5, bold = T) %>%
  column_spec(3, width = "5cm") %>%
  save_kable(file = "../data_scAgeCom/figures/fig2_LRI_data.png", zoom = 3)

## Fig.2.3 Potential CCI table ####

fig2_cci_table_pot <- copy(readRDS( "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds")$CCI_table)[Dataset == "TMS FACS (female)" & Tissue == "Mammary_Gland"]
fig2_cci_table_pot[
  data.table(
    old_ct = unique(fig2_cci_table_pot$`Emitter Cell Type`),
    new_ct = paste0("Cell Type ", 1:length( unique(fig2_cci_table_pot$`Emitter Cell Type`)))
  ),
  on = "`Emitter Cell Type`==old_ct",
  emitter_cell_type := new_ct
]
fig2_cci_table_pot[
  data.table(
    old_ct = unique(fig2_cci_table_pot$`Receiver Cell Type`),
    new_ct = paste0("Cell Type ", 1:length( unique(fig2_cci_table_pot$`Receiver Cell Type`)))
  ),
  on = "`Receiver Cell Type`==old_ct",
  receiver_cell_type := new_ct
]
fig2_cci_table_pot <- fig2_cci_table_pot[sample(1:nrow(fig2_cci_table_pot), 10), c(24, 25, 3, 9, 10, 6)]
setorder(
  fig2_cci_table_pot,
  emitter_cell_type,
  receiver_cell_type
)

fig2_cci_table_pot <- fig2_cci_table_pot[c(1,2,6,7,10)]
fig2_cci_table_pot$`Log2 FC` <- fig2_cci_table_pot$`Log2 FC`/log2(exp(1))
fig2_cci_table_pot[, 4:6] <- round(fig2_cci_table_pot[, 4:6] , 1)
fig2_cci_table_pot <- data.frame(lapply(fig2_cci_table_pot, as.character), stringsAsFactors = FALSE)
fig2_cci_table_pot[c(1,3,5), ] <- "..."
colnames(fig2_cci_table_pot) <- c("Emitter", "Receiver", "LRI", "Score Cond1", "Score Cond2", "Log FC")

fig2_cci_table_pot

fig2_cci_table_pot %>%
  kbl(
    caption = "<span style = 'font-size: 70px;font-weight: bold;'>Hypothetic cell-cell interactions</span>",
    align = rep("c", 6)
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  kable_styling(font_size = 55) %>%
  column_spec(1:6, bold = T) %>%
  save_kable(file = "../data_scAgeCom/figures/fig2_cci_table_pot.png", zoom = 2, vwidth = 2200)


## Fig.2.1 tsne plots by cell types and age ####

seurat_mg <- readRDS("../data_scAgeCom/analysis/inputs_data/seurat_testing_tms_facs_mammary_gland.rds")
seurat_mg$age_group <- ifelse(
  seurat_mg$age_group == "YOUNG",
  "Cond1",
  "Cond2"
)
seurat_mg$cell_ontology_final <- ifelse(
  seurat_mg$cell_ontology_final == "basal cell",
  "Cell Type 1",
  ifelse(
    seurat_mg$cell_ontology_final == "stromal cell",
    "Cell Type 3",
    ifelse(
      seurat_mg$cell_ontology_final == "endothelial cell",
      "Cell Type 4",
      "Cell Type 2"
    )
  )
)


seurat_mg <- Seurat::RunTSNE(seurat_mg, reduction = "PCA")

Idents(seurat_mg) <- seurat_mg$cell_ontology_final

Seurat::DimPlot(
  seurat_mg,
  reduction = "tsne",
  pt.size = 1
) + theme(text=element_text(size=40))

fig2_tsne_celltype <- Seurat::DimPlot(
  seurat_mg,
  reduction = "tsne",
  pt.size = 3
) + theme(text=element_text(size=40))

ggsave("../data_scAgeCom/figures/fig2_tsne_celltype.png", fig2_tsne_celltype, scale = 1)

Idents(seurat_mg) <- seurat_mg$age_group

fig2_tsne_cond <- Seurat::DimPlot(
  seurat_mg,
  reduction = "tsne",
  pt.size = 1
) + theme(text=element_text(size=40))

ggsave("../data_scAgeCom/figures/fig2_tsne_cond.png", fig2_tsne_cond)

## Toy model simulation for distributions #####

aging_example <- run_interaction_analysis(
  seurat_object = seurat_sample_tms_liver,
  LRI_species = "mouse",
  seurat_celltype_id = "cell_type",
  seurat_condition_id = list(
    column_name = "age_group",
    cond1_name = "YOUNG",
    cond2_name = "OLD"
  ),
  return_distributions = TRUE
)

## Fig.2.4 Differential permutations ####

distr_de <- data.frame(
  counts = GetDistributions(aging_example)$DISTRIBUTIONS_DE[235, ]
)
hist(distr_de$counts)

distr_de$is_above <- ifelse(
  abs(distr_de$counts) >= GetDistributions(aging_example)$DISTRIBUTIONS_DE[235, 1001],
  TRUE,
  FALSE
)

fig2_distr_de <- ggplot(
  distr_de,
  aes(
    x = counts,
    fill = is_above
  )
) + geom_histogram(
  bins = 40,
  alpha = 0.7,
  color = "black",
  show.legend = FALSE
) + scale_fill_manual(
  values = c("blue", "red")
) + annotate(
  "segment",
  x = GetDistributions(aging_example)$DISTRIBUTIONS_DE[235, 1001],
  xend = GetDistributions(aging_example)$DISTRIBUTIONS_DE[235, 1001],
  y = 38,
  yend = 0.5,
  size = 3,
  arrow = arrow(length = unit(0.5, "cm"))
) + annotate(
  "text",
  x = GetDistributions(aging_example)$DISTRIBUTIONS_DE[235, 1001]+0.8,
  y = 42,
  label = "True Difference",
  size = 22
) + xlab(
  "Score(Cond2) - Score(Cond1)"
) + ylab(
  "Frequency"
) + theme(
  text=element_text(size=80),
  axis.text.x = element_text(size = 80),
  axis.text.y = element_blank()
) 
fig2_distr_de
#manual save 2000x1400
#ggsave("../data_scAgeCom/figures/fig2_distr_de.png", fig2_distr_de, scale = 1.5)

## Fig.2.4 Cond1 permutations ####

distr_cond1 <- data.frame(
  counts = GetDistributions(aging_example)$DISTRIBUTIONS_YOUNG[235, ]
)
hist(distr_cond1$counts)
distr_cond1$counts[1001]

distr_cond1$is_above <- ifelse(
  abs(distr_cond1$counts) >= GetDistributions(aging_example)$DISTRIBUTIONS_YOUNG[235, 1001],
  TRUE,
  FALSE
)

fig2_distr_cond1 <- ggplot(
  distr_cond1,
  aes(
    x = counts,
    fill = is_above
  )
) + geom_histogram(
  bins = 50,
  alpha = 0.7,
  color = "black",
  show.legend = FALSE
) + scale_fill_manual(
  values = c("blue", "red")
) + annotate(
  "segment",
  x = GetDistributions(aging_example)$DISTRIBUTIONS_YOUNG[235, 1001],
  xend = GetDistributions(aging_example)$DISTRIBUTIONS_YOUNG[235, 1001],
  y = 22,
  yend = 0.5,
  size = 3,
  arrow = arrow(length = unit(0.5, "cm"))
) + annotate(
  "text",
  x = GetDistributions(aging_example)$DISTRIBUTIONS_YOUNG[235, 1001],
  y = 25,
  label = "True Score 1",
  size = 22
) + xlab(
  "Score(Cond1)"
) + ylab(
  "Frequency"
) + theme(
  text=element_text(size=80),
  axis.text.x = element_text(size = 80),
  axis.text.y = element_blank()
) + xlim(-0.5, 10.5)
fig2_distr_cond1
#manual save 2000x1400

#ggsave("../data_scAgeCom/figures/fig2_distr_cond1.png", fig2_distr_cond1, scale = 1.5)


## Fig.2.4 Cond2 permutations ####

distr_cond2 <- data.frame(
  counts = GetDistributions(aging_example)$DISTRIBUTIONS_OLD[235, ]
)

distr_cond2$is_above <- ifelse(
  abs(distr_cond2$counts) >= GetDistributions(aging_example)$DISTRIBUTIONS_OLD[235, 1001],
  TRUE,
  FALSE
)

fig2_distr_cond2 <- ggplot(
  distr_cond2,
  aes(
    x = counts,
    fill = is_above
  )
) + geom_histogram(
  bins = 50,
  alpha = 0.7,
  color = "black",
  show.legend = FALSE
) + scale_fill_manual(
  values = c("blue", "red")
) + annotate(
  "segment",
  x = GetDistributions(aging_example)$DISTRIBUTIONS_OLD[235, 1001],
  xend = GetDistributions(aging_example)$DISTRIBUTIONS_OLD[235, 1001],
  y = 22,
  yend = 0.5,
  size = 3,
  arrow = arrow(length = unit(0.5, "cm"))
) + annotate(
  "text",
  x = GetDistributions(aging_example)$DISTRIBUTIONS_OLD[235, 1001]-0.3,
  y = 25,
  label = "True Score 2",
  size = 22
) + xlab(
  "Score(Cond2)"
) + ylab(
  "Frequency"
) + theme(
  text=element_text(size=80),
  axis.text.x = element_text(size = 80),
  axis.text.y = element_blank()
) + xlim(-0.5, 13)
fig2_distr_cond2

ggsave("../data_scAgeCom/figures/fig2_distr_cond2.png", fig2_distr_cond2, scale = 1.5)


## Fig.2.6 Final CCI table ####

fig2_cci_table_final <- copy(readRDS( "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds")$CCI_table)[Dataset == "TMS FACS (female)" & Tissue == "Mammary_Gland"]
fig2_cci_table_final[
  data.table(
    old_ct = unique(fig2_cci_table_final$`Emitter Cell Type`),
    new_ct = paste0("Cell Type ", 1:length( unique(fig2_cci_table_final$`Emitter Cell Type`)))
  ),
  on = "`Emitter Cell Type`==old_ct",
  emitter_cell_type := new_ct
]
fig2_cci_table_final[
  data.table(
    old_ct = unique(fig2_cci_table_final$`Receiver Cell Type`),
    new_ct = paste0("Cell Type ", 1:length( unique(fig2_cci_table_final$`Receiver Cell Type`)))
  ),
  on = "`Receiver Cell Type`==old_ct",
  receiver_cell_type := new_ct
]
fig2_cci_table_final <- fig2_cci_table_final[sample(1:nrow(fig2_cci_table_final), 10), c(24, 25, 3, 6, 7, 8)]
setorder(
  fig2_cci_table_final,
  emitter_cell_type,
  receiver_cell_type
)
fig2_cci_table_final$`Log2 FC` <- fig2_cci_table_final$`Log2 FC`/log2(exp(1))
fig2_cci_table_final[, 4] <- round(fig2_cci_table_final[, 4] , 1)
fig2_cci_table_final[, 5] <- round(fig2_cci_table_final[, 5] , 2)
fig2_cci_table_final <- data.frame(lapply(fig2_cci_table_final, as.character), stringsAsFactors = FALSE)

fig2_cci_table_final <- fig2_cci_table_final[c(1,4,5,9,10),]

fig2_cci_table_final[c(1,3,5), ] <- "..."
colnames(fig2_cci_table_final) <- c("Emitter", "Receiver", "LRI", "Log FC", "Adj. p-value", "Regulation")

fig2_cci_table_final

fig2_cci_table_final %>%
  kbl(
    caption = "<span style = 'font-size: 70px;font-weight: bold;'>Detected cell-cell interactions</span>",
    align = rep("c", 6),
    row.names = FALSE
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  kable_styling(font_size = 55) %>%
  column_spec(1:6, bold = T) %>%
  save_kable(file = "../data_scAgeCom/figures/fig2_cci_table_final.png", zoom = 2, vwidth = 2100)





## Figure 4 ####

CCI_table <- readRDS( "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds")$CCI_table
ORA_table <- readRDS( "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds")$ORA_table
ABBR_CELLTYPE <- readRDS( "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds")$ABBR_CELLTYPE

fig4_cci <- CCI_table[
  Dataset == "TMS Droplet (male)" &
    Tissue == "Bladder"
]

plot_volcano_CCI <- function(
  CCI_table
) {
  dt <- CCI_table[
    ,
    3:12
  ]
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
    l = 10,
    r = 10,
    b = 10,
    t = 30,
    pad = 10
  )
  plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~`Log2 FC`,
    y = ~-log10(`Adj. p-value` + 1E-4),
    text = ~paste(
      "LRI: ",
      `Ligand-Receptor Interaction`, 
      '<br>Emitter:',
      `Emitter Cell Type`,
      '<br>Receiver:',
      `Receiver Cell Type`
    ),
    color = ~`Age Regulation`,
    colors = stats::setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    ),
    marker = list(size = 15)
  ) %>% plotly::layout(
    title = list(
      text = "",
      font = list(size = 20),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = "Log2(FC)",
        font = list(size = 36)
      ),
      tickfont = list(
        size = 24
      )
    ),
    yaxis = list(
      title = list(
        text = "-Log10(Adj. p-value)",
        font = list(size = 36),
        standoff = 40
      ),
      tickfont = list(
        size = 28
      )
    ),
    shapes = list(
      vline(log2(1.5)),
      vline(-log2(1.5)),
      hline(-log10(0.05))
    ),
    legend = list(
      orientation = "h",
      xanchor = "center",
      x = 0.5,
      y = 1.1,
      font = list(size = 34)
    ),
    margin = m
  ) 
}

plot_volcano_CCI(fig4_cci)

plot_scores_CCI <- function(
  CCI_table
) {
  `Young CCI Score` <- `Old CCI Score` <- NULL
  dt <- CCI_table[
    ,
    3:12
  ]
  dt$`Age Regulation` <- factor(
    dt$`Age Regulation`,
    levels = c("UP", "DOWN", "FLAT", "NSC")
  )
  min_score <-  10^(floor(
    log10(
      min(
        min(dt[`Young CCI Score` > 0]$`Young CCI Score`),
        min(dt[`Old CCI Score` > 0]$`Old CCI Score`)
      )
    )
  ))
  m <- list(
    l = 10,
    r = 10,
    b = 10,
    t = 30,
    pad = 10
  )
  plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~log10(`Young CCI Score` + min_score),
    y = ~log10(`Old CCI Score` + min_score),
    text = ~paste(
      "LRI: ",
      `Ligand-Receptor Interaction`, 
      '<br>Emitter:',
      `Emitter Cell Type`,
      '<br>Receiver:',
      `Receiver Cell Type`
    ),
    color = ~`Age Regulation`,
    colors = stats::setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    ),
    marker = list(size = 15)
  ) %>% plotly::layout(
    title = list(
      text = "Interactive Score Plot",
      font = list(size = 16),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = "Log10(Young CCI Score)",
        font = list(size = 32)
      ),
      tickfont = list(
        size = 24
      )
    ),
    yaxis = list(
      title = list(
        text = "Log10(Old CCI Score)",
        font = list(size = 32)
      ),
      tickfont = list(
        size = 24
      )
    ),
    legend = list(
      orientation = "h",
      xanchor = "center",
      x = 0.5,
      y = 1.02,
      font = list(size = 28)
    ),
    margin = m
  )
}

plot_scores_CCI(fig4_cci)

plot_lrfc_CCI <- function(
  CCI_table
) {
  dt <- CCI_table[
    ,
    3:12
  ]
  dt$`Age Regulation` <- factor(
    dt$`Age Regulation`,
    levels = c("UP", "DOWN", "FLAT", "NSC")
  )
  m <- list(
    l = 10,
    r = 10,
    b = 10,
    t = 30,
    pad = 10
  )
  plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~`Ligand Log2 FC`,
    y = ~`Receptor Log2 FC`,
    text = ~paste(
      "LRI: ",
      `Ligand-Receptor Interaction`, 
      '<br>Emitter:',
      `Emitter Cell Type`,
      '<br>Receiver:',
      `Receiver Cell Type`
    ),
    color = ~`Age Regulation`,
    colors = stats::setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    ),
    marker = list(size = 15)
  )  %>% plotly::layout(
    title = list(
      text = "Interactive 'Ligand-FC vs Receptor-FC' Plot",
      font = list(size = 16),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = "Ligand Log2(FC)",
        font = list(size = 32)
      ),
      tickfont = list(
        size = 24
      )
    ),
    yaxis = list(
      title = list(
        text = "Receptor Log2(FC)",
        font = list(size = 32)
      ),
      tickfont = list(
        size = 24
      )
    ),
    legend = list(
      orientation = "h",
      xanchor = "center",
      x = 0.5,
      y = 1.02,
      font = list(size = 28)
    ),
    margin = m
  )
}

plot_lrfc_CCI(fig4_cci)

display_CCI_table <- function(
  CCI_table
) {
  dt <- CCI_table[
    ,
    3:12
  ]
  CCI_DT <- DT::datatable(
    data = dt[, -c(9, 10)],
    class = "display compact",
    options =list(
      pageLength = 10,
      dom = '<"top"f>rt<"bottom"lip><"clear">'
    ),
    caption = tags$caption(
      style = paste0(
        'caption-side: top; text-align: center; ',
        'color:black; font-size:120% ;'
      ),
      "Table of Cell-Cell Interactions"
    ),
    rownames = rownames,
    extensions = c("Buttons")
  ) %>% DT::formatStyle(
    colnames(dt[, -c(9, 10)])[4:8],
    `text-align` = 'center'
  )
}

fig4_cci_table <- display_CCI_table(fig4_cci)

plot_ORA_visnetwork <- function(
  CCI_table,
  ORA_table,
  tissue_choice,
  dataset_choice,
  abbr_celltype
) {
  Dataset <- Tissue <- ORA_CATEGORY <- ORIGINAL_CELLTYPE <-
    EMITTER_CELLTYPE <- RECEIVER_CELLTYPE <- i.ABBR_CELLTYPE <- NULL
  CCI_dt <- copy(CCI_table)
  data.table::setnames(
    CCI_dt,
    old = c("Emitter Cell Type", "Receiver Cell Type", "Age Regulation"),
    new = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION")
  )
  ora_table_ER <- ORA_table[
    Dataset == dataset_choice &
      Tissue == tissue_choice &
      ORA_CATEGORY == "ER_CELLTYPES"
  ]
  ora_table_EMITTER <- ORA_table[
    Dataset == dataset_choice &
      Tissue == tissue_choice &
      ORA_CATEGORY == "EMITTER_CELLTYPE"
  ]
  ora_table_RECEIVER <- ORA_table[
    Dataset == dataset_choice &
      Tissue == tissue_choice &
      ORA_CATEGORY == "RECEIVER_CELLTYPE"
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
      ,
      "EMITTER_CELLTYPE_ORIGINAL" := EMITTER_CELLTYPE
    ]
    cci_table_detected[
      ,
      "RECEIVER_CELLTYPE_ORIGINAL" := RECEIVER_CELLTYPE
    ]
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
    ora_table_EMITTER[
      abbreviation_table,
      on = "VALUE==ORIGINAL_CELLTYPE",
      "VALUE" := i.ABBR_CELLTYPE]
    ora_table_RECEIVER[
      abbreviation_table,
      on = "VALUE==ORIGINAL_CELLTYPE",
      "VALUE" := i.ABBR_CELLTYPE]
  }
  scDiffCom:::interactive_from_igraph(
    cci_table_detected = cci_table_detected,
    conds = c("YOUNG", "OLD"),
    ora_table_ER = ora_table_ER,
    ora_table_EMITTER = ora_table_EMITTER,
    ora_table_RECEIVER = ora_table_RECEIVER,
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

plot_ORA_visnetwork(
  CCI_table,
  ORA_table,
  "Bladder",
  "TMS Droplet (male)",
  ABBR_CELLTYPE
)

plot_ora_local <- function(
  ora_dt,
  category,
  regulation,
  max_terms_show,
  GO_aspect,
  OR_threshold,
  bh_p_value_threshold
) {
  VALUE <- ASPECT <- LEVEL <-
    OR <- OR_UP <- OR_DOWN <- OR_FLAT <-
    BH_PVAL <- BH_P_VALUE_UP <- BH_P_VALUE_DOWN <- BH_P_VALUE_FLAT <-
    ORA_SCORE <- ORA_SCORE_UP <- ORA_SCORE_DOWN <- ORA_SCORE_FLAT <- NULL
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"ggplot2\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (OR_threshold < 1) {
    stop(
      "'OR_thtreshold' muste be bigger than 1"
    )
  }
  if (bh_p_value_threshold > 0.05) {
    stop(
      "'bh_p_value_threshold' must be smaller than 0.05"
    )
  }
  dt <- copy(ora_dt)
  if(regulation == "UP") {
    dt[, OR := OR_UP]
    dt[, BH_PVAL := BH_P_VALUE_UP]
    dt[, ORA_SCORE := ORA_SCORE_UP ]
  } else if(regulation == "DOWN") {
    dt[, OR := OR_DOWN]
    dt[, BH_PVAL := BH_P_VALUE_DOWN]
    dt[, ORA_SCORE := ORA_SCORE_DOWN ]
  } else if(regulation == "FLAT") {
    dt[, OR := OR_FLAT]
    dt[, BH_PVAL := BH_P_VALUE_FLAT]
    dt[, ORA_SCORE := ORA_SCORE_FLAT ]
  } else {
    stop("Can't find 'regulation' type")
  }
  if (category == "GO_TERMS") {
    dt <- dt[
      ASPECT == GO_aspect
    ]
    dt[, VALUE := paste0(
      "(L",
      LEVEL,
      ") ",
      VALUE
    )]
  }
  dt <- dt[
    OR > 1 &
      BH_PVAL <= 0.05
  ]
  if (any(is.infinite(dt$OR))) {
    extra_label_annotation <- " (* : infinite odds ratios are normalized)"
    dt[
      ,
      VALUE := ifelse(
        is.infinite(OR),
        paste0("* ", VALUE),
        VALUE
      )
    ]
    dt_finite <- dt[is.finite(OR)]
    if (nrow(dt_finite) > 0) {
      dt[
        ,
        OR := ifelse(
          is.infinite(OR),
          1 + max(dt_finite$OR),
          OR
        )
      ]
    } else {
      dt[, OR := 100]
    }
    dt[
      ,
      ORA_SCORE := -log10(BH_PVAL) * log2(OR)
    ]
  } else {
    extra_label_annotation <- NULL
  }
  dt <- dt[
    OR > OR_threshold &
      BH_PVAL <= bh_p_value_threshold
  ][order(-ORA_SCORE)]
  n_row_tokeep <- min(max_terms_show, nrow(dt))
  dt <- dt[1:n_row_tokeep]
  dt$VALUE <- sapply(
    dt$VALUE,
    function(i) {
      words <- strsplit(i, " ")[[1]]
      n_words <- length(words)
      if (n_words >= 5) {
        if (n_words%%2 == 0) {
          mid <- n_words/2
        } else {
          mid <- (n_words+1)/2
        }
        res <- paste0(
          paste0(words[1:mid], collapse = " "),
          "\n",
          paste0(words[(mid+1):length(words)], collapse = " ")
        )
      } else {
        res <- i
      }
      res
    }
  )
  category_label <- ifelse(
    category == "LRI",
    "Ligand-Receptor Interactions",
    ifelse(
      category == "LIGAND_COMPLEX",
      "Ligand Genes",
      ifelse(
        category == "RECEPTOR_COMPLEX",
        "Receptor Genes",
        ifelse(
          category == "ER_CELLTYPES",
          "Emitter-Receiver Cell Types",
          ifelse(
            category == "EMITTER_CELLTYPE",
            "Emitter Cell Types",
            ifelse(
              category == "RECEIVER_CELLTYPE",
              "Receiver Cell Types",
              ifelse(
                category == "GO_TERMS",
                ifelse(
                  GO_aspect == "biological_process",
                  "GO Biological Processes",
                  ifelse(
                    GO_aspect == "molecular_function",
                    "GO Molecular Functions",
                    "GO Cellular Components"
                  )
                ),
                ifelse(
                  category == "KEGG_PWS",
                  "KEGG Pathways",
                  category
                )
              )
            )
          )
        )
      )
    )
  )
  ggplot2::ggplot(
    dt,
    ggplot2::aes(
      ORA_SCORE,
      stats::reorder(VALUE, ORA_SCORE)
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = -log10(BH_PVAL),
        color = log2(OR)
      )
    ) +
    ggplot2::scale_color_gradient(low = "orange", high = "red") +
    ggplot2::scale_size_continuous(range = c(5, 8)) +
    ggplot2::xlab(paste0("ORA score ", regulation)) +
    ggplot2::ylab("") +
    ggplot2::labs(
      size = "-log10(Adj. P-Value)",
      color = "log2(Odds Ratio)",
      caption = extra_label_annotation
    ) +
    ggplot2::theme(text = ggplot2::element_text(size = 40)) +
    ggplot2::theme(legend.position = c(0.8, 0.3)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 34)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 26)) +
    ggplot2::theme(plot.title.position = "plot") +
    ggplot2::ggtitle(
      paste0(
        "Top ",
        n_row_tokeep,
        " over-represented ",
        regulation,
        "-regulated ",
        category_label
      )
    )
}

plot_ora_local(
  ora_dt = ORA_table[
    ORA_CATEGORY == "LRI" &
      Dataset == "TMS Droplet (male)" &
      Tissue == "Bladder"
  ],
  category = "LRI",
  regulation = "UP",
  max_terms_show = 20,
  GO_aspect = NULL,
  OR_threshold = 1,
  bh_p_value_threshold = 0.045
)

plot_ORA_GO_treemap_local <- function(
  GO_REDUCED_table,
  tissue_choice,
  dataset_choice,
  type_choice,
  go_aspect_choice,
  title_text,
  domain = NULL
) {
  Dataset <- Tissue <- ASPECT <- REGULATION <- 
    new_parent <- term <- parentTerm <- ids <- 
    parents <- score <- text <- NULL
  ex_data <- GO_REDUCED_table[
    Dataset == dataset_choice &
      Tissue == tissue_choice &
      ASPECT == go_aspect_choice &
      REGULATION == type_choice
  ][, c("score", "term", "parentTerm")]
  if (nrow(ex_data) == 0) return(NULL)
  ex_data[, new_parent := ifelse(
    term %in% parentTerm,
    "",
    parentTerm
  )]
  new_data <- data.table(
    labels = c(ex_data$term, ex_data[new_parent == ""]$term),
    parents = c(ex_data$parentTerm, rep("", length(ex_data[new_parent == ""]$term)))
  )
  new_data[
    ,
    ids := sapply(
      1:nrow(.SD),
      function(i) {
        if (labels[[i]] == parents[[i]]) {
          res <- paste(labels[[i]], parents[[i]], sep = " - ")
        } else {
          res <- labels[[i]]
        }
        res
      }
    )
  ]
  new_data[
    ,
    score := sapply(
      1:nrow(.SD),
      function(i) {
        if (parents[[i]] == "") {
          res <- sum(ex_data[parentTerm == labels[[i]]]$score)
        } else {
          res <- ex_data[term == labels[[i]]]$score
        }
        res
      }
    )
  ]
  new_data[
    ,
    text := gsub(" ", "\n", labels)
  ]
  m <- list(
    l = 5,
    r = 5,
    b = 5,
    t = 30,
    pad = 0
  )
  plotly::plot_ly(
    new_data,
    type = "treemap",
    opacity = 1,
    ids = ~ids,
    parents = ~parents,
    values = ~score,
    labels = ~labels,
    text = ~text,
    textposition = "middle center",
    branchvalues = "total",
    hoverinfo = "label+value",
    marker = list(
      line = list(color = "black")
    ),
    textinfo = "text",
    domain = domain
  ) %>% plotly::layout(
    title = list(
      text = title_text,
      font = list(size = 16),
      xanchor = "left",
      x = 0.0
    ),
    #uniformtext = list(
    #  minsize = 16,
    #  mode = "hide"
    #),
    margin = m
  )
}

GO_REDUCED_table_local <- readRDS("../data_scAgeCom/analysis/outputs_data/scAgeCom_shiny_data.rds")$GO_REDUCED_table

p_treemap <- plot_ORA_GO_treemap_local(
  GO_REDUCED_table = GO_REDUCED_table_local,
  tissue_choice = "Bladder",
  dataset_choice = "TMS Droplet (male)",
  type_choice = "UP",
  go_aspect_choice = "biological_process",
  title_text = paste0(
    "GO Biological Processes - ",
    "UP"
  )
)

Sys.setenv(
  "PATH" = paste(
    Sys.getenv("PATH"),
    "C:\\Users\\clagger\\AppData\\Local\\Programs\\orca",
    sep = .Platform$path.sep)
)
orca(p_treemap, "../data_scAgeCom/figures/fig4_ora_go.png", scale = 3)

## Figure 5 ####

fig5_data <- readRDS("../data_scAgeCom/analysis/outputs_data/scAgeCom_shiny_data.rds")

plot_KEYWORD_summary <- function(
  ora_keyword_summary,
  ora_keyword_template,
  category,
  keyword
) {
  ORA_CATEGORY <- VALUE <- Regulation <- i.ORA_REGULATION <-
    Dataset <- Tissue <- NULL
  dt <- copy(ora_keyword_summary)[
    ORA_CATEGORY == category
  ]
  if (!(keyword %in% dt$VALUE)) {
    stop("`keyword` not found")
  }
  dt <- dt[
    VALUE == keyword
  ]
  dt <- copy(ora_keyword_template)[
    dt,
    on = c("Tissue", "Dataset"),
    Regulation := i.ORA_REGULATION
  ]
  dt[is.na(dt)] <- "Not Detected"
  print(dt)
  p <- ggplot2::ggplot(dt) +
    ggplot2::geom_tile(
      ggplot2::aes(
        Dataset,
        Tissue,
        fill = Regulation,
        width = 0.9,
        height = 0.9
      ),
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c(
        #"No Data" = "transparent",
        "Not Over-represented" = "white",
        "Not Detected" = "gray",
        "UP" = "red"#,
        #"DOWN" = "blue",
        #"FLAT" = "green",
        #"UP:DOWN" = "yellow"
      ),
      drop = FALSE
    ) +
    ggplot2::ggtitle(
      stringr::str_trunc(
        paste0(
          "Over-representation of ",
          keyword
        ),
        70, 
        "right"
      )
    ) +
    ggplot2::scale_x_discrete(
      limits = c(
        "TMS FACS (male)",
        "TMS FACS (female)" ,
        "TMS Droplet (male)",
        "TMS Droplet (female)",
        "Calico Droplet (male)"
      ),
      labels = c(
        "TMS\nFACS\n(male)",
        "TMS\nFACS\n(female)",
        "TMS\nDroplet\n(male)",
        "TMS\nDroplet\n(female)",
        "Calico\nDroplet\n(male)"
      )
    ) +
    ggplot2::scale_y_discrete(
      limits = sort(
        unique(dt$Tissue),
        decreasing = TRUE
      )
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme(text=ggplot2::element_text(size = 32)) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 32, face = "bold")) +
    ggplot2::theme(legend.position = c(0.8, 0.8))
  p
}


plot_KEYWORD_summary(
  fig5_data$ORA_KEYWORD_SUMMARY,
  fig5_data$ORA_KEYWORD_TEMPLATE,
  "Ligand-Receptor Interaction",
  "Lgals3:Lag3"
)
