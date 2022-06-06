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
## Create main figures for manuscript
##
####################################################
##

## Add libraries ####

library(ggplot2)
library(kableExtra)
library(webshot2)
library(htmlwidgets)
library(shiny)

## Prepare Figure LRI distribution (Fig.1a) ####

fig_dt_lri_mouse <- copy(dt_lri_mouse)
fig_dt_lri_mouse <- fig_dt_lri_mouse[
  ,
  c(
    "LIGAND_1", "LIGAND_2",
    "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
    "DATABASE"
  )
]
fig_dt_lri_mouse[
  ,
  Type := ifelse(
    !is.na(LIGAND_2) | !is.na(RECEPTOR_2),
    "Complex",
    "Simple"
  )
]
fig_dt_lri_mouse[
  ,
  c(lri_dbs) := lapply(
    lri_dbs,
    function(i) {
      ifelse(grepl(i, DATABASE), TRUE, FALSE)
    }
  )
]

figp_lri_upset_mouse <- ComplexUpset::upset(
  as.data.frame(fig_dt_lri_mouse),
  lri_dbs,
  name = "Database Groupings by Frequency",
  set_sizes = ComplexUpset::upset_set_size()
  + ylab("Database Size"),
  base_annotations = list(
    "Intersection Size" = ComplexUpset::intersection_size(
      mapping = aes(fill = Type),
      counts = TRUE,
      bar_number_threshold = 100,
      text = list(size = 8)
    ) + scale_fill_manual(
      values = c("purple", "coral")
    ) + theme(
      axis.title.y = element_text(
        margin = margin(t = 0, r = -200, b = 0, l = 0)
      ),
      legend.position = c(0.8, 0.85)
    )
  ),
  themes = ComplexUpset::upset_default_themes(
    text = element_text(size = 40),
    plot.title = element_text(size = 34)
  ),
  min_size = 40
) + ggtitle(
  "Origin of curated mouse ligand-receptor interactions"
)
figp_lri_upset_mouse
ggsave(
  paste0(
    path_scagecom_output,
    "fig_lri_upset_mouse.png"
  ),
  plot = figp_lri_upset_mouse,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 3
)
#manual save: 2100x1200

## Prepare Figure LRI functional annotation (Fig.1b) ####

#in BioRender directly

## Prepare Figure Workflow tSNE/UMAP (Fig.2.1) ####

fig_seurat_mammary <- readRDS(
  paste0(
    path_scagecom_input,
    "seurat_testing_tms_facs_mammary_gland.rds"
  )
)
fig_seurat_mammary$age_group <- ifelse(
  fig_seurat_mammary$age_group == "YOUNG",
  "Cond1",
  "Cond2"
)
fig_seurat_mammary$cell_ontology_final <- ifelse(
  fig_seurat_mammary$cell_ontology_final == "basal cell",
  "Cell Type 1",
  ifelse(
    fig_seurat_mammary$cell_ontology_final == "stromal cell",
    "Cell Type 3",
    ifelse(
      fig_seurat_mammary$cell_ontology_final == "endothelial cell",
      "Cell Type 4",
      "Cell Type 2"
    )
  )
)

fig_seurat_mammary <- RunTSNE(fig_seurat_mammary, reduction = "PCA")
Idents(fig_seurat_mammary) <- fig_seurat_mammary$cell_ontology_final

fig_tsne_celltype <- DimPlot(
  fig_seurat_mammary,
  reduction = "tsne",
  pt.size = 3,
  cols = c("grey", "blue", "brown", "red"),
) + theme(text = element_text(size = 40))

ggsave(
  paste0(
    path_scagecom_output,
    "fig_tsne_celltype.png"
  ),
  fig_tsne_celltype,
  scale = 1
)

Idents(fig_seurat_mammary) <- fig_seurat_mammary$age_group

fig_tsne_cond <- Seurat::DimPlot(
  fig_seurat_mammary,
  reduction = "tsne",
  pt.size = 1
) + theme(text = element_text(size = 40))

ggsave(
  paste0(
    path_scagecom_output,
    "fig_tsne_cond.png"
  ),
  fig_tsne_cond,
  scale = 1
)

## Prepare Figure Workflow pre-processing (Fig.2.1) ####

fig_sc_data <- expm1(
  scDiffCom::seurat_sample_tms_liver$RNA@data[1:4, 1:4]
)
fig_sc_data <- as.data.frame(
  as.matrix(
    fig_sc_data
  )
)
fig_sc_data <- round(
  fig_sc_data,
  1
)
fig_sc_data <- data.frame(
  lapply(
    fig_sc_data,
    as.character
  ),
  stringsAsFactors = FALSE
)
colnames(fig_sc_data) <- c(
  paste(
    "Cell",
    1:(ncol(fig_sc_data) - 2)
  ),
  "...", "Cell N"
)
rownames(fig_sc_data) <- c(
  paste(
    "Gene",
    1:(nrow(fig_sc_data) - 2)
  ),
  "...", "Gene M"
)
fig_sc_data$... <- "..."
fig_sc_data[3, 1:4] <- "..."

fig_aggr_group <- paste(
  seurat_sample_tms_liver$cell_type,
  seurat_sample_tms_liver$age_group,
  sep = "_"
)

fig_sc_aggr <- t(
  DelayedArray::rowsum(
    Matrix::t(scDiffCom::seurat_sample_tms_liver$RNA@data),
    group = fig_aggr_group
  ) / as.vector(table(fig_aggr_group))
)
fig_sc_aggr <- fig_sc_aggr[1:4, c(1, 2, 3, 7, 8)]
fig_sc_aggr <- round(
  fig_sc_aggr,
  1
)
colnames(fig_sc_aggr) <- c("Cond1", "Cond2", "...", "Cond1", "Cond2")
rownames(fig_sc_aggr) <- c(
  paste(
    "Gene",
    1:(nrow(fig_sc_aggr) - 2)
  ),
  "...",
  "Gene M"
)
fig_sc_aggr[3, ] <- "..."
fig_sc_aggr[, 3] <- "..."

cbind(
  fig_sc_data,
  fig_sc_aggr
) %>%
  kbl(
    caption = paste0(
      "<span style = 'font-size: 50px;",
      "font-weight: bold;'>Averaged expression</span>"
    ),
    align = rep("c", 5)
  ) %>%
  kable_styling(
    "striped",
    full_width = FALSE
  ) %>%
  add_header_above(
    c(" " = 5, "Cell Type 1" = 2, " " = 1, "Cell Type 4" = 2)
  ) %>%
  kable_styling(
    font_size = 50
  )

fig_sc_aggr %>%
  kbl(
    caption = paste0(
      "<span style = 'font-size: 70px;",
      "font-weight: bold;'>Averaged expression</span>"
    ),
    align = rep("c", 5)
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(
    c(" " = 1, "Cell Type 1" = 2, " " = 1, "Cell Type 4" = 2)
  ) %>%
  kable_styling(font_size = 60) %>%
  column_spec(1:6, bold = T) %>%
  save_kable(
    file = paste0(
      path_scagecom_output,
      "fig_sc_aggr.png"
    ),
    zoom = 2,
    vwidth = 1200
  )

fig_sc_data %>%
  kbl(
    caption = paste0(
      "<span style = 'font-size: 70px; ",
      "font-weight: bold;'>Normalized counts</span>"
    ),
    align = rep("c", 4)
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  kable_styling(font_size = 60) %>%
  column_spec(1:5, bold = T) %>%
  save_kable(
    file = paste0(
      path_scagecom_output,
      "fig_sc_data.png"
    ),
    zoom = 2,
    vwidth = 1200
  )

## Prepare Figure Workflow LRI (Fig.2.2) ####

fig_LRI_rd_data <- as.data.frame(
  scDiffCom::LRI_mouse$LRI_curated[c(1, 643, 2, 3573), c(2, 3, 1, 4, 5, 6)]
)
fig_LRI_rd_data[is.na(fig_LRI_rd_data)] <- ""
fig_LRI_rd_data[c(1, 3), ] <- "..."
fig_LRI_rd_data[, 3] <- "\u00A0 \u00A0"
colnames(fig_LRI_rd_data) <- c("L1", "L2", "\u00A0 \u00A0  ", "R1", "R2", "R3")

fig_LRI_rd_data %>%
  kbl(
    caption = paste0(
      "<span style = 'font-size: 50px; ",
      "font-weight: bold'>Ligand-Receptor Interactions</span>"
    ),
    align = rep("l", 5)
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(c("Ligand" = 2, "\u00A0" = 1, "Receptor" = 3)) %>%
  kable_styling(font_size = 50) %>%
  column_spec(1:5, bold = T) %>%
  column_spec(3, width = "5cm") %>%
  save_kable(
    file = paste0(
      path_scagecom_output,
      "fig_lri_rd_data.png"
    ),
    zoom = 3
  )

## Prepare Figure Workflow potential CCIs (Fig.2.3) ####

fig_cci_pot <- copy(
  dt_cci_full
)[dataset == "TMS FACS (female)" & tissue == "Mammary_Gland"]
fig_cci_pot[
  data.table(
    old_ct = unique(fig_cci_pot$EMITTER_CELLTYPE),
    new_ct = paste0(
      "Cell Type ",
      seq_along(unique(fig_cci_pot$EMITTER_CELLTYPE))
    )
  ),
  on = "EMITTER_CELLTYPE==old_ct",
  emitter_cell_type := new_ct
]
fig_cci_pot[
  data.table(
    old_ct = unique(fig_cci_pot$RECEIVER_CELLTYPE),
    new_ct = paste0(
      "Cell Type ",
      seq_along(unique(fig_cci_pot$RECEIVER_CELLTYPE))
    )
  ),
  on = "RECEIVER_CELLTYPE==old_ct",
  receiver_cell_type := new_ct
]
fig_cci_pot <- fig_cci_pot[
  sample(seq_len(nrow(fig_cci_pot)), 10),
  c(
    "emitter_cell_type",
    "receiver_cell_type",
    "LRI",
    "CCI_SCORE_YOUNG",
    "CCI_SCORE_OLD",
    "LOGFC"
  )
]
setorder(
  fig_cci_pot,
  emitter_cell_type,
  receiver_cell_type
)

fig_cci_pot <- fig_cci_pot[c(1, 2, 6, 7, 10)]
fig_cci_pot$LOGFC <- fig_cci_pot$LOGFC / log2(exp(1))
fig_cci_pot[, 4:6] <- round(fig_cci_pot[, 4:6], 1)
fig_cci_pot <- data.frame(
  lapply(fig_cci_pot, as.character),
  stringsAsFactors = FALSE
)
fig_cci_pot[c(1, 3, 5), ] <- "..."
colnames(fig_cci_pot) <- c(
  "Emitter", "Receiver", "LRI",
  "Score Cond1", "Score Cond2", "Log FC"
)
fig_cci_pot

fig_cci_pot %>%
  kbl(
    caption = paste0(
      "<span style = 'font-size: 70px;font-weight: bold;",
      "'>Hypothetic cell-cell interactions</span>"
    ),
    align = rep("c", 6)
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  kable_styling(font_size = 55) %>%
  column_spec(1:6, bold = T) %>%
  save_kable(
    file = paste0(
      path_scagecom_output,
      "fig_cci_pot.png"
    ),
    zoom = 2,
    vwidth = 2200
  )

## Prepare Figure Workflow permutations (Fig.2.4) ####

fig_aging_example <- run_interaction_analysis(
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

fig_distr_de <- data.frame(
  counts = GetDistributions(fig_aging_example)$DISTRIBUTIONS_DE[55, ]
)
hist(fig_distr_de$counts)

fig_distr_de$is_above <- ifelse(
  abs(fig_distr_de$counts) >=
    GetDistributions(fig_aging_example)$DISTRIBUTIONS_DE[55, 1001],
  TRUE,
  FALSE
)

figp_distr_de <- ggplot(
  fig_distr_de,
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
  x = GetDistributions(fig_aging_example)$DISTRIBUTIONS_DE[55, 1001],
  xend = GetDistributions(fig_aging_example)$DISTRIBUTIONS_DE[55, 1001],
  y = 38,
  yend = 0.5,
  size = 3,
  arrow = arrow(length = unit(0.5, "cm"))
) + annotate(
  "text",
  x = GetDistributions(fig_aging_example)$DISTRIBUTIONS_DE[55, 1001] + 0.8,
  y = 42,
  label = "True Difference",
  size = 22
) + xlab(
  "Score(Cond2) - Score(Cond1)"
) + ylab(
  "Frequency"
) + theme(
  text = element_text(size = 80),
  axis.text.x = element_text(size = 80),
  axis.text.y = element_blank()
)
figp_distr_de
#manual save 2000x1400
ggsave(
  paste0(
    path_scagecom_output,
    "fig_distr_de.png"
  ),
  figp_distr_de,
  scale = 0.7
)

fig_distr_cond1 <- data.frame(
  counts = GetDistributions(fig_aging_example)$DISTRIBUTIONS_YOUNG[55, ]
)
hist(fig_distr_cond1$counts)
fig_distr_cond1$counts[1001]

fig_distr_cond1$is_above <- ifelse(
  abs(fig_distr_cond1$counts) >=
    GetDistributions(fig_aging_example)$DISTRIBUTIONS_YOUNG[55, 1001],
  TRUE,
  FALSE
)

figp_distr_cond1 <- ggplot(
  fig_distr_cond1,
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
  x = GetDistributions(fig_aging_example)$DISTRIBUTIONS_YOUNG[55, 1001],
  xend = GetDistributions(fig_aging_example)$DISTRIBUTIONS_YOUNG[55, 1001],
  y = 22,
  yend = 0.5,
  size = 3,
  arrow = arrow(length = unit(0.5, "cm"))
) + annotate(
  "text",
  x = GetDistributions(fig_aging_example)$DISTRIBUTIONS_YOUNG[55, 1001],
  y = 25,
  label = "True Score 1",
  size = 22
) + xlab(
  "Score(Cond1)"
) + ylab(
  "Frequency"
) + theme(
  text = element_text(size = 80),
  axis.text.x = element_text(size = 80),
  axis.text.y = element_blank()
) + xlim(-0.5, 2.5)
figp_distr_cond1
#manual save 2000x1400
ggsave(
  paste0(
    path_scagecom_output,
    "fig_distr_cond1.png"
  ),
  figp_distr_cond1,
  scale = 0.7
)

fig_distr_cond2 <- data.frame(
  counts = GetDistributions(fig_aging_example)$DISTRIBUTIONS_OLD[55, ]
)
hist(fig_distr_cond2$counts)

fig_distr_cond2$is_above <- ifelse(
  abs(fig_distr_cond2$counts) >=
    GetDistributions(fig_aging_example)$DISTRIBUTIONS_OLD[55, 1001],
  TRUE,
  FALSE
)

figp_distr_cond2 <- ggplot(
  fig_distr_cond2,
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
  x = GetDistributions(fig_aging_example)$DISTRIBUTIONS_OLD[55, 1001],
  xend = GetDistributions(fig_aging_example)$DISTRIBUTIONS_OLD[55, 1001],
  y = 22,
  yend = 0.5,
  size = 3,
  arrow = arrow(length = unit(0.5, "cm"))
) + annotate(
  "text",
  x = GetDistributions(fig_aging_example)$DISTRIBUTIONS_OLD[55, 1001] - 0.3,
  y = 25,
  label = "True Score 2",
  size = 22
) + xlab(
  "Score(Cond2)"
) + ylab(
  "Frequency"
) + theme(
  text = element_text(size = 80),
  axis.text.x = element_text(size = 80),
  axis.text.y = element_blank()
) + xlim(-0.5, 3)
figp_distr_cond2

ggsave(
  paste0(
    path_scagecom_output,
    "fig_distr_cond2.png"
  ),
  figp_distr_cond2,
  scale = 0.7
)

## Prepare Figure Workflow CCI final (Fig.2.6) ####

fig_cci_table_final <- copy(dt_cci_full)[
  dataset == "TMS FACS (female)" & tissue == "Mammary_Gland"
]
fig_cci_table_final[
  data.table(
    old_ct = unique(fig_cci_table_final$EMITTER_CELLTYPE),
    new_ct = paste0(
      "Cell Type ",
      seq_along(unique(fig_cci_table_final$EMITTER_CELLTYPE))
    )
  ),
  on = "EMITTER_CELLTYPE==old_ct",
  emitter_cell_type := new_ct
]
fig_cci_table_final[
  data.table(
    old_ct = unique(fig_cci_table_final$RECEIVER_CELLTYPE),
    new_ct = paste0(
      "Cell Type ",
      seq_along(unique(fig_cci_table_final$RECEIVER_CELLTYPE))
    )
  ),
  on = "RECEIVER_CELLTYPE==old_ct",
  receiver_cell_type := new_ct
]
fig_cci_table_final <- fig_cci_table_final[
  sample(seq_len(nrow(fig_cci_table_final)), 10),
  c(
    "emitter_cell_type",
    "receiver_cell_type",
    "LRI",
    "LOGFC",
    "BH_P_VALUE_DE",
    "REGULATION"
  )
]
setorder(
  fig_cci_table_final,
  emitter_cell_type,
  receiver_cell_type
)
fig_cci_table_final$LOGFC <- fig_cci_table_final$LOGFC / log2(exp(1))
fig_cci_table_final[, 4] <- round(fig_cci_table_final[, 4], 1)
fig_cci_table_final[, 5] <- round(fig_cci_table_final[, 5], 2)
fig_cci_table_final <- data.frame(
  lapply(fig_cci_table_final, as.character), stringsAsFactors = FALSE
)

fig_cci_table_final <- fig_cci_table_final[c(1, 4, 5, 9, 10), ]

fig_cci_table_final[c(1, 3, 5), ] <- "..."
colnames(fig_cci_table_final) <- c(
  "Emitter", "Receiver", "LRI", "Log FC", "Adj. p-value", "Regulation"
)

fig_cci_table_final

fig_cci_table_final %>%
  kbl(
    caption = paste0(
      "<span style = 'font-size: 70px;font-weight: bold;",
      "'>Detected cell-cell interactions</span>"
    ),
    align = rep("c", 6),
    row.names = FALSE
  ) %>%
  kable_styling("striped", full_width = FALSE) %>%
  kable_styling(font_size = 55) %>%
  column_spec(1:6, bold = T) %>%
  save_kable(
    file = paste0(
      path_scagecom_output,
      "fig_cci_table_final.png"
    ),
    zoom = 2,
    vwidth = 2100
  )

## Prepare Figure Dataset summary (Fig.3) ####

fun_process_md <- function(
  md_path
) {
  # retrieve and load each object in the dataset
  tissues <- gsub(".*md_(.+)\\.rds.*", "\\1", list.files(md_path))
  tissues <- tissues[!grepl("scdiffcom_", tissues)]
  print(tissues)
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

mds_processed <- lapply(
  paths_scd_results,
  function(path) {
    fun_process_md(path)
  }
)
names(mds_processed) <- scd_dataset_names

mds_age <- rbindlist(
  lapply(
    mds_processed,
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
        idcol = "tissue"
      )
    }
  ),
  idcol = "dataset"
)

dt_datasets_summary <- dt_cci_full[!grepl("mixed", dataset)][
  ,
  c("dataset", "tissue", "EMITTER_CELLTYPE")
][
  ,
  uniqueN(EMITTER_CELLTYPE),
  by = c("dataset", "tissue")
][
  dcast.data.table(
    mds_age,
    dataset + tissue ~ age_group,
    value.var = "age_group"
  )[
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
  ][
    ,
    age_group2 := ifelse(
      dataset == "TMS FACS (female)" & tissue == "Mammary_Gland",
      "(3m vs 18/21m)",
      ifelse(
        dataset == "TMS Droplet (male)" & tissue %in% c("Liver", "Spleen"),
        "(3m vs 24m)",
        ifelse(
          dataset == "TMS Droplet (male)" & tissue %in% c("Lung"),
          "(3m vs 18m)",
          ""
        )
      )
    )
  ],
  on = c("dataset", "tissue"),
  age_group2 := i.age_group2
]

fig_datasets_summary <- ggplot(
  dt_datasets_summary,
  aes(
    dataset,
    tissue
  )
) + geom_tile(
  aes(
    width = 0.9,
    height = 0.9
  ),
  colour = "black",
  fill = "coral",
  alpha = 1
) + ggtitle(
  "Number of cell types per dataset (additional age information)"
) + geom_text(
  #aes(label = paste(V1, "Cell Types", age_group2)),
  aes(label = paste(V1, age_group2)),
  size = 8,
  fontface = "bold"
) + scale_x_discrete(
  limits = c(
    "TMS FACS (male)",
    "TMS FACS (female)",
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
) + scale_y_discrete(
  limits = sort(
    unique(dt_datasets_summary$tissue),
    decreasing = TRUE
  )
) + xlab(
  ""
) + ylab(
  ""
) + theme(
  text = element_text(size = 28),
  axis.text = element_text(size = 28, colour = "black")
)
fig_datasets_summary
ggsave(
  paste0(
    path_scagecom_output,
    "fig_datasets_summary.png"
  ),
  plot = fig_datasets_summary,
  width = 2100,
  height = 1400,
  units = "px",
  scale = 3
)
#manual save 2000x1400


## Prepare Figure Tissue Specific volcano (Fig. 4.a) ####

fun_plot_volcano_cci <- function(
  cci_table
) {
  dt <- copy(
    cci_table[
      ,
      c(
        "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",
        "LRI", "REGULATION", "LOG2FC", "BH_P_VALUE_DE"
      )
    ]
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
  dt$REGULATION <- factor(
    dt$REGULATION,
    levels = c("UP", "DOWN", "FLAT", "NSC")
  )
  m <- list(
    l = 10,
    r = 10,
    b = 30,
    t = 30,
    pad = 10
  )
  plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~LOG2FC,
    y = ~-log10(BH_P_VALUE_DE + 1E-4),
    text = ~paste(
      "LRI: ",
      LRI,
      "<br>Emitter:",
      EMITTER_CELLTYPE,
      "<br>Receiver:",
      RECEIVER_CELLTYPE
    ),
    color = ~REGULATION,
    colors = stats::setNames(
      c("red", "blue", "beige", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    ),
    marker = list(size = 10)
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
        font = list(size = 36),
        standoff = 30
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

fun_plot_volcano_cci(
  dt_cci_full[
    dataset == "TMS Droplet (male)" &
      tissue == "Bladder"
  ]
)

## Prepare Figure Tissue Specific scores (Fig. 4.b) ####

fun_plot_scores_cci <- function(
  cci_table
) {
  dt <- copy(
    cci_table[
      ,
      c(
        "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",
        "LRI", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"
      )
    ]
  )
  dt$REGULATION <- factor(
    dt$REGULATION,
    levels = c("UP", "DOWN", "FLAT", "NSC")
  )
  min_score <-  10 ^ (floor(
    log10(
      min(
        min(dt[CCI_SCORE_YOUNG > 0]$CCI_SCORE_YOUNG),
        min(dt[CCI_SCORE_OLD > 0]$CCI_SCORE_OLD)
      )
    )
  ))
  m <- list(
    l = 10,
    r = 10,
    b = 30,
    t = 30,
    pad = 10
  )
  plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~log10(CCI_SCORE_YOUNG + min_score),
    y = ~log10(CCI_SCORE_OLD + min_score),
    text = ~paste(
      "LRI: ",
      LRI,
      "<br>Emitter:",
      EMITTER_CELLTYPE,
      "<br>Receiver:",
      RECEIVER_CELLTYPE
    ),
    color = ~REGULATION,
    colors = stats::setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    ),
    marker = list(size = 10)
  ) %>% plotly::layout(
    title = list(
      text = "Interactive Score Plot",
      font = list(size = 20),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = "Log10(Young CCI Score)",
        font = list(size = 36),
        standoff = 30
      ),
      tickfont = list(
        size = 24
      )
    ),
    yaxis = list(
      title = list(
        text = "Log10(Old CCI Score)",
        font = list(size = 36),
        standoff = 40
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
      font = list(size = 34)
    ),
    margin = m
  )
}

fun_plot_scores_cci(
  dt_cci_full[
    dataset == "TMS Droplet (male)" &
      tissue == "Bladder"
  ]
)

## Prepare Figure Tissue Specific lrfc (Fig. 4.c) ####

fun_plot_lrfc_cci <- function(
  cci_table
) {
  dt <- copy(
    cci_table[
      ,
      c(
        "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",
        "LRI", "REGULATION", "LOG2FC_L", "LOG2FC_R"
      )
    ]
  )
  dt$REGULATION <- factor(
    dt$REGULATION,
    levels = c("UP", "DOWN", "FLAT", "NSC")
  )
  m <- list(
    l = 10,
    r = 10,
    b = 30,
    t = 30,
    pad = 10
  )
  plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~LOG2FC_L,
    y = ~LOG2FC_R,
    text = ~paste(
      "LRI: ",
      LRI,
      "<br>Emitter:",
      EMITTER_CELLTYPE,
      "<br>Receiver:",
      RECEIVER_CELLTYPE
    ),
    color = ~REGULATION,
    colors = stats::setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    ),
    marker = list(size = 10)
  )  %>% plotly::layout(
    title = list(
      text = "",
      font = list(size = 20),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = "Ligand Log2(FC)",
        font = list(size = 36),
        standoff = 30
      ),
      tickfont = list(
        size = 24
      )
    ),
    yaxis = list(
      title = list(
        text = "Receptor Log2(FC)",
        font = list(size = 36),
        standoff = 40
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
      font = list(size = 34)
    ),
    margin = m
  )
}

fun_plot_lrfc_cci(
  dt_cci_full[
    dataset == "TMS Droplet (male)" &
      tissue == "Bladder"
  ]
)

## Prepare Figure Tissue Specific visnetwork (Fig. 4.d) ####

fun_plot_ora_visnetwork <- function(
  cci_table,
  ora_table,
  tissue_choice,
  dataset_choice,
  abbr_celltype
) {
  cci_dt <- data.table::copy(cci_table)
  ora_table_er <- ora_table[
    dataset == dataset_choice &
      tissue == tissue_choice &
      ORA_CATEGORY == "ER_CELLTYPES"
  ]
  ora_table_emitter <- ora_table[
    dataset == dataset_choice &
      tissue == tissue_choice &
      ORA_CATEGORY == "EMITTER_CELLTYPE"
  ]
  ora_table_receiver <- ora_table[
    dataset == dataset_choice &
      tissue == tissue_choice &
      ORA_CATEGORY == "RECEIVER_CELLTYPE"
  ]
  cci_table_detected <- cci_dt[
    dataset == dataset_choice &
      tissue == tissue_choice
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
    ora_table_er[
      abbreviation_table,
      on = "EMITTER_CELLTYPE==ORIGINAL_CELLTYPE",
      "EMITTER_CELLTYPE" := i.ABBR_CELLTYPE]
    ora_table_er[
      abbreviation_table,
      on = "RECEIVER_CELLTYPE==ORIGINAL_CELLTYPE",
      "RECEIVER_CELLTYPE" := i.ABBR_CELLTYPE]
    ora_table_emitter[
      abbreviation_table,
      on = "VALUE==ORIGINAL_CELLTYPE",
      "VALUE" := i.ABBR_CELLTYPE]
    ora_table_receiver[
      abbreviation_table,
      on = "VALUE==ORIGINAL_CELLTYPE",
      "VALUE" := i.ABBR_CELLTYPE]
  }
  scDiffCom:::interactive_from_igraph(
    cci_table_detected = cci_table_detected,
    conds = c("YOUNG", "OLD"),
    ora_table_ER = ora_table_er,
    ora_table_EMITTER = ora_table_emitter,
    ora_table_RECEIVER = ora_table_receiver,
    ora_table_LR = ORA_table[
     dataset == dataset_choice &
        tissue == tissue_choice &
        ORA_CATEGORY == "LRI"
    ],
    network_type = "ORA_network",
    layout_type = "bipartite",
    object_name = tissue_choice
  )
}

fun_plot_ora_visnetwork(
  dt_cci_full,
  dt_ora_full,
  "Bladder",
  "TMS Droplet (male)",
  shiny_abbr_celltype
)

## Prepare Figure Tissue Specific go treemap (Fig. 4.e) ####

fun_plot_ora_go_treemap <- function(
  go_reduced_table,
  tissue_choice,
  dataset_choice,
  type_choice,
  go_aspect_choice,
  title_text,
  domain = NULL,
  min_size
) {
  ex_data <- go_reduced_table[
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
    parents = c(
      ex_data$parentTerm,
      rep("", length(ex_data[new_parent == ""]$term))
    )
  )
  new_data[
    ,
    ids := sapply(
      seq_len(nrow(.SD)),
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
      seq_len(nrow(.SD)),
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
    uniformtext = list(
      minsize = min_size,
      mode = "hide"
    ),
    margin = m
  )
}

fun_plot_ora_go_treemap(
  go_reduced_table = shiny_dt_go_reduced[
  Dataset ==  "TMS Droplet (male)" &
    Tissue == "Bladder" &
    ASPECT == "biological_process" &
    REGULATION == "UP" &
    cluster %in% c(1, 2, 3, 4)
],
  tissue_choice = "Bladder",
  dataset_choice = "TMS Droplet (male)",
  type_choice = "UP",
  go_aspect_choice = "biological_process",
  title_text = paste0(
    "Top GO Biological Processes - ",
    "UP"
  ),
  min_size = 12
)

## Prepare Figure Tissue Specific gene ORA (Fig. 4.f) ####

#TODO

## Prepare Figure Cross Tissue GO table (Fig. 5.a) ####

fun_display_keyword_counts <- function(
  ora_keyword_counts,
  category,
  regulation,
  go_aspect = NULL
) {
  dt <- ora_keyword_counts[
    ORA_CATEGORY == category &
      ORA_REGULATION == regulation
  ]
  data.table::setnames(dt, old = "VALUE", new = category)
  if (category == "GO Term") {
    temp_aspect <- ifelse(
      go_aspect == "Biological Process",
      "biological_process",
      ifelse(
        go_aspect == "Molecular Function",
        "molecular_function",
        "cellular_component"
      )
    )
    dt <- dt[ASPECT == temp_aspect]
    dt <- dt[, -c(1,2,10)]
    if (go_aspect == "Biological Process") {
      category_label <- paste0("GO ", go_aspect, "es")
    } else {
      category_label <- paste0("GO ", go_aspect, "s")
    }
  } else {
    dt <- dt[, -c(1,2,10,11)]
    if(grepl("Family", category)){
      category_label <- sub("Family", "Families", category)
    } else {
      category_label <- paste0(category, "s")
    }
  }
  DT <- DT::datatable(
    data = dt,
    filter = list(
      position ="top",
      clear = FALSE,
      plain = FALSE
    ),
    class = "display compact",
    options =list(
      pageLength = 5,
      dom = '<"top"f>rt<"bottom"lip><"clear">',
      columnDefs = list(
        list(width = '300px', targets = c(1))
      )
    ),
    caption = tags$caption(
      style = paste0(
        "caption-side: top; ",
        "text-align: center; ",
        "color: black; ",
        "font-size: 120%;"
      ),
      paste0(
        "Number of tissues in which ",
        category_label,
        " are over-represented among ",
        regulation,
        "-regulated cell-cell interactions"
      )
    )
  ) %>% DT::formatStyle(
    colnames(dt)[-1],
    `text-align` = "center"
  )
  if (category == "GO Term") {
    DT <- DT %>% DT::formatStyle(c(7), `border-right` = "solid 2px")
  }
  DT
}

fig_table_bp_up <- fun_display_keyword_counts(
  shiny_dt_ora_key_counts,
  category = "GO Term",
  regulation = "UP",
  go_aspect = "Biological Process"
)
fig_table_bp_up$width <- "650px"

fig_table_bp_up_html <- "fig_table_BP_up.html"
saveWidget(fig_table_bp_up, fig_table_bp_up_html)
webshot(
  fig_table_bp_up_html,
  paste0(
    path_scagecom_output,
    "fig_table_BP_up.png"
  ),
  zoom = 3
)




table_LRI_down <- display_KEYWORD_counts(
  shiny_results$ORA_KEYWORD_COUNTS,
  category = "Ligand-Receptor Interaction",
  regulation = "DOWN"
)
table_LRI_down$width <- "600px"

table_LRI_down_html <- "table_LRI_down.html"
saveWidget(table_LRI_down, table_LRI_down_html)
webshot(
  table_LRI_down_html,
  "../../../../../table_LRI_down.png",
  zoom = 3
)

table_L_down <- display_KEYWORD_counts(
  shiny_results$ORA_KEYWORD_COUNTS,
  category = "Receptor",
  regulation = "DOWN"
)
table_L_down$width <- "600px"

table_LRI_down_html <- "table_LRI_down.html"
saveWidget(table_LRI_down, table_LRI_down_html)
webshot(
  table_LRI_down_html,
  "../../../../../table_LRI_down.png",
  zoom = 3
)








plot_keyword_summary <- function(
  ora_keyword_summary,
  ora_keyword_template,
  category,
  keyword
) {
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
  p <- ggplot(dt) +
    geom_tile(
      aes(
        Dataset,
        Tissue,
        fill = Regulation,
        width = 0.9,
        height = 0.9
      ),
      colour = "black"
    ) +
    scale_fill_manual(
      name = NULL,
      values = c(
        "Not Over-represented" = "white",
        "Not Detected" = "gray",
        "UP" = "red"
      ),
      drop = FALSE
    ) +
    ggtitle(
      stringr::str_trunc(
        paste0(
          "Over-representation of ",
          keyword
        ),
        70,
        "right"
      )
    ) +
    scale_x_discrete(
      limits = c(
        "TMS FACS (male)",
        "TMS FACS (female)",
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
    scale_y_discrete(
      limits = sort(
        unique(dt$Tissue),
        decreasing = TRUE
      )
    ) +
    xlab("") +
    ylab("") +
    theme(text = element_text(size = 32)) +
    theme(
      axis.text = element_text(size = 32, face = "bold")
    ) +
    theme(legend.position = c(0.8, 0.8))
  p
}

plot_keyword_summary(
  shiny_dt_ora_key_summary,
  shiny_dt_ora_key_template,
  "Ligand-Receptor Interaction",
  "Lgals3:Lag3"
)

## Figure 5 ####




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
        "UP" = "red",
        "DOWN" = "blue",
        "FLAT" = "green"#,
        #"UP:DOWN" = "yellow"
      )
    ) +
    ggplot2::ggtitle(
      stringr::str_trunc(
        paste0(
          "Over-representation with age of ",
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
    #theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::theme(text=ggplot2::element_text(size = 30)) +
    ggplot2::theme(
      axis.text=ggplot2::element_text(size = 30),
      legend.position = c(0.85, 0.8)
      )# +
  #ggplot2::theme(legend.position = c(0.8, 0.8))
  p
  #plotly::ggplotly(
  #  p,
  #  source = "TCA_PLOT_KEYWORD_SUMMARY",
  #  tooltip = c("Dataset", "Tissue", "Regulation")
  #) #%>% plotly::layout(
  #  legend = list(
  #    title = list(text = "")
  #  ) 
  # )
}

plot_KEYWORD_summary(
  ora_keyword_summary = shiny_results$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_results$ORA_KEYWORD_TEMPLATE,
  category = "GO Term",
  keyword = "T cell differentiation"
)
#2000x1400

plot_KEYWORD_summary(
  ora_keyword_summary = shiny_results$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_results$ORA_KEYWORD_TEMPLATE,
  category = "Ligand-Receptor Interaction",
  keyword = "B2m:Cd3g"
)

plot_KEYWORD_summary(
  ora_keyword_summary = shiny_results$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_results$ORA_KEYWORD_TEMPLATE,
  category = "Ligand-Receptor Interaction",
  keyword = "Gpi1:Amfr"
)

## Figure 6 ####

ggplot(
  data = regulation_distr_long,
  aes(
    y = Tissue,
    x = pct,
    fill = REGULATION
  )
) + geom_bar(
  stat = "identity",
  position = "fill",
  alpha = 0.8
) + scale_fill_manual(
  "Age-Regulation",
  values = c("UP" = "red", "DOWN" = "blue", "FLAT" = "green", "NSC" = "grey")
) + facet_wrap(
  ~ Dataset,
  ncol = 5
) + ggplot2::scale_y_discrete(
    limits = sort(
      unique(regulation_distr_long$Tissue),
      decreasing = TRUE
    )
) + xlab(
  "Fraction of CCIs per regulation group"
) + theme(
  text = element_text(size = 40, face = "bold"),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 36, face = "bold"),
  axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
  panel.spacing = unit(2, "lines"),
  legend.position = c(0.905, 0.8)
)
#manual save 3000x1800

## To remove? #####

display_cci_table <- function(
  cci_table
) {
  dt <- cci_table[
    ,
    3:12
  ]
  CCI_DT <- DT::datatable(
    data = dt[, -c(9, 10)],
    class = "display compact",
    options = list(
      pageLength = 10,
      dom = '<"top"f>rt<"bottom"lip><"clear">'
    ),
    caption = tags$caption(
      style = paste0(
        "caption-side: top; text-align: center; ",
        "color:black; font-size:120% ;"
      ),
      "Table of Cell-Cell Interactions"
    ),
    rownames = rownames,
    extensions = c("Buttons")
  ) %>% DT::formatStyle(
    colnames(dt[, -c(9, 10)])[4:8],
    `text-align` = "center"
  )
}

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
  dt <- data.table::copy(ora_dt)
  if(regulation == "UP") {
    dt[, OR := OR_UP]
    dt[, BH_PVAL := BH_P_VALUE_UP]
    dt[, ORA_SCORE := ORA_SCORE_UP]
  } else if (regulation == "DOWN") {
    dt[, OR := OR_DOWN]
    dt[, BH_PVAL := BH_P_VALUE_DOWN]
    dt[, ORA_SCORE := ORA_SCORE_DOWN]
  } else if (regulation == "FLAT") {
    dt[, OR := OR_FLAT]
    dt[, BH_PVAL := BH_P_VALUE_FLAT]
    dt[, ORA_SCORE := ORA_SCORE_FLAT]
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
        if (n_words %% 2 == 0) {
          mid <- n_words / 2
        } else {
          mid <- (n_words + 1) / 2
        }
        res <- paste0(
          paste0(words[1:mid], collapse = " "),
          "\n",
          paste0(words[(mid + 1):length(words)], collapse = " ")
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
  ggplot(
    dt,
    aes(
      ORA_SCORE,
      stats::reorder(VALUE, ORA_SCORE)
    )
  ) +
    geom_point(
      aes(
        size = -log10(BH_PVAL),
        color = log2(OR)
      )
    ) +
    scale_color_gradient(low = "orange", high = "red") +
    scale_size_continuous(range = c(8, 10)) +
    xlab(paste0("ORA score ", regulation)) +
    ylab("") +
    labs(
      size = "-log10(Adj. P-Value)",
      color = "log2(Odds Ratio)",
      caption = extra_label_annotation
    ) +
    theme(text = element_text(size = 40)) +
    theme(legend.position = c(0.8, 0.3)) +
    theme(legend.title = element_text(size = 34)) +
    theme(legend.text = element_text(size = 26)) +
    theme(plot.title.position = "plot") +
    ggtitle(
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