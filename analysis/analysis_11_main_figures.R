####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## Create (6) main figures for manuscript
##
####################################################
##

## Add libraries ####

library(ggplot2)
library(kableExtra)
library(webshot2)
library(htmlwidgets)
library(shiny)
library(ComplexUpset)
library(plotly)

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
  name = "",
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
    "fig_lri_upset_mouse3.pdf"
  ),
  plot = figp_lri_upset_mouse,
  width = 88,
  height = 50,
  units = "mm",
  scale = 6
)

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
fig_seurat_mammary <- RunUMAP(fig_seurat_mammary, dims=1:15, reduction = "PCA")
Idents(fig_seurat_mammary) <- fig_seurat_mammary$cell_ontology_final

fig_tsne_celltype <- DimPlot(
  fig_seurat_mammary,
  reduction = "tsne",
  pt.size = 0.5,
  cols = c("blue", "brown", "red", "grey"),
) + theme(text = element_text(size = 26))

ggsave(
  paste0(
    path_scagecom_output,
    "fig_tsne_celltype2.pdf"
  ),
  fig_tsne_celltype,
  width = 30,
  height = 25,
  units = "mm",
  scale = 6
)

Idents(fig_seurat_mammary) <- fig_seurat_mammary$age_group

fig_tsne_cond <- Seurat::DimPlot(
  fig_seurat_mammary,
  reduction = "tsne",
  pt.size = 0.5
) + theme(text = element_text(size = 26))

ggsave(
  paste0(
    path_scagecom_output,
    "fig_tsne_cond2.pdf"
  ),
  fig_tsne_cond,
  width = 30,
  height = 25,
  units = "mm",
  scale = 6
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
  "....", "Cell N"
)
rownames(fig_sc_data) <- c(
  paste(
    "Gene",
    1:(nrow(fig_sc_data) - 2)
  ),
  "....", "Gene M"
)
fig_sc_data$.... <- "...."
fig_sc_data[3, 1:4] <- "...."

fig_sc_data_gt <- gt(fig_sc_data, rownames_to_stub = TRUE) |>
  fmt_auto() |>
  cols_align(
    align = "center",
    columns = 1:5
  ) |>
  opt_table_font(
    font = list("Arial")
  ) |> 
  tab_options(
    table.font.size = 10
  ) |>
  gtsave(paste0(path_scagecom_output, "sfig_workflow_ncounts.pdf"))

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
fig_sc_aggr <- as.data.frame(fig_sc_aggr)

colnames(fig_sc_aggr) <- c("Cell Type 1 Cond1", "Cell Type 1 Cond2", "....", "Cell Type 4 Cond1", "Cell Type 4 Cond2")
rownames(fig_sc_aggr) <- c(
  paste(
    "Gene",
    1:(nrow(fig_sc_aggr) - 2)
  ),
  "....",
  "Gene M"
)
fig_sc_aggr[3, ] <- "...."
fig_sc_aggr[, 3] <- "...."


fig_sc_agg_gt <- gt(fig_sc_aggr, rownames_to_stub = TRUE) |>
  fmt_auto() |>
  cols_align(
    align = "center",
    columns = 1:6
  ) |>
  opt_table_font(
    font = list("Arial")
  ) |>
  cols_width(
    starts_with("Cell") ~ px(60)
  ) |> 
  tab_options(
    table.font.size = 10
  ) |>
  gtsave(paste0(path_scagecom_output, "sfig_workflow_agg.pdf"))

## Prepare Figure Workflow LRI (Fig.2.2) ####

fig_LRI_rd_data <- as.data.frame(
  scDiffCom::LRI_mouse$LRI_curated[c(1, 643, 2, 3573), c(2, 3, 1, 4, 5, 6)]
)
fig_LRI_rd_data[is.na(fig_LRI_rd_data)] <- ""
fig_LRI_rd_data[c(1, 3), ] <- "...."
fig_LRI_rd_data[, 3] <- "\u00A0 \u00A0"
colnames(fig_LRI_rd_data) <- c("Ligand 1", "Ligand 2", "\u00A0 \u00A0  ", "Receptor 1", "Receptor 2", "Receptor 3")

fig_LRI_rd_data_gt <- gt(fig_LRI_rd_data) |>
  fmt_auto() |>
  cols_align(
    align = "center",
    columns = 1:6
  ) |>
  opt_table_font(
    font = list("Arial")
  ) |>
  tab_options(
    table.font.size = 10
  ) |>
  gtsave(paste0(path_scagecom_output, "sfig_workflow_lri_rd_data.pdf"))


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
fig_cci_pot[c(1, 3, 5), ] <- "...."
colnames(fig_cci_pot) <- c(
  "Emitter", "Receiver", "LRI",
  "Score Cond1", "Score Cond2", "Log FC"
)
fig_cci_pot

fig_cci_pot_gt <- gt(fig_cci_pot) |>
  fmt_auto() |>
  cols_align(
    align = "center",
    columns = 1:6
  ) |>
  opt_table_font(
    font = list("Arial")
  ) |>
  tab_options(
    table.font.size = 10
  ) |>
  gtsave(paste0(path_scagecom_output, "sfig_workflow_cci_pot.pdf"))

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
  x = GetDistributions(fig_aging_example)$DISTRIBUTIONS_DE[55, 1001] + 0.35,
  y = 42,
  label = "True Difference",
  size = 32
) + xlab(
  "Score(Cond2) - Score(Cond1)"
) + ylab(
  "Frequency"
) + theme(
  text = element_text(size = 100),
  axis.text.x = element_text(size = 100),
  axis.text.y = element_blank()
)
figp_distr_de

ggsave(
  paste0(
    path_scagecom_output,
    "fig_distr_de2.pdf"
  ),
  plot = figp_distr_de,
  width = 1580,
  height = 1400,
  units = "px",
  scale = 6
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
  y = 35,
  yend = 0.5,
  size = 3,
  arrow = arrow(length = unit(0.5, "cm"))
) + annotate(
  "text",
  x = 0.33 + GetDistributions(fig_aging_example)$DISTRIBUTIONS_YOUNG[55, 1001],
  y = 38,
  label = "True Score 1",
  size = 32
) + xlab(
  "Score(Cond1)"
) + ylab(
  "Frequency"
) + theme(
  text = element_text(size = 100),
  axis.text.x = element_text(size = 100),
  axis.text.y = element_blank()
) + xlim(-0.5, 2.5)
figp_distr_cond1

ggsave(
  paste0(
    path_scagecom_output,
    "fig_distr_cond1_b.pdf"
  ),
  plot = figp_distr_cond1,
  width = 1580,
  height = 1400,
  units = "px",
  scale = 6
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
  size = 32
) + xlab(
  "Score(Cond2)"
) + ylab(
  "Frequency"
) + theme(
  text = element_text(size = 100),
  axis.text.x = element_text(size = 100),
  axis.text.y = element_blank()
) + xlim(-0.5, 3)
figp_distr_cond2

ggsave(
  paste0(
    path_scagecom_output,
    "fig_distr_cond2_b.pdf"
  ),
  plot = figp_distr_cond2,
  width = 1580,
  height = 1400,
  units = "px",
  scale = 6
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

fig_cci_table_final[c(1, 3, 5), ] <- "...."
colnames(fig_cci_table_final) <- c(
  "Emitter", "Receiver", "LRI", "Log FC", "Adj. p-value", "Regulation"
)

fig_cci_table_final

fig_cci_table_final_dt <- gt(fig_cci_table_final) |>
  fmt_auto() |>
  cols_align(
    align = "center",
    columns = 1:6
  ) |>
  opt_table_font(
    font = list("Arial")
  ) |>
  tab_options(
    table.font.size = 10
  ) |>
  gtsave(paste0(path_scagecom_output, "sfig_workflow_final_dt.pdf"))

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
  paths_scd_results[grepl("diffage", paths_scd_results)],
  function(path) {
    fun_process_md(path)
  }
)
names(mds_processed) <- scd_dataset_names_diffage

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
  text = element_text(size = 32),
  axis.text = element_text(size = 32, colour = "black")
)
fig_datasets_summary
ggsave(
  paste0(
    path_scagecom_output,
    "fig_datasets_summary3.pdf"
  ),
  plot = fig_datasets_summary,
  width = 2100,
  height = 1400,
  units = "px",
  scale = 3.2
)

## Prepare Figure Tissue Specific volcano (Fig. 4.a) ####

fun_plot_volcano_cci <- function(
  cci_table,
  title_volcano
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
    t = 90,
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
      c("red", "blue", "black", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    ),
    marker = list(size = 8)
  ) %>% plotly::layout(
    title = list(
      text = title_volcano,
      font = list(size = 30),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = "Log2(FC)",
        font = list(size = 30),
        standoff = 30
      ),
      tickfont = list(
        size = 24
      )
    ),
    yaxis = list(
      title = list(
        text = "-Log10(Adj. p-value)",
        font = list(size = 30),
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
      font = list(size = 24)
    ),
    margin = m
  )
}

fig4a_volcano <- fun_plot_volcano_cci(
  dt_cci_full[
    dataset == "TMS Droplet (male)" &
      tissue == "Bladder"
  ],
  "Bladder - TMS Droplet (male)"
)

save_image(
  fig4a_volcano, 
  paste0(
   path_scagecom_output,
   "fig4a.svg"
  ),
  scale = 1
)

## Prepare Figure Tissue Specific lrfc (Fig. 4.b) ####

fun_plot_lrfc_cci <- function(
  cci_table,
  title_lrfc,
  x_axis_label = "Ligand Log2(FC)"
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
    t = 90,
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
      c("red", "blue", "black", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    ),
    marker = list(size = 8)
  )  %>% plotly::layout(
    title = list(
      text = title_lrfc,
      font = list(size = 30),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = x_axis_label,
        font = list(size = 30),
        standoff = 30
      ),
      tickfont = list(
        size = 24
      )
    ),
    yaxis = list(
      title = list(
        text = "Receptor Log2(FC)",
        font = list(size = 30),
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
      y = 1.1,
      font = list(size = 24)
    ),
    margin = m
  )
}

fig4b_lrfc <- fun_plot_lrfc_cci(
  dt_cci_full[
    dataset == "TMS Droplet (male)" &
      tissue == "Bladder"
  ],
  "Bladder - TMS Droplet (male)"
)

save_image(
  fig4b_lrfc, 
  paste0(
    path_scagecom_output,
    "fig4b.svg"
  ),
  scale = 1
)

## Prepare Figure Tissue Specific visnetwork (Fig. 4.c) ####

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

## Prepare Figure Tissue Specific LRI ORA (Fig. 4.d) ####

fig4d_data <- readRDS(
  "../data_scAgeCom/output/scDiffCom_14_05_2022/scdiffcom_droplet_diffage_male/scdiffcom_Bladder.rds"
)
fig4d_fig <- PlotORA(
  fig4d_data, "LRI", "UP", 10
) +  theme(text = element_text(size = 20))
ggsave(
  paste0(
    path_scagecom_output,
    "fig4d.pdf"
  ),
  plot = fig4d_fig,
  width = 1580,
  height = 1400,
  units = "px",
  scale = 1.4
)


## Prepare Figure Cross Tissue T cell summary (Fig 4e and 4f) ####

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
        "FLAT" = "black"#,
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

fig4e_tcell <- plot_KEYWORD_summary(
  ora_keyword_summary = shiny_list_full$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_list_full$ORA_KEYWORD_TEMPLATE,
  category = "GO Term",
  keyword = "T cell differentiation"
)
#2000x1400
ggsave(
  paste0(
    path_scagecom_output,
    "fig4e_tcell.pdf"
  ),
  plot = fig4e_tcell,
  width = 2000,
  height = 1400,
  units = "px",
  scale = 2.7
)

fig4f_b2m <- plot_KEYWORD_summary(
  ora_keyword_summary = shiny_list_full$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_list_full$ORA_KEYWORD_TEMPLATE,
  category = "Ligand-Receptor Interaction",
  keyword = "B2m:Cd3g"
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig4f.pdf"
  ),
  plot = fig4f_b2m,
  width = 2000,
  height = 1400,
  units = "px",
  scale = 2.7
)

plot_KEYWORD_summary(
  ora_keyword_summary = shiny_list_full$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_list_full$ORA_KEYWORD_TEMPLATE,
  category = "Ligand-Receptor Interaction",
  keyword = "Gpi1:Amfr"
)

## Prepare Figure ORA hUVEC family (Fig 5a) ####

fig_ora_huvec_fam <- ggplot(
  dt_val_huvec_fam,
  aes(
    x = OR,
    y = reorder(category, OR),
    size = -log10(
      BH + min(dt_val_huvec_fam[BH != 0]$BH)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
      dt_val_huvec_fam$BH + min(dt_val_huvec_fam[BH != 0]$BH)
    )))[1:4],
  labels = round(fivenum(-log10(
    dt_val_huvec_fam$BH + min(dt_val_huvec_fam[BH != 0]$BH)
  )))[1:4],
  limit = c(0, 500)
) + xlab(
  "Odds Ratio"
) + ylab(
  "scAgeCom emitter cell-type family"
#) + theme_minimal(
) + ggtitle(
  "Association with the hUVEC secretome"
)  + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 18)
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_huvec_fam2.pdf"
  ),
  plot = fig_ora_huvec_fam,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)

## Prepare Figure ORA mBMMM family (Fig 5b) ####

fig_ora_mbmm_fam <- ggplot(
  dt_val_macro_fam,
  aes(
    x = OR,
    y = reorder(category, OR),
    size = -log10(
      BH + min(dt_val_macro_fam[BH != 0]$BH)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
      dt_val_macro_fam$BH + min(dt_val_macro_fam[BH != 0]$BH)
    )))[c(1, 2, 3, 5)],
  labels = round(fivenum(-log10(
    dt_val_macro_fam$BH + min(dt_val_macro_fam[BH != 0]$BH)
  )))[c(1, 2, 3, 5)],
  limit = c(0, 500)
) + xlab(
  "Odds Ratio"
) + ylab(
  "scAgeCom emitter cell-type family"
#) + theme_minimal(
) + ggtitle(
  "Association with the mBMM secretome"
) + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 18)
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_mbmm_fam2.pdf"
  ),
  plot = fig_ora_mbmm_fam,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)

## Prepare Figure ORA neuron family (Fig 5c) ####

fig_ora_neuron_fam <- ggplot(
  dt_val_neuron_fam,
  aes(
    x = OR,
    y = reorder(category, OR),
    size = -log10(
      BH + min(dt_val_neuron_fam[BH != 0]$BH)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
    dt_val_neuron_fam$BH + min(dt_val_neuron_fam[BH != 0]$BH)
  ))),
  labels = round(fivenum(-log10(
    dt_val_neuron_fam$BH + min(dt_val_neuron_fam[BH != 0]$BH)
  ))),
  limit = c(0, 500)
) + xlab(
  "Odds Ratio"
) + ylab(
  "scAgeCom emitter cell-type family"
  #) + theme_minimal(
) + ggtitle(
  "Association with the mNeuron secretome"
)  + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 18)
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_neuron_fam2.pdf"
  ),
  plot = fig_ora_neuron_fam,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)


## Prepare Figure ORA mMSC-AT family (Fig 5d) ####

fig_ora_mmscat_fam <- ggplot(
  dt_val_mscat_fam,
  aes(
    x = OR,
    y = reorder(category, OR),
    size = -log10(
      BH + min(dt_val_mscat_fam[BH != 0]$BH)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
      dt_val_mscat_fam$BH + min(dt_val_mscat_fam[BH != 0]$BH)
    )))[c(1, 2, 3, 5)],
  labels = round(fivenum(-log10(
    dt_val_mscat_fam$BH + min(dt_val_mscat_fam[BH != 0]$BH)
  )))[c(1, 2, 3, 5)],
  limit = c(0, 500)
) + xlab(
  "Odds Ratio"
) + ylab(
  "scAgeCom emitter cell-type family"
#) + theme_minimal(
) + ggtitle(
  "Association with the mMSC-AT secretome"
) + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 18)
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_mmscat_fam2.pdf"
  ),
  plot = fig_ora_mmscat_fam,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)

## Prepare Figure ORA cardio family (Fig 5e) ####

fig_ora_cardio_fam <- ggplot(
  dt_val_cardio_fam,
  aes(
    x = OR,
    y = reorder(category, OR),
    size = -log10(
      BH + min(dt_val_cardio_fam[BH != 0]$BH)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
      dt_val_cardio_fam$BH + min(dt_val_cardio_fam[BH != 0]$BH)
    ))),
  labels = round(fivenum(-log10(
    dt_val_cardio_fam$BH + min(dt_val_cardio_fam[BH != 0]$BH)
  ))),
  limit = c(0, 500)
) + xlab(
  "Odds Ratio"
) + ylab(
  "scAgeCom emitter cell-type family"
#) + theme_minimal(
) + ggtitle(
  "Association with the rCM secretome"
)  + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 18)
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_cardio_fam2.pdf"
  ),
  plot = fig_ora_cardio_fam,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)

## Prepare Figure ORA hpde family (Fig 5f) ####

fig_ora_hpde_fam <- ggplot(
  dt_val_pancreas_fam,
  aes(
    x = OR,
    y = reorder(category, OR),
    size = -log10(
      BH + min(dt_val_pancreas_fam[BH != 0]$BH)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
      dt_val_pancreas_fam$BH + min(dt_val_pancreas_fam[BH != 0]$BH)
    ))),
  labels = round(fivenum(-log10(
    dt_val_pancreas_fam$BH + min(dt_val_pancreas_fam[BH != 0]$BH)
  ))),
  limit = c(0, 500)
) + xlab(
  "Odds Ratio"
) + ylab(
  "scAgeCom emitter cell-type family"
#) + theme_minimal(
) + ggtitle(
  "Association with the hPDE secretome"
) + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 18)
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_hpde_fam2.pdf"
  ),
  plot = fig_ora_hpde_fam,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)

## Prepare Figure 5 at once ####

fig5_complete <- cowplot::plot_grid(
  plotlist = list(
    fig_ora_huvec_fam,
    fig_ora_mbmm_fam,
    fig_ora_neuron_fam,
    fig_ora_mmscat_fam,
    fig_ora_cardio_fam,
    fig_ora_hpde_fam
  ),
  ncol = 2
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_complete_fam.png"
  ),
  plot = fig5_complete,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 3.2
)

## Prepare Figure 6 ####

regulation_distr_long <- copy(shiny_tissue_counts_summary)
setnames(
  regulation_distr_long,
  old = c("Flat CCIs", "Down CCIs", "UP CCIs", "NSC CCIs"),
  new= c("FLAT", "DOWN", "UP", "NSC")
)

regulation_distr_long <- melt.data.table(
  regulation_distr_long,
  id.vars = c("Tissue", "Dataset"),
  measure.vars = c("FLAT", "DOWN", "UP", "NSC"),
  variable.name = "REGULATION",
  value.name = "N"
)

regulation_distr_long <- regulation_distr_long[
  ,
  {
  totwt = sum(N)
  .SD[,.(pct=N/totwt), by = REGULATION ]
  },
  by = c("Tissue", "Dataset")
]

regulation_distr_long[, Dataset := factor(Dataset, levels = c(
  "TMS FACS (male)", "TMS FACS (female)", "TMS Droplet (male)", "TMS Droplet (female)", "Calico Droplet (male)"
  ))
]

regulation_distr_long[, REGULATION := factor(REGULATION, levels = c(
  "UP", "DOWN", "FLAT", "NSC"
))
]

fig6_complete <- ggplot(
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
  values = c("UP" = "red", "DOWN" = "blue", "FLAT" = "black", "NSC" = "grey")
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
ggsave(
  paste0(
    path_scagecom_output,
    "fig_6_complete.pdf"
  ),
  plot = fig6_complete,
  width = 3000,
  height = 1800,
  units = "px",
  scale = 3.2
)

fwrite(
  regulation_distr_long,
  paste0(
    path_scagecom_output,
    "dt_source_data_regulation_distr_fig6.csv"
  )
)

