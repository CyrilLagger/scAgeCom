####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - September 2020
##
## Provide some statistacal information about tissue,
## cell-types, genes, LR-pair and sex.
##
####################################################
##

## Libraries ####

library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)

## Analysis data path ####
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Read metadata from the seurat files ####

seurat_md <- readRDS(paste0(dir_data_analysis, "analysis_1_data_seurat_md.rds"))

## Add and modify metadata ####

# add age groups
anyNA(seurat_md$tms_facs$age)
unique(seurat_md$tms_facs$age)
seurat_md$tms_facs$age_group <- ifelse(seurat_md$tms_facs$age %in% c('1m', '3m'), 'YOUNG', 'OLD')
anyNA(seurat_md$tms_droplet$age)
unique(seurat_md$tms_droplet$age)
seurat_md$tms_droplet$age_group <- ifelse(seurat_md$tms_droplet$age %in% c('1m', '3m'), 'YOUNG', 'OLD')
seurat_md$calico_kidney$age_group <- ifelse(seurat_md$calico_kidney$age == "young", 'YOUNG', 'OLD')
seurat_md$calico_lung$age_group <- ifelse(seurat_md$calico_lung$age == "young", 'YOUNG', 'OLD')
seurat_md$calico_spleen$age_group <- ifelse(seurat_md$calico_spleen$age == "young", 'YOUNG', 'OLD')

# add tissue-cell types
seurat_md[1:2] <- lapply(seurat_md[1:2], function(x) {
  x$tissue_cell_type <- paste(x$tissue, x$cell_ontology_class, sep = "_")
  return(x)
})
seurat_md$calico_kidney$tissue <- "Kidney"
seurat_md$calico_lung$tissue <- "Lung"
seurat_md$calico_spleen$tissue <- "Spleen"
seurat_md[3:5] <- lapply(seurat_md[3:5], function(x) {
  x$tissue_cell_type <- paste(x$tissue, x$cell_type, sep = "_")
  return(x)
})

# change some categories to character
seurat_md <- lapply(seurat_md, function(x) {
  x$tissue <- as.character(x$tissue)
  x$age_group <- as.character(x$age_group)
  x$tissue_cell_type <- as.character(x$tissue_cell_type)
  return(x)
})
seurat_md$tms_facs$sex <- as.character(seurat_md$tms_facs$sex)
seurat_md$tms_droplet$sex <- as.character(seurat_md$tms_droplet$sex)

## Create tables of age vs sex for each cell type ####

tables_age_vs_sex <- lapply(
  seurat_md[1:2],
  function(md) {
    sapply(
      unique(md$tissue_cell_type),
      function(tct) {
        temp <- md[md$tissue_cell_type == tct,]
        table(temp$age_group, temp$sex)
      },
      simplify = FALSE,
      USE.NAMES = TRUE)
  }
)

odds_age_vs_sex <- lapply(
  seurat_md[1:2],
  function(md) {
    sapply(
      unique(md$tissue_cell_type),
      function(tct) {
        temp <- md[md$tissue_cell_type == tct,]
        tb <- table(temp$age_group, temp$sex) 
        if(ncol(tb) == 2 & nrow(tb) == 2) {
          res <- fisher.test(tb + 1)$estimate
        } else {
          res <- NA
        }
        return(res)
      },
      simplify = TRUE)
  }
)

tables_age_vs_sex$tms_facs$`Diaphragm_skeletal muscle satellite cell` + 1
odds_age_vs_sex$tms_facs[["Diaphragm_skeletal muscle satellite cell.odds ratio"]]

hist(log10(odds_age_vs_sex$tms_facs), breaks = 20)

## Filter tissues and cell-types with less than a given amount of cells ####

#function for filtering
filter_by_cells <- function(
  md,
  min_cells,
  group_1,
  group_2,
  all_or_any = "all"
) {
  if(all_or_any == "all") {
    res <- apply(
      table(md[[group_1]], md[[group_2]]) >= min_cells,
      MARGIN = 1,
      FUN = all
    )
  } else {
    res <- apply(
      table(md[[group_1]], md[[group_2]]) >= min_cells,
      MARGIN = 1,
      FUN = any
    )
  }
  names(res[res])
}

#tissues that remain after filtering
tissues_more_5cells <- lapply(
  seurat_md,
  filter_by_cells,
  min_cells = 5,
  group_1 = "tissue",
  group_2 = "age_group"
)
tissues_more_10cells <- lapply(
  seurat_md,
  filter_by_cells,
  min_cells = 10,
  group_1 = "tissue",
  group_2 = "age_group"
)

tissues_female_more_10cells <- lapply(
  seurat_md[1:2],
  function(md) {
    filter_by_cells(md[md$sex == "female",], 10, "tissue", "age_group")
  }
)
tissues_male_more_10cells <- lapply(
  seurat_md[1:2],
  function(md) {
    filter_by_cells(md[md$sex == "male",], 10, "tissue", "age_group")
  }
)

tissues_sex_more_10cells <- lapply(
  seurat_md[1:2],
  filter_by_cells,
  min_cells = 10,
  group_1 = "tissue",
  group_2 = "sex"
)

#cell types that remains after filtering
tct_more_5cells <- lapply(
  seurat_md,
  filter_by_cells,
  min_cells = 5,
  group_1 = "tissue_cell_type",
  group_2 = "age_group",
  all_or_any = "any"
)
tct_more_10cells <- lapply(
  seurat_md,
  filter_by_cells,
  min_cells = 10,
  group_1 = "tissue_cell_type",
  group_2 = "age_group",
  all_or_any = "any"
)
tct_more_10cells_all <- lapply(
  seurat_md,
  filter_by_cells,
  min_cells = 10,
  group_1 = "tissue_cell_type",
  group_2 = "age_group",
  all_or_any = "all"
)

tct_female_more_10cells <- lapply(
  seurat_md[1:2],
  function(md) {
    filter_by_cells(md[md$sex == "female",], 10, "tissue_cell_type", "age_group")
  }
)
tct_male_more_10cells <- lapply(
  seurat_md[1:2],
  function(md) {
    filter_by_cells(md[md$sex == "male",], 10, "tissue_cell_type", "age_group")
  }
)
tct_sex_more_10cells <- lapply(
  seurat_md[1:2],
  filter_by_cells,
  min_cells = 10,
  group_1 = "tissue_cell_type",
  group_2 = "sex"
)

## Create a summary data.frame for presentation ####

create_summary_df <- function(
  tissue_distr,
  tct_distr,
  is_sex
) {
  if(!is_sex) {
    number_tct <- sapply(tct_distr, length)
    df <- data.frame(
    name = c("Tabula Muris Senis - FACS", "Tabula Muris Senis - Droplet", "Calico"),
    tissue_number = c(sapply(tissue_distr, length)[1:2], 3),
    cell_type_number = c(number_tct[1:2], sum(number_tct[3:5]))
  )
  } else {
    df <- data.frame(
      name = c("Tabula Muris Senis - FACS", "Tabula Muris Senis - Droplet"),
      tissue_number = sapply(tissue_distr, length),
      cell_type_number = sapply(tct_distr, length)
    )
  }
  df <- transpose(df, keep.names = "")
  colnames(df) <- df[1,]
  df <- df[-c(1),]
  df$name <- c(
    "Number of tissues",
    "Number of cell types"
  )
  colnames(df)[[1]] <- ""
  g_df <- tableGrob(df, rows = NULL)
  return(g_df)
}


plot_seurat_md_number_tct_mixed <- create_summary_df(tissues_more_10cells, tct_more_10cells, FALSE)
grid.newpage()
grid.draw(plot_seurat_md_number_tct_mixed)
ggsave(filename = paste0(dir_data_analysis, "analysis_1_plot_seurat_md_number_tct_mixed.png"),
       plot = plot_seurat_md_number_tct_mixed, scale = 1.5)

plot_seurat_md_number_tct_female <- create_summary_df(tissues_female_more_10cells, tct_female_more_10cells, TRUE)
grid.newpage()
grid.draw(plot_seurat_md_number_tct_female)
ggsave(filename = paste0(dir_data_analysis, "analysis_1_plot_seurat_md_number_tct_female.png"),
       plot = plot_seurat_md_number_tct_female, scale = 1.5)

plot_seurat_md_number_tct_male <- create_summary_df(tissues_male_more_10cells, tct_male_more_10cells, TRUE)
grid.newpage()
grid.draw(plot_seurat_md_number_tct_male)
ggsave(filename = paste0(dir_data_analysis, "analysis_1_plot_seurat_md_number_tct_male.png"),
       plot = plot_seurat_md_number_tct_male, scale = 1.5)

plot_seurat_md_number_tct_sex <- create_summary_df(tissues_sex_more_10cells, tct_sex_more_10cells, TRUE)
grid.newpage()
grid.draw(plot_seurat_md_number_tct_sex)
ggsave(filename = paste0(dir_data_analysis, "analysis_1_plot_seurat_md_number_tct_sex.png"),
       plot = plot_seurat_md_number_tct_sex, scale = 1.5)


## Order tissue alphabetically ####
seurat_md <- lapply(seurat_md, function(md) {
  md$tissue <- factor(md$tissue, levels = sort(unique(md$tissue), decreasing = TRUE))
  return(md)
})

## Plot number of cell-types ####

plot_cells_per_tissue <- cowplot::plot_grid(
  plotlist = lapply(seurat_md, function(md) {
    ggplot(md, aes(x = tissue, fill = age_group)) + geom_bar(position = "dodge") +
      theme(text=element_text(size=28)) +
      scale_y_log10() +
      xlab("Tissue") +
      ylab("Number of cells") +
      labs(fill = "Age") +
      coord_flip()
  }),
  ncol = 2,
  labels = c("TMS FACS", "TMS Droplet", "Calico"),
  align = "v",
  rel_heights = c(2, 1.2)#,
  #label_x = 0, label_y = 0,
  #hjust = -0.5, vjust = -0.5
)
plot_cells_per_tissue
ggsave(paste0(dir_data_analysis, "analysis_1_plot_cells_per_tissue.png"),
       plot = plot_cells_per_tissue, scale = 2.5)

# number of cell-types per tissues
celltype_distr <- lapply(seurat_md, function(md) {
  setDT(md)
  unique(md[,c("tissue", "tissue_cell_type")])
})

plot_celltypes_per_tissue <- cowplot::plot_grid(
  plotlist = lapply(celltype_distr, function(md) {
    ggplot(md, aes(x = tissue)) + geom_bar() +
      geom_text(stat='count', aes(label=..count..), hjust = -0.1) +
      theme(text=element_text(size=28)) +
      xlab("Tissue") +
      ylab("Number of cell-types") +
      coord_flip()
  }),
  ncol = 2,
  labels = c("TMS FACS", "TMS Droplet", "Calico"),
  align = "v",
  rel_heights = c(2, 1.2)#,
  #label_x = 0, label_y = 0,
  #hjust = -0.5, vjust = -0.5
)
plot_celltypes_per_tissue
ggsave(paste0(dir_data_analysis, "analysis_1_plot_celltypes_per_tissue.png"),
       plot = plot_celltypes_per_tissue, scale = 2.4)


## Comparison of common tissues ####

seurat_md_2 <- seurat_md[1:2]
seurat_md_2$calico <- rbindlist(
  l = list(seurat_md$calico_kidney, seurat_md$calico_lung, seurat_md$calico_spleen)
)

seurat_common_md <- lapply(seurat_md_2, function(x) {
  x <- x[x$tissue %in% c("Kidney", "Lung", "Spleen"),]
  x$tissue <- factor(x$tissue, levels = sort(unique(as.character(x$tissue))))
  return(x)
})

plot_cells_per_celltype <- cowplot::plot_grid(
  plotlist = lapply(seurat_common_md, function(md) {
    ggplot(md, aes(x = tissue_cell_type, fill = age_group)) + geom_bar(position = "dodge") +
      theme(text=element_text(size=16),
            axis.text.y=element_blank()) +
      xlab("Cell type (Kidney, Lung, Spleen)") +
      scale_y_log10() +
      ylab("Number of cells") +
      labs(fill = "Age") +
      coord_flip()
  })
  ,
  ncol = 2,
  labels = c("TMS FACS", "TMS Droplet", "Calico"),
  align = "v",
  rel_heights = c(1.3, 1)#,
  #label_x = 0, label_y = 0,
  #hjust = -0.5, vjust = -0.5
)
plot_cells_per_celltype
ggsave(paste0(dir_data_analysis, "analysis_1_plot_cells_per_celltypes.png"),
       plot = plot_cells_per_celltype, scale = 2.4)


kidney_celltypes <- list(
  tms_facs = unique(seurat_common_md$tms_facs[tissue == "Kidney",]$tissue_cell_type),
  tms_droplet = unique(seurat_common_md$tms_droplet[tissue == "Kidney",]$tissue_cell_type),
  calico = unique(seurat_common_md$calico[tissue == "Kidney",]$tissue_cell_type)
)
kidney_celltypes

lung_celltypes <- list(
  tms_facs = unique(seurat_common_md$tms_facs[tissue == "Lung",]$tissue_cell_type),
  tms_droplet = unique(seurat_common_md$tms_droplet[tissue == "Lung",]$tissue_cell_type),
  calico = unique(seurat_common_md$calico[tissue == "Lung",]$tissue_cell_type)
)
lung_celltypes

spleen_celltypes <- list(
  tms_facs = sort(unique(as.character(seurat_common_md$tms_facs[tissue == "Spleen",]$cell_ontology_class))),
  tms_droplet = sort(unique(as.character(seurat_common_md$tms_droplet[tissue == "Spleen",]$cell_ontology_class))),
  calico = sort(unique(as.character(seurat_common_md$calico[tissue == "Spleen",]$cell_type))),
  calico_sub = sort(unique(as.character(seurat_common_md$calico[tissue == "Spleen",]$subtype)))
)
spleen_celltypes

# spleen data.frame for presentation

spleen_ct_df <- data.frame(
  tms_facs = c(spleen_celltypes$tms_facs, rep("",6)),
  tms_droplet = spleen_celltypes$tms_droplet,
  calico = c(spleen_celltypes$calico, rep("", 8)),
  calico_sub = c(spleen_celltypes$calico_sub, rep("", 4))
)
colnames(spleen_ct_df) <- c("TMS FACS", "TMS Droplet", "Calico", "Calico (subtypes)")
g_spleen <- tableGrob(spleen_ct_df, rows = NULL)
grid.newpage()
grid.draw(g_spleen)

ggsave(paste0(dir_data_analysis, "analysis_1_plot_spleen_celltypes.png"),
       plot = g_spleen, scale = 1.5)

## Create a table of all tissues and cell-types accross with numbers accross the three datasets ####

cell_types_full_table <- dcast.data.table(
  rbindlist(
    l = lapply(
      1:3,
      function(i) {
        temp <- seurat_md_2[[i]][, c("tissue", "tissue_cell_type", "age_group")]
        temp[, dataset := names(seurat_md_2)[[i]]]
        temp <- temp[, .N, by = c("dataset", "tissue", "tissue_cell_type", "age_group")]
        return(temp)
      }
    ),
    use.names = TRUE
  ),
  formula = tissue + tissue_cell_type ~ dataset + age_group ,
  value.var = "N"
)

fwrite(cell_types_full_table, paste0(dir_data_analysis, "analysis_1_cell_types_full_table.csv"))
