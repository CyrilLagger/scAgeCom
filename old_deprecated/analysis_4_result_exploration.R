
## Old stuff ####

genage_mouse <- setDT(read.csv(paste0(path_output_directory, "inputs_data/genage_mouse.tsv"), sep = "\t", header = TRUE))
dt[, GENAGE := ifelse(
  LIGAND_1 %in% genage_mouse$Gene.Symbol | LIGAND_2 %in% genage_mouse$Gene.Symbol |
    RECEPTOR_1 %in% genage_mouse$Gene.Symbol | RECEPTOR_2 %in% genage_mouse$Gene.Symbol,
  "YES",
  "NO"
)
]

# DATASETS_PROCESSED <- lapply(
#   DATASETS_PROCESSED,
#   function(i) {
#     dt <- i@cci_table_detected
#     dt[, NCELLS_MIN :=
#          pmin(
#            get(paste0("NCELLS_EMITTER_", i@parameters$seurat_condition_id$cond1_name)),
#            get(paste0("NCELLS_EMITTER_", i@parameters$seurat_condition_id$cond2_name)),
#            get(paste0("NCELLS_RECEIVER_", i@parameters$seurat_condition_id$cond1_name)),
#            get(paste0("NCELLS_RECEIVER_", i@parameters$seurat_condition_id$cond2_name))
#          )
#     ]
#     i@cci_table_detected <- dt
#     return(i)
#   }
# )


lung_100k <- readRDS("../data_scAgeCom/scDiffCom_results_18_03_2021/scdiffcom_facs_male_size_factor_logFALSE_1e+05iter/scdiffcom_Lung.rds")
lung_10k <- readRDS("../data_scAgeCom/scDiffCom_results_17_03_2021/scdiffcom_facs_male_size_factor_logFALSE_10000iter/scdiffcom_Lung.rds")

BuildNetwork(
  lung_10k,
  network_type = "ORA_network",
  layout_type = "bipartite"
)

BuildNetwork(
  lung_100k,
  network_type = "ORA_network",
  layout_type = "bipartite"
)

lung_100k@cci_table_detected[, .N, by = REGULATION]
lung_10k@cci_table_detected[, .N, by = REGULATION]

intersect(
  lung_100k@cci_table_detected[REGULATION == "DOWN"]$CCI,
  lung_10k@cci_table_detected[REGULATION == "DOWN"]$CCI
)

plot(lung_10k@cci_table_raw$P_VALUE_DE, lung_100k@cci_table_raw$P_VALUE_DE)
identical(lung_10k@cci_table_raw$LOGFC, lung_100k@cci_table_raw$LOGFC)

# quantile threshold behaviour

lung_10k_3 <- FilterCCI(
  lung_10k,
  new_threshold_quantile_score = 0.3
)

lung_10k@cci_table_detected[, .N, by = REGULATION]
lung_10k_2@cci_table_detected[, .N, by = REGULATION]
lung_10k_3@cci_table_detected[, .N, by = REGULATION]

PlotORA(
  lung_10k,
  category = "LRI",
  regulation = "DOWN",
  max_terms_show = 10
)

tissues_of_interest <- tissues_of_interest[!(tissues_of_interest %in% c("Brain_Myeloid", "Brain_Non-Myeloid"))]

DATASETS_PROCESSED_CLEAN <- lapply(
  DATASETS_PROCESSED,
  function(i) {
    temp <- i$dataset
    temp@ora_combined_default <- list()
    temp@cci_table_detected <- temp@cci_table_detected[ID %in% tissues_of_interest]
    temp_ora <- temp@ora_table
    temp_ora <- lapply(
      temp_ora,
      function(j) {
        j[ID %in% tissues_of_interest]
        j
      }
    )
    temp@ora_table <- temp_ora
    temp
  }
)





DATASETS_PROCESSED_old <- readRDS("../data_scAgeCom/analysis/outputs_data/scAgeCom_results_processed_clean.rds")

tissues_of_interest <- sort(unique(unlist(lapply(
  DATASETS_PROCESSED,
  function(i) {
    unique(i@cci_table_detected$ID)
  }
))))


## Visualize general results and behaviour ####

#distribution of regulations by dataset

REGULATION_DISTRIBUTION <- dcast.data.table(
  rbindlist(
    lapply(
      DATASETS_PROCESSED,
      function(i) {
        dt <- i$dataset@cci_table_detected[, .N, by = REGULATION]
        dt[, PCT := N/sum(N)*100]
      }
    ),
    use.names = TRUE,
    idcol = "dataset"),
  formula = dataset ~ REGULATION,
  value.var = c("N", "PCT")
)

#regulation vs NCELLS_MIN
ggplot(DATASETS_PROCESSED$calico$dataset@cci_table_detected, aes(REGULATION, NCELLS_MIN)) + geom_boxplot() + scale_y_log10()
ggplot(DATASETS_PROCESSED$facs_male$dataset@cci_table_detected, aes(REGULATION, NCELLS_MIN)) + geom_boxplot() + scale_y_log10()


#regulation vs logfc_threshold

LOGFC_behaviour <- rbindlist(
  lapply(
    DATASETS_PROCESSED,
    function(i) {
      i$logfc_behaviour
    }
  ),
  idcol = "dataset"
)

ggplot(LOGFC_behaviour, aes(x = as.character(logfc_threshold), y = pct, fill = dataset)) + 
  geom_boxplot() +
  facet_wrap(~REGULATION) +
  xlab("Fold Change Threshold") +
  ylab("Percentage of CCIs") +
  theme(text=element_text(size=20)) +
  ggtitle("Summary over all tissues")

ggplot(LOGFC_behaviour[ID %in% c("Kidney", "Spleen", "Lung")], aes(x = as.character(logfc_threshold), y = pct, fill = dataset)) + 
  geom_boxplot() +
  facet_wrap(~REGULATION) +
  xlab("Fold Change Threshold") +
  ylab("Percentage of CCIs") +
  theme(text=element_text(size=20)) +
  ggtitle("Summary over Kidney, Lung and Spleen")

#regulation vs dataset/sex

GENDER_behaviour <- rbindlist(
  lapply(
    DATASETS_PROCESSED,
    function(i) {
      dt <- i$dataset@cci_table_detected[, .N, by = c("REGULATION", "ID")]
      dt[, PCT := 100*N/sum(N), by = "ID"]
    }
  ),
  idcol = "dataset"
)

ggplot(GENDER_behaviour, aes(x = as.character(REGULATION), y = PCT, fill = dataset)) + 
  geom_boxplot() +
  xlab("Regulation") +
  ylab("Percentage of CCIs") +
  theme(text=element_text(size=20)) +
  ggtitle("Summary over all tissues")


## Similarity analysis per tissue from upset plots #####

build_upsets <- function(
  category,
  regulation,
  min_size = 2
) {
  dt <- rbindlist(
    lapply(
      DATASETS_PROCESSED,
      function(i) {
        i$dataset@ora_table[[category]][
          get(paste0("OR_", regulation)) >= 1 & 
            get(paste0("BH_P_VALUE_", regulation)) <= 0.05, c("ID", "VALUE")
        ]
      }
    ),
    idcol = "dataset"
  )
  upset_list <- sapply(
    tissues_of_interest,
    function(tiss) {
      temp <- dt[ID == tiss]
      if (nrow(temp) == 0) return(NULL)
      temp <- dcast.data.table(
        temp,
        VALUE ~ dataset,
        fun.aggregate = length
      )
      if (ncol(temp) < 3) return(NULL)
      ComplexUpset::upset(
        as.data.frame(temp),
        colnames(temp)[-1],
        min_size = min_size,
        max_size = 500,
        set_sizes = FALSE,
        min_degree = 1,
        height_ratio = 0.7
      )
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )
  upset_list<- upset_list[lengths(upset_list) != 0]
  upset_list
}

upsets_GO_UP <- build_upsets("GO_TERMS", "UP")
plot_upsets_GO_UP <- cowplot::plot_grid(
  plotlist = upsets_GO_UP,
  labels = names(upsets_GO_UP)
)

upsets_GO_DOWN <- build_upsets("GO_TERMS", "DOWN")
plot_upsets_GO_DOWN <- cowplot::plot_grid(
  plotlist = upsets_GO_DOWN,
  labels = names(upsets_GO_DOWN)
)

upsets_KEGG_UP <- build_upsets("KEGG_PWS", "UP")
plot_upsets_KEGG_UP <- cowplot::plot_grid(
  plotlist = upsets_KEGG_UP,
  labels = names(upsets_KEGG_UP)
)

upsets_KEGG_DOWN <- build_upsets("KEGG_PWS", "DOWN")
plot_upsets_KEGG_DOWN <- cowplot::plot_grid(
  plotlist = upsets_KEGG_DOWN,
  labels = names(upsets_KEGG_DOWN)
)

upsets_LRI_UP <- build_upsets("LRI", "UP")
plot_upsets_LRI_UP <- cowplot::plot_grid(
  plotlist = upsets_LRI_UP,
  labels = names(upsets_LRI_UP)
)

upsets_LRI_DOWN <- build_upsets("LRI", "DOWN")
plot_upsets_LRI_DOWN <- cowplot::plot_grid(
  plotlist = upsets_LRI_DOWN,
  labels = names(upsets_LRI_DOWN)
)

upsets_CFAM_UP <- build_upsets("ER_CELLFAMILIES", "UP", 1)
plot_upsets_CFAM_UP <- cowplot::plot_grid(
  plotlist = upsets_CFAM_UP,
  labels = names(upsets_CFAM_UP)
)

upsets_CFAM_DOWN <- build_upsets("ER_CELLFAMILIES", "DOWN", 1)
plot_upsets_CFAM_DOWN <- cowplot::plot_grid(
  plotlist = upsets_CFAM_DOWN,
  labels = names(upsets_CFAM_DOWN)
)





## Networks of GO terms ####

library(ontologyPlot)
library(igraph)
library(visNetwork)

ontoGO <- ontoProc::getGeneOnto()
go_names <- ontoGO$name
go_id <- ontoGO$id

go_of_interest <- ora_keyword_counts[
  TYPE == "GO_TERMS" &
    REGULATION == "DOWN" & 
    GENDER == "either" & 
    DATASET == "either" & 
    N > 8
]$VALUE

go_of_interest <- ora_keyword_summary[
  TYPE == "GO_TERMS" &
    REGULATION == "UP" & 
    GENDER == "male" & 
    DATASET == "facs" & 
    ID == "Lung"
]$VALUE

go_of_interest <- test4[REGULATION == "DOWN" & N > 9]$VALUE

go_of_interest <- go_names[go_names %in% go_of_interest]

go_of_interest_bp <- exclude_descendants(
  ontoGO,
  names(go_names[go_names %in% c("molecular_function", "cellular_component")]),
  names(go_of_interest)
)


p <- onto_plot(
  ontoGO,
  remove_links(ontoGO, go_of_interest_bp),
  #frequencies = all_ora_counts_go[either_all > 12 & REGULATION == "UP"]$either_all,
  #fillcolor = "green",
  #fillcolor = colour_by_population_frequency,
  width = 0.75
)
p <- onto_plot(
  ontoGO,
  go_of_interest_bp,
  #frequencies = all_ora_counts_go[either_all > 12 & REGULATION == "UP"]$either_all,
  #fillcolor = "green",
  #fillcolor = colour_by_population_frequency,
  width = 0.75
)

p
g <- graph_from_adjacency_matrix(p$adjacency_matrix)
V(g)$name <- as.vector(sapply(
  V(g)$name,
  function(i) {
    go_names[names(go_names) == i]
  }
))
lay <- layout_as_tree(g)
lay[, 2] <- -lay[, 2]
visIgraph(g, layout = 'layout.norm', layoutMatrix = lay) 

ggraph(g, "tree") + geom_edge_diagonal()

visIgraph(g, layout = 'layout.norm', layoutMatrix = as.matrix(-g2$data[, 1:2]))


g4 <- plot(graph::graphAM(p$adjacency_matrix, edgemode = "directed"), "dot")

lay3 <- t(sapply(
  g4@AgNode,
  function(i) {
    c(i@center@x, i@center@y)
  }
))

visIgraph(g, layout = 'layout.norm', layoutMatrix = -lay3) 

process_dataset <- function(
  dataset_path,
  dataset_name,
  fc_behaviour
) {
  tissues <- gsub(".*scdiffcom_(.+)\\.rds.*", "\\1", list.files(dataset_path))
  dataset <- lapply(
    X = tissues,
    FUN = function(
      tiss
    ) {
      res <- readRDS(paste0(dataset_path, "/scdiffcom_", tiss, ".rds"))
    }
  )
  tissues[tissues == "BAT"] <- "Adipose_Tissue_(Brown)"
  tissues[tissues == "GAT"] <- "Adipose_Tissue_(Gonadal)"
  tissues[tissues == "MAT"] <- "Adipose_Tissue_(Mesenteric)"
  tissues[tissues == "SCAT"] <- "Adipose_Tissue_(Subcutaneous)"
  #tissues[tissues == "brain_combined"] <- "Brain"
  names(dataset) <- tissues
  dataset <- lapply(
    seq_along(dataset), 
    function (i) {
      #names(dataset[[i]]@parameters)[names(dataset[[i]]@parameters) == "LRdb_species"] <- "LRI_species"
      dataset[[i]]@parameters$object_name <- tissues[[i]]
      dataset[[i]]
    }
  )
  #create an scDiffComCombined object
  dataset_combined <- Combine_scDiffCom(
    l = dataset,
    object_name = dataset_name,
    verbose = TRUE
  )
  rm(dataset)
  #study dataset in function of fc
  if(fc_behaviour) {
    temp_seq <- seq(1.1, 2, 0.1)
    temp_list <- lapply(
      temp_seq,
      function(thr) {
        temp_obj <- FilterCCI(
          object = dataset_combined,
          new_threshold_logfc = log(thr),
          skip_ora = TRUE
        )
        get_cci_table_detected(temp_obj)[, {
          temp = .N
          .SD[, .(pct=.N/temp*100), by = REGULATION]
        }, by = ID ]
      }
    )
    names(temp_list) <- as.character(temp_seq)
    byfc_dt <- rbindlist(
      temp_list,
      use.names = TRUE,
      idcol = "logfc_threshold"
    )
    byfc_dt$logfc_threshold <- as.numeric(byfc_dt$logfc_threshold)
  } else {
    byfc_dt <- NULL
  }
  #keep object with logfc(1.5) and remove cci_table_raw
  if(dataset_combined@parameters$threshold_logfc != log(1.5)) {
    stop("Error")
  }
  dataset_combined@cci_table_raw <- list()
  gc()
  # add information
  cell_types_dt <- setDT(read.csv(paste0(path_output_directory, "inputs_data/scDiffCom_cell_types_clean.csv"), stringsAsFactors = FALSE))
  genage_mouse <- setDT(read.csv(paste0(path_output_directory, "inputs_data/genage_mouse.tsv"), sep = "\t", header = TRUE))
  dt <- copy(dataset_combined@cci_table_detected)
  dt[cell_types_dt, on = "EMITTER_CELLTYPE==Final_annotation", EMITTER_CELLFAMILY := i.Family_broad]
  dt[cell_types_dt, on = "RECEIVER_CELLTYPE==Final_annotation", RECEIVER_CELLFAMILY := i.Family_broad]
  dt[cell_types_dt, on = "EMITTER_CELLTYPE==Final_annotation", EMITTER_CELL_ABR := i.Abbreviation]
  dt[cell_types_dt, on = "RECEIVER_CELLTYPE==Final_annotation", RECEIVER_CELL_ABR := i.Abbreviation]
  dt[, ER_CELLFAMILIES := paste(EMITTER_CELLFAMILY, RECEIVER_CELLFAMILY, sep = "_")]
  dt[, GENAGE := ifelse(
    LIGAND_1 %in% genage_mouse$Gene.Symbol | LIGAND_2 %in% genage_mouse$Gene.Symbol |
      RECEPTOR_1 %in% genage_mouse$Gene.Symbol | RECEPTOR_2 %in% genage_mouse$Gene.Symbol,
    "YES",
    "NO"
  )
  ]
  dataset_combined@cci_table_detected <- dt
  dataset_combined <- RunORA(
    object = dataset_combined,
    categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS", "ER_CELLFAMILIES", "GENAGE"),
    overwrite = TRUE,
    #stringent_or_default = "default",
    #stringent_logfc_threshold = NULL,
    #global = FALSE
  )
  # dataset_combined <- RunORA(
  #   object = dataset_combined,
  #   categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS", "ID", "ER_CELLFAMILIES", "GENAGE"),
  #   overwrite = FALSE,
  #   stringent_or_default = "default",
  #   stringent_logfc_threshold = NULL,
  #   global = TRUE
  # )
  return(list(dataset = dataset_combined, logfc_behaviour = byfc_dt))
}
