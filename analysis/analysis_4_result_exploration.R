####################################################
##
## Project: scAgeCom
##
## Last update - February 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## Process scDiffCom results and check filtering
## parameters.
##
####################################################
##

## Libraries ####
library(scDiffCom)
library(data.table)
library(ggplot2)

## Specify the directory with all scDiffCom results ####
path_input_directory <- "../data_scAgeCom/scdiffcom_results_22_02_2021"
path_output_directory <- "../data_scAgeCom/analysis/"

## Retrieve all the results for each tissue and dataset ####
RESULT_PATHS <- list.dirs(path_input_directory, recursive = FALSE)
RESULT_PATHS2 <- RESULT_PATHS[c(4,6)]


RESULT_NAMES <- c("calico", "droplet_female", "droplet_male", "droplet_mixed", "droplet_sex",
                     "facs_female", "facs_male", "facs_mixed", "facs_sex")
RESULT_NAMES2 <- c("droplet_male_w30", "droplet_mixed_w30")

## Pipeline for each dataset independently (too heavy to load them all on a laptop) ####

process_dataset <- function(dataset_path, dataset_name, fc_behaviour) {
  #retrieve scDiffCom objects
  tissues <- gsub(".*scdiffcom_(.+)\\.rds.*", "\\1", list.files(dataset_path))
  dataset <- lapply(
    X = tissues,
    FUN = function(
      tiss
    ) {
      res <- readRDS(paste0(dataset_path, "/scdiffcom_", tiss, ".rds"))
    }
  )
  names(dataset) <- tissues
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
        get_cci_detected(temp_obj)[, {
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
  #keep object with logfc(1.5) and remove cci_raw
  if(dataset_combined@parameters$threshold_logfc != log(1.5)) {
    stop("Error")
  }
  dataset_combined@cci_raw <- list()
  gc()
  # add information
  cell_types_dt <- setDT(read.csv(paste0(path_output_directory, "inputs_data/scDiffCom_cell_types_clean.csv"), stringsAsFactors = FALSE))
  genage_mouse <- setDT(read.csv(paste0(path_output_directory, "inputs_data/genage_mouse.tsv"), sep = "\t", header = TRUE))
  dt <- copy(dataset_combined@cci_detected)
  dt[cell_types_dt, on = "EMITTER_CELLTYPE==Final_annotation", EMITTER_CELL_FAMILY := i.Family_broad]
  dt[cell_types_dt, on = "RECEIVER_CELLTYPE==Final_annotation", RECEIVER_CELL_FAMILY := i.Family_broad]
  dt[cell_types_dt, on = "EMITTER_CELLTYPE==Final_annotation", EMITTER_CELL_ABR := i.Abbreviation]
  dt[cell_types_dt, on = "RECEIVER_CELLTYPE==Final_annotation", RECEIVER_CELL_ABR := i.Abbreviation]
  dt[, ER_CELL_FAMILY := paste(EMITTER_CELL_FAMILY, RECEIVER_CELL_FAMILY, sep = "_")]
  dt[, GENAGE := ifelse(
    LIGAND_1 %in% genage_mouse$Gene.Symbol | LIGAND_2 %in% genage_mouse$Gene.Symbol |
      RECEPTOR_1 %in% genage_mouse$Gene.Symbol | RECEPTOR_2 %in% genage_mouse$Gene.Symbol,
    "YES",
    "NO"
  )
  ]
  dataset_combined@cci_detected <- dt
  dataset_combined <- RunORA(
    object = dataset_combined,
    categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "KEGG_PWS", "ER_CELL_FAMILY", "GENAGE"),
    overwrite = FALSE,
    stringent_or_default = "default",
    stringent_logfc_threshold = NULL,
    global = FALSE
  )
  dataset_combined <- RunORA(
    object = dataset_combined,
    categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "KEGG_PWS", "ID", "ER_CELL_FAMILY", "GENAGE"),
    overwrite = FALSE,
    stringent_or_default = "default",
    stringent_logfc_threshold = NULL,
    global = TRUE
  )
  return(list(dataset = dataset_combined, logfc_behaviour = byfc_dt))
}

DATASETS_PROCESSED <- lapply(
  seq_along(RESULT_PATHS),
  function(i) {
    process_dataset(RESULT_PATHS[[i]], RESULT_NAMES[[i]], TRUE)
  }
)
names(DATASETS_PROCESSED) <- RESULT_NAMES

DATASETS_PROCESSED2 <- lapply(
  seq_along(RESULT_PATHS2),
  function(i) {
    process_dataset(RESULT_PATHS2[[i]], RESULT_NAMES2[[i]], TRUE)
  }
)
names(DATASETS_PROCESSED2) <- RESULT_NAMES2

## Add extra information such as minimum number of cells in each CCI

DATASETS_PROCESSED2 <- lapply(
  DATASETS_PROCESSED2,
  function(i) {
    dt <- i$dataset@cci_detected
    dt[, NCELLS_MIN :=
                pmin(
                  get(paste0("NCELLS_EMITTER_", i$dataset@parameters$seurat_condition_id$cond1_name)),
                  get(paste0("NCELLS_EMITTER_", i$dataset@parameters$seurat_condition_id$cond2_name)),
                  get(paste0("NCELLS_RECEIVER_", i$dataset@parameters$seurat_condition_id$cond1_name)),
                  get(paste0("NCELLS_RECEIVER_", i$dataset@parameters$seurat_condition_id$cond2_name))
                )
              ]
    i$dataset@cci_detected <- dt
    return(i)
  }
)

## Save datasets ####
saveRDS(DATASETS_PROCESSED, "../data_scAgeCom/analysis/outputs_data/TMS_scAgeCom_processed.rds")
saveRDS(DATASETS_PROCESSED2, "../data_scAgeCom/analysis/outputs_data/TMS_scAgeCom_processed_bestORA_w30m.rds")
#DATASETS_PROCESSED <- readRDS("../data_scAgeCom/analysis/outputs_data/TMS_scAgeCom_processed.rds")

## Visualize general results and behaviour ####

#distribution of regulations by dataset

REGULATION_DISTRIBUTION <- dcast.data.table(
  rbindlist(
    lapply(
      DATASETS_PROCESSED,
      function(i) {
        dt <- i$dataset@cci_detected[, .N, by = REGULATION]
        dt[, PCT := N/sum(N)*100]
      }
    ),
    use.names = TRUE,
    idcol = "dataset"),
  formula = dataset ~ REGULATION,
  value.var = c("N", "PCT")
)

#regulation vs NCELLS_MIN
ggplot(DATASETS_PROCESSED$calico$dataset@cci_detected, aes(REGULATION, NCELLS_MIN)) + geom_boxplot() + scale_y_log10()
ggplot(DATASETS_PROCESSED$droplet_mixed$dataset@cci_detected, aes(REGULATION, NCELLS_MIN)) + geom_boxplot() + scale_y_log10()


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
      dt <- i$dataset@cci_detected[, .N, by = c("REGULATION", "ID")]
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

## LOGFC coherence between ligand-receptor ####

coherence_plot <- function(
  dataset,
  axis_lim
) {
  dt_filt <- copy(dataset@cci_detected)
  dt_filt[, LOG2FC_L := log2(
    pmin(
      get(paste0("L1_EXPRESSION_", dataset@parameters$seurat_condition_id$cond2_name)),
      get(paste0("L2_EXPRESSION_", dataset@parameters$seurat_condition_id$cond2_name)),
      na.rm = TRUE
    )
    /
      pmin(
        get(paste0("L1_EXPRESSION_", dataset@parameters$seurat_condition_id$cond1_name)),
        get(paste0("L2_EXPRESSION_", dataset@parameters$seurat_condition_id$cond1_name)),
        na.rm = TRUE
      )
  ) ]
  dt_filt[, LOG2FC_R := log2(
    pmin(
      get(paste0("R1_EXPRESSION_", dataset@parameters$seurat_condition_id$cond2_name)),
      get(paste0("R2_EXPRESSION_", dataset@parameters$seurat_condition_id$cond2_name)),
      get(paste0("R3_EXPRESSION_", dataset@parameters$seurat_condition_id$cond2_name)),
      na.rm = TRUE
    )
    /
      pmin(
        get(paste0("R1_EXPRESSION_", dataset@parameters$seurat_condition_id$cond1_name)),
        get(paste0("R2_EXPRESSION_", dataset@parameters$seurat_condition_id$cond1_name)),
        get(paste0("R3_EXPRESSION_", dataset@parameters$seurat_condition_id$cond1_name)),
        na.rm = TRUE
      )
  ) ]
  dt_filt[, LOG2FC_L := ifelse(is.infinite(LOG2FC_L) & LOG2FC_L > 0, axis_lim,
                               ifelse(is.infinite(LOG2FC_L) & LOG2FC_L < 0, -axis_lim, LOG2FC_L))]
  dt_filt[, LOG2FC_R := ifelse(is.infinite(LOG2FC_R) & LOG2FC_R > 0, axis_lim,
                               ifelse(is.infinite(LOG2FC_R) & LOG2FC_R < 0, -axis_lim, LOG2FC_R))]
  ggplot(dt_filt, aes(x = LOG2FC_L, y = LOG2FC_R, color = REGULATION)) +
    geom_point() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlim(c(-axis_lim-1, axis_lim+1)) +
    ylim(c(-axis_lim-1,axis_lim+1))
}

coherence_plot(DATASETS_PROCESSED$calico$dataset, 5)
coherence_plot(DATASETS_PROCESSED$facs_mixed$dataset, 10)
coherence_plot(DATASETS_PROCESSED$facs_sex$dataset, 10)
coherence_plot(DATASETS_PROCESSED$facs_male$dataset, 10)
coherence_plot(DATASETS_PROCESSED$facs_female$dataset, 10)
coherence_plot(DATASETS_PROCESSED$droplet_mixed$dataset, 10)
coherence_plot(DATASETS_PROCESSED$droplet_sex$dataset, 10)



## Number of cell-types per tissue

NCT_TISSUE <- rbindlist(
  lapply(
    DATASETS_PROCESSED,
    function(i) {
      i$dataset@cci_detected[, uniqueN(EMITTER_CELLTYPE), by = "ID"]
    }
  ),
  idcol = "dataset"
)

##Temporaryl rerun ORA ####

DATASETS_PROCESSED2 <- lapply(
  DATASETS_PROCESSED2,
  function(i) {
    new_object <- RunORA(
      object = i$dataset,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "KEGG_PWS", "ER_CELL_FAMILY", "GENAGE"),
      overwrite = TRUE,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      global = FALSE
    )
    new_object <- RunORA(
      object = new_object,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "KEGG_PWS", "ID", "ER_CELL_FAMILY", "GENAGE"),
      overwrite = TRUE,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      global = TRUE
    )
    i$dataset <- new_object
    return(i)
  }
)

##GET ORA results by rank #####

ORA_RANK_COMP <- rbindlist(
  lapply(
    DATASETS_PROCESSED2,
    function(i) {
      dt <- copy(i$dataset@ora_default$GO_TERMS[, c("ID", "VALUE", "ORA_SCORE_UP", "ORA_SCORE_DOWN")])
      dt[, ORA_RANK_UP := rank(-ORA_SCORE_UP, ties.method = "min"), by = "ID"]
      dt[, ORA_RANK_DOWN := rank(-ORA_SCORE_DOWN, ties.method = "min"), by = "ID"]
    }
  ),
  idcol = "dataset"
)

ORA_RANK_COMP2 <- ORA_RANK_COMP[ORA_SCORE_UP != 0 | ORA_SCORE_DOWN != 0]

test_lung_up <- dcast.data.table(
  ORA_RANK_COMP2[ID == "Lung", c("dataset", "VALUE", "ORA_RANK_UP")], 
  formula = VALUE ~ dataset,
  value.var = "ORA_RANK_UP"
)

ggplot(test_lung_up, aes(facs_mixed, facs_male)) + geom_point() +
  geom_smooth(method = "lm")

ggplot(test_lung_up, aes(droplet_female, droplet_mixed)) + geom_point() +
  geom_smooth(method = "lm")

##########
PlotORA(
  DATASETS_PROCESSED$droplet_female$dataset,
  subID = "Lung",
  category = "GO_TERMS",
  regulation = "DOWN",
  max_terms_show = 40,
  global = FALSE
)

identical(DATASETS_PROCESSED$droplet_mixed$dataset@cci_detected,
          DATASETS_PROCESSED2$droplet_mixed$dataset@cci_detected)

BuildNetwork(
  DATASETS_PROCESSED$droplet_mixed$dataset,
  network_type = "ORA",
  network_layout = "celltypes",
  ID = "Spleen"
)

test_calico <- copy(DATASETS_PROCESSED$calico$dataset)

test_calico <- RunORA(
  object = test_calico,
  categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "KEGG_PWS", "ER_CELL_FAMILY", "GENAGE"),
  overwrite = TRUE,
  stringent_or_default = "default",
  stringent_logfc_threshold = NULL,
  global = FALSE
)


## Do stringent ORA with logfc > log(2) ####

DATASETS_COMBINED_log13 <- lapply(
  DATASETS_COMBINED_log13,
  function(i) {
    RunORA(
      object = i,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ER_CELL_FAMILY", "GENAGE"),
      stringent = "stringent",
      stringent_logfc_threshold = log(2),
      verbose = TRUE,
      global = FALSE
    )
  }
)

DATASETS_COMBINED_log13 <- lapply(
  DATASETS_COMBINED_log13,
  function(i) {
    RunORA(
      object = i,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ID", "ER_CELL_FAMILY", "GENAGE"),
      stringent = "stringent",
      stringent_logfc_threshold = log(2),
      verbose = TRUE,
      global = TRUE
    )
  }
)

DATASETS_COMBINED_log15 <- lapply(
  DATASETS_COMBINED_log15,
  function(i) {
    RunORA(
      object = i,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ER_CELL_FAMILY", "GENAGE"),
      stringent = "stringent",
      stringent_logfc_threshold = log(2),
      verbose = TRUE,
      global = FALSE
    )
  }
)

DATASETS_COMBINED_log15 <- lapply(
  DATASETS_COMBINED_log15,
  function(i) {
    RunORA(
      object = i,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ID", "ER_CELL_FAMILY", "GENAGE"),
      stringent = "stringent",
      stringent_logfc_threshold = log(2),
      verbose = TRUE,
      global = TRUE
    )
  }
)





