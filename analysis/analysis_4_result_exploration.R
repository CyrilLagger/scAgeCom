####################################################
##
## Project: scAgeCom
##
## Last update - December 2020
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
path_input_directory <- "../data_scAgeCom/scdiffcom_results_30_11_2020"
path_output_directory <- "../data_scAgeCom/analysis/"

## Retrieve all the results for each tissue and dataset ####
RESULT_PATHS <- list.dirs(path_input_directory, recursive = FALSE)
RESULT_PATHS

#focus on mixed datasets for now (postpone sex analysis to later on) 
RESULT_PATHS <- RESULT_PATHS[c(1, 4, 8)]

DATASETS <- lapply(
  RESULT_PATHS, 
  function(path) {
    tissues <- gsub(".*scdiffcom_(.+)\\.rds.*", "\\1", list.files(path))
    temp <- lapply(
      X = tissues,
      FUN = function(
        tiss
      ) {
        res <- readRDS(paste0(path, "/scdiffcom_", tiss, ".rds"))
      }
    )
    names(temp) <- tissues
    return(temp)
  }
)

# Set the names manually, be careful to check what you are doing
RESULT_PATHS
#names(DATASETS) <- c("calico")
names(DATASETS) <- c("calico", "droplet_mixed", "facs_mixed")
#names(DATASETS) <- c("calico", "droplet_female", "droplet_male", "droplet_mixed", "droplet_sex",
#                     "facs_female", "facs_male", "facs_mixed", "facs_sex")

## Combine the results in a scDiffComCombined object ####

DATASETS_COMBINED <- lapply(
  seq_along(DATASETS),
  function(i) {
    Combine_scDiffCom(
      l = DATASETS[[i]],
      object_name = names(DATASETS)[[i]],
      verbose = TRUE
    )
  }
)
names(DATASETS_COMBINED) <- names(DATASETS)

rm(DATASETS)

## See the parameters ####
lapply(
  DATASETS_COMBINED,
  function(i) {
    parameters(i)
  }
)

## Inspect the behaviour in functions of the logfc threshold ####

LOGFC_behaviour <- lapply(
  DATASETS_COMBINED,
  function(data_comb) {
    temp_seq <- seq(1.1, 2, 0.1)
    temp_list <- lapply(
      temp_seq,
      function(thr) {
        temp_obj <- FilterCCI(
          object = data_comb,
          new_threshold_logfc = log(thr),
          skip_ora = TRUE
        )
        get_cci_detected(temp_obj)[, {
          temp = .N
          .SD[,.(pct=.N/temp*100),by=REGULATION_SIMPLE]
        }, by = ID ]
      }
    )
    names(temp_list) <- as.character(temp_seq)
    temp_dt <- rbindlist(
    temp_list,
    use.names = TRUE,
    idcol = "logfc_threshold"
    )
    temp_dt$logfc_threshold <- as.numeric(temp_dt$logfc_threshold)
    return(temp_dt)
  }
)

#save for later use
saveRDS(LOGFC_behaviour, "../data_scAgeCom/analysis/outputs_data/analysis4_logfc_behaviour.rds")

LOGFC_behaviour <- readRDS("../data_scAgeCom/analysis/outputs_data/analysis4_logfc_behaviour.rds")

# visualize results
ggplot(LOGFC_behaviour$calico[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, color = ID)) + geom_point()
ggplot(LOGFC_behaviour$droplet_mixed[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, color = ID)) + geom_point()
ggplot(LOGFC_behaviour$facs_mixed[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, color = ID)) + geom_point()
ggplot(LOGFC_behaviour$calico[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, group = logfc_threshold )) + geom_boxplot()
ggplot(LOGFC_behaviour$droplet_mixed[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, group = logfc_threshold )) + geom_boxplot()
ggplot(LOGFC_behaviour$facs_mixed[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, group = logfc_threshold )) + geom_boxplot()

## Prepare two analysis results with log(1.3) and log(1.5) ####

# Filter with new logfc thresholds
DATASETS_COMBINED_log13 <- lapply(
  DATASETS_COMBINED,
  function(i) {
    FilterCCI(
      object = i,
      new_threshold_logfc = log(1.3),
      skip_ora = FALSE,
      verbose = TRUE
    )
  }
)

rm(DATASETS_COMBINED)

DATASETS_COMBINED_log15 <- lapply(
  DATASETS_COMBINED_log13,
  function(i) {
    FilterCCI(
      object = i,
      new_threshold_logfc = log(1.5),
      skip_ora = FALSE,
      verbose = TRUE
    )
  }
)

## Add specific terms to cell-types and genes (only to CCI raw) ####

#cell-type families
cell_types_dt <- setDT(read.csv(paste0(path_output_directory, "inputs_data/scDiffCom_cell_types.csv"), stringsAsFactors = FALSE))

#genage longevity genes
genage_mouse <- setDT(read.csv(paste0(path_output_directory, "inputs_data/genage_mouse.tsv"), sep = "\t", header = TRUE))

#add the terms to the detected CCIs
DATASETS_COMBINED_log13 <- lapply(
  DATASETS_COMBINED_log13,
  function(dataset) {
    dt <- dataset@cci_detected
    dt[cell_types_dt, on = "EMITTER_CELLTYPE==scDiffCom.cell.type", EMITTER_CELL_FAMILY := i.Family...broad]
    dt[cell_types_dt, on = "RECEIVER_CELLTYPE==scDiffCom.cell.type", RECEIVER_CELL_FAMILY := i.Family...broad]
    dt[, ER_CELL_FAMILY := paste(EMITTER_CELL_FAMILY, RECEIVER_CELL_FAMILY, sep = "_")]
    dt[, GENAGE := ifelse(
      LIGAND_1 %in% genage_mouse$Gene.Symbol | LIGAND_2 %in% genage_mouse$Gene.Symbol |
        RECEPTOR_1 %in% genage_mouse$Gene.Symbol | RECEPTOR_2 %in% genage_mouse$Gene.Symbol,
      "YES",
      "NO"
    )
    ]
    dataset@cci_detected <- dt
    return(dataset)
  }
)

DATASETS_COMBINED_log15 <- lapply(
  DATASETS_COMBINED_log15,
  function(dataset) {
    dt <- dataset@cci_detected
    dt[cell_types_dt, on = "EMITTER_CELLTYPE==scDiffCom.cell.type", EMITTER_CELL_FAMILY := i.Family...broad]
    dt[cell_types_dt, on = "RECEIVER_CELLTYPE==scDiffCom.cell.type", RECEIVER_CELL_FAMILY := i.Family...broad]
    dt[, ER_CELL_FAMILY := paste(EMITTER_CELL_FAMILY, RECEIVER_CELL_FAMILY, sep = "_")]
    dt[,GENAGE := ifelse(
      LIGAND_1 %in% genage_mouse$Gene.Symbol | LIGAND_2 %in% genage_mouse$Gene.Symbol |
        RECEPTOR_1 %in% genage_mouse$Gene.Symbol | RECEPTOR_2 %in% genage_mouse$Gene.Symbol,
      "YES",
      "NO"
    )
    ]
    dataset@cci_detected <- dt
    return(dataset)
  }
)

## Redo ORA with new categories ####

DATASETS_COMBINED_log13 <- lapply(
  DATASETS_COMBINED_log13,
  function(dataset) {
    RunORA(
      object = dataset,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ER_CELL_FAMILY", "GENAGE"),
      overwrite = FALSE,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      global = FALSE
    )
  }
)

DATASETS_COMBINED_log13 <- lapply(
  DATASETS_COMBINED_log13,
  function(dataset) {
    RunORA(
      object = dataset,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ID", "ER_CELL_FAMILY", "GENAGE"),
      overwrite = FALSE,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      global = TRUE
    )
  }
)

DATASETS_COMBINED_log15 <- lapply(
  DATASETS_COMBINED_log15,
  function(dataset) {
    RunORA(
      object = dataset,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ER_CELL_FAMILY", "GENAGE"),
      overwrite = FALSE,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      global = FALSE
    )
  }
)

DATASETS_COMBINED_log15 <- lapply(
  DATASETS_COMBINED_log15,
  function(dataset) {
    RunORA(
      object = dataset,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ID", "ER_CELL_FAMILY", "GENAGE"),
      overwrite = FALSE,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      global = TRUE
    )
  }
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

## save the full objects with the raw CCIs ####
saveRDS(DATASETS_COMBINED_log13, "../data_scAgeCom/analysis/outputs_data/analysis4_DATASETS_COMBINED_log13.rds")
saveRDS(DATASETS_COMBINED_log15, "../data_scAgeCom/analysis/outputs_data/analysis4_DATASETS_COMBINED_log15.rds")

DATASETS_COMBINED_log13 <- readRDS("../data_scAgeCom/analysis/outputs_data/analysis4_DATASETS_COMBINED_log13.rds")
DATASETS_COMBINED_log15 <- readRDS("../data_scAgeCom/analysis/outputs_data/analysis4_DATASETS_COMBINED_log15.rds")

## Coherence between ligand vs receptor logfc ####

coherence_plot <- function(
  dataset,
  axis_lim
) {
  dt_filt <- copy(dataset@cci_detected)
  dt_filt[, LOG2FC_L := log2(
    pmin(
      L1_EXPRESSION_OLD,
      L2_EXPRESSION_OLD,
      na.rm = TRUE
    )
    /
      pmin(
        L1_EXPRESSION_YOUNG,
        L2_EXPRESSION_YOUNG,
        na.rm = TRUE
      )
  ) ]
  dt_filt[, LOG2FC_R := log2(
    pmin(
      R1_EXPRESSION_OLD,
      R2_EXPRESSION_OLD,
      R3_EXPRESSION_OLD,
      na.rm = TRUE
    )
    /
      pmin(
        R1_EXPRESSION_YOUNG,
        R2_EXPRESSION_YOUNG,
        R3_EXPRESSION_YOUNG,
        na.rm = TRUE
      )
  ) ]
  dt_filt[, LOG2FC_L := ifelse(is.infinite(LOG2FC_L) & LOG2FC_L > 0, axis_lim,
                               ifelse(is.infinite(LOG2FC_L) & LOG2FC_L < 0, -axis_lim, LOG2FC_L))]
  dt_filt[, LOG2FC_R := ifelse(is.infinite(LOG2FC_R) & LOG2FC_R > 0, axis_lim,
                               ifelse(is.infinite(LOG2FC_R) & LOG2FC_R < 0, -axis_lim, LOG2FC_R))]
  ggplot(dt_filt, aes(x = LOG2FC_L, y = LOG2FC_R, color = REGULATION_SIMPLE)) +
    geom_point() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlim(c(-axis_lim-1, axis_lim+1)) +
    ylim(c(-axis_lim-1,axis_lim+1))
}

coherence_plot(DATASETS_COMBINED_log15$calico, 5)
coherence_plot(DATASETS_COMBINED_log15$droplet_mixed, 9)
coherence_plot(DATASETS_COMBINED_log15$facs_mixed, 12)

## Check re-classified CCIs ####

#to do later on
#test <- DATASETS_COMBINED_log15$facs_mixed@cci_detected
#test_raw <- DATASETS_COMBINED_log15$facs_mixed@cci_raw
#test2 <- test[, .N, by = c("REGULATION", "DIFFERENTIAL_DIRECTION", "DIFFERENTIALLY_EXPRESSED", "CCI_DETECTED_AND_SIGNIFICANT_IN_YOUNG", "CCI_DETECTED_AND_SIGNIFICANT_IN_OLD" )]

## create light objects (removing the raw CCIs), cannot perform filtering on them, but can do ORA ####
DATASETS_COMBINED_log13 <- lapply(
  DATASETS_COMBINED_log13,
  function(i) {
    i@cci_raw <- list()
    return(i)
  }
)

DATASETS_COMBINED_log15 <- lapply(
  DATASETS_COMBINED_log15,
  function(i) {
    i@cci_raw <- list()
    return(i)
  }
)

## save the ligth objects without the raw CCIs ####
saveRDS(DATASETS_COMBINED_log13, "../data_scAgeCom/analysis/outputs_data/analysis4_DATASETS_COMBINED_log13_light.rds")
saveRDS(DATASETS_COMBINED_log15, "../data_scAgeCom/analysis/outputs_data/analysis4_DATASETS_COMBINED_log15_light.rds")


