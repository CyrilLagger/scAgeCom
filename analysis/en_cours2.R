
DATASETS_PROCESSED1 <- readRDS("../data_scAgeCom/analysis/outputs_data/TMS_scAgeCom_processed_bestORA.rds")
DATASETS_PROCESSED2 <- readRDS("../data_scAgeCom/analysis/outputs_data/TMS_scAgeCom_processed_bestORA_w30m.rds")

DATASETS_PROCESSED <- c(DATASETS_PROCESSED1, DATASETS_PROCESSED2)
#

ALL_CELLTYPES <- dcast.data.table(
  rbindlist(
    lapply(
      DATASETS_PROCESSED,
      function(i) {
        i$dataset@cci_detected[, uniqueN(EMITTER_CELLTYPE), by = "ID"]
      }
    ),
    idcol = "dataset"
  ),
  formula = ID ~ dataset,
  value.var = "V1"
)
ALL_CELLTYPES[is.na(ALL_CELLTYPES)] <- 0

## min cells by CCI

NCELLS_MIN_dt <- rbindlist(
  lapply(
    DATASETS_PROCESSED,
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
      dt[, c("NCELLS_MIN")]
    }
  ),
  idcol = "dataset"
)

NCELLS_MEAN_dt <- rbindlist(
  lapply(
    DATASETS_PROCESSED,
    function(i) {
      dt <- i$dataset@cci_detected
      dt[, NCELLS_MEAN :=
           (get(paste0("NCELLS_EMITTER_", i$dataset@parameters$seurat_condition_id$cond1_name)) +
                  get(paste0("NCELLS_EMITTER_", i$dataset@parameters$seurat_condition_id$cond2_name)) +
                  get(paste0("NCELLS_RECEIVER_", i$dataset@parameters$seurat_condition_id$cond1_name)) +
                  get(paste0("NCELLS_RECEIVER_", i$dataset@parameters$seurat_condition_id$cond2_name)))/4
         ]
      dt[, c("NCELLS_MEAN")]
    }
  ),
  idcol = "dataset"
)

ggplot(NCELLS_MIN_dt, aes(x = dataset, y = log10(NCELLS_MIN))) + geom_boxplot()
ggplot(NCELLS_MEAN_dt, aes(x = dataset, y = log10(NCELLS_MEAN))) + geom_boxplot()

##All ORA up-down per tissues

ALL_ORA_GO_TERMS <- rbindlist(
  lapply(
    DATASETS_PROCESSED[-c(5,9)],
    function(dataset) {
      dt_up <- dataset$dataset@ora_default$GO_TERMS[OR_UP > 1 & BH_P_VALUE_UP <= 0.05, c("ID", "VALUE")]
      dt_down <- dataset$dataset@ora_default$GO_TERMS[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05, c("ID", "VALUE")]
      dt <- rbindlist(
        list(UP = dt_up, DOWN = dt_down),
        idcol = "REGULATION"
      )
    }
  ),
  idcol = "DATASET"
)

ALL_ORA_KEGG_PWS <- rbindlist(
  lapply(
    DATASETS_PROCESSED[-c(5,9)],
    function(dataset) {
      dt_up <- dataset$dataset@ora_default$KEGG_PWS[OR_UP > 1 & BH_P_VALUE_UP <= 0.05, c("ID", "VALUE")]
      dt_down <- dataset$dataset@ora_default$KEGG_PWS[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05, c("ID", "VALUE")]
      dt <- rbindlist(
        list(UP = dt_up, DOWN = dt_down),
        idcol = "REGULATION"
      )
    }
  ),
  idcol = "DATASET"
)

ALL_ORA_ER_CELLTYPES <- rbindlist(
  lapply(
    DATASETS_PROCESSED[-c(5,9)],
    function(dataset) {
      dt_up <- dataset$dataset@ora_default$ER_CELLTYPES[OR_UP > 1 & BH_P_VALUE_UP <= 0.05, c("ID", "VALUE")]
      dt_down <- dataset$dataset@ora_default$ER_CELLTYPES[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05, c("ID", "VALUE")]
      dt <- rbindlist(
        list(UP = dt_up, DOWN = dt_down),
        idcol = "REGULATION"
      )
    }
  ),
  idcol = "DATASET"
)

ALL_ORA_ER_FAMILY <- rbindlist(
  lapply(
    DATASETS_PROCESSED[-c(5,9)],
    function(dataset) {
      dt_up <- dataset$dataset@ora_default$ER_CELL_FAMILY[OR_UP > 1 & BH_P_VALUE_UP <= 0.05, c("ID", "VALUE")]
      dt_down <- dataset$dataset@ora_default$ER_CELL_FAMILY[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05, c("ID", "VALUE")]
      dt <- rbindlist(
        list(UP = dt_up, DOWN = dt_down),
        idcol = "REGULATION"
      )
    }
  ),
  idcol = "DATASET"
)

ALL_ORA_LR_GENES <- rbindlist(
  lapply(
    DATASETS_PROCESSED[-c(5,9)],
    function(dataset) {
      dt_up <- dataset$dataset@ora_default$LR_GENES[OR_UP > 1 & BH_P_VALUE_UP <= 0.05, c("ID", "VALUE")]
      dt_down <- dataset$dataset@ora_default$LR_GENES[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05, c("ID", "VALUE")]
      dt <- rbindlist(
        list(UP = dt_up, DOWN = dt_down),
        idcol = "REGULATION"
      )
    }
  ),
  idcol = "DATASET"
)
ALL_ORA_LR_GENES[, REG := ifelse(REGULATION == "UP", 1, -1)]

test <- dcast.data.table(
  ALL_ORA_ER_FAMILY[DATASET %in% c("calico", "facs_mixed", "droplet_mixed"), .N, by = c("REGULATION", "VALUE")],
  formula = VALUE ~ REGULATION,
  value.var = "N"
)

test <- dcast.data.table(
  ALL_ORA_LR_GENES[DATASET %in% c("droplet_mixed", "droplet_mixed_w30"), .N, by = c("REGULATION", "VALUE")],
  formula = VALUE ~ REGULATION,
  value.var = "N"
)

test2 <- ALL_ORA_LR_GENES[DATASET %in% c("droplet_mixed", "droplet_mixed_w30")]

#Lung celltypes
LUNG_CELLTYPES <- lapply(
  DATASETS_PROCESSED,
  function(i) {
    unique(i$dataset@cci_detected[ID == "Lung", ]$EMITTER_CELLTYPE)
  }
)
LUNG_CELLTYPES_UPSET <- data.table(
  CELLTYPE = sort(unique(unlist(LUNG_CELLTYPES)))
)
LUNG_CELLTYPES_UPSET[, c(names(LUNG_CELLTYPES)) := lapply(LUNG_CELLTYPES, function(i) {
  ifelse(CELLTYPE %in% i, 1, 0)
})]

UpSetR::upset(
  LUNG_CELLTYPES_UPSET,
  nsets = 9,
  nintersects = 60,
  order.by = "freq"
)

## Lung ORA-UP

LUNG_ORA_GO_TERMS_UP <- lapply(
  DATASETS_PROCESSED,
  function(i) {
    dt <- i$dataset@ora_default$GO_TERMS[ID == "Lung", c("ID", "VALUE", "OR_UP", "BH_P_VALUE_UP")]
    dt <- dt[OR_UP > 1 & BH_P_VALUE_UP <= 0.05]
    dt$VALUE
  }
)

LUNG_ORA_GO_TERMS_UP_UPSET <- data.table(
  GO_TERM = sort(unique(unlist(LUNG_ORA_GO_TERMS_UP)))
)
LUNG_ORA_GO_TERMS_UP_UPSET[, c(names(LUNG_ORA_GO_TERMS_UP)) := lapply(LUNG_ORA_GO_TERMS_UP, function(i) {
  ifelse(GO_TERM %in% i, 1, 0)
})]

UpSetR::upset(
  LUNG_ORA_GO_TERMS_UP_UPSET[, -c(6, 10)],
  nsets = 7,
  nintersects = 60,
  order.by = "freq"
)

## Lung ORA-DOWN

LUNG_ORA_GO_TERMS_DOWN <- lapply(
  DATASETS_PROCESSED,
  function(i) {
    dt <- i$dataset@ora_default$GO_TERMS[ID == "Kidney", c("ID", "VALUE", "OR_DOWN", "BH_P_VALUE_DOWN")]
    dt <- dt[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05]
    dt$VALUE
  }
)

LUNG_ORA_GO_TERMS_DOWN <- lapply(
  DATASETS_PROCESSED,
  function(i) {
    dt <- i$dataset@ora_combined_default$LR_GENES[, c("VALUE", "OR_DOWN", "BH_P_VALUE_DOWN")]
    dt <- dt[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05]
    dt$VALUE
  }
)

LUNG_ORA_GO_TERMS_DOWN_UPSET <- data.table(
  GO_TERM = sort(unique(unlist(LUNG_ORA_GO_TERMS_DOWN)))
)
LUNG_ORA_GO_TERMS_DOWN_UPSET[, c(names(LUNG_ORA_GO_TERMS_DOWN)) := lapply(LUNG_ORA_GO_TERMS_DOWN, function(i) {
  ifelse(GO_TERM %in% i, 1, 0)
})]


UpSetR::upset(
  LUNG_ORA_GO_TERMS_DOWN_UPSET[, -c(6, 10)],
  nsets = 7,
  nintersects = 60,
  order.by = "freq"
)


UpSetR::upset(
  LUNG_ORA_GO_TERMS_UP_UPSET[, c(2,5,9)],
  nsets = 4,
  nintersects = 60,
  order.by = "freq"
)


UpSetR::upset(
  LUNG_ORA_GO_TERMS_UP_UPSET[, c(2,7,8,9)],
  nsets = 4,
  nintersects = 60,
  order.by = "freq"
)

UpSetR::upset(
  LUNG_ORA_GO_TERMS_UP_UPSET[, c(2,8,9)],
  nsets = 8,
  nintersects = 60
)

UpSetR::upset(
  LUNG_ORA_GO_TERMS_UP_UPSET[, 2:9],
  nsets = 8,
  nintersects = 60
)

