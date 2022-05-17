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
## Process tissue-specific results
##
####################################################
##

## Retrieve scDiffCom results ####

paths_scd_results <- list.dirs(
  path_scagecom_output_scdiffcom,
  recursive = FALSE
)
paths_scd_results

## Rename datasets ####

scd_dataset_names <- c(
  "Calico Droplet (male)",
  "TMS Droplet (female)",
  "TMS Droplet (male)",
  "TMS Droplet (mixed)",
  "TMS FACS (female)",
  "TMS FACS (male)",
  "TMS FACS (mixed)"
)
names(scd_dataset_names) <- scd_dataset_names

## Change some tissue names #####

dt_celltype_conversion[
  ,
  new_tissue := ifelse(
    tissue == "BAT", "Adipose_Brown",
    ifelse(
      tissue == "GAT", "Adipose_Gonadal",
      ifelse(
        tissue == "MAT", "Adipose_Mesenteric",
        ifelse(
          tissue == "SCAT", "Adipose_Subcutaneous",
          ifelse(
            tissue == "Brain_Myeloid", "Brain",
            ifelse(
              tissue == "Brain_Non-Myeloid", "Brain",
              tissue
            )
          )
        )
      )
    )
  )
]

## Clean scDiffCom results and add cell type families ####

fun_process_scd <- function(
  dataset_path
) {
  # retrieve and load each object in the dataset
  tissues <- gsub(".*scdiffcom_(.+)\\.rds.*", "\\1", list.files(dataset_path))
  tissues <- tissues[!grepl("md_", tissues)]
  print(tissues)
  dataset <- lapply(
    X = tissues,
    FUN = function(
      tiss
    ) {
      readRDS(paste0(dataset_path, "/scdiffcom_", tiss, ".rds"))
    }
  )
  # change some tissue names
  tissues[tissues == "BAT"] <- "Adipose_Brown"
  tissues[tissues == "GAT"] <- "Adipose_Gonadal"
  tissues[tissues == "MAT"] <- "Adipose_Mesenteric"
  tissues[tissues == "SCAT"] <- "Adipose_Subcutaneous"
  dataset <- lapply(
    seq_along(dataset),
    function(i) {
      dataset[[i]]@parameters$object_name <- tissues[[i]]
      dataset[[i]]
    }
  )
  names(dataset) <- tissues
  # erase raw CCI to make the object lighter
  dataset <- lapply(
    dataset,
    function(i) {
      scDiffCom::EraseRawCCI(i)
    }
  )
  # add cell type families ORA
  dataset <- lapply(
    dataset,
    function(i) {
      print(scDiffCom::GetParameters(i)$object_name)
      temp_cell_types <- unique(
        GetTableCCI(
          object = i,
          type = "detected",
          simplified = FALSE
        )$EMITTER_CELLTYPE
      )
      temp_meta_family <- unique(
        dt_celltype_conversion[
          new_tissue == GetParameters(i)$object_name &
            cell_ontology_final %in% temp_cell_types
        ][, c("cell_ontology_final", "cell_family")]
      )
      EMITTER_dt <- data.table::copy(temp_meta_family)
      setnames(
        EMITTER_dt,
        old = colnames(EMITTER_dt),
        new = c("EMITTER_CELLTYPE", "EMITTER_CELLFAMILY")
      )
      RECEIVER_dt <- data.table::copy(temp_meta_family)
      setnames(
        RECEIVER_dt,
        old = colnames(RECEIVER_dt),
        new = c("RECEIVER_CELLTYPE", "RECEIVER_CELLFAMILY")
      )
      ER_dt <- CJ(
        EMITTER_TYPE = unique(EMITTER_dt$EMITTER_CELLTYPE),
        RECEIVER_TYPE = unique(RECEIVER_dt$RECEIVER_CELLTYPE)
      )
      ER_dt[
        EMITTER_dt,
        on = "EMITTER_TYPE==EMITTER_CELLTYPE",
        EMITTER_FAMILY := i.EMITTER_CELLFAMILY
      ]
      ER_dt[
        RECEIVER_dt,
        on = "RECEIVER_TYPE==RECEIVER_CELLTYPE",
        RECEIVER_FAMILY := i.RECEIVER_CELLFAMILY
      ]
      ER_dt[, ER_CELLTYPES := paste(
        EMITTER_TYPE,
        RECEIVER_TYPE,
        sep = "_"
      )]
      ER_dt[, ER_CELLFAMILIES := paste(
        EMITTER_FAMILY,
        RECEIVER_FAMILY,
        sep = "_"
      )]
      ER_dt <- ER_dt[, c("ER_CELLTYPES", "ER_CELLFAMILIES")]
      print(EMITTER_dt)
      new_object <- RunORA(
        object = i,
        categories = c(
          "LRI", "LIGAND_COMPLEX", "RECEPTOR_COMPLEX", "ER_CELLTYPES",
          "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "GO_TERMS", "KEGG_PWS"
        ),
        extra_annotations = list(
          EMITTER_dt,
          RECEIVER_dt,
          ER_dt
        ),
        overwrite = FALSE,
        verbose = TRUE
      )
    }
  )
  return(dataset)
}

scds_datasets <- lapply(
  paths_scd_results,
  function(path) {
    fun_process_scd(path)
  }
)
names(scds_datasets) <- scd_dataset_names

## Create full CCI tables and add information ####

dt_cci_full <- rbindlist(
  lapply(
    scds_datasets,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            GetTableCCI(
              object = tissue,
              type = "detected",
              simplified = FALSE
            )
          }
        ),
        fill = TRUE,
        idcol = "tissue"
      )
    }
  ),
  idcol = "dataset"
)

# add log2fc on top of logfc and deal with infinite values
dt_cci_full[, LOG2FC_BASE := LOGFC * log2(exp(1))]
dt_cci_full[
  ,
  LOG2FC := {
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
    )},
  by = c(
    "dataset",
    "tissue"
  )
]

# add ligand and receptor log2fc
dt_cci_full[
  ,
  LOG2FC_L := log2(
    pmin(L1_EXPRESSION_OLD, L2_EXPRESSION_OLD, na.rm = TRUE)
    /
      pmin(L1_EXPRESSION_YOUNG, L2_EXPRESSION_YOUNG, na.rm = TRUE)
  )
]
dt_cci_full[
  ,
  LOG2FC_R := log2(
    pmin(R1_EXPRESSION_OLD, R2_EXPRESSION_OLD, R3_EXPRESSION_OLD, na.rm = TRUE)
    /
      pmin(
        R1_EXPRESSION_YOUNG, R2_EXPRESSION_YOUNG,
        R3_EXPRESSION_YOUNG,
        na.rm = TRUE
      )
  )
]
dt_cci_full[
  ,
  LOG2FC_L := {
    max_L <- ceiling(max(.SD[is.finite(LOG2FC_L)][["LOG2FC_L"]]))
    min_L <- floor(min(.SD[is.finite(LOG2FC_L)][["LOG2FC_L"]]))
    max_L <- max(max_L, -min_L)
    min_L <- min(-max_L, min_L)
    ifelse(
      is.infinite(LOG2FC_L) & LOG2FC_L > 0,
      max_L,
      ifelse(
        is.infinite(LOG2FC_L) & LOG2FC_L < 0,
        min_L,
        LOG2FC_L
      )
    )
  },
  by = c(
    "dataset",
    "tissue"
  )
]
dt_cci_full[
  ,
  LOG2FC_R := {
    max_R <- ceiling(max(.SD[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
    min_R <- floor(min(.SD[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
    max_R <- max(max_R, -min_R)
    min_R <- min(-max_R, min_R)
    ifelse(
      is.infinite(LOG2FC_R) & LOG2FC_R > 0,
      max_R,
      ifelse(
        is.infinite(LOG2FC_R) & LOG2FC_R < 0,
        min_R,
        LOG2FC_R
      )
    )
  },
  by = c(
    "dataset",
    "tissue"
  )
]

# add all GO names attached to each LRI

table(
  scDiffCom::LRI_mouse$LRI_curated_GO$GO_ID %in%
    scDiffCom::gene_ontology_level$ID
)
any(duplicated(scDiffCom::gene_ontology_level$ID))

dt_cci_full[
  dcast.data.table(
    copy(scDiffCom::LRI_mouse$LRI_curated_GO)[
      scDiffCom::gene_ontology_level,
      on = "GO_ID==ID",
      GO_NAME := i.NAME
    ][, c(1, 3)],
    LRI ~ .,
    value.var = "GO_NAME",
    fun.aggregate = paste0,
    collapse = ";"
  ),
  on = "LRI",
  GO_NAMES := i..
][
  ,
  GO_NAMES := paste0(";", GO_NAMES, ";")
]

dt_cci_full$GO_NAMES[1:3]

# add all KEGG names attached to each LRI

dt_cci_full[
  dcast.data.table(
    scDiffCom::LRI_mouse$LRI_curated_KEGG[, c(1, 3)],
    LRI ~ .,
    value.var = "KEGG_NAME",
    fun.aggregate = paste0,
    collapse = ";"
  ),
  on = "LRI",
  KEGG_NAMES := i..
][
  ,
  KEGG_NAMES := paste0(";", KEGG_NAMES, ";")
]

dt_cci_full$KEGG_NAMES[1:10]

## Create a full ORA table for analysis and add information ####

dt_ora_full <- rbindlist(
  lapply(
    scds_datasets,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            rbindlist(
              GetTableORA(
                object = tissue,
                categories = "all",
                simplified = FALSE
              ),
              idcol = "ORA_CATEGORY",
              fill = TRUE
            )
          }
        ),
        fill = TRUE,
        idcol = "tissue"
      )
    }
  ),
  idcol = "dataset"
)

# add cell types for ERI

dt_ora_full[
  ORA_CATEGORY == "ER_CELLTYPES",
  c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE") := list(
    sub("_.*", "", VALUE),
    sub(".*_", "", VALUE)
  )]

# add ORA regulation annotations

dt_ora_full[
  ,
  ':='(
    IS_UP = ifelse(
      OR_UP >= 1 & BH_P_VALUE_UP <= 0.05,
      TRUE,
      FALSE
    ),
    IS_DOWN = ifelse(
      OR_DOWN >= 1 & BH_P_VALUE_DOWN <= 0.05,
      TRUE,
      FALSE
    ),
    IS_FLAT = ifelse(
      OR_FLAT >= 1 & BH_P_VALUE_FLAT <= 0.05,
      TRUE,
      FALSE
    )
  )
]

dt_ora_full[
  ,
  ORA_REGULATION := ifelse(
    !IS_UP & !IS_DOWN & !IS_FLAT,
    "Not Over-represented",
    ifelse(
      !IS_UP & !IS_DOWN & IS_FLAT,
      "FLAT",
      ifelse(
        !IS_UP & IS_DOWN & !IS_FLAT,
        "DOWN",
        ifelse(
          IS_UP & !IS_DOWN & !IS_FLAT,
          "UP",
          ifelse(
            IS_UP & !IS_DOWN & IS_FLAT,
            "UP",
            ifelse(
              !IS_UP & IS_DOWN & IS_FLAT,
              "DOWN",
              ifelse(
                IS_UP & IS_DOWN & !IS_FLAT,
                "UP:DOWN",
                "UP:DOWN"
              )
            )
          )
        )
      )
    )
  )
]

## Prepare Figure "dataset summary" for the manuscript ####

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

## Prepare Supplementary Data (CCI classification) ####

dt_cci_classification <- lapply(
  paths_scd_results,
  function(path) {
    tissues <- gsub(".*scdiffcom_(.+)\\.rds.*", "\\1", list.files(path))
    tissues <- tissues[!grepl("md_", tissues)]
    print(tissues)
    dataset <- lapply(
      X = tissues,
      FUN = function(
        tiss
      ) {
        res <- readRDS(paste0(path, "/scdiffcom_", tiss, ".rds"))
      }
    )
    tissues[tissues == "BAT"] <- "Adipose_Brown"
    tissues[tissues == "GAT"] <- "Adipose_Gonadal"
    tissues[tissues == "MAT"] <- "Adipose_Mesenteric"
    tissues[tissues == "SCAT"] <- "Adipose_Subcutaneous"
    dataset <- lapply(
      seq_along(dataset),
      function(i) {
        dataset[[i]]@parameters$object_name <- tissues[[i]]
        dataset[[i]]
      }
    )
    names(dataset) <- tissues
    regulation_dt <- rbindlist(
      lapply(
        dataset,
        function(i) {
          cci_dt <- copy(i@cci_table_raw)
          dt <- cci_dt[
            IS_CCI_EXPRESSED_YOUNG == TRUE | IS_CCI_EXPRESSED_OLD == TRUE
          ]
          dt[
            ,
            c(
              "BH_P_VALUE_YOUNG",
              "BH_P_VALUE_OLD",
              "BH_P_VALUE_DE"
            ) := list(
              stats::p.adjust(P_VALUE_YOUNG, method = "BH"),
              stats::p.adjust(P_VALUE_OLD, method = "BH"),
              stats::p.adjust(P_VALUE_DE, method = "BH")
            )
          ]
          dt[
            ,
            c(
              "IS_CCI_SCORE_YOUNG",
              "IS_CCI_SCORE_OLD",
              "IS_CCI_SPECIFIC_YOUNG",
              "IS_CCI_SPECIFIC_OLD",
              "IS_DE_LOGFC",
              "IS_DE_SIGNIFICANT",
              "DE_DIRECTION",
              "IS_CCI_DETECTED_YOUNG",
              "IS_CCI_DETECTED_OLD",
              "IS_CCI_DE"
            ) := {
              threshold_score_temp <- stats::quantile(
                x = c(
                  .SD[IS_CCI_EXPRESSED_YOUNG == TRUE][["CCI_SCORE_YOUNG"]],
                  .SD[IS_CCI_EXPRESSED_OLD == TRUE][["CCI_SCORE_OLD"]]
                ),
                probs = 0.2
              )
              is_cci_score_1 <- (CCI_SCORE_YOUNG >= threshold_score_temp)
              is_cci_score_2 <- (CCI_SCORE_OLD >= threshold_score_temp)
              is_cci_specific_1 <- BH_P_VALUE_YOUNG <= 0.05
              is_cci_specific_2 <- BH_P_VALUE_OLD <= 0.05
              is_de_logfc <- LOGFC_ABS >= log(1.5)
              is_de_significant <- BH_P_VALUE_DE <= 0.05
              de_direction <- fifelse(LOGFC > 0, "UP", "DOWN")
              is_cci_detected_1 <- (IS_CCI_EXPRESSED_YOUNG == TRUE) &
                is_cci_score_1 & is_cci_specific_1
              is_cci_detected_2 <- (IS_CCI_EXPRESSED_OLD == TRUE) &
                is_cci_score_2 & is_cci_specific_2
              is_cci_de <- is_de_logfc & is_de_significant
              list(
                is_cci_score_1, is_cci_score_2,
                is_cci_specific_1, is_cci_specific_2,
                is_de_logfc, is_de_significant, de_direction,
                is_cci_detected_1, is_cci_detected_2,
                is_cci_de
              )
            }
          ]
          dt[
            ,
            REGULATION :=
              ifelse(
                !IS_CCI_DETECTED_YOUNG & !IS_CCI_DETECTED_OLD,
                "NOT_DETECTED",
                ifelse(
                  IS_DE_LOGFC & IS_DE_SIGNIFICANT & DE_DIRECTION == "UP",
                  "UP",
                  ifelse(
                    IS_DE_LOGFC & IS_DE_SIGNIFICANT & DE_DIRECTION == "DOWN",
                    "DOWN",
                    ifelse(
                      !IS_DE_LOGFC,
                      "FLAT",
                      ifelse(
                        IS_DE_LOGFC & !IS_DE_SIGNIFICANT,
                        "NSC",
                        "There is a problem here!"
                      )
                    )
                  )
                )
              )
          ]
          dt <- dt[
            ,
            c(
              "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LRI",
              "IS_CCI_EXPRESSED_YOUNG",
              "IS_CCI_EXPRESSED_OLD",
              "IS_CCI_SCORE_YOUNG",
              "IS_CCI_SCORE_OLD",
              "IS_CCI_SPECIFIC_YOUNG",
              "IS_CCI_SPECIFIC_OLD",
              "IS_DE_LOGFC",
              "IS_DE_SIGNIFICANT",
              "DE_DIRECTION",
              "IS_CCI_DETECTED_YOUNG",
              "IS_CCI_DETECTED_OLD",
              "IS_CCI_DE",
              "REGULATION"
            )
          ]
          cci_dt <- cci_dt[
            ,
            c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LRI")
          ]
          cci_dt <- merge.data.table(
            cci_dt,
            dt,
            by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LRI"),
            all.x = TRUE,
            sort = FALSE
          )
          cci_dt[
            ,
            ER_CELLTYPES := paste(
              EMITTER_CELLTYPE,
              RECEIVER_CELLTYPE,
              sep = "_")]
          cci_dt[, CCI := paste(ER_CELLTYPES, LRI, sep = "_")]
          cci_dt[
            ,
            DE_DIRECTION := ifelse(is.na(DE_DIRECTION), "NONE", DE_DIRECTION)
          ]
          cci_dt[
            ,
            REGULATION := ifelse(is.na(REGULATION), "NOT_DETECTED", REGULATION)
          ]
          cci_dt[is.na(cci_dt)] <- FALSE
          cci_dt[
            ,
            .N,
            by =  c(
              "IS_CCI_EXPRESSED_YOUNG",
              "IS_CCI_SCORE_YOUNG",
              "IS_CCI_SPECIFIC_YOUNG",
              "IS_CCI_EXPRESSED_OLD",
              "IS_CCI_SCORE_OLD",
              "IS_CCI_SPECIFIC_OLD",
              "IS_DE_LOGFC",
              "IS_DE_SIGNIFICANT",
              "DE_DIRECTION",
              "IS_CCI_DETECTED_YOUNG",
              "IS_CCI_DETECTED_OLD",
              "IS_CCI_DE",
              "REGULATION"
            )
          ]
        }
      ),
      idcol = "tissue"
    )
  }
)
names(dt_cci_classification) <- scd_dataset_names
dt_cci_classification <- rbindlist(
  dt_cci_classification,
  idcol = "dataset"
)
dt_cci_classification[REGULATION == "NOT_DETECTED"]
dt_cci_classification <- dt_cci_classification[REGULATION != "NOT_DETECTED"]
dt_cci_classification[
  ,
  PCT := N / sum(N) * 100,
  by = c("dataset", "tissue")
]
dt_cci_classification <- dt_cci_classification[
  !grepl("mixed", dataset)
]

fwrite(
  dt_cci_classification,
  paste0(
    path_scagecom_output,
    "Supplementary_Data_cci_classification.csv"
  )
)