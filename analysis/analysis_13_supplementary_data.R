####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## Prepare Supplementary Data
##
####################################################
##

## LRI database - general (Sup. Data 1-6) ####

supd_lri_human <- scDiffCom::LRI_human[c(1, 2, 3)]
names(supd_lri_human) <- paste0(names(supd_lri_human), "_human")
supd_lri_mouse <- scDiffCom::LRI_mouse[c(1, 2, 3)]
names(supd_lri_mouse) <- paste0(names(supd_lri_mouse), "_mouse")

openxlsx::write.xlsx(
  c(supd_lri_mouse, supd_lri_human),
  file = paste0(
    path_scagecom_output,
    "Supplementary_Data_LRI_db.xlsx"
  )
)

## Ligand/receptor - aging (Sup. Data 7) ####

supd_lr_aging <- copy(pmid_aging_lri_n)
supd_lr_aging[
  data.table(hagr_gene = hagr_genes, in_hagr = TRUE),
  on = "gene==hagr_gene",
  in_hagr := i.in_hagr
]
supd_lr_aging[is.na(supd_lr_aging)] <- FALSE
supd_lr_aging <- supd_lr_aging[order(-count)]

table(supd_lr_aging$in_hagr)
summary(supd_lr_aging$count)

openxlsx::write.xlsx(
  supd_lr_aging,
  file = paste0(
    path_scagecom_output,
    "Supplementary_Data_LR_aging.xlsx"
  )
)

## Cell number distribution (Sup. Data 8) #####

supd_cell_number_distr <- copy(dt_md_tms_age_sex_dc)
setnames(
  supd_cell_number_distr,
  old = c("dataset", "tissue", "cell_ontology_final",
          "female_OLD", "female_YOUNG", "male_OLD", "male_YOUNG", 
          "male_scd", "female_scd"),
  new = c("Dataset", "Tissue", "Cell type",
          "Cell Number (female old)", "Cell Number (female young)",
          "Cell Number (male old)", "Cell Number (male young)", 
          "Available in scAgeCom (male)", "Available in scAgeCom (female)"),
)

supd_cell_number_distr <- supd_cell_number_distr[
  ,
  c("Dataset", "Tissue", "Cell type",
    "Cell Number (female old)", "Cell Number (female young)",
    "Cell Number (male old)", "Cell Number (male young)", 
    "Available in scAgeCom (male)", "Available in scAgeCom (female)")
]

setcolorder(
  supd_cell_number_distr,
  c("Dataset", "Tissue", "Cell type",
    "Cell Number (female young)", "Cell Number (female old)", 
    "Cell Number (male young)", "Cell Number (male old)",  
    "Available in scAgeCom (female)", "Available in scAgeCom (male)")
)

openxlsx::write.xlsx(
  supd_cell_number_distr,
  file = paste0(
    path_scagecom_output,
    "Supplementary_Data_Cell_Number_Distr.xlsx"
  )
)

## LRI not detected ####

openxlsx::write.xlsx(
  LRI_not_detected,
  file = paste0(
    path_scagecom_output,
    "Supplementary_Data_LRI_not_detected.xlsx"
  )
)

## SDEA vs scDiffCom comparison ####

openxlsx::write.xlsx(
  dt_sdea_comp_count_list,
  file = paste0(
    path_scagecom_output,
    "Supplementary_Data_SDEA_comparison.xlsx"
  )
)

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

dt_cci_classification <- fread(
  paste0(
    path_scagecom_output,
    "Supplementary_Data_cci_classification.csv"
  )
)

######### Old code ############



## Do the DEG analysis ####

for (analysis in analysis_list) {
  message(paste0("Starting DEG analysis of ", analysis))
  output_dir <- paste0(
    dir_deg_output,
    "/deg_",
    analysis
  )
  if (!dir.exists(output_dir)) {
    message("Create output directory.")
    dir.create(output_dir)
  }
  message("Reading seurat object.")
  if (analysis == "calico_male") {
    seurat_kidney <- readRDS(dataset_path[["calico_kidney"]])
    seurat_kidney$age_group <- ifelse(
      seurat_kidney$age == "young",
      "YOUNG",
      "OLD"
    )
    seurat_lung <- readRDS(dataset_path[["calico_lung"]])
    seurat_lung$age_group <- ifelse(
      seurat_lung$age == "young",
      "YOUNG",
      "OLD"
    )
    seurat_spleen <- readRDS(dataset_path[["calico_spleen"]])
    seurat_spleen$age_group <- ifelse(
      seurat_spleen$age == "young",
      "YOUNG",
      "OLD"
    )
    tissue_list <- c("Kidney", "Lung", "Spleen")
    seurat_list <- list(
      Kidney = seurat_kidney,
      Lung = seurat_lung,
      Spleen = seurat_spleen
    )
    n_tissue <- 3
  } else if (analysis == "facs_female") {
    seurat_obj <- readRDS(dataset_path[["tms_facs"]])
    seurat_obj <- subset(seurat_obj, subset = sex == "female")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "BAT", "Bladder", "Brain", "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle",
      "Lung", "MAT", "Mammary_Gland", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "facs_male") {
    seurat_obj <- readRDS(dataset_path[["tms_facs"]])
    seurat_obj <- subset(seurat_obj, subset = sex == "male")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Aorta", "BAT", "Bladder", "Brain", "Diaphragm",
      "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle", "Liver",
      "Lung", "MAT", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "facs_mixed") {
    seurat_obj <- readRDS(dataset_path[["tms_facs"]])
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Aorta", "BAT", "Bladder", "Brain", "Diaphragm",
      "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle", "Liver",
      "Lung", "Mammary_Gland", "MAT", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_female") {
    seurat_obj <- readRDS(dataset_path[["tms_droplet"]])
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj <- subset(seurat_obj, subset = sex == "female")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Heart_and_Aorta", "Kidney", "Limb_Muscle",
      "Liver", "Lung", "Mammary_Gland", "Marrow",
      "Spleen", "Thymus"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_male") {
    seurat_obj <- readRDS(dataset_path[["tms_droplet"]])
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj <- subset(seurat_obj, subset = sex == "male")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Bladder", "Kidney", "Liver", "Lung",
      "Spleen", "Tongue"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_mixed") {
    seurat_obj <- readRDS(dataset_path[["tms_droplet"]])
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Bladder", "Heart_and_Aorta", "Kidney", "Liver", "Limb_Muscle",
      "Liver", "Lung", "Mammary_Gland", "Marrow",
      "Spleen", "Tongue", "Thymus"
    )
    n_tissue <- length(tissue_list)
  }
  for (i in 1:n_tissue) {
    tiss <- tissue_list[i]
    message(paste0(
      "DEG analysis of the ",
      tiss,
      ". Tissue ",
      i,
      " out of ",
      n_tissue,
      "."
    ))
    if (analysis == "calico_male") {
      seurat_tiss <- seurat_list[[tiss]]
    } else {
      if (tiss == "Brain") {
        seurat_tiss <- subset(
          seurat_obj,
          subset = tissue %in% c(
            "Brain_Myeloid",
            "Brain_Non-Myeloid"
          )
        )
      } else {
        cells_tiss <- colnames(seurat_obj)[which(seurat_obj$tissue == tiss)]
        seurat_tiss <- subset(seurat_obj, cells = cells_tiss)
      }
    }
    cell_type_id <- "cell_ontology_final"
    n_ct <- length(unique(seurat_tiss[[cell_type_id]][[1]]))
    DefaultAssay(seurat_tiss) <- "RNA"
    if (n_ct > 1) {
      message("Size-factor normalization:")
      seurat_tiss <- NormalizeData(seurat_tiss, assay = "RNA")
      Idents(seurat_tiss) <- seurat_tiss[[cell_type_id]]
      saveRDS(seurat_tiss[[]], file = paste0(output_dir, "/md_", tiss, ".rds"))
      future::plan(multicore, workers = 24)
      deg_res <- lapply(
        levels(seurat_tiss),
        function(ct_temp) {
          tryCatch(
            FindMarkers(
              object = seurat_tiss,
              ident.1 = "OLD",
              ident.2 = "YOUNG",
              group.by = "age_group",
              subset.ident = ct_temp,
              logfc.threshold = log2(1.2),
              test.use = "wilcox",
              min.pct = 0.1,
              min.cells.feature = 5,
              min.cells.group = 5,
              pseudocount.use = 0.001
            ),
            error = function(e) e,
            warning = function(w) w
          )
        }
      )
      names(deg_res) <- levels(seurat_tiss)
      message(paste0("Saving results for the ", tiss, "."))
      saveRDS(deg_res, file = paste0(output_dir, "/deg_", tiss, ".rds"))
    } else {
      message("Not enough cell types, not performing the analysis.")
    }
    future::plan(sequential)
  }
}

## Path of MAST/wilcox results ####

MAST_path <- "data_deg_mast_aging"
MAST_dataset_paths <- list.dirs(MAST_path, recursive = FALSE)
MAST_dataset_names <- c(
  "Calico Droplet (male)",
  "TMS Droplet (female)",
  "TMS Droplet (male)",
  "TMS Droplet (mixed)",
  "TMS FACS (female)",
  "TMS FACS (male)",
  "TMS FACS (mixed)"
)
names(MAST_dataset_names) <- MAST_dataset_names

wilcox_path <- "data_deg_wilcox_aging"
wilcox_dataset_paths <- list.dirs(wilcox_path, recursive = FALSE)
wilcox_dataset_names <- c(
  "Calico Droplet (male)",
  "TMS Droplet (female)",
  "TMS Droplet (male)",
  "TMS Droplet (mixed)",
  "TMS FACS (female)",
  "TMS FACS (male)",
  "TMS FACS (mixed)"
)
names(wilcox_dataset_names) <- wilcox_dataset_names

## Clean MAST/wilcox results in a single data.table ####

process_deg_dataset <- function(
  deg_dataset_path
) {
  # retrieve and load each object in the dataset
  tissues <- gsub(".*deg_(.+)\\.rds.*", "\\1", list.files(deg_dataset_path))
  tissues <- tissues[!grepl("md_", tissues)]
  dataset <- lapply(
    X = tissues,
    FUN = function(
      tiss
    ) {
      readRDS(paste0(deg_dataset_path, "/deg_", tiss, ".rds"))
    }
  )
  tissues[tissues == "BAT"] <- "Adipose_Brown"
  tissues[tissues == "GAT"] <- "Adipose_Gonadal"
  tissues[tissues == "MAT"] <- "Adipose_Mesenteric"
  tissues[tissues == "SCAT"] <- "Adipose_Subcutaneous"
  names(dataset) <- tissues
  # remove potential errors and transform to data.table
  dataset <- rbindlist(
    l = lapply(
      dataset,
      function(tiss) {
        rbindlist(
          l = lapply(
            tiss,
            function(ct) {
              if (!is.data.frame(ct)) {
                return(NULL)
              } else {
                as.data.table(ct, keep.rownames = "gene")
              }
            }
          ),
          idcol = "cell_type"
        )
      }
    ),
    idcol = "tissue"
  )
  return(dataset)
}

MAST_results <- rbindlist(
  l = lapply(
    setNames(
      MAST_dataset_paths,
      MAST_dataset_names
    ),
    process_deg_dataset
  ),
  idcol = "dataset"
)
setnames(
  MAST_results,
  old = c("pct.1", "pct.2"),
  new = c("pct_old", "pct_young")
)
saveRDS(MAST_results, "data_deg_mast_aging/MAST_results_dt.rds")

wilcox_results <- rbindlist(
  l = lapply(
    setNames(
      wilcox_dataset_paths,
      wilcox_dataset_names
    ),
    process_deg_dataset
  ),
  idcol = "dataset"
)
setnames(
  wilcox_results,
  old = c("pct.1", "pct.2"),
  new = c("pct_old", "pct_young")
)
saveRDS(wilcox_results, "data_deg_wilcox_aging/wilcox_results_dt.rds")

## Perform FDR correction by tissue-cell-type ####

wilcox_results[
  ,
  p_val_bh := p.adjust(p_val, method = "BH"),
  by = c("dataset", "tissue", "cell_type")
]

## DEG regulation ####

MAST_results[
  ,
  regulation := ifelse(
    p_val <= 0.05 & avg_log2FC > log2(1.5),
    "UP",
    ifelse(
      p_val <= 0.05 & avg_log2FC < -log2(1.5),
      "DOWN",
      "NO"
    )
  )
]
table(MAST_results$regulation)

wilcox_results[
  ,
  regulation := ifelse(
    p_val_bh <= 0.05 & avg_log2FC > log2(1.5),
    "UP",
    ifelse(
      p_val_bh <= 0.05 & avg_log2FC < -log2(1.5),
      "DOWN",
      "NO"
    )
  )
]
table(wilcox_results$regulation)