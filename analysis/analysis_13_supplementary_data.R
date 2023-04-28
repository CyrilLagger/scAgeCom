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

## LRI database - general ####

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

## Ligand/receptor - aging ####

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

## Cell number distribution #####

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

## CCI classification ####

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
