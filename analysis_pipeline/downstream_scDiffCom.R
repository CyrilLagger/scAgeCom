####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - June 2020
##
## Combine the results from scDiffCom for the
## various datasets considered here. Then load the
## downstream parameters and proceed to filtering
## and overrepresenation analysis.
##
####################################################
##

library(Seurat)
library(scDiffCom)
library(data.table)
library(config)
library(conflicted)
library(circlize)

source("downstream_script/src/utils.R")

# `merge` from config is in conflict with base merge used in downstream.R
conflicted::conflict_prefer(name = "merge", winner = "base")

#load parameters
load_globals_from_config(R_CONFIG_ACTIVE = "default")

#Function to add useful column (and remove others) and rbind all tissues in a single data.table
bind_tissues <- function(
  path,
  list_of_tissues
) {
  data.table::rbindlist(
    l = lapply(
      X = list_of_tissues,
      FUN = function(
        tiss
      ) 
      {
        dt <- readRDS(paste0(path, "/diffcom_", tiss, ".rds"))
        dt[, TISSUE := tiss]
        dt[, LR_LOGFC := LR_SCORE_old - LR_SCORE_young]
        dt[, LR_LOGFC_ABS := abs(LR_LOGFC)]
        dt[, LR_CELLTYPES := paste(L_CELLTYPE, R_CELLTYPE, sep = "_")]
        #dt[, c("L_EXPRESSION_young",
        #       "L_EXPRESSION_old",
        #       "R_EXPRESSION_young",
        #       "R_EXPRESSION_old",
        #       "L_DETECTED_young",
        #       "L_DETECTED_old",
        #       "R_DETECTED_young",
        #       "R_DETECTED_old",
        #       "PVAL_old", 
        #       "PVAL_young") := NULL]
      }),
    use.names = TRUE
  )
}

#datasets of interest
datasets <- c("calico", "tms_facs", "tms_droplet")
#path of datasets
data_path <- list(
  calico = "scDiffCom_results/diffcom_calico_size_factor_log_10000iter_mixed",
  tms_facs = "scDiffCom_results/diffcom_tms_facs_size_factor_log_10000iter_mixed",
  tms_droplet = "scDiffCom_results/diffcom_tms_droplet_size_factor_log_10000iter_mixed"
)
#tissues of interest
tissue_list <- list(
  calico = c("kidney", "lung", "spleen"),
  tms_facs = c("Aorta", "BAT", "Bladder", "Brain_Myeloid",
               "Brain_Non-Myeloid", "Diaphragm", "GAT",
               "Heart", "Kidney", "Large_Intestine",
               "Limb_Muscle", "Liver", "Lung", "Mammary_Gland",
               "Marrow", "MAT", "Pancreas", "SCAT", "Skin",
               "Spleen", "Thymus", "Tongue", "Trachea"),
  tms_droplet = c("Bladder", "Heart_and_Aorta", "Kidney",
                  "Limb_Muscle", "Liver", "Lung",
                  "Mammary_Gland", "Marrow", "Spleen",
                  "Thymus", "Tongue")
)


#list of all results
diffcom_results <- mapply(bind_tissues, data_path, tissue_list, SIMPLIFY = FALSE)
#diffcom_results <- lapply(diffcom_results, setDF)

#check that the format and colnames are OK
check_representation_inv(diffcom_results$calico)

#Cutoff exploration
df_detec_in_young_and_old <- lapply(
  diffcom_results,
  remove_undetected_interactions
)
#calico
explore_cutoffs_dir_calico <- "test_and_comparison/explore_cutoffs/calico/"
create_dir(explore_cutoffs_dir_calico)
explore_calico <- explore_filter_cutoffs(df_detec_in_young_and_old$calico, explore_cutoffs_dir_calico) 
#tms FACS
explore_cutoffs_dir_tms_facs <- "test_and_comparison/explore_cutoffs/tms_facs/"
create_dir(explore_cutoffs_dir_tms_facs)
explore_tms_facs <- explore_filter_cutoffs(df_detec_in_young_and_old$tms_facs, explore_cutoffs_dir_tms_facs) 
#tms Droplet
explore_cutoffs_dir_tms_droplet <- "test_and_comparison/explore_cutoffs/tms_droplet/"
create_dir(explore_cutoffs_dir_tms_droplet)
explore_tms_droplet <- explore_filter_cutoffs(df_detec_in_young_and_old$tms_droplet, explore_cutoffs_dir_tms_droplet)

explore_calico
explore_tms_facs
explore_tms_droplet

# Filtering process
filter_CCI <- function(
  dt,
  CUTOFF_SCORE_YOUNG,
  CUTOFF_SCORE_OLD,
  CUTOFF_CPDB_PVAL_YOUNG,
  CUTOFF_CPDB_PVAL_OLD,
  CUTOFF_PVAL_ADJ,
  CUTOFF_LOGFC
  ) {
  dt[, LR_KEEP_young := ifelse(
    (LR_DETECTED_young == TRUE & LR_SCORE_young >= CUTOFF_SCORE_YOUNG & BH_PVAL_young <= CUTOFF_CPDB_PVAL_YOUNG),
    TRUE,
    FALSE
  )]
  dt[, LR_KEEP_old := ifelse(
    (LR_DETECTED_old == TRUE & LR_SCORE_old >= CUTOFF_SCORE_OLD & BH_PVAL_old <= CUTOFF_CPDB_PVAL_OLD),
    TRUE,
    FALSE
  )]
  dt[, LR_KEEP_DIFF := ifelse(
    (LR_KEEP_young | LR_KEEP_old) & BH_PVAL_DIFF <= CUTOFF_PVAL_ADJ & LR_LOGFC_ABS >= CUTOFF_LOGFC,
    TRUE,
    FALSE
  )]
  dt[, DIRECTION := ifelse(LR_LOGFC > 0, "UP", "DOWN")]
  # dt[, LR_UP := ifelse(
  #   LR_KEEP_DIFF == TRUE & LR_KEEP_old == TRUE & LR_LOGFC > 0, 
  #   TRUE, 
  #   FALSE
  # )]
  # dt[, LR_DOWN := ifelse(
  #   LR_KEEP_DIFF == TRUE & LR_KEEP_young == TRUE & LR_LOGFC < 0, 
  #   TRUE, 
  #   FALSE
  # )]
  return(dt)
}

reassign_CCI <- function(
  dt,
  CUTOFF_SCORE_YOUNG,
  CUTOFF_SCORE_OLD,
  CUTOFF_CPDB_PVAL_YOUNG,
  CUTOFF_CPDB_PVAL_OLD,
  CUTOFF_PVAL_ADJ,
  CUTOFF_LOGFC
) {
  dt[, SIG_TYPE := ifelse(
    LR_KEEP_young == TRUE & LR_KEEP_old == TRUE &  LR_KEEP_DIFF == TRUE,
    ifelse(
      DIRECTION == "UP",
      "TTTU",
      "TTTD"
    ),
    ifelse(
      LR_KEEP_young == TRUE & LR_KEEP_old == TRUE &  LR_KEEP_DIFF == FALSE,
      "TTF",
      ifelse(
        LR_KEEP_young == TRUE & LR_KEEP_old == FALSE &  LR_KEEP_DIFF == TRUE,
        ifelse(
          DIRECTION == "DOWN",
          "TFT",
          ifelse(
            sum(c(BH_PVAL_old > CUTOFF_CPDB_PVAL_OLD, LR_DETECTED_old == FALSE, LR_SCORE_old < CUTOFF_SCORE_OLD)) == 1,
            "TTTU",
            "FFF"
          )
        ),
        ifelse(
          LR_KEEP_young == TRUE & LR_KEEP_old == FALSE &  LR_KEEP_DIFF == FALSE,
          ifelse(
            sum(c(BH_PVAL_old > CUTOFF_CPDB_PVAL_OLD, LR_DETECTED_old == FALSE, LR_SCORE_old < CUTOFF_SCORE_OLD)) == 1,
            "TTF",
            "FFF"
          ),
          ifelse(
            LR_KEEP_young == FALSE & LR_KEEP_old == TRUE &  LR_KEEP_DIFF == TRUE,
            ifelse(
              DIRECTION == "UP",
              "FTT",
              ifelse(
                sum(c(BH_PVAL_young > CUTOFF_CPDB_PVAL_YOUNG, LR_DETECTED_young == FALSE, LR_SCORE_young < CUTOFF_SCORE_YOUNG)) == 1,
                "TTTD",
                "FFF"
              )
            ),
            ifelse(
              LR_KEEP_young == FALSE & LR_KEEP_old == TRUE &  LR_KEEP_DIFF == FALSE,
              ifelse(
                sum(c(BH_PVAL_young > CUTOFF_CPDB_PVAL_YOUNG, LR_DETECTED_young == FALSE, LR_SCORE_young < CUTOFF_SCORE_YOUNG)) == 1,
                "TTF",
                "FFF"
              ),
              "FFF"
            )
          )
        )
      )
    )
  )]
}

quantile(diffcom_results$tms_facs[LR_DETECTED_young == TRUE,]$LR_SCORE_young, probs = seq(0, 1, 0.1))
quantile(diffcom_results$tms_facs[LR_DETECTED_old == TRUE,]$LR_SCORE_old, probs = seq(0, 1, 0.1))
quantile(diffcom_results$tms_droplet[LR_DETECTED_young == TRUE,]$LR_SCORE_young, probs = seq(0, 1, 0.1))
quantile(diffcom_results$tms_droplet[LR_DETECTED_old == TRUE,]$LR_SCORE_old, probs = seq(0, 1, 0.1))
quantile(diffcom_results$calico[LR_DETECTED_young == TRUE,]$LR_SCORE_young, probs = seq(0, 1, 0.1))
quantile(diffcom_results$calico[LR_DETECTED_old == TRUE ,]$LR_SCORE_old, probs = seq(0, 1, 0.1))

filter_CCI(
  dt = diffcom_results$tms_facs,
  CUTOFF_SCORE_YOUNG = 0.85, #case dependent
  CUTOFF_SCORE_OLD = 0.85, # case dependent
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05,
  CUTOFF_LOGFC = log(1.5)
)

filter_CCI(
  dt = diffcom_results$tms_droplet,
  CUTOFF_SCORE_YOUNG = 1.0, #case dependent
  CUTOFF_SCORE_OLD = 1.0, # case dependent
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05,
  CUTOFF_LOGFC = log(1.1)
)

filter_CCI(
  dt = diffcom_results$calico,
  CUTOFF_SCORE_YOUNG = 0.8, #case dependent
  CUTOFF_SCORE_OLD = 0.8, # case dependent
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05,
  CUTOFF_LOGFC = log(1.1)
)

diffcom_results_detected <- list(
  tms_facs = reassign_CCI(
    dt = diffcom_results$tms_facs[LR_KEEP_young == TRUE | LR_KEEP_old == TRUE,],
    CUTOFF_SCORE_YOUNG = 0.85, #case dependent
    CUTOFF_SCORE_OLD = 0.85, # case dependent
    CUTOFF_CPDB_PVAL_YOUNG = 0.05,
    CUTOFF_CPDB_PVAL_OLD = 0.05,
    CUTOFF_PVAL_ADJ = 0.05,
    CUTOFF_LOGFC = log(1.5)
  ),
  tms_droplet = reassign_CCI(
    dt = diffcom_results$tms_droplet[LR_KEEP_young == TRUE | LR_KEEP_old == TRUE,],
    CUTOFF_SCORE_YOUNG = 1.0, #case dependent
    CUTOFF_SCORE_OLD = 1.0, # case dependent
    CUTOFF_CPDB_PVAL_YOUNG = 0.05,
    CUTOFF_CPDB_PVAL_OLD = 0.05,
    CUTOFF_PVAL_ADJ = 0.05,
    CUTOFF_LOGFC = log(1.1)
  ),
  calico = reassign_CCI(
    dt = diffcom_results$calico[LR_KEEP_young == TRUE | LR_KEEP_old == TRUE,],
    CUTOFF_SCORE_YOUNG = 0.8, #case dependent
    CUTOFF_SCORE_OLD = 0.8, # case dependent
    CUTOFF_CPDB_PVAL_YOUNG = 0.05,
    CUTOFF_CPDB_PVAL_OLD = 0.05,
    CUTOFF_PVAL_ADJ = 0.05,
    CUTOFF_LOGFC = log(1.1)
  )
)

#detected interaction by LR_CELLTYPES
table(diffcom_results_detected$tms_facs[TISSUE == "Liver",]$SIG_TYPE)

diffcom_results_detected$tms_facs[, LR_T_CELLTYPES := paste(TISSUE, L_CELLTYPE, TISSUE, R_CELLTYPE, sep = "_" ) ]
distr_CCI_tms_facs <- dcast(diffcom_results_detected$tms_facs[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                            TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_CCI_tms_facs[, N := rowSums(.SD), .SDcols = colnames(distr_CCI_tms_facs)[-c(1,2)]]
distr_CCI_tms_facs[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_facs[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_facs[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_CCI_tms_facs[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_CCI_tms_facs[, appear := old - young]
distr_CCI_tms_facs[, net_reg := up - down]
hist(distr_CCI_tms_facs$appear, breaks = 50)
hist(distr_CCI_tms_facs$net_reg, breaks = 50)
hist(distr_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$appear, breaks = 50)
hist(distr_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$net_reg, breaks = 50)
ggplot(distr_CCI_tms_facs, aes(x = appear, y = net_reg)) + geom_point()

het_tms_facs <- merge.data.table(
  diffcom_results_detected$tms_facs[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                       .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                       by = "LR_T_CELLTYPES" ],
  diffcom_results_detected$tms_facs[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                       .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                       by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_tms_facs[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_tms_facs$diff_LR, breaks = 100)
hist(het_tms_facs$diff_L, breaks = 50)
hist(het_tms_facs$diff_R, breaks = 50)



diffcom_results_detected$tms_droplet[, LR_T_CELLTYPES := paste(TISSUE, L_CELLTYPE, TISSUE, R_CELLTYPE, sep = "_" ) ]
distr_CCI_tms_droplet <- dcast(diffcom_results_detected$tms_droplet[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                            TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_CCI_tms_droplet[, N := rowSums(.SD), .SDcols = colnames(distr_CCI_tms_droplet)[-c(1,2)]]
distr_CCI_tms_droplet[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_droplet[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_droplet[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_CCI_tms_droplet[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_CCI_tms_droplet[, appear := old - young]
distr_CCI_tms_droplet[, net_reg := up - down]
hist(distr_CCI_tms_droplet[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$appear, breaks = 50)
hist(distr_CCI_tms_droplet[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$net_reg, breaks = 50)
hist(distr_CCI_tms_droplet$appear, breaks = 50)
hist(distr_CCI_tms_droplet$net_reg, breaks = 50)
ggplot(distr_CCI_tms_droplet, aes(x = appear, y = net_reg)) + geom_point()

het_tms_droplet <- merge.data.table(
  diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_tms_droplet[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_tms_droplet$diff_LR, breaks = 100)
hist(het_tms_droplet$diff_L, breaks = 50)
hist(het_tms_droplet$diff_R, breaks = 50)

diffcom_results_detected$calico[, LR_T_CELLTYPES := paste(TISSUE, L_CELLTYPE, TISSUE, R_CELLTYPE, sep = "_" ) ]
distr_CCI_calico <- dcast(diffcom_results_detected$calico[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                            TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_CCI_calico[, N := rowSums(.SD), .SDcols = colnames(distr_CCI_calico)[-c(1,2)]]
distr_CCI_calico[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_CCI_calico[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_CCI_calico[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_CCI_calico[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_CCI_calico[, appear := old - young]
distr_CCI_calico[, net_reg := up - down]
hist(distr_CCI_calico$appear, breaks = 50)
hist(distr_CCI_calico$net_reg, breaks = 50)
ggplot(distr_CCI_calico, aes(x = appear, y = net_reg)) + geom_point()

het_calico <- merge.data.table(
  diffcom_results_detected$calico[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  diffcom_results_detected$calico[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_calico[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_calico$diff_LR, breaks = 100)
hist(het_calico$diff_L, breaks = 50)
hist(het_calico$diff_R, breaks = 50)


hist(distr_CCI_tms_facs$N, breaks = 50)
hist(distr_CCI_tms_droplet$N, breaks = 50)
hist(distr_CCI_calico$N, breaks = 30)

#maybe not all genes in calico!!!
hist(distr_CCI_tms_facs[TISSUE == "Kidney", ]$N, breaks = 50)
hist(distr_CCI_tms_droplet[TISSUE == "Kidney", ]$N, breaks = 50)
hist(distr_CCI_calico[TISSUE == "kidney", ]$N, breaks = 50)

#difference in detection
test_lung_sig[SIG_TYPE %in% c("TTTU", "TTTD", "TTF", "TFT"),]
test_lung_sig[SIG_TYPE %in% c("TTTU", "TTTD", "TTF", "FTT"),]


distr_LR_tms_facs <- dcast(diffcom_results_detected$tms_facs[!(SIG_TYPE %in% "FFF"), c("LR_GENES", "SIG_TYPE", "TISSUE") ],
                            LR_GENES ~ SIG_TYPE, value.var = "SIG_TYPE")

distr_LR_tiss_tms_facs <- merge.data.table(
  dcast(diffcom_results_detected$tms_facs[!(SIG_TYPE %in% "FFF"), c("LR_GENES", "SIG_TYPE", "TISSUE") ],
        TISSUE + LR_GENES ~ SIG_TYPE, value.var = "SIG_TYPE"),
  diffcom_results$tms_facs[ , length(unique(LR_CELLTYPES)), by = "TISSUE"],
  by = "TISSUE",
  all.x =  TRUE
)
cols <- c("FTT", "TFT", "TTF", "TTTD", "TTTU")
distr_LR_tiss_tms_facs[, N := rowSums(.SD) , .SDcols = cols]

distr_LR_tiss_tms_facs[, paste0("pct_", c(cols, "N")) := .SD/V1 , .SDcols = c(cols, "N")]

####################



ggplot(diffcom_results$tms_facs, aes(x = LR_SCORE_young)) + geom_histogram(bins = 100) + scale_y_log10()
ggplot(diffcom_results$tms_facs, aes(x = LR_SCORE_old)) + geom_histogram(bins = 150) + scale_y_log10()
ggplot(diffcom_results$tms_facs, aes(x = LR_LOGFC)) + geom_histogram(bins = 150) + scale_y_log10()

ggplot(diffcom_results$tms_facs[LR_DETECTED_young == TRUE,], aes(x = LR_SCORE_young)) + 
  geom_histogram(bins = 100) + scale_y_log10()
ggplot(diffcom_results$tms_facs[LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05,], aes(x = LR_SCORE_old)) + 
  geom_histogram(bins = 100) + scale_y_log10()
ggplot(diffcom_results$tms_facs[(LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05) |
                                  (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05),], aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_facs[(LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 0.85) |
                                  (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 0.85),],
       aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_facs[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 0.85) |
                                  (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 0.85)) &
                                  BH_PVAL_DIFF <= 0.05,],
       aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$calico[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 1) |
                                   (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 1)) &
                                  BH_PVAL_DIFF <= 0.05,],
       aes(x = LR_LOGFC)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_droplet[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 1) |
                                 (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 1)) &
                                BH_PVAL_DIFF <= 0.05,],
       aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_facs[LR_DETECTED_old == TRUE & BH_PVAL_old | LR_DETECTED_young == TRUE ,], aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()



quantile(diffcom_results$tms_facs[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 0.85) |
                                     (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 0.85)) &
                                    BH_PVAL_DIFF <= 0.05,]$LR_LOGFC_ABS, probs = seq(0,1, 0.1))


#non immune cells
diffcom_results_detected$tms_facs[, c("L_T_CELLTYPE", "R_T_CELLTYPE") := .(paste(TISSUE, L_CELLTYPE, sep = "_"), paste(TISSUE, R_CELLTYPE, sep = "_"))]

write.csv(unique(diffcom_results_detected$tms_facs[, .(L_T_CELLTYPE)]),
          file = "../../../../../facs_t_cell_type.csv",
          row.names = FALSE
)

identical(
  sort(unique(diffcom_results_detected$tms_facs[, L_T_CELLTYPE])),
  sort(unique(diffcom_results_detected$tms_facs[, R_T_CELLTYPE]))
)

cell_types_immune_tms_facs <- read.csv("../../../../../facs_t_cell_type_sorted.csv",
                                       header = TRUE, 
                                       stringsAsFactors = FALSE)
table(cell_types_immune_tms_facs$Immune)
anyNA(cell_types_immune_tms_facs)

setDT(cell_types_immune_tms_facs)
ct_nonIm_tms_facs <- cell_types_immune_tms_facs[Immune == FALSE, L_T_CELLTYPE]

tms_facs_nonIm <- diffcom_results_detected$tms_facs[L_T_CELLTYPE %in% ct_nonIm_tms_facs & R_T_CELLTYPE %in% ct_nonIm_tms_facs, ]

distr_nonIm_CCI_tms_facs <- dcast(tms_facs_nonIm[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                            TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_nonIm_CCI_tms_facs[, N := rowSums(.SD), .SDcols = colnames(distr_nonIm_CCI_tms_facs)[-c(1,2)]]
distr_nonIm_CCI_tms_facs[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_nonIm_CCI_tms_facs[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_nonIm_CCI_tms_facs[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_nonIm_CCI_tms_facs[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_nonIm_CCI_tms_facs[, appear := old - young]
distr_nonIm_CCI_tms_facs[, net_reg := up - down]
hist(distr_nonIm_CCI_tms_facs$appear, breaks = 50)
hist(distr_nonIm_CCI_tms_facs$net_reg, breaks = 50)
hist(distr_nonIm_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$appear, breaks = 50)
hist(distr_nonIm_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$net_reg, breaks = 50)
ggplot(distr_nonIm_CCI_tms_facs, aes(x = appear, y = net_reg)) + geom_point()

het_nonIm_tms_facs <- merge.data.table(
  tms_facs_nonIm[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  tms_facs_nonIm[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_nonIm_tms_facs[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_nonIm_tms_facs$diff_LR, breaks = 100)
hist(het_nonIm_tms_facs$diff_L, breaks = 50)
hist(het_nonIm_tms_facs$diff_R, breaks = 50)

################

test_lung <- test_d[TISSUE == "Lung",]
test_lung <- diffcom_results_detected$tms_facs[TISSUE == "Liver",]
test_lung <- tms_facs_nonIm[TISSUE == "Lung",]

test_lung[SIG_TYPE %in% c("TFT", "FTT", "TTTU", "TTTD") & L_GENE == "Ccl5" & R_GENE == "Cxcr6",]
test_lung[SIG_TYPE %in% c("TFT", "FTT", "TTTU", "TTTD") & L_GENE == "Mrc1" & R_GENE == "Ptprc" & L_CELLTYPE == "hepatocyte",]


test_lung <- test_lung[, c("L_GENE", "R_GENE", "L_CELLTYPE", "R_CELLTYPE",
                               "LR_KEEP_young", "LR_KEEP_old", "LR_LOGFC", "LR_CELLTYPES", "LR_GENES", "SIG_TYPE")]
test_lung[, L := paste("L", L_CELLTYPE, L_GENE, sep = "_")]
test_lung[, R := paste("R", R_CELLTYPE, R_GENE, sep = "_")]

test_lung_d <- test_lung[ SIG_TYPE != "FFF", ]
test_lung_up <- test_lung[SIG_TYPE %in% c("FTT", "TTTU"),  ]
test_lung_down <- test_lung[SIG_TYPE %in% c("TFT", "TTTD"),  ]

test_lung_d_up <- merge.data.table(test_lung_d[, .N, by = LR_CELLTYPES],
                 test_lung_up[, .N, by = LR_CELLTYPES],
                 by = "LR_CELLTYPES", 
                 all.x = TRUE
)
test_lung_d_up[is.na(test_lung_d_up)] <- 0
test_lung_d_up[, ratio := N.y/N.x]

test_lung_d_down <- merge.data.table(test_lung_d[, .N, by = LR_CELLTYPES],
                                   test_lung_down[, .N, by = LR_CELLTYPES],
                                   by = "LR_CELLTYPES", 
                                   all.x = TRUE
)
test_lung_d_down[is.na(test_lung_d_down)] <- 0
test_lung_d_down[, ratio := N.y/N.x]



test_lung_chord <- data.table::merge.data.table(y = test_lung_d[, .N, by = LR_CELLTYPES],
                                                   x = unique(test_lung_d[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
                                                   by = "LR_CELLTYPES",
                                                   all = FALSE)[, LR_CELLTYPES := NULL]
test_lung_chord[, L_CELLTYPE := paste0("L_", L_CELLTYPE)]
test_lung_chord[, R_CELLTYPE := paste0("R_", R_CELLTYPE)]
chordDiagram(test_lung_chord, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()


test_lung_up_chord <- data.table::merge.data.table(y = test_lung_d_up[, c("LR_CELLTYPES", "ratio")],
                             x = unique(test_lung_up[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
                             by = "LR_CELLTYPES",
                             all = FALSE)[, LR_CELLTYPES := NULL]
test_lung_up_chord[, L_CELLTYPE := paste0("L_", L_CELLTYPE)]
test_lung_up_chord[, R_CELLTYPE := paste0("R_", R_CELLTYPE)]
chordDiagram(test_lung_up_chord, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_up_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

test_lung_down_chord <- data.table::merge.data.table(y = test_lung_d_down[, c("LR_CELLTYPES", "ratio")],
                                                   x = unique(test_lung_down[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
                                                   by = "LR_CELLTYPES",
                                                   all = FALSE)[, LR_CELLTYPES := NULL]
test_lung_down_chord[, L_CELLTYPE := paste0("L_", L_CELLTYPE)]
test_lung_down_chord[, R_CELLTYPE := paste0("R_", R_CELLTYPE)]
chordDiagram(test_lung_down_chord, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_down_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

test_lung_up_chord_g <- data.table::merge.data.table(y = test_lung_up[, .N, by = LR_GENES],
                                                           x = unique(test_lung_up[, c("L_GENE", "R_GENE", "LR_GENES")]),
                                                           by = "LR_GENES",
                                                           all = FALSE)[, LR_GENES := NULL]
test_lung_up_chord_g[, L_GENE := paste0("L_", L_GENE)]
test_lung_up_chord_g[, R_GENE := paste0("R_", R_GENE)]
chordDiagram(test_lung_up_chord_g, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_up_chord_g))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

test_lung_down_chord_g <- data.table::merge.data.table(y = test_lung_down[, .N, by = LR_GENES][N > 7, ],
                                                     x = unique(test_lung_down[, c("L_GENE", "R_GENE", "LR_GENES")]),
                                                     by = "LR_GENES",
                                                     all = FALSE)[, LR_GENES := NULL]
test_lung_down_chord_g[, L_GENE := paste0("L_", L_GENE)]
test_lung_down_chord_g[, R_GENE := paste0("R_", R_GENE)]
chordDiagram(test_lung_down_chord_g, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_down_chord_g))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()






#try with circos
#test_lung_d_up <- test_lung_d[LR_UP == TRUE, c("L", "R") ]
#test_lung_d_down <- test_lung_d[LR_DOWN == TRUE, c("L", "R") ]

factor_total <- factor(c(sort(unique(test_lung_up$L)), sort(unique(test_lung_up$R))))
factor_LR <- factor(c(rep("L", length(unique(test_lung_up$L))), rep("R", length(unique(test_lung_up$R))) ))
factor_col_LR <- c(rep(rand_color(1), length(unique(test_lung_up$L))), rep(rand_color(1), length(unique(test_lung_up$R))) )

circos.par(
  cell.padding = c(0.00, 0, 0.00, 0),
  gap.degree = 0.,
  start.degree = 270,
  track.height = 0.01)
circos.initialize(factors = factor_total, xlim = c(0,1) )
#circos.track(ylim = c(0,0.1))
#circos.track(ylim = c(0,0.1))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ycenter,
              #CELL_META$sector.index,
              factor_LR[CELL_META$sector.numeric.index],
              col = factor_col_LR[CELL_META$sector.numeric.index]
              #facing = "inside", niceFacing = TRUE
  )
})
for(i in 1:nrow(test_lung_up)) {
  circos.link(test_lung_up$L[[i]], 0, test_lung_up$R[[i]], 0,
              col = "red",
              directional = 0)
}
for(i in 1:nrow(test_lung_down)) {
  circos.link(test_lung_down$L[[i]], 0, test_lung_down$R[[i]], 0,
              col = "blue",
              directional = 1)
}
circos.clear()

######ORA

test_lung_detected <- test_lung[SIG_TYPE %in% c("FTT", "TFT", "TTF", "TTTD", "TTTU"), ]
test_lung_sig <- test_lung[SIG_TYPE %in% c("FTT", "TFT", "TTTD", "TTTU"),]
colnames(test_lung_detected)

df_d <- setDF(test_lung_detected)
df_sig <- setDF(test_lung_sig)


ora_test <- analyze_overrepresentation(df_d, df_sig,
                           exclude_tissue_overrepresentation=TRUE,
                           adjust_pvals=FALSE)

orat_test_LR_GENES <- ora_test$L_GENE



########enrichment####
setDT(test_lung_sig)
up_genes <- c(unique(test_lung_sig[SIG_TYPE %in% c("FTT", "TTTU"),]$L_GENE),
  unique(test_lung_sig[SIG_TYPE %in% c("FTT", "TTTU"),]$R_GENE))
down_genes <- c(unique(test_lung_sig[SIG_TYPE %in% c("TFT", "TTTD"),]$L_GENE),
              unique(test_lung_sig[SIG_TYPE %in% c("TFT", "TTTD"),]$R_GENE))
intersect(up_genes, down_genes)

library(clusterProfiler)
library(org.Mm.eg.db)




test_lung_diffnet <- test_lung[, lapply(.SD, sum), by = LR_CELLTYPES, .SDcols = c("LR_UP", "LR_DOWN")]

test_lung_net <- test_lung[, lapply(.SD, sum), by = LR_CELLTYPES, .SDcols = c("LR_KEEP_young", "LR_KEEP_old", "LR_KEEP_DIFF")]


test_lung[LR_CELLTYPES == "B cell_dendritic cell" & (LR_KEEP_old == TRUE | LR_KEEP_young == TRUE),]


table(test$LR_KEEP_DIFF & test$LR_LOGFC > 0)
table(test$LR_KEEP_DIFF & test$LR_LOGFC < 0)
table(test$LR_KEEP_young)
table(test$LR_KEEP_old)





sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC > 0, ]$L_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC > 0, ]$R_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC > 0, ]$LR_GENES), decreasing = TRUE)

sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC < 0, ]$L_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC < 0, ]$R_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC < 0, ]$LR_GENES), decreasing = TRUE)

#Apoe
test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]
test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]

sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]$TISSUE), decreasing = TRUE)
sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]$TISSUE), decreasing = TRUE)

sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]$L_CELLTYPE), decreasing = TRUE)
sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]$L_CELLTYPE), decreasing = TRUE)

sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]$R_GENE), decreasing = TRUE)
sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]$R_GENE), decreasing = TRUE)


#############

which(test$LR_KEEP_old & !test$LR_KEEP_young & test$LR_KEEP_DIFF & (test$LR_LOGFC < 0))

diffcom_results$tms_facs[16265744, ]

df_detected <- run_filtering(df, filter_detected, 
                            'raw_to_detected', res_filtering_dir)



df_significant = run_filtering(df_detected, filter_significant, 
                               'detected_to_significant', res_filtering_dir)


filter_detected(df)
#summarize


df_detected = run_filtering(df, filter_detected, 
                            'raw_to_detected', res_filtering_dir)

###################################################
create_dir(RESULTS_DIR)
copy_config(RESULTS_DIR)

# logging by redirecting output to file
file_log = paste(RESULTS_DIR, "log.txt", sep='')
sink(file=file_log)


# FILTERING
res_filtering_dir = paste0(RESULTS_DIR, "/", "FILTERING")
create_dir(res_filtering_dir)
df_detected = run_filtering(df, filter_detected, 
                            'raw_to_detected', res_filtering_dir)
df_significant = run_filtering(df_detected, filter_significant, 
                               'detected_to_significant', res_filtering_dir)

ftable(test_d$LR_KEEP_young, test_d$LR_KEEP_old, test_d$LR_KEEP_DIFF, test_d$DIRECTION)

test_inc_TFTU <- test_d[
  LR_KEEP_young == TRUE &
    LR_KEEP_old == FALSE &
    LR_KEEP_DIFF == TRUE &
    DIRECTION == "UP",]
ftable(test_inc_TFTU$BH_PVAL_old > 0.05,
       test_inc_TFTU$LR_DETECTED_old == FALSE,
       test_inc_TFTU$LR_SCORE_old < 0.85)

test_inc_FTTD <- test_d[
  LR_KEEP_young == FALSE &
    LR_KEEP_old == TRUE &
    LR_KEEP_DIFF == TRUE &
    DIRECTION == "DOWN",]
ftable(test_inc_FTTD$BH_PVAL_young > 0.05,
       test_inc_FTTD$LR_DETECTED_young == FALSE,
       test_inc_FTTD$LR_SCORE_young < 0.85)


test_inc_TFFU <- test_d[
  LR_KEEP_young == TRUE &
    LR_KEEP_old == FALSE &
    LR_KEEP_DIFF == FALSE&
    DIRECTION == "UP",]
ftable(test_inc_TFFU$BH_PVAL_old > 0.05,
       test_inc_TFFU$LR_DETECTED_old == FALSE,
       test_inc_TFFU$LR_SCORE_old < 0.85)

test_inc_TFFD <- test_d[
  LR_KEEP_young == TRUE &
    LR_KEEP_old == FALSE &
    LR_KEEP_DIFF == FALSE&
    DIRECTION == "DOWN",]
ftable(test_inc_TFFD$BH_PVAL_old > 0.05,
       test_inc_TFFD$LR_DETECTED_old == FALSE,
       test_inc_TFFD$LR_SCORE_old < 0.85)

ftable(test_inc_TFFD$BH_PVAL_old > 0.05,
       test_inc_TFFD$LR_DETECTED_old == TRUE,
       test_inc_TFFD$LR_SCORE_old < 0.85)

#
test_inc_FTFD <- test_d[
  LR_KEEP_young == FALSE &
    LR_KEEP_old == TRUE &
    LR_KEEP_DIFF == FALSE&
    DIRECTION == "DOWN",]
ftable(test_inc_FTFD$BH_PVAL_young > 0.05,
       test_inc_FTFD$LR_DETECTED_young == FALSE,
       test_inc_FTFD$LR_SCORE_young < 0.85)

test_inc_FTFU <- test_d[
  LR_KEEP_young == FALSE &
    LR_KEEP_old == TRUE &
    LR_KEEP_DIFF == FALSE&
    DIRECTION == "UP",]
ftable(test_inc_FTFU$BH_PVAL_young > 0.05,
       test_inc_FTFU$LR_DETECTED_young == FALSE,
       test_inc_FTFU$LR_SCORE_young < 0.85)

