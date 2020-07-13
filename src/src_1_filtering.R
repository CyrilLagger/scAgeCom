####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - July 2020
##
## Utility functions to assign the CCI to their
## correct category based on various cutoffs.
##
####################################################
##

bind_tissues <- function(
  path,
  list_of_tissues,
  is_log
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
        if(is_log) {
          dt[, LR_LOGFC := LR_SCORE_old - LR_SCORE_young]
        } else {
          dt[, LR_LOGFC := log(LR_SCORE_old/LR_SCORE_young)]
        }
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

# Filtering process
filter_CCI <- function(
  dt,
  CUTOFF_SCORE_YOUNG,
  CUTOFF_SCORE_OLD,
  CUTOFF_LOGFC,
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05
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
  CUTOFF_LOGFC,
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05
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

filter_and_reassign_CCI <- function(
  dt,
  CUTOFF_SCORE_YOUNG,
  CUTOFF_SCORE_OLD,
  CUTOFF_LOGFC,
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05
) {
  filter_CCI(
    dt = dt,
    CUTOFF_SCORE_YOUNG = CUTOFF_SCORE_YOUNG,
    CUTOFF_SCORE_OLD = CUTOFF_SCORE_OLD,
    CUTOFF_LOGFC = CUTOFF_LOGFC,
    CUTOFF_CPDB_PVAL_YOUNG = CUTOFF_CPDB_PVAL_YOUNG,
    CUTOFF_CPDB_PVAL_OLD = CUTOFF_CPDB_PVAL_OLD,
    CUTOFF_PVAL_ADJ = CUTOFF_PVAL_ADJ
  )
  res_dt <- reassign_CCI(
    dt = dt[LR_KEEP_young == TRUE | LR_KEEP_old == TRUE],
    CUTOFF_SCORE_YOUNG = CUTOFF_SCORE_YOUNG,
    CUTOFF_SCORE_OLD = CUTOFF_SCORE_OLD,
    CUTOFF_LOGFC = CUTOFF_LOGFC,
    CUTOFF_CPDB_PVAL_YOUNG = CUTOFF_CPDB_PVAL_YOUNG,
    CUTOFF_CPDB_PVAL_OLD = CUTOFF_CPDB_PVAL_OLD,
    CUTOFF_PVAL_ADJ = CUTOFF_PVAL_ADJ
  )
  return(res_dt)
}
