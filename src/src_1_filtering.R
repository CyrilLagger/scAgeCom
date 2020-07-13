####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - July 2020
##
## Two utility functions to assign the CCI to their
## correct category based on various cutoffs.
##
####################################################
##

# Filtering process
filter_CCI <- function(
  dt,
  CUTOFF_SCORE_YOUNG,
  CUTOFF_SCORE_OLD,
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05,
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
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05,
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
