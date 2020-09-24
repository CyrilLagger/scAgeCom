####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - July 2020
## ursu_eugen@hotmail.com
##
## Utility functions to assign the CCI to their
## correct category based on various cutoffs.
##
####################################################
##

bind_tissues <- function(
  path,
  list_of_tissues,
  pre_filtering = TRUE
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
        if(pre_filtering) {
          dt <- dt[LR_DETECTED_young == TRUE | LR_DETECTED_old == TRUE]
        }
        return(dt)
        # if(is_log) {
        #   dt[, LR_LOGFC := LR_SCORE_old - LR_SCORE_young]
        # } else {
        #   dt[, LR_LOGFC := log(LR_SCORE_old/LR_SCORE_young)]
        # }
        # dt[, LR_LOGFC_ABS := abs(LR_LOGFC)]
        # dt[, LR_CELLTYPES := paste(L_CELLTYPE, R_CELLTYPE, sep = "_")]
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
    use.names = TRUE,
    fill = TRUE
  )
}

analyze_CCI_per_tissue <- function(
  data, 
  cutoff_quantile,
  cutoff_specificity_young = 0.05,
  cutoff_specificity_old = 0.05,
  cutoff_pval = 0.05,
  cutoff_logFC_abs = log(1.1),
  reassignment = NULL,
  is_log = TRUE
) {
  message("Analyzing CCI...")
  data <- copy(data)
  data <- preprocess_results(data, is_log)
  tissues <- unique(data$TISSUE)
  res <- rbindlist(
    l = lapply(
      tissues,
      function(tiss) {
        temp <- data[TISSUE == tiss]
        cutoff_score <- quantile(c(temp$LR_SCORE_old, temp$LR_SCORE_young), cutoff_quantile)
        temp <- analyze_detected(
          temp, 
          cutoff_score,
          cutoff_specificity_young,
          cutoff_specificity_old
        )
        temp <- analyze_significant(
          temp,
          cutoff_pval,
          cutoff_logFC_abs
        )
        temp <- add_case_type(
          temp,
          cutoff_score,
          cutoff_specificity_young,
          cutoff_specificity_old
        )
        return(temp)
      }
    ),
    use.names = TRUE
  )
  return(res)
}


analyze_CCI <- function(
  data, 
  cutoff_score,
  cutoff_specificity_young = 0.05,
  cutoff_specificity_old = 0.05,
  cutoff_pval = 0.05,
  cutoff_logFC_abs = log(1.1),
  reassignment = NULL,
  is_log = TRUE
) {
  message("Analyzing CCI...")
  data <- copy(data)
  data <- preprocess_results(data, is_log)
  data <- analyze_detected(
    data, 
    cutoff_score,
    cutoff_specificity_young,
    cutoff_specificity_old
  )
  data <- analyze_significant(
    data,
    cutoff_pval,
    cutoff_logFC_abs
  )
  data <- add_case_type(
    data,
    cutoff_score,
    cutoff_specificity_young,
    cutoff_specificity_old
  )
  return(data)
}


preprocess_results <- function(
  data,
  is_log = TRUE
) {
  message("Preprocessing results...")
  if (!all(class(data) == c("data.table", "data.frame"))) {
    stop("data should be data.table")
  }
  if(is_log) {
    data[, LOGFC := LR_SCORE_old - LR_SCORE_young]
  } else {
    data[, LOGFC := log(LR_SCORE_old / LR_SCORE_young)]
  }
  data[, LOGFC_ABS := abs(LOGFC)]
  data[, LR_CELLTYPE := paste(L_CELLTYPE, R_CELLTYPE, sep = "_")]
  data[, LR_NAME := list(sapply(1:nrow(.SD), function(i) {
    temp1 <- c(LIGAND_1[[i]], LIGAND_2[[i]])
    temp1 <- temp1[!is.na(temp1)]
    temp1 <- paste0(temp1, collapse = "_")
    temp2 <- c(RECEPTOR_1[[i]], RECEPTOR_2[[i]], RECEPTOR_3[[i]])
    temp2 <- temp2[!is.na(temp2)]
    temp2 <- paste0(temp2, collapse = "_")
    return(paste(temp1, temp2, sep = ":"))
  }))]
  return(data)
}

analyze_detected <- function(
  data,
  cutoff_score,
  cutoff_specificity_young,
  cutoff_specificity_old
) {
  message("Analyzing detected...")
  # Can have separate column for specificity to analyze later, since it
  #  can reject interesting interactions that we may want to detect.
  data[, LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG := 
         (LR_DETECTED_young == TRUE) 
       & (LR_SCORE_young >=  cutoff_score)
       & (BH_PVAL_young <= cutoff_specificity_young)
       ]
  data[, LR_DETECTED_AND_SIGNIFICANT_IN_OLD := 
         (LR_DETECTED_old == TRUE) 
       & (LR_SCORE_old >=  cutoff_score)
       & (BH_PVAL_old <= cutoff_specificity_old)
       ]
  return(data)
}

analyze_significant <- function(
  data, 
  cutoff_pval,
  cutoff_logFC_abs
) {
  message("Analyzing significant...")
  data[, DIFFERENTIAL_EXPRESSED := 
         (LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG
          | LR_DETECTED_AND_SIGNIFICANT_IN_OLD)
       & (BH_PVAL_DIFF <= cutoff_pval)
       & (LOGFC_ABS >= cutoff_logFC_abs)
       ]
  data[, DIFFERENTIAL_DIRECTION := fifelse(LOGFC > 0, "UP", "DOWN")]
  return(data)
}

add_case_type <- function(
  data,
  cutoff_score,
  cutoff_specificity_young,
  cutoff_specificity_old
) {
  message("Analyzing type of interactions...")
  data[, CASE_TYPE := ifelse(
    LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG &
      LR_DETECTED_AND_SIGNIFICANT_IN_OLD &
      DIFFERENTIAL_EXPRESSED,
    ifelse(
      DIFFERENTIAL_DIRECTION == "UP",
      "TTTU",
      "TTTD"
    ),
    ifelse(
      LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG &
        LR_DETECTED_AND_SIGNIFICANT_IN_OLD &
        !DIFFERENTIAL_EXPRESSED,
      ifelse(
        DIFFERENTIAL_DIRECTION == "UP",
        "TTFU",
        "TTFD"
      ),
      ifelse(
        LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG &
          !LR_DETECTED_AND_SIGNIFICANT_IN_OLD &
          DIFFERENTIAL_EXPRESSED,
        ifelse(
          DIFFERENTIAL_DIRECTION == "DOWN",
          "TFTD",
          ifelse(
            sum(c(BH_PVAL_old > cutoff_specificity_old,
                  !LR_DETECTED_old,
                  LR_SCORE_old < cutoff_score
            )
            ) == 1,
            "TTTU",
            "FFF"
          )
        ),
        ifelse(
          LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG & 
            !LR_DETECTED_AND_SIGNIFICANT_IN_OLD &  
            !DIFFERENTIAL_EXPRESSED,
          ifelse(
            sum(c(BH_PVAL_old > cutoff_specificity_old,
                  !LR_DETECTED_old,
                  LR_SCORE_old < cutoff_score
            )
            ) == 1,
            ifelse(
              DIFFERENTIAL_DIRECTION == "UP",
              "TTFU",
              "TTFD"
            ),
            "FFF"
          ),
          ifelse(
            !LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG &
              LR_DETECTED_AND_SIGNIFICANT_IN_OLD & 
              DIFFERENTIAL_EXPRESSED,
            ifelse(
              DIFFERENTIAL_DIRECTION == "UP",
              "FTTU",
              ifelse(
                sum(c(BH_PVAL_young > cutoff_specificity_young,
                      !LR_DETECTED_young,
                      LR_SCORE_young < cutoff_score
                )
                ) == 1,
                "TTTD",
                "FFF"
              )
            ),
            ifelse(
              !LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG &
                LR_DETECTED_AND_SIGNIFICANT_IN_OLD &
                !DIFFERENTIAL_EXPRESSED,
              ifelse(
                sum(c(BH_PVAL_young > cutoff_specificity_young,
                      !LR_DETECTED_young,
                      LR_SCORE_young < cutoff_score
                )
                ) == 1,
                ifelse(
                  DIFFERENTIAL_DIRECTION == "UP",
                  "TTFU",
                  "TTFD"
                ),
                "FFF"
              ),
              "FFF"
            )
          )
        )
      )
    )
  )
  ]
  
  return(data)
}



#################################


#' Creating a list of column names.
#' @return list
#' @export
get_default_colnames <- function() {
  
  cols = list(
    "TISSUES" = "TISSUE",
    "LR_SORTED" = "LR_SORTED",
    
    #"L_GENE" = "L_GENE",
    #"R_GENE" = "R_GENE",
    
    "LIGAND_CELLTYPE" = "L_CELLTYPE",
    "RECEPTOR_CELLTYPE" = "R_CELLTYPE",
    "LR_SCORE_YOUNG" = "LR_SCORE_young",
    "LR_SCORE_OLD" = "LR_SCORE_old",
    "BH_PVAL" = "BH_PVAL_DIFF",
    "RAW_PVAL" = "PVAL_DIFF",
    "LR_DETECTED_YOUNG" = "LR_DETECTED_young",
    "LR_DETECTED_OLD" = "LR_DETECTED_old",
    "BH_PVAL_SPECIFICITY_YOUNG" = "BH_PVAL_young",
    "BH_PVAL_SPECIFICITY_OLD" = "BH_PVAL_old",
    
    
    "LIGAND_RECEPTOR_CELLTYPES" = "LR_CELLTYPES",
    "LOGFC" = "LR_LOGFC",
    "LOGFC_ABS" = "LR_LOGFC_ABS",
    
    
    #"LIGAND_EXPRESSION_YOUNG" = "L_EXPRESSION_young",
    #"LIGAND_EXPRESSION_OLD" = "L_EXPRESSION_old",
    #"LIGAND_DETECTED_YOUNG" = "L_DETECTED_young",
    #"LIGAND_DETECTED_OLD" = "L_DETECTED_old",
    #"RECEPTOR_EXPRESSION_YOUNG" = "R_EXPRESSION_young",
    #"RECEPTOR_EXPRESSION_OLD" = "R_EXPRESSION_old",
    #"RECEPTOR_DETECTED_YOUNG" = "R_DETECTED_young",
    #"RECEPTOR_DETECTED_OLD" = "R_DETECTED_old",
    
    # New
    "LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG" = "LR_DETECTED_AND_SIGN_young",
    "LR_DETECTED_AND_SIGNIFICANT_IN_OLD" = "LR_DETECTED_AND_SIGN_old",
    "DIFFERENTIAL_EXPRESSED" = "DIFF_EXPRESSION",
    "DIFFERENTIAL_DIRECTION" = "DIRECTION",
    "CASE_TYPE" = "CASE_TYPE"
  )
  return(cols)
}

#' Check columns.
#' Stops execution if essential columns not found.
#' 
#' @param data
#' @param cols
#' @return void
#' 
check_columns <- function(
  data,
  cols
) {
  message("Checking column names...")
  
  essential_cols = c(
    cols$TISSUES,
    cols$LR_SORTED,
    cols$LIGAND_CELLTYPE,
    cols$RECEPTOR_CELLTYPE,
    cols$LR_SCORE_YOUNG,
    cols$LR_SCORE_OLD,
    cols$BH_PVAL,
    cols$RAW_PVAL,
    cols$LR_DETECTED_YOUNG,
    cols$LR_DETECTED_OLD#,
    #cols$LIGAND_EXPRESSION_YOUNG,
    #cols$LIGAND_EXPRESSION_OLD,
    #cols$LIGAND_DETECTED_YOUNG,
    #cols$LIGAND_DETECTED_OLD,
    #cols$RECEPTOR_EXPRESSION_YOUNG,
    #cols$RECEPTOR_EXPRESSION_OLD,
    #cols$RECEPTOR_DETECTED_YOUNG,
    #cols$RECEPTOR_DETECTED_OLD,
    #cols$RAW_PVAL_SPECIFICITY_YOUNG,
    #cols$RAW_PVAL_SPECIFICITY_OLD
  )
  
  essential_columns_are_detected = all(
    essential_cols %in% colnames(data)
  )
  
  if (!essential_columns_are_detected) {
    stop("At least one essential column not detected.")
  }
  
  return(0)
}

#' Extract first letter.
#' @param s, character vector or boolean
#' @return first letter as a character vector with length 1
extract_first_letter <- Vectorize(function(s) {
  strsplit(as.character(s), split="")[[1]][1]
}
)






