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
  conds,
  pre_filtering = TRUE,
  min_cells = NULL
) {
  data.table::rbindlist(
    l = lapply(
      X = list_of_tissues,
      FUN = function(
        tiss
      ) 
      {
        dt <- readRDS(paste0(path, "/scdiffcom_", tiss, ".rds"))
        dt[, TISSUE := tiss]
        if(pre_filtering) {
          dt <- dt[get(paste0("LR_DETECTED_", conds[[1]])) == TRUE | get(paste0("LR_DETECTED_", conds[[2]])) == TRUE]
        }
        if(!is.null(min_cells)) {
          # Colnames
          L_NCELLS_COND1 = paste0("L_NCELLS_", conds[[1]])
          R_NCELLS_COND1 = paste0("R_NCELLS_", conds[[1]])
          L_NCELLS_COND2 = paste0("L_NCELLS_", conds[[2]])
          R_NCELLS_COND2 = paste0("R_NCELLS_", conds[[2]])
          
          dt <- dt[get(L_NCELLS_COND1) >= min_cells & 
                     get(R_NCELLS_COND1) >= min_cells &
                     get(L_NCELLS_COND2) >= min_cells &
                     get(R_NCELLS_COND2) >= min_cells]
        }
        return(dt)
      }),
    use.names = TRUE,
    fill = TRUE
  )
}

analyze_CCI_per_tissue <- function(
  data,
  conds,
  cutoff_quantile,
  cutoff_logFC_abs,
  cutoff_specificity = 0.05,
  cutoff_pval = 0.05,
  reassignment = NULL,
  is_log = TRUE,
  recompute_BH = FALSE
) {
  message("Analyzing CCI...")
  data <- copy(data)
  data <- preprocess_results(data, conds, is_log)
  tissues <- unique(data$TISSUE)
  
  # Colnames
  BH_PVAL_COND1 = paste0("BH_PVAL_", conds[[1]])
  BH_PVAL_COND2 = paste0("BH_PVAL_", conds[[2]])
  PVAL_COND1 = paste0("PVAL_", conds[[1]])
  PVAL_COND2 = paste0("PVAL_", conds[[2]])
  LR_SCORE_COND1 = paste0("LR_SCORE_", conds[[1]])
  LR_SCORE_COND2 = paste0("LR_SCORE_", conds[[2]])
  
  res <- rbindlist(
    l = lapply(
      tissues,
      function(tiss) {
        temp <- data[TISSUE == tiss]
        if(recompute_BH) {
          temp[, (BH_PVAL_DIFF) := p.adjust(PVAL_DIFF, method = "BH")]
          temp[, (BH_PVAL_COND1) := p.adjust(get(PVAL_COND1), method = "BH")]
          temp[, (BH_PVAL_COND2) := p.adjust(get(PVAL_COND2), method = "BH")]
        }
        cutoff_score <- quantile(c(temp[[LR_SCORE_COND1]], temp[[LR_SCORE_COND2]]), cutoff_quantile)
        temp <- analyze_detected(
          temp,
          conds,
          cutoff_score,
          cutoff_specificity
        )
        temp <- analyze_significant(
          temp,
          conds,
          cutoff_pval,
          cutoff_logFC_abs
        )
        temp <- add_case_type(
          temp,
          conds,
          cutoff_score,
          cutoff_specificity
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
  conds,
  cutoff_score,
  cutoff_specificity = 0.05,
  cutoff_pval = 0.05,
  cutoff_logFC_abs = log(1.1),
  reassignment = NULL,
  is_log = TRUE
) {
  message("Analyzing CCI...")
  data <- copy(data)
  data <- preprocess_results(data, conds, is_log)
  data <- analyze_detected(
    data, 
    conds,
    cutoff_score,
    cutoff_specificity
  )
  data <- analyze_significant(
    data,
    conds,
    cutoff_pval,
    cutoff_logFC_abs
  )
  data <- add_case_type(
    data,
    conds,
    cutoff_score,
    cutoff_specificity
  )
  return(data)
}


preprocess_results <- function(
  data,
  conds,
  is_log = TRUE
) {
  message("Preprocessing results...")
  
  # Colnames
  LR_SCORE_COND1 = paste0("LR_SCORE_", conds[[1]])
  LR_SCORE_COND2 = paste0("LR_SCORE_", conds[[2]])
  
  if (!all(class(data) == c("data.table", "data.frame"))) {
    stop("data should be data.table")
  }
  if(is_log) {
    data[, LOGFC := get(LR_SCORE_COND2) - get(LR_SCORE_COND1)]
  } else {
    data[, LOGFC := log(get(LR_SCORE_COND2) / get(LR_SCORE_COND1))]
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
  conds,
  cutoff_score,
  cutoff_specificity
) {
  message("Analyzing detected...")
  
  # Colnames
  LR_DETECT_AND_SIGN_COND1 = paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", conds[[1]])
  LR_DETECT_AND_SIGN_COND2 = paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", conds[[2]])
  LR_DETECT_COND1 = paste0("LR_DETECTED_", conds[[1]])
  LR_DETECT_COND2 = paste0("LR_DETECTED_", conds[[2]])
  BH_PVAL_COND1 = paste0("BH_PVAL_", conds[[1]])
  BH_PVAL_COND2 = paste0("BH_PVAL_", conds[[2]])
  LR_SCORE_COND1 = paste0("LR_SCORE_", conds[[1]])
  LR_SCORE_COND2 = paste0("LR_SCORE_", conds[[2]])
  
  data[, (LR_DETECT_AND_SIGN_COND1) := 
         (get(LR_DETECT_COND1) == TRUE) 
       & (get(LR_SCORE_COND1) >=  cutoff_score)
       & (get(BH_PVAL_COND1) <= cutoff_specificity)
       ]
  data[, (LR_DETECT_AND_SIGN_COND2) := 
         (get(LR_DETECT_COND2) == TRUE) 
       & (get(LR_SCORE_COND2) >=  cutoff_score)
       & (get(BH_PVAL_COND2) <= cutoff_specificity)
       ]
  return(data)
}

analyze_significant <- function(
  data,
  conds,
  cutoff_pval,
  cutoff_logFC_abs
) {
  message("Analyzing significant...")
  
  # Columns
  LR_DETECT_AND_SIGN_COND1 = paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", conds[[1]])
  LR_DETECT_AND_SIGN_COND2 = paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", conds[[2]])
  
  data[, DIFFERENTIAL_EXPRESSED := 
         (get(LR_DETECT_AND_SIGN_COND1)
          | get(LR_DETECT_AND_SIGN_COND2))
       & (BH_PVAL_DIFF <= cutoff_pval)
       & (LOGFC_ABS >= cutoff_logFC_abs)
       ]
  data[, DIFFERENTIAL_DIRECTION := fifelse(LOGFC > 0, "UP", "DOWN")]
  return(data)
}

add_case_type <- function(
  data,
  conds,
  cutoff_score,
  cutoff_specificity
) {
  message("Analyzing type of interactions...")
  
  # Colnames
  LR_DETECT_AND_SIGN_COND1 = paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", conds[[1]])
  LR_DETECT_AND_SIGN_COND2 = paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", conds[[2]])
  BH_PVAL_COND1 = paste0("BH_PVAL_", conds[[1]])
  BH_PVAL_COND2 = paste0("BH_PVAL_", conds[[2]])
  LR_DETECT_COND1 = paste0("LR_DETECTED_", conds[[1]])
  LR_DETECT_COND2 = paste0("LR_DETECTED_", conds[[2]])
  LR_SCORE_COND1 = paste0("LR_SCORE_", conds[[1]])
  LR_SCORE_COND2 = paste0("LR_SCORE_", conds[[2]])
  
  data[, CASE_TYPE := ifelse(
    get(LR_DETECT_AND_SIGN_COND1) &
      get(LR_DETECT_AND_SIGN_COND2) &
      DIFFERENTIAL_EXPRESSED,
    ifelse(
      DIFFERENTIAL_DIRECTION == "UP",
      "TTTU",
      "TTTD"
    ),
    ifelse(
      get(LR_DETECT_AND_SIGN_COND1) &
        get(LR_DETECT_AND_SIGN_COND2) &
        !DIFFERENTIAL_EXPRESSED,
      ifelse(
        DIFFERENTIAL_DIRECTION == "UP",
        "TTFU",
        "TTFD"
      ),
      ifelse(
        get(LR_DETECT_AND_SIGN_COND1) &
          !get(LR_DETECT_AND_SIGN_COND2) &
          DIFFERENTIAL_EXPRESSED,
        ifelse(
          DIFFERENTIAL_DIRECTION == "DOWN",
          "TFTD",
          ifelse(
            sum(c(get(BH_PVAL_COND2) > cutoff_specificity,
                  !get(LR_DETECT_COND2),
                  get(LR_SCORE_COND2) < cutoff_score
            )
            ) == 1,
            "TTTU",
            "FFF"
          )
        ),
        ifelse(
          get(LR_DETECT_AND_SIGN_COND1) & 
            !get(LR_DETECT_AND_SIGN_COND2) &  
            !DIFFERENTIAL_EXPRESSED,
          ifelse(
            sum(c(get(BH_PVAL_COND2) > cutoff_specificity,
                  !get(LR_DETECT_COND2),
                  get(LR_SCORE_COND2) < cutoff_score
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
            !get(LR_DETECT_AND_SIGN_COND1) &
              get(LR_DETECT_AND_SIGN_COND2) & 
              DIFFERENTIAL_EXPRESSED,
            ifelse(
              DIFFERENTIAL_DIRECTION == "UP",
              "FTTU",
              ifelse(
                sum(c(get(BH_PVAL_COND1) > cutoff_specificity,
                      !get(LR_DETECT_COND1),
                      get(LR_SCORE_COND1) < cutoff_score
                )
                ) == 1,
                "TTTD",
                "FFF"
              )
            ),
            ifelse(
              !get(LR_DETECT_AND_SIGN_COND1) &
                get(LR_DETECT_AND_SIGN_COND2) &
                !DIFFERENTIAL_EXPRESSED,
              ifelse(
                sum(c(get(BH_PVAL_COND1) > cutoff_specificity,
                      !get(LR_DETECT_COND1),
                      get(LR_SCORE_COND1) < cutoff_score
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
    "LR_SCORE_YOUNG" = "LR_SCORE_YOUNG",
    "LR_SCORE_OLD" = "LR_SCORE_OLD",
    "BH_PVAL" = "BH_PVAL_DIFF",
    "RAW_PVAL" = "PVAL_DIFF",
    "LR_DETECTED_YOUNG" = "LR_DETECTED_YOUNG",
    "LR_DETECTED_OLD" = "LR_DETECTED_OLD",
    "BH_PVAL_SPECIFICITY_YOUNG" = "BH_PVAL_YOUNG",
    "BH_PVAL_SPECIFICITY_OLD" = "BH_PVAL_OLD",
    
    
    "LIGAND_RECEPTOR_CELLTYPES" = "LR_CELLTYPES",
    "LOGFC" = "LR_LOGFC",
    "LOGFC_ABS" = "LR_LOGFC_ABS",
    
    
    #"LIGAND_EXPRESSION_YOUNG" = "L_EXPRESSION_YOUNG",
    #"LIGAND_EXPRESSION_OLD" = "L_EXPRESSION_OLD",
    #"LIGAND_DETECTED_YOUNG" = "L_DETECTED_YOUNG",
    #"LIGAND_DETECTED_OLD" = "L_DETECTED_OLD",
    #"RECEPTOR_EXPRESSION_YOUNG" = "R_EXPRESSION_YOUNG",
    #"RECEPTOR_EXPRESSION_OLD" = "R_EXPRESSION_OLD",
    #"RECEPTOR_DETECTED_YOUNG" = "R_DETECTED_YOUNG",
    #"RECEPTOR_DETECTED_OLD" = "R_DETECTED_OLD",
    
    # New
    "LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG" = "LR_DETECTED_AND_SIGN_YOUNG",
    "LR_DETECTED_AND_SIGNIFICANT_IN_OLD" = "LR_DETECTED_AND_SIGN_OLD",
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






