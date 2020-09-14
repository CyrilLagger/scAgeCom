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


# Filtering --------------------------------------------------------------------

#' Filtering of results data for statistically and scientifically significant
#'  interactions. The function acts as a mutator for data by adding relevant
#'  columns which can be acted upon for desired filtering.
#' 
#' @param data data.table with the ligand receptor interaction results
#' @param cols list of columnnames to be used.
#' @param cutoff_score_young numeric indicating the cutoff on ligand-receptor 
#' interaction score on young samples. Default is 0.85.
#' @param cutoff_score_old numeric indicating the cutoff on ligand-receptor 
#' interaction score on old samples. Default is 0.85.
#' @param cutoff_specificity_young numeric indicating the cutoff on specificity 
#' p-value obtained from cellphonedb for the young samples. Default is 0.05.
#' @param cutoff_specificity_old numeric indicating the cutoff on specificity 
#' p-value obtained from cellphonedb for the old samples. Default is 0.05.
#' @param cutoff_pval numeric indicating the cutoff on the significance p-value
#' obtained during results generation. Default is 0.05.
#' @param cutoff_logFC_abs numeric indicating the cutoff on absolute logFC which
#' is an indicator of scientific significance and interest. Default is log(1.1).
#'
#' @return mutated data.table with columns used for filtering.
#' @export
analyze_CCI <- function(
  data, 
  cols = get_default_colnames(),
  cutoff_score,
  cutoff_specificity_young = 0.05,
  cutoff_specificity_old = 0.05,
  cutoff_pval = 0.05,
  cutoff_logFC_abs = log(1.1),
  reassignment = NULL,
  is_log = TRUE
) {
  message("Analyzing CCI...")
  #data <- copy(data)
  
  data <- preprocess_results(data, cols, is_log)
  
  data <- analyze_detected(
    data, 
    cols,
    cutoff_score,
    cutoff_specificity_young,
    cutoff_specificity_old
  )
  
  data <- analyze_significant(
    data,
    cols,
    cutoff_pval,
    cutoff_logFC_abs
  )
  
  data <- add_case_type(
    data,
    cols,
    cutoff_score,
    cutoff_specificity_young,
    cutoff_specificity_old
  )
  
  return(data)
}

#' Preprocess the results data.table. Namely, the function mutates the 
#' data.table by introducing additionally derived columns for convenience: 
#'   cols$L_GENE, cols$R_GENE, cols$LIGAND_RECEPTOR_CELLTYPES, cols$LOGFC
#'   cols$LOGFC_ABS
#' 
#' @param data data.table
#' @return data.table the mutated data.table
#' @export
preprocess_results <- function(
  data,
  cols,
  is_log = TRUE
) {
  message("Preprocessing results...")
  
  if (!all(class(data) == c("data.table", "data.frame"))) {
    stop("data should be data.table")
  }
  
  # check all necessary columns exist
  check_columns(data, cols)
  
  # Define colnames to fit into data.table syntax. Haven't figured out how to
  #  nicely use the referenced values from the list with data.table.
  COL_LR_SORTED = cols$LR_SORTED
  #COL_L_GENE = cols$L_GENE
  #COL_R_GENE = cols$R_GENE
  #COL_LIGAND_RECEPTOR_CELLTYPES = cols$LIGAND_RECEPTOR_CELLTYPES
  COL_LIGAND_CELLTYPE = cols$LIGAND_CELLTYPE
  COL_RECEPTOR_CELLTYPE = cols$RECEPTOR_CELLTYPE
  COL_LR_SCORE_OLD = cols$LR_SCORE_OLD
  COL_LR_SCORE_YOUNG = cols$LR_SCORE_YOUNG
  COL_LOGFC = cols$LOGFC
  COL_LOGFC_ABS = cols$LOGFC_ABS
  
  # Add separate columns for ligand and receptor gene name
  #LR_sorted = data[[COL_LR_SORTED]]
  #splitter = function(string) strsplit(string, "_")[[1]]
  #LR_genes_split = t(sapply(LR_genes, splitter, USE.NAMES=FALSE))
  #data[, (COL_L_GENE) := LR_genes_split[, 1]]
  #data[, (COL_R_GENE) := LR_genes_split[, 2]]
  
  # Add column with joined ligand to receptor cell type
  #data[, (COL_LIGAND_RECEPTOR_CELLTYPES) := paste(
  #  get(COL_LIGAND_CELLTYPE), 
  #  get(COL_RECEPTOR_CELLTYPE), 
  #  sep = " : "  # didn't check if "_" can be found in celltypes names
  #)
  #]
  
  if(is_log) {
    data[, (COL_LOGFC) := get(COL_LR_SCORE_OLD) - get(COL_LR_SCORE_YOUNG)]
  } else {
    data[, (COL_LOGFC) := log(get(COL_LR_SCORE_OLD) / get(COL_LR_SCORE_YOUNG))]
  }
  data[, (COL_LOGFC_ABS) := abs(get(COL_LOGFC))]
  
  return(data)
}

#' Mutates data by adding columns indicating detected and significant 
#'  interactions; differentially expressed interactions between young and old
#'  samples; direction of difference (up or down-regulation).
#' 
#' @param data data.table with the ligand receptor interaction results
#' @param cols ...
#' @param cutoff_score_young numeric indicating the cutoff on ligand-receptor 
#' interaction score on young samples.
#' @param cutoff_score_old numeric indicating the cutoff on ligand-receptor 
#' interaction score on old samples.
#' @param cutoff_specificity_young numeric indicating the cutoff on specificity 
#' p-value obtained from cellphonedb for the young samples.
#' @param cutoff_specificity_old numeric indicating the cutoff on specificity 
#' p-value obtained from cellphonedb for the old samples.
#' 
#' @return data.table mutated data.table
#' @export
analyze_detected <- function(
  data,
  cols,
  cutoff_score,
  cutoff_specificity_young,
  cutoff_specificity_old
) {
  message("Analyzing detected...")
  
  # Define column variables
  COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG = cols$LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG
  COL_LR_DETECTED_YOUNG = cols$LR_DETECTED_YOUNG
  COL_LR_SCORE_YOUNG = cols$LR_SCORE_YOUNG
  COL_BH_PVAL_SPECIFICITY_YOUNG = cols$BH_PVAL_SPECIFICITY_YOUNG
  COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD = cols$LR_DETECTED_AND_SIGNIFICANT_IN_OLD
  COL_LR_DETECTED_OLD = cols$LR_DETECTED_OLD
  COL_LR_SCORE_OLD = cols$LR_SCORE_OLD
  COL_BH_PVAL_SPECIFICITY_OLD = cols$BH_PVAL_SPECIFICITY_OLD
  
  
  # Can have separate column for specificity to analyze later, since it
  #  can reject interesting interactions that we may want to detect.
  data[, (COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG) := 
         (get(COL_LR_DETECTED_YOUNG) == TRUE) 
       & (get(COL_LR_SCORE_YOUNG) >=  cutoff_score)
       & (get(COL_BH_PVAL_SPECIFICITY_YOUNG) <= cutoff_specificity_young)
       ]
  
  data[, (COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD) := 
         (get(COL_LR_DETECTED_OLD) == TRUE) 
       & (get(COL_LR_SCORE_OLD) >=  cutoff_score)
       & (get(COL_BH_PVAL_SPECIFICITY_OLD) <= cutoff_specificity_old)
       ]
  
  return(data)
}

#' Filtering of results data for statistically and scientifically significant
#'  interactions. Acts as mutator by adding relevant columns.
#' 
#' @param data data.table with the ligand receptor interaction results
#' @param cutoff_pval numeric indicating the cutoff on the significance p-value
#' obtained during results generation.
#' @param cutoff_logFC_abs numeric indicating the cutoff on absolute logFC which
#' is an indicator of scientific significance and interest.
#'
#' @return mutated data.table
#' @export
analyze_significant <- function(
  data, 
  cols,
  cutoff_pval,
  cutoff_logFC_abs
) {
  message("Analyzing significant...")
  
  # Define colnames
  COL_DIFFERENTIAL_EXPRESSED = cols$DIFFERENTIAL_EXPRESSED
  COL_DIFFERENTIAL_DIRECTION = cols$DIFFERENTIAL_DIRECTION
  COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG = cols$LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG
  COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD = cols$LR_DETECTED_AND_SIGNIFICANT_IN_OLD
  COL_BH_PVAL = cols$BH_PVAL
  COL_LOGFC = cols$LOGFC
  COL_LOGFC_ABS = cols$LOGFC_ABS
  
  data[, (COL_DIFFERENTIAL_EXPRESSED) := 
         (get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG)
          | get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD))
       & (get(COL_BH_PVAL) <= cutoff_pval)
       & (get(COL_LOGFC_ABS) >= cutoff_logFC_abs)
       ]
  
  data[, (COL_DIFFERENTIAL_DIRECTION) := fifelse(get(COL_LOGFC) > 0, "UP", "DOWN")]
  
  return(data)
}

#' ...
#' @param data
#' @param cols
#' @param reassignment
#' @return data.table with case type column
add_case_type <- function(
  data,
  cols,
  cutoff_score,
  cutoff_specificity_young,
  cutoff_specificity_old
  #reassignment = NULL
) {
  COL_CASE_TYPE = cols$CASE_TYPE
  COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG = cols$LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG
  COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD = cols$LR_DETECTED_AND_SIGNIFICANT_IN_OLD
  COL_DIFFERENTIAL_EXPRESSED = cols$DIFFERENTIAL_EXPRESSED
  COL_DIFFERENTIAL_DIRECTION = cols$DIFFERENTIAL_DIRECTION
  COL_LR_DETECTED_YOUNG = cols$LR_DETECTED_YOUNG
  COL_LR_SCORE_YOUNG = cols$LR_SCORE_YOUNG
  COL_BH_PVAL_SPECIFICITY_YOUNG = cols$BH_PVAL_SPECIFICITY_YOUNG
  COL_LR_DETECTED_OLD = cols$LR_DETECTED_OLD
  COL_LR_SCORE_OLD = cols$LR_SCORE_OLD
  COL_BH_PVAL_SPECIFICITY_OLD = cols$BH_PVAL_SPECIFICITY_OLD
  
  # data[, (COL_CASE_TYPE) := paste0(
  #   extract_first_letter(get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG)),
  #   extract_first_letter(get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD)),
  #   extract_first_letter(get(COL_DIFFERENTIAL_EXPRESSED)),
  #   extract_first_letter(get(COL_DIFFERENTIAL_DIRECTION))
  # )
  # ]
  
  # if (!is.null(reassignment)) {
  #   for (i in 1:length(reassignment)) {
  #     
  #     target_val = names(reassignment)[i]
  #     reassignment_val = reassignment[i]
  #     
  #     data[get(COL_CASE_TYPE) == target_val, 
  #          (COL_CASE_TYPE) := reassignment_val]
  #   }
  # }
  
  data[, (COL_CASE_TYPE) := ifelse(
    get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG) &
      get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD) &
      get(COL_DIFFERENTIAL_EXPRESSED),
    ifelse(
      get(COL_DIFFERENTIAL_DIRECTION) == "UP",
      "TTTU",
      "TTTD"
    ),
    ifelse(
      get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG) &
        get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD) &
        !get(COL_DIFFERENTIAL_EXPRESSED),
      ifelse(
        get(COL_DIFFERENTIAL_DIRECTION) == "UP",
        "TTFU",
        "TTFD"
      ),
      ifelse(
        get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG) &
          !get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD) &
          get(COL_DIFFERENTIAL_EXPRESSED),
        ifelse(
          get(COL_DIFFERENTIAL_DIRECTION) == "DOWN",
          "TFTD",
          ifelse(
            sum(c(get(COL_BH_PVAL_SPECIFICITY_OLD) > cutoff_specificity_old,
                  !get(COL_LR_DETECTED_OLD),
                  get(COL_LR_SCORE_OLD) < cutoff_score
                  )
                ) == 1,
            "TTTU",
            "FFF"
          )
        ),
        ifelse(
          get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG) & 
            !get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD) &  
            !get(COL_DIFFERENTIAL_EXPRESSED),
          ifelse(
            sum(c(get(COL_BH_PVAL_SPECIFICITY_OLD) > cutoff_specificity_old,
                  !get(COL_LR_DETECTED_OLD),
                  get(COL_LR_SCORE_OLD) < cutoff_score
                  )
                ) == 1,
            ifelse(
              get(COL_DIFFERENTIAL_DIRECTION) == "UP",
              "TTFU",
              "TTFD"
            ),
            "FFF"
          ),
          ifelse(
            !get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG) &
              get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD) & 
              get(COL_DIFFERENTIAL_EXPRESSED),
            ifelse(
              get(COL_DIFFERENTIAL_DIRECTION) == "UP",
              "FTTU",
              ifelse(
                sum(c(get(COL_BH_PVAL_SPECIFICITY_YOUNG) > cutoff_specificity_young,
                      !get(COL_LR_DETECTED_YOUNG),
                      get(COL_LR_SCORE_YOUNG) < cutoff_score
                      )
                    ) == 1,
                "TTTD",
                "FFF"
              )
            ),
            ifelse(
              !get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG) &
                get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD) &
                !get(COL_DIFFERENTIAL_EXPRESSED),
              ifelse(
                sum(c(get(COL_BH_PVAL_SPECIFICITY_YOUNG) > cutoff_specificity_young,
                      !get(COL_LR_DETECTED_YOUNG),
                      get(COL_LR_SCORE_YOUNG) < cutoff_score
                      )
                    ) == 1,
                ifelse(
                  get(COL_DIFFERENTIAL_DIRECTION) == "UP",
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
  )]
  
  return(data)
}

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






