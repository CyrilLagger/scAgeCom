#' Downstream processing of intercellular interactions results ------------------
#'
#' It consists of several functions corresponding to steps:
#'   1) Filtering
#'   2) Overrepresentation analysis
#'   3) Frequent itemsets analysis
#'   4) Plots

# #' @importFrom utils globalVariables
# #' @importFrom ggplot2 ggproto GeomViolin
# #' @import jsonlite
# #' @include utilities.R
# #'
# NULL


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
analyze_CCI <- function(data, 
                       cols = get_default_colnames(),
                       cutoff_score_young = 0.85,
                       cutoff_score_old = 0.85,
                       cutoff_specificity_young = 0.05,
                       cutoff_specificity_old = 0.05,
                       cutoff_pval = 0.05,
                       cutoff_logFC_abs = log(1.1),
                       reassignment = NULL) {
    
    message("Analyzing CCI...")
    
    data <- preprocess_results(data, cols)
    
    data <- analyze_detected(
        data, 
        cols,
        cutoff_score_young,
        cutoff_score_old,
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
        reassignment
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
preprocess_results <- function(data, cols) {
    message("Preprocessing results...")
    
    if (!all(class(data) == c("data.table", "data.frame"))) {
        stop("data should data.table")
    }
    
    # check all necessary columns exist
    check_columns(data, cols)
    
    # Define colnames to fit into data.table syntax. Haven't figured out how to
    #  nicely use the referenced values from the list with data.table.
    COL_LR_GENES = cols$LR_GENES
    COL_L_GENE = cols$L_GENE
    COL_R_GENE = cols$R_GENE
    COL_LIGAND_RECEPTOR_CELLTYPES = cols$LIGAND_RECEPTOR_CELLTYPES
    COL_LIGAND_CELLTYPE = cols$LIGAND_CELLTYPE
    COL_RECEPTOR_CELLTYPE = cols$RECEPTOR_CELLTYPE
    COL_LR_SCORE_OLD = cols$LR_SCORE_OLD
    COL_LR_SCORE_YOUNG = cols$LR_SCORE_YOUNG
    COL_LOGFC = cols$LOGFC
    COL_LOGFC_ABS = cols$LOGFC_ABS

    # Add separate columns for ligand and receptor gene name
    LR_genes = data[[COL_LR_GENES]]
    splitter = function(string) strsplit(string, "_")[[1]]
    LR_genes_split = t(sapply(LR_genes, splitter, USE.NAMES=FALSE))
    data[, (COL_L_GENE) := LR_genes_split[, 1]]
    data[, (COL_R_GENE) := LR_genes_split[, 2]]
    
    # Add column with joined ligand to receptor cell type
    data[, (COL_LIGAND_RECEPTOR_CELLTYPES) := paste(
        get(COL_LIGAND_CELLTYPE), 
        get(COL_RECEPTOR_CELLTYPE), 
        sep = " : "  # didn't check if "_" can be found in celltypes names
    )
    ]
    
    data[, (COL_LOGFC) := get(COL_LR_SCORE_OLD) - get(COL_LR_SCORE_YOUNG)]
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
analyze_detected <- function(data, cols,
                           cutoff_score_young,
                           cutoff_score_old,
                           cutoff_specificity_young,
                           cutoff_specificity_old) {
    message("Analyzing detected...")
    
    # Define column variables
    COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG = cols$LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG
    COL_LR_DETECTED_YOUNG = cols$LR_DETECTED_YOUNG
    COL_LR_SCORE_YOUNG = cols$LR_SCORE_YOUNG
    COL_RAW_PVAL_SPECIFICITY_YOUNG = cols$RAW_PVAL_SPECIFICITY_YOUNG
    COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD = cols$LR_DETECTED_AND_SIGNIFICANT_IN_OLD
    COL_LR_DETECTED_OLD = cols$LR_DETECTED_OLD
    COL_LR_SCORE_OLD = cols$LR_SCORE_OLD
    COL_RAW_PVAL_SPECIFICITY_OLD = cols$RAW_PVAL_SPECIFICITY_OLD
    
    
    # Can have separate column for specificity to analyze later, since it
    #  can reject interesting interactions that we may want to detect.
    data[, (COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG) := 
             (get(COL_LR_DETECTED_YOUNG) == TRUE) 
             & (get(COL_LR_SCORE_YOUNG) >=  cutoff_score_young)
             & (get(COL_RAW_PVAL_SPECIFICITY_YOUNG) <= cutoff_specificity_young)
         ]
    
    data[, (COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD) := 
             (get(COL_LR_DETECTED_OLD) == TRUE) 
             & (get(COL_LR_SCORE_OLD) >=  cutoff_score_old)
             & (get(COL_RAW_PVAL_SPECIFICITY_OLD) <= cutoff_specificity_old)
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
analyze_significant <- function(data, 
                           cols,
                           cutoff_pval,
                           cutoff_logFC_abs) {
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
add_case_type <- function(data, cols,
                            reassignment = NULL) {
    
    COL_CASE_TYPE = cols$CASE_TYPE
    COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG = cols$LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG
    COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD = cols$LR_DETECTED_AND_SIGNIFICANT_IN_OLD
    COL_DIFFERENTIAL_EXPRESSED = cols$DIFFERENTIAL_EXPRESSED
    COL_DIFFERENTIAL_DIRECTION = cols$DIFFERENTIAL_DIRECTION
    
    data[, (COL_CASE_TYPE) := paste0(
        extract_first_letter(get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG)),
        extract_first_letter(get(COL_LR_DETECTED_AND_SIGNIFICANT_IN_OLD)),
        extract_first_letter(get(COL_DIFFERENTIAL_EXPRESSED)),
        extract_first_letter(get(COL_DIFFERENTIAL_DIRECTION))
        )
    ]
    
    if (!is.null(reassignment)) {
        for (i in 1:length(reassignment)) {
            
            target_val = names(reassignment)[i]
            reassignment_val = reassignment[i]
            
            data[get(COL_CASE_TYPE) == target_val, 
                 (COL_CASE_TYPE) := reassignment_val]
        }
    }

    return(data)
}

#' Creating a list of column names.
#' @return list
#' @export
get_default_colnames <- function() {
    
    cols = list(
        "TISSUES" = "TISSUE",
        "LR_GENES" = "LR_GENES",
        "L_GENE" = "L_GENE",
        "R_GENE" = "R_GENE",
        "LIGAND_CELLTYPE" = "L_CELLTYPE",
        "RECEPTOR_CELLTYPE" = "R_CELLTYPE",
        "LIGAND_RECEPTOR_CELLTYPES" = "LR_CELLTYPES",
        "LR_SCORE_YOUNG" = "LR_SCORE_young",
        "LR_SCORE_OLD" = "LR_SCORE_old",
        "LOGFC" = "LR_LOGFC",
        "LOGFC_ABS" = "LR_LOGFC_ABS",
        "BH_PVAL" = "BH_PVAL_DIFF",
        "RAW_PVAL" = "PVAL_DIFF",
        "LR_DETECTED_YOUNG" = "LR_DETECTED_young",
        "LR_DETECTED_OLD" = "LR_DETECTED_old",
        "LIGAND_EXPRESSION_YOUNG" = "L_EXPRESSION_young",
        "LIGAND_EXPRESSION_OLD" = "L_EXPRESSION_old",
        "LIGAND_DETECTED_YOUNG" = "L_DETECTED_young",
        "LIGAND_DETECTED_OLD" = "L_DETECTED_old",
        "RECEPTOR_EXPRESSION_YOUNG" = "R_EXPRESSION_young",
        "RECEPTOR_EXPRESSION_OLD" = "R_EXPRESSION_old",
        "RECEPTOR_DETECTED_YOUNG" = "R_DETECTED_young",
        "RECEPTOR_DETECTED_OLD" = "R_DETECTED_old",
        "RAW_PVAL_SPECIFICITY_YOUNG" = "BH_PVAL_young",
        "RAW_PVAL_SPECIFICITY_OLD" = "BH_PVAL_old",
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
check_columns <- function(data, cols) {
    message("Checking column names...")
    
    essential_cols = c(
       cols$TISSUES,
       cols$LR_GENES,
       cols$LIGAND_CELLTYPE,
       cols$RECEPTOR_CELLTYPE,
       cols$LR_SCORE_YOUNG,
       cols$LR_SCORE_OLD,
       cols$BH_PVAL,
       cols$RAW_PVAL,
       cols$LR_DETECTED_YOUNG,
       cols$LR_DETECTED_OLD,
       cols$LIGAND_EXPRESSION_YOUNG,
       cols$LIGAND_EXPRESSION_OLD,
       cols$LIGAND_DETECTED_YOUNG,
       cols$LIGAND_DETECTED_OLD,
       cols$RECEPTOR_EXPRESSION_YOUNG,
       cols$RECEPTOR_EXPRESSION_OLD,
       cols$RECEPTOR_DETECTED_YOUNG,
       cols$RECEPTOR_DETECTED_OLD,
       cols$RAW_PVAL_SPECIFICITY_YOUNG,
       cols$RAW_PVAL_SPECIFICITY_OLD
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


# Overrepresentation analysis --------------------------------------------------

#' Title
#'
#' @param data 
#' @param cols 
#'
#' @return
#' @export
#'
#' @examples
analyze_ORA <- function(data, 
                        cols = get_default_colnames()) {
    
    COL_DIFFERENTIAL_DIRECTION = cols$DIFFERENTIAL_DIRECTION
    
    dt_ora = ora(data, cols)
    dt_ora_up = ora(data[get(COL_DIFFERENTIAL_DIRECTION) == "UP"], cols)
    dt_ora_down = ora(data[get(COL_DIFFERENTIAL_DIRECTION) == "DOWN"], cols)
    
    dt_complete = merge(dt_ora, dt_ora_up,
                        by = c("Tissue", "Category", "Value"), 
                        all = TRUE, 
                        suffixes = c("", "_UP"))
    
    dt_complete = merge(dt_complete, dt_ora_down,
                        by = c("Tissue", "Category", "Value"),
                        all = TRUE,
                        suffixes = c("", "_DOWN"))
    
    return(dt_complete)
}


#' Title
#'
#' @param data 
#' @param cols 
#'
#' @return
#' @export
#'
#' @examples
ora <- function(data, cols) {
    
    # tissues = c(unique(data[[cols$TISSUES]]))
    categories = c(
        cols$TISSUES, cols$LIGAND_CELLTYPE, cols$RECEPTOR_CELLTYPE, 
        cols$LIGAND_RECEPTOR_CELLTYPES, cols$L_GENE, cols$R_GENE, cols$LR_GENES
    )
    
    dt_counts = fast_counts(data, cols, categories)
    dt_ora = perform_ora_from_counts(dt_counts)
    
    return(dt_ora)
}


#' Title
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
perform_ora_from_counts <- function(data) {
    
    # https://stackoverflow.com/questions/11680579/assign-multiple-columns-using-in-data-table-by-group
    
    data[, c("OR", "pval") := 
             vfisher_2sided(Counts_value_significant,
                           Counts_value_notsignificant,
                           Counts_notvalue_significant, 
                           Counts_notvalue_notsignificant)
    ]
    
    data[, c("Kulc_distance", "Imbalance_ratio") := list(
                kulc(Counts_value_significant, 
                     Counts_value_notsignificant, 
                     Counts_notvalue_significant, 
                     Counts_notvalue_notsignificant),
                imbalance_ratio(Counts_value_significant, 
                                Counts_value_notsignificant, 
                                Counts_notvalue_significant,
                                Counts_notvalue_notsignificant)
    )]
    
    # Add adjusted pval
    data[, "pval_adjusted" := .(p.adjust(pval, method="BH"))]
    
    return(data)
}

#' Title
#'
#' @param data 
#' @param cols 
#' @param categories 
#'
#' @return
#' @export
#'
#' @examples
fast_counts <- function(data, cols, categories) {
    
    COL_TISSUE = cols$TISSUES
    
    dt_categories = count_on_all_columns(data, cols, categories)
    dt_tissues = count_significant_on_tissues(data, cols)
    
    dt_counts = merge(dt_categories, dt_tissues, by=COL_TISSUE, all=TRUE)
    
    dt_counts = dt_counts[, .(
        Tissue = get(COL_TISSUE),
        Category = Category,
        Value = Value,
        
        Counts_value_significant = Counts_value_significant,
        Counts_value_notsignificant = Counts_value_notsignificant,
        Counts_notvalue_significant = Counts_significant - Counts_value_significant,
        Counts_notvalue_notsignificant = Counts_notsignificant - Counts_value_notsignificant
        )
    ]
    
    return(dt_counts)
}


#' Title
#'
#' @param data 
#' @param cols 
#' @param categories 
#'
#' @return
#' @export
#'
#' @examples
count_on_all_columns <- function(data, cols, categories) {
    
    dt_categories = list()
    
    for (category in categories) {
        dt_cat = count_on_one_column(data, cols, category)
        dt_categories[[category]] = dt_cat
    }
    
    dt_counts = rbindlist(dt_categories, use.names = TRUE)
    return(dt_counts)
}


#' Title
#'
#' @param data 
#' @param cols 
#' @param column 
#'
#' @return
#' @export
#'
#' @examples
count_on_one_column <- function(data, cols, column) {
    
    COL_TISSUE = cols$TISSUES
    COL_DIFFERENTIAL_EXPRESSED = cols$DIFFERENTIAL_EXPRESSED
    COL_DIFFERENTIAL_DIRECTION = cols$DIFFERENTIAL_DIRECTION
    
    dt_val_sig = data[get(COL_DIFFERENTIAL_EXPRESSED) == TRUE,
                      .(Category=column, Counts_value_significant=.N),
                      by=.(get(COL_TISSUE), Value=get(column))]
    
    dt_val_notsig = data[get(COL_DIFFERENTIAL_EXPRESSED) == FALSE,
                         .(Category=column, Counts_value_notsignificant=.N),
                         by=.(get(COL_TISSUE), Value=get(column))]
    
    dt_list = list(dt_val_sig, dt_val_notsig)
    
    lapply(dt_list, function(dt) {
        setnames(dt, new = COL_TISSUE, old = "get")
    })
    

    dt_cat_final = merge(dt_list[[1]], dt_list[[2]], 
                         by=c(COL_TISSUE, "Value", "Category"),
                         all=TRUE)

    setnafill(dt_cat_final, 
              "const",
              fill=0,
              cols=4:5)
    
    # Add entries for all organs
    dt_all = dt_cat_final[, lapply(.SD, sum), by=.(Value, Category), .SDcols = !c(COL_TISSUE)]
    dt_all[, (COL_TISSUE) := "All"]
    dt_cat_final = rbindlist(list(dt_cat_final, dt_all), use.names=TRUE)
    
    return(dt_cat_final)
    
}

#' Title
#'
#' @param data 
#' @param cols 
#'
#' @return
#' @export
#'
#' @examples
count_significant_on_tissues <- function(data, cols) {
    
    COL_TISSUE = cols$TISSUES
    COL_DIFFERENTIAL_EXPRESSED = cols$DIFFERENTIAL_EXPRESSED
    COL_DIFFERENTIAL_DIRECTION = cols$DIFFERENTIAL_DIRECTION
    
    # To be merged at last on tissue level + to create "All" entry
    dt_sig = data[get(COL_DIFFERENTIAL_EXPRESSED) == TRUE,
                  .(Counts_significant = .N),
                  by = .(get(COL_TISSUE))]
    dt_notsig = data[get(COL_DIFFERENTIAL_EXPRESSED) == FALSE,
                     .(Counts_notsignificant = .N),
                     by = .(get(COL_TISSUE))]
    
    dt_list = list(dt_sig, dt_notsig)
    
    lapply(dt_list, function(dt) {
        setnames(dt, new = COL_TISSUE, old = "get")
    })
    
    dt_tissues_counts = merge(dt_list[[1]], dt_list[[2]], 
                              by=c(COL_TISSUE),
                              all=TRUE)
    
    setnafill(dt_tissues_counts, "const", fill=0, cols=2:3)
    
    dt_all = dt_tissues_counts[, lapply(.SD, sum), .SDcols = !c(COL_TISSUE)]
    dt_all[, (COL_TISSUE) := "All"]
    dt_tissues_counts = rbindlist(list(dt_tissues_counts, dt_all), use.names=TRUE)
    
    return(dt_tissues_counts)
}


#' Title
#'
#' @param Counts_value_significant 
#' @param Counts_value_notsignificant 
#' @param Counts_notvalue_significant 
#' @param Counts_notvalue_notsignificant 
#'
#' @return
#' @export
#'
#' @examples
odds_ratio <- function(Counts_value_significant,
                       Counts_value_notsignificant,
                       Counts_notvalue_significant,
                       Counts_notvalue_notsignificant) {
    
    # There are functions that have better methods for estimating
    #  odds ratio rather than the simple formula below.
    # uses only basic arithmetic operations which are already vectorized
    # or = ((Counts_value_significant * Counts_notvalue_notsignificant)
    #       / (Counts_value_notsignificant * Counts_notvalue_significant))
    
    or = 
    
    return(or)
}


#' Title
#'
#' @param Counts_value_significant 
#' @param Counts_value_notsignificant 
#' @param Counts_notvalue_significant 
#' @param Counts_notvalue_notsignificant 
#'
#' @return
#' @export
#'
#' @examples
fisher_2sided <- function(Counts_value_significant,
                               Counts_value_notsignificant,
                               Counts_notvalue_significant,
                               Counts_notvalue_notsignificant) {
    
    m = matrix(nrow = 2, ncol = 2, byrow=TRUE,
               data = c(Counts_value_significant, Counts_value_notsignificant,
                        Counts_notvalue_significant, Counts_notvalue_notsignificant)
    )
    
    test = fisher.test(m, alternative="two.sided")
    res = c(test$estimate, test$p.value)
    return(res)
}

#' Title
#'
#' @param Counts_value_significant 
#' @param Counts_value_notsignificant 
#' @param Counts_notvalue_significant 
#' @param Counts_notvalue_notsignificant 
#'
#' @return
#' @export
#'
#' @examples
vfisher_2sided <- function(Counts_value_significant,
                           Counts_value_notsignificant,
                           Counts_notvalue_significant,
                           Counts_notvalue_notsignificant) {

    v = mapply(fisher_2sided, 
               Counts_value_significant,
               Counts_value_notsignificant,
               Counts_notvalue_significant,
               Counts_notvalue_notsignificant)
    
    l = list(OR = v[1, ], pval = v[2, ])
    return(l)
}

#' Title
#'
#' @param Counts_value_significant 
#' @param Counts_value_notsignificant 
#' @param Counts_notvalue_significant 
#' @param Counts_notvalue_notsignificant 
#'
#' @return
#' @export
#'
#' @examples
kulc <- function(Counts_value_significant,
                 Counts_value_notsignificant,
                 Counts_notvalue_significant,
                 Counts_notvalue_notsignificant) {
    
    # vectorized by using only arithmetic operations
    P_AB = Counts_value_significant / (Counts_value_significant + Counts_notvalue_significant)
    P_BA = Counts_value_significant / (Counts_value_significant + Counts_value_notsignificant)
    avg = (P_AB + P_BA) / 2
    return(avg)
}


#' Title
#'
#' @param Counts_value_significant 
#' @param Counts_value_notsignificant 
#' @param Counts_notvalue_significant 
#' @param Counts_notvalue_notsignificant 
#'
#' @return
#' @export
#'
#' @examples
imbalance_ratio <- function(Counts_value_significant,
                            Counts_value_notsignificant,
                            Counts_notvalue_significant,
                            Counts_notvalue_notsignificant) {
    # Vectorized
    numerator = abs(Counts_value_notsignificant - Counts_notvalue_significant)
    denominator = (Counts_value_significant 
                   + Counts_value_notsignificant
                   + Counts_notvalue_significant)
    return(numerator / denominator)
    
}


# Frequent itemsets analysis ---------------------------------------------------

#' analyze_FreqItemSets
#'
#' @param data 
#' @param cols 
#' @param items_to_include 
#' @param target 
#' @param support 
#' @param confidence 
#'
#' @return
#' @export
#'
#' @examples
analyze_FreqItemSets <- function(data, cols, 
                                 items_to_include = NULL,
                                 target = "closed frequent itemsets",
                                 support = 0.05,
                                 confidence = 0.1) {
    
    transactions = convert_data_to_transactions(data, cols)
    
    res = compute_freqitemsets_and_rules(transactions, 
                                         target,
                                         support,
                                         confidence)
    
    freqsets = res[["freqsets"]]
    rules = res[["rules"]]
    
    interesting_rules = get_interesting_rules(rules, transactions, cols)
    
    return(interesting_rules)
}

#' get_interesting_rules
#'
#' @param rules 
#' @param transactions 
#'
#' @return
#' @export
#'
#' @examples
get_interesting_rules <- function(rules, transactions, cols) {
    
    STRING_sign_diff_expression = glue("{s}=Sign", s=cols$DIFFERENTIAL_EXPRESSED)
    STRING_tissue = glue("{s}=", s=cols$TISSUES)
    STRING_ligand_gene = glue("{s}=", s=cols$L_GENE)
    STRING_receptor_gene = glue("{s}=", s=cols$R_GENE)
    STRING_ligand_celltype = glue("{s}=", s=cols$LIGAND_CELLTYPE)
    STRING_receptor_celltype = glue("{s}=", s=cols$RECEPTOR_CELLTYPE)
    
    sub1 = subset(
        rules, 
        subset = rhs %in% STRING_sign_diff_expression
        & lhs %pin% STRING_tissue
        & lhs %pin% STRING_ligand_gene
        & lhs %pin% STRING_receptor_celltype
        # & !(lhs %pin% "Direction=")
    )
    sub1_measure = interestMeasure(sub1,
                                   measure=c("fishersExactTest", "oddsRatio"),
                                   transactions=transactions,
                                   reuse=TRUE)
    sub1_dt = data.table(
        lhs=labels(lhs(sub1)),
        rhs=labels(rhs(sub1)), 
        sub1@quality,
        sub1_measure
    )
    
    
    sub2 = subset(
        rules, 
        subset = rhs %in% STRING_sign_diff_expression
        & lhs %pin% STRING_tissue
        & lhs %pin% STRING_receptor_gene
        & lhs %pin% STRING_ligand_celltype
        # & !(lhs %pin% "Direction=")
    )
    sub2_measure = interestMeasure(sub2,
                                   measure=c("fishersExactTest", "oddsRatio"),
                                   transactions=transactions,
                                   reuse=TRUE)
    sub2_dt = data.table(
        lhs=labels(lhs(sub2)),
        rhs=labels(rhs(sub2)), 
        sub2@quality,
        sub2_measure
    )
    
    subsets = list()
    subsets[["sub1"]] = sub1_dt
    subsets[["sub2"]] = sub2_dt
    
    return(subsets)
}

#' compute_freqitemsets_and_rules
#'
#' @param transactions 
#' @param support 
#' @param confidence 
#'
#' @return
#' @export
#'
#' @examples
compute_freqitemsets_and_rules <- function(transactions,
                                           target,
                                           support,
                                           confidence) {
        
    freqsets = apriori(
        transactions, 
        parameter = list(support = support, 
                        minlen = 1,
                        maxlen = 5,
                        target = target, 
                        ext = TRUE,
                        smax = 1,
                        minval = 0, 
                        maxtime = 0),
        control = list(verbose=FALSE)
    )
        
    rules = ruleInduction(
        freqsets, 
        transactions = transactions,
        confidence = confidence,
        control = list(method = "ptree", verbose = FALSE)
    )
        
    return(list(freqsets = freqsets, rules = rules))
}
    
#' convert_data_to_transactions
#'
#' @param data 
#' @param cols 
#'
#' @return
#' @export
#'
#' @examples
convert_data_to_transactions <- function(data, cols) {
    
    COL_DIFFERENTIAL_EXPRESSED = cols$DIFFERENTIAL_EXPRESSED
    
    cols_for_items = c(
        cols$TISSUES, 
        cols$LIGAND_CELLTYPE, cols$RECEPTOR_CELLTYPE, 
        cols$L_GENE, cols$R_GENE, 
        cols$DIFFERENTIAL_EXPRESSED, cols$DIFFERENTIAL_DIRECTION
    )
    
    dt_p = data[, ..cols_for_items]
    
    dt_p[, (COL_DIFFERENTIAL_EXPRESSED) := 
             fifelse(get(COL_DIFFERENTIAL_EXPRESSED), "Sign", "Notsign")]
    
    
    transactions = as(dt_p, "transactions")
    return(transactions)
}


# Plots ------------------------------------------------------------------------

plot_heatmaps <- function(data, cols) {
    
    
    
}

build_heatmaps <- function(dt, cols, dir_path) {
    
    # add "/" if not last char of dir_path
    if (!substr(dir_path, length(dir_path), length(dir_path)) == "/") {
        dir = paste0(dir_path, "/")
    } else {
        dir = dir_path
    }
    
    dir.create(dir, showWarnings = FALSE) 
    
    OR_max = 10
    OR_min = 0.1
    tissues = unique(dt[["Tissue"]])
    ncols = dim(dt)[2]
    
    for (tissue in tissues) {
        
        dt_edge = dt[Tissue == tissue, 2:ncols]
        nrows = dim(dt_edge)[1]
        if (nrows == 1) {message(paste0("One edge for: ", tissue))}
        
        G = graph_from_data_frame(d = dt_edge,
                                  directed = TRUE,
                                  vertices = NULL)
        
        matrix_all = as_adjacency_matrix(G, attr = "OR", sparse = FALSE)
        matrix_all = clip_matrix(matrix_all, 1, OR_min, OR_max)
        
        matrix_up = as_adjacency_matrix(G, attr = "OR_UP", sparse = FALSE)
        matrix_up = clip_matrix(matrix_up, 1, OR_min, OR_max)
        
        matrix_down = as_adjacency_matrix(G, attr = "OR_DOWN", sparse = FALSE)
        matrix_down = clip_matrix(matrix_down, 1, OR_min, OR_max)
        
        dir.create(paste0(dir, tissue, '/'), showWarnings = FALSE)
        
        explanatory_string = paste0("Row-transmitter, col-receiver; ",
                                    "Values are clipped log(OR)")
        eps = 10**(-5)
        generate_heatmap(log(matrix_all + eps),
                         #matrix_all,
                         title = paste0(tissue, ": ", "ALL ", explanatory_string),
                         filename = paste0(dir, tissue, '/', 'ALL.png'))
        generate_heatmap(log(matrix_up + eps), 
                         #matrix_up,
                         title = paste0(tissue, ": ", "UP ", explanatory_string),
                         filename = paste0(dir, tissue, '/', 'UP.png'))
        generate_heatmap(log(matrix_down + eps),
                         #matrix_down,
                         title = paste0(tissue, ": ", "DOWN ", explanatory_string),
                         filename = paste0(dir, tissue, '/', 'DOWN.png'))
    }
}



clip_matrix <- function(matrix, zero_val, min_val, max_val) {
    matrix[matrix == 0] = zero_val
    matrix[matrix < min_val] = min_val
    matrix[matrix > max_val] = max_val
    return(matrix)
}


generate_heatmap <- function(matrix, title, filename) {
    
    message(paste0("Preparing figure: ", title, " ", filename))
    
    nrows = dim(matrix)[1]
    ncols = dim(matrix)[2]
    
    if (nrows < 2 | ncols < 2) {
        message(paste0("Matrix size is < 2x2: ", title))
        return()
    }
    
    breaks = c(seq(-3, -0.5, 0.5), 0, seq(0.5, 3, 0.5))
    pheatmap(matrix, #cellwidth=10, cellheight=10,
             width=11, height=11,
             color = colorRampPalette(brewer.pal(8, "RdYlBu"))(length(breaks)-1),
             breaks = breaks,
             na_col = "green",
             scale = "none",
             cluster_rows = nrows >= 2,
             cluster_cols = ncols >= 2,
             cutree_rows = ifelse(nrows >= 5, min(nrows, 3), 1),
             cutree_cols = ifelse(ncols >= 5, min(ncols, 3), 1),
             display_numbers=FALSE, 
             fontsize=10, 
             main = title,
             ylab = "Transmitter cell", 
             xlab = "Receiver cell",
             legend=TRUE,
             filename = filename)
}


get_celltypes_enrichment <- function(dt, cols) {
    
    COL_LIGAND_RECEPTOR_CELLTYPES = cols$LIGAND_RECEPTOR_CELLTYPES
    
    dt_ctypes = dt[Category == COL_LIGAND_RECEPTOR_CELLTYPES
                   & pval < 0.05]
    dt_ctypes = dt_ctypes[, c("Ligand_cell", "Receptor_cell") := list(
        sub("(.*) : ", "", Value),
        sub(" : (.*)", "", Value)
    )
    ]
    
    cols_to_select = c(
        "Tissue", "Ligand_cell", "Receptor_cell", "OR", "OR_UP", "OR_DOWN"
    )
    
    dt_ctypes = dt_ctypes[, ..cols_to_select]
    
    return(dt_ctypes)
}






