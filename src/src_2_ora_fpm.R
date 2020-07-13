
load_globals_from_config <- function(R_CONFIG_ACTIVE) {
  
  Sys.setenv(R_CONFIG_ACTIVE = R_CONFIG_ACTIVE)
  
  #assign("DATA_PATH", config::get("DATA_PATH"), envir=.GlobalEnv)
  #assign("RESULTS_DIR", config::get("RESULTS_DIR"), envir=.GlobalEnv)
  # CUTOFFS
  #assign("CUTOFF_PVAL_ADJ", config::get("CUTOFF_PVAL_ADJ"), envir=.GlobalEnv)
  #assign("CUTOFF_SCORE_YOUNG", config::get("CUTOFF_SCORE_YOUNG"), envir=.GlobalEnv)
  #assign("CUTOFF_SCORE_OLD", config::get("CUTOFF_SCORE_OLD"), envir=.GlobalEnv)
  #assign("CUTOFF_LOGFC", config::get("CUTOFF_LOGFC"), envir=.GlobalEnv)
  #assign("CUTOFF_CPDB_PVAL_YOUNG", config::get("CUTOFF_CPDB_PVAL_YOUNG"), envir=.GlobalEnv)
  #assign("CUTOFF_CPDB_PVAL_OLD", config::get("CUTOFF_CPDB_PVAL_OLD"), envir=.GlobalEnv)
  # COLUMNS' NAMES
  assign("COL_TISSUES", config::get("COL_TISSUES"), envir=.GlobalEnv)
  assign("COL_LR_GENES", config::get("COL_LR_GENES"), envir=.GlobalEnv)
  assign("COL_L_GENE", config::get("COL_L_GENE"), envir=.GlobalEnv)
  assign("COL_R_GENE", config::get("COL_R_GENE"), envir=.GlobalEnv)
  assign("COL_LIGAND_CELLTYPE", config::get("COL_LIGAND_CELLTYPE"), envir=.GlobalEnv)
  assign("COL_RECEPTOR_CELLTYPE", config::get("COL_RECEPTOR_CELLTYPE"), envir=.GlobalEnv)
  assign("COL_LIGAND_RECEPTOR_CELLTYPES", config::get("COL_LIGAND_RECEPTOR_CELLTYPES"), envir=.GlobalEnv)
  assign("COL_LR_SCORE_YOUNG", config::get("COL_LR_SCORE_YOUNG"), envir=.GlobalEnv)
  assign("COL_LR_SCORE_OLD", config::get("COL_LR_SCORE_OLD"), envir=.GlobalEnv)
  assign("COL_LOGFC", config::get("COL_LOGFC"), envir=.GlobalEnv)
  assign("COL_LOGFC_ABS", config::get("COL_LOGFC_ABS"), envir=.GlobalEnv)
  assign("COL_BH_PVAL", config::get("COL_BH_PVAL"), envir=.GlobalEnv)
  assign("COL_RAW_PVAL", config::get("COL_RAW_PVAL"), envir=.GlobalEnv)
  assign("COL_LR_DETECTED_YOUNG", config::get("COL_LR_DETECTED_YOUNG"), envir=.GlobalEnv)
  assign("COL_LR_DETECTED_OLD", config::get("COL_LR_DETECTED_OLD"), envir=.GlobalEnv)
  #assign("COL_LIGAND_EXPRESSION_YOUNG", config::get("COL_LIGAND_EXPRESSION_YOUNG"), envir=.GlobalEnv)
  #assign("COL_LIGAND_EXPRESSION_OLD", config::get("COL_LIGAND_EXPRESSION_OLD"), envir=.GlobalEnv)
  #assign("COL_LIGAND_DETECTED_YOUNG", config::get("COL_LIGAND_DETECTED_YOUNG"), envir=.GlobalEnv)
  #assign("COL_LIGAND_DETECTED_OLD", config::get("COL_LIGAND_DETECTED_OLD"), envir=.GlobalEnv)
  #assign("COL_RECEPTOR_EXPRESSION_YOUNG", config::get("Receptor_expr_young"), envir=.GlobalEnv)
  #assign("COL_RECEPTOR_EXPRESSION_OLD", config::get("COL_RECEPTOR_EXPRESSION_OLD"), envir=.GlobalEnv)
  #assign("COL_RECEPTOR_DETECTED_YOUNG", config::get("COL_RECEPTOR_DETECTED_YOUNG"), envir=.GlobalEnv)
  #assign("COL_RECEPTOR_DETECTED_OLD", config::get("COL_RECEPTOR_DETECTED_OLD"), envir=.GlobalEnv)
  assign("COL_BH_PVAL_SPECIFICITY_YOUNG", config::get("COL_BH_PVAL_SPECIFICITY_YOUNG"), envir=.GlobalEnv)
  assign("COL_BH_PVAL_SPECIFICITY_OLD", config::get("COL_BH_PVAL_SPECIFICITY_OLD"), envir=.GlobalEnv)
}

# funcs to perform overrepresentation analysis


build_contingency_table <- function(df, df_filt, col, value) {
  
  # Compute overrepresentation of the value from col in df_filt
  #  in comparison to df. 
  #  E.g. Overrepresentation of value "Kidney" from column "Tissue".
  
  have_same_cols = all(colnames(df) == colnames(df_filt))
  have_col = (col %in% colnames(df))
  
  if ( !(have_same_cols & have_col) ) {
    stop("Invalid input")
  }
  
  # Count
  counts_filtered_interest = sum(df_filt[, col] == value)
  counts_filtered_not_interest = sum(df_filt[, col] != value)
  counts_not_filtered_interest = sum(df[, col] == value) - counts_filtered_interest
  counts_not_filtered_not_interest = sum(df[, col] != value) - counts_filtered_not_interest
  
  contingency_df = data.frame(value.not.interest=c(counts_filtered_not_interest, 
                                                   counts_not_filtered_not_interest), 
                              value.of.interest=c(counts_filtered_interest, 
                                                  counts_not_filtered_interest))
  colnames(contingency_df) <- c(paste("not_", value, sep=""), value)
  row.names(contingency_df) <- c("filtered", "not_in_filtered")
  
  return(contingency_df)
}


compute_overrepr_pval <- function(df) {
  
  # Compute overrepresentation from a contingency table
  # df has to be 2x2 and is assumed to follow structure:
  #   - columns: (not_in_filtered_df, in_filtered_df)
  #   - rows: (in_category, not_in_category), e.g. counts_tissue_interest and counts_not_tissue_interest
  
  if ( !(all(dim(df) == c(2,2))) ) {
    stop("Invalid contingency table dimension.")
  }
  
  # Alternative hypothesis is "less" due to construction of the contingency table.
  #  The H1 is OR < 1, because we search for features where the odds that
  #  a value of interest is in filtered is > odds that a value NOT of interest
  #  is in filtered. 
  #  E.g. Odds(Not kidney is in filtered) < Odds(Kidney is in filtered) => OR < 1.
  test = fisher.test(df, alternative="less")
  return(test$p.value)
}


compute_column_overrepr <- function(df, df_filt, col, adjust_pvals=FALSE) {
  
  # Computes overrepresentation of all terms from
  #  a selected column.
  
  is_col_char = (class(col) == "character")
  is_single_element = (length(col) == 1)
  
  if ( !(is_col_char & is_single_element) ) {
    stop("Invalid column as input. Has to be character vector of length 1.")
  }
  
  unique_values = unique(df[, col])
  
  overrepr_df = data.frame()
  for (val in unique_values) {
    contingency = build_contingency_table(df, df_filt, col, val)
    pval = compute_overrepr_pval(contingency)
    
    if (adjust_pvals) {
      adjpval = p.adjust(pval, method='BH')
      row = data.frame(overrepr_pval=pval,
                       overrepr_pval_adj=adjpval,
                       in_filtered_of_interest=contingency[1,2],
                       in_filtered_not_interest=contingency[1,1], 
                       not_in_filtered_of_interest=contingency[2,2],
                       not_in_filtered_not_interest=contingency[2,1],
                       row.names=val)
    } else {
      row = data.frame(overrepr_pval=pval,
                       in_filtered_of_interest=contingency[1,2],
                       in_filtered_not_interest=contingency[1,1], 
                       not_in_filtered_of_interest=contingency[2,2],
                       not_in_filtered_not_interest=contingency[2,1],
                       row.names=val)
    }
    
    overrepr_df = rbind(overrepr_df, row)
  }
  
  # add column explicitely, reorder row index
  overrepr_df[, col] = rownames(overrepr_df)
  rownames(overrepr_df) = NULL
  
  
  # Sort if not empty
  
  if ( dim(overrepr_df)[1] == 0) {
    return(overrepr_df)
  } else {
    overrepr_df = overrepr_df[order(overrepr_df[, "overrepr_pval"]), ]
    return(overrepr_df)
  }
}


analyze_overrepresentation <- function(df, df_filt, exclude_tissue_overrepresentation=FALSE, adjust_pvals=FALSE) {
  
  # Perform overrepresentation analysis over predefined
  #  columns.
  
  # Replaced by global column names variables such as COL_TISSUE  
  #  cols_to_analyze = c("Tissue", "LR_genes", "L_gene", "R_gene",
  #                      "Ligand_cell_type", "Receptor_cell_type",
  #                      "Ligand_Receptor_cell_types")
  
  col_pval = "overrepr_pval"
  cols_to_analyze = c(COL_TISSUES, COL_LR_GENES, COL_L_GENE, COL_R_GENE,
                      COL_LIGAND_CELLTYPE, COL_RECEPTOR_CELLTYPE, 
                      COL_LIGAND_RECEPTOR_CELLTYPES)
  
  if (exclude_tissue_overrepresentation) {
    cols_to_analyze = cols_to_analyze[2:length(cols_to_analyze)]
  }
  
  df_up = df[df[, COL_LOGFC] > 0, ]
  df_filt_up = df_filt[df_filt[, COL_LOGFC] > 0, ]
  
  df_down = df[df[, COL_LOGFC] < 0, ]
  df_filt_down = df_filt[df_filt[, COL_LOGFC] < 0, ]
  
  overrepr_dfs = list()
  for (col in cols_to_analyze) {
    print(col)
    
    overrepr_all = compute_column_overrepr(df, df_filt, col, adjust_pvals)
    overrepr_up = compute_column_overrepr(df_up, df_filt_up, col, adjust_pvals)
    overrepr_down = compute_column_overrepr(df_down, df_filt_down, col, adjust_pvals)
    
    merged = merge(overrepr_all, overrepr_up, by=col, all=TRUE, suffixes=c("_ALL", ""))  # outer join
    merged = merge(merged, overrepr_down, by=col, all=TRUE, suffixes=c("_UP", "_DOWN"))
    #    colnames(merged) = c( col, "ORA_all_interactions", 
    #                         "ORA_up_interactions", "ORA_down_interactions")
    
    overrepr_dfs[[col]] = merged
  }
  
  return(overrepr_dfs)
}

