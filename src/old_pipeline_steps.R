filter_detected <- function(df) {
  
  # Filter for getting detected interactions.
  #  It also wrapps each filter to account for NULL cutoffs.
  
  res = remove_undetected_interactions(df)
  
  
  mask <- ((df[, COL_LR_SCORE_YOUNG] >= CUTOFF_SCORE_YOUNG) | (df[, COL_LR_SCORE_OLD] >= CUTOFF_SCORE_OLD))
  res <- df[mask, ]
  
  
  
  
  if ( !(is.null(CUTOFF_CPDB_PVAL_YOUNG)) & !(is.null(CUTOFF_CPDB_PVAL_OLD)) ) {
    res = remove_nonspecific_interactions(res, CUTOFF_CPDB_PVAL_YOUNG, CUTOFF_CPDB_PVAL_OLD)
  }
  
  return(res)
}


filter_significant <- function(df) {
  
  # Composite filter on the results dataframe.
  #  It also wrapps each filter to account for NULL cutoffs.
  
  if ( !(is.null(CUTOFF_PVAL_ADJ)) ) {
    res = remove_unsign_adj_pval(df, CUTOFF_PVAL_ADJ) 
  } else {
    res = df
  }
  
  if ( !(is.null(CUTOFF_LOGFC)) ) {
    res = remove_low_logFC(res, CUTOFF_LOGFC) 
  }

  return(res)
}


run_filtering <- function(df, filter, name, dir) {
  
  # Runs filtering: applyies filter, summarizes and saves summaries.
  
  #message = paste0("Running filtering: ", name, " in ", dir)
  #print(message)
  
  df_filt = filter(df)
  
  res_dir = paste(dir, "/", name, "/", sep='')
  create_dir(res_dir)
  
  summ_filt = summarize_filtering(df, df_filt)
  save_summary(summ_filt, "Removed_summary", res_dir)
  
  summ_res = summarize_dataframe(df_filt)
  save_summary(summ_res, "Output_dataframe_summary", res_dir)
  
  return(df_filt)
}


analyze_overrepr_per_tissue <- function(df_detected, df_significant, dir_, adjust_pvals=FALSE) {
  
  tissues = unique(df_detected[, COL_TISSUES])
  for (tissue in tissues) {
    
    print(paste0("Computing overrepresentation for: ", tissue))
    
    df_tiss_detected = filter_tissues(df_detected, tissue)
    df_tiss_significant = filter_tissues(df_significant, tissue)
    
    overrepr_tiss = analyze_overrepresentation(df_tiss_detected, df_tiss_significant,
                                               exclude_tissue_overrepresentation=TRUE,
                                               adjust_pvals=adjust_pvals)
    # overrepr_tiss = adjust_pvals_of_list_dfs(overrepr_tiss, col_pval="overrepr_pval")
    
    dir_res_over_tissue = paste0(dir_, tissue, '/')
    save_list_dfs(overrepr_tiss, dir_res_over_tissue)
  }
}


adjust_pvals_of_list_dfs <- function(dfs, col_pval) {
  # Operates on list of dataframes dfs and adjust all the
  #  pvals from the selected col_pval column.
  
  stop("Not implemented")
}