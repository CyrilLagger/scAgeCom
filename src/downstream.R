# Primitive procedures for downstream analysis of intercellular communications
#  analysis results. Procedures operate on the intercellular communications
#  dataframe, which comes out directly from the analysis pipeline. Functions have
#  to preserve the results representation invariant, defined in a separate function below.


###############################################################
                    # METHODS FILTERING #
###############################################################


# preprocess_input_df <- function(df) {
#   
#   # Adds some columns for convenience
#   
#   # Add separate columns for ligand and receptor gene name
#   LR_genes = df[, COL_LR_GENES]
#   splitter = function(string) strsplit(string, "_")[[1]]
#   LR_genes_split = t(sapply(LR_genes, splitter, USE.NAMES=FALSE))
#   df[, COL_L_GENE] = LR_genes_split[, 1]
#   df[, COL_R_GENE] = LR_genes_split[, 2]
#   
#   # Add column with joined ligand to receptor cell type
#   df[, COL_LIGAND_RECEPTOR_CELLTYPES] = paste(
#     df[, COL_LIGAND_CELLTYPE], 
#     df[, COL_RECEPTOR_CELLTYPE], 
#     sep=" : ")
#   
#   # Add logFC
#   df[, COL_LOGFC] = log2(df[, COL_LR_SCORE_OLD] / (df[, COL_LR_SCORE_YOUNG] + 1))
#   df[, COL_LOGFC_ABS] = abs( df[, COL_LOGFC] )
#   
#   return(df)
# }


check_representation_inv <- function(df) {
  
  # Checks representation invariant of the dataframe
  #   - class is data.table/data.frame
  #   - includes the required columns
  
  cols = c(COL_TISSUES, COL_LR_GENES, COL_L_GENE, COL_R_GENE,
           COL_LIGAND_CELLTYPE, COL_RECEPTOR_CELLTYPE,  COL_LIGAND_RECEPTOR_CELLTYPES,
           COL_LR_SCORE_YOUNG, COL_LR_SCORE_OLD, COL_BH_PVAL, COL_RAW_PVAL,
           COL_LR_DETECTED_YOUNG, COL_LR_DETECTED_OLD,
           #COL_LIGAND_EXPRESSION_YOUNG, COL_LIGAND_EXPRESSION_OLD,
           #COL_LIGAND_DETECTED_YOUNG, COL_LIGAND_DETECTED_OLD, COL_RECEPTOR_EXPRESSION_YOUNG,
           #COL_RECEPTOR_EXPRESSION_OLD, COL_RECEPTOR_DETECTED_YOUNG, COL_RECEPTOR_DETECTED_OLD,
           COL_BH_PVAL_SPECIFICITY_YOUNG, COL_BH_PVAL_SPECIFICITY_OLD)

  
  is_correct_type = all(class(df) == c("data.table", "data.frame")) | 
    all(class(df) == c("data.frame"))
  has_reqd_cols = all(cols %in% colnames(df))
  rep_inv = (is_correct_type & has_reqd_cols)
  
  if (!rep_inv) { 
    stop("Representation invariant violated.") 
  }
  
  return(rep_inv)
}

remove_undetected_interactions <- function(df) {
  
  # Removes undetected interactions, defined as those for which expression of
  #  BOTH ligand and receptor is not detected.
  
  check_representation_inv(df)
  
  col_is_young_detected = COL_LR_DETECTED_YOUNG
  col_is_old_detected = COL_LR_DETECTED_OLD
  
  mask = (df[, col_is_old_detected] != FALSE) | (df[, col_is_young_detected] != FALSE)
  res = df[mask, ]
  
  check_representation_inv(res)
  return(res)
}

explore_filter_cutoffs <- function(df, dir) {
  
  # to_explore
  colnames = c(COL_LR_SCORE_YOUNG, COL_LR_SCORE_OLD, COL_LOGFC_ABS,
               COL_BH_PVAL_SPECIFICITY_YOUNG, COL_BH_PVAL_SPECIFICITY_OLD)
  df_cols = df[, colnames]
  
  # Build quantile df
  probs = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
  apply_quantile = function(vec) { quantile(vec, probs) }
  quantiles = t( sapply(df_cols, apply_quantile) )
  write.table(quantiles, file=paste0(dir, "quantiles.tsv"), sep='\t',
              quote=FALSE, row.names=TRUE, col.names=TRUE)
  
  # Build the correlogram
  correlogram = cor(df_cols, method="spearman")
  # to use pearson, need to convert INF2NA:
  #  function(x) { x[is.infinite(x)] <- NA; x }
  write.table(correlogram, file=paste0(dir, "correlogram.tsv"), sep='\t',
              quote=FALSE, row.names=TRUE, col.names=TRUE)
  # Dropped building the samples
  # for all combination of selected quantiles
  return(list(quantiles = quantiles,
              correlogram = correlogram))
}

##

remove_unsign_adj_pval <- function(df, cutoff) {
  
  # Filters dataframe by applying cutoff on adjusted p-value.
  
  check_representation_inv(df)
  
  adj_pval_col = COL_BH_PVAL
  res = df[df[, adj_pval_col] <= cutoff, ]
  
  check_representation_inv(res)
  return(res)
}





remove_lowScore <- function(df, age, cutoff) {
  
  # Remove LR interactions of low intensity.
  # Args:
  #     age: str in {"old", "young"}
  #     cutoff: float in (0, inf)
  
  if (!(cutoff>0)) {
    stop("Invalid input: cutoff is float > 0.")
  }
  
  check_representation_inv(df)
  
  if (age=='old') {
    col = COL_LR_SCORE_OLD
    
  } else if (age=='young') {
    col = COL_LR_SCORE_YOUNG
    
  } else {
    stop("Invalid input: age is str in {'old', 'young'}")
  }
  
  mask = df[, col] > cutoff
  res = df[mask, ]
  
  check_representation_inv(res)
  return(res)
}


remove_low_logFC <- function(df, cutoff) {
  
  # Removes interactions with low logFC between old and 
  #  young samples.
  # Low logFC is defined as abs(logFC) < cutoff
  
  check_representation_inv(df)
  
  col_lr_old = COL_LR_SCORE_OLD
  col_lr_young = COL_LR_SCORE_YOUNG
  
  logFC = log2(df[, col_lr_old] / df[, col_lr_young])
  mask = (logFC <= cutoff)
  res = df[mask, ]
  
  check_representation_inv(res)
  return(res)
}

remove_nonspecific_interactions <- function(df, pval_cutoff_young, pval_cutoff_old) {
  # Removes nonspecific interactions that are detected as those with a
  #  p-value from CPDB > pval_cutoff.
  
  col_pval_spec_young = COL_RAW_PVAL_SPECIFICITY_YOUNG
  col_pval_spec_old = COL_RAW_PVAL_SPECIFICITY_OLD
  
  mask = ( (df[, col_pval_spec_young] <= pval_cutoff_young)
           | (df[, col_pval_spec_old] <= pval_cutoff_old) )

  res = df[mask, ]
  return(res)
}


#build_monad_filter <- function(filter, ...) {
#  # Build filter that operates only on the dataframe,
#  #  hard-setting the arguments. Useful for composing
#  #  filters into a pipeline.
#  
#  monad <- function(df) {
#    return(filter(df, ...))
#  }
#  return(monad)
#}


### A nice way of composing more complex filters and pipelines, but
###  suffers from the limits of the low C stack size allocated by default
###  to R. Can be configured, but may introduce platform-specific issues,
###  so currently abandoned.

#compose_filters <- function(...) {
#  
#  # Composes sequentially the input filters and builds
#  #  a function for the composite filter.
#  
#  stop("NotImplemented. Exceeds standard allocated stack of 8Mb.")
#  
#  filters <- c(...)
#  base_filter <- function(df) { return(df) }
#  
#  composite_filter <- base_filter
#  
#  for (f in filters) {
#    composite_filter <- function(df) {
#      return(f(composite_filter(df)))
#    }
#  }
#  
#  return(composite_filter)
#}


###############################################################
            # METHODS FOR SUMMARIZING THE DATAFRAME #
###############################################################


filter_tissues <- function(df, tissues) {
  
  # Returns dataframe with entries corresponding to tissues.
  
  col_tissue = COL_TISSUES
  if (!all(tissues %in% df[, col_tissue]) | (class(tissues) != "character")) {
    stop("Invalid input: tissues not found in dataframe OR are not 'character' vector.")
  }
  
  mask = df[, col_tissue] %in% tissues
  res = df[mask, ]
  return(res)
}


filter_celltypes <- function(df, ctypes_ligand=NULL, ctypes_receptor=NULL) {
  
  # Return dataframe with entries in the input ligand cell types AND receptor cell types.
  
  col_ctype_ligand = COL_LIGAND_CELLTYPE
  col_ctype_receptor = COL_RECEPTOR_CELLTYPE
  
  if ( (is.null(ctypes_ligand)) & (is.null(ctypes_receptor)) ) {
    
    stop("Invalid input: Ligand and receptors cell types are None.")
    
  } else if ( (!is.null(ctypes_ligand)) & (is.null(ctypes_receptor)) ) {
    
    ctypes_ligand_ = ctypes_ligand
    ctypes_receptor_ = unique(df[, col_ctype_receptor])
    
  } else if ( (is.null(ctypes_ligand)) & (!is.null(ctypes_receptor)) ) {
    
    ctypes_ligand_ = unique(df[, col_ctype_ligand])
    ctypes_receptor_ = ctypes_receptor
    
  } else {
    
    ctypes_ligand_ = ctypes_ligand
    ctypes_receptor_ = ctypes_receptor
    
  }
  
  if ( !all(ctypes_ligand_ %in% df[, col_ctype_ligand])
       | !all(ctypes_receptor_ %in% df[, col_ctype_receptor]) ) {
    stop("Invalid input: at least one celltype in ligand or receptor not found in dataframe.")
  }
  
  mask_l = df[, col_ctype_ligand] %in% ctypes_ligand_
  mask_r = df[, col_ctype_receptor] %in% ctypes_receptor_
  mask = (mask_l & mask_r)
  
  res = df[mask, ]
  return(res)
}


#filter_tissues_celltypes <- function(df, tissues, ctypes_ligand=NULL, ctypes_receptor=NULL) {
  
#  # Filter by tissue and ligand/receptor cell types
  
#  df_t = filter_tissues(df, tissues)
#  res = filter_celltypes(df_t, ctypes_ligand, ctypes_receptor)
#  return(res)
#}


filter_lr_genes <- function(df, lr_genes) {
  
  # Return dataframe 
  stop("NotImplemented")
}


filter_interactions <- function(df, tissues=NULL, ctypes_ligand=NULL, ctypes_receptor=NULL,
                                lr_genes=NULL) {
  
  # General-purpose filter on interactions dataframe
  # Returns dataframe with records satysfing 
  #  input interactions (tissue, l cell, r cell and lr genes)
  stop("NotImplemented")
  return(-1)
}


count_interactions <- function(df) {
  
  # Returns a one-column dataframe with interactions' counts.

  # Counts:
  #   - total number of interactions
  #   - number of up-regulated interactions
  #   - number of down-regulated interactions
  
  col_score_old = COL_LR_SCORE_OLD
  col_score_young = COL_LR_SCORE_YOUNG
  
  score_diff = df[, col_score_old] - df[, col_score_young]
  
  count_total = length(score_diff)
  count_up = sum(score_diff > 0)
  count_down = sum(score_diff < 0)

  counts <- data.frame(num_interactions=c(count_up, count_down, count_total),
                       row.names=c("up", "down", "total"))
  return(counts)
}


summarize_dataframe <- function(df) {
  
  # Summarizes the analysis dataframe with useful statistics.
  
  check_representation_inv(df)
  
  col_all = "Interactions_all"
  col_up = "Interactions_up"
  col_down = "Interactions_down"
  
  dt = data.table(df)
  
  # Get counts overall
  dt_all = dt[, .N, by=""]
  dt_all_up = dt[log(LR_SCORE_old/LR_SCORE_young) > 0, .N, by=""]
  dt_all_down = dt[log(LR_SCORE_old/LR_SCORE_young) < 0, .N, by=""]
  dt_all_res = cbind(dt_all, dt_all_up, dt_all_down)
  colnames(dt_all_res) = c(col_all, col_up, col_down)
  
  # Get counts at tissue level
  dt_tissues = dt[, .N, by=tissue]
  setkey(dt_tissues, tissue)
  dt_tissues_up = dt[log(LR_SCORE_old/LR_SCORE_young) > 0, .N, by=tissue]
  dt_tissues_down = dt[log(LR_SCORE_old/LR_SCORE_young) < 0, .N, by=tissue]
  dt_tissues_res = merge(dt_tissues, dt_tissues_up, all=TRUE)
  dt_tissues_res = merge(dt_tissues_res, dt_tissues_down, all=TRUE)
  colnames(dt_tissues_res) = c("Tissues", col_all, col_up, col_down)
  
  # Get counts at cell types level
  dt_ctypes = dt[, .N, by=.(L_CELLTYPE, R_CELLTYPE)]
  setkey(dt_ctypes, L_CELLTYPE, R_CELLTYPE)
  dt_ctypes_up = dt[log(LR_SCORE_old/LR_SCORE_young) > 0, .N, by=.(L_CELLTYPE, R_CELLTYPE)]
  dt_ctypes_down = dt[log(LR_SCORE_old/LR_SCORE_young) < 0, .N, by=.(L_CELLTYPE, R_CELLTYPE)]
  dt_ctypes_res = merge(dt_ctypes, dt_ctypes_up, all=TRUE)
  dt_ctypes_res = merge(dt_ctypes_res, dt_ctypes_down, all=TRUE)
  colnames(dt_ctypes_res) = c(COL_LIGAND_CELLTYPE, COL_RECEPTOR_CELLTYPE, col_all, col_up, col_down)
  
  # Aggregate stats in list
  setDF(dt_all_res)
  setDF(dt_tissues_res)
  setDF(dt_ctypes_res)
  res = list(totals=dt_all_res, tissues=dt_tissues_res, celltypes=dt_ctypes_res)
  
  return(res)
}


summarize_filtering <- function(df_in, df_out) {
  
  # Generates a summary for the filtering operation(s) outcome.
  # Operates on the input and output of the entire filtering pipeline.
  # Args:
  #   df_in: input dataframe to post-processing
  #   df_out: output dataframe of post-processing
  
  check_representation_inv(df_in)
  check_representation_inv(df_out)
  
  mask = !(rownames(df_in) %in% rownames(df_out))
  
  df_removed = df_in[mask, ]
  
  summarize_removed <- summarize_dataframe(df_removed)
  return(summarize_removed)
}


save_summary <- function(summ, name, base_dir) {
  
  dir = paste(base_dir, name, "/", sep="")
  create_dir(dir)  # introduces indirect dependencies, to mv to upstream code
  
  fn_totals = paste(dir, "totals", ".csv", sep="")
  fn_tissues = paste(dir, "tissues", ".csv", sep="")
  fn_celltypes = paste(dir, "celltypes", ".csv", sep="")
  
  fwrite(summ$totals, fn_totals, row.names=FALSE, col.names=TRUE, compress="none")
  fwrite(summ$tissues, fn_tissues, row.names=FALSE, col.names=TRUE, compress="none")
  fwrite(summ$celltypes, fn_celltypes, row.names=FALSE, col.names=TRUE, compress="none")
}


