# Utils

create_dir <- function(dir) {
  if (!(dir.exists(dir))) {
    dir.create(dir)
  } else {
    stop("Directory already exists.")
  }
}


copy_config <- function(dir) {
  
  if ( !file.exists("config.yml") ) {
    stop("config.yml not found.")
  }
  
  if ( !(dir.exists(dir)) ) {
    stop("Directory already exists.")
  } else {
    file.copy("config.yml", dir)
  }
}


read_explorecutoffs_arg <- function() {
  
  option_list = list(
    make_option(c("--explorecutoffs"), type="character", default=NULL, 
                help="directory where to explore-cutoff", metavar="character")
  ); 
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  
  if (is.null(opt$explorecutoffs)){
    return(NULL)
  } else {
    return(opt$explorecutoffs)
  }
}


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


save_list_dfs <- function(l, dir) {
  
  is_not_list = !(class(l) == "list")
  is_first_not_df = (class(l[[1]]) != "data.frame")
  
  if ( is_not_list | is_first_not_df) {
    stop("Input must be list of dfs.")
  }
  
  create_dir(dir)
  
  for (name in names(l)) {
    
    df = l[[name]]
    
    # The issue with below commented code is that some cell types have ',' in
    #  their names. One solution is to save as .tsv
    
    #    fname = paste(dir, name, ".csv", sep="")
    #    write.csv(df, fname, quote=FALSE, row.names=TRUE, col.names=TRUE)
    
    fname = paste(dir, name, ".tsv", sep="")
    write.table(df, fname, sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
  }
  
}
