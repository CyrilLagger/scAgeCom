#!/usr/bin/env Rscript


###   Downstream processing script.
###
###   USAGE:
###
###  1. For exploring cutoffs:
###     Rscript downstream_processing.R --explorecutoffs <path-to-dir>/  
###     # for now don't forget the forwardslash, will correct in future.
###
###  2. For filtering and overrepresentation analysis (ORA):
###     Rscript downstream_processing.R
###
###   Make sure the config file is in the same folder and appropriately setup.
###   Also make sure that no directory exists with same path as results path for either
###    use case 1) or 2). The script will raise an exception in such case. The file to be
###    configured is "config.yml", the "default" section. Other sections can be deleted in 
###    that file, I use them for testing, experimenting.
###   To shut off a filter, just insert NULL at the coresponding CUTOFF in config.yml. Most
###    filters can be disabled this way, except the old/young detected filter, which stays
###    always active.



# TODOS:
# Aggregated adjpval per category.


library(optparse)
library(config)
library(conflicted)
library(data.table)  # loading it changes behavior of existing code, reqs refactoring
source("src/downstream.R")
source("src/overrepresentation.R")
source("src/utils.R")
source("src/pipeline_steps.R")


# process cli args, will embbed in a function
explore_cutoffs_dir = read_explorecutoffs_arg()
if (is.null(explore_cutoffs_dir)) {
  TO_EXPLORE_CUTOFFS = FALSE
} else {
  TO_EXPLORE_CUTOFFS = TRUE
}


# `merge` from config is in conflict with base merge used in downstream.R
conflicted::conflict_prefer(name = "merge", winner = "base")


# load constants
load_globals_from_config(R_CONFIG_ACTIVE = "default")


main <- function() {

  df = readRDS(DATA_PATH)
  setDF(df)

  df = preprocess_input_df(df)
  check_representation_inv(df)

  
  # FILTER CUTOFFS EXPLORATION
  if (TO_EXPLORE_CUTOFFS) {
    create_dir(explore_cutoffs_dir)
    df_detec_in_young_and_old = remove_undetected_interactions(df)
    exlpore_filter_cutoffs(df_detec_in_young_and_old, explore_cutoffs_dir)
  } 
  
  else {
    
    create_dir(RESULTS_DIR)
    copy_config(RESULTS_DIR)
    
    # logging by redirecting output to file
    file_log = paste(RESULTS_DIR, "log.txt", sep='')
    sink(file=file_log)
    
    
    # FILTERING
    res_filtering_dir = paste0(RESULTS_DIR, "/", "FILTERING")
    create_dir(res_filtering_dir)
    df_detected = run_filtering(df, filter_detected, 
                                'raw_to_detected', res_filtering_dir)
    df_significant = run_filtering(df_detected, filter_significant, 
                                   'detected_to_significant', res_filtering_dir)
    
    # OVERREPRESENTATION
    print("Performing ORA")
    overrepr_res_dir = paste(RESULTS_DIR, "ORA", "/", sep='')
    create_dir(overrepr_res_dir)
    
    overrepr = analyze_overrepresentation(df_detected, df_significant, adjust_pvals=TRUE)
    # overrepr = adjust_pvals_of_list_dfs(overrepr, col_pval="overrepr_pval")

    dir_res_over_all_interacts = paste(overrepr_res_dir,
                                      'ALL/', sep='')
    save_list_dfs(overrepr, dir_res_over_all_interacts)
    
    
    analyze_overrepr_per_tissue(df_detected, df_significant, overrepr_res_dir, adjust_pvals=TRUE)

    sink()  # back to stdout
  }
}


main()
