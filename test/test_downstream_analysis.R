####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - June 2020
##
## Test the downstream pipeline. 
##
####################################################
##

library(Seurat)
library(scDiffCom)
library(data.table)

#Function to add useful column and rbind all tissues in a single data.table
bind_tissues <- function(
  path,
  list_of_tissues
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
       dt[, LR_diff := LR_SCORE_old - LR_SCORE_young]
      }),
    use.names = TRUE
  )
}

#datasets of interest
datasets <- c("calico", "tms_facs", "tms_droplet")
#path of datasets
data_path <- list(
  calico = "scDiffCom_results/diffcom_calico_size_factor_log_10000iter_mixed",
  tms_facs = "scDiffCom_results/diffcom_tms_facs_size_factor_log_10000iter_mixed",
  tms_droplet = "scDiffCom_results/diffcom_tms_droplet_size_factor_log_10000iter_mixed"
)
#tissues of interest
tissue_list <- list(
  calico = c("kidney", "lung", "spleen"),
  tms_facs = c("Aorta", "BAT", "Bladder", "Brain_Myeloid",
               "Brain_Non-Myeloid", "Diaphragm", "GAT",
               "Heart", "Kidney", "Large_Intestine",
               "Limb_Muscle", "Liver", "Lung", "Mammary_Gland",
               "Marrow", "MAT", "Pancreas", "SCAT", "Skin",
               "Spleen", "Thymus", "Tongue", "Trachea"),
  tms_droplet = c("Bladder", "Heart_and_Aorta", "Kidney",
                  "Limb_Muscle", "Liver", "Lung",
                  "Mammary_Gland", "Marrow", "Spleen",
                  "Thymus", "Tongue")
)


#list of all results
diffcom_results <- mapply(bind_tissues, data_path, tissue_list, SIMPLIFY = FALSE)

#a test file as well
diffcom_test <- readRDS("test_and_comparison/data_results_diffcom.rds")
diffcom_test[, TISSUE := "Liver"]
diffcom_test[, LR_diff := LR_SCORE_old - LR_SCORE_young]

### Eugen code adapated from command line to in-script



#Explore cutoffs
TO_EXPLORE_CUTOFFS = TRUE
explore_cutoffs_dir = "test_and_comparison/explore_cutoffs/"





######################

# Add relevant columns
# correct downstream logfc

diffcom_test$LR_CELLTYPES <- paste(diffcom_test$L_CELLTYPE, diffcom_test$R_CELLTYPE, sep = "_")





#
library(optparse)
library(config)
library(conflicted)
library(data.table)  # loading it changes behavior of existing code, reqs refactoring
#
source("downstream_script/src/downstream.R")
source("downstream_script/src/overrepresentation.R")
source("downstream_script/src/utils.R")
source("downstream_script/src/pipeline_steps.R")

#
# `merge` from config is in conflict with base merge used in downstream.R
conflicted::conflict_prefer(name = "merge", winner = "base")
#
# load constants
load_globals_from_config(R_CONFIG_ACTIVE = "default")
RESULTS_DIR <- paste0("test_and_comparison/", RESULTS_DIR)
#
main <- function() {
  df = diffcom_test
  #df = readRDS(DATA_PATH)
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

calico_ora_lung <- read.csv("test_and_comparison/downstream_test_calico/ORA/lung/LR_GENES.tsv", header = TRUE, sep =  "\t")

