####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - July 2020
##
## Combine scDiffCom results and filtering.
##
####################################################
##

## Libraries ####
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)

## Sources ####
source("src/src_1_filtering.R")

## Data information ####

# general data path
dir_data <- "../data_scAgeCom/"

# path of scDiffCom results
diffcom_res_path <- list(
  calico = "../data_scAgeCom/scDiffCom_results/diffcom_calico_size_factor_log_10000iter_mixed",
  calico_sub = "../data_scAgeCom/scDiffCom_results/diffcom_calico_subtype_size_factor_log_10000iter_mixed",
  tms_facs = "../data_scAgeCom/scDiffCom_results/diffcom_tms_facs_size_factor_log_10000iter_mixed",
  tms_droplet = "../data_scAgeCom/scDiffCom_results/diffcom_tms_droplet_size_factor_log_10000iter_mixed"
)

# tissues of interest
tissue_list <- list(
  calico = c("Kidney", "Lung", "Spleen"),
  calico_sub = c("Kidney", "Lung", "Spleen"),
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

## Load scDiffCom results and bind them in single data.tables ####

diffcom_results <- mapply(
  bind_tissues,
  diffcom_res_path,
  tissue_list,
  MoreArgs = list(is_log = TRUE),
  SIMPLIFY = FALSE
)

## Explore the cutoffs ####

#Note: this part can be modified as proposed by Eugen with the functions "explore_cutoff"

# cutoff based on top 20% CCIs
cutoff_score_young_25pct <- lapply(
  diffcom_results,
  function(x) {
    quantile(x[LR_DETECTED_young == TRUE]$LR_SCORE_young, 0.75)
  }
)
cutoff_score_old_25pct <- lapply(
  diffcom_results,
  function(x) {
    quantile(x[LR_DETECTED_old == TRUE]$LR_SCORE_old, 0.75)
  }
)

cutoff_score_25_pct <- mapply(
  FUN = function(x,y) {min(c(x,y))},
  cutoff_score_young_25pct,
  cutoff_score_young_25pct,
  SIMPLIFY = FALSE
)


# some histograms
cowplot::plot_grid(
  plotlist = lapply(
    diffcom_results,
    function(x) {
      ggplot(x[LR_DETECTED_young == TRUE], aes(x = LR_SCORE_young)) + 
        geom_histogram(bins = 100) #+
      #scale_y_log10()
    }
  ),
  nrow = 2,
  align = "v"
)

## Filtering process and save the results ####

diffcom_filter <- mapply(
  FUN = function(x,y) {
    filter_and_reassign_CCI(
      dt = x,
      CUTOFF_SCORE_YOUNG = y,
      CUTOFF_SCORE_OLD = y,
      CUTOFF_LOGFC = log(1.1)
    )
  },
  diffcom_results,
  cutoff_score_25_pct,
  SIMPLIFY = FALSE
)

saveRDS(diffcom_filter, file = paste0(dir_data, "analysis/analysis_4_data_diffcom_filter.rds"))

## Overall detected CCIs ####

# data.frame summarizing the 6 possible scenarios
six_scenarios <- data.frame(
  Young = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  Old = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
  Diff = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE),
  Direction = c("Up", "Down", "", "Down", "UP", "")
)
g_six <- tableGrob(six_scenarios, rows = NULL)
grid.newpage()
grid.draw(g_six)


# detection per dataset and per scenario
distr_scenar <- merge.data.table(
  merge.data.table(
    diffcom_filter$tms_facs[,.N, by = SIG_TYPE],
    diffcom_filter$tms_droplet[,.N, by = SIG_TYPE],
    by = "SIG_TYPE"
  ),
  diffcom_filter$calico[,.N, by = SIG_TYPE],
  by = "SIG_TYPE"
)
colnames(distr_scenar) <- c("SIG_TYPE", "TMS FACS", "TMS Droplet", "Calico")
setDT(six_scenarios)
six_scenarios[, SIG_TYPE := c("TTTU", "TTTD", "TTF", "TFT", "FTT", "FFF")]
distr_scenar <- merge.data.table(
  six_scenarios,
  distr_scenar,
  by = "SIG_TYPE",
  sort = FALSE
)
distr_scenar[, SIG_TYPE := NULL]
setcolorder(distr_scenar, neworder = c("Young", "Old", "Diff", "Direction", "TMS FACS", "TMS Droplet", "Calico"))
g_distr_all <- tableGrob(distr_scenar, rows = NULL)
grid.newpage()
grid.draw(g_distr_all)

ggsave(filename = paste0(dir_data, "analysis/analysis_4_plot_detection_distr.png") , plot = g_distr_all, scale = 1.5)



