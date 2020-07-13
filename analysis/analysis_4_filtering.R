####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - June 2020
##
## Combine the results from scDiffCom for the
## various datasets considered here. Then load the
## downstream parameters and proceed to filtering
## and overrepresenation analysis.
##
####################################################
##

library(Seurat)
library(scDiffCom)
library(data.table)
library(config)
library(conflicted)
library(circlize)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
library(clusterProfiler)
library(org.Mm.eg.db)

LRone2one <- LRall$LRall_one2one

#Function to add useful column (and remove others) and rbind all tissues in a single data.table
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
        dt[, LR_LOGFC := LR_SCORE_old - LR_SCORE_young]
        dt[, LR_LOGFC_ABS := abs(LR_LOGFC)]
        dt[, LR_CELLTYPES := paste(L_CELLTYPE, R_CELLTYPE, sep = "_")]
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
    use.names = TRUE
  )
}

#####
#datasets of interest
datasets <- c("calico", "calico_sub", "tms_facs", "tms_droplet")
#path of datasets
data_path <- list(
  calico = "../data_scAgeCom/scDiffCom_results/diffcom_calico_size_factor_log_10000iter_mixed",
  calico_sub = "../data_scAgeCom/scDiffCom_results/diffcom_calico_subtype_size_factor_log_10000iter_mixed",
  tms_facs = "../data_scAgeCom/scDiffCom_results/diffcom_tms_facs_size_factor_log_10000iter_mixed",
  tms_droplet = "../data_scAgeCom/scDiffCom_results/diffcom_tms_droplet_size_factor_log_10000iter_mixed"
)
#tissues of interest
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


#list of all results
diffcom_results <- mapply(bind_tissues, data_path, tissue_list, SIMPLIFY = FALSE)
#diffcom_results <- lapply(diffcom_results, setDF)

#check that the format and colnames are OK
#check_representation_inv(diffcom_results$calico)

#Cutoff exploration
#df_detec_in_young_and_old <- lapply(
#  diffcom_results,
#  remove_undetected_interactions
#)
#calico
#explore_cutoffs_dir_calico <- "test_and_comparison/explore_cutoffs/calico/"
#create_dir(explore_cutoffs_dir_calico)
#explore_calico <- explore_filter_cutoffs(df_detec_in_young_and_old$calico, explore_cutoffs_dir_calico) 
#tms FACS
#explore_cutoffs_dir_tms_facs <- "test_and_comparison/explore_cutoffs/tms_facs/"
#create_dir(explore_cutoffs_dir_tms_facs)
#explore_tms_facs <- explore_filter_cutoffs(df_detec_in_young_and_old$tms_facs, explore_cutoffs_dir_tms_facs) 
#tms Droplet
#explore_cutoffs_dir_tms_droplet <- "test_and_comparison/explore_cutoffs/tms_droplet/"
#create_dir(explore_cutoffs_dir_tms_droplet)
#explore_tms_droplet <- explore_filter_cutoffs(df_detec_in_young_and_old$tms_droplet, explore_cutoffs_dir_tms_droplet)

#explore_calico
#explore_tms_facs
#explore_tms_droplet

# Filtering process
filter_CCI <- function(
  dt,
  CUTOFF_SCORE_YOUNG,
  CUTOFF_SCORE_OLD,
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ= 0.05,
  CUTOFF_LOGFC
  ) {
  dt[, LR_KEEP_young := ifelse(
    (LR_DETECTED_young == TRUE & LR_SCORE_young >= CUTOFF_SCORE_YOUNG & BH_PVAL_young <= CUTOFF_CPDB_PVAL_YOUNG),
    TRUE,
    FALSE
  )]
  dt[, LR_KEEP_old := ifelse(
    (LR_DETECTED_old == TRUE & LR_SCORE_old >= CUTOFF_SCORE_OLD & BH_PVAL_old <= CUTOFF_CPDB_PVAL_OLD),
    TRUE,
    FALSE
  )]
  dt[, LR_KEEP_DIFF := ifelse(
    (LR_KEEP_young | LR_KEEP_old) & BH_PVAL_DIFF <= CUTOFF_PVAL_ADJ & LR_LOGFC_ABS >= CUTOFF_LOGFC,
    TRUE,
    FALSE
  )]
  dt[, DIRECTION := ifelse(LR_LOGFC > 0, "UP", "DOWN")]
  # dt[, LR_UP := ifelse(
  #   LR_KEEP_DIFF == TRUE & LR_KEEP_old == TRUE & LR_LOGFC > 0, 
  #   TRUE, 
  #   FALSE
  # )]
  # dt[, LR_DOWN := ifelse(
  #   LR_KEEP_DIFF == TRUE & LR_KEEP_young == TRUE & LR_LOGFC < 0, 
  #   TRUE, 
  #   FALSE
  # )]
  return(dt)
}

reassign_CCI <- function(
  dt,
  CUTOFF_SCORE_YOUNG,
  CUTOFF_SCORE_OLD,
  CUTOFF_CPDB_PVAL_YOUNG= 0.05,
  CUTOFF_CPDB_PVAL_OLD= 0.05,
  CUTOFF_PVAL_ADJ= 0.05,
  CUTOFF_LOGFC
) {
  dt[, SIG_TYPE := ifelse(
    LR_KEEP_young == TRUE & LR_KEEP_old == TRUE &  LR_KEEP_DIFF == TRUE,
    ifelse(
      DIRECTION == "UP",
      "TTTU",
      "TTTD"
    ),
    ifelse(
      LR_KEEP_young == TRUE & LR_KEEP_old == TRUE &  LR_KEEP_DIFF == FALSE,
      "TTF",
      ifelse(
        LR_KEEP_young == TRUE & LR_KEEP_old == FALSE &  LR_KEEP_DIFF == TRUE,
        ifelse(
          DIRECTION == "DOWN",
          "TFT",
          ifelse(
            sum(c(BH_PVAL_old > CUTOFF_CPDB_PVAL_OLD, LR_DETECTED_old == FALSE, LR_SCORE_old < CUTOFF_SCORE_OLD)) == 1,
            "TTTU",
            "FFF"
          )
        ),
        ifelse(
          LR_KEEP_young == TRUE & LR_KEEP_old == FALSE &  LR_KEEP_DIFF == FALSE,
          ifelse(
            sum(c(BH_PVAL_old > CUTOFF_CPDB_PVAL_OLD, LR_DETECTED_old == FALSE, LR_SCORE_old < CUTOFF_SCORE_OLD)) == 1,
            "TTF",
            "FFF"
          ),
          ifelse(
            LR_KEEP_young == FALSE & LR_KEEP_old == TRUE &  LR_KEEP_DIFF == TRUE,
            ifelse(
              DIRECTION == "UP",
              "FTT",
              ifelse(
                sum(c(BH_PVAL_young > CUTOFF_CPDB_PVAL_YOUNG, LR_DETECTED_young == FALSE, LR_SCORE_young < CUTOFF_SCORE_YOUNG)) == 1,
                "TTTD",
                "FFF"
              )
            ),
            ifelse(
              LR_KEEP_young == FALSE & LR_KEEP_old == TRUE &  LR_KEEP_DIFF == FALSE,
              ifelse(
                sum(c(BH_PVAL_young > CUTOFF_CPDB_PVAL_YOUNG, LR_DETECTED_young == FALSE, LR_SCORE_young < CUTOFF_SCORE_YOUNG)) == 1,
                "TTF",
                "FFF"
              ),
              "FFF"
            )
          )
        )
      )
    )
  )]
}

hist(diffcom_results$tms_facs[LR_DETECTED_young == TRUE,]$LR_SCORE_young, breaks = 50)
quantile(diffcom_results$tms_facs[LR_DETECTED_young == TRUE,]$LR_SCORE_young, probs = seq(0, 1, 0.1))
quantile(diffcom_results$tms_facs[LR_DETECTED_old == TRUE,]$LR_SCORE_old, probs = seq(0, 1, 0.1))
filter_CCI(
  dt = diffcom_results$tms_facs,
  CUTOFF_SCORE_YOUNG = 0.9, #case dependent
  CUTOFF_SCORE_OLD = 0.9, # case dependent
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05,
  CUTOFF_LOGFC = log(1.1)
)

hist(diffcom_results$tms_droplet[LR_DETECTED_young == TRUE,]$LR_SCORE_young, breaks = 50)
quantile(diffcom_results$tms_droplet[LR_DETECTED_young == TRUE,]$LR_SCORE_young, probs = seq(0, 1, 0.1))
quantile(diffcom_results$tms_droplet[LR_DETECTED_old == TRUE,]$LR_SCORE_old, probs = seq(0, 1, 0.1))
filter_CCI(
  dt = diffcom_results$tms_droplet,
  CUTOFF_SCORE_YOUNG = 1.5, #case dependent
  CUTOFF_SCORE_OLD = 1.5, # case dependent
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05,
  CUTOFF_LOGFC = log(1.1)
)

hist(diffcom_results$calico[LR_DETECTED_young == TRUE,]$LR_SCORE_young, breaks = 50)
quantile(diffcom_results$calico[LR_DETECTED_young == TRUE,]$LR_SCORE_young, probs = seq(0, 1, 0.1))
quantile(diffcom_results$calico[LR_DETECTED_old == TRUE ,]$LR_SCORE_old, probs = seq(0, 1, 0.1))
filter_CCI(
  dt = diffcom_results$calico,
  CUTOFF_SCORE_YOUNG = 1.25, #case dependent
  CUTOFF_SCORE_OLD = 1.25, # case dependent
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05,
  CUTOFF_LOGFC = log(1.1)
)


hist(diffcom_results$calico_sub[LR_DETECTED_young == TRUE,]$LR_SCORE_young, breaks = 50)
quantile(diffcom_results$calico_sub[LR_DETECTED_young == TRUE,]$LR_SCORE_young, probs = seq(0, 1, 0.1))
quantile(diffcom_results$calico_sub[LR_DETECTED_old == TRUE ,]$LR_SCORE_old, probs = seq(0, 1, 0.1))
filter_CCI(
  dt = diffcom_results$calico_sub,
  CUTOFF_SCORE_YOUNG = 1.25, #case dependent
  CUTOFF_SCORE_OLD = 1.25, # case dependent
  CUTOFF_CPDB_PVAL_YOUNG = 0.05,
  CUTOFF_CPDB_PVAL_OLD = 0.05,
  CUTOFF_PVAL_ADJ = 0.05,
  CUTOFF_LOGFC = log(1.1)
)

#reassign
diffcom_results_detected <- list(
  tms_facs = reassign_CCI(
    dt = diffcom_results$tms_facs[LR_KEEP_young == TRUE | LR_KEEP_old == TRUE,],
    CUTOFF_SCORE_YOUNG = 0.9, #case dependent
    CUTOFF_SCORE_OLD = 0.9, # case dependent
    CUTOFF_CPDB_PVAL_YOUNG = 0.05,
    CUTOFF_CPDB_PVAL_OLD = 0.05,
    CUTOFF_PVAL_ADJ = 0.05,
    CUTOFF_LOGFC = log(1.1)
  ),
  tms_droplet = reassign_CCI(
    dt = diffcom_results$tms_droplet[LR_KEEP_young == TRUE | LR_KEEP_old == TRUE,],
    CUTOFF_SCORE_YOUNG = 1.5, #case dependent
    CUTOFF_SCORE_OLD = 1.5, # case dependent
    CUTOFF_CPDB_PVAL_YOUNG = 0.05,
    CUTOFF_CPDB_PVAL_OLD = 0.05,
    CUTOFF_PVAL_ADJ = 0.05,
    CUTOFF_LOGFC = log(1.1)
  ),
  calico = reassign_CCI(
    dt = diffcom_results$calico[LR_KEEP_young == TRUE | LR_KEEP_old == TRUE,],
    CUTOFF_SCORE_YOUNG = 1.25, #case dependent
    CUTOFF_SCORE_OLD = 1.25, # case dependent
    CUTOFF_CPDB_PVAL_YOUNG = 0.05,
    CUTOFF_CPDB_PVAL_OLD = 0.05,
    CUTOFF_PVAL_ADJ = 0.05,
    CUTOFF_LOGFC = log(1.1)
  ),
  calico_sub = reassign_CCI(
    dt = diffcom_results$calico_sub[LR_KEEP_young == TRUE | LR_KEEP_old == TRUE,],
    CUTOFF_SCORE_YOUNG = 1.25, #case dependent
    CUTOFF_SCORE_OLD = 1.25, # case dependent
    CUTOFF_CPDB_PVAL_YOUNG = 0.05,
    CUTOFF_CPDB_PVAL_OLD = 0.05,
    CUTOFF_PVAL_ADJ = 0.05,
    CUTOFF_LOGFC = log(1.1)
  )
)



#####
#a data.frame summarizing the 6 possible scenarios
six_scenarios <- data.frame(
  Young = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  Old = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
  Diff = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE),
  Direction = c("Up", "Down", "", "Down", "UP", "")
)
g_six <- tableGrob(six_scenarios, rows = NULL)
grid.newpage()
grid.draw(g_six)
#ggsave(filename = "../data_scAgeCom/scenario_tables.png", plot = g_six, scale = 1.5)

#distribution of detection per dataset
distr_scenar <- merge.data.table(
  merge.data.table(
  diffcom_results_detected$tms_facs[,.N, by = SIG_TYPE],
  diffcom_results_detected$tms_droplet[,.N, by = SIG_TYPE],
  by = "SIG_TYPE"
  ),
  diffcom_results_detected$calico[,.N, by = SIG_TYPE],
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
#ggsave(filename = "../data_scAgeCom/distr_scenarios_all_data.png", plot = g_distr_all, scale = 1.5)

######
#Some sub-datatables
diffcom_results_detected_strong <- lapply(
  diffcom_results_detected,
  function(x) {
    dt <- x[!(SIG_TYPE == "FFF"),]
  }
)

diffcom_results_significant <- lapply(
  diffcom_results_detected_strong,
  function(x){
    dt <- x[SIG_TYPE %in% c("FTT", "TFT", "TTTD", "TTTU"),]
  }
)

##enrichment on all results
L_Up <- list(
  tms_facs = unique(diffcom_results_detected$tms_facs[SIG_TYPE %in% c("TTTU", "FTT"),]$L_GENE),
  tms_droplet = unique(diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("TTTU", "FTT"),]$L_GENE),
  calico = unique(diffcom_results_detected$calico[SIG_TYPE %in% c("TTTU", "FTT"),]$L_GENE)
)
R_Up <- list(
  tms_facs = unique(diffcom_results_detected$tms_facs[SIG_TYPE %in% c("TTTU", "FTT"),]$R_GENE),
  tms_droplet = unique(diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("TTTU", "FTT"),]$R_GENE),
  calico = unique(diffcom_results_detected$calico[SIG_TYPE %in% c("TTTU", "FTT"),]$R_GENE)
)
LR_Up <- list(
  tms_facs = unique(c(L_Up$tms_facs, R_Up$tms_facs) ),
  tms_droplet = unique(c(L_Up$tms_droplet, R_Up$tms_droplet) ),
  calico = unique(c(L_Up$calico, R_Up$calico) )
)

L_Down <- list(
  tms_facs = unique(diffcom_results_detected$tms_facs[SIG_TYPE %in% c("TTTD", "TFT"),]$L_GENE),
  tms_droplet = unique(diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("TTTD", "TFT"),]$L_GENE),
  calico = unique(diffcom_results_detected$calico[SIG_TYPE %in% c("TTTD", "TFT"),]$L_GENE)
)
R_Down <- list(
  tms_facs = unique(diffcom_results_detected$tms_facs[SIG_TYPE %in% c("TTTD", "TFT"),]$R_GENE),
  tms_droplet = unique(diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("TTTD", "TFT"),]$R_GENE),
  calico = unique(diffcom_results_detected$calico[SIG_TYPE %in% c("TTTD", "TFT"),]$R_GENE)
)
LR_Down <- list(
  tms_facs = unique(c(L_Down$tms_facs, R_Down$tms_facs) ),
  tms_droplet = unique(c(L_Down$tms_droplet, R_Down$tms_droplet) ),
  calico = unique(c(L_Down$calico, R_Down$calico) )
)

universe_enrich_L <- unique(LRone2one[scsr == TRUE | cpdb == TRUE, ]$GENESYMB_L)
universe_enrich_R <- unique(LRone2one[scsr == TRUE | cpdb == TRUE, ]$GENESYMB_R)
universe_enrich_LR <- unique(c(universe_enrich_L, universe_enrich_R))

ego_analysis_list <- list(
  LR_Up$tms_facs, LR_Down$tms_facs,
  LR_Up$tms_droplet, LR_Down$tms_droplet,
  LR_Up$calico, LR_Down$calico
)

ego_BP_overall <- lapply(
  ego_analysis_list,
  function(x) {
    enrichGO(
      gene = x,
      OrgDb = org.Mm.eg.db,
      keyType = 'SYMBOL',
      ont = "BP",
      universe = universe_enrich_LR,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff  = 0.05
    )
  }
)

ego_BP_overall_plots <- lapply(
  ego_BP_overall,
  function(x) {
    dotplot(
      x,
      showCategory = 10,
      font.size = 10
    )
  }
)

g_ego_BP_overall <- cowplot::plot_grid(
  plotlist = ego_BP_overall_plots,
  ncol = 2,
  align = "v",
  labels =  c("TMS FACS Up", "TMS FACS Down",
              "TMS Droplet Up", "TMS Droplet Down",
              "Calico Up", "Calico Down"
              )
)

ggsave(filename = "../data_scAgeCom/ego_BP_overall.png", plot = g_ego_BP_overall, scale = 1.8)



#######
# OVERREPRESENTATION

#load parameters
load_globals_from_config(R_CONFIG_ACTIVE = "default")
conflicted::conflict_prefer(name = "merge", winner = "base")

lapply(diffcom_results_detected_strong, function(x) {
  setDF(x)
})

lapply(diffcom_results_significant, function(x) {
  setDF(x)
})
#
overrep_calico <- analyze_overrepresentation(
  diffcom_results_detected_strong$calico,
  diffcom_results_significant$calico,
  adjust_pvals = TRUE
)
overrep_calico_tiss <- sapply(unique(diffcom_results_detected_strong$calico$TISSUE), function(tiss) {
  analyze_overrepresentation(
    diffcom_results_detected_strong$calico[diffcom_results_detected_strong$calico$TISSUE == tiss,],
    diffcom_results_significant$calico[diffcom_results_significant$calico$TISSUE == tiss,],
    exclude_tissue_overrepresentation=TRUE,
    adjust_pvals=TRUE)
}, USE.NAMES = TRUE, simplify = FALSE)
##
overrep_tms_facs <- analyze_overrepresentation(
  diffcom_results_detected_strong$tms_facs,
  diffcom_results_significant$tms_facs,
  adjust_pvals = TRUE
)
overrep_tms_facs_tiss <- sapply(unique(diffcom_results_detected_strong$tms_facs$TISSUE), function(tiss) {
  analyze_overrepresentation(
    diffcom_results_detected_strong$tms_facs[diffcom_results_detected_strong$tms_facs$TISSUE == tiss,],
    diffcom_results_significant$tms_facs[diffcom_results_significant$tms_facs$TISSUE == tiss,],
    exclude_tissue_overrepresentation=TRUE,
    adjust_pvals=TRUE)
}, USE.NAMES = TRUE, simplify = FALSE)
##
overrep_tms_droplet <- analyze_overrepresentation(
  diffcom_results_detected_strong$tms_droplet,
  diffcom_results_significant$tms_droplet,
  adjust_pvals = TRUE
)
overrep_tms_droplet_tiss <- sapply(unique(diffcom_results_detected_strong$tms_droplet$TISSUE), function(tiss) {
  analyze_overrepresentation(
    diffcom_results_detected_strong$tms_droplet[diffcom_results_detected_strong$tms_droplet$TISSUE == tiss,],
    diffcom_results_significant$tms_droplet[diffcom_results_significant$tms_droplet$TISSUE == tiss,],
    exclude_tissue_overrepresentation=TRUE,
    adjust_pvals=TRUE)
}, USE.NAMES = TRUE, simplify = FALSE)

#overall K overrep
overrep_L_overall <- list(
  tms_facs = overrep_tms_facs$L_GENE,
  tms_droplet = overrep_tms_droplet$L_GENE,
  calico = overrep_calico$L_GENE
)
lapply(overrep_L_overall, setDT)
overrep_L_sig <- list(
  overrep_L_overall$tms_facs[,c("L_GENE", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_L_overall$tms_droplet[,c("L_GENE", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_L_overall$calico[,c("L_GENE", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_L_overall$tms_facs[,c("L_GENE", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,],
  overrep_L_overall$tms_droplet[,c("L_GENE", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,],
  overrep_L_overall$calico[,c("L_GENE", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,]
)
overrep_L_plot <- lapply(overrep_L_sig, function(x)
  tableGrob(x, rows = NULL)
)
g_over_L <- cowplot::plot_grid(
  plotlist = overrep_L_plot,
  nrow = 2,
  align = "v",
  labels = c("FACS", "Droplet", "Calico", "FACS", "Droplet", "Calico")
)

#oerall R overrep
overrep_R_overall <- list(
  tms_facs = overrep_tms_facs$R_GENE,
  tms_droplet = overrep_tms_droplet$R_GENE,
  calico = overrep_calico$R_GENE
)
lapply(overrep_R_overall, setDT)
overrep_R_sig <- list(
  overrep_R_overall$tms_facs[,c("R_GENE", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_R_overall$tms_droplet[,c("R_GENE", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_R_overall$calico[,c("R_GENE", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_R_overall$tms_facs[,c("R_GENE", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,],
  overrep_R_overall$tms_droplet[,c("R_GENE", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,],
  overrep_R_overall$calico[,c("R_GENE", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,]
)
overrep_R_plot <- lapply(overrep_R_sig, function(x)
  tableGrob(x, rows = NULL)
)
g_over_R <- cowplot::plot_grid(
  plotlist = overrep_R_plot,
  nrow = 2,
  align = "v",
  labels = c("FACS", "Droplet", "Calico", "FACS", "Droplet", "Calico")
)


#overall LR overrep
overrep_LR_overall <- list(
  tms_facs = overrep_tms_facs$LR_GENES,
  tms_droplet = overrep_tms_droplet$LR_GENES,
  calico = overrep_calico$LR_GENES
)
lapply(overrep_LR_overall, setDT)
overrep_LR_sig <- list(
  overrep_LR_overall$tms_facs[,c("LR_GENES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_LR_overall$tms_droplet[,c("LR_GENES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_LR_overall$calico[,c("LR_GENES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_LR_overall$tms_facs[,c("LR_GENES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,],
  overrep_LR_overall$tms_droplet[,c("LR_GENES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,],
  overrep_LR_overall$calico[,c("LR_GENES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,]
)
overrep_LR_plot <- lapply(overrep_LR_sig, function(x)
  tableGrob(x, rows = NULL)
)
g_over_LR <- cowplot::plot_grid(
  plotlist = overrep_LR_plot,
  nrow = 2,
  align = "v",
  labels = c("FACS", "Droplet", "Calico", "FACS", "Droplet", "Calico")
)

#ggsave(filename = "../data_scAgeCom/overrep_LR_overall.png", plot = g_over_LR, scale = 1.7)


#tissue LR overrep
#facs
lapply(overrep_tms_facs_tiss, function(x) {
  setDT(x$LR_GENES)
})
overrep_tiss_facs_up <- lapply(overrep_tms_facs_tiss, function(x) {
  tableGrob(x$LR_GENES[,c("LR_GENES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:5,], rows = NULL)
})
overrep_tiss_facs_down <- lapply(overrep_tms_facs_tiss, function(x) {
  tableGrob(x$LR_GENES[,c("LR_GENES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:5,], rows = NULL)
})
g_over_LR_facs_tiss_up <- cowplot::plot_grid(
  plotlist = overrep_tiss_facs_up,
  ncol = 5,
  align = "v",
  labels = names(overrep_tms_facs_tiss)
)
g_over_LR_facs_tiss_down <- cowplot::plot_grid(
  plotlist = overrep_tiss_facs_down,
  ncol = 5,
  align = "v",
  labels = names(overrep_tms_facs_tiss)
)
#ggsave(filename = "../data_scAgeCom/overrep_LR_facs_tiss_up.png", plot = g_over_LR_facs_tiss_up, scale = 2)
#ggsave(filename = "../data_scAgeCom/overrep_LR_facs_tiss_down.png", plot = g_over_LR_facs_tiss_down, scale = 2.1)
#droplet
lapply(overrep_tms_droplet_tiss, function(x) {
  setDT(x$LR_GENES)
})
overrep_tiss_droplet_up <- lapply(overrep_tms_droplet_tiss, function(x) {
  tableGrob(x$LR_GENES[,c("LR_GENES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:5,], rows = NULL)
})
overrep_tiss_droplet_down <- lapply(overrep_tms_droplet_tiss, function(x) {
  tableGrob(x$LR_GENES[,c("LR_GENES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:5,], rows = NULL)
})
g_over_LR_droplet_tiss_up <- cowplot::plot_grid(
  plotlist = overrep_tiss_droplet_up,
  ncol = 4,
  align = "v",
  labels = names(overrep_tms_droplet_tiss)
)
g_over_LR_droplet_tiss_down <- cowplot::plot_grid(
  plotlist = overrep_tiss_droplet_down,
  ncol = 4,
  align = "v",
  labels = names(overrep_tms_droplet_tiss)
)
#ggsave(filename = "../data_scAgeCom/overrep_LR_droplet_tiss_up.png", plot = g_over_LR_droplet_tiss_up, scale = 1.5)
#ggsave(filename = "../data_scAgeCom/overrep_LR_droplet_tiss_down.png", plot = g_over_LR_droplet_tiss_down, scale = 1.7)
#calico
lapply(overrep_calico_tiss, function(x) {
  setDT(x$LR_GENES)
})
overrep_tiss_calico_up <- lapply(overrep_calico_tiss, function(x) {
  tableGrob(x$LR_GENES[,c("LR_GENES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:5,], rows = NULL)
})
overrep_tiss_calico_down <- lapply(overrep_calico_tiss, function(x) {
  tableGrob(x$LR_GENES[,c("LR_GENES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:5,], rows = NULL)
})
g_over_LR_calico_tiss_up <- cowplot::plot_grid(
  plotlist = overrep_tiss_calico_up,
  ncol = 2,
  align = "v",
  labels = names(overrep_calico_tiss)
)
g_over_LR_calico_tiss_down <- cowplot::plot_grid(
  plotlist = overrep_tiss_calico_down,
  ncol = 2,
  align = "v",
  labels = names(overrep_calico_tiss)
)
#ggsave(filename = "../data_scAgeCom/overrep_LR_calico_tiss_up.png", plot = g_over_LR_calico_tiss_up, scale = 1.1)
#ggsave(filename = "../data_scAgeCom/overrep_LR_calico_tiss_down.png", plot = g_over_LR_calico_tiss_down, scale = 1.1)

#overall CC overrep
overrep_CC_overall <- list(
  tms_facs = overrep_tms_facs$LR_CELLTYPES,
  tms_droplet = overrep_tms_droplet$LR_CELLTYPES,
  calico = overrep_calico$LR_CELLTYPES
)
lapply(overrep_CC_overall, setDT)
overrep_CC_sig <- list(
  overrep_CC_overall$tms_facs[,c("LR_CELLTYPES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_CC_overall$tms_droplet[,c("LR_CELLTYPES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_CC_overall$calico[,c("LR_CELLTYPES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:15,],
  overrep_CC_overall$tms_facs[,c("LR_CELLTYPES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,],
  overrep_CC_overall$tms_droplet[,c("LR_CELLTYPES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,],
  overrep_CC_overall$calico[,c("LR_CELLTYPES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:15,]
)
overrep_CC_plot <- lapply(overrep_CC_sig, function(x)
  tableGrob(x, rows = NULL)
)
g_over_CC <- cowplot::plot_grid(
  plotlist = overrep_CC_plot[c(1,2,4,5)],
  nrow = 2,
  align = "v",
  labels = c("FACS", "Droplet", "FACS", "Droplet")
)
#ggsave(filename = "../data_scAgeCom/overrep_CC_overall.png", plot = g_over_CC, scale = 2.2)

#tissue CC overrep
#facs
lapply(overrep_tms_facs_tiss, function(x) {
  setDT(x$LR_CELLTYPES)
})
overrep_CC_tiss_facs_up <- lapply(overrep_tms_facs_tiss, function(x) {
  tableGrob(x$LR_CELLTYPES[,c("LR_CELLTYPES", "overrepr_pval_UP")][order(overrepr_pval_UP)][1:5,], rows = NULL)
})
overrep_CC_tiss_facs_down <- lapply(overrep_tms_facs_tiss, function(x) {
  tableGrob(x$LR_CELLTYPES[,c("LR_CELLTYPES", "overrepr_pval_DOWN")][order(overrepr_pval_DOWN)][1:5,], rows = NULL)
})
g_over_CC_facs_tiss_up <- cowplot::plot_grid(
  plotlist = overrep_CC_tiss_facs_up,
  ncol = 5,
  align = "v",
  labels = names(overrep_tms_facs_tiss)
)
g_over_CC_facs_tiss_down <- cowplot::plot_grid(
  plotlist = overrep_CC_tiss_facs_down,
  ncol = 5,
  align = "v",
  labels = names(overrep_tms_facs_tiss)
)



##
SCAT_facs <- diffcom_results_detected$tms_facs[TISSUE == "SCAT",]

CCI_chord(SCAT_facs, "flat", TRUE, 10)
CCI_chord(SCAT_facs, "Up", TRUE, 5)
CCI_chord(SCAT_facs, "Down", TRUE, 5)

CCI_chord(SCAT_facs, "flat", FALSE, 10)
CCI_chord(SCAT_facs, "Up", FALSE, 5)
CCI_chord(SCAT_facs, "Down", FALSE, 5)

CCI_chord(SCAT_facs, "Down", TRUE, 10)
CCI_chord(SCAT_facs, "Down", FALSE, 10, "../data_scAgeCom/test.png")

gtest <- cowplot::as_grob(CCI_chord(SCAT_facs, "Down", TRUE, 10))

CCI_chord <- function(
  dt,
  direction,
  is_LR,
  filter,
  dir = NULL
) {
  dt_clean <- dt[SIG_TYPE != "FFF", ]
  if(direction == "Up") {
    dt_clean <- dt_clean[SIG_TYPE %in% c("FTT", "TTTU"),  ]
  } else if(direction == "Down") {
    dt_clean <- dt_clean[SIG_TYPE %in% c("TFT", "TTTD"),  ]
  }
  if(is_LR) {
    dt_clean <- merge.data.table(
      x = unique(dt_clean[, c("L_GENE", "R_GENE", "LR_GENES")]),
      y = dt_clean[,.N, by = "LR_GENES"][N >= filter, ],
      by = "LR_GENES",
      all = FALSE
    )[, LR_GENES := NULL]
    dt_clean[, L_GENE := paste0("L-", L_GENE)]
    dt_clean[, R_GENE := paste0("R-", R_GENE)]
  } else{
    dt_clean <- merge.data.table(
      x = unique(dt_clean[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
      y = dt_clean[,.N, by = "LR_CELLTYPES"][N >= filter, ],
      by = "LR_CELLTYPES",
      all = FALSE
    )[, LR_CELLTYPES := NULL]
    dt_clean[, L_CELLTYPE := paste0("L-", L_CELLTYPE)]
    dt_clean[, R_CELLTYPE := paste0("R-", R_CELLTYPE)]
  }
  if(!is.null(dir)) {
    png(dir, width = 2000, height = 2000, res = 200)
  }
  circos.clear()
  chordDiagram(dt_clean, annotationTrack = "grid", 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(dt_clean))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.2)
  }, bg.border = NA)
  
  if(!is.null(dir)) {
    dev.off()
  }
}


###############
#############
#########

test_lung_detected <- test_lung[SIG_TYPE %in% c("FTT", "TFT", "TTF", "TTTD", "TTTU"), ]
test_lung_sig <- test_lung[SIG_TYPE %in% c("FTT", "TFT", "TTTD", "TTTU"),]
colnames(test_lung_detected)

df_d <- setDF(test_lung_detected)
df_sig <- setDF(test_lung_sig)


ora_test <- analyze_overrepresentation(df_d, df_sig,
                                       exclude_tissue_overrepresentation=TRUE,
                                       adjust_pvals=FALSE)

orat_test_LR_GENES <- ora_test$L_GENE


#########

as.data.frame.matrix(table(diffcom_results_detected$tms_facs$SIG_TYPE))

#detected interaction by LR_CELLTYPES
table(diffcom_results_detected$tms_facs[TISSUE == "Liver",]$SIG_TYPE)

diffcom_results_detected$tms_facs[, LR_T_CELLTYPES := paste(TISSUE, L_CELLTYPE, TISSUE, R_CELLTYPE, sep = "_" ) ]
distr_CCI_tms_facs <- dcast(diffcom_results_detected$tms_facs[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                            TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_CCI_tms_facs[, N := rowSums(.SD), .SDcols = colnames(distr_CCI_tms_facs)[-c(1,2)]]
distr_CCI_tms_facs[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_facs[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_facs[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_CCI_tms_facs[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_CCI_tms_facs[, appear := old - young]
distr_CCI_tms_facs[, net_reg := up - down]
hist(distr_CCI_tms_facs$appear, breaks = 50)
hist(distr_CCI_tms_facs$net_reg, breaks = 50)
hist(distr_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$appear, breaks = 50)
hist(distr_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$net_reg, breaks = 50)
ggplot(distr_CCI_tms_facs, aes(x = appear, y = net_reg)) + geom_point()

het_tms_facs <- merge.data.table(
  diffcom_results_detected$tms_facs[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                       .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                       by = "LR_T_CELLTYPES" ],
  diffcom_results_detected$tms_facs[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                       .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                       by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_tms_facs[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_tms_facs$diff_LR, breaks = 100)
hist(het_tms_facs$diff_L, breaks = 50)
hist(het_tms_facs$diff_R, breaks = 50)



diffcom_results_detected$tms_droplet[, LR_T_CELLTYPES := paste(TISSUE, L_CELLTYPE, TISSUE, R_CELLTYPE, sep = "_" ) ]
distr_CCI_tms_droplet <- dcast(diffcom_results_detected$tms_droplet[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                            TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_CCI_tms_droplet[, N := rowSums(.SD), .SDcols = colnames(distr_CCI_tms_droplet)[-c(1,2)]]
distr_CCI_tms_droplet[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_droplet[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_droplet[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_CCI_tms_droplet[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_CCI_tms_droplet[, appear := old - young]
distr_CCI_tms_droplet[, net_reg := up - down]
hist(distr_CCI_tms_droplet[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$appear, breaks = 50)
hist(distr_CCI_tms_droplet[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$net_reg, breaks = 50)
hist(distr_CCI_tms_droplet$appear, breaks = 50)
hist(distr_CCI_tms_droplet$net_reg, breaks = 50)
ggplot(distr_CCI_tms_droplet, aes(x = appear, y = net_reg)) + geom_point()

het_tms_droplet <- merge.data.table(
  diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_tms_droplet[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_tms_droplet$diff_LR, breaks = 100)
hist(het_tms_droplet$diff_L, breaks = 50)
hist(het_tms_droplet$diff_R, breaks = 50)

diffcom_results_detected$calico[, LR_T_CELLTYPES := paste(TISSUE, L_CELLTYPE, TISSUE, R_CELLTYPE, sep = "_" ) ]
distr_CCI_calico <- dcast(diffcom_results_detected$calico[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                            TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_CCI_calico[, N := rowSums(.SD), .SDcols = colnames(distr_CCI_calico)[-c(1,2)]]
distr_CCI_calico[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_CCI_calico[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_CCI_calico[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_CCI_calico[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_CCI_calico[, appear := old - young]
distr_CCI_calico[, net_reg := up - down]
hist(distr_CCI_calico$appear, breaks = 50)
hist(distr_CCI_calico$net_reg, breaks = 50)
ggplot(distr_CCI_calico, aes(x = appear, y = net_reg)) + geom_point()

het_calico <- merge.data.table(
  diffcom_results_detected$calico[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  diffcom_results_detected$calico[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_calico[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_calico$diff_LR, breaks = 100)
hist(het_calico$diff_L, breaks = 50)
hist(het_calico$diff_R, breaks = 50)


hist(distr_CCI_tms_facs$N, breaks = 50)
hist(distr_CCI_tms_droplet$N, breaks = 50)
hist(distr_CCI_calico$N, breaks = 30)

#maybe not all genes in calico!!!
hist(distr_CCI_tms_facs[TISSUE == "Kidney", ]$N, breaks = 50)
hist(distr_CCI_tms_droplet[TISSUE == "Kidney", ]$N, breaks = 50)
hist(distr_CCI_calico[TISSUE == "kidney", ]$N, breaks = 50)

#difference in detection
test_lung_sig[SIG_TYPE %in% c("TTTU", "TTTD", "TTF", "TFT"),]
test_lung_sig[SIG_TYPE %in% c("TTTU", "TTTD", "TTF", "FTT"),]


distr_LR_tms_facs <- dcast(diffcom_results_detected$tms_facs[!(SIG_TYPE %in% "FFF"), c("LR_GENES", "SIG_TYPE", "TISSUE") ],
                            LR_GENES ~ SIG_TYPE, value.var = "SIG_TYPE")

distr_LR_tiss_tms_facs <- merge.data.table(
  dcast(diffcom_results_detected$tms_facs[!(SIG_TYPE %in% "FFF"), c("LR_GENES", "SIG_TYPE", "TISSUE") ],
        TISSUE + LR_GENES ~ SIG_TYPE, value.var = "SIG_TYPE"),
  diffcom_results$tms_facs[ , length(unique(LR_CELLTYPES)), by = "TISSUE"],
  by = "TISSUE",
  all.x =  TRUE
)
cols <- c("FTT", "TFT", "TTF", "TTTD", "TTTU")
distr_LR_tiss_tms_facs[, N := rowSums(.SD) , .SDcols = cols]

distr_LR_tiss_tms_facs[, paste0("pct_", c(cols, "N")) := .SD/V1 , .SDcols = c(cols, "N")]

####################



ggplot(diffcom_results$tms_facs, aes(x = LR_SCORE_young)) + geom_histogram(bins = 100) + scale_y_log10()
ggplot(diffcom_results$tms_facs, aes(x = LR_SCORE_old)) + geom_histogram(bins = 150) + scale_y_log10()
ggplot(diffcom_results$tms_facs, aes(x = LR_LOGFC)) + geom_histogram(bins = 150) + scale_y_log10()

ggplot(diffcom_results$tms_facs[LR_DETECTED_young == TRUE,], aes(x = LR_SCORE_young)) + 
  geom_histogram(bins = 100) + scale_y_log10()
ggplot(diffcom_results$tms_facs[LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05,], aes(x = LR_SCORE_old)) + 
  geom_histogram(bins = 100) + scale_y_log10()
ggplot(diffcom_results$tms_facs[(LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05) |
                                  (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05),], aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_facs[(LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 0.85) |
                                  (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 0.85),],
       aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_facs[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 0.85) |
                                  (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 0.85)) &
                                  BH_PVAL_DIFF <= 0.05,],
       aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$calico[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 1) |
                                   (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 1)) &
                                  BH_PVAL_DIFF <= 0.05,],
       aes(x = LR_LOGFC)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_droplet[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 1) |
                                 (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 1)) &
                                BH_PVAL_DIFF <= 0.05,],
       aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_facs[LR_DETECTED_old == TRUE & BH_PVAL_old | LR_DETECTED_young == TRUE ,], aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()



quantile(diffcom_results$tms_facs[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 0.85) |
                                     (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 0.85)) &
                                    BH_PVAL_DIFF <= 0.05,]$LR_LOGFC_ABS, probs = seq(0,1, 0.1))


########
##non immune cells
diffcom_results_detected$tms_facs[, c("L_T_CELLTYPE", "R_T_CELLTYPE") := .(paste(TISSUE, L_CELLTYPE, sep = "_"), paste(TISSUE, R_CELLTYPE, sep = "_"))]

write.csv(unique(diffcom_results_detected$tms_facs[, .(L_T_CELLTYPE)]),
          file = "../../../../../facs_t_cell_type.csv",
          row.names = FALSE
)

identical(
  sort(unique(diffcom_results_detected$tms_facs[, L_T_CELLTYPE])),
  sort(unique(diffcom_results_detected$tms_facs[, R_T_CELLTYPE]))
)

cell_types_immune_tms_facs <- read.csv("../../../../../facs_t_cell_type_sorted.csv",
                                       header = TRUE, 
                                       stringsAsFactors = FALSE)
table(cell_types_immune_tms_facs$Immune)
anyNA(cell_types_immune_tms_facs)

setDT(cell_types_immune_tms_facs)
ct_nonIm_tms_facs <- cell_types_immune_tms_facs[Immune == FALSE, L_T_CELLTYPE]

tms_facs_nonIm <- diffcom_results_detected$tms_facs[L_T_CELLTYPE %in% ct_nonIm_tms_facs & R_T_CELLTYPE %in% ct_nonIm_tms_facs, ]

distr_nonIm_CCI_tms_facs <- dcast(tms_facs_nonIm[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                            TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_nonIm_CCI_tms_facs[, N := rowSums(.SD), .SDcols = colnames(distr_nonIm_CCI_tms_facs)[-c(1,2)]]
distr_nonIm_CCI_tms_facs[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_nonIm_CCI_tms_facs[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_nonIm_CCI_tms_facs[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_nonIm_CCI_tms_facs[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_nonIm_CCI_tms_facs[, appear := old - young]
distr_nonIm_CCI_tms_facs[, net_reg := up - down]
hist(distr_nonIm_CCI_tms_facs$appear, breaks = 50)
hist(distr_nonIm_CCI_tms_facs$net_reg, breaks = 50)
hist(distr_nonIm_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$appear, breaks = 50)
hist(distr_nonIm_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$net_reg, breaks = 50)
ggplot(distr_nonIm_CCI_tms_facs, aes(x = appear, y = net_reg)) + geom_point()

het_nonIm_tms_facs <- merge.data.table(
  tms_facs_nonIm[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  tms_facs_nonIm[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_nonIm_tms_facs[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_nonIm_tms_facs$diff_LR, breaks = 100)
hist(het_nonIm_tms_facs$diff_L, breaks = 50)
hist(het_nonIm_tms_facs$diff_R, breaks = 50)

################

test_lung <- test_d[TISSUE == "Lung",]
test_lung <- diffcom_results_detected$tms_facs[TISSUE == "Liver",]
test_lung <- tms_facs_nonIm[TISSUE == "Lung",]

test_lung[SIG_TYPE %in% c("TFT", "FTT", "TTTU", "TTTD") & L_GENE == "Ccl5" & R_GENE == "Cxcr6",]
test_lung[SIG_TYPE %in% c("TFT", "FTT", "TTTU", "TTTD") & L_GENE == "Mrc1" & R_GENE == "Ptprc" & L_CELLTYPE == "hepatocyte",]

diffcom_results_detected$tms_facs[LR_GENES == "Ltb_Ltbr" & LR_CELLTYPES ==  "B cell_Kupffer cell"]
test_lung[LR_GENES == "Mrc1_Ptprc" & LR_CELLTYPES == "endothelial cell of hepatic sinusoid_B cell"]

test_lung <- test_lung[, c("L_GENE", "R_GENE", "L_CELLTYPE", "R_CELLTYPE",
                               "LR_KEEP_young", "LR_KEEP_old", "LR_LOGFC", "LR_CELLTYPES", "LR_GENES", "SIG_TYPE")]
test_lung[, L := paste("L", L_CELLTYPE, L_GENE, sep = "_")]
test_lung[, R := paste("R", R_CELLTYPE, R_GENE, sep = "_")]

test_lung_d <- test_lung[ SIG_TYPE != "FFF", ]
test_lung_up <- test_lung[SIG_TYPE %in% c("FTT", "TTTU"),  ]
test_lung_down <- test_lung[SIG_TYPE %in% c("TFT", "TTTD"),  ]

test_lung_d_up <- merge.data.table(test_lung_d[, .N, by = LR_CELLTYPES],
                 test_lung_up[, .N, by = LR_CELLTYPES],
                 by = "LR_CELLTYPES", 
                 all.x = TRUE
)
test_lung_d_up[is.na(test_lung_d_up)] <- 0
test_lung_d_up[, ratio := N.y/N.x]

test_lung_d_down <- merge.data.table(test_lung_d[, .N, by = LR_CELLTYPES],
                                   test_lung_down[, .N, by = LR_CELLTYPES],
                                   by = "LR_CELLTYPES", 
                                   all.x = TRUE
)
test_lung_d_down[is.na(test_lung_d_down)] <- 0
test_lung_d_down[, ratio := N.y/N.x]



test_lung_chord <- data.table::merge.data.table(y = test_lung_d[, .N, by = LR_CELLTYPES],
                                                   x = unique(test_lung_d[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
                                                   by = "LR_CELLTYPES",
                                                   all = FALSE)[, LR_CELLTYPES := NULL]
test_lung_chord[, L_CELLTYPE := paste0("L_", L_CELLTYPE)]
test_lung_chord[, R_CELLTYPE := paste0("R_", R_CELLTYPE)]
chordDiagram(test_lung_chord, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()




test_lung_up_chord <- data.table::merge.data.table(y = test_lung_d_up[, c("LR_CELLTYPES", "ratio")],
                             x = unique(test_lung_up[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
                             by = "LR_CELLTYPES",
                             all = FALSE)[, LR_CELLTYPES := NULL]
test_lung_up_chord[, L_CELLTYPE := paste0("L_", L_CELLTYPE)]
test_lung_up_chord[, R_CELLTYPE := paste0("R_", R_CELLTYPE)]
chordDiagram(test_lung_up_chord, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_up_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

test_lung_down_chord <- data.table::merge.data.table(y = test_lung_d_down[, c("LR_CELLTYPES", "ratio")],
                                                   x = unique(test_lung_down[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
                                                   by = "LR_CELLTYPES",
                                                   all = FALSE)[, LR_CELLTYPES := NULL]
test_lung_down_chord[, L_CELLTYPE := paste0("L_", L_CELLTYPE)]
test_lung_down_chord[, R_CELLTYPE := paste0("R_", R_CELLTYPE)]
chordDiagram(test_lung_down_chord, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_down_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

test_lung_up_chord_g <- data.table::merge.data.table(y = test_lung_up[, .N, by = LR_GENES],
                                                           x = unique(test_lung_up[, c("L_GENE", "R_GENE", "LR_GENES")]),
                                                           by = "LR_GENES",
                                                           all = FALSE)[, LR_GENES := NULL]
test_lung_up_chord_g[, L_GENE := paste0("L_", L_GENE)]
test_lung_up_chord_g[, R_GENE := paste0("R_", R_GENE)]
chordDiagram(test_lung_up_chord_g, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_up_chord_g))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

test_lung_down_chord_g <- data.table::merge.data.table(y = test_lung_down[, .N, by = LR_GENES][N > 7, ],
                                                     x = unique(test_lung_down[, c("L_GENE", "R_GENE", "LR_GENES")]),
                                                     by = "LR_GENES",
                                                     all = FALSE)[, LR_GENES := NULL]
test_lung_down_chord_g[, L_GENE := paste0("L_", L_GENE)]
test_lung_down_chord_g[, R_GENE := paste0("R_", R_GENE)]
chordDiagram(test_lung_down_chord_g, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_down_chord_g))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()






#try with circos
#test_lung_d_up <- test_lung_d[LR_UP == TRUE, c("L", "R") ]
#test_lung_d_down <- test_lung_d[LR_DOWN == TRUE, c("L", "R") ]

factor_total <- factor(c(sort(unique(test_lung_up$L)), sort(unique(test_lung_up$R))))
factor_LR <- factor(c(rep("L", length(unique(test_lung_up$L))), rep("R", length(unique(test_lung_up$R))) ))
factor_col_LR <- c(rep(rand_color(1), length(unique(test_lung_up$L))), rep(rand_color(1), length(unique(test_lung_up$R))) )

circos.par(
  cell.padding = c(0.00, 0, 0.00, 0),
  gap.degree = 0.,
  start.degree = 270,
  track.height = 0.01)
circos.initialize(factors = factor_total, xlim = c(0,1) )
#circos.track(ylim = c(0,0.1))
#circos.track(ylim = c(0,0.1))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ycenter,
              #CELL_META$sector.index,
              factor_LR[CELL_META$sector.numeric.index],
              col = factor_col_LR[CELL_META$sector.numeric.index]
              #facing = "inside", niceFacing = TRUE
  )
})
for(i in 1:nrow(test_lung_up)) {
  circos.link(test_lung_up$L[[i]], 0, test_lung_up$R[[i]], 0,
              col = "red",
              directional = 0)
}
for(i in 1:nrow(test_lung_down)) {
  circos.link(test_lung_down$L[[i]], 0, test_lung_down$R[[i]], 0,
              col = "blue",
              directional = 1)
}
circos.clear()





########enrichment####
setDT(test_lung_sig)
up_genes <- c(unique(test_lung_sig[SIG_TYPE %in% c("FTT", "TTTU"),]$L_GENE),
  unique(test_lung_sig[SIG_TYPE %in% c("FTT", "TTTU"),]$R_GENE))
down_genes <- c(unique(test_lung_sig[SIG_TYPE %in% c("TFT", "TTTD"),]$L_GENE),
              unique(test_lung_sig[SIG_TYPE %in% c("TFT", "TTTD"),]$R_GENE))
intersect(up_genes, down_genes)

library(clusterProfiler)
library(org.Mm.eg.db)




test_lung_diffnet <- test_lung[, lapply(.SD, sum), by = LR_CELLTYPES, .SDcols = c("LR_UP", "LR_DOWN")]

test_lung_net <- test_lung[, lapply(.SD, sum), by = LR_CELLTYPES, .SDcols = c("LR_KEEP_young", "LR_KEEP_old", "LR_KEEP_DIFF")]


test_lung[LR_CELLTYPES == "B cell_dendritic cell" & (LR_KEEP_old == TRUE | LR_KEEP_young == TRUE),]


table(test$LR_KEEP_DIFF & test$LR_LOGFC > 0)
table(test$LR_KEEP_DIFF & test$LR_LOGFC < 0)
table(test$LR_KEEP_young)
table(test$LR_KEEP_old)





sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC > 0, ]$L_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC > 0, ]$R_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC > 0, ]$LR_GENES), decreasing = TRUE)

sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC < 0, ]$L_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC < 0, ]$R_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC < 0, ]$LR_GENES), decreasing = TRUE)

#Apoe
test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]
test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]

sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]$TISSUE), decreasing = TRUE)
sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]$TISSUE), decreasing = TRUE)

sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]$L_CELLTYPE), decreasing = TRUE)
sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]$L_CELLTYPE), decreasing = TRUE)

sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]$R_GENE), decreasing = TRUE)
sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]$R_GENE), decreasing = TRUE)


#############
