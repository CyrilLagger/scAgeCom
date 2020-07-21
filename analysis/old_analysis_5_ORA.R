
##old stuff below to modify and update!!

library(config)

#######
# OVERREPRESENTATION

#load parameters
load_globals_from_config(R_CONFIG_ACTIVE = "default")
conflicted::conflict_prefer(name = "merge", winner = "base")

diffcom_detected <- lapply(
  diffcom_filter,
  function(x) {
    return(setDF(x[SIG_TYPE != "FFF"]))
  }
)
diffcom_significant <- lapply(
  diffcom_filter,
  function(x) {
    return(setDF(x[SIG_TYPE %in% c("FTT", "TFT", "TTTU", "TTTD")]))
  }
)


overrep_calico <- analyze_overrepresentation(
  diffcom_detected$calico,
  diffcom_significant$calico,
  adjust_pvals = TRUE
)
overrep_calico_tiss <- sapply(unique(diffcom_results$calico$TISSUE), function(tiss) {
  analyze_overrepresentation(
    diffcom_detected$calico[diffcom_detected$calico$TISSUE == tiss,],
    diffcom_significant$calico[diffcom_significant$calico$TISSUE == tiss,],
    exclude_tissue_overrepresentation=TRUE,
    adjust_pvals=TRUE)
}, USE.NAMES = TRUE, simplify = FALSE)

overrep_tms_facs <- analyze_overrepresentation(
  diffcom_detected$tms_facs,
  diffcom_significant$tms_facs,
  adjust_pvals = TRUE
)
overrep_tms_facs_tiss <- sapply(unique(diffcom_detected$tms_facs$TISSUE), function(tiss) {
  analyze_overrepresentation(
    diffcom_detected$tms_facs[diffcom_detected$tms_facs$TISSUE == tiss,],
    diffcom_significant$tms_facs[diffcom_significant$tms_facs$TISSUE == tiss,],
    exclude_tissue_overrepresentation=TRUE,
    adjust_pvals=TRUE)
}, USE.NAMES = TRUE, simplify = FALSE)

overrep_tms_droplet <- analyze_overrepresentation(
  diffcom_detected$tms_droplet,
  diffcom_significant$tms_droplet,
  adjust_pvals = TRUE
)
overrep_tms_droplet_tiss <- sapply(unique(diffcom_detected$tms_droplet$TISSUE), function(tiss) {
  analyze_overrepresentation(
    diffcom_detected$tms_droplet[diffcom_detected$tms_droplet$TISSUE == tiss,],
    diffcom_significant$tms_droplet[diffcom_significant$tms_droplet$TISSUE == tiss,],
    exclude_tissue_overrepresentation=TRUE,
    adjust_pvals=TRUE)
}, USE.NAMES = TRUE, simplify = FALSE)

########################

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



##################################################

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

