library(clusterProfiler)
library(org.Mm.eg.db)
library(scDiffCom)

diffcom_results <- readRDS("../../data_scAgeCom/analysis/analysis_4_data_diffcom_filter_new.rds")

setorder(diffcom_results$tms_facs, -LR_LOGFC)


head(diffcom_results$tms_facs)
test <- unique(diffcom_results$tms_facs[BH_PVAL_DIFF <= 0.05 & LR_LOGFC >= log(1.1)]$L_GENE)

universe_L <- unique(LRall[scsr == TRUE | cpdb == TRUE]$GENESYMB_L)

ego_test <- enrichGO(
  gene = test,
  OrgDb = org.Mm.eg.db,
  keyType = 'SYMBOL',
  ont = "MF",
  universe = universe_L,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff  = 0.05
)

dotplot(ego_test)

##old stuff below, not working!!

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
