
## Messy file with some ideas
## will eventually be removed#####

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


B2M <- LRone2one[GENESYMB_L == "B2m",]


#setDT(diffcom_results_detected_strong$tms_facs)
B2M_facs <- diffcom_results_detected_strong$tms_facs[L_GENE == "B2m", ]
B2M_droplet <- diffcom_results_detected_strong$tms_droplet[L_GENE == "B2m", ]
table(B2M_facs$SIG_TYPE)

B2M_facs_up <- B2M_facs[SIG_TYPE %in% c("FTT", "TTTU")]
B2M_facs_up[, TISSUE_LR_CELLTYPES := paste(TISSUE, LR_CELLTYPES, sep = "_")]
B2M_droplet_up <- B2M_droplet[SIG_TYPE %in% c("FTT", "TTTU")]


sort(table(B2M_facs_up$TISSUE), decreasing = TRUE)
sort(table(B2M_droplet_up$TISSUE), decreasing = TRUE)
head(sort(table(B2M_facs_up$TISSUE_LR_CELLTYPES), decreasing = TRUE),10)
sort(table(B2M_facs_up$L_CELLTYPE), decreasing = TRUE)
sort(table(B2M_facs_up$R_GENE), decreasing = TRUE)
ggplot(B2M_facs_up, aes(x= LR_SCORE_old)) + geom_histogram(bins = 50)
ggplot(B2M_facs_up, aes(x= LR_SCORE_young)) + geom_histogram(bins = 50)
ggplot(B2M_facs_up, aes(x= LR_LOGFC)) + geom_histogram(bins = 50)
    


App_nc <- diffcom_results_significant$tms_facs[LR_GENES == "App_Ncstn",]
sort(table(diffcom_results_significant$tms_facs[L_GENE == "App",]$R_GENE))


diffcom_results_detected_strong$tms_facs[L_GENE == "Fgf2",]
test <- diffcom_results_significant$tms_facs[R_GENE == "Fgfr2",]

Itgb1 <- diffcom_results_significant$tms_facs[R_GENE == "Itgb1",]
table(Itgb1$SIG_TYPE)
sort(table(Itgb1$TISSUE), decreasing = TRUE)
sort(table(Itgb1$L_GENE), decreasing = TRUE)
sort(table(Itgb1$L_CELLTYPE), decreasing = TRUE)

Itgb1_drop <- diffcom_results_significant$tms_droplet[R_GENE == "Itgb1",]
table(Itgb1_drop$SIG_TYPE)
sort(table(Itgb1_drop$TISSUE), decreasing = TRUE)
sort(table(Itgb1$L_GENE), decreasing = TRUE)
sort(table(Itgb1$L_CELLTYPE), decreasing = TRUE)


Spp1_facs <- diffcom_results_significant$tms_facs[L_GENE == "Spp1",]
Spp1_drop <- diffcom_results_significant$tms_droplet[L_GENE == "Spp1",]

Itgav_drop <- diffcom_results_significant$tms_droplet[R_GENE == "Itgav",]
table(Itgav_drop$SIG_TYPE)

"Il6" %in% diffcom_results$tms_facs$L_GENE

"Il6" %in% LRone2one$GENESYMB_L

Il6 <- LRone2one[GENESYMB_L == "Il6",]
LRone2one[GENESYMB_R == "Il6ra",]

library(SingleCellSignalR)
data("LRdb")
LRdb
"IL6" %in% LRdb$ligand

LRdb[LRdb$ligand == "IL6",]

data()

IL6sc <- LRdb[LRdb$ligand ==  "IL6",]

LRall$LRall_many2many[LRall$LRall_many2many$GENESYMB_R == "Il6ra", ]


mm2Hs[mm2Hs$`Gene name` == "IL6R",]

test <- diffcom_results_significant$tms_facs
test2 <- diffcom_results_significant$tms_droplet



##volcano plot #####
library(ggrepel)
diffcom_1000iter <- readRDS("test_and_comparison/data_results_diffcom_10000iter_log.rds")
diffcom_1000iter[, LR_DIFF := LR_SCORE_old - LR_SCORE_young]
ggplot(diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old), ], 
       aes(x = LR_DIFF, y = -log10(BH_PVAL_DIFF + 1E-4))) +
  geom_point() +
  geom_vline(xintercept = log(1.1)) +
  geom_vline(xintercept = -log(1.1)) +
  geom_hline(yintercept = -log10(0.05))

+
  geom_text(
    data = diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old) & BH_PVAL_DIFF < 0.01 & abs(LR_DIFF) > 0.5, ],
    aes(label = LR_GENES),
    size = 5
    #box.padding = unit(0.35, "lines"),
    #point.padding = unit(0.3, "lines")
  )

hist(log10(diffcom_1000iter$LR_SCORE_old), breaks = 100)
hist(log10(diffcom_1000iter$LR_SCORE_young), breaks = 100)

ggplot(diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old) & (BH_PVAL_young <= 0.05 | BH_PVAL_old <= 0.05) & 
                          (LR_SCORE_young > 0.1 | LR_SCORE_old > 0.1), ], 
       aes(x = LR_DIFF, y = -log10(BH_PVAL_DIFF + 1E-4))) +
  geom_point() +
  geom_vline(xintercept = log(1.1)) +
  geom_vline(xintercept = -log(1.1)) +
  geom_hline(yintercept = -log10(0.05))


+
  geom_text(
    data = diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old) & BH_PVAL_DIFF < 0.01 & abs(LR_DIFF) > 0.3, ],
    aes(label = L_CELLTYPE),
    size = 5
    #box.padding = unit(0.35, "lines"),
    #point.padding = unit(0.3, "lines")
  )






calico_previous <- bind_tissues(data_path$calico, tissue_list$calico)

calico_new <- bind_tissues("../../../../../scDiffCom_results/diffcom_calico_size_factor_log_10000iter_mixed", tissue_list$calico)

table(calico_new$LR_DETECTED_old, calico_new$LR_DETECTED_young)
table(calico_previous$LR_DETECTED_old, calico_previous$LR_DETECTED_young)


diffcom_results <- mapply(bind_tissues, data_path, tissue_list, SIMPLIFY = FALSE)

facs_previous <- bind_tissues(data_path$tms_facs, tissue_list$tms_facs)
facs_news <- bind_tissues("../../../../../scDiffCom_results/diffcom_tms_facs_size_factor_log_10000iter_mixed", tissue_list$tms_facs)

facs_previous[, id_col := paste(LR_GENES, L_CELLTYPE, R_CELLTYPE, TISSUE, sep = "_")]
facs_news[, id_col := paste(LR_GENES, L_CELLTYPE, R_CELLTYPE, TISSUE, sep = "_")]


table(facs_news$LR_DETECTED_old, facs_news$LR_DETECTED_young)
table(facs_previous$LR_DETECTED_old, facs_previous$LR_DETECTED_young)


facs_comp <- facs_news[id_col %in% facs_previous$id_col,]
facs_previous <- unique(facs_previous)


table(facs_comp$LR_DETECTED_old, facs_comp$LR_DETECTED_young)
table(facs_previous$LR_DETECTED_old, facs_previous$LR_DETECTED_young)



library(config)
library(conflicted)
library(circlize)
library(ggplot2)

library(clusterProfiler)
library(org.Mm.eg.db)



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



#####


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



ora_facs <- ora_results$tms_facs
ora_LR_facs_all <- ora_facs[Tissue == "All" & Category == "LR_GENES"]
ora_LR_facs_all_up <- ora_LR_facs_all[pval_adjusted_UP <= 0.05 & OR_UP >= 1, c("Value", "pval_adjusted_UP", "OR_UP")]
ora_LR_facs_all_down <- ora_LR_facs_all[pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1, c("Value", "pval_adjusted_DOWN", "OR_DOWN")]


ora_LR_facs_tiss <- ora_facs[Tissue != "All" & Category == "LR_GENES"]
ora_LR_facs_tiss_up <- ora_LR_facs_tiss[pval_adjusted_UP <= 0.5 & OR_UP >= 1, c("Tissue", "Value", "pval_adjusted_UP", "OR_UP")]
ora_LR_facs_tiss_down <- ora_LR_facs_tiss[pval_adjusted_DOWN <= 0.5 & OR_DOWN >= 1, c("Tissue", "Value", "pval_adjusted_DOWN", "OR_DOWN")]


test2 <- test[Tissue == "All" & Category == "LR_GENES" & pval_adjusted_UP <= 0.05 & OR_UP >= 1]
test3 <- test[Tissue == "All" & Category == "LR_GENES" & pval_adjusted_UP <= 0.05 & OR_UP <= 1]

test4 <- test[Category == "LR_GENES" & pval_adjusted_UP <= 0.05 & OR_UP >= 1 & Tissue != "All"]

test5 <- test[Category == "LR_CELLTYPES"  & Tissue == "Liver"]


fpm_test <- fpm_results$tms_facs$sub1
fpm_test2 <- fpm_results$tms_facs$sub2


table(datasets_filtered$tms_facs$CASE_TYPE)
table(datasets_filtered$tms_droplet$CASE_TYPE)
table(datasets_filtered$calico$CASE_TYPE)

cluster_test <- datasets_filtered$tms_droplet[ TISSUE == "Thymus"]
cluster_test[, SIG_VAL := ifelse(CASE_TYPE %in% c("TTFU", "TTFD"),
                                 0,
                                 ifelse(CASE_TYPE %in% c("TTTD", "TFTD"),
                                        -1,
                                        1))]
cluster_test <- dcast(cluster_test[, c("LR_CELLTYPES", "LR_GENES", "SIG_VAL")],
                      formula = LR_GENES ~ LR_CELLTYPES,
                      value.var = "SIG_VAL")
cluster_test[is.na(cluster_test)] <- 0

mat_cluster <- as.matrix(cluster_test, rownames = 1)
mat_cluster <- mat_cluster[rowSums(abs(mat_cluster)) > 4,]
mat_cluster <- mat_cluster[, colSums(abs(mat_cluster)) > 3]

pheatmap(
  mat_cluster,
  show_rownames = TRUE,
  show_colnames = FALSE
)

pheatmap(mat_cluster, #cellwidth=10, cellheight=10,
         width=11, height=11,
         #color = colorRampPalette(brewer.pal(8, "RdYlBu"))(length(breaks)-1),
         #breaks = breaks,
         na_col = "green",
         scale = "none",
         #cluster_rows = nrows >= 2,
         #cluster_cols = ncols >= 2,
         #cutree_rows = ifelse(nrows >= 5, min(nrows, 3), 1),
         #cutree_cols = ifelse(ncols >= 5, min(ncols, 3), 1),
         display_numbers=FALSE, 
         fontsize=10, 
         #main = title,
         #ylab = "Transmitter cell", 
         #xlab = "Receiver cell",
         legend=TRUE,
         #filename = filename
)

LRall[SYMB_LR == "B2m_Cd3g"]


b2m_facs <- datasets_filtered$tms_facs[L_GENE == "B2m"]

table(datasets_filtered$tms_facs$CASE_TYPE)


ggplot(datasets_filtered$tms_facs, aes(x = LR_LOGFC*log2(exp(1)), y = -log10(BH_PVAL_DIFF + 1E-4), color = TISSUE)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = log2(1.1)) +
  geom_vline(xintercept = -log2(1.1)) +
  geom_vline(xintercept = log2(1.5)) +
  geom_vline(xintercept = -log2(1.5)) +
  xlab(expression(paste(Log[2], "FC"))) +
  ylab(expression(paste(-Log[10], " ", p[BH])))
#geom_density_2d() +

volcano_facs <- cowplot::plot_grid(
  plotlist = lapply(
    unique(datasets_filtered$tms_facs$TISSUE),
    function(tiss) {
      ggplot(datasets_filtered$tms_facs[TISSUE == tiss],
             aes(x = LR_LOGFC*log2(exp(1)), y = -log10(BH_PVAL_DIFF + 1E-4)), color = CASE_TYPE) +
        geom_point() +
        geom_hline(yintercept = -log10(0.05)) +
        geom_vline(xintercept = log2(1.1)) +
        geom_vline(xintercept = -log2(1.1)) +
        geom_vline(xintercept = log2(1.5)) +
        geom_vline(xintercept = -log2(1.5)) +
        xlab(expression(paste(Log[2], "FC"))) +
        ylab(expression(paste(-Log[10], " ", p[BH]))) +
        theme(text=element_text(size=20))
    }
  ),
  ncol = 5,
  align = "v",
  labels = unique(datasets_filtered$tms_facs$TISSUE)
)

volcano_droplet <- cowplot::plot_grid(
  plotlist = lapply(
    unique(datasets_filtered$tms_droplet$TISSUE),
    function(tiss) {
      ggplot(datasets_filtered$tms_droplet[TISSUE == tiss],
             aes(x = LR_LOGFC*log2(exp(1)), y = -log10(BH_PVAL_DIFF + 1E-4)), color = CASE_TYPE) +
        geom_point() +
        geom_hline(yintercept = -log10(0.05)) +
        geom_vline(xintercept = log2(1.1)) +
        geom_vline(xintercept = -log2(1.1)) +
        geom_vline(xintercept = log2(1.5)) +
        geom_vline(xintercept = -log2(1.5)) +
        xlab(expression(paste(Log[2], "FC"))) +
        ylab(expression(paste(-Log[10], " ", p[BH]))) +
        theme(text=element_text(size=20))
    }
  ),
  ncol = 4,
  align = "v",
  labels = unique(datasets_filtered$tms_droplet$TISSUE)
)

ggsave("../../../../../volcano_test_facs.png", plot = volcano_facs, scale =  2.5)

ggplot(datasets_filtered$tms_droplet, aes(x = LR_LOGFC*log10(exp(1)), y = -log10(BH_PVAL_DIFF + 1E-4))) + geom_point() +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = log(1.1)) +
  geom_vline(xintercept = -log(1.1)) +
  geom_vline(xintercept = log(1.5)) +
  geom_vline(xintercept = -log(1.5)) +
  geom_density_2d()



