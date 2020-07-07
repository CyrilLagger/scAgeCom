
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
