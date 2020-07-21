

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

