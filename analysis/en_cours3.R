#######

seurat_facs <- readRDS("seurat_final_tms_facs.rds")
data_facs <- GetAssayData(seurat_facs, assay = "RNA", slot = "data")
data_facs <- expm1(data_facs)
data_tr_facs <- Matrix::t(data_facs)
colsum_facs <- Matrix::colSums(data_facs)

seurat_facs$age_group <- ifelse(seurat_facs$age %in% c("1m", "3m"), "YOUNG", "OLD")
md_facs <- seurat_facs[[]]
setDT(md_facs)
md_facs[, grouping := paste(tissue, cell_ontology_final, age_group, sex)]

gene_sum_facs <- DelayedArray::rowsum(data_tr_facs, md_facs$grouping)
gene_mean_facs <- gene_sum_facs/as.vector(table(md_facs$grouping))

data_young_facs <- GetAssayData(subset(seurat_facs, subset = age_group == "YOUNG"), assay = "RNA", slot = "data")
data_old_facs <- GetAssayData(subset(seurat_facs, subset = age_group == "OLD"), assay = "RNA", slot = "data")

data_tr_young_facs <- Matrix::t(data_young_facs)
data_tr_old_facs <- Matrix::t(data_old_facs)

md_young_facs <- subset(seurat_facs, subset = age_group == "YOUNG")[[]]
md_old_facs <- subset(seurat_facs, subset = age_group == "OLD")[[]]
setDT(md_young_facs)
setDT(md_old_facs)

md_young_facs[, tissue_celltype := paste(tissue, cell_ontology_final, sep = "_")]
md_old_facs[, tissue_celltype := paste(tissue, cell_ontology_final, sep = "_")]

gene_sum_young_facs <- DelayedArray::rowsum(data_tr_young_facs, md_young_facs$tissue_celltype)
gene_sum_old_facs <- DelayedArray::rowsum(data_tr_old_facs, md_old_facs$tissue_celltype)

gene_mean_young_facs <- gene_sum_young_facs/as.vector(table(md_young_facs$tissue_celltype))
gene_mean_old_facs <- gene_sum_old_facs/as.vector(table(md_old_facs$tissue_celltype))








#######

facs_mean_expr <- readRDS("../data_scAgeCom/analysis/outputs_data/facs_mean_expression.rds")
facs_mean_expr_young_female <- facs_mean_expr[grepl(":YOUNG:female", rownames(facs_mean_expr), fixed = TRUE),]
facs_mean_expr_old_female <- facs_mean_expr[grepl(":OLD:female", rownames(facs_mean_expr), fixed = TRUE),]
facs_mean_expr_young_male <- facs_mean_expr[grepl(":YOUNG:male", rownames(facs_mean_expr), fixed = TRUE),]
facs_mean_expr_old_male <- facs_mean_expr[grepl(":OLD:male", rownames(facs_mean_expr), fixed = TRUE),]

rownames(facs_mean_expr_young_female) <- sub(":YOUNG:female", "", rownames(facs_mean_expr_young_female), fixed = TRUE)
rownames(facs_mean_expr_old_female) <- sub(":OLD:female", "", rownames(facs_mean_expr_old_female), fixed = TRUE)
rownames(facs_mean_expr_young_male) <- sub(":YOUNG:male", "", rownames(facs_mean_expr_young_male), fixed = TRUE)
rownames(facs_mean_expr_old_male) <- sub(":OLD:male", "", rownames(facs_mean_expr_old_male), fixed = TRUE)

keep_female <- intersect(rownames(facs_mean_expr_young_female), rownames(facs_mean_expr_old_female))
keep_male <- intersect(rownames(facs_mean_expr_young_male), rownames(facs_mean_expr_old_male))

keep_mixed <- intersect(keep_female, keep_male)

logfc_female <- log(facs_mean_expr_old_female[keep_female,]/facs_mean_expr_young_female[keep_female,])
logfc_male <- log(facs_mean_expr_old_male[keep_male,]/facs_mean_expr_young_male[keep_male,])

dt_logfc_female <- as.data.table(as.table(logfc_female))
dt_logfc_male <- as.data.table(as.table(logfc_male))

dt_logfc_female[, N := ifelse(is.nan(N), 0, N)]
dt_logfc_female[, N := ifelse(is.infinite(N) & N > 0, 20, N)]
dt_logfc_female[, N := ifelse(is.infinite(N) & N < 0, -20, N)]

dt_logfc_male[, N := ifelse(is.nan(N), 0, N)]
dt_logfc_male[, N := ifelse(is.infinite(N) & N > 0, 20, N)]
dt_logfc_male[, N := ifelse(is.infinite(N) & N < 0, -20, N)]

summary(dt_logfc_female$N)
summary(dt_logfc_male$N)

dt_logfc <- rbindlist(
  list(female = dt_logfc_female, male = dt_logfc_male),
  idcol = "gender"
)

ggplot(dt_logfc, aes(N, color = gender, fill = gender)) +
  geom_histogram(bins = 100, alpha = 0.4, position = "identity") +
  scale_y_log10()

LRdb_genes <- unique(unlist(LRdb_mouse$LRdb_curated[, c(2,3,4,5,6)]))
dt_logfc_LR <- dt_logfc[V2 %in% LRdb_genes]

dt_logfc_LR_wide <- dcast.data.table(
  dt_logfc_LR,
  formula = V1 + V2 ~ gender,
  value.var = "N"
)

ggplot(dt_logfc_LR_wide, aes(female, male)) + geom_point() + geom_smooth(method = "lm") +
  geom_density2d()

ggplot(dt_logfc_LR, aes(N, color = gender, fill = gender)) +
  geom_histogram(aes(y = ..density..), bins = 100, alpha = 0.4, position = "identity")

logfc_female_2 <- log(facs_mean_expr_old_female[keep_mixed,]/facs_mean_expr_young_female[keep_mixed,])
logfc_male_2 <- log(facs_mean_expr_old_male[keep_mixed,]/facs_mean_expr_young_male[keep_mixed,])

hist(as.vector(logfc_female), breaks = 100)
hist(as.vector(logfc_male), breaks = 100)



plot(as.vector(logfc_female_2), as.vector(logfc_male_2))

plot(logfc_female_2[grepl("Lung|Bladder", rownames(logfc_female_2)), "App"], logfc_male_2[grepl("Lung|Bladder", rownames(logfc_female_2)), "App"])
abline(h = 0)

identical(rownames(logfc_female_2), rownames(logfc_male_2))

summary(as.vector(logfc_female))
summary(as.vector(logfc_male))

facs_ids <- strsplit(rownames(facs_mean_expr), ":")

dt_template <- data.table(
  id = rownames(facs_mean_expr),
  tissue = sapply(facs_ids, function(i) i[[1]]),
  celltype = sapply(facs_ids, function(i) i[[2]]),
  age_group = sapply(facs_ids, function(i) i[[3]]),
  gender = sapply(facs_ids, function(i) i[[4]])
)




