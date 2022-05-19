library(Seurat)
library(scDiffCom)
library(data.table)

library(gridExtra)

library(future)
future::plan(future::multisession, workers = 8)

library(ggplot2)

source("utils_random_pairs.R")
source("coexpression.R")

### Input data ###

INPUT = "droplet"  # facs, sample
TISSUE = "Liver"
AGE_GROUP = list(
  YOUNG = c("3m"),
  OLD = c("21m", "24m", "30m")
)
SEX = c("male")

##################


if (INPUT == "droplet") {
  seurat_sample_obj = readRDS("../data/seurat_shared_tms_droplet.rds")
  seurat_sample_obj = subset(x = seurat_sample_obj, subset = tissue == TISSUE)
  seurat_sample_obj[["age_group"]] = ifelse(
    seurat_sample_obj[[]]$age %in% AGE_GROUP$YOUNG, 
    "YOUNG", 
    ifelse(
      seurat_sample_obj[[]]$age %in% AGE_GROUP$OLD, 
      "OLD", 
      "IGNORE"
    )
  )
  seurat_sample_obj = subset(x = seurat_sample_obj, subset = age_group != "IGNORE")
  seurat_sample_obj = subset(x = seurat_sample_obj, subset = sex %in% SEX)
  seurat_sample_obj[["cell_type"]] = seurat_sample_obj$cell_type_scagecom
  seurat_sample_obj[["cell_abbreviation"]] = seurat_sample_obj$cell_abbreviation_scagecom
} else if (INPUT == "facs") {
  stop("Not implemented.")
  seurat_sample_obj = readRDS("../data/seurat_shared_tms_facs.rds")
  seurat_sample_obj = subset(x = seurat_sample_obj, subset = tissue == TISSUE)
} else if (INPUT == "sample") {
  load("../data/seurat_sample_tms_liver.rda")
  seurat_sample_obj = seurat_sample_tms_liver
} else {
  stop("INPUT not recognized.")
}

# seurat_sample_obj
# metadata = seurat_sample_obj[[]]
# View distribution cell types ~ sex ~ age group
# table(metadata[, c("cell_type", "age_group", "sex")])

## Genes in seurat object and LRI
seurat_g = get_genes_from_seurat(seurat_sample_obj)
lri = scDiffCom::LRI_mouse$LRI_curated
lri_simple = subset_simple_lri(lri)
lri_simple_in_seurat = lri_simple[ (LIGAND_1 %in% seurat_g) & (RECEPTOR_1 %in% seurat_g) ]
print_("Seurat genes: ", length(seurat_g))
print_("LRI LR_simple genes: ", length(get_genes_lri(lri_simple_in_seurat, "LR_simple")))
print_("Seurat genes in LRI LR: ", sum(seurat_g %in% get_genes_lri(lri, "LR")))
print_("Seurat genes in LRI LR_simple: ", sum(seurat_g %in% get_genes_lri(lri_simple_in_seurat, "LR_simple")))


results = list()

## Compute results by scDiffCom methodology
name = "scDiffCom"
res_scDiffCom = run_internal_analysis(
  seurat_obj = seurat_sample_obj,
  lri_table = lri_simple_in_seurat
)
res_scDiffCom = FilterCCI(res_scDiffCom, skip_ora = TRUE)
results[[name]] = res_scDiffCom


## Compute results on Seurat object with the gene names shuffled
name = "Seurat shuffled"
res_shuffled = run_internal_analysis(
  seurat_obj = seurat_shuffle(seurat_sample_obj),
  lri_table = lri_simple_in_seurat
)
res_shuffled = FilterCCI(res_shuffled, skip_ora = TRUE)
results[[name]] = res_shuffled


## Compute results on original Seurat object with randomly synthesized LRI
name = "Random LRI"
res_random = run_internal_analysis(
  seurat_obj = seurat_sample_obj,
  lri_table = create_random_LRIs(
    seurat_sample_obj,
    lri_simple_in_seurat,
    random_lri_nrow = NULL)
)
res_random = FilterCCI(res_random, skip_ora = TRUE)
results[[name]] = res_random


## Compute results on pairs made of the same gene in different cell types (L, R or non-LRI gene)
# res_single_gene = run_internal_analysis(
#   seurat_obj = seurat_sample_obj,
#   lri_table = create_random_LRIs(
#     seurat_sample_obj,
#     scDiffCom::LRI_mouse$LRI_curated,
#     random_lri_nrow = 100,
#     build_pairs_from_single_genes = TRUE)
# )
# res_single_gene = FilterCCI(res_single_gene, skip_ora = TRUE)
# results = append(results, res_single_gene)

# res_single_gene_L = run_internal_analysis(
#   seurat_obj = seurat_sample_obj,
#   lri_table = create_random_LRIs(
#     seurat_sample_obj,
#     scDiffCom::LRI_mouse$LRI_curated,
#     random_lri_nrow = 100,
#     build_pairs_from_single_genes = TRUE,
#     build_pairs_from_single_genes_type = "Ligands")
# )
# res_single_gene_L = FilterCCI(res_single_gene_L, skip_ora = TRUE)
# results = append(results, res_single_gene_L)

# On p-value histograms: http://varianceexplained.org/statistics/interpreting-pvalue-histogram/

# qvalue::pi0est() -> percent of null and alternative hypotheses

pdf(file="naive_exploration_plots.pdf")
#
# # Histograms: Raw P_VALUE_DE
for (i in 1:length(results)) {
  name = names(results)[i]
  res = results[[i]]
  hist(res@cci_table_raw$P_VALUE_DE, freq = FALSE, main = paste0(name, ": RAW P_VALUE_DE"))
}

# # Histograms: Detected P_VALUE_DE
for (i in 1:length(results)) {
  name = names(results)[i]
  res = results[[i]]
  hist(res@cci_table_detected$P_VALUE_DE, freq = FALSE, main = paste0(name, ": DETECTED P_VALUE_DE"))
}

# # Histograms: Raw LogFC
for (i in 1:length(results)) {
  name = names(results)[i]
  res = results[[i]]
  hist(res@cci_table_raw$LOGFC, freq = FALSE, main = paste0(name, ": RAW LOGFC"))
}


# # Histograms: Detected LogFC
for (i in 1:length(results)) {
  name = names(results)[i]
  res = results[[i]]
  hist(res@cci_table_detected$LOGFC, freq = FALSE, main = paste0(name, ": DETECTED LOGFC"))
}

# # Histogram: Specificity P_VALUE_YOUNG
for (i in 1:length(results)) {
  name = names(results)[i]
  res = results[[i]]
  barplot(table(res@cci_table_detected$P_VALUE_YOUNG), freq = FALSE, main = paste0(name, ": P_VALUE_YOUNG"))
}

# # Histogram: Specificity P_VALUE_OLD
for (i in 1:length(results)) {
  name = names(results)[i]
  res = results[[i]]
  barplot(table(res@cci_table_detected$P_VALUE_YOUNG), freq = FALSE, main = paste0(name, ": P_VALUE_YOUNG"))
}

# # Counts: Regulation
for (i in 1:length(results)) {
  name = names(results)[i]
  res = results[[i]]
  barplot(table(res@cci_table_detected$REGULATION), freq = FALSE, main = paste0(name, ": counts REGULATION"))
}

dev.off()


## DGE Analysis
dge_markers = dge_by_celltypes(seurat_sample_obj)
dge_markers_sign = dge_markers[p_val_adj < 0.05]
print_("Number of DGE events: ", dim(dge_markers_sign)[1])
num_celltypes = length(unique(dge_markers[, cell_type]))
# table( dge_markers_sign[, .N, by=gene][, N] )

dge_genes = unique(dge_markers_sign$gene)
num_dge_genes = length(dge_genes)
stable_genes = dge_markers[order(abs(avg_log2FC))]$gene[1:num_dge_genes]

# # Simulate drop in dge genes [shuffling]
t = simulate_random_gene_drop(seurat_sample_obj, dge_genes, 50)
dt = data.table(t)
write.csv(dt, "regulation_vs_dge_drop.csv")

# # Simulate drop in non-dge genes [shuffling]
# t = simulate_random_gene_drop(seurat_sample_obj, stable_genes, 50)
# dt = data.table(t)
# write.csv(dt, "regulation_vs_dge_stable_drop.csv")

# # Simulate drop in dge genes [no shuffling]
# t = simulate_random_gene_drop(seurat_sample_obj, dge_genes, 50, with_shuffling = FALSE)
# dt = data.table(t)
# write.csv(dt, "regulation_vs_dge_drop_no_shuffle.csv")

# # Simulate drop in non-dge genes [no shuffling]
# t = simulate_random_gene_drop(seurat_sample_obj, stable_genes, 50, with_shuffling = FALSE)
# dt = data.table(t)
# write.csv(dt, "regulation_vs_dge_stable_drop_no_shuffle.csv")

df_dge = read_drop_simulation_results("dge_drop")
m_dge = lm(SC ~ DROP, data=df_dge)
print(summary(m_dge))
# plot(df_dge$DROP, df_dge$SC)
# ggplot(df_dge,aes(DROP, SC)) +
#   geom_point() +
#   geom_smooth(method='lm')

df_stable = read_drop_simulation_results("stable_drop")
m_stable = lm(SC ~ DROP, data=df_stable)
print(summary(m_stable))
# plot(df_stable$DROP, df_stable$SC)
# ggplot(df_stable,aes(DROP, SC)) +
#   geom_point() +
#   geom_smooth(method='lm')

# Multiple reg
df_dge$DROP_TYPE = "DGE_SHUFFLE"
df_stable$DROP_TYPE = "STABLE_SHUFFLE"
df_comb = rbind(df_dge, df_stable)
m_multi = lm(SC ~ DROP + DROP_TYPE, data=df_comb)
print(summary(m_multi))

df_dge_noshuffle = read_drop_simulation_results("dge_drop_no_shuffle")
# m_dge_noshuffle = lm(SC ~ DROP, data=df_dge_noshuffle)
# print(summary(m_dge_noshuffle))

df_stable_noshuffle = read_drop_simulation_results("stable_drop_no_shuffle")
# m_stable_noshuffle = lm(SC ~ DROP, data=df_stable_noshuffle)
# print(summary(m_stable_noshuffle))


df_dge_noshuffle$DROP_TYPE = "DGE"
df_stable_noshuffle$DROP_TYPE = "STABLE"

df_comb = rbind(df_dge, df_stable, df_dge_noshuffle, df_stable_noshuffle)

ggplot(df_comb, aes(x=DROP, y=SC, color=DROP_TYPE)) +
  geom_point() + 
  geom_smooth(method=lm) +
  ggtitle("# significant changes ~ dropping DGE vs stable genes")