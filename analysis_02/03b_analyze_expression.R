source("utils_random_pairs.R")
library(ggplot2)
library(data.table)

## Check the full data - why only 8 datasets?
df = data.table::data.table(read.csv("data/expression/markers_full.csv"))
df$dataset = paste0(df$input, "_", df$tissue, "_", df$sex)

## Plot volcano (aggregated data)
lri = scDiffCom::LRI_mouse$LRI_curated
lri_genes = get_genes_lri(lri)

df$gene_in_lri = df$gene %in% lri_genes
df$avg_logFC = log(2) * df$avg_log2FC

df_dset = df[
  (p_val_adj < 0.05) & (abs(avg_log2FC) > 1),
  list(
    num_genes_in_lri = sum(gene_in_lri),
    num_genes = .N
  ), 
  by=dataset]
df_dset$lri_ratio = df_dset$num_genes_in_lri / df_dset$num_genes
df_dset = df_dset[order(-lri_ratio)]

# df for stacked barplot
df_lri = copy(df_dset)
df_lri[, num_genes := num_genes_in_lri]
df_lri[, gene_type := "LRI"]
df_nonlri = copy(df_dset)
df_nonlri[, gene_type := "NONLRI"]
df_g = rbind(df_lri, df_nonlri)
rm(df_lri)
rm(df_nonlri)

# g = ggplot2::ggplot(data=df_dset, aes(x=num_genes, y=lri_ratio)) +
#   geom_point()
selected_datasets = df_g[gene_type == "NONLRI"][order(-num_genes)]$dataset[1:30]
g = ggplot(data=df_g[dataset %in% selected_datasets], aes(x=reorder(dataset, num_genes), y=num_genes, fill=gene_type)) +
  geom_bar(position="stack", stat='identity') +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 1000, 100)) +
  xlab("Dataset") +
  ylab("Number of DGE genes") +
  ggtitle("DGE: LRI and non-LRI genes by dataset") +
  theme_pubr()
g

g = ggplot(data=df_dset, aes(x=lri_ratio)) +
  geom_histogram(binwidth = 0.05) +
  geom_vline(xintercept=1822/20000, color="red") +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("Fraction of LRI in DGE genes") +
  ylab("Number of datasets") +
  ggtitle("Distribution of LRI fraction in DGE across datasets") +
  theme_pubr()
g
