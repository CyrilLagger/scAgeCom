library(ggplot2)

# Renaming issue
# files = list.files(path="./data/shuffled_regulation", pattern="*pvals*")
# files = paste0("data/shuffled_regulation/", files)
# renamed = paste0(sapply(strsplit(files, "csv"), function (x) {x}), ".csv")
# file.rename(files, renamed)

files = list.files(path="./data/shuffled_regulation", pattern="*pvals*")
files = paste0("data/shuffled_regulation/", files)

dt_stat = data.table::data.table()
categories = c("UP", "DOWN", "FLAT", "NSC", "TOTAL")
for (category in categories) {
  for (file in files) {
    dt = data.table::data.table(read.csv(file))
    dt$category = dt$X
    dt$X = NULL
    dt$file = file
    # dataset from file
    dataset = strsplit(file, ".csv")[[1]]
    dataset = strsplit(dataset, "pvals_")[[1]][2]
    dt$dataset = dataset
    
    mask = dt$category == category
    dt = dt[mask]
    dt_stat = rbind(dt_stat, dt)
  }
}
dt_stat[, "-log_pval" := -log(p_val, base=10)]
dt_stat[, "logFC" := log(fc_obs_vs_random, base=2)]

dt_sim = data.table::data.table()
for (dataset in unique(dt_stat$dataset)) {
  simulation_res_filepath = glue::glue("data/shuffled_regulation/shuffled_regulation_{dataset}.csv")
  dt = data.table::data.table(read.csv(simulation_res_filepath))
  dt$dataset = dataset
  dt_sim = rbind(dt_sim, dt)
}
dt_sim[is.na(dt_sim)] = 0
dt_sim[, TOTAL := DOWN+FLAT+NSC+UP]
dt_sim = dt_sim[!(UP==1 & DOWN==1 & FLAT==1 & NSC==1)]
dt_sim = dt_sim[order(TOTAL)]


## Plots

# Total detected - comparative boxplot
selected_datasets = dt_stat[category == "TOTAL"][order(-mean_random_vals)]$dataset[1:30]
g = ggplot(data=dt_sim[dataset %in% selected_datasets], aes(y=reorder(dataset, TOTAL), x=TOTAL)) +
  geom_boxplot() +
  geom_point(
    data=dt_stat[dataset %in% selected_datasets & category == "TOTAL"],
    aes(x=observed_val, y=dataset),
    color = "red"
  ) +
  scale_x_continuous(breaks = seq(0, 30000, 5000)) +
  xlab("Total detected CCI") +
  ylab("TMS dataset") +
  ggtitle("Detection events in TMS datasets with LRI and random genes") +
  theme_pubr()
g

# Total detected - histogram (not useful) and table
# ggplot(data=dt_stat[category == "TOTAL"], aes(x=p_val)) +
#   geom_histogram() +
#   scale_x_continuous(breaks = seq(0, 1, 0.05))
dt_stat[, is_pval_sig := p_val < 0.05]
dt_stat[, direction := ifelse(logFC > 0, "INCREASE", "DECREASE")]
table(dt_stat[category == "TOTAL", list(is_pval_sig, direction)])

table(dt_stat[category == "UP", list(is_pval_sig, direction)])

# logFC histogram per group
# g = ggplot(data=dt_stat[category != "TOTAL"], aes(y=logFC, x=category)) +
#   geom_boxplot() +
#   theme_pubr()
# g

p <- ggboxplot(dt_stat[category != "TOTAL"], x = "category", y = "logFC",
               add = "jitter")
my_comparisons <- list( c("DOWN", "NSC"), c("DOWN", "FLAT"), c("UP", "FLAT"), c("UP", "NSC") )
p = p + stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 5) +
  xlab("Category of regulation") +
  ylab("log2FC") +
  scale_y_continuous(breaks = seq(-1, 5, 0.5)) +
  ggtitle("Increased detection events across TMS datasets by regulation category by using LRI or random genes")
p

# Bingo
t.test(dt_stat[category == "UP"]$fc_obs_vs_random, dt_stat[category == "NSC"]$fc_obs_vs_random)
t.test(dt_stat[category == "UP"]$fc_obs_vs_random, dt_stat[category == "FLAT"]$fc_obs_vs_random)
t.test(dt_stat[category == "DOWN"]$fc_obs_vs_random, dt_stat[category == "NSC"]$fc_obs_vs_random)
t.test(dt_stat[category == "DOWN"]$fc_obs_vs_random, dt_stat[category == "FLAT"]$fc_obs_vs_random)
