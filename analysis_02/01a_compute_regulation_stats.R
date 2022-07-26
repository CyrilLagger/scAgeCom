library(Seurat)
library(scDiffCom)
library(data.table)

library(gridExtra)

library(future)
future::plan(future::multisession, workers = 4)

library(glue)
library(ggplot2)

source("utils_random_pairs.R")
source("coexpression.R")

### Parameters ###

NUM_SAMPLES = 50
READ_ONLY = TRUE

##################

### Input data ###

DATASETS = NULL
# DATASETS = list(
#   
#   dataset1 = list(
#     INPUT = "droplet",
#     TISSUE = "Liver",
#     AGE_GROUP = list(
#       YOUNG = c("3m"),
#       OLD = c("18m", "21m", "24m")
#     ),
#     SEX = c("male")
#   )
# 
# dataset2 = list(
#   INPUT = "droplet",
#   TISSUE = "Liver",
#   AGE_GROUP = list(
#     YOUNG = c("3m"),
#     OLD = c("18m", "21m", "24m")
#   ),
#   SEX = c("female")
# ),
# 
# dataset3 = list(
#   INPUT = "droplet",
#   TISSUE = "Kidney",
#   AGE_GROUP = list(
#     YOUNG = c("3m"),
#     OLD = c("18m", "21m", "24m")
#   ),
#   SEX = c("male")
# ),
# 
# dataset4 = list(
#   INPUT = "droplet",
#   TISSUE = "Kidney",
#   AGE_GROUP = list(
#     YOUNG = c("3m"),
#     OLD = c("18m", "21m", "24m")
#   ),
#   SEX = c("female")
# ),
# 
# dataset5 = list(
#   INPUT = "droplet",
#   TISSUE = "Spleen",
#   AGE_GROUP = list(
#     YOUNG = c("3m"),
#     OLD = c("18m", "21m", "24m")
#   ),
#   SEX = c("male")
# ),
# 
# dataset6 = list(
#   INPUT = "droplet",
#   TISSUE = "Spleen",
#   AGE_GROUP = list(
#     YOUNG = c("3m"),
#     OLD = c("18m", "21m", "24m")
#   ),
#   SEX = c("female")
# )
# )

# INPUT = "droplet"  # facs, sample
# TISSUE = "Liver"
# AGE_GROUP = list(
#   YOUNG = c("3m"),
#   OLD = c("18m", "21m", "24m")
# )
# SEX = c("male")

##################

if (is.null(DATASETS)) {
  DATASETS = list()
  
  counter = 0
  for (input in c("droplet", "facs")) {
    seurat_obj = readRDS(glue("../data/seurat_shared_tms_{input}.rds"))
    tissues = unique(seurat_obj[[]]$tissue)
    for (tissue in tissues) {
      for (sex in c("male", "female")) {
        counter = counter + 1
        DATASETS[[counter]] = list(
          INPUT=input,
          TISSUE=tissue,
          AGE_GROUP=list(
              YOUNG = c("3m"),
              OLD = c("18m", "21m", "24m")
            ),
          SEX=sex
        )
      }
    }
  }
}
rm(seurat_obj)

for (dataset in DATASETS) {
  print(dataset)
  
  INPUT = dataset$INPUT
  TISSUE = dataset$TISSUE
  AGE_GROUP = dataset$AGE_GROUP
  SEX = dataset$SEX

  SAVE_PATH = paste0("data/shuffled_regulation/shuffled_regulation_", INPUT, "_", TISSUE, "_", SEX, ".csv")
  SVG_SAVE = paste0("data/shuffled_regulation/shuffled_regulation_", INPUT, "_", TISSUE, "_", SEX, ".svg")
  PNG_SAVE = paste0("data/shuffled_regulation/shuffled_regulation_", INPUT, "_", TISSUE, "_", SEX, ".png")
  STATS_SAVE = paste0("data/shuffled_regulation/shuffled_regulation_pvals_", INPUT, "_", TISSUE, "_", SEX, ".csv")

  tryCatch(
    {
      seurat_sample_obj = get_seurat(INPUT, TISSUE, AGE_GROUP, SEX)
      seurat_g = get_genes_from_seurat(seurat_sample_obj)
      lri = scDiffCom::LRI_mouse$LRI_curated
      lri_simple = subset_simple_lri(lri)
      lri_simple_in_seurat = lri_simple[ (LIGAND_1 %in% seurat_g) & (RECEPTOR_1 %in% seurat_g) ]
      
      
      res_scDiffCom = run_internal_analysis(
        seurat_obj = seurat_sample_obj,
        lri_table = lri_simple_in_seurat
      )
      res_scDiffCom = FilterCCI(res_scDiffCom, skip_ora = TRUE)
      
      
      t = simulate_random_seurat_regulation(
        seurat_sample_obj, 
        lri_simple_in_seurat, 
        SAVE_PATH,
        NUM_SAMPLES,
        read_only=READ_ONLY
      )
      t[is.na(t)] = 0
      
      detected_reshaped = reshape2::melt(table(res_scDiffCom@cci_table_detected$REGULATION))
      rownames(detected_reshaped) = detected_reshaped$Var1
      
      tests = list()
      fc_average = c()
      observed_vals = c()
      mean_random_vals = c()
      for (regulation_type in c("DOWN", "UP", "NSC", "FLAT")) {
        random_vals = t[,regulation_type, with=FALSE][[1]]
        observed_val = detected_reshaped[regulation_type, "value"]
        tests[[regulation_type]] = t.test(random_vals, mu = observed_val, alternative = "two.sided")
        fc_average[regulation_type] = observed_val / mean(random_vals)
        observed_vals[regulation_type] = observed_val
        mean_random_vals[regulation_type] = mean(random_vals)
      }
      test_pvals = sapply(tests, function(test) {test$p.value})
      
      random_totals = t$UP + t$DOWN + t$NSC + t$FLAT
      observed_total = length(res_scDiffCom@cci_table_detected$REGULATION)
      total = c(
        p_val=t.test(random_totals, mu = observed_total, alternative = "two.sided")$p.value,
        fc_obs_vs_random=observed_total / mean(random_totals),
        observed_val=observed_total,
        mean_random_vals=mean(random_totals)
      )
      
      stats = cbind(test_pvals, fc_average, observed_vals, mean_random_vals)
      colnames(stats) = c("p_val", "fc_obs_vs_random", "observed_val", "mean_random_vals")
      stats = rbind(stats, total)
      rownames(stats) = c("DOWN", "UP", "NSC", "FLAT", "TOTAL")
      write.csv(stats, STATS_SAVE)
      
      g = ggplot(reshape2::melt(t), aes(x=variable, y=value)) + 
        geom_boxplot() +
        geom_point(
          data=detected_reshaped,
          aes(x=Var1, y=value),
          color = "red"
        ) +
        xlab("REGULATION type") +
        ylab("COUNTS") +
        ggtitle(glue("{INPUT}, {TISSUE}, {SEX}\nBoxplot of `REGULATION` for random LRI"))

      ggsave(SVG_SAVE, 
             plot = g)
      ggsave(PNG_SAVE, 
             plot = g)  
    },

    error=function(cond) {
      message("An Error occured")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )
}
