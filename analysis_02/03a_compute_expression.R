library(Seurat)
library(scDiffCom)
library(data.table)

library(gridExtra)

library(future)
future::plan(future::multisession, workers = 4)

library(glue)
library(ggplot2)

source("utils_random_pairs.R")

### Parameters ###

MARKERS_EXHAUSTIVE_FP = "data/expression/markers_full.csv"
MARKERS_AGG_FP = "data/expression/markers_agg.csv"

##################

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
rm(seurat_obj)


# Main

lri = scDiffCom::LRI_mouse$LRI_curated
lri = subset_simple_lri(lri)
lri_genes = get_genes_lri(lri)

if (file.exists("data/expression/markers_full.csv")) {
  df = data.table::data.table(read.csv("data/expression/markers_full.csv"))
  df$dataset = paste0(df$input, "_", df$tissue, "_", df$sex)
  computed_datasets = unique(df$dataset)
  rm(df)
} else {
  computed_datasets = c()
}

markers_exhaustive = data.table()
markers_agg = data.table()
counter = 0
for (dataset in DATASETS) {
  print(dataset)
  counter = counter + 1
  
  INPUT = dataset$INPUT
  TISSUE = dataset$TISSUE
  AGE_GROUP = dataset$AGE_GROUP
  SEX = dataset$SEX
  
  dataset = paste0(INPUT, "_", TISSUE, "_", SEX)
  if (dataset %in% computed_datasets) {
    print(paste0("Dataset ", dataset, " already processed."))
    next()
  }
  
  tryCatch(
    {
      seurat_sample_obj = get_seurat(INPUT, TISSUE, AGE_GROUP, SEX)
      markers = dge_by_celltypes(seurat_sample_obj)
      
      if (nrow(markers) == 0) {
        print(glue("No markers for dataset # {counter}"))
        next()
      }
      
      markers[, "input"] = INPUT
      markers[, "tissue"] = TISSUE
      markers[, "sex"] = SEX
      
      markers_exhaustive = rbind(markers_exhaustive, markers)
      write.csv(markers_exhaustive, MARKERS_EXHAUSTIVE_FP)
      
      
      num_dge = nrow(markers)
      num_dge_lri = sum(markers$gene %in% lri_genes)
      markers_agg_dataset = data.table(
        input=INPUT,
        tissue=TISSUE,
        sex=SEX,
        num_dge=num_dge,
        num_dge_lri=num_dge_lri
      )
      markers_agg = rbind(markers_agg, markers_agg_dataset)
      write.csv(markers_agg, MARKERS_AGG_FP)
    },
    
    error=function(cond) {
      message("An Error occured")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )
  
  if (counter == 2) {
    stop()
  }

}
