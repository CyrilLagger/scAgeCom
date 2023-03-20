library(Seurat)
library(scDiffCom)
library(data.table)

library(gridExtra)

library(purrr)
library(furrr)
library(future)
future::plan(future::multisession, workers = 8)

library(glue)
library(ggplot2)

source("utils_random_pairs.R")
source("coexpression.R")


lri = scDiffCom::LRI_mouse$LRI_curated
lri = subset_simple_lri(lri)


## Mouse LRI to Ensembl IDs
gene_ids = map_ensembl_ids(
  get_genes_lri(lri)
)$ensembl_gene_id


coexp_adj = read_network()
if (sum(is.na(coexp_adj)) > 0) {
  stop("NAs in coexpression matrix.")
}

lri = add_lri_ensembl_ids(lri)
lri = add_lri_coexpression_adj(lri, coexp_adj)
lri = lri[!is.na(coexpression)]
rm(coexp_adj)

# Drop lri NA coexpression (just not found)
# lri_coexpression = lri$coexpression
# lri_coexpression = lri_coexpression[!is.na(lri_coexpression)]

DATASETS = NULL
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
  
  SAVE_PATH = paste0("data/datasets_detections/detected_", INPUT, "_", TISSUE, "_", SEX, ".csv")
  
  tryCatch(
    {
      s = get_seurat(
        INPUT = INPUT,
        TISSUE = TISSUE,
        AGE_GROUP = AGE_GROUP,
        SEX = SEX
      )
      res_scDiffCom = run_internal_analysis(
        seurat_obj = s,
        lri_table = lri
      )
      res_scDiffCom = FilterCCI(res_scDiffCom, skip_ora = TRUE)
      det = res_scDiffCom@cci_table_detected
      det = merge(det, lri[, list(LRI, coexpression)], by="LRI", all.x = TRUE)
      write.csv(det, SAVE_PATH)
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

