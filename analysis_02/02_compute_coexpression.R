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
library(ggpubr)

source("utils_random_pairs.R")
source("coexpression.R")


lri = scDiffCom::LRI_mouse$LRI_curated
lri = subset_simple_lri(lri)
 

## Mouse LRI to Ensembl IDs
gene_ids = map_ensembl_ids(
  get_genes_lri(lri)
)$ensembl_gene_id


## Main

DF_COEXP_FP = "data/coexpression_plots/df_coexp.csv"

if (!file.exists(DF_COEXP_FP)) {
  coexp_adj = read_network()
  if (sum(is.na(coexp_adj)) > 0) {
    stop("NAs in coexpression matrix.")
  }
  
  lri = add_lri_ensembl_ids(lri)
  lri = add_lri_coexpression_adj(lri, coexp_adj)
  
  # Drop lri NA coexpression (just not found)
  lri_coexpression = lri$coexpression
  lri_coexpression = lri_coexpression[!is.na(lri_coexpression)]
  
  
  # Fetch non LRI co-expressions
  COEXP_NO_LRI_PATH = "data/coexp_adj_no_lri.Rdata"
  
  if (!file.exists(COEXP_NO_LRI_PATH)) {
    coexp_adj_no_lri = copy(coexp_adj)
    ids_in_coexp = colnames(coexp_adj_no_lri)
    L1_coexp_idx = purrr:::map_int(lri$LIGAND_1_ID, function(id) {
      match(id, ids_in_coexp)
    })
    R1_coexp_idx = purrr:::map_int(lri$RECEPTOR_1_ID, function(id) {
      match(id, ids_in_coexp)
    })
    lri$LIGAND_1_COEXP_IDX = L1_coexp_idx
    lri$RECEPTOR_1_COEXP_IDX = R1_coexp_idx
    
    sum((!is.na(L1_coexp_idx) & !is.na(R1_coexp_idx)))
    
    # SLOW
    counter = 0
    r = map2(lri$LIGAND_1_COEXP_IDX, lri$RECEPTOR_1_COEXP_IDX, function(L_idx, R_idx) {
      if (
        ( (!is.na(L_idx)) & (!is.na(R_idx)) )
      ) {
        coexp_adj_no_lri[L_idx, R_idx] <<- NA  # global assignment
        coexp_adj_no_lri[R_idx, L_idx] <<- NA
        counter = counter + 1
        if (counter %% 50 == 0) {
          print(glue("COUNTER = {counter}"))
        }
      }
    })
    save(coexp_adj_no_lri, file=COEXP_NO_LRI_PATH)  
  } else {
    load(COEXP_NO_LRI_PATH)
  }
  
  uptri = upper.tri(coexp_adj_no_lri, diag = FALSE)
  non_lri_coexp = coexp_adj_no_lri[uptri]
  non_lri_coexp = non_lri_coexp[!is.na(non_lri_coexp)]
  rm(coexp_adj_no_lri)
  rm(uptri)
  
  df_coexp = rbind(
    data.frame(weights = lri_coexpression, type = "LRI"),
    data.frame(weights = non_lri_coexp, type = "NON-LRI")
  )
  write.csv(df_coexp, DF_COEXP_FP)
  rm(lri_coexpression)
  rm(non_lri_coexp)
  
} else {
  df_coexp = read.csv(DF_COEXP_FP)
}


g_hist_lri = ggplot(data=df_coexp[df_coexp$type=="LRI",], aes(x=weights)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  geom_vline(xintercept=median(df_coexp[df_coexp$type=="LRI",]$weights), color="red") +
  xlab("Normalized coexpression rank") +
  ylab("Count") +
  ggtitle("Distribution of coexpression ranks in gene pairs from the LRI dataset") +
  theme_pubr()

g_hist_nonlri = ggplot(data=df_coexp[df_coexp$type=="NON-LRI",], aes(x=weights)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', boundary=0) +
  xlab("Normalized coexpression rank") +
  ylab("Count") +
  ggtitle("Distribution of coexpression ranks in random gene pairs") +
  theme_pubr()

ggsave(glue("data/coexpression_plots/coexp_lri_hist.png"),
       plot = g_hist_lri)
ggsave(glue("data/coexpression_plots/coexp_lri_hist.svg"),
       plot = g_hist_lri)
ggsave(glue("data/coexpression_plots/coexp_nonlri_hist.png"),
       plot = g_hist_nonlri)
ggsave(glue("data/coexpression_plots/coexp_nonlri_hist.svg"),
       plot = g_hist_nonlri)
