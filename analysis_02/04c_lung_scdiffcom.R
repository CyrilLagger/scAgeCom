library(data.table)

load("data/proteins/GSE124872_raw_counts_single_cell.RData")
dim(raw_counts)

meta = data.table(read.csv("data/proteins/GSE124872_Angelidis_2018_metadata.csv"))
meta$age_group = ifelse(meta$grouping == "3m", "young", "old")
selected_celltypes = meta[
  , 
  list(
    N_y = sum(age_group=="young"), 
    N_o = sum(age_group=="old")), 
  by = celltype][
    N_y >= 5 & N_o >= 5
    ]$celltype
meta = meta[celltype %in% selected_celltypes]
