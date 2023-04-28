library(glue)
library(data.table)
library(scDiffCom)

dt_cci_age = readRDS("data/data_02-08-2022/dt_cci_age.rds")
dt_cci_sex = readRDS("data/data_02-08-2022/dt_cci_sex.rds")
# dt_secr_validation_final = readRDS("data/data_02-08-2022/dt_secr_validation_final.rds")
# scAgeComShiny_data_02_08_2022 = readRDS("data/data_02-08-2022/scAgeComShiny_data_02_08_2022.rds")

dt_cci_age_2 <- copy(dt_cci_age)
dt_cci_sex_2 <- copy(dt_cci_sex)

## What dataset:tissue pairs are age and sex analyzable?
dt_cci_age_2$dataset_src = sapply(strsplit(dt_cci_age_2$dataset, split=' \\('), function(elem) {elem[1]})
dt_cci_age_2$dataset_sample_type = ifelse(
  sapply(strsplit(dt_cci_age_2$dataset, split=' \\('), function(elem) {elem[2]}) == "male)",
  "male",
  "female"
)
dt_cci_age_2$dataset_src_tissue = paste0(dt_cci_age_2$dataset_src, "_", dt_cci_age_2$tissue)

dt_cci_sex_2$dataset_src = sapply(strsplit(dt_cci_sex_2$dataset, split=' \\('), function(elem) {elem[1]})
dt_cci_sex_2$dataset_sample_type = ifelse(
  sapply(strsplit(dt_cci_sex_2$dataset, split=' \\('), function(elem) {elem[2]}) == "old)",
  "old",
  "young"
)
dt_cci_sex_2$dataset_src_tissue = paste0(dt_cci_sex_2$dataset_src, "_", dt_cci_sex_2$tissue)

ds_tissue_pairs = union(unique(dt_cci_age_2$dataset_src_tissue), unique(dt_cci_sex_2$dataset_src_tissue))

dt_ds_tissue_counts = data.table()
for (ds_tissue in ds_tissue_pairs) {
  ds_src = strsplit(ds_tissue, "_")[[1]][1]
  ts = strsplit(ds_tissue, "_")[[1]][2]

  ds_sex_old = glue("{ds_src} (old)")
  ds_sex_young = glue("{ds_src} (young)")
  ds_age_male = glue("{ds_src} (male)")
  ds_age_female = glue("{ds_src} (female)")
    
  N_sex_old = dt_cci_sex_2[dataset == ds_sex_old & tissue == ts, .N]
  N_sex_young = dt_cci_sex_2[dataset == ds_sex_young & tissue == ts, .N]
  N_age_male = dt_cci_age_2[dataset == ds_age_male & tissue == ts, .N]
  N_age_female = dt_cci_age_2[dataset == ds_age_female & tissue == ts, .N]
  
  v = data.table(
    dataset = ds_src, 
    tissue = ts, 
    N_sex_old = N_sex_old, 
    N_sex_young = N_sex_young, 
    N_age_male = N_age_male, 
    N_age_female = N_age_female
  )
  dt_ds_tissue_counts = rbind(dt_ds_tissue_counts, v)
}
dt_ds_tissue_counts = dt_ds_tissue_counts[N_sex_old > 0 & N_sex_young > 0 & N_age_male > 0 & N_age_female]
dt_ds_tissue_counts[, N_total := N_sex_old + N_sex_young + N_age_male + N_age_female]
dt_ds_tissue_counts = dt_ds_tissue_counts[order(-N_total)]


## Finding the direction in sex analysis: F -> M
# dt_cci_sex_2[, list(
#   L1_EXPRESSION_female, L2_EXPRESSION_female, R1_EXPRESSION_female, R2_EXPRESSION_female, R3_EXPRESSION_female, CCI_SCORE_female,
#   L1_EXPRESSION_male, L2_EXPRESSION_male, R1_EXPRESSION_male, R2_EXPRESSION_male, R3_EXPRESSION_male, CCI_SCORE_male,
#   LOGFC)]


## Evaluate the relation: Sy + Am = Af + So

dt_age_male = dt_cci_age_2[dataset_sample_type == "male"]
dt_age_female = dt_cci_age_2[dataset_sample_type == "female"]
dt_sex_young = dt_cci_sex_2[dataset_sample_type == "young"]
dt_sex_old = dt_cci_sex_2[dataset_sample_type == "old"]

dt_age = merge(
  x = dt_age_male,
  y = dt_age_female,
  by = c("dataset_src", "tissue", "CCI"),
  all = FALSE,
  suffixes = c("_Am", "_Af")
)

dt_sex = merge(
  x = dt_sex_young,
  y = dt_sex_old,
  by = c("dataset_src", "tissue", "CCI"),
  all = FALSE,
  suffixes = c("_Sy", "_So")
)

dt = merge(
  x = dt_age,
  y = dt_sex,
  by = c("dataset_src", "tissue", "CCI"),
  all = FALSE,
  suffixes = c("", "")
)

## Histogram the differences: zero.
Sy = dt$LOGFC_Sy
So = dt$LOGFC_So
Am = dt$LOGFC_Am
Af = dt$LOGFC_Af
hist(Sy + Am - (Af + So))

dim(dt)
## 72.3k CCI (from all dataset:tissue)

delta = Sy + Am  # or = Af + So
delta_clean = delta[!is.na(delta) & abs(delta) != Inf]
sd(delta_clean)
# SD = 0.5 --> delta is 0 +/- 1

## Ignore REGULATION_So, as 3 of the 4 are enough to characterize the results
dt_regulation_counts = dt[, list(REGULATION_Sy, REGULATION_Am, REGULATION_Af)]
dt_regulation_counts[dt_regulation_counts == "NSC"] = "NSC/FLAT"
dt_regulation_counts[dt_regulation_counts == "FLAT"] = "NSC/FLAT"
dt_regulation_counts = dt_regulation_counts[, .N, by=list(REGULATION_Sy, REGULATION_Am, REGULATION_Af)]
dt_regulation_counts = dt_regulation_counts[order(-N)]
dt_regulation_counts_updown = dt_regulation_counts[
  REGULATION_Am %in% c("UP", "DOWN") 
  | REGULATION_Af  %in% c("UP", "DOWN")
  | REGULATION_Sy %in% c("UP", "DOWN")
  # | REGULATION_So %in% c("UP", "DOWN")
]
dt_regulation_counts_updown$N_perc = dt_regulation_counts_updown$N / sum(dt_regulation_counts_updown$N)
dt_regulation_counts_updown$N_perc_cumsum = cumsum(dt_regulation_counts_updown$N_perc)
dt_regulation_counts_updown$SIGNATURE = paste0(
  substr(dt_regulation_counts_updown$REGULATION_Sy, 1, 1),
  substr(dt_regulation_counts_updown$REGULATION_Am, 1, 1),
  substr(dt_regulation_counts_updown$REGULATION_Af, 1, 1)
)
plot(dt_regulation_counts_updown$N_perc_cumsum)

dt_updown = dt[
  REGULATION_Am %in% c("UP", "DOWN") 
  | REGULATION_Af  %in% c("UP", "DOWN")
  | REGULATION_Sy %in% c("UP", "DOWN")
  # | REGULATION_So %in% c("UP", "DOWN")
]

dt_updown[dt_updown == "FLAT"] = "NSC/FLAT"
dt_updown[dt_updown == "NSC"] = "NSC/FLAT"
dt_updown$SIGNATURE = paste0(
  substr(dt_updown$REGULATION_Sy, 1, 1),
  substr(dt_updown$REGULATION_Am, 1, 1),
  substr(dt_updown$REGULATION_Af, 1, 1)
)

fwrite(
  dt_regulation_counts_updown,
  paste0(
    path_scagecom_output,
    "dt_regulation_counts_updown.csv"
  )
)

## Patterns by Dataset_src : Tissue
dt_tab = data.table()
for (dataset_src_tissue_i in unique(dt$dataset_src_tissue_Am)) {
  tab = table(dt_updown[dataset_src_tissue_Af == dataset_src_tissue_i]$SIGNATURE)
  tab = tab / sum(tab)
  dt_tab_i = data.table(tab)
  colnames(dt_tab_i) = c("SIGNATURE", dataset_src_tissue_i)
  
  if (nrow(dt_tab) == 0) {
    dt_tab = dt_tab_i
    next()
  } else {
    dt_tab = merge(
      dt_tab,
      dt_tab_i,
      by = "SIGNATURE",
      all.x = TRUE,
      all.y = TRUE
    )
  }
}

dt_tab = merge(dt_regulation_counts_updown, dt_tab, by = "SIGNATURE", all=TRUE)
dt_tab = dt_tab[order(N_perc_cumsum)]


## Some particular analyses
### Cross-tissue
tab = table(dt_updown[SIGNATURE == "UDN"]$LIGAND_1_Sy)
tab = tab[order(tab)]
tab = tab / sum(tab)

### Adipose
dt_updown[dataset_src_tissue_Sy == "TMS FACS_Adipose_Subcutaneous" & SIGNATURE == "NNU"]

### Lung
dt_lung = dt_updown[dataset_src_tissue_Sy == "TMS FACS_Lung"]
dt_lung_udu = dt_updown[dataset_src_tissue_Sy == "TMS FACS_Lung" & SIGNATURE == "UDU"]
tab = table(dt_lung_udu$LIGAND_1_Sy)
tab = tab[order(tab)]
tab = tab / sum(tab)

tab = table(dt_lung_udu$RECEPTOR_1_Sy)
tab = tab[order(tab)]
tab = tab / sum(tab)

tab = table(dt_lung_udu$ER_CELLTYPES_Am)
tab = tab[order(tab)]
tab = tab / sum(tab)

perc_bronchial = sum(
  dt_lung_udu$EMITTER_CELLTYPE_Sy == "bronchial smooth muscle cell" 
  | dt_lung_udu$RECEIVER_CELLTYPE_Sy == "bronchial smooth muscle cell"
) / dim(dt_lung_udu)[1]

# perc_bronchial_lung = sum(
#   dt_lung$EMITTER_CELLTYPE_Sy == "bronchial smooth muscle cell" 
#   | dt_lung$RECEIVER_CELLTYPE_Sy == "bronchial smooth muscle cell"
# ) / dim(dt_lung)[1]

## App:Lrp10
dt[LIGAND_1_Af == "App" & RECEPTOR_1_Af == "Lrp10"]$SIGNATURE
tab = table(dt_updown[LIGAND_1_Af == "App" & RECEPTOR_1_Af == "Lrp10"]$SIGNATURE)
tab = tab[order(tab)]
tab = tab / sum(tab)
