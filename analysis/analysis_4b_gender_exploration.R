

## Load the processed datasets ####

DATASETS_PROCESSED <- c(
  readRDS("../data_scAgeCom/analysis/outputs_data/TMS_scAgeCom_processed_bestORA.rds"),
  readRDS("../data_scAgeCom/analysis/outputs_data/TMS_scAgeCom_processed_bestORA_w30m.rds")
)


## Load and process general md ####
seurat_md <- readRDS("../data_scAgeCom/analysis/outputs_data/seurat_md.rds")
seurat_md <- lapply(
  seurat_md,
  setDT
)
seurat_md$calico_kidney[, tissue :=  "Kidney"]
seurat_md$calico_lung[, tissue := "Lung"]
seurat_md$calico_spleen[, tissue := "Spleen"]
seurat_md$calico <- rbindlist(
  l = list(seurat_md$calico_kidney, seurat_md$calico_lung, seurat_md$calico_spleen)
)
seurat_md$calico_kidney <- NULL
seurat_md$calico_lung <- NULL
seurat_md$calico_spleen <- NULL
seurat_md$calico[, sex := "male"]

lapply(
  seurat_md, 
  function(i) {
    i[, age_group := ifelse(age %in% c('1m', '3m', "young"), 'YOUNG', 'OLD')]
  }
)

seurat_md <- lapply(seurat_md, function(x) {
  x$tissue_cell_type <- paste(x$tissue, x$cell_ontology_final)
  return(x)
})

seurat_md$tms_facs[, .N, by = c("sex", "mouse.id", "tissue") ][, .N, by = "tissue"]
seurat_md$tms_facs[, .N, by = c("mouse.id", "tissue") ][, .N, by = "tissue"]
counts_mice_ct <- seurat_md$tms_facs[, .N, by = c("mouse.id", "tissue_cell_type") ][, .N, by = "tissue_cell_type"]

facs_summary <- merge.data.table(
  seurat_md$tms_facs[, .N, by = c("tissue", "tissue_cell_type", "age_group", "sex")],
  seurat_md$tms_facs[, .N, by = c("tissue", "mouse.id", "tissue_cell_type", "age_group", "sex")][
    , .N, by = c("tissue", "tissue_cell_type", "age_group", "sex")],
  by =  c("tissue", "tissue_cell_type", "age_group", "sex"),
  suffixes = c("_cells", "_mice")
)

facs_ct_by_genderAge <- dcast.data.table(
  facs_summary,
  tissue_cell_type ~ sex + age_group,
  value.var = "N_mice"
)

droplet_summary <- merge.data.table(
  seurat_md$tms_droplet[, .N, by = c("tissue", "tissue_cell_type", "age_group", "sex")],
  seurat_md$tms_droplet[, .N, by = c("tissue", "mouse.id", "tissue_cell_type", "age_group", "sex")][
    , .N, by = c("tissue", "tissue_cell_type", "age_group", "sex")],
  by =  c("tissue", "tissue_cell_type", "age_group", "sex"),
  suffixes = c("_cells", "_mice")
)

droplet_ct_by_genderAge <- dcast.data.table(
  droplet_summary,
  tissue_cell_type ~ sex + age_group,
  value.var = "N_mice"
)

ggplot(facs_summary, aes(N_mice, log10(N_cells), color = sex)) + geom_point()

ggplot(facs_summary, aes(sex, N_mice, color = age_group)) + geom_boxplot()
ggplot(droplet_summary, aes(sex, N_mice, color = age_group)) + geom_boxplot()

ggplot(facs_summary[, .N, by = c("tissue", "sex", "age_group")], 
       aes(age_group, N, color = sex)) +
  geom_boxplot() +
  facet_wrap(~tissue)


seurat_md$tms_facs$n_counts

ggplot(seurat_md$tms_facs, aes(age_group, log10(n_counts), color = sex)) + geom_boxplot() + 
  facet_wrap(~tissue)

ggplot(seurat_md$tms_droplet, aes(age_group, log10(n_counts), color = sex)) + geom_boxplot() + 
  facet_wrap(~tissue)

ggplot(seurat_md$calico, aes(age_group, log10(n_counts), color = sex)) + geom_boxplot() + 
  facet_wrap(~tissue)

table(seurat_md$tms_facs$mouse.id)
table(seurat_md$tms_droplet$mouse.id)

sort(unique(seurat_md$tms_facs$mouse.id))
table(seurat_md$tms_droplet$mouse.id)

## Table of number of cells by age and by gender ####

cells_by_AgeGender <- rbindlist(
  lapply(
    seurat_md,
    function(dataset) {
      dt <- dcast.data.table(
        dataset[age != "30m", .N, by = c("age_group", "sex", "tissue_cell_type")],
        formula = tissue_cell_type ~ age_group + sex,
        value.var = "N"
      )
      dt[is.na(dt)] <- 0
      #dt[, OR := OLD_male*YOUNG_female/(YOUNG_male*OLD_female)]
      #dt[, JACCARD := (OLD_male + OLD_female)/((OLD_male + YOUNG_female)+(YOUNG_male + OLD_female)+(OLD_male + OLD_female))]
      dt
    }
  ),
  idcol = "dataset",
  fill = TRUE
)
cells_by_AgeGender[is.na(cells_by_AgeGender)] <- 0

## Common CCI between mixed and gender ####

names(DATASETS_PROCESSED)

droplet_all_cci <- rbindlist(
  lapply(
    DATASETS_PROCESSED[c(2, 10, 11)],
    function(i) {
      i$dataset@cci_detected
    }
  ),
  idcol = "dataset"
)
droplet_all_cci[, TCCI := paste(ID, CCI, sep = "_")]

facs_tcci_distr <- dcast.data.table(
  facs_all_cci[, c("dataset", "REGULATION", "TCCI")],
  formula = TCCI ~ dataset + REGULATION,
  fun.aggregate = length
)
facs_tcci_counts <- facs_tcci_distr[, .N, by = c(colnames(facs_tcci_distr)[-1])]
facs_updown_female_male <- facs_tcci_distr[facs_female_UP == 1 & facs_male_DOWN == 1, "TCCI"]
facs_updown_male_female <- facs_tcci_distr[facs_male_UP == 1 & facs_female_DOWN == 1, "TCCI"]

facs_updown_female_male[facs_all_cci, on = "TCCI", LR_GENES := i.LR_GENES]
facs_updown_female_male[facs_all_cci, on = "TCCI", LIGAND_1 := i.LIGAND_1]
facs_updown_female_male[facs_all_cci, on = "TCCI", LIGAND_2 := i.LIGAND_2]
facs_updown_female_male[facs_all_cci, on = "TCCI", RECEPTOR_1 := i.RECEPTOR_1]
facs_updown_female_male[facs_all_cci, on = "TCCI", RECEPTOR_2 := i.RECEPTOR_2]
test4 <- facs_updown_female_male[, .N, by = LIGAND_1][order(-N)]

facs_all_cci[LIGAND_1 == "Itgb1"]$L1_EXPRESSION_OLD

droplet_tcci_distr <- dcast.data.table(
  droplet_all_cci[, c("dataset", "REGULATION", "TCCI")],
  formula = TCCI ~ dataset + REGULATION,
  fun.aggregate = length
)
droplet_tcci_counts <- droplet_tcci_distr[, .N, by = c(colnames(droplet_tcci_distr)[-1])]
droplet_updown_female_male <- droplet_tcci_distr[droplet_female_UP == 1 & droplet_male_w30_DOWN == 1, "TCCI"]
droplet_updown_male_female <- droplet_tcci_distr[droplet_male_w30_UP == 1 & droplet_female_DOWN == 1, "TCCI"]


ggplot(facs_all_cci[RECEPTOR_1 == "Itgb1"],
       aes(R1_EXPRESSION_YOUNG, R1_EXPRESSION_OLD, color = dataset)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)

ggplot(facs_all_cci[RECEPTOR_1 == "Itgb1"],
       aes(CCI_SCORE_YOUNG, CCI_SCORE_OLD, color = dataset)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)

ggplot(facs_all_cci[LIGAND_1 == "App"],
       aes(L1_EXPRESSION_YOUNG, L1_EXPRESSION_OLD, color = dataset)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)

ggplot(facs_all_cci[LIGAND_1 == "App"],
       aes(CCI_SCORE_YOUNG, CCI_SCORE_OLD, color = dataset)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)

ggplot(facs_all_cci[LIGAND_1 == "App"],
       aes(LOGFC, color = dataset, fill = dataset)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.4)

ggplot(facs_all_cci[RECEPTOR_1 == "Itgb1"],
       aes(LOGFC, color = dataset, fill = dataset)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.4)

ggplot(facs_all_cci,
       aes(LOGFC, color = dataset, fill = dataset)) +
  geom_histogram(bins = 100, aes(y = ..density..), position = "identity", alpha = 0.4) +
  geom_density(alpha = 0.4)

ggplot(facs_all_cci,
       aes(LOGFC, color = dataset, fill = dataset)) +
  geom_histogram(bins = 100, aes(y = ..density..), position = "identity", alpha = 0.4) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ID) +
  xlim(c(-3, 3))

ggplot(facs_all_cci,
       aes(LOGFC, color = dataset, fill = dataset)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.4) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ID) +
  xlim(c(-3, 3))

ggplot(droplet_all_cci,
       aes(LOGFC, color = dataset, fill = dataset)) +
  geom_histogram(bins = 100, aes(y = ..density..), position = "identity", alpha = 0.4) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ID)

ggplot(facs_all_cci,
       aes(LOGFC, color = ID, fill = ID)) +
  geom_histogram(bins = 100, aes(y = ..density..), position = "identity", alpha = 0.4) +
  geom_density(alpha = 0.4)

ggplot(droplet_all_cci,
       aes(LOGFC, color = dataset, fill = dataset)) +
  geom_histogram(bins = 100, aes(y = ..density..), position = "identity", alpha = 0.4) +
  geom_density(alpha = 0.4)

summary(facs_all_cci[dataset == "facs_male"]$LOGFC)
summary(facs_all_cci[dataset == "facs_female"]$LOGFC)

droplet_pct_comp <- droplet_all_cci[, .N, by = c("dataset", "ID", "REGULATION")]
droplet_pct_comp[, PCT := N/sum(N)*100, by = c("dataset", "ID") ]
droplet_pct_comp <- dcast.data.table(
  droplet_pct_comp,
  formula = dataset + ID ~ REGULATION,
  value.var = "PCT"
)


TCCI_by_REGULATION <- dcast.data.table(
  droplet_all_cci[ID == "Heart_and_Aorta", .N, by = c("REGULATION", "TCCI")],
  formula = TCCI ~ REGULATION,
  value.var = "N"
)
TCCI_by_REGULATION[is.na(TCCI_by_REGULATION)] <- 0
TCCI_by_REGULATION[, TOTAL := UP + DOWN + NON_SIGNIFICANT_CHANGE + FLAT]
TCCI_by_REGULATION[, .N, by = c("TOTAL")]

TCCI_by_REGULATION[TOTAL == 2, .N, by = c("DOWN", "FLAT", "UP", "NON_SIGNIFICANT_CHANGE")][order(-N)]
TCCI_by_REGULATION[TOTAL == 1, .N, by = c("DOWN", "FLAT", "UP", "NON_SIGNIFICANT_CHANGE")][order(-N)]

droplet_all_cci[ID == "Heart_and_Aorta"]

TCCI_by_REGULATION <- dcast.data.table(
  droplet_all_cci[ID == "Heart_and_Aorta", .N, by = c("REGULATION", "TCCI", "dataset")],
  formula = TCCI ~ REGULATION + dataset,
  value.var = "N"
)
test <- TCCI_by_REGULATION[, .N, by = c(colnames(TCCI_by_REGULATION)[-1])]

facs_all_cci <- rbindlist(
  lapply(
    DATASETS_PROCESSED[c(6, 7, 8)],
    function(i) {
      i$dataset@cci_detected
    }
  ),
  idcol = "dataset"
)
facs_all_cci[, TCCI := paste(ID, CCI, sep = "_")]
facs_pct_comp <- facs_pct_comp[, .N, by = c("dataset", "ID", "REGULATION")]
facs_pct_comp[, PCT := N/sum(N)*100, by = c("dataset", "ID") ]
facs_pct_comp <- dcast.data.table(
  facs_pct_comp,
  formula = dataset + ID ~ REGULATION,
  value.var = "PCT"
)

##########

cci_droplet_mixed <- copy(DATASETS_PROCESSED$droplet_mixed$dataset@cci_detected)
cci_droplet_mixed[, TCCI := paste(ID, CCI, sep = "_")]
cci_droplet_sex <- copy(DATASETS_PROCESSED$droplet_sex$dataset@cci_detected)
cci_droplet_sex[, TCCI := paste(ID, CCI, sep = "_")]
cci_droplet_female <- copy(DATASETS_PROCESSED$droplet_female$dataset@cci_detected)
cci_droplet_female[, TCCI := paste(ID, CCI, sep = "_")]
cci_droplet_male <- copy(DATASETS_PROCESSED$droplet_male$dataset@cci_detected)
cci_droplet_male[, TCCI := paste(ID, CCI, sep = "_")]


cci_droplet_comp <- merge.data.table(
  x = cci_droplet_mixed[, c("TCCI", "REGULATION")],
  y = cci_droplet_sex[, c("TCCI", "REGULATION")],
  by = "TCCI",
  all.x = TRUE,
  all.y = TRUE,
  suffixes = c("_mixed", "_sex")
)
cci_droplet_comp[is.na(cci_droplet_comp)] <- "NON_DETECTED"

cci_droplet_comp_mf <- merge.data.table(
  x = cci_droplet_female[, c("TCCI", "REGULATION")],
  y = cci_droplet_male[, c("TCCI", "REGULATION")],
  by = "TCCI",
  all.x = TRUE,
  all.y = TRUE,
  suffixes = c("_female", "_male")
)
cci_droplet_comp_mf[is.na(cci_droplet_comp_mf)] <- "NON_DETECTED"

sex_overlap_droplet <- cci_droplet_comp[, .N, by = c("REGULATION_mixed", "REGULATION_sex")]
mf_overlap_droplet <- cci_droplet_comp_mf[, .N, by = c("REGULATION_female", "REGULATION_male")]

cci_facs_mixed <- copy(DATASETS_PROCESSED$facs_mixed$dataset@cci_detected)
cci_facs_mixed[, TCCI := paste(ID, CCI, sep = "_")]
cci_facs_sex <- copy(DATASETS_PROCESSED$facs_sex$dataset@cci_detected)
cci_facs_sex[, TCCI := paste(ID, CCI, sep = "_")]
cci_facs_female <- copy(DATASETS_PROCESSED$facs_female$dataset@cci_detected)#[ID == "Lung"])
cci_facs_female[, TCCI := paste(ID, CCI, sep = "_")]
cci_facs_male <- copy(DATASETS_PROCESSED$facs_male$dataset@cci_detected)#[ID == "Lung"])
cci_facs_male[, TCCI := paste(ID, CCI, sep = "_")]


cci_facs_comp <- merge.data.table(
  x = cci_facs_mixed[, c("TCCI", "REGULATION")],
  y = cci_facs_sex[, c("TCCI", "REGULATION")],
  by = "TCCI",
  all.x = TRUE,
  all.y = TRUE,
  suffixes = c("_mixed", "_sex")
)
cci_facs_comp[is.na(cci_facs_comp)] <- "NON_DETECTED"

cci_facs_comp_MF <- merge.data.table(
  x = cci_facs_mixed[, c("TCCI", "REGULATION")],
  y = cci_facs_female[, c("TCCI", "REGULATION")],
  by = "TCCI",
  all.x = TRUE,
  all.y = TRUE,
  suffixes = c("_mixed", "_female")
)
cci_facs_comp_MF[is.na(cci_facs_comp_MF)] <- "NON_DETECTED"

sex_overlap_facs <- cci_facs_comp[, .N, by = c("REGULATION_mixed", "REGULATION_sex")]
MF_overlap_facs <- cci_facs_comp_MF[, .N, by = c("REGULATION_mixed", "REGULATION_female")]