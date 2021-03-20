

## Part on the server logfc per cell-type in function of sex ####

#actually it could be done from the raw data directly

library(Seurat)
library(scDiffCom)
library(data.table)


seurat_facs <- readRDS("seurat_final_tms_facs.rds")
seurat_facs$tissue_celltype <- paste(seurat_facs$tissue, seurat_facs$cell_ontology_final, sep = "_")
all_tissue_celltype <- sort(unique(seurat_facs$tissue_celltype))
LRI_genes <- sort(unique(unlist(LRI_mouse$LRI_curated[, 2:6])))
LRI_genes <- LRI_genes[!is.na(LRI_genes)]
LRI_genes <- LRI_genes[LRI_genes %in% rownames(seurat_facs)]
seurat_facs_sub <- subset(seurat_facs, features = LRI_genes)
seurat_facs_sub$age_group <- ifelse(
  seurat_facs_sub$age == "3m",
  "YOUNG",
  "OLD"
)
seurat_data <- expm1(GetAssayData(seurat_facs_sub, assay = "RNA", slot = "data"))
logfc_celltype_facs <- rbindlist(
  sapply(
    all_tissue_celltype,
    function(ct) {
      cells_tokeep_young_male <- colnames(seurat_data)[
        seurat_facs_sub$tissue_celltype == ct &
          seurat_facs_sub$age_group == "YOUNG" &
          seurat_facs_sub$sex == "male"
      ]
      cells_tokeep_young_female <- colnames(seurat_data)[
        seurat_facs_sub$tissue_celltype == ct &
          seurat_facs_sub$age_group == "YOUNG" &
          seurat_facs_sub$sex == "female"
      ]
      cells_tokeep_old_male <- colnames(seurat_data)[
        seurat_facs_sub$tissue_celltype == ct &
          seurat_facs_sub$age_group == "OLD" &
          seurat_facs_sub$sex == "male"
      ]
      cells_tokeep_old_female <- colnames(seurat_data)[
        seurat_facs_sub$tissue_celltype == ct &
          seurat_facs_sub$age_group == "OLD" &
          seurat_facs_sub$sex == "female"
      ]
      if (length(cells_tokeep_old_male) < 5 | length(cells_tokeep_young_male) < 5) {
        res_male <- NULL
      } else {
        mat_young_male <- seurat_data[, cells_tokeep_young_male]
        mat_old_male <- seurat_data[, cells_tokeep_old_male]
        dr_young_male <- rowMeans(mat_young_male > 0)
        dr_old_male <- rowMeans(mat_old_male > 0)
        logfc_male <- log(rowMeans(mat_old_male)/rowMeans(mat_young_male))
        logfc_male[dr_young_male < 0.1 & dr_old_male < 0.1] <- 0
        res_male <- data.table(
          GENE = names(logfc_male),
          LOGFC = logfc_male,
          GENDER = "male")
      }
      if (length(cells_tokeep_old_female) < 5 | length(cells_tokeep_young_female) < 5) {
        res_female <- NULL
      } else {
        mat_young_female <- seurat_data[, cells_tokeep_young_female]
        mat_old_female <- seurat_data[, cells_tokeep_old_female]
        dr_young_female <- rowMeans(mat_young_female > 0)
        dr_old_female <- rowMeans(mat_old_female > 0)
        logfc_female <- log(rowMeans(mat_old_female)/rowMeans(mat_young_female))
        logfc_female[dr_young_female < 0.1 & dr_old_female < 0.1] <- 0
        res_female <- data.table(
          GENE = names(logfc_female),
          LOGFC = logfc_female,
          GENDER = "female")
      }
      return(rbind(res_male, res_female))
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  ),
  idcol = "CELL_TYPE"
)

seurat_droplet <- readRDS("seurat_final_tms_droplet.rds")
seurat_droplet$tissue_celltype <- paste(seurat_droplet$tissue, seurat_droplet$cell_ontology_final, sep = "_")
all_tissue_celltype <- sort(unique(seurat_droplet$tissue_celltype))
LRI_genes <- sort(unique(unlist(LRI_mouse$LRI_curated[, 2:6])))
LRI_genes <- LRI_genes[!is.na(LRI_genes)]
LRI_genes <- LRI_genes[LRI_genes %in% rownames(seurat_droplet)]
seurat_droplet_sub <- subset(seurat_droplet, features = LRI_genes)
seurat_droplet_sub$age_group <- ifelse(
  seurat_droplet_sub$age == "3m",
  "YOUNG",
  ifelse (
    seurat_droplet_sub$age %in% c("18m", "21m", "24m"),
    "OLD",
    "OTHER"
  )
)
seurat_data <- expm1(GetAssayData(seurat_droplet_sub, assay = "RNA", slot = "data"))

logfc_celltype_droplet <- rbindlist(
  sapply(
    all_tissue_celltype,
    function(ct) {
      cells_tokeep_young_male <- colnames(seurat_data)[
        seurat_droplet_sub$tissue_celltype == ct &
          seurat_droplet_sub$age_group == "YOUNG" &
          seurat_droplet_sub$sex == "male"
      ]
      cells_tokeep_young_female <- colnames(seurat_data)[
        seurat_droplet_sub$tissue_celltype == ct &
          seurat_droplet_sub$age_group == "YOUNG" &
          seurat_droplet_sub$sex == "female"
      ]
      cells_tokeep_old_male <- colnames(seurat_data)[
        seurat_droplet_sub$tissue_celltype == ct &
          seurat_droplet_sub$age_group == "OLD" &
          seurat_droplet_sub$sex == "male"
      ]
      cells_tokeep_old_female <- colnames(seurat_data)[
        seurat_droplet_sub$tissue_celltype == ct &
          seurat_droplet_sub$age_group == "OLD" &
          seurat_droplet_sub$sex == "female"
      ]
      if (length(cells_tokeep_old_male) < 5 | length(cells_tokeep_young_male) < 5) {
        res_male <- NULL
      } else {
        mat_young_male <- seurat_data[, cells_tokeep_young_male]
        mat_old_male <- seurat_data[, cells_tokeep_old_male]
        dr_young_male <- rowMeans(mat_young_male > 0)
        dr_old_male <- rowMeans(mat_old_male > 0)
        logfc_male <- log(rowMeans(mat_old_male)/rowMeans(mat_young_male))
        logfc_male[dr_young_male < 0.1 & dr_old_male < 0.1] <- 0
        res_male <- data.table(
          GENE = names(logfc_male),
          LOGFC = logfc_male,
          GENDER = "male")
      }
      if (length(cells_tokeep_old_female) < 5 | length(cells_tokeep_young_female) < 5) {
        res_female <- NULL
      } else {
        mat_young_female <- seurat_data[, cells_tokeep_young_female]
        mat_old_female <- seurat_data[, cells_tokeep_old_female]
        dr_young_female <- rowMeans(mat_young_female > 0)
        dr_old_female <- rowMeans(mat_old_female > 0)
        logfc_female <- log(rowMeans(mat_old_female)/rowMeans(mat_young_female))
        logfc_female[dr_young_female < 0.1 & dr_old_female < 0.1] <- 0
        res_female <- data.table(
          GENE = names(logfc_female),
          LOGFC = logfc_female,
          GENDER = "female")
      }
      return(rbind(res_male, res_female))
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  ),
  idcol = "CELL_TYPE"
)


## read genes logfc computed on the server ####

logfc_celltype_facs <- readRDS("../data_scAgeCom/analysis/outputs_data/tms_facs_logfc_LRIGenes_by_celltypes.rds")
logfc_celltype_droplet <- readRDS("../data_scAgeCom/analysis/outputs_data/tms_droplet_logfc_LRIGenes_by_celltypes.rds")

logfc_ct_facs <- dcast.data.table(
  logfc_celltype_facs,
  GENE + CELL_TYPE ~ GENDER,
  value.var = "LOGFC"
)
logfc_ct_facs[is.na(logfc_ct_facs)] <- 0
logfc_ct_facs_detected <- logfc_ct_facs[female !=0 | male != 0]
logfc_ct_facs_detected2 <- logfc_ct_facs[female !=0 & male != 0]
logfc_ct_facs_detected3 <- logfc_ct_facs_detected2[is.finite(male) & is.finite(female)]

logfc_ct_droplet <- dcast.data.table(
  logfc_celltype_droplet,
  GENE + CELL_TYPE ~ GENDER,
  value.var = "LOGFC"
)
logfc_ct_droplet[is.na(logfc_ct_droplet)] <- 0
logfc_ct_droplet_detected <- logfc_ct_droplet[female !=0 | male != 0]
logfc_ct_droplet_detected2 <- logfc_ct_droplet[female !=0 & male != 0]
logfc_ct_droplet_detected3 <- logfc_ct_droplet_detected2[is.finite(male) & is.finite(female)]

cor(logfc_ct_facs_detected3$female, logfc_ct_facs_detected3$male, method = "pearson")
table(logfc_ct_facs_detected2$male < - log(1.5), logfc_ct_facs_detected2$female > log(1.5))
ggplot(logfc_ct_facs_detected2, aes(male, female)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
ggplot(logfc_ct_facs_detected2[grepl("App", GENE)], aes(male, female)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

cor(logfc_ct_droplet_detected3$female, logfc_ct_droplet_detected3$male, method = "pearson")
table(logfc_ct_droplet_detected2$male > 0, logfc_ct_droplet_detected2$female < 0)
ggplot(logfc_ct_droplet_detected2, aes(male, female)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
ggplot(logfc_ct_droplet_detected2[grepl("App", GENE)], aes(male, female)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

## Look at LRI genes fc based on detected CCI only #####

DATASETS_PROCESSED_GENDER <- readRDS("../data_scAgeCom/analysis/outputs_data/scAgeCom_results_processed.rds")

get_genes_logfc <- function(
  cci_dt
) {
  L1 <- unique(cci_dt[, c("ID", "EMITTER_CELLTYPE", "LIGAND_1", "L1_EXPRESSION_YOUNG", "L1_EXPRESSION_OLD")])
  L1 <- na.omit(L1)
  L1[, TISSUE_CELLTYPE := paste(ID, EMITTER_CELLTYPE, sep = "_")]
  L1[, GENE := LIGAND_1]
  L1[, LOGFC := log(L1_EXPRESSION_OLD/L1_EXPRESSION_YOUNG)]
  L1 <- L1[, c("TISSUE_CELLTYPE", "GENE", "LOGFC")]
  
  L2 <- unique(cci_dt[, c("ID", "EMITTER_CELLTYPE", "LIGAND_2", "L2_EXPRESSION_YOUNG", "L2_EXPRESSION_OLD")])
  L2 <- na.omit(L2)
  L2[, TISSUE_CELLTYPE := paste(ID, EMITTER_CELLTYPE, sep = "_")]
  L2[, GENE := LIGAND_2]
  L2[, LOGFC := log(L2_EXPRESSION_OLD/L2_EXPRESSION_YOUNG)]
  L2 <- L2[, c("TISSUE_CELLTYPE", "GENE", "LOGFC")]
  
  R1 <- unique(cci_dt[, c("ID", "RECEIVER_CELLTYPE", "RECEPTOR_1", "R1_EXPRESSION_YOUNG", "R1_EXPRESSION_OLD")])
  R1 <- na.omit(R1)
  R1[, TISSUE_CELLTYPE := paste(ID, RECEIVER_CELLTYPE, sep = "_")]
  R1[, GENE := RECEPTOR_1]
  R1[, LOGFC := log(R1_EXPRESSION_OLD/R1_EXPRESSION_YOUNG)]
  R1 <- R1[, c("TISSUE_CELLTYPE", "GENE", "LOGFC")]
  
  R2 <- unique(cci_dt[, c("ID", "RECEIVER_CELLTYPE", "RECEPTOR_2", "R2_EXPRESSION_YOUNG", "R2_EXPRESSION_OLD")])
  R2 <- na.omit(R2)
  R2[, TISSUE_CELLTYPE := paste(ID, RECEIVER_CELLTYPE, sep = "_")]
  R2[, GENE := RECEPTOR_2]
  R2[, LOGFC := log(R2_EXPRESSION_OLD/R2_EXPRESSION_YOUNG)]
  R2 <- R2[, c("TISSUE_CELLTYPE", "GENE", "LOGFC")]
  
  R3 <- unique(cci_dt[, c("ID", "RECEIVER_CELLTYPE", "RECEPTOR_3", "R3_EXPRESSION_YOUNG", "R3_EXPRESSION_OLD")])
  R3 <- na.omit(R3)
  R3[, TISSUE_CELLTYPE := paste(ID, RECEIVER_CELLTYPE, sep = "_")]
  R3[, GENE := RECEPTOR_3]
  R3[, LOGFC := log(R3_EXPRESSION_OLD/R3_EXPRESSION_YOUNG)]
  R3 <- R3[, c("TISSUE_CELLTYPE", "GENE", "LOGFC")]
  
  dt <- unique(rbind(L1, L2, R1, R2, R3))
  
}

fc_df <- get_genes_logfc(DATASETS_PROCESSED_GENDER$droplet_female$dataset@cci_detected)
fc_dm <- get_genes_logfc(DATASETS_PROCESSED_GENDER$droplet_male$dataset@cci_detected)
fc_d <- merge.data.table(
  fc_df,
  fc_dm,
  all.x = TRUE,
  all.y = TRUE,
  by = c("TISSUE_CELLTYPE", "GENE"),
  suffixes = c("_female", "_male")
)
fc_d <- na.omit(fc_d)
ggplot(fc_d, aes(LOGFC_male, LOGFC_female)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
cor(
  fc_d[is.finite(LOGFC_female) & is.finite(LOGFC_male)]$LOGFC_female,
  fc_d[is.finite(LOGFC_female) & is.finite(LOGFC_male)]$LOGFC_male
)

fc_ff <- get_genes_logfc(DATASETS_PROCESSED_GENDER$facs_female$dataset@cci_detected)
fc_fm <- get_genes_logfc(DATASETS_PROCESSED_GENDER$facs_male$dataset@cci_detected)
fc_f <- merge.data.table(
  fc_ff,
  fc_fm,
  all.x = TRUE,
  all.y = TRUE,
  by = c("TISSUE_CELLTYPE", "GENE"),
  suffixes = c("_female", "_male")
)
fc_f <- na.omit(fc_f)
ggplot(fc_f, aes(LOGFC_male, LOGFC_female)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
cor(
  fc_f[is.finite(LOGFC_female) & is.finite(LOGFC_male)]$LOGFC_female,
  fc_f[is.finite(LOGFC_female) & is.finite(LOGFC_male)]$LOGFC_male
)

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

seurat_md$tms_droplet <- seurat_md$tms_droplet[age != "1m"]

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