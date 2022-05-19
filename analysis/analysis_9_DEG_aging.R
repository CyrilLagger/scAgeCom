####################################################
##
## Project: scAgeCom
##
## Last update - April 2022
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## Perform DEG analysis on TMS data
## that match scAgeCom cell types and conditions
##
####################################################
##

## Libraries ####

library(Seurat)
library(scDiffCom)
library(future)
library(future.apply)
library(data.table)
library(ggplot2)

options(future.globals.maxSize = 15 * 1024^3)

## Dataset path on Server ####

server_path <- "/workspace/lcyril_data/scRNA_seq/seurat_processed/"

dataset_path <- c(
  tms_facs = paste0(server_path, "seurat_final_tms_facs.rds"),
  tms_droplet = paste0(server_path, "seurat_final_tms_droplet.rds"),
  calico_kidney = paste0(server_path, "seurat_final_calico_kidney.rds"),
  calico_lung = paste0(server_path, "seurat_final_calico_lung.rds"),
  calico_spleen = paste0(server_path, "seurat_final_calico_spleen.rds")
)

## Output path ####

dir_deg_output <- paste0(getwd(), "/data_deg_wilcox_aging")

## List of analyses to do over datasets and sex ####

analysis_list <- list(
  "calico_male",
  "facs_female",
  "facs_male",
  "facs_mixed",
  "droplet_female",
  "droplet_male",
  "droplet_mixed"
)

## Do the DEG analysis ####

for (analysis in analysis_list) {
  message(paste0("Starting DEG analysis of ", analysis))
  output_dir <- paste0(
    dir_deg_output,
    "/deg_",
    analysis
  )
  if (!dir.exists(output_dir)) {
    message("Create output directory.")
    dir.create(output_dir)
  }
  message("Reading seurat object.")
  if (analysis == "calico_male") {
    seurat_kidney <- readRDS(dataset_path[["calico_kidney"]])
    seurat_kidney$age_group <- ifelse(
      seurat_kidney$age == "young",
      "YOUNG",
      "OLD"
    )
    seurat_lung <- readRDS(dataset_path[["calico_lung"]])
    seurat_lung$age_group <- ifelse(
      seurat_lung$age == "young",
      "YOUNG",
      "OLD"
    )
    seurat_spleen <- readRDS(dataset_path[["calico_spleen"]])
    seurat_spleen$age_group <- ifelse(
      seurat_spleen$age == "young",
      "YOUNG",
      "OLD"
    )
    tissue_list <- c("Kidney", "Lung", "Spleen")
    seurat_list <- list(
      Kidney = seurat_kidney,
      Lung = seurat_lung,
      Spleen = seurat_spleen
    )
    n_tissue <- 3
  } else if (analysis == "facs_female") {
    seurat_obj <- readRDS(dataset_path[["tms_facs"]])
    seurat_obj <- subset(seurat_obj, subset = sex == "female")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "BAT", "Bladder", "Brain", "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle",
      "Lung", "MAT", "Mammary_Gland", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "facs_male") {
    seurat_obj <- readRDS(dataset_path[["tms_facs"]])
    seurat_obj <- subset(seurat_obj, subset = sex == "male")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Aorta", "BAT", "Bladder", "Brain", "Diaphragm",
      "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle", "Liver",
      "Lung", "MAT", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "facs_mixed") {
    seurat_obj <- readRDS(dataset_path[["tms_facs"]])
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Aorta", "BAT", "Bladder", "Brain", "Diaphragm",
      "GAT", "Heart", "Kidney",
      "Large_Intestine", "Limb_Muscle", "Liver",
      "Lung", "Mammary_Gland", "MAT", "Marrow",
      "Pancreas", "SCAT", "Skin", "Spleen", "Thymus",
      "Tongue", "Trachea"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_female") {
    seurat_obj <- readRDS(dataset_path[["tms_droplet"]])
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj <- subset(seurat_obj, subset = sex == "female")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Heart_and_Aorta", "Kidney", "Limb_Muscle",
      "Liver", "Lung", "Mammary_Gland", "Marrow",
      "Spleen", "Thymus"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_male") {
    seurat_obj <- readRDS(dataset_path[["tms_droplet"]])
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj <- subset(seurat_obj, subset = sex == "male")
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Bladder", "Kidney", "Liver", "Lung",
      "Spleen", "Tongue"
    )
    n_tissue <- length(tissue_list)
  } else if (analysis == "droplet_mixed") {
    seurat_obj <- readRDS(dataset_path[["tms_droplet"]])
    seurat_obj <- subset(
      seurat_obj,
      subset = age %in% c("3m", "18m", "21m", "24m")
    )
    seurat_obj$age_group <- ifelse(
      seurat_obj$age %in% c("3m"),
      "YOUNG",
      "OLD"
    )
    tissue_list <- c(
      "Bladder", "Heart_and_Aorta", "Kidney", "Liver", "Limb_Muscle",
      "Liver", "Lung", "Mammary_Gland", "Marrow",
      "Spleen", "Tongue", "Thymus"
    )
    n_tissue <- length(tissue_list)
  }
  for (i in 1:n_tissue) {
    tiss <- tissue_list[i]
    message(paste0(
      "DEG analysis of the ",
      tiss,
      ". Tissue ",
      i,
      " out of ",
      n_tissue,
      "."
    ))
    if (analysis == "calico_male") {
      seurat_tiss <- seurat_list[[tiss]]
    } else {
      if (tiss == "Brain") {
        seurat_tiss <- subset(
          seurat_obj,
          subset = tissue %in% c(
            "Brain_Myeloid",
            "Brain_Non-Myeloid"
          )
        )
      } else {
        cells_tiss <- colnames(seurat_obj)[which(seurat_obj$tissue == tiss)]
        seurat_tiss <- subset(seurat_obj, cells = cells_tiss)
      }
    }
    cell_type_id <- "cell_ontology_final"
    n_ct <- length(unique(seurat_tiss[[cell_type_id]][[1]]))
    DefaultAssay(seurat_tiss) <- "RNA"
    if (n_ct > 1) {
      message("Size-factor normalization:")
      seurat_tiss <- NormalizeData(seurat_tiss, assay = "RNA")
      Idents(seurat_tiss) <- seurat_tiss[[cell_type_id]]
      saveRDS(seurat_tiss[[]], file = paste0(output_dir, "/md_", tiss, ".rds"))
      future::plan(multicore, workers = 24)
      deg_res <- lapply(
        levels(seurat_tiss),
        function(ct_temp) {
          tryCatch(
            FindMarkers(
              object = seurat_tiss,
              ident.1 = "OLD",
              ident.2 = "YOUNG",
              group.by = "age_group",
              subset.ident = ct_temp,
              logfc.threshold = log2(1.2),
              test.use = "wilcox",
              min.pct = 0.1,
              min.cells.feature = 5,
              min.cells.group = 5,
              pseudocount.use = 0.001
            ),
            error = function(e) e,
            warning = function(w) w
          )
        }
      )
      names(deg_res) <- levels(seurat_tiss)
      message(paste0("Saving results for the ", tiss, "."))
      saveRDS(deg_res, file = paste0(output_dir, "/deg_", tiss, ".rds"))
    } else {
      message("Not enough cell types, not performing the analysis.")
    }
    future::plan(sequential)
  }
}

## Path of MAST/wilcox results ####

MAST_path <- "data_deg_mast_aging"
MAST_dataset_paths <- list.dirs(MAST_path, recursive = FALSE)
MAST_dataset_names <- c(
  "Calico Droplet (male)",
  "TMS Droplet (female)",
  "TMS Droplet (male)",
  "TMS Droplet (mixed)",
  "TMS FACS (female)",
  "TMS FACS (male)",
  "TMS FACS (mixed)"
)
names(MAST_dataset_names) <- MAST_dataset_names

wilcox_path <- "data_deg_wilcox_aging"
wilcox_dataset_paths <- list.dirs(wilcox_path, recursive = FALSE)
wilcox_dataset_names <- c(
  "Calico Droplet (male)",
  "TMS Droplet (female)",
  "TMS Droplet (male)",
  "TMS Droplet (mixed)",
  "TMS FACS (female)",
  "TMS FACS (male)",
  "TMS FACS (mixed)"
)
names(wilcox_dataset_names) <- wilcox_dataset_names

## Clean MAST/wilcox results in a single data.table ####

process_deg_dataset <- function(
  deg_dataset_path
) {
  # retrieve and load each object in the dataset
  tissues <- gsub(".*deg_(.+)\\.rds.*", "\\1", list.files(deg_dataset_path))
  tissues <- tissues[!grepl("md_", tissues)]
  dataset <- lapply(
    X = tissues,
    FUN = function(
      tiss
    ) {
      readRDS(paste0(deg_dataset_path, "/deg_", tiss, ".rds"))
    }
  )
  tissues[tissues == "BAT"] <- "Adipose_Brown"
  tissues[tissues == "GAT"] <- "Adipose_Gonadal"
  tissues[tissues == "MAT"] <- "Adipose_Mesenteric"
  tissues[tissues == "SCAT"] <- "Adipose_Subcutaneous"
  names(dataset) <- tissues
  # remove potential errors and transform to data.table
  dataset <- rbindlist(
    l = lapply(
      dataset,
      function(tiss) {
        rbindlist(
          l = lapply(
            tiss,
            function(ct) {
              if (!is.data.frame(ct)) {
                return(NULL)
              } else {
                as.data.table(ct, keep.rownames = "gene")
              }
            }
          ),
          idcol = "cell_type"
        )
      }
    ),
    idcol = "tissue"
  )
  return(dataset)
}

MAST_results <- rbindlist(
  l = lapply(
    setNames(
      MAST_dataset_paths,
      MAST_dataset_names
    ),
    process_deg_dataset
  ),
  idcol = "dataset"
)
setnames(
  MAST_results,
  old = c("pct.1", "pct.2"),
  new = c("pct_old", "pct_young")
)
saveRDS(MAST_results, "data_deg_mast_aging/MAST_results_dt.rds")

wilcox_results <- rbindlist(
  l = lapply(
    setNames(
      wilcox_dataset_paths,
      wilcox_dataset_names
    ),
    process_deg_dataset
  ),
  idcol = "dataset"
)
setnames(
  wilcox_results,
  old = c("pct.1", "pct.2"),
  new = c("pct_old", "pct_young")
)
saveRDS(wilcox_results, "data_deg_wilcox_aging/wilcox_results_dt.rds")

## Perform FDR correction by tissue-cell-type ####

wilcox_results[
  ,
  p_val_bh := p.adjust(p_val, method = "BH"),
  by = c("dataset", "tissue", "cell_type")
]

## DEG regulation ####

MAST_results[
  ,
  regulation := ifelse(
    p_val <= 0.05 & avg_log2FC > log2(1.5),
    "UP",
    ifelse(
      p_val <= 0.05 & avg_log2FC < -log2(1.5),
      "DOWN",
      "NO"
    )
  )
]
table(MAST_results$regulation)

wilcox_results[
  ,
  regulation := ifelse(
    p_val_bh <= 0.05 & avg_log2FC > log2(1.5),
    "UP",
    ifelse(
      p_val_bh <= 0.05 & avg_log2FC < -log2(1.5),
      "DOWN",
      "NO"
    )
  )
]
table(wilcox_results$regulation)

## Load scAgeCom results ####

scAgeCom_results <- readRDS(
  "data_scAgeCom_11_04_2022_processed/scAgeCom_results_processed.rds"
)

CCI_results <- rbindlist(
  lapply(
    scAgeCom_results,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            GetTableCCI(
              object = tissue,
              type = "detected",
              simplified = FALSE
            )
          }
        ),
        fill = TRUE,
        idcol = "Tissue"
      )
    }
  ),
  idcol = "Dataset"
)

## Comparing wilcox DEG with scAgeCom ####

CCI_results[
  wilcox_results,
  on = c(
    "Dataset==dataset", "Tissue==tissue",
    "EMITTER_CELLTYPE==cell_type", "LIGAND_1==gene"
  ),
  l1_deg_reg := i.regulation
]
CCI_results[
  ,
  l1_deg_reg := ifelse(
    is.na(l1_deg_reg),
    "NO",
    l1_deg_reg
  )
]
CCI_results[
  wilcox_results,
  on = c(
    "Dataset==dataset", "Tissue==tissue",
    "EMITTER_CELLTYPE==cell_type", "LIGAND_2==gene"
  ),
  l2_deg_reg := i.regulation
]
CCI_results[
  ,
  l2_deg_reg := ifelse(
    is.na(l2_deg_reg) & !is.na(LIGAND_2),
    "NO",
    l2_deg_reg
  )
]
CCI_results[
  wilcox_results,
  on = c(
    "Dataset==dataset", "Tissue==tissue",
    "RECEIVER_CELLTYPE==cell_type", "RECEPTOR_1==gene"
  ),
  r1_deg_reg := i.regulation
]
CCI_results[
  ,
  r1_deg_reg := ifelse(
    is.na(r1_deg_reg),
    "NO",
    r1_deg_reg
  )
]
CCI_results[
  wilcox_results,
  on = c(
    "Dataset==dataset", "Tissue==tissue",
    "RECEIVER_CELLTYPE==cell_type", "RECEPTOR_2==gene"
  ),
  r2_deg_reg := i.regulation
]
CCI_results[
  ,
  r2_deg_reg := ifelse(
    is.na(r2_deg_reg) & !is.na(RECEPTOR_2),
    "NO",
    r2_deg_reg
  )
]
CCI_results[
  wilcox_results,
  on = c(
    "Dataset==dataset", "Tissue==tissue",
    "RECEIVER_CELLTYPE==cell_type", "RECEPTOR_3==gene"
  ),
  r3_deg_reg := i.regulation
]
CCI_results[
  ,
  r3_deg_reg := ifelse(
    is.na(r3_deg_reg) & !is.na(RECEPTOR_3),
    "NO",
    r3_deg_reg
  )
]

table(CCI_results$l1_deg_reg)
table(CCI_results$l2_deg_reg)
table(CCI_results$r1_deg_reg)
table(CCI_results$r2_deg_reg)
table(CCI_results$r3_deg_reg)
anyNA(CCI_results$r3_deg_reg)

## Summary Wilcox comparison table ####

CCI_results[
  ,
  scDiffCom_regulation := ifelse(
    !REGULATION %in% c("UP", "DOWN"),
    "NO",
    REGULATION
  )
]

deg_comp_wilcox <- CCI_results[, .N, by = c(
  "scDiffCom_regulation",
  "l1_deg_reg",
  "l2_deg_reg",
  "r1_deg_reg",
  "r2_deg_reg",
  "r3_deg_reg"
)
]
deg_comp_wilcox[
  ,
  pct := N / sum(N) * 100
]

fwrite(deg_comp_wilcox, "data_deg_wilcox_aging/deg_comparison.csv")

## Compare permutation DEG with CCIs ####

CCI_results[
  ,
  l1_deg_perm :=  ifelse(
    L1_BH_P_VALUE_DE <= 0.05 & L1_LOGFC > log(1.5),
    "UP",
    ifelse(
      L1_BH_P_VALUE_DE <= 0.05 & L1_LOGFC < -log(1.5),
      "DOWN",
      "NO"
    )
  )
]
table(CCI_results$l1_deg_perm)
table(CCI_results$l1_deg_perm, CCI_results$l1_deg_reg)

CCI_results[
  ,
  l2_deg_perm :=  ifelse(
    is.na(L2_BH_P_VALUE_DE),
    NA,
    ifelse(
      L2_BH_P_VALUE_DE <= 0.05 & L2_LOGFC > log(1.5),
      "UP",
      ifelse(
        L2_BH_P_VALUE_DE <= 0.05 & L2_LOGFC < -log(1.5),
        "DOWN",
        "NO"
      )
    )
  )
]
table(CCI_results$l2_deg_perm)
table(CCI_results$l2_deg_perm, CCI_results$l2_deg_reg)

CCI_results[
  ,
  r1_deg_perm :=  ifelse(
    R1_BH_P_VALUE_DE <= 0.05 & R1_LOGFC > log(1.5),
    "UP",
    ifelse(
      R1_BH_P_VALUE_DE <= 0.05 & R1_LOGFC < -log(1.5),
      "DOWN",
      "NO"
    )
  )
]
table(CCI_results$r1_deg_perm)
table(CCI_results$r1_deg_perm, CCI_results$r1_deg_reg)

CCI_results[
  ,
  r2_deg_perm :=  ifelse(
    is.na(R2_BH_P_VALUE_DE),
    NA,
    ifelse(
      R2_BH_P_VALUE_DE <= 0.05 & R2_LOGFC > log(1.5),
      "UP",
      ifelse(
        R2_BH_P_VALUE_DE <= 0.05 & R2_LOGFC < -log(1.5),
        "DOWN",
        "NO"
      )
    )
  )
]
table(CCI_results$r2_deg_perm)
table(CCI_results$r2_deg_perm, CCI_results$r2_deg_reg)

CCI_results[
  ,
  r3_deg_perm :=  ifelse(
    is.na(R3_BH_P_VALUE_DE),
    NA,
    ifelse(
      R3_BH_P_VALUE_DE <= 0.05 & R3_LOGFC > log(1.5),
      "UP",
      ifelse(
        R3_BH_P_VALUE_DE <= 0.05 & R3_LOGFC < -log(1.5),
        "DOWN",
        "NO"
      )
    )
  )
]
table(CCI_results$r3_deg_perm)
table(CCI_results$r3_deg_perm, CCI_results$r3_deg_reg)

## DEG permutation comparison ####

cci_deg_comp <- CCI_results[
  ,
  c(
  "Dataset",
  "scDiffCom_regulation",
  "l1_deg_perm",
  "l2_deg_perm",
  "r1_deg_perm",
  "r2_deg_perm",
  "r3_deg_perm"
  )
]

cci_deg_comp[
  ,
  ligand_deg_perm := ifelse(
    is.na(l2_deg_perm),
    l1_deg_perm,
    ifelse(
      l1_deg_perm == l2_deg_perm,
      l1_deg_perm,
      ifelse(
        l1_deg_perm == "NO",
        l2_deg_perm,
        ifelse(
          l2_deg_perm == "NO",
          l1_deg_perm,
          "NO"
        )
      )
    )
  )
]

cci_deg_comp[
  ,
  receptor_deg_perm := ifelse(
    is.na(r2_deg_perm),
    r1_deg_perm,
    ifelse(
      r1_deg_perm == r2_deg_perm,
      r1_deg_perm,
      ifelse(
        r1_deg_perm == "NO",
        r2_deg_perm,
        ifelse(
          r2_deg_perm == "NO",
          r1_deg_perm,
          "NO"
        )
      )
    )
  )
]

table(cci_deg_comp$ligand_deg_perm)
table(cci_deg_comp$receptor_deg_perm)
anyNA(cci_deg_comp$receptor_deg_perm)

deg_comp_perm <- cci_deg_comp[, .N, by = c(
  "Dataset",
  "scDiffCom_regulation",
  "ligand_deg_perm",
  "receptor_deg_perm"
)
]
deg_comp_perm[
  ,
  pct := N / sum(N) * 100,
  by = c("Dataset")
]

fwrite(
  deg_comp_perm,
   "data_scAgeCom_11_04_2022_processed/deg_comparison_perm.csv"
)
deg_comp_perm <- fread("data_scAgeCom_11_04_2022_processed/deg_comparison_perm.csv")

## Comparison plot SDEA vs scDiffCom ####

deg_comp_perm_test <- deg_comp_perm[
  Dataset == "TMS FACS (male)"
]

cci_deg_comp[, scdiffcom_up := scDiffCom_regulation == "UP"]
cci_deg_comp[, scdiffcom_down := scDiffCom_regulation == "DOWN"]
cci_deg_comp[, scdiffcom_no := scDiffCom_regulation == "NO"]
cci_deg_comp[, ligand_up := ligand_deg_perm == "UP"]
cci_deg_comp[, ligand_down := ligand_deg_perm == "DOWN"]
cci_deg_comp[, ligand_no := ligand_deg_perm == "NO"]
cci_deg_comp[, receptor_up := receptor_deg_perm == "UP"]
cci_deg_comp[, receptor_down := receptor_deg_perm == "DOWN"]
cci_deg_comp[, receptor_no := receptor_deg_perm == "NO"]

ComplexUpset::upset(
  cci_deg_comp[
    Dataset == "TMS Droplet (male)"
  ],
  c("scdiffcom_up", "scdiffcom_down", "scdiffcom_no",
  "ligand_up", "ligand_down", "ligand_no",
  "receptor_up", "receptor_down", "receptor_no")
)






## Extract specific examples to discuss ####

## LRI not detected at all and counts of detection for each LRI per dataset ####

seurat_genes <- lapply(
  dataset_path,
  function(i) {
    rownames(readRDS(i))
  }
)
seurat_genes$calico <- unique(
  c(seurat_genes$calico_kidney,
  seurat_genes$calico_lung,
  seurat_genes$calico_spleen
  )
)

LRI_template <- LRI_mouse$LRI_curated[, 1:6]


LRI_in_seurat <- lapply(
  seurat_genes,
  function (i) {
    LRI_mouse$LRI_curated[
      LIGAND_1 %in% i &
        RECEPTOR_1 %in% i &
        LIGAND_2 %in% c(i, NA) &
        RECEPTOR_2 %in% c(i, NA) &
        RECEPTOR_3 %in% c(i, NA)
    ][]
  }
)
LRI_in_seurat <- unique(
  rbindlist(LRI_in_seurat)
)

nrow(LRI_template) - nrow(LRI_in_seurat)
(nrow(LRI_template) - nrow(LRI_in_seurat))/nrow(LRI_template)*100

(nrow(LRI_template) - nrow(LRI_in_seurat)) +
  sum(!LRI_in_seurat$LRI %in% unique(CCI_results$LRI))

((nrow(LRI_template) - nrow(LRI_in_seurat)) +
  sum(!LRI_in_seurat$LRI %in% unique(CCI_results$LRI))) /
  nrow(LRI_template)*100

table(LRI_in_seurat$LRI %in% unique(CCI_results$LRI))
table(LRI_in_seurat$LRI %in% unique(CCI_results$LRI))/
  nrow(LRI_template)*100



LRI_not_detected <- LRI_in_seurat[
  !LRI %in% unique(CCI_results$LRI)
]

LRI_not_detected_genes <- sort(unique(
  unlist(
    LRI_not_detected[
      ,
      c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")
    ]
  )
))
LRI_not_detected_genes <- LRI_not_detected_genes[!is.na(LRI_not_detected_genes)]

# example of LRI not detected that could have been DEG

wilcox_results[
  gene == "Agtr2"
]

LRI_not_detected_deg_up <- rbindlist(
  lapply(
    unique(CCI_results$Dataset),
    function(i) {
      rbindlist(
        lapply(
          unique(CCI_results[Dataset == i]$Tissue),
          function(j) {
            temp <- wilcox_results[
              dataset == i &
                tissue == j
            ]
            temp2 <- copy(LRI_not_detected)[, c("LIGAND_1", "RECEPTOR_1")]
            temp2[
              unique(
                temp[regulation == "UP"][, c("gene", "cell_type")]
              ),
              on = "LIGAND_1==gene",
              emitter := i.cell_type
            ]
            temp2[
              unique(
                temp[regulation == "UP"][, c("gene", "cell_type")]
              ),
              on = "RECEPTOR_1==gene",
              receiver := i.cell_type
            ]
            na.omit(temp2)
          }
        ),
        idcol = "Tissue"
      )
    }
  ),
  idcol = "Dataset"
)

LRI_not_detected_deg_down <- rbindlist(
  lapply(
    unique(CCI_results$Dataset),
    function(i) {
      rbindlist(
        lapply(
          unique(CCI_results[Dataset == i]$Tissue),
          function(j) {
            temp <- CCI_results[
              Dataset == i &
                Tissue == j
            ]
            temp2 <- copy(LRI_not_detected)[, c("LIGAND_1", "RECEPTOR_1")]
            temp2[
              unique(
                temp[l1_deg_reg == "DOWN"][, c("LIGAND_1", "EMITTER_CELLTYPE")]
              ),
              on = "LIGAND_1",
              emitter := i.EMITTER_CELLTYPE
            ]
            temp2[
              unique(
                temp[r1_deg_reg == "DOWN"][, c("RECEPTOR_1", "RECEIVER_CELLTYPE")]
              ),
              on = "RECEPTOR_1",
              receiver := i.RECEIVER_CELLTYPE
            ]
            na.omit(temp2)
          }
        ),
        idcol = "Tissue"
      )
    }
  ),
  idcol = "Dataset"
)

# Stable LRI

LRI_regulation_table <- dcast.data.table(
  CCI_results[, .N, by = c("LRI", "REGULATION")],
  formula = LRI ~ REGULATION,
  value.var = "N",
  fill = 0
)

LRI_regulation_table_notreg <- LRI_regulation_table[
  DOWN == 0 & UP == 0
][order(-FLAT)]
