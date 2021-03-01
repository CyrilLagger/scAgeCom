####################################################
##
## Project: scAgeCom
##
## Last update - December 2020
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## Perform a global analysis (e.g. on all tissues)
## on the scDiffCom results.
##
####################################################
##

## Libraries ####
library(scDiffCom)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

## Read light scDiffComCombined objects ####

DATASETS_COMBINED <- readRDS("../data_scAgeCom/analysis/outputs_data/TMS_scAgeCom_processed_bestORA.rds")

## Read heavy data for subsetting facs to make it comparable to droplet ####

DATASETS_COMBINED_log15_heavy <- readRDS("../data_scAgeCom/analysis/outputs_data/analysis4_DATASETS_COMBINED_log15.rds")
facs_tissue_to_keep <- c(
  intersect(
    unique(DATASETS_COMBINED_log15_heavy$droplet_mixed@cci_raw$ID),
    unique(DATASETS_COMBINED_log15_heavy$facs_mixed@cci_raw$ID)
    ),
  "Aorta", "Heart"
)
facs_subobject <- DATASETS_COMBINED_log15_heavy$facs_mixed
facs_subobject@cci_raw <- facs_subobject@cci_raw[ID %in% facs_tissue_to_keep]
facs_subobject@cci_detected <- list()
facs_subobject@ora_combined_default <- list()
facs_subobject@ora_combined_stringent <- list()
facs_subobject@ora_default <- list()
facs_subobject@ora_stringent <- list()

facs_subobject <- FilterCCI(
  facs_subobject
)


cell_types_dt <- setDT(read.csv(paste0(path_output_directory, "inputs_data/scDiffCom_cell_types.csv"), stringsAsFactors = FALSE))
genage_mouse <- setDT(read.csv(paste0(path_output_directory, "inputs_data/genage_mouse.tsv"), sep = "\t", header = TRUE))

temp_dt <- facs_subobject@cci_detected
temp_dt[cell_types_dt, on = "EMITTER_CELLTYPE==scDiffCom.cell.type", EMITTER_CELL_FAMILY := i.Family...broad]
temp_dt[cell_types_dt, on = "RECEIVER_CELLTYPE==scDiffCom.cell.type", RECEIVER_CELL_FAMILY := i.Family...broad]
temp_dt[, ER_CELL_FAMILY := paste(EMITTER_CELL_FAMILY, RECEIVER_CELL_FAMILY, sep = "_")]
temp_dt[,GENAGE := ifelse(
  LIGAND_1 %in% genage_mouse$Gene.Symbol | LIGAND_2 %in% genage_mouse$Gene.Symbol |
    RECEPTOR_1 %in% genage_mouse$Gene.Symbol | RECEPTOR_2 %in% genage_mouse$Gene.Symbol,
  "YES",
  "NO"
)]
temp_dt[cell_types_dt, on = "EMITTER_CELLTYPE==scDiffCom.cell.type", EMITTER_CELLTYPE_ABR := i.Abbreviation]
temp_dt[cell_types_dt, on = "RECEIVER_CELLTYPE==scDiffCom.cell.type", RECEIVER_CELLTYPE_ABR := i.Abbreviation]

facs_subobject@cci_detected <- temp_dt

facs_subobject <- RunORA(
  object = facs_subobject,
  categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ER_CELL_FAMILY", "GENAGE"),
  overwrite = FALSE,
  stringent_or_default = "default",
  stringent_logfc_threshold = NULL,
  global = FALSE
)

facs_subobject <- RunORA(
  object = facs_subobject,
  categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ID", "ER_CELL_FAMILY", "GENAGE"),
  overwrite = FALSE,
  stringent_or_default = "default",
  stringent_logfc_threshold = NULL,
  global = TRUE
)

facs_subobject@cci_raw <- list()

rm(DATASETS_COMBINED_log15_heavy)

saveRDS(facs_subobject, "../data_scAgeCom/analysis/outputs_data/analysis5_DATASET_FACS_SUBSETTED.rds")

facs_subobject <- readRDS("../data_scAgeCom/analysis/outputs_data/analysis5_DATASET_FACS_SUBSETTED.rds")

DATASETS_COMBINED_log15$facs_mixed_subsetted <- facs_subobject


## Create a function that summarize results by category ####

get_summary_cci <- function(
  object,
  category
) {
  cci_dt_full <- copy(object@cci_detected)
  if (category == "ID") {
    cci_dt_full[, N_CELLTYPES := max(uniqueN(EMITTER_CELLTYPE), uniqueN(RECEIVER_CELLTYPE)), by = "ID"]
    cci_dt <- cci_dt_full[, c("ID", "N_CELLTYPES", "REGULATION")]
    cols_key <- c("ID", "REGULATION")
    ora_dt <- copy(object@ora_combined_default$ID)
    category_name <- category
    col_id <- "N_CCI"
    ora_on <- "ID==VALUE"
    final_col_order <- c(
      "ID", "N_CELLTYPES", "N_CCI", "AVG_LR_PER_ER", "N_CCI_DOWN",
      "N_CCI_FLAT", "N_CCI_UP", "N_CCI_NON_SIGNIFICANT_CHANGE",
      "PCT_CCI_STABLE", "ORA_SCORE_UP", "ORA_SCORE_DOWN",
      "ORA_OR_UP", "ORA_OR_DOWN",  "ORA_P_VALUE_UP" , "ORA_P_VALUE_DOWN"
    )
  } else if (category == "ER_CELLTYPES") {
    cci_dt <- cci_dt_full[, c("ID", "ER_CELLTYPES", "REGULATION")]
    cci_dt[, ID_ER_CELLTYPES := paste(ID, ER_CELLTYPES, sep = "_")]
    cols_key <- c("ID_ER_CELLTYPES", "REGULATION")
    ora_dt <- object@ora_default$ER_CELLTYPES
    category_name <- "ID_ER_CELLTYPES"
    col_id <- "N_LR"
    ora_on <- c("ID==ID", "ER_CELLTYPES==VALUE")
    final_col_order <- c(
      "ER_CELLTYPES", "ID", "N_LR", "N_LR_DOWN", "N_LR_FLAT", "N_LR_UP",
      "N_LR_NON_SIGNIFICANT_CHANGE", "ORA_SCORE_UP", "ORA_SCORE_DOWN",
      "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
    )
  } else if (category == "ER_CELL_FAMILY") {
    cci_dt <- cci_dt_full[, c("ID", "ER_CELL_FAMILY", "REGULATION")]
    cci_dt[, ID_ER_CELL_FAMILY := paste(ID, ER_CELL_FAMILY, sep = "_")]
    cols_key <- c("ID_ER_CELL_FAMILY", "REGULATION")
    ora_dt <- object@ora_default$ER_CELL_FAMILY
    category_name <- "ID_ER_CELL_FAMILY"
    col_id <- "N_LR"
    ora_on <- c("ID==ID", "ER_CELL_FAMILY==VALUE")
    final_col_order <- c(
      "ER_CELL_FAMILY", "ID", "N_LR", "N_LR_DOWN", "N_LR_FLAT", "N_LR_UP",
      "N_LR_NON_SIGNIFICANT_CHANGE", "ORA_SCORE_UP", "ORA_SCORE_DOWN",
      "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
    )
  } else if (category == "ER_CELL_FAMILY_GLOBAL") {
    cci_dt <- cci_dt_full[, c("ER_CELL_FAMILY", "REGULATION")]
    cols_key <- c("ER_CELL_FAMILY", "REGULATION")
    ora_dt <- object@ora_combined_default$ER_CELL_FAMILY
    category_name <- "ER_CELL_FAMILY"
    col_id <- "N_LR"
    ora_on <- "ER_CELL_FAMILY==VALUE"
    final_col_order <- c(
      "ER_CELL_FAMILY", "N_LR", "N_LR_DOWN", "N_LR_FLAT", "N_LR_UP",
      "N_LR_NON_SIGNIFICANT_CHANGE", "ORA_SCORE_UP", "ORA_SCORE_DOWN",
      "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
    )
  } else if (category == "LR_GENES") {
    cci_dt <- cci_dt_full[, c("ID", "LR_GENES", "REGULATION")]
    cci_dt[, ID_LR_GENES := paste(ID, LR_GENES, sep = "_")]
    cols_key <- c("ID_LR_GENES", "REGULATION")
    ora_dt <- object@ora_default$LR_GENES
    category_name <- "ID_LR_GENES"
    col_id <- "N_ER"
    ora_on <- c("ID==ID", "LR_GENES==VALUE")
    final_col_order <- c(
      "LR_GENES", "ID", "N_ER", "N_ER_DOWN", "N_ER_FLAT", "N_ER_UP",
      "N_ER_NON_SIGNIFICANT_CHANGE", "ORA_SCORE_UP", "ORA_SCORE_DOWN",
      "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
    )
  } else if (category %in% c("LR_GENES_GLOBAL", "GO_TERMS_GLOBAL", "KEGG_PWS_GLOBAL")) {
    cci_dt <- cci_dt_full[, c("LR_GENES", "REGULATION")]
    cols_key <- c("LR_GENES", "REGULATION")
    ora_dt<- object@ora_combined_default$LR_GENES
    category_name <- "LR_GENES"
    col_id <- "N_ER"
    ora_on <-  "LR_GENES==VALUE"
    final_col_order <- c(
      "LR_GENES",
      "ORA_SCORE_UP", "ORA_SCORE_DOWN", "N_ID", "N_ID_UP", "N_ID_DOWN",
      "N_ORA_ID_UP", "N_ORA_ID_DOWN",
      "N_ER", "N_ER_DOWN", "N_ER_FLAT", "N_ER_UP", "N_ER_NON_SIGNIFICANT_CHANGE",
      "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
    )
  }
  setkeyv(cci_dt, cols = cols_key)
  res <- cci_dt[CJ(get(cols_key[[1]]), get(cols_key[[2]]), unique = TRUE), .(N_temp = .N), by = .EACHI]
  res <- dcast.data.table(
    data = res,
    formula = get(cols_key[[1]]) ~ get(cols_key[[2]]),
    value.var = "N_temp"
  )
  cci_temp_N <- cci_dt[, .N, by = get(cols_key[[1]])]
  res[cci_temp_N, on = "cols_key==get", (col_id) := i.N]
  #res[, (paste0(col_id, "_YOUNG")) := DOWN + DOWN_DISAPPEARS + FLAT + UP]
  #res[, (paste0(col_id, "_OLD")):= DOWN + FLAT + UP + UP_APPEARS]
  setnames(
    res,
    old = c("cols_key", colnames(res)[2:5]),
    new = c(category_name, paste(col_id, colnames(res)[2:5], sep = "_"))
  )
  if (category == "ID") {
    res[cci_dt, on = "ID", N_CELLTYPES := i.N_CELLTYPES]
    res[, AVG_LR_PER_ER := N_CCI/(N_CELLTYPES^2)]
    res[, PCT_CCI_STABLE := N_CCI_FLAT/N_CCI*100]
  } else if (category == "ER_CELLTYPES") {
    res[cci_dt, on = "ID_ER_CELLTYPES", `:=` (ID = i.ID, ER_CELLTYPES = i.ER_CELLTYPES)]
    res[, ID_ER_CELLTYPES := NULL]
  } else if (category == "ER_CELL_FAMILY") {
    res[cci_dt, on = "ID_ER_CELL_FAMILY", `:=` (ID = i.ID, ER_CELL_FAMILY = i.ER_CELL_FAMILY)]
    res[, ID_ER_CELL_FAMILY := NULL]
  } else if (category == "LR_GENES") {
    res[cci_dt, on = "ID_LR_GENES", `:=` (ID = i.ID, LR_GENES = i.LR_GENES)]
    res[, ID_LR_GENES := NULL]
  } else if (category %in% c("LR_GENES_GLOBAL", "GO_TERMS_GLOBAL", "KEGG_PWS_GLOBAL")) {
    res[cci_dt_full[, uniqueN(ID), by = "LR_GENES"], on = "LR_GENES", N_ID := i.V1]
    res[cci_dt_full[REGULATION == "UP", uniqueN(ID), by = "LR_GENES"], on = "LR_GENES", N_ID_UP := i.V1 ]
    res[cci_dt_full[REGULATION == "DOWN", uniqueN(ID), by = "LR_GENES"], on = "LR_GENES", N_ID_DOWN := i.V1 ]
    ora_dt_2 <- object@ora_default$LR_GENES
    res[ora_dt_2[OR_UP > 1 & BH_P_VALUE_UP <= 0.05, .N, by = c("VALUE")], on = "LR_GENES==VALUE", N_ORA_ID_UP := i.N]
    res[ora_dt_2[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05, .N, by = c("VALUE")], on = "LR_GENES==VALUE", N_ORA_ID_DOWN := i.N]
    setnafill(res, fill = 0, cols = c("N_ID_UP", "N_ID_DOWN", "N_ORA_ID_UP", "N_ORA_ID_DOWN"))
  }

  res[ora_dt, on = ora_on,
      `:=` (ORA_OR_UP = i.OR_UP, ORA_P_VALUE_UP = i.P_VALUE_UP, ORA_SCORE_UP = i.ORA_SCORE_UP,
            ORA_OR_DOWN = i.OR_DOWN, ORA_P_VALUE_DOWN = i.P_VALUE_DOWN, ORA_SCORE_DOWN = i.ORA_SCORE_DOWN) ]
  setcolorder(
    x = res,
    neworder = final_col_order
  )
  if (category %in% c("GO_TERMS_GLOBAL", "KEGG_PWS_GLOBAL")) {
    if (category == "GO_TERMS_GLOBAL") {
      LRdb_temp <- LRdb_mouse$LRdb_curated_GO
      id_name <- "GO_NAME"
      id_category <- "GO_TERMS"
    }
    if (category == "KEGG_PWS_GLOBAL") {
      LRdb_temp <- LRdb_mouse$LRdb_curated_KEGG
      id_name <- "KEGG_NAME"
      id_category <- "KEGG_PWS"
    }
    res <- merge.data.table(
      LRdb_temp,
      res[, c("LR_GENES", "N_ER", "N_ER_DOWN", "N_ER_FLAT", "N_ER_UP",
              "N_ER_NON_SIGNIFICANT_CHANGE")],
      all.x = TRUE,
      by.x = "LR_GENES",
      by.y = "LR_GENES",
      sort = FALSE
    )
    res <- na.omit(res, cols = "LR_GENES")
    res[is.na(res)] <- 0
    res <- res[, lapply(.SD, sum), by = c(id_name), .SDcols = c(
      "N_ER", "N_ER_DOWN", "N_ER_FLAT",
      "N_ER_UP", "N_ER_NON_SIGNIFICANT_CHANGE"  
    )]
    res <- res[N_ER > 0]
    ora_dt_3 <- object@ora_combined_default[[id_category]]
    res[ora_dt_3, on = c(paste0(id_name, "==VALUE")),
        `:=` (ORA_OR_UP = i.OR_UP, ORA_P_VALUE_UP = i.P_VALUE_UP, ORA_SCORE_UP = i.ORA_SCORE_UP,
              ORA_OR_DOWN = i.OR_DOWN, ORA_P_VALUE_DOWN = i.P_VALUE_DOWN, ORA_SCORE_DOWN = i.ORA_SCORE_DOWN) ]
    ora_dt_4 <- object@ora_default[[id_category]]
    res[ora_dt_4[OR_UP > 1 & BH_P_VALUE_UP <= 0.05, .N, by = c("VALUE")], on = paste0(id_name, "==VALUE"), N_ORA_ID_UP := i.N]
    res[ora_dt_4[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05, .N, by = c("VALUE")], on = paste0(id_name, "==VALUE"), N_ORA_ID_DOWN := i.N]
    setnafill(res, fill = 0, cols = c("N_ORA_ID_UP", "N_ORA_ID_DOWN"))
    setcolorder(
      x = res,
      neworder = c(id_name,
                   "ORA_SCORE_UP", "ORA_SCORE_DOWN", "N_ORA_ID_UP", "N_ORA_ID_DOWN",
                   "N_ER", "N_ER_DOWN", "N_ER_FLAT", "N_ER_UP",
                   "N_ER_NON_SIGNIFICANT_CHANGE",
                   "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
      )
    )
  }
  res[, ORA_SCORE_UP := ifelse(ORA_P_VALUE_UP <= 0.05 & ORA_OR_UP > 1, ORA_SCORE_UP, NA)]
  res[, ORA_SCORE_DOWN := ifelse(ORA_P_VALUE_DOWN <= 0.05 & ORA_OR_DOWN > 1, ORA_SCORE_DOWN, NA)]
  return(res)
}

## Apply summary function on all categories ####
categories <- c("ID", "ER_CELLTYPES", "ER_CELL_FAMILY", "ER_CELL_FAMILY_GLOBAL", "LR_GENES",
                "LR_GENES_GLOBAL", "GO_TERMS_GLOBAL", "KEGG_PWS_GLOBAL")


SUMMARY_DATA <- lapply(
  DATASETS_COMBINED,
  function(dataset) {
    sapply(
      categories,
      function(i) {
        print(i)
        get_summary_cci(object = dataset$dataset, category = i)
      },
      USE.NAMES = TRUE,
      simplify = FALSE
    )
  }
)


## Save summary information ####

saveRDS(SUMMARY_DATA, "../data_scAgeCom/analysis/outputs_data/analysis5_SUMMARY_DATA.rds")

test <- readRDS("../data_scAgeCom/analysis/outputs_data/an")

## Function to do Disease Ontology analysis (temporarily here) ####
get_DO_interactions <- function(
  species,
  LR_db
) {
  DO_NAME <- i.DO_name <- NULL
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  if (!requireNamespace("ontoProc", quietly = TRUE)) {
    stop("Package \"ontoProc\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!requireNamespace("ontologyIndex", quietly = TRUE)) {
    stop("Package \"ontologyIndex\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  LR_genes <- unique(unlist(LR_db[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
  LR_genes <- LR_genes[!is.na(LR_genes)]
  DO_mgi_info <- readxl::read_xlsx("../../../../../geneTabReport.xlsx")
  setDT(DO_mgi_info)
  DO_mgi_info <- na.omit(DO_mgi_info, cols = "Mouse Homologs")
  DO_mgi_info[, mgi_symbol := gsub("*", "", `Mouse Homologs`, fixed = TRUE)]
  DO_mgi_info <- DO_mgi_info[mgi_symbol %in% LR_genes]
  onto_do_terms <- ontologyIndex::get_OBO("../../../../../doid.obo")
  do_names <- onto_do_terms$name
  LR_genes_do <- sapply(
    LR_genes,
    function(gene) {
      temp_go <- unique(DO_mgi_info[mgi_symbol == gene]$`Disease Term`)
      ontologyIndex::get_ancestors(onto_do_terms, names(do_names[do_names %in% temp_go]))
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )
  LR_interactions_do_intersection <- rbindlist(
    apply(
      LR_db,
      MARGIN = 1,
      function(row) {
        LIGAND_DO <- unique(c(
          LR_genes_do[[row[["LIGAND_1"]]]],
          LR_genes_do[[row[["LIGAND_2"]]]]
        ))
        RECEPTOR_DO <- unique(c(
          LR_genes_do[[row[["RECEPTOR_1"]]]],
          LR_genes_do[[row[["RECEPTOR_2"]]]],
          LR_genes_do[[row[["RECEPTOR_3"]]]]
        ))
        res_inter <- intersect(LIGAND_DO, RECEPTOR_DO)
        if (length(res_inter) > 0) {
          res_inter <- data.table(
            LR_SORTED = rep(row[["LR_SORTED"]], length(res_inter)),
            DO_ID = res_inter
          )
        } else {
          res_inter <- NULL
        }
        return(res_inter)
      }
    )
  )
  do_id_name_dt <- data.table(
    DO_name = do_names,
    ID = names(do_names)
  )
  LR_interactions_do_intersection[
    do_id_name_dt,
    on = "DO_ID==ID",
    DO_NAME := i.DO_name
    ]
  return(list(
    LR_DO_intersection = LR_interactions_do_intersection
  ))
}

## Study of interesting global LR_GENES signals #####

global_LR_GENES_facs <- SUMMARY_log15$facs_mixed$LR_GENES_GLOBAL
global_LR_GENES_facs_sub <- SUMMARY_log15$facs_mixed_subsetted$LR_GENES_GLOBAL
global_LR_GENES_droplet <- SUMMARY_log15$droplet_mixed$LR_GENES_GLOBAL
global_LR_GENES_calico <- SUMMARY_log15$calico$LR_GENES_GLOBAL

sort(intersect(global_LR_GENES_droplet[ORA_SCORE_UP > 10]$LR_GENES, global_LR_GENES_facs_sub[ORA_SCORE_UP > 10]$LR_GENES))
sort(intersect(global_LR_GENES_droplet[ORA_SCORE_DOWN > 10]$LR_GENES, global_LR_GENES_facs_sub[ORA_SCORE_DOWN > 10]$LR_GENES))

sort(intersect(global_LR_GENES_droplet[ORA_SCORE_UP > 1]$LR_GENES, global_LR_GENES_facs_sub[ORA_SCORE_UP > 1]$LR_GENES))
sort(intersect(global_LR_GENES_droplet[ORA_SCORE_DOWN > 1]$LR_GENES, global_LR_GENES_facs_sub[ORA_SCORE_DOWN > 1]$LR_GENES))

sort(global_LR_GENES_facs[ORA_SCORE_UP > 50 & N_ID_UP > 10 & N_ORA_ID_UP > 5]$LR_GENES)
sort(global_LR_GENES_droplet[ORA_SCORE_UP > 50 & N_ID_UP > 5 & N_ORA_ID_UP > 3]$LR_GENES)

sort(intersect(global_LR_GENES_facs_sub[ORA_SCORE_UP > 0]$LR_GENES, global_LR_GENES_facs[ORA_SCORE_UP > 0]$LR_GENES))

sort(setdiff(global_LR_GENES_facs_sub[ORA_SCORE_UP > 10]$LR_GENES, global_LR_GENES_facs[ORA_SCORE_UP > 10]$LR_GENES))
sort(setdiff(global_LR_GENES_facs[ORA_SCORE_UP > 10]$LR_GENES, global_LR_GENES_facs_sub[ORA_SCORE_UP > 10]$LR_GENES))

#full_cci_table_droplet[LR_NAME == "Lgals1:Itgb1"][, .N, by = L_CELLTYPE][order(-N)][1:20]
#full_cci_table_facs[LR_NAME == "Hspa8:Adrb2"][, .N, by = L_CELLTYPE][order(-N)][1:20]
#heat-shock protein related to aging, what's its role in IC, with adrenergic receptors...
#see https://pubmed.ncbi.nlm.nih.gov/23153586/
#postulate new mechanism? e.g. development, see heart, marrow
#The present data provide new evidence that serum concentration of Hsp70 decreases with age in a normal population.
#Our study also shows that higher levels of Hsp70 are associated with inflammation and frailty in elderly patients.
#full_cci_table_facs[LR_NAME == "Psen1:Notch2"][, .N, by = L_CELLTYPE][order(-N)][1:20]
#full_cci_table_facs[LR_NAME == "Col3a1:Itgb1"][, .N, by = L_CELLTYPE][order(-N)][1:20]
#colagen, integrin stuff ###
# The bulk of the current literature suggests that collagen-binding integrins only have a limited role
#in adult connective tissue homeostasis, partly due to a limited availability of cell-binding sites in the 
#mature fibrillar collagen matrices. However, some recent data suggest that, instead, they are more crucial for
#dynamic connective tissue remodeling events 
#– such as wound healing – where they might act specifically to remodel and restore the tissue architecture. 
#full_cci_table_facs[LR_NAME == "Uba52:Tgfbr2"][, .N, by = L_CELLTYPE][order(-N)][1:20]
#full_cci_table_facs[LR_NAME == "Uba52:Tgfbr2"][, .N, by = R_CELLTYPE][order(-N)][1:20]
#Ubb up but Uba52 and Ubc down!!!!
#full_cci_table_facs[LR_NAME == "Lgals1:Ptprc"][, .N, by = L_CELLTYPE][order(-N)][1:20]
#link with apoptosis ## T cell - T cell death is going down!!!!!
#Although not absolutely required for susceptibility to galectin-1,
#CD45 is a major receptor for galectin-1 on T cells, acts as a negative and positive regulator of galectin-1 death,
#and enhances phagocytic clearance of cells killed by galectin-1 (11,–15).
#full_cci_table_facs[LR_NAME == "App:Ncstn"][, .N, by = L_CELLTYPE][order(-N)][1:20]
#related to AD, but why in so many other tissues! see related article! new role of App!
#APP function and Aβ pathology in adipose tissue
#APP function and Aβ pathology in muscle
#full_cci_table_facs[LR_NAME == "Apoe:Sorl1"][, .N, by = L_CELLTYPE][order(-N)][1:20]
#full_cci_table_facs[LR_NAME == "Ltb:Tnfrsf1a"][, .N, by = L_CELLTYPE][order(-N)][1:20]
#full_cci_table_facs[LR_NAME == "Calm1:Ptpra"][, .N, by = L_CELLTYPE][order(-N)]
#full_cci_table_facs[LR_NAME == "Ubb:Ripk1"][, .N, by = REGULATION_SIMPLE]
#full_cci_table_facs[LR_NAME == "B2m:Hfe"][, .N, by = REGULATION_SIMPLE]
## look at iron stuff
#full_cci_table_facs[LR_NAME == "Calm2:Insr"][, .N, by = REGULATION_SIMPLE]
#full_cci_table_facs[LR_NAME == "Hmgb1:Thbd"][, .N, by = REGULATION_SIMPLE]
#anti-inflammatory signal !!!
#full_cci_table_facs[LR_NAME == "Ubb:Ripk1", c("Tissue", "LR_CELLTYPE", "LOGFC", "REGULATION")][order(-LOGFC)]
#summary(full_cci_table_facs[LR_NAME == "Ubb:Ripk1"]$LOGFC)


## Tissue comparison ####

test_UP <- test[REGULATION_SIMPLE == "UP"]

test_2 <- sapply(
  unique(test_UP$ID),
  function(tiss) {
    unique(test_UP[ID == tiss]$LR_GENES)
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

test_3 <- sapply(
  seq_along(test_2),
  function(i) {
    sapply(
      seq_along(test_2),
      function(j) {
        if (i == j) {
          0
        } else {
          length(intersect(test_2[[i]], test_2[[j]]))/length(unique(c(test_2[[i]], test_2[[j]])))
        }
      }
    )
  }
)
colnames(test_3) <- names(test_2)
rownames(test_3) <- names(test_2)


hist(test_UP$LOGFC, breaks = 50)


## Clustering analysis ####

cci_dt <- DATASETS_COMBINED_log15$facs_mixed@cci_detected
cci_dt[, ID_ER := paste(ID, ER_CELLTYPES, sep = "_")]
cci_dt[, PI_VALUE := -log(P_VALUE_DE+1E-4)*LOGFC]

cci_pi_matrix_ER <- as.matrix(
  dcast.data.table(
    cci_dt[, c("ID_ER", "LR_GENES", "PI_VALUE") ],
    formula = ID_ER ~ LR_GENES,
    value.var = "PI_VALUE"
  ),
  rownames = "ID_ER"
)
cci_pi_matrix_ER[is.na(cci_pi_matrix_ER)] <- 0

cci_pi_matrix_LR <- t(cci_pi_matrix_ER)

#scaling
cci_pi_matrix_ER_scaled <- scale(cci_pi_matrix_ER)
cci_pi_matrix_LR_scaled <- scale(cci_pi_matrix_LR)


#try sparcl
library(sparcl)
sparcl_LR <- KMeansSparseCluster(
  x = cci_pi_matrix_LR,
  K = 10,
  wbounds = 1.5
)
sparcl_LR_cl <- sparcl_LR[[1]]$Cs
table(sparcl_LR_cl)
sparcl_LR_cl[sparcl_LR_cl == 10]


sparcl_ER <- KMeansSparseCluster(
  x = cci_pi_matrix_ER,
  K = 10,
  wbounds = 1.5
)
sparcl_ER_cl <- sparcl_ER[[1]]$Cs
table(sparcl_ER_cl)

Heatmap(
  matrix = cci_pi_matrix_LR,
  #col = hm_col,
  name = "test",
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_split = sparcl_ER_cl,
  row_split = sparcl_LR_cl,
  use_raster = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE
)

library(skmeans)
skm_LR <- skmeans(cci_pi_matrix_LR, k = 20)
skm_ER <- skmeans(cci_pi_matrix_ER, k = 20)
table(skm_LR$cluster)
table(skm_ER$cluster)

km_LR <- kmeans(cci_pi_matrix_LR,  10)
km_ER <- kmeans(cci_pi_matrix_ER, 10)
table(km_LR$cluster)
table(km_ER$cluster)

Heatmap(
  matrix = cci_pi_matrix_LR,
  #row_km = 20,
  #column_km = 10,
  #col = hm_col,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = skm_ER$cluster,
  row_split = skm_LR$cluster,
  use_raster = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE
)

#use binary data
#cci_pi_vals <- as.vector(cci_pi_matrix_LR)
#cci_pi_vals <- cci_pi_vals[cci_pi_vals != 0]
#summary(cci_pi_vals)

-log(0.05)*log(1.5)

cci_pi_matrix_binary_LR <- cci_pi_matrix_LR
cci_pi_matrix_binary_LR[abs(cci_pi_matrix_binary_LR) < -log(0.05)*log(1.5)] <- 0
cci_pi_matrix_binary_LR[cci_pi_matrix_binary_LR >= -log(0.05)*log(1.5)] <- 1
cci_pi_matrix_binary_LR[cci_pi_matrix_binary_LR <= log(0.05)*log(1.5)] <- -1
unique(as.vector(cci_pi_matrix_binary_LR))

table(rowSums(cci_pi_matrix_binary_LR))
table(colSums(cci_pi_matrix_binary_LR))
which(rowSums(cci_pi_matrix_binary_LR) == -124)

cci_pi_matrix_binary_LR <- cci_pi_matrix_binary_LR[abs(rowSums(cci_pi_matrix_binary_LR)) > 2, abs(colSums(cci_pi_matrix_binary_LR)) > 1]
cci_pi_matrix_binary_LR <- cci_pi_matrix_binary_LR[abs(rowSums(cci_pi_matrix_binary_LR)) > 0, abs(colSums(cci_pi_matrix_binary_LR)) > 0]

skm_bin_LR <- skmeans(cci_pi_matrix_binary_LR, k = 25)
skm_bin_ER <- skmeans(t(cci_pi_matrix_binary_LR), k = 25)
table(skm_bin_LR$cluster)
table(skm_bin_ER$cluster)

Heatmap(
  matrix = cci_pi_matrix_binary_LR,
  #row_km = 15,
  #column_km = 10,
  #col = hm_col,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = skm_bin_ER$cluster,
  row_split = skm_bin_LR$cluster,
  use_raster = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  row_names_max_width = unit(6, "cm"),
  row_names_centered = FALSE,
  row_names_gp = grid::gpar(fontsize = 4),
  column_names_gp = grid::gpar(fontsize = 4)
)

hist(binary_groups_sd, breaks = 70)

skm_bin_ER$cluster[skm_bin_ER$cluster == 9]
skm_bin_LR$cluster[skm_bin_LR$cluster == 6]

binary_groups_means <- t(sapply(
  sort(unique(skm_bin_LR$cluster)),
  function(i) {
    sapply(
      sort(unique(skm_bin_ER$cluster)),
      function(j) {
        rows_to_keep <- names(skm_bin_LR$cluster[skm_bin_LR$cluster == i])
        cols_to_keep <- names(skm_bin_ER$cluster[skm_bin_ER$cluster == j])
        mean(cci_pi_matrix_binary_LR[rows_to_keep, cols_to_keep])
      }
    )
  }
))

binary_groups_abs_means <- t(sapply(
  sort(unique(skm_bin_LR$cluster)),
  function(i) {
    sapply(
      sort(unique(skm_bin_ER$cluster)),
      function(j) {
        rows_to_keep <- names(skm_bin_LR$cluster[skm_bin_LR$cluster == i])
        cols_to_keep <- names(skm_bin_ER$cluster[skm_bin_ER$cluster == j])
        mean(abs(cci_pi_matrix_binary_LR[rows_to_keep, cols_to_keep]))
      }
    )
  }
))

binary_groups_sd <- t(sapply(
  sort(unique(skm_bin_LR$cluster)),
  function(i) {
    sapply(
      sort(unique(skm_bin_ER$cluster)),
      function(j) {
        rows_to_keep <- names(skm_bin_LR$cluster[skm_bin_LR$cluster == i])
        cols_to_keep <- names(skm_bin_ER$cluster[skm_bin_ER$cluster == j])
        sd(cci_pi_matrix_binary_LR[rows_to_keep, cols_to_keep])
      }
    )
  }
))
binary_groups_dt <- data.table(
  avg = as.vector(binary_groups_means),
  sd = as.vector(binary_groups_sd),
  avg_abs = as.vector(binary_groups_abs_means)
)

ggplot(binary_groups_dt, aes(avg, sd)) + geom_point()
ggplot(binary_groups_dt, aes(avg_abs, sd)) + geom_point()
ggplot(binary_groups_dt, aes(avg_abs, avg)) + geom_point()

bin_sd_sorted <- sort(as.vector(binary_groups_sd), decreasing = TRUE)

rank(binary_groups_sd)

which(binary_groups_sd == bin_sd_sorted[[1]], arr.ind = TRUE)

order(c(-1, 1, 4,7, 5))

hist(binary_groups_sd, breaks = 100)

skm_bin_LR$cluster[skm_bin_LR$cluster == 5]
skm_bin_ER$cluster[skm_bin_ER$cluster == 19]
test <- cci_pi_matrix_binary_LR[names(skm_bin_LR$cluster[skm_bin_LR$cluster == 5]), names(skm_bin_ER$cluster[skm_bin_ER$cluster == 19])]

Heatmap(
  test
)


######
library(biclust)
bi_clust <- biclust(cci_pi_matrix_binary_LR, method=BCCC())
plot(bi_clust)
bi_clust
drawHeatmap(cci_pi_matrix_binary_LR, bicResult = bi_clust)
library(fabia)

res <- fabia(cci_pi_matrix_binary_LR, 5, 0.01, 500)
show(res)
plot(res)
extractPlot(res,ti="FABIA")

rb <- extractBic(res)

library(isa2)

isa.result <- isa(cci_pi_matrix_binary_LR)

ncol(isa.result$rows)
layout(cbind(1:5))
sapply(1:5, function(x) {
  par(mar = c(3, 5, 1, 1))
  plot(isa.result$rows[, x], ylim = c(-1, 1), ylab = "scores",
       xlab = NA, cex.lab = 2, cex.axis = 1.5)
})

mymodules <- lapply(seq(ncol(isa.result$rows)), function(x) {
  list(rows = which(isa.result$rows[, x] != 0),
       columns = which(isa.result$columns[, x] !=
                         0))
})

sort(mymodules[[1]]$rows)

my_rows <- lapply(mymodules, function(i) {i[["rows"]]})

sort(table(unlist(my_rows)))

sub_mat <- cci_pi_matrix_binary_LR[mymodules[[30]]$rows, mymodules[[30]]$columns]
Heatmap(sub_mat)

Bc <- isa.biclust(isa.result)

drawHeatmap(cci_pi_matrix_binary_LR, Bc, 1, local = FALSE)

data <- isa.in.silico(noise=0.1)
isa.result <- isa(data[[1]])
## Find the best bicluster for each block in the input
best <- apply(cor(isa.result$rows, data[[2]]), 2, which.max)
## Check correlation
sapply(seq_along(best),
       function(x) cor(isa.result$rows[,best[x]], data[[2]][,x]))
## The same for the columns
sapply(seq_along(best),
       function(x) cor(isa.result$columns[,best[x]], data[[3]][,x]))
## Plot the data and the modules found
if (interactive()) {
  layout(rbind(1:2,3:4))
  image(data[[1]], main="In-silico data")
  sapply(best, function(b) image(outer(isa.result$rows[,b],
                                       isa.result$columns[,b]),
                                 main=paste("Module", b)))
}

plotModules(res_isa, cci_pi_matrix_binary_LR)

image(cci_pi_matrix_binary_LR)

heatmap(cci_pi_matrix_binary_LR)

plotclust(bi_clust, cci_pi_matrix_binary_LR)






mean(test)
sd(as.vector(test))

# 5,19 / 3,6 / 16, 17 / 1,15

binary_groups_means[5, 19]

hist(binary_groups_means, breaks = 50)


#use skmean and then another method inside interesting blocks


#try pca, knn, louvain etc

colMeans(cci_pi_matrix_LR_scaled)

cci_pca <- prcomp(cci_pi_matrix_LR_scaled)
cci_pca
library(factoextra)
fviz_eig(cci_pca, ncp = 20)
library("mstknnclust")

library(PCAtools)

cci_pca_LR <- pca(t(cci_pi_matrix_LR))
screeplot(cci_pca_LR)
biplot(cci_pca_LR)
findElbowPoint(cci_pca_LR$variance)
cci_direct_dist_LR <- dist(cci_pi_matrix_LR)
cci_pca_dist_LR <- dist(cci_pca_LR$rotated[, 1:3])
cci_clust_LR <- mst.knn(as.matrix(cci_pca_dist_LR))$cluster
cci_pca_h_clust_tree_LR <- hclust(cci_pca_dist_LR)
plot(cci_h_clust_tree_LR, labels = FALSE)
cci_pca_h_clust_LR <- cutree(cci_pca_h_clust_tree_LR, k = 100)
table(cci_pca_h_clust_LR)
cci_pca_h_clust_LR[cci_pca_h_clust_LR == 2]

cci_LR_means <- rowMeans(cci_pi_matrix_LR)
cci_LR_sd <- apply(
  cci_pi_matrix_LR,
  MARGIN = 1,
  sd
)
cci_LR_max <- apply(
  abs(cci_pi_matrix_LR),
  MARGIN = 1,
  max
)
cci_LR_dt <- data.table(
  id = names(cci_LR_means),
  avg = cci_LR_means,
  sd = cci_LR_sd,
  max = cci_LR_max
)


ggplot(cci_LR_dt, aes(x = avg, y = sd)) + geom_point() + geom_density2d()
ggplot(cci_LR_dt, aes(log10(abs(avg)), log10(sd))) + geom_point()+ geom_density2d()
ggplot(cci_LR_dt, aes(avg)) + geom_histogram(bins = 100)
ggplot(cci_LR_dt, aes(sd)) + geom_histogram(bins = 100)

LR_keep <- cci_LR_dt[abs(avg) > 0.02 | sd > 1]$id

cci_ER_means <- rowMeans(cci_pi_matrix_ER)
cci_ER_sd <- apply(
  cci_pi_matrix_ER,
  MARGIN = 1,
  sd
)
cci_ER_max <- apply(
  abs(cci_pi_matrix_ER),
  MARGIN = 1,
  max
)
cci_ER_dt <- data.table(
  id = names(cci_ER_means),
  avg = cci_ER_means,
  sd = cci_ER_sd,
  max = cci_ER_max
)
ggplot(cci_ER_dt, aes(x = avg, y = sd)) + geom_point() + geom_density2d()
ggplot(cci_ER_dt, aes(log10(abs(avg)), log10(sd))) + geom_point()+ geom_density2d()
ggplot(cci_ER_dt, aes(avg)) + geom_histogram(bins = 100)
ggplot(cci_ER_dt, aes(sd)) + geom_histogram(bins = 100)


ER_keep <- cci_ER_dt[abs(avg) > 0.05 | sd > 0.5]$id

cci_pi_matrix_LR_sub <- cci_pi_matrix_LR[LR_keep, ER_keep]

ggplot(cci_LR_dt, aes(x = sd, y = max)) + geom_point() + geom_density2d()

plot(log10(abs(cci_LR_means)), log10(cci_LR_sd))


cci_direct_h_clust_tree_LR <- hclust(cci_direct_dist_LR)
plot(cci_direct_h_clust_tree_LR, labels = FALSE)
cci_direct_h_clust_LR <- cutree(cci_direct_h_clust_tree_LR, k = 10)
table(cci_direct_h_clust_LR)
cci_direct_h_clust_LR[cci_direct_h_clust_LR == 2]

cci_pca_ER <- pca((cci_pi_matrix_LR))
screeplot(cci_pca_ER)
biplot(cci_pca_ER)
findElbowPoint(cci_pca_ER$variance)
cci_pca_dist_ER <- dist(cci_pca_ER$rotated[, 1:3])
cci_clust_ER <- mst.knn(as.matrix(cci_pca_dist_ER))$cluster
cci_h_clust_ER <- hclust(cci_pca_dist_ER)
plot(cci_h_clust_ER, labels = FALSE)

hm_col_val <- abs(quantile(as.vector(cci_pi_matrix_LR), 0.01))
hm_col <- colorRamp2(c(-hm_col_val, 0, hm_col_val), c("blue", "green", "red"))

Heatmap(
  matrix = cci_pi_matrix_LR_sub,
  #col = hm_col,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE,
  #column_split = cci_clust_ER,
  #row_split = cci_clust_LR,
  use_raster = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE
)

cci_dt_dc <- dcast.data.table(
  cci_dt[, c("ID_ER", "LR_GENES", "PI_VALUE") ],
  formula = ID_ER ~ LR_GENES,
  value.var = "PI_VALUE"
)
cci_dt_dc[is.na(cci_dt_dc)] <- 0
cci_dt_dc[cci_dt, on = "ID_ER", ID := i.ID]
cci_dt_dc$ID
cci_dt_dc[cci_dt, on = "ID_ER", CF := i.ER_CELL_FAMILY]

cci_pi_Tissue <- rowsum(t(cci_pi_matrix_LR), group = cci_dt_dc$ID)/as.vector(table(cci_dt_dc$ID))
cci_pi_Tissue_sub <- cci_pi_Tissue[, LR_keep]

cci_pi_CF <- rowsum(t(cci_pi_matrix_LR), group = cci_dt_dc$CF)/as.vector(table(cci_dt_dc$CF))
cci_pi_CF_sub <- cci_pi_CF[, LR_keep]


hm_col <- colorRamp2(c(-4, 0, 4), c("blue", "gray", "red"))
Heatmap(
  matrix = t(cci_pi_Tissue_sub),
  #col = hm_col,
  name = "test",
  show_row_names = FALSE,
  show_column_names = TRUE,
  #column_split = cci_clust_ER,
  #row_split = cci_clust_LR,
  use_raster = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE
)

# the matrix are very sparse, kmeans does not work well on it...

sum(as.vector(cci_pi_matrix_LR == 0))/length(as.vector(cci_pi_matrix_LR))




# Heatmap


Heatmap(
  matrix = cci_pi_matrix_LR,
  col = hm_col,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE
)



table(rs$cluster)
table(cs$cluster)

cs_cl <- cs$cluster
rs_cl <- rs$cluster

test <- sapply(
  unique(rs_cl),
  function(i) {
    sapply(
      unique(cs_cl),
      function(j) {
        rows_to_keep <- names(rs_cl[rs_cl == i])
        cols_to_keep <- names(cs_cl[cs_cl == j])
        mean(cci_pi_matrix_LR[rows_to_keep, cols_to_keep])
      }
    )
  }
)

which(test == max(test), arr.ind = TRUE)

rs$cluster[rs$cluster == 16]
cs$cluster[cs$cluster == 2]






#distance
library(factoextra)
library(pheatmap)
ER_dist <- get_dist(facs_pi_matrix_ER_scaled, stand = TRUE, method = "euclidean")
fviz_dist(ER_dist, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

LR_dist <- get_dist(facs_pi_matrix_LR_scaled, stand = TRUE, method = "euclidean")
fviz_dist(LR_dist, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

fviz_nbclust(facs_pi_matrix_ER_scaled, cluster::pam, diss = ,  method = "gap_stat", k.max = 35)

pheatmap(facs_pi_matrix_ER_scaled, scale = "none", show_rownames = FALSE, show_colnames = FALSE, kmeans_k = 10)
pheatmap(facs_pi_matrix_LR_scaled, scale = "none", show_rownames = FALSE, show_colnames = FALSE, kmeans_k = 10)

cl <- kmeans(facs_pi_matrix_ER, centers = 20)$cluster
rw <- kmeans(t(facs_pi_matrix_ER), centers = 20)$cluster

table(cl)
table(rw)

names(cl[cl == 13])

sub_pi_mat <- facs_pi_matrix_ER[! (rownames(facs_pi_matrix_ER) %in% names(cl[cl == 8])),]
sub_pi_mat <- sub_pi_mat[, !(colnames(sub_pi_mat) %in% names(rw[rw == 8]))]

cl2 <- cl[cl != 8]
rw2 <- rw[rw != 8]

max(sub_pi_mat)
min(sub_pi_mat)

hist(as.vector(sub_pi_mat), breaks = 100)

Heatmap(
  sub_pi_mat,
  row_split = cl2,
  column_split = rw2,
  show_row_names = FALSE,
  show_column_names = FALSE
)

cl[cl == 13]
rw[rw == 13]

Heatmap(
  facs_pi_matrix_ER,
  row_km = 15,
  column_km = 15,
  show_row_names = FALSE,
  show_column_names = FALSE
)

LR_km <- kmeans(facs_pi_matrix_LR_scaled, 10)
table(LR_km$cluster)
LR_km2 <- kmeans(facs_pi_matrix_LR, 15)
table(LR_km2$cluster)


ER_km <- kmeans(facs_pi_matrix_ER_scaled, 10)
table(ER_km$cluster)

ER_km$cluster[ER_km$cluster == 1]
ER_km$cluster[ER_km$cluster == 2]
ER_km$cluster[ER_km$cluster == 9]
ER_km$cluster[ER_km$cluster == 10]

LR_km2$cluster[LR_km$cluster == 3]

test <- dcast.data.table(
  facs_cci[ID == "Kidney", c("ID_ER", "LR_GENES", "PI_VALUE") ],
  formula = ID_ER ~ LR_GENES,
  value.var = "PI_VALUE"
)
rownames(test) <- test$ID_ER
test[, ID_ER := NULL]

test2 <- CA(test, ncp = 5, graph = TRUE)

db <- fpc::dbscan(facs_pi_matrix_LR_scaled, eps = 0.15, MinPts = 5)

fviz_cluster(db, data = facs_pi_matrix_LR_scaled, stand = FALSE,
             ellipse = FALSE, show.clust.cent = FALSE,
             geom = "point",palette = "jco", ggtheme = theme_classic())

pca_1 <- prcomp(t(test_mat), scale. = TRUE)
library(factoextra)
fviz_eig(pca_1)
fviz_pca_ind(pca_1,
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             label = "none"
)



hist(colSums(abs(test_mat) > 0), breaks = 10)
head(table(colSums(abs(test_mat) > 0)))
plot(table(colSums(abs(test_mat) > 0)))

sort(colSums(abs(test_mat)), decreasing = TRUE)

hist(colSums(abs(test_mat)), breaks = 50)

test_mat2 <- test_mat[, colSums(abs(test_mat)) > 50]

colSums(abs(test_mat) > 0)["App:Notch1"]
test_mat2[, "App:Notch1"]

Heatmap(test_mat2)
pheatmap(test_mat2, show_rownames = FALSE, show_colnames = FALSE, scale = "column")
library(pheatmap)
library(ComplexHeatmap)

pheatmap(test_3, scale = "none")







## Do some heatmaps??? #####

testxx <- as.matrix(dcast.data.table(
  testx,
  formula = LR_NAME ~ LR_CELL_FAMILY,
  value.var = "frac"
),
rownames = "LR_NAME"
)
testxx[is.na(testxx)] <- 0

heatmap(testxx)

library(circlize)
ComplexHeatmap::Heatmap(
  testxx,
  #clustering_distance_rows = "pearson",
  show_column_names = FALSE,
  show_row_names = FALSE
)

ComplexHeatmap::Heatmap(
  testxx,
  row_km = 10,
  column_km = 10,
  show_column_names = FALSE,
  show_row_names = FALSE
)

cl = kmeans(testxx, centers = 10)
cr = kmeans(t(testxx), centers = 10)

cl_list <- lapply(
  unique(cl$cluster),
  function(i) {
    names(cl$cluster[cl$cluster == i])
  }
)

cr_list <- lapply(
  unique(cr$cluster),
  function(i) {
    names(cr$cluster[cr$cluster == i])
  }
)

cl$cluster

library(seriation)
o = seriate(testxx, method = "BEA_TSP")
Heatmap(testxx, name = "mat", 
        row_order = get_order(o, 1), column_order = get_order(o, 2),
        column_title = "seriation by BEA_TSP method")


pheatmap::pheatmap(testxx, scale = "row")
pheatmap::pheatmap(testxx, scale = "column")

test4 <- test3[, list(text = paste(LR_CELL_FAMILY, collapse = "_")), by = "LR_NAME"]
table(test4$text)



ggplot(test3, aes(x= LR_CELL_FAMILY, y = LR_NAME)) + geom_point(aes(size = N))

facs_LR_GENES_summary



