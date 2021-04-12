####################################################
##
## Project: scAgeCom
##
## Last update - March 2021
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

DATASETS_COMBINED <- readRDS("../data_scAgeCom/analysis/outputs_data/scAgeCom_results_processed.rds")

## Create a function that summarize results by category ####

get_summary_cci <- function(
  object,
  category
) {
  cci_dt_full <- copy(object@cci_table_detected)
  if (category == "ID") {
    cci_dt_full[, N_CELLTYPES := max(uniqueN(EMITTER_CELLTYPE), uniqueN(RECEIVER_CELLTYPE)), by = "ID"]
    cci_dt <- cci_dt_full[, c("ID", "N_CELLTYPES", "REGULATION")]
    cols_key <- c("ID", "REGULATION")
    #ora_dt <- copy(object@ora_combined_default$ID)
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
    ora_dt <- object@ora_table$ER_CELLTYPES
    category_name <- "ID_ER_CELLTYPES"
    col_id <- "N_LR"
    ora_on <- c("ID==ID", "ER_CELLTYPES==VALUE")
    final_col_order <- c(
      "ER_CELLTYPES", "ID", "N_LR", "N_LR_DOWN", "N_LR_FLAT", "N_LR_UP",
      "N_LR_NON_SIGNIFICANT_CHANGE", "ORA_SCORE_UP", "ORA_SCORE_DOWN",
      "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
    )
  } else if (category == "ER_CELLFAMILIES") {
    cci_dt <- cci_dt_full[, c("ID", "ER_CELLFAMILIES", "REGULATION")]
    cci_dt[, ID_ER_CELLFAMILIES := paste(ID, ER_CELLFAMILIES, sep = "_")]
    cols_key <- c("ID_ER_CELLFAMILIES", "REGULATION")
    ora_dt <- object@ora_table$ER_CELLFAMILIES
    category_name <- "ID_ER_CELLFAMILIES"
    col_id <- "N_LR"
    ora_on <- c("ID==ID", "ER_CELLFAMILIES==VALUE")
    final_col_order <- c(
      "ER_CELLFAMILIES", "ID", "N_LR", "N_LR_DOWN", "N_LR_FLAT", "N_LR_UP",
      "N_LR_NON_SIGNIFICANT_CHANGE", "ORA_SCORE_UP", "ORA_SCORE_DOWN",
      "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
    )
  } else if (category == "ER_CELLFAMILIES_GLOBAL") {
    cci_dt <- cci_dt_full[, c("ER_CELLFAMILIES", "REGULATION")]
    cols_key <- c("ER_CELLFAMILIES", "REGULATION")
    #ora_dt <- object@ora_combined_default$ER_CELLFAMILIES
    category_name <- "ER_CELLFAMILIES"
    col_id <- "N_LR"
    ora_on <- "ER_CELLFAMILIES==VALUE"
    final_col_order <- c(
      "ER_CELLFAMILIES", "N_LR", "N_LR_DOWN", "N_LR_FLAT", "N_LR_UP",
      "N_LR_NON_SIGNIFICANT_CHANGE", "ORA_SCORE_UP", "ORA_SCORE_DOWN",
      "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
    )
  } else if (category == "LRI") {
    cci_dt <- cci_dt_full[, c("ID", "LRI", "REGULATION")]
    cci_dt[, ID_LRI := paste(ID, LRI, sep = "_")]
    cols_key <- c("ID_LRI", "REGULATION")
    ora_dt <- object@ora_table$LRI
    category_name <- "ID_LRI"
    col_id <- "N_ER"
    ora_on <- c("ID==ID", "LRI==VALUE")
    final_col_order <- c(
      "LRI", "ID", "N_ER", "N_ER_DOWN", "N_ER_FLAT", "N_ER_UP",
      "N_ER_NON_SIGNIFICANT_CHANGE", "ORA_SCORE_UP", "ORA_SCORE_DOWN",
      "ORA_OR_UP", "ORA_OR_DOWN", "ORA_P_VALUE_UP", "ORA_P_VALUE_DOWN"
    )
  } else if (category %in% c("LRI_GLOBAL", "GO_TERMS_GLOBAL", "KEGG_PWS_GLOBAL")) {
    cci_dt <- cci_dt_full[, c("LRI", "REGULATION")]
    cols_key <- c("LRI", "REGULATION")
    #ora_dt<- object@ora_combined_default$LRI
    category_name <- "LRI"
    col_id <- "N_ER"
    ora_on <-  "LRI==VALUE"
    final_col_order <- c(
      "LRI",
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
  } else if (category == "ER_CELLFAMILIES") {
    res[cci_dt, on = "ID_ER_CELLFAMILIES", `:=` (ID = i.ID, ER_CELLFAMILIES = i.ER_CELLFAMILIES)]
    res[, ID_ER_CELLFAMILIES := NULL]
  } else if (category == "LRI") {
    res[cci_dt, on = "ID_LRI", `:=` (ID = i.ID, LRI = i.LRI)]
    res[, ID_LRI := NULL]
  } else if (category %in% c("LRI_GLOBAL", "GO_TERMS_GLOBAL", "KEGG_PWS_GLOBAL")) {
    res[cci_dt_full[, uniqueN(ID), by = "LRI"], on = "LRI", N_ID := i.V1]
    res[cci_dt_full[REGULATION == "UP", uniqueN(ID), by = "LRI"], on = "LRI", N_ID_UP := i.V1 ]
    res[cci_dt_full[REGULATION == "DOWN", uniqueN(ID), by = "LRI"], on = "LRI", N_ID_DOWN := i.V1 ]
    ora_dt_2 <- object@ora_table$LRI
    res[ora_dt_2[OR_UP > 1 & BH_P_VALUE_UP <= 0.05, .N, by = c("VALUE")], on = "LRI==VALUE", N_ORA_ID_UP := i.N]
    res[ora_dt_2[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05, .N, by = c("VALUE")], on = "LRI==VALUE", N_ORA_ID_DOWN := i.N]
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
      LRI_temp <- LRI_mouse$LRI_curated_GO
      id_name <- "GO_NAME"
      id_category <- "GO_TERMS"
    }
    if (category == "KEGG_PWS_GLOBAL") {
      LRI_temp <- LRI_mouse$LRI_curated_KEGG
      id_name <- "KEGG_NAME"
      id_category <- "KEGG_PWS"
    }
    res <- merge.data.table(
      LRI_temp,
      res[, c("LRI", "N_ER", "N_ER_DOWN", "N_ER_FLAT", "N_ER_UP",
              "N_ER_NON_SIGNIFICANT_CHANGE")],
      all.x = TRUE,
      by.x = "LRI",
      by.y = "LRI",
      sort = FALSE
    )
    res <- na.omit(res, cols = "LRI")
    res[is.na(res)] <- 0
    res <- res[, lapply(.SD, sum), by = c(id_name), .SDcols = c(
      "N_ER", "N_ER_DOWN", "N_ER_FLAT",
      "N_ER_UP", "N_ER_NON_SIGNIFICANT_CHANGE"  
    )]
    res <- res[N_ER > 0]
    #ora_dt_3 <- object@ora_combined_default[[id_category]]
    res[ora_dt_3, on = c(paste0(id_name, "==VALUE")),
        `:=` (ORA_OR_UP = i.OR_UP, ORA_P_VALUE_UP = i.P_VALUE_UP, ORA_SCORE_UP = i.ORA_SCORE_UP,
              ORA_OR_DOWN = i.OR_DOWN, ORA_P_VALUE_DOWN = i.P_VALUE_DOWN, ORA_SCORE_DOWN = i.ORA_SCORE_DOWN) ]
    ora_dt_4 <- object@ora_table[[id_category]]
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
categories <- c("ID", "ER_CELLTYPES", "ER_CELLFAMILIES", "ER_CELLFAMILIES_GLOBAL", "LRI",
                "LRI_GLOBAL", "GO_TERMS_GLOBAL", "KEGG_PWS_GLOBAL")

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

saveRDS(SUMMARY_DATA, "../data_scAgeCom/analysis/outputs_data/scAgeCom_summary.rds")

test <- readRDS("../data_scAgeCom/analysis/outputs_data/analysis5_SUMMARY_DATA.rds")

## Create summary information for each LRI ####

all_cci <- rbindlist(
  lapply(
  DATASETS_COMBINED[c(1,2,3,6,7)],
  function(i) {
    i$dataset@cci_table_detected[, c("ID", "ER_CELLTYPES", "LRI", "ER_CELL_FAMILY", "REGULATION")]
  }
  ),
  use.names = TRUE,
  idcol = "DATASET"
)

all_ORA_LRI  <- rbindlist(
  lapply(
    DATASETS_COMBINED[c(1,2,3,6,7)],
    function(i) {
      i$dataset@ora_table$LRI
    }
  ),
  use.names = TRUE,
  idcol = "DATASET"
)

all_cci[LRI == unique(all_cci$LRI)[[1]]]




## Upset of intersection ####

SUMMARY_DATA <- readRDS("../data_scAgeCom/analysis/outputs_data/analysis5_SUMMARY_DATA.rds")

test <- dcast.data.table(
 DATASETS_COMBINED$droplet_male$dataset@ora_table$GO_TERMS[OR_UP >=1 & BH_P_VALUE_UP <= 0.05, c("ID", "VALUE")],
  VALUE ~ ID,
  fun.aggregate = length
)
colnames(test) <- paste0(colnames(test), "_male")


test2 <- dcast.data.table(
  DATASETS_COMBINED$droplet_female$dataset@ora_table$GO_TERMS[OR_UP >=1 & BH_P_VALUE_UP <= 0.05, c("ID", "VALUE")],
  VALUE ~ ID,
  fun.aggregate = length
)
colnames(test2) <- paste0(colnames(test2), "_female")

test3 <- merge.data.table(
  test,
  test2,
  all = TRUE,
  by.x = "VALUE_male",
  by.y = "VALUE_female"
)
test3[is.na(test3)] <- 0


test2 <- test[Lung == 1 & Marrow == 1]

ComplexUpset::upset(
  data = as.data.frame(test3),
  intersect = colnames(test3)[-1],
  min_size = 5
)


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
  LRI <- unique(unlist(LR_db[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
  LRI <- LRI[!is.na(LRI)]
  DO_mgi_info <- readxl::read_xlsx("../../../../../geneTabReport.xlsx")
  setDT(DO_mgi_info)
  DO_mgi_info <- na.omit(DO_mgi_info, cols = "Mouse Homologs")
  DO_mgi_info[, mgi_symbol := gsub("*", "", `Mouse Homologs`, fixed = TRUE)]
  DO_mgi_info <- DO_mgi_info[mgi_symbol %in% LRI]
  onto_do_terms <- ontologyIndex::get_OBO("../../../../../doid.obo")
  do_names <- onto_do_terms$name
  LRI_do <- sapply(
    LRI,
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
          LRI_do[[row[["LIGAND_1"]]]],
          LRI_do[[row[["LIGAND_2"]]]]
        ))
        RECEPTOR_DO <- unique(c(
          LRI_do[[row[["RECEPTOR_1"]]]],
          LRI_do[[row[["RECEPTOR_2"]]]],
          LRI_do[[row[["RECEPTOR_3"]]]]
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


## Tissue comparison ####

test_UP <- test[REGULATION_SIMPLE == "UP"]

test_2 <- sapply(
  unique(test_UP$ID),
  function(tiss) {
    unique(test_UP[ID == tiss]$LRI)
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

