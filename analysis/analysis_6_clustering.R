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
## Perform clustering analysis on each tissue
## and on the global data to find patterns.
##
####################################################
##

## Libraries ####
library(scDiffCom)
library(data.table)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(iBBiG)
library(blockcluster)

#library(blockcluster)

## Load data ####
DATASETS_COMBINED <- readRDS("../data_scAgeCom/analysis/outputs_data/TMS_scAgeCom_processed_bestORA.rds")

## Extract detected CCIs #####

DATASETS_CCI <- lapply(
  DATASETS_COMBINED,
  function(i) {
    temp_dt <- i$dataset@cci_detected
    temp_dt[, ID_ER := paste(ID, ER_CELLTYPES, sep = "_")]
    temp_dt[, REGULATION_SCORE := ifelse(
      REGULATION == "UP",
      1,
      ifelse(
        REGULATION == "DOWN",
        -1,
        0
      )
    )
    ]
    #temp_dt[, PI_VALUE := -log(BH_P_VALUE_DE + 1E-4)*LOGFC]
    temp_dt[, EMITTER_CELL_GENES := paste(EMITTER_CELLTYPE, LIGAND_1, LIGAND_2, sep = "_")]
    temp_dt[, RECEIVER_CELL_GENES := paste(RECEIVER_CELLTYPE, RECEPTOR_1, RECEPTOR_2, RECEPTOR_3, sep = "_")]
    return(temp_dt)
  }
)

## Create Cell-types vs Genes matrices globally ####

DATASETS_MATRIX_CELLS_vs_GENES <- lapply(
  DATASETS_CCI,
  function(i) {
    temp_mat <- as.matrix(
      dcast.data.table(
        i[, c("LR_GENES", "ID_ER", "REGULATION_SCORE") ],
        formula = LR_GENES ~ ID_ER,
        value.var = "REGULATION_SCORE"
      ),
      rownames = "LR_GENES"
    )
    temp_mat[is.na(temp_mat)] <- 0
    temp_mat_up <- temp_mat
    temp_mat_up[temp_mat_up == -1] <- 0
    temp_mat_down <- temp_mat
    temp_mat_down[temp_mat_down == 1] <- 0
    temp_mat_down[temp_mat_down == -1] <- 1
    temp_mat <- temp_mat[rowSums(abs(temp_mat)) >= 1, ]
    temp_mat <- temp_mat[, colSums(abs(temp_mat)) >= 1]
    temp_mat_up <- temp_mat_up[rowSums(temp_mat_up) >= 1, ]
    temp_mat_up <- temp_mat_up[, colSums(temp_mat_up) >= 1]
    temp_mat_down <- temp_mat_down[rowSums(temp_mat_down) >= 1, ]
    temp_mat_down <- temp_mat_down[, colSums(temp_mat_down) >= 1]
    return(
      list(
        matrix_up_down = temp_mat,
        matrix_up = temp_mat_up,
        matrix_down = temp_mat_down
      )
    )
  }
)

## Create Cell-types vs Genes matrices for each tissue ####

DATASETS_MATRIX_CELLS_vs_GENES_PER_ID <- lapply(
  DATASETS_CCI,
  function(i) {
    ids <- sort(unique(i$ID))
    temp_all <- lapply(
      ids,
      function(id) {
        temp_mat <- as.matrix(
          dcast.data.table(
            i[ID == id, c("LR_GENES", "ID_ER", "REGULATION_SCORE") ],
            formula = LR_GENES ~ ID_ER,
            value.var = "REGULATION_SCORE"
          ),
          rownames = "LR_GENES"
        )
        temp_mat[is.na(temp_mat)] <- 0
        temp_mat_up <- temp_mat
        temp_mat_up[temp_mat_up == -1] <- 0
        temp_mat_down <- temp_mat
        temp_mat_down[temp_mat_down == 1] <- 0
        temp_mat_down[temp_mat_down == -1] <- 1
        temp_mat <- temp_mat[rowSums(abs(temp_mat)) >= 1, , drop = FALSE]
        temp_mat <- temp_mat[, colSums(abs(temp_mat)) >= 1, drop = FALSE]
        temp_mat_up <- temp_mat_up[rowSums(temp_mat_up) >= 1, , drop = FALSE]
        temp_mat_up <- temp_mat_up[, colSums(temp_mat_up) >= 1, drop = FALSE]
        temp_mat_down <- temp_mat_down[rowSums(temp_mat_down) >= 1, , drop = FALSE]
        temp_mat_down <- temp_mat_down[, colSums(temp_mat_down) >= 1, drop = FALSE]
        return(
          list(
            matrix_up_down = temp_mat,
            matrix_up = temp_mat_up,
            matrix_down = temp_mat_down
          )
        )
      }
    )
    names(temp_all) <- ids
    return(temp_all)
  }
)

## Create EMITTERS vs RECEIVERS matrices for each tissue ####

DATASETS_MATRIX_EMITTERS_vs_RECEIVERS <- lapply(
  DATASETS_CCI,
  function(i) {
    ids <- sort(unique(i$ID))
    temp_all <- lapply(
      ids,
      function(id) {
        print(id)
        temp_mat <- as.matrix(
          dcast.data.table(
            i[ID == id, c("RECEIVER_CELL_GENES", "EMITTER_CELL_GENES", "REGULATION_SCORE") ],
            formula = RECEIVER_CELL_GENES ~ EMITTER_CELL_GENES,
            value.var = "REGULATION_SCORE"
          ),
          rownames = "RECEIVER_CELL_GENES"
        )
        temp_mat[is.na(temp_mat)] <- 0
        temp_mat_up <- temp_mat
        temp_mat_up[temp_mat_up == -1] <- 0
        temp_mat_down <- temp_mat
        temp_mat_down[temp_mat_down == 1] <- 0
        temp_mat_down[temp_mat_down == -1] <- 1
        temp_mat <- temp_mat[rowSums(abs(temp_mat)) >= 1, , drop = FALSE]
        temp_mat <- temp_mat[, colSums(abs(temp_mat)) >= 1, drop = FALSE]
        temp_mat_up <- temp_mat_up[rowSums(temp_mat_up) >= 1, , drop = FALSE]
        temp_mat_up <- temp_mat_up[, colSums(temp_mat_up) >= 1, drop = FALSE]
        temp_mat_down <- temp_mat_down[rowSums(temp_mat_down) >= 1, , drop = FALSE]
        temp_mat_down <- temp_mat_down[, colSums(temp_mat_down) >= 1, drop = FALSE]
        return(
          list(
            matrix_up_down = temp_mat,
            matrix_up = temp_mat_up,
            matrix_down = temp_mat_down
          )
        )
      }
    )
    names(temp_all) <- ids
    return(temp_all)
  }
)

## Set seed ####
set.seed(42)

## We find modules with iBBiG ####

#see how to obtain weighed scores
#see how to remove bad clusters
#can choose a large number of clusters to start with
#how to rank stuff in a cluster

test_obj <- DATASETS_MATRIX_CELLS_vs_GENES_PER_ID$facs_female$SCAT$matrix_down

test <- iBBiG(
 test_obj,
  nModules = 20,
  alpha = 0.5
)

sapply(1:20,
       function(i) {
         dim( test_obj[RowxNumber(test)[,i], NumberxCol(test)[i,], drop = FALSE])}
)

test
sort(test@Clusterscores, decreasing = TRUE)
test@RowScorexNumber
RowScorexNumber(test)[1:2,]

i <- 2
Heatmap(
  (test_obj[RowxNumber(test)[,i], NumberxCol(test)[i,], drop = FALSE]),
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 9),
  col = colorRamp2(c(-1, 0, 1), c("red", "gray", "blue")),
  height = unit(4, "cm")
)


test <- iBBiG(
  binaryMatrix = DATASETS_MATRIX_EMITTERS_vs_RECEIVERS$facs_mixed$Diaphragm$matrix_up,
  nModules = 15,
  alpha = 0.4
)

test_subobj <- test_obj[RowxNumber(test)[,i], NumberxCol(test)[i,], drop = FALSE]
test_subLR <- data.table(LR_GENES = rownames(test_subobj))
test_subLR[, c("LIGAND", "RECEPTOR") := tstrsplit(LR_GENES, ":", fixed=TRUE)]
test_subLR2 <- dcast.data.table(test_subLR,
                                formula = LIGAND ~ RECEPTOR,
                                fun.aggregate = length)
test_subLR2 <- as.matrix(test_subLR2, rownames = "LIGAND")

test_subER <- data.table(ER = colnames(test_subobj))
test_subER[, c("TISSUE", "EM", "RE") := tstrsplit(ER, "_", fixed=TRUE)]
test_subER2 <- dcast.data.table(test_subER,
                                formula = EM ~ RE,
                                fun.aggregate = length)
test_subER2 <- as.matrix(test_subER2, rownames = "EM")

library(igraph)

LO <- layout_as_bipartite(test_G_LR)

test_G_LR <- graph_from_incidence_matrix(test_subLR2)
test_G_ER <- graph_from_incidence_matrix(test_subER2)

test_G
plot(test_G_LR, layout = LO[, 2:1], vertex.size = 6)

is_bipartite(test_G)


test2 <- iBBiG::Clusterscores(test)
iBBiG::summary(test)




sort(table(unlist(sapply(
  1:15,
  function(i) {
    names(RowxNumber(res)[,i][RowxNumber(res)[,i] == TRUE])
  }
))), decreasing = TRUE)
sort(table(unlist(sapply(
  1:15,
  function(i) {
    names(NumberxCol(res)[i,][NumberxCol(res)[i,] == TRUE])
  }
))), decreasing = TRUE)


names(RowxNumber(res)[,i][RowxNumber(res)[,i] == TRUE])

## We block-organize the heatmap with blockcluster ####

cci_blok_up_1 <- coclusterCategorical(data = DATASETS_MATRIX_CELLS_vs_GENES_PER_ID$facs_female$SCAT$matrix_up_down, nbcocluster = c(10,6))

Heatmap(
  matrix =DATASETS_MATRIX_CELLS_vs_GENES_PER_ID$facs_mixed$SCAT$matrix_up_down,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = cci_blok_up_1@colclass,
  row_split = cci_blok_up_1@rowclass,
  use_raster = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = TRUE,
  cluster_column_slices = TRUE,
  row_names_centered = FALSE,
  row_names_gp = grid::gpar(fontsize = 6),
  column_names_gp = grid::gpar(fontsize = 6)
)




### 

LRdb_kegg <- LRdb_mouse$LRdb_curated_KEGG

kegg_cholesterol <- LRdb_kegg[KEGG_ID == "path:mmu04150"]
LRG_kegg_cholesterol <- unique(kegg_cholesterol$LR_GENES)
LRG_kegg_cholesterol <- LRG_kegg_cholesterol[LRG_kegg_cholesterol %in% 
                                               rownames(DATASETS_MATRIX_CELLS_vs_GENES$facs_mixed$matrix_up_down)]

cholesterol_matrix <- DATASETS_MATRIX_CELLS_vs_GENES$facs_mixed$matrix_up_down[LRG_kegg_cholesterol, ]
cholesterol_matrix <- cholesterol_matrix[,colSums(abs(cholesterol_matrix)) > 0]
cholesterol_matrix <- cholesterol_matrix[rowSums(abs(cholesterol_matrix)) > 0, ]

cholesterol_matrix_up <- cholesterol_matrix
cholesterol_matrix_up[cholesterol_matrix_up == -1] <- 0
cholesterol_matrix_up <- cholesterol_matrix_up[, colSums(cholesterol_matrix_up) > 0]
cholesterol_matrix_up <- cholesterol_matrix_up[rowSums(cholesterol_matrix_up) > 0, ]

cholesterol_matrix_down <- cholesterol_matrix
cholesterol_matrix_down[cholesterol_matrix_down == 1] <- 0
cholesterol_matrix_down[cholesterol_matrix_down == -1] <- 1
cholesterol_matrix_down <- cholesterol_matrix_down[, colSums(cholesterol_matrix_down) > 0]
cholesterol_matrix_down <- cholesterol_matrix_down[rowSums(cholesterol_matrix_down) > 0, ]

sum(cholesterol_matrix == 1)
sum(cholesterol_matrix == -1)

Heatmap(cholesterol_matrix)
Heatmap(
  cholesterol_matrix_down,
  col = colorRamp2(c(-1, 0, 1), c("red", "gray", "blue"))
)
Heatmap(
  cholesterol_matrix_up,
  col = colorRamp2(c(-1, 0, 1), c("blue", "gray", "red"))
)

test2 <- blockcluster::coclusterCategorical(cholesterol_matrix, nbcocluster = c(5, 6))

Heatmap(
  matrix = cholesterol_matrix,
  name = "test2",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = test2@colclass,
  row_split = test2@rowclass,
  use_raster = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = TRUE,
  cluster_column_slices = TRUE,
  row_names_centered = FALSE,
  row_names_gp = grid::gpar(fontsize = 6),
  column_names_gp = grid::gpar(fontsize = 6)
)

colnames(cholesterol_matrix[, test2@colclass == 3 ])

cholesterol_ibbig_up <- iBBiG(
  cholesterol_matrix_up,
  nModules = 10,
  alpha = 0.5
)
plot(cholesterol_ibbig_up)
sort(cholesterol_ibbig_up@Clusterscores, decreasing = TRUE)
i <- 4
Heatmap(
  cholesterol_matrix_up[RowxNumber(cholesterol_ibbig_up)[,i], NumberxCol(cholesterol_ibbig_up)[i,], drop = FALSE],
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 9),
  col = colorRamp2(c(-1, 0, 1), c("blue", "gray", "red")),
  height = unit(4, "cm")
)

cholesterol_ibbig_down <- iBBiG(
  cholesterol_matrix_down,
  nModules = 10,
  alpha = 0.5
)
plot(cholesterol_ibbig_down)
sort(cholesterol_ibbig_down@Clusterscores, decreasing = TRUE)
i <- 1
Heatmap(
  cholesterol_matrix_down[RowxNumber(cholesterol_ibbig_down)[,i], NumberxCol(cholesterol_ibbig_down)[i,], drop = FALSE],
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 9),
  col = colorRamp2(c(0, 1), c("gray", "blue")),
  height = unit(4, "cm")
)

## sub network by GO or KEGG ####


LRdb_kegg <- LRdb_mouse$LRdb_curated_KEGG
LRdb_go <- LRdb_mouse$LRdb_curated_GO

cci_term <- DATASETS_CCI$facs_mixed[LR_GENES %in% LRdb_kegg[KEGG_ID == "path:mmu04150"]$LR_GENES ]
cci_term_up <- cci_term[REGULATION_SIMPLE == "UP"]
cci_term_down <- cci_term[REGULATION_SIMPLE == "DOWN"]

cci_term_up[, .N, by = c("ID", "ER_CELLTYPES")]
cci_term_down[, .N, by = c("ID", "ER_CELLTYPES")]

cci_term_down[ID == "Aorta", .N, by = "EMITTER_CELLTYPE"]
cci_term_down[ID == "Aorta", .N, by = "RECEIVER_CELLTYPE"]
