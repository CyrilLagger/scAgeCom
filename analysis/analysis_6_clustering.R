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
library(skmeans)
library(igraph)
library(signnet)
library(blockcluster)

## Data loading ####
DATASETS_COMBINED_log15 <- readRDS("../data_scAgeCom/analysis/outputs_data/analysis4_DATASETS_COMBINED_log15_light.rds")

## Prepare matrix ####
cci_dt <- DATASETS_COMBINED_log15$facs_mixed@cci_detected
cci_dt[, ID_ER := paste(ID, ER_CELLTYPES, sep = "_")]
cci_dt[, PI_VALUE := -log(P_VALUE_DE+1E-4)*LOGFC]

cci_pi_matrix_ER <- as.matrix(
  dcast.data.table(
    cci_dt[ID == "Aorta", c("ID_ER", "LR_GENES", "PI_VALUE") ],
    formula = ID_ER ~ LR_GENES,
    value.var = "PI_VALUE"
  ),
  rownames = "ID_ER"
)
cci_pi_matrix_ER[is.na(cci_pi_matrix_ER)] <- 0
cci_pi_matrix_LR <- t(cci_pi_matrix_ER)

## Binary matrix ####
cci_pi_matrix_binary_LR <- cci_pi_matrix_LR
cci_pi_matrix_binary_LR[abs(cci_pi_matrix_binary_LR) < -log(0.05)*log(1.5)] <- 0
cci_pi_matrix_binary_LR[cci_pi_matrix_binary_LR >= -log(0.05)*log(1.5)] <- 1
cci_pi_matrix_binary_LR[cci_pi_matrix_binary_LR <= log(0.05)*log(1.5)] <- -1

## Remove rows/cols ####
table(rowSums(abs(cci_pi_matrix_binary_LR)))
cci_pi_matrix_binary_LR <- cci_pi_matrix_binary_LR[rowSums(abs(cci_pi_matrix_binary_LR)) >= 2, ]

table(colSums(abs(cci_pi_matrix_binary_LR)))
cci_pi_matrix_binary_LR <- cci_pi_matrix_binary_LR[, colSums(abs(cci_pi_matrix_binary_LR)) >= 1]

cci_pi_matrix_LR <- cci_pi_matrix_LR[rownames(cci_pi_matrix_binary_LR), colnames(cci_pi_matrix_binary_LR)]

## Set seed ###
set.seed(42)

## Standard kmeans ####

cci_km_LR <- kmeans(cci_pi_matrix_LR, centers = 15, iter.max = 20, nstart = 25)$cluster
cci_km_ER <- kmeans(t(cci_pi_matrix_LR), centers = 5, iter.max = 20, nstart = 25)$cluster
table(cci_km_LR)
table(cci_km_ER)
cci_km_bin_LR <- kmeans(cci_pi_matrix_binary_LR, centers = 15, iter.max = 20, nstart = 25)$cluster
cci_km_bin_ER <- kmeans(t(cci_pi_matrix_binary_LR), centers = 15, iter.max = 20, nstart = 25)$cluster
table(cci_km_bin_LR)
table(cci_km_bin_ER)

Heatmap(
  matrix = cci_pi_matrix_binary_LR,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = cci_km_bin_ER,
  row_split = cci_km_bin_LR,
  use_raster = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = TRUE,
  cluster_column_slices = TRUE,
  row_names_centered = FALSE,
  row_names_gp = grid::gpar(fontsize = 4),
  column_names_gp = grid::gpar(fontsize = 4)
)

## Spherical kmeans ####

cci_skm_LR <- skmeans(cci_pi_matrix_LR, k = 15)$cluster
cci_skm_ER <- skmeans(t(cci_pi_matrix_LR), k = 5)$cluster

cci_skm_bin_LR <- skmeans(cci_pi_matrix_binary_LR, k = 15)$cluster
cci_skm_bin_ER <- skmeans(t(cci_pi_matrix_binary_LR), k = 8)$cluster


Heatmap(
  matrix = cci_pi_matrix_binary_LR,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = cci_skm_bin_ER,
  row_split = cci_skm_bin_LR,
  use_raster = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = TRUE,
  cluster_column_slices = TRUE,
  row_names_centered = FALSE,
  row_names_gp = grid::gpar(fontsize = 4),
  column_names_gp = grid::gpar(fontsize = 4)
)


## blockcluster way ####

cci_blok <- coclusterCategorical(data = cci_pi_matrix_binary_LR, nbcocluster = c(10,6))
summary(cci_blok)
plot(cci_blok)
table(cci_blok@colclass)
table(cci_blok@rowclass)

cci_blok@rowclass

table(cci_blok@rowclass)

Heatmap(
  matrix = cci_pi_matrix_binary_LR,
  name = "test",
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_split = cci_blok@colclass,
  row_split = cci_blok@rowclass,
  use_raster = TRUE,
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


## igraph way ####

cci_g <- graph.incidence(abs(cci_pi_matrix_binary_LR), weighted = TRUE)

test <- E(cci_g)


get.incidence(cci_g)
plot(cci_g, layout=layout_as_bipartite)


skm_bin_LR <- skmeans(cci_pi_matrix_LR, k = 12)
skm_bin_ER <- skmeans(t(cci_pi_matrix_LR), k = 7)
table(skm_bin_LR$cluster)
table(skm_bin_ER$cluster)

Heatmap(
  matrix = cci_pi_matrix_binary_LR,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = skm_bin_ER$cluster,
  row_split = skm_bin_LR$cluster,
  use_raster = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = TRUE,
  cluster_column_slices = TRUE,
  row_names_centered = FALSE,
  row_names_gp = grid::gpar(fontsize = 4),
  column_names_gp = grid::gpar(fontsize = 4)
)

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

which(binary_groups_sd == bin_sd_sorted[[8]], arr.ind = TRUE)

hist(binary_groups_sd, breaks = 100)

skm_bin_LR$cluster[skm_bin_LR$cluster == 20]
skm_bin_ER$cluster[skm_bin_ER$cluster == 20]
test <- cci_pi_matrix_binary_LR[names(skm_bin_LR$cluster[skm_bin_LR$cluster == 20]), names(skm_bin_ER$cluster[skm_bin_ER$cluster == 10])]

Heatmap(
  test,
  row_km = 2,
  column_km = 2
)

table(sapply(names(skm_bin_ER$cluster[skm_bin_ER$cluster == 22]),
       function(i) {
         unique(cci_dt[ID_ER == i]$ER_CELL_FAMILY)
       }
)
)
table(sapply(names(skm_bin_ER$cluster[skm_bin_ER$cluster == 22]),
             function(i) {
               unique(cci_dt[ID_ER == i]$ID)
             }
)
)

#biclustering

library(isa2)
library(biclust)
thr.row <- 2.7
thr.col <- 1.4
set.seed(42) # to get the same results, always
modules <- isa(test, thr.col = 0.4, thr.row = 0.4)
modules
length(modules)
Bc <- isa.biclust(modules)
Bc

mymodules <- lapply(seq(ncol(modules$rows)), function(x) {
  list(rows = which(modules$rows[, x] != 0),
       columns = which(modules$columns[, x] !=
                         0))
})

Heatmap(
  test[mymodules[[2]]$rows, mymodules[[2]]$columns]
)

bin_fabia <- fabia(test, 2, 0.1, 500)
summary(bin_fabia)
show(bin_fabia)
extractPlot(bin_fabia ,ti="FABIA")
plot(bin_fabia)
rb <- extractBic(bin_fabia)
bin_fabia@avini
rb$bic[1,]
plotBicluster(rb,1)

