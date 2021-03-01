
#library(skmeans)
#library(igraph)
#library(signnet)







###



sum(cci_pi_matrix_up==1)/100

test <- res@Seeddata
unique(as.vector(res@Seeddata))

cci_blok_up <- coclusterBinary(data = cci_pi_matrix_up, nbcocluster = c(10,4))

Heatmap(
  matrix = cci_pi_matrix_up,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE,
  #column_split = cci_blok_up@colclass,
  #row_split = cci_blok_up@rowclass,
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

table(rowSums(abs(cci_pi_matrix_binary_LR)))
cci_pi_matrix_binary_LR <- cci_pi_matrix_binary_LR[rowSums(abs(cci_pi_matrix_binary_LR)) >= 2, ]
table(colSums(abs(cci_pi_matrix_binary_LR)))
cci_pi_matrix_binary_LR <- cci_pi_matrix_binary_LR[, colSums(abs(cci_pi_matrix_binary_LR)) >= 1]
cci_pi_matrix_LR <- cci_pi_matrix_LR[rownames(cci_pi_matrix_binary_LR), colnames(cci_pi_matrix_binary_LR)]

strat <- coclusterStrategy (algo = "BEM")


cci_blok <- coclusterCategorical(data = cci_pi_matrix_binary_LR, nbcocluster = c(10,4))
#cci_blok_ct <- coclusterContinuous(data = cci_pi_matrix_LR, nbcocluster = c(7,4))
#cci_blok_ct

summary(cci_blok)
plot(cci_blok)
plot(cci_blok, type = "distribution")
table(cci_blok@colclass)
table(cci_blok@rowclass)

cci_blok@rowclass

cci_blok<- readRDS("../data_scAgeCom/analysis/outputs_data/blok_facs_50_30.rds")

Heatmap(
  matrix = cci_pi_matrix_binary_LR,
  name = "test",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = cci_blok@colclass,
  row_split = cci_blok@rowclass,
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

Heatmap(cci_pi_matrix_LR)
names(cci_blok@colclass) <- colnames(cci_pi_matrix_binary_LR)
names(cci_blok@rowclass) <- rownames(cci_pi_matrix_binary_LR)

binary_groups_means <- t(sapply(
  sort(unique(cci_blok@rowclass)),
  function(i) {
    sapply(
      sort(unique(cci_blok@colclass)),
      function(j) {
        rows_to_keep <- names(cci_blok@rowclass[cci_blok@rowclass == i])
        cols_to_keep <- names(cci_blok@colclass[cci_blok@colclass == j])
        mean(cci_pi_matrix_binary_LR[rows_to_keep, cols_to_keep])
      }
    )
  }
))

binary_groups_abs_means <- t(sapply(
  sort(unique(cci_blok@rowclass)),
  function(i) {
    sapply(
      sort(unique(cci_blok@colclass)),
      function(j) {
        rows_to_keep <- names(cci_blok@rowclass[cci_blok@rowclass == i])
        cols_to_keep <- names(cci_blok@colclass[cci_blok@colclass == j])
        mean(abs(cci_pi_matrix_binary_LR[rows_to_keep, cols_to_keep]))
      }
    )
  }
))

binary_groups_sd <- t(sapply(
  sort(unique(cci_blok@rowclass)),
  function(i) {
    sapply(
      sort(unique(cci_blok@colclass)),
      function(j) {
        rows_to_keep <- names(cci_blok@rowclass[cci_blok@rowclass == i])
        cols_to_keep <- names(cci_blok@colclass[cci_blok@colclass == j])
        sd(cci_pi_matrix_binary_LR[rows_to_keep, cols_to_keep])
      }
    )
  }
))
binary_groups_up<- t(sapply(
  sort(unique(cci_blok@rowclass)),
  function(i) {
    sapply(
      sort(unique(cci_blok@colclass)),
      function(j) {
        rows_to_keep <- names(cci_blok@rowclass[cci_blok@rowclass == i])
        cols_to_keep <- names(cci_blok@colclass[cci_blok@colclass == j])
        mat_temp <- cci_pi_matrix_binary_LR[rows_to_keep, cols_to_keep]
        sum(mat_temp == 1)/(nrow(mat_temp)*ncol(mat_temp))
      }
    )
  }
))
binary_groups_down<- t(sapply(
  sort(unique(cci_blok@rowclass)),
  function(i) {
    sapply(
      sort(unique(cci_blok@colclass)),
      function(j) {
        rows_to_keep <- names(cci_blok@rowclass[cci_blok@rowclass == i])
        cols_to_keep <- names(cci_blok@colclass[cci_blok@colclass == j])
        mat_temp <- cci_pi_matrix_binary_LR[rows_to_keep, cols_to_keep]
        sum(mat_temp == -1)/(nrow(mat_temp)*ncol(mat_temp))
      }
    )
  }
))

binary_groups_dt <- data.table(
  avg = as.vector(binary_groups_means),
  sd = as.vector(binary_groups_sd),
  avg_abs = as.vector(binary_groups_abs_means),
  up_pct = as.vector(binary_groups_up),
  down_pct = as.vector(binary_groups_down)
)

ggplot(binary_groups_dt, aes(avg, sd)) + geom_point()
ggplot(binary_groups_dt, aes(avg_abs, sd)) + geom_point()
ggplot(binary_groups_dt, aes(avg_abs, avg)) + geom_point()
ggplot(binary_groups_dt, aes(sd, up_pct)) + geom_point()
ggplot(binary_groups_dt, aes(down_pct, up_pct)) + geom_point()

bin_sd_sorted <- sort(as.vector(binary_groups_sd), decreasing = TRUE)
bin_up_sorted <- sort(as.vector(binary_groups_up), decreasing = TRUE)
bin_down_sorted <- sort(as.vector(binary_groups_down), decreasing = TRUE)

which(binary_groups_sd == bin_sd_sorted[[2]], arr.ind = TRUE)
which(binary_groups_up == bin_up_sorted[[1]], arr.ind = TRUE)
which(binary_groups_down == bin_down_sorted[[1]], arr.ind = TRUE)

hist(bin_up_sorted, breaks = 50)

hist(binary_groups_sd, breaks = 100)

cci_blok@rowclass[cci_blok@rowclass == 1]
cci_blok@colclass[cci_blok@colclass == 0]
test <- cci_pi_matrix_binary_LR[names(cci_blok@rowclass[cci_blok@rowclass == 0]),
                                names(cci_blok@colclass[cci_blok@colclass == 3])]

Heatmap(
  test,
  col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "gray", "red"))
)

sort(table(sapply(names(cci_blok@colclass[cci_blok@colclass == 6]),
             function(i) {
               unique(cci_dt[ID_ER == i]$ER_CELL_FAMILY)
             }
)
), decreasing = TRUE)
sort(table(sapply(names(cci_blok@colclass[cci_blok@colclass == 6]),
             function(i) {
               unique(cci_dt[ID_ER == i]$ID)
             }
)
), decreasing = TRUE)

sum(cci_pi_matrix_binary_LR == 1)
sum(cci_pi_matrix_binary_LR == -1)
sum(cci_pi_matrix_binary_LR == 0)
data(categoricaldata)





ratio_up_down <- rowSums(cci_pi_matrix_binary_LR == 1)/rowSums(cci_pi_matrix_binary_LR == -1)
n_up_down <- rowSums(abs(cci_pi_matrix_binary_LR) == 1)

head(sort(rowSums(cci_pi_matrix_binary_LR == -1), decreasing = TRUE))

plot(log(n_up_down), log(ratio_up_down))

sum(is.infinite(ratio_up_down))

hist(log(test), breaks = 100)
plot(rowSums(cci_pi_matrix_binary_LR == 1), rowSums(cci_pi_matrix_binary_LR == -1))

library(fabia)

fab_modules <- fabia(cci_pi_matrix_binary_LR, 5, norm = 0, center = 0, alpha = 0.001)
summary(fab_modules)
show(fab_modules)

rb <- extractBic(fab_modules)
names(rb$bic[1,]$bixv)
rb$bic[4,]$bixn
rb$bic[4,]$biypn

i <- 2

Heatmap(
  cci_pi_matrix_binary_LR[rb$bic[i,]$bixn, rb$bic[i,]$biypn],
  col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "gray", "red"))
)


library(isa2)

modules <- isa(cci_pi_matrix_binary_LR)
summary(modules)

modules <- isa(cci_pi_matrix_LR)

mymodules <- lapply(seq(ncol(modules$rows)), function(x) {
  list(rows = which(modules$rows[, x] != 0),
       columns = which(modules$columns[, x] !=
                         0))
})

i <- 5
Heatmap(
  cci_pi_matrix_LR[mymodules[[i]]$rows, mymodules[[i]]$columns]
)


Heatmap(
  cci_pi_matrix_binary_LR[mymodules[[i]]$rows, mymodules[[i]]$columns],
  col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "gray", "red"))
)

isa_summary <- t(sapply(
  1:length(mymodules),
  function(i) {
    mat_temp <- cci_pi_matrix_binary_LR[mymodules[[i]]$rows, mymodules[[i]]$columns]
    return(
      c(
        mean(mat_temp),
        mean(abs(mat_temp)),
        sd(mat_temp),
        sum(mat_temp == 1),
        sum(mat_temp == -1)
      )
    )
  }
)
)

test_norm <- isa.normalize(cci_pi_matrix_binary_LR)

Heatmap(test_norm[[1]], show_column_names = FALSE, show_row_names = FALSE, use_raster = TRUE)

#without normalization

data_isa <- cci_pi_matrix_binary_LR + 1
normed.data <- list(Er = t(data_isa), Ec = data_isa)

row.seeds <- generate.seeds(length = nrow(data_isa),
                            count = 1000, sparsity = c(2, 5, 10, 100))

cols.seeds <- generate.seeds(length = nrow(t(data_isa)),
                            count = 1000, sparsity = c(2, 5, 10))

modules <- isa.iterate(normed.data = normed.data, row.seeds = row.seeds,
                             col.seeds = cols.seeds,
                             thr.row = 1, thr.col = 1, direction = "updown")
ncol(modules$rows)
sum(apply(modules$rows == 0, 2, all))
modules2 <- isa.unique(normed.data, modules, cor.limit = 0.9)
ncol(modules2$rows)
rob <- robustness(normed.data, modules2$rows,
                  modules2$columns)
summary(rob)
par(cex.lab = 1.5, cex.axis = 1.5)
boxplot(rob, ylab = "Robustness")
modules3 <- isa.filter.robust(data_isa, normed.data,
                              modules2, perms = 2, row.seeds = row.seeds)
ncol(modules3$rows)

mymodules2 <- lapply(seq(ncol(modules3$rows)), function(x) {
  list(rows = which(modules3$rows[, x] != 0),
       columns = which(modules3$columns[, x] !=
                         0))
})

plot(modules3$seeddata$rob, ylab = "Robustness",
     cex.lab = 1.5)

i <- 56
Heatmap(
  cci_pi_matrix_binary_LR[mymodules2[[i]]$rows, mymodules2[[i]]$columns],
  col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "gray", "red"))
)

isa_summary <- t(sapply(
  1:length(mymodules2),
  function(i) {
    mat_temp <- cci_pi_matrix_binary_LR[mymodules2[[i]]$rows, mymodules2[[i]]$columns]
    return(
      c(
        mean(mat_temp),
        mean(abs(mat_temp)),
        sd(mat_temp),
        sum(mat_temp == 1),
        sum(mat_temp == -1)
      )
    )
  }
)
)

Heatmap(cci_pi_matrix_binary_LR)

#######
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

cci_blok <- coclusterCategorical(data = cci_pi_matrix_binary_LR, nbcocluster = c(10,5))
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

#biclust


LE_test <- names(NumberxCol(res)[3,][ NumberxCol(res)[3,] == TRUE])
LE_test_split <- strsplit(LE_test, "_")
E_test <- paste0(sapply(LE_test_split, `[[`, 2), "_L")
L_test <- sapply(LE_test_split, function(i) {paste0(i[c(3,4)], collapse = "_")})

LE_test_mat <- matrix(c(E_test, L_test), nrow = length(E_test), ncol = 2)


RR_test <- names(RowxNumber(res)[,3][RowxNumber(res)[,3] == TRUE])
RR_test_split <- strsplit(RR_test, "_")
Rc_test <- paste0(sapply(RR_test_split, `[[`, 2), "_R")
Rg_test <- sapply(RR_test_split, function(i) {paste0(i[c(3,4, 5)], collapse = "_")})

RR_test_mat <- matrix(c(Rc_test, Rg_test), nrow = length(Rc_test), ncol = 2)

test_mat <- rbind(LE_test_mat, RR_test_mat)

test_hm <- cci_pi_matrix_up[RowxNumber(res)[,i], NumberxCol(res)[i,]]

apply(
  test_hm,
  MARGIN = c(1,2),
  function(i,j) {
    print(test_hm[i,j])
  }
)

library(igraph)

test_g <- graph_from_edgelist(test_mat, directed = TRUE)

test_g
plot(test_g)



unique(cci_dt$LR_SORTED)
test <- LRdb_human$LRdb_curated



my_db <- LRdb_mouse$LRdb_curated
my_db <- my_db[LR_SORTED %in% unique(cci_dt$LR_SORTED) ]



library(skmeans)
cci_skm_bin_LR <- skmeans(cci_pi_matrix_up, k = 15)$cluster
cci_skm_bin_ER <- skmeans(t(cci_pi_matrix_up), k = 8)$cluster


Heatmap(
  matrix = cci_pi_matrix_up,
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
