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





#from other notebook #################




## Clustering analysis ####

cci_dt <- DATASETS_COMBINED$droplet_male$dataset@cci_detected
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

hist(as.vector(cci_pi_matrix_LR), breaks = 100)
summary(as.vector(cci_pi_matrix_LR))

cci_pi_matrix_LR[is.infinite(cci_pi_matrix_LR) & cci_pi_matrix_LR > 0] <- max(cci_pi_matrix_LR[!is.infinite(cci_pi_matrix_LR)])
cci_pi_matrix_LR[is.infinite(cci_pi_matrix_LR) & cci_pi_matrix_LR < 0] <- min(cci_pi_matrix_LR[!is.infinite(cci_pi_matrix_LR)])

cci_pi_matrix_ER <- t(cci_pi_matrix_LR)

skm_LR <- skmeans(cci_pi_matrix_LR, k = 20)
skm_ER <- skmeans(cci_pi_matrix_ER, k = 20)

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

cci_pi_matrix_LR[1:2, 1:2]

names(skm_LR$cluster[skm_LR$cluster == 17])
names(skm_ER$cluster[skm_ER$cluster == 17])

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

cci_pi_matrix_LR[1:5, 1:5]

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


names(skm_bin_ER$cluster[skm_bin_ER$cluster == 4])
names(skm_bin_LR$cluster[skm_bin_LR$cluster == 1])

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





###############################

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
